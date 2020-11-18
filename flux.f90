module flux_mod
  !
!#define PROFILE_ON
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer => debug_timer
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx, only : ncs,ng
  !
  use ovar,    only : force_nd_subroutines
  !
  implicit none
  !
  private
  !
  public :: interp_cv_to_flxpts
  public :: interp_dusp_to_flxpts
  public :: flux_gradient_at_solpts
  public :: discontinuous_gradient_at_solpts
  public :: correct_gradients_at_solpts_gpu
  public :: correct_gradients_br2_gpu
  public :: get_invs_common_fluxes_interior
  public :: get_visc_common_fluxes_interior
  public :: get_invs_common_fluxes_bc
  public :: get_visc_common_fluxes_bc
  public :: init_cbias
  !
#ifdef ENABLE_OLD_FLUX_SUBROUTINES
!!public :: flux_gradient_at_solpts_visc
!!public :: flux_gradient_at_solpts_visc_gpu
!!public :: correct_interface_flux_visc
!!public :: correct_gradients_at_solpts
!!public :: correct_gradients_br2
#endif
  !
  real(wp), save, allocatable :: cbias(:)
  !$ACC DECLARE CREATE(cbias)
  real(wp), save, allocatable :: face_jump(:,:)
  !$ACC DECLARE CREATE(face_jump)
  !
  real(wp), save, allocatable :: faceuqp(:,:,:)
  real(wp), save, allocatable :: facegpqp(:,:,:)
  real(wp), save, allocatable :: faceduqp(:,:,:,:)
  real(wp), save, allocatable :: flux_qdtpts(:,:)
  !
  integer, parameter :: Separate_All_Projections = 0
  integer, parameter :: Separate_BC_Projections = 1
  integer, parameter :: Combine_All_Projections = 2
  !
  interface write_array
    module procedure write_2d_array, write_3d_array
  end interface write_array
  !
!#define WRITE_DATA_FOR_PROFILING
#ifdef WRITE_DATA_FOR_PROFILING
  integer, save :: num_prof_vflx_outputs = 0
  integer, save :: num_prof_br2_outputs = 0
#endif
#ifdef TIME_MINIAPPS_FOR_GPU_COMPARISON
  public :: write_br2_timings
  integer,  save :: timings_counter = 0
  integer,  save :: timings_points(1:1000) = 0
  real(wp), save :: timings_beg(1:1000) = zero
  real(wp), save :: timings_end(1:1000) = zero
  integer,  save :: timings_counter1 = 0
  integer,  save :: timings_points1(1:1000) = 0
  real(wp), save :: timings_beg1(1:1000) = zero
  real(wp), save :: timings_end1(1:1000) = zero
#endif
  !
contains
!
!###############################################################################
!
#ifdef TIME_MINIAPPS_FOR_GPU_COMPARISON
subroutine write_br2_timings()
  !
  integer :: i,l,n,io,ierr,ncalls
  character(len=100) :: fname
  !
  write (fname,2) 1,mypnum
  open (newunit=io,file=fname,status="replace",action="write",iostat=ierr)
  l = maxval(timings_points1(1:timings_counter1))
  ncalls = count(timings_points1(1:timings_counter1) == l)
  write (io,10) 1,mypnum+1,ncpu,ncalls,l
  i = 0
  do n = 1,timings_counter1
    if (timings_points1(n) == l) then
      i = i + 1
      write (io,1) i,timings_end1(n)-timings_beg1(n)
    end if
  end do
  close (io,iostat=ierr)
  !
  write (fname,2) 2,mypnum
  open (newunit=io,file=fname,status="replace",action="write",iostat=ierr)
  l = maxval(timings_points(1:timings_counter))
  ncalls = count(timings_points(1:timings_counter) == l)
  write (io,10) 2,mypnum+1,ncpu,ncalls,l
  i = 0
  do n = 1,timings_counter
    if (timings_points(n) == l) then
      i = i + 1
      write (io,1) i,timings_end(n)-timings_beg(n)
    end if
  end do
  close (io,iostat=ierr)
  !
  1 format (i4,2x,es23.15)
  2 format ("br2_timings.k-",i1,".",i2.2,".dat")
  10 format ('VARIABLES = "Call #", "Wall Time"',/, &
             'ZONE T="br2, k-',i1,': CPU ',i2,' of ',i2, &
             '", I=',i0,', DATAPACKING=POINT',/, &
             'AUXDATA loop_size="'i0,'"')
end subroutine write_br2_timings
#endif
!
!###############################################################################
!
subroutine interp_cv_to_flxpts(iter,nst)
  !
  !.. Use Statements ..
  use geovar,          only : nr,n_bnd_flxpts,n_phys_bnd_flxpts,n_flxpts
  use flowvar,         only : faceusp_array,facegpsp
  use thermo_mod,      only : update_thermo_properties
  use thermo_mod,      only : update_mixture_density
  use transport_mod,   only : update_transport_properties
  use module_limiters, only : check_face_thermo_properties
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: iter,nst
  !
  !.. Local Scalars ..
  integer :: n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_flxpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (force_nd_subroutines) then
    !
    call interp_cv_to_flxpts_nd(1,nq)
    !
  else
    !
    if (nr == 3) then
      call interp_cv_to_flxpts_3d(iter,nst)
    else if (nr == 2) then
      call interp_cv_to_flxpts_2d
    else
      call interp_cv_to_flxpts_1d
    end if
    !
    if (ncs > 0) then
      call interp_cv_to_flxpts_nd(1,ncs)
    end if
    !
  end if
  !
  ! Update the thermodynamic and transport gas properties
  ! at all the flux points of the interior faces
  !
  n1 = n_bnd_flxpts+1
  n2 = n_flxpts
  !
  call update_thermo_properties(n1,n2,faceusp_array,facegpsp,caller_name=pname)
  !
  ! Check for any negative thermodynamic properties at the flux points
  !
 !call check_face_thermo_properties()
  !
!block
!  use flowvar, only : facexyz_array
!  use ovar, only : grid_scaling_factor,rhoref,aref,tref
!  use eqn_idx, only : npres,ntemp
!  integer :: l,m,n
!  real(wp) :: rho_inv
!  call mpi_barrier(MPI_COMM_WORLD,mpierr)
!  cpu_loop: &
!  do m = 0,ncpu-1
!    if (mypnum == m) then
!      if (any(facegpsp(n1:n2,ntemp,:) <= zero) .or. &
!          any(faceusp_array(n1:n2,nec,:) <= zero) .or. &
!          any(faceusp_array(n1:n2,nee,:) <= zero)) then
!        pts_loop: &
!        do n = n1,n2
!          if (any(facegpsp(n,ntemp,:) <= zero) .or. &
!              any(faceusp_array(n,nec,:) <= zero) .or. &
!              any(faceusp_array(n,nec,:) <= zero)) then
!            write (iout,1) mypnum,n, &
!                           facexyz_array(1,n)/grid_scaling_factor, &
!                           facexyz_array(2,n)/grid_scaling_factor, &
!                           facexyz_array(3,n)/grid_scaling_factor
!            flush (iout)
!            sides_loop: &
!            do l = 1,2
!              if (faceusp_array(n,nec,l) /= zero) then
!                rho_inv = one/abs(faceusp_array(n,nec,l))
!              else
!                rho_inv = one
!              end if
!              write (iout,2) mypnum,n,l, &
!                             faceusp_array(n,nec,l)*rhoref, &
!                             faceusp_array(n,nmx,l)*rho_inv*aref, &
!                             faceusp_array(n,nmy,l)*rho_inv*aref, &
!                             faceusp_array(n,nmz,l)*rho_inv*aref, &
!                             faceusp_array(n,nee,l)*rhoref*aref*aref, &
!                             facegpsp(n,npres,l)*rhoref*aref*aref, &
!                             facegpsp(n,ntemp,l)*tref
!              flush (iout)
!            end do sides_loop
!            flush (iout)
!          end if
!          flush (iout)
!        end do pts_loop
!        flush (iout)
!      end if
!    end if
!    flush (iout)
!    call mpi_barrier(MPI_COMM_WORLD,mpierr)
!    flush (iout)
!    call mpi_barrier(MPI_COMM_WORLD,mpierr)
!  end do cpu_loop
!  flush (iout)
!  call mpi_barrier(MPI_COMM_WORLD,mpierr)
!  flush (iout)
!  call mpi_barrier(MPI_COMM_WORLD,mpierr)
!  1 format (/,"[",i5,"] invalid T @ fp=",i0, &
!            " : [x,y,z] = [",f13.8,",",f13.8,",",f13.8,"]")
!  2 format ("[",i5,"] invalid T @ fp=",i0,", side=",i0,/, &
!            8x,"[r,u,v,w,rE] = [",f9.5,3(",",f10.4),",",f12.2,"]",/, &
!            8x,"[P,T] = [",f11.2,",",f9.3,"]")
!end block
  call update_transport_properties(n1,n2,faceusp_array,facegpsp)
  !
  ! We need to update the mixture density from interpolating the species
  ! densities to the flux points in case the interpolation has caused
  ! any of the species densities to become negative. This only needs to be
  ! done for the communication boundary faces, since these interpolated
  ! values will be sent to the adjacent partitions before they would otherwise
  ! be updated by calling update_thermo_properties in subroutine get_faceusp.
  ! This doesnt need to be done for the flux points of the interior faces
  ! since the above call to update_thermo_properties also calls
  ! update_mixture_density.
  !
  n1 = n_phys_bnd_flxpts+1
  n2 = n_bnd_flxpts
  !
  call update_mixture_density(n1,n2,faceusp_array(:,:,1))
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine interp_cv_to_flxpts
!
!###############################################################################
!
subroutine interp_cv_to_flxpts_3d(iter,nst)
  !
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use flowvar,           only : usp
  use geovar,            only : face
  use flowvar,           only : faceusp_array
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: iter,nst
  !
  !.. Local Scalars ..
  integer :: n,nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,face_ec,face_mx,face_my,face_mz,face_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_flxpts_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp_array) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      !
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      !$ACC LOOP VECTOR
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        face_ec = zero
        face_mx = zero
        face_my = zero
        face_mz = zero
        face_ee = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          face_ec = face_ec + a * usp(n0+j,nec)
          face_mx = face_mx + a * usp(n0+j,nmx)
          face_my = face_my + a * usp(n0+j,nmy)
          face_mz = face_mz + a * usp(n0+j,nmz)
          face_ee = face_ee + a * usp(n0+j,nee)
          !
        end do
        !
        faceusp_array(f0+kf,nec,sf) = face_ec
        faceusp_array(f0+kf,nmx,sf) = face_mx
        faceusp_array(f0+kf,nmy,sf) = face_my
        faceusp_array(f0+kf,nmz,sf) = face_mz
        faceusp_array(f0+kf,nee,sf) = face_ee
        !
      end do
      !
    end do
    !
  end do
  !
 !if (present(iter) .and. present(nst)) then
 !  if (iter > 0 .and. nst > 0) then
 !    call write_array(100,iter,nst,faceusp_array(:,:,1))
 !  end if
 !end if
  !
end subroutine interp_cv_to_flxpts_3d
!
!###############################################################################
!
subroutine interp_cv_to_flxpts_2d()
  !
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use flowvar,           only : usp
  use geovar,            only : face
  use flowvar,           only : faceusp_array
  !
  !.. Local Scalars ..
  integer :: n,nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,face_ec,face_mx,face_my,face_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_flxpts_2d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp_array) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      !
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      !$ACC LOOP VECTOR
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        face_ec = zero
        face_mx = zero
        face_my = zero
        face_ee = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          face_ec = face_ec + a * usp(n0+j,nec)
          face_mx = face_mx + a * usp(n0+j,nmx)
          face_my = face_my + a * usp(n0+j,nmy)
          face_ee = face_ee + a * usp(n0+j,nee)
          !
        end do
        !
        faceusp_array(f0+kf,nec,sf) = face_ec
        faceusp_array(f0+kf,nmx,sf) = face_mx
        faceusp_array(f0+kf,nmy,sf) = face_my
        faceusp_array(f0+kf,nee,sf) = face_ee
        !
      end do
      !
    end do
    !
  end do
  !
end subroutine interp_cv_to_flxpts_2d
!
!###############################################################################
!
subroutine interp_cv_to_flxpts_1d()
  !
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use flowvar,           only : usp
  use geovar,            only : face
  use flowvar,           only : faceusp_array
  !
  !.. Local Scalars ..
  integer :: n,nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,face_ec,face_mx,face_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_flxpts_1d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp_array) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      !
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      !$ACC LOOP VECTOR
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        face_ec = zero
        face_mx = zero
        face_ee = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          face_ec = face_ec + a * usp(n0+j,nec)
          face_mx = face_mx + a * usp(n0+j,nmx)
          face_ee = face_ee + a * usp(n0+j,nee)
          !
        end do
        !
        faceusp_array(f0+kf,nec,sf) = face_ec
        faceusp_array(f0+kf,nmx,sf) = face_mx
        faceusp_array(f0+kf,nee,sf) = face_ee
        !
      end do
      !
    end do
    !
  end do
  !
end subroutine interp_cv_to_flxpts_1d
!
!###############################################################################
!
subroutine interp_cv_to_flxpts_nd(mbeg,mend)
  !
  use order_mod,         only : maxSP,maxFP
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use flowvar,           only : usp
  use geovar,            only : face
  use flowvar,           only : faceusp_array
  !
  !.. Formal Arguments ..
  integer, intent(in) :: mbeg,mend
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,k,m,f0
  integer  :: nf,kf,lf,sf,i,j,nfp
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,face_rot
  integer  :: num_cell_faces
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  integer, dimension(1:maxFP) :: fpidx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_flxpts_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp_array) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRIVATE(fpidx)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1 ! beginning index of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do n = 1,num_cell_faces
      !
      nf = cell(nc)%face(n)%idx  ! index of the current face
      sf = cell(nc)%face(n)%side ! side of face that cell nc is on
      !
      face_geom = cell(nc)%face(n)%geom   ! geom  of current face
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      fp_offset = cell(nc)%face(n)%fp_offset
      face_rot = cell(nc)%face(n)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
      !$ACC LOOP VECTOR COLLAPSE(2)
      do kf = 1,nfp
        do m = mbeg,mend
          !
          ! kf: index of the flx-pt local to this face
          ! lf: index of flx-pt within interp%mat
          !
          lf = fpidx(kf)
          val = zero
          !
          do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1,&
                 interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
            !
            j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
            a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
            !
            val = val + a * usp(n0+j,m)
            !
          end do
          !
          faceusp_array(f0+kf,m,sf) = val
          !
        end do
        !
      end do
      !
    end do
    !
  end do
  !
end subroutine interp_cv_to_flxpts_nd
!
!###############################################################################
!
subroutine interp_dusp_to_flxpts()
  !
  use geovar, only : nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_dusp_to_flxpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call interp_dusp_to_flxpts_3d
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call interp_dusp_to_flxpts_2d
  else
    call interp_dusp_to_flxpts_nd
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine interp_dusp_to_flxpts
!
!###############################################################################
!
subroutine interp_dusp_to_flxpts_3d()
  !
  !.. Use Statements ..
  use geovar,            only : ncell,cell,face,fp
  use interpolation_mod, only : interp
  use flowvar,           only : dusp,facedusp_array
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,k,f0
  integer  :: nf,kf,lf,sf,i,j,nfp
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,face_rot
  integer  :: num_cell_faces
  real(wp) :: a
  real(wp) :: decdx,dmxdx,dmydx,dmzdx,deedx
  real(wp) :: decdy,dmxdy,dmydy,dmzdy,deedy
  real(wp) :: decdz,dmxdz,dmydz,dmzdz,deedz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_dusp_to_flxpts_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(facedusp_array) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(dusp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1 ! beginning index of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do n = 1,num_cell_faces
      !
      nf = cell(nc)%face(n)%idx  ! index of the current face
      sf = cell(nc)%face(n)%side ! side of face that cell nc is on
      !
      face_geom = cell(nc)%face(n)%geom   ! geom  of current face
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      fp_offset = cell(nc)%face(n)%fp_offset
      face_rot = cell(nc)%face(n)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      !$ACC LOOP VECTOR
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        decdx = zero ; decdy = zero ; decdz = zero
        dmxdx = zero ; dmxdy = zero ; dmxdz = zero
        dmydx = zero ; dmydy = zero ; dmydz = zero
        dmzdx = zero ; dmzdy = zero ; dmzdz = zero
        deedx = zero ; deedy = zero ; deedz = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          decdx = decdx + a * dusp(n0+j,1,nec)
          decdy = decdy + a * dusp(n0+j,2,nec)
          decdz = decdz + a * dusp(n0+j,3,nec)
          !
          dmxdx = dmxdx + a * dusp(n0+j,1,nmx)
          dmxdy = dmxdy + a * dusp(n0+j,2,nmx)
          dmxdz = dmxdz + a * dusp(n0+j,3,nmx)
          !
          dmydx = dmydx + a * dusp(n0+j,1,nmy)
          dmydy = dmydy + a * dusp(n0+j,2,nmy)
          dmydz = dmydz + a * dusp(n0+j,3,nmy)
          !
          dmzdx = dmzdx + a * dusp(n0+j,1,nmz)
          dmzdy = dmzdy + a * dusp(n0+j,2,nmz)
          dmzdz = dmzdz + a * dusp(n0+j,3,nmz)
          !
          deedx = deedx + a * dusp(n0+j,1,nee)
          deedy = deedy + a * dusp(n0+j,2,nee)
          deedz = deedz + a * dusp(n0+j,3,nee)
          !
        end do
        !
        facedusp_array(f0+kf,1,nec,sf) = decdx
        facedusp_array(f0+kf,2,nec,sf) = decdy
        facedusp_array(f0+kf,3,nec,sf) = decdz
        !
        facedusp_array(f0+kf,1,nmx,sf) = dmxdx
        facedusp_array(f0+kf,2,nmx,sf) = dmxdy
        facedusp_array(f0+kf,3,nmx,sf) = dmxdz
        !
        facedusp_array(f0+kf,1,nmy,sf) = dmydx
        facedusp_array(f0+kf,2,nmy,sf) = dmydy
        facedusp_array(f0+kf,3,nmy,sf) = dmydz
        !
        facedusp_array(f0+kf,1,nmz,sf) = dmzdx
        facedusp_array(f0+kf,2,nmz,sf) = dmzdy
        facedusp_array(f0+kf,3,nmz,sf) = dmzdz
        !
        facedusp_array(f0+kf,1,nee,sf) = deedx
        facedusp_array(f0+kf,2,nee,sf) = deedy
        facedusp_array(f0+kf,3,nee,sf) = deedz
        !
      end do
      !
    end do
    !
  end do
  !
end subroutine interp_dusp_to_flxpts_3d
!
!###############################################################################
!
subroutine interp_dusp_to_flxpts_2d()
  !
  !.. Use Statements ..
  use geovar,            only : ncell,cell,face,fp
  use interpolation_mod, only : interp
  use flowvar,           only : dusp,facedusp_array
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,k,f0
  integer  :: nf,kf,lf,sf,i,j,nfp
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,face_rot
  integer  :: num_cell_faces
  real(wp) :: a
  real(wp) :: decdx,dmxdx,dmydx,deedx
  real(wp) :: decdy,dmxdy,dmydy,deedy
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_dusp_to_flxpts_2d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(facedusp_array) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(dusp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1 ! beginning index of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do n = 1,num_cell_faces
      !
      nf = cell(nc)%face(n)%idx  ! index of the current face
      sf = cell(nc)%face(n)%side ! side of face that cell nc is on
      !
      face_geom = cell(nc)%face(n)%geom   ! geom  of current face
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      fp_offset = cell(nc)%face(n)%fp_offset
      face_rot = cell(nc)%face(n)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      !$ACC LOOP VECTOR
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        decdx = zero ; decdy = zero
        dmxdx = zero ; dmxdy = zero
        dmydx = zero ; dmydy = zero
        deedx = zero ; deedy = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          decdx = decdx + a * dusp(n0+j,1,nec)
          decdy = decdy + a * dusp(n0+j,2,nec)
          !
          dmxdx = dmxdx + a * dusp(n0+j,1,nmx)
          dmxdy = dmxdy + a * dusp(n0+j,2,nmx)
          !
          dmydx = dmydx + a * dusp(n0+j,1,nmy)
          dmydy = dmydy + a * dusp(n0+j,2,nmy)
          !
          deedx = deedx + a * dusp(n0+j,1,nee)
          deedy = deedy + a * dusp(n0+j,2,nee)
          !
        end do
        !
        facedusp_array(f0+kf,1,nec,sf) = decdx
        facedusp_array(f0+kf,2,nec,sf) = decdy
        !
        facedusp_array(f0+kf,1,nmx,sf) = dmxdx
        facedusp_array(f0+kf,2,nmx,sf) = dmxdy
        !
        facedusp_array(f0+kf,1,nmy,sf) = dmydx
        facedusp_array(f0+kf,2,nmy,sf) = dmydy
        !
        facedusp_array(f0+kf,1,nee,sf) = deedx
        facedusp_array(f0+kf,2,nee,sf) = deedy
        !
      end do
      !
    end do
    !
  end do
  !
end subroutine interp_dusp_to_flxpts_2d
!
!###############################################################################
!
subroutine interp_dusp_to_flxpts_nd()
  !
  !.. Use Statements ..
  use order_mod,         only : maxFP
  use geovar,            only : nr,ncell,cell,face,fp
  use interpolation_mod, only : interp
  use flowvar,           only : dusp,facedusp_array
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,k,l,m,f0
  integer  :: nf,kf,lf,sf,i,j,nfp
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,face_rot
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  integer, dimension(1:maxFP) :: fpidx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_dusp_to_flxpts_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(facedusp_array) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(dusp) &
  !$ACC&                   PRIVATE(fpidx)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1 ! beginning index of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    do n = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(n)%idx  ! index of the current face
      sf = cell(nc)%face(n)%side ! side of face that cell nc is on
      !
      face_geom = cell(nc)%face(n)%geom   ! geom  of current face
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      fp_offset = cell(nc)%face(n)%fp_offset
      face_rot = cell(nc)%face(n)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
      !$ACC LOOP VECTOR COLLAPSE(3)
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        do m = 1,nq
          do l = 1,nr
            !
            lf = fpidx(kf)
            val = zero
            !
            do i = interp(this_geom,this_order)% &
                             toFace(face_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)% &
                             toFace(face_order)%row_ptr(lf+1)
              !
              j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
              a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
              !
              val = val + a * dusp(n0+j,l,m)
              !
            end do
            !
            facedusp_array(f0+kf,l,m,sf) = val
            !
          end do
        end do
       !facedusp_array(f0+kf,1:nr,1:nq,sf) = interp(this_geom,this_order)% &
       !                                            toFace(face_order)% &
       !                                  dot_mat2( lf , tdusp(1:np,1:nr,1:nq) )
      end do
      !
    end do
    !
  end do
  !
end subroutine interp_dusp_to_flxpts_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Procedures for computing the divergence of inviscid and viscous fluxes
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine flux_gradient_at_solpts(iter,nst)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts
  use geovar,  only : using_quadrature
  use flowvar, only : usp,gpsp,dusp
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iter,nst
  !
  !.. Local Scalars ..
  integer :: l,m,n
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: ifluxes(:,:,:)
  real(wp), allocatable :: vfluxes(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "flux_gradient_at_solpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( ifluxes(1:n_solpts,1:nr,1:nq) )
  allocate ( vfluxes(1:n_solpts,1:nr,1:nq) )
  !$ACC DATA CREATE(ifluxes,vfluxes)
  !
  if (using_quadrature) then
    call compute_fluxes_at_qdtpts(usp,dusp,ifluxes)
  else
    call compute_fluxes_at_solpts(usp,gpsp,dusp,ifluxes,vfluxes)
  end if
  !
  ! Transform the fluxes from physical to computation space. The viscous
  ! fluxes currently stored in vfluxes are no longer needed, so store the
  ! transformed fluxes in computational space in the vfluxes array.
  !
  call transform_fluxes(ifluxes,vfluxes)
  !
  call compute_flux_divergence(vfluxes,iter,nst)
  !
  call add_flux_correction(ifluxes)
  !
  !$ACC END DATA
  deallocate ( vfluxes )
  deallocate ( ifluxes )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine flux_gradient_at_solpts
!
!###############################################################################
!
subroutine compute_fluxes_at_solpts(usp,gpsp,dusp,ifluxes,vfluxes)
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts
  use ovar,    only : governing_equations
  !
  use visc_flux_mod, only : visc_flx_cv_gpu
  use invs_flux_mod, only : invs_flx_cv_gpu
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: gpsp(:,:)
  real(wp),    intent(in) :: dusp(:,:,:)
  real(wp), intent(inout) :: ifluxes(:,:,:)
  real(wp), intent(inout) :: vfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer :: l,m,n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_fluxes_at_solpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call invs_flx_cv_gpu( usp , gpsp , ifluxes )
  !
  if (governing_equations == NavierStokes_Eqns) then
    !
#ifdef WRITE_DATA_FOR_PROFILING
    if (num_prof_vflx_outputs == 0) then
      if (mypnum == glb_root) then
        call output_vflx_data(n_solpts,usp,dusp,gpsp)
      end if
    end if
#endif
    !
!#define DEBUG_FLUXES_AND_RESIDUAL
    call visc_flx_cv_gpu( 1 , n_solpts , usp , gpsp , dusp , vfluxes &
#ifdef DEBUG_FLUXES_AND_RESIDUAL
                          , output_visc_fluxes=true &
#endif
                          )
    !
#ifdef WRITE_DATA_FOR_PROFILING
    if (num_prof_vflx_outputs == 0) then
      if (mypnum == glb_root) then
        call output_gold_vflx_data(n_solpts,vfluxes)
      end if
      num_prof_vflx_outputs = 1
    end if
#endif
    !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,l,m,n
    character(len=30) :: fname
    !
    write (fname,1) "invs_fluxes",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,n_solpts
      write (io,2) n,((ifluxes(n,l,m), l=1,2), m=1,nq)
    end do
    close (io)
    !
    write (fname,1) "visc_fluxes",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,n_solpts
      write (io,2) n,((vfluxes(n,l,m), l=1,2), m=1,nq)
    end do
    close (io)
    !
    1 format (a,".",i0,".dat")
    2 format (i7,2x,20es14.6)
    !
  end block
#endif
    !
    !$ACC PARALLEL LOOP COLLAPSE(3) &
    !$ACC&              PRESENT(ifluxes,vfluxes)
    do m = 1,nq
      do l = 1,nr
        do n = 1,n_solpts
          ifluxes(n,l,m) = ifluxes(n,l,m) - vfluxes(n,l,m)
        end do
      end do
    end do
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_fluxes_at_solpts
!
!###############################################################################
!
subroutine output_vflx_data(n_solpts,cv,dcvdx,gp)
  !
  !.. Use Statements ..
  use geovar,  only : nr
  use ovar,    only : gam,rgasref_nd
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n_solpts
  real(wp), intent(in) :: cv(:,:)
  real(wp), intent(in) :: gp(:,:)
  real(wp), intent(in) :: dcvdx(:,:,:)
  !
  !.. Local Scalars ..
  integer :: io,l,m,n,ierr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dpvdx(:,:,:)
  !
continue
  !
  !
  !
  open (newunit=io,file="vflx_cv.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) ((cv(n,m), n=1,n_solpts), m=1,nq)
  close (io)
  !
  dpvdx = init_array(n_solpts,nr,nq,zero) ! F2003 AUTO-REALLOC
  call grad_cv_to_grad_pv_gpu(1,n_solpts,cv,gp,dcvdx,dpvdx)
  open (newunit=io,file="vflx_dpvdx.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) (((dpvdx(n,l,m), n=1,n_solpts), l=1,nr), m=1,nq)
  close (io)
  !
  open (newunit=io,file="vflx_gp.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) ((gp(n,m), n=1,n_solpts), m=1,ng)
  close (io)
  !
  !
  !
  open (newunit=io,file="grad_ovar.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) gam,rgasref_nd
  close (io)
  !
  open (newunit=io,file="grad_cv.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) ((cv(n,m), n=1,n_solpts), m=1,nq)
  close (io)
  !
  open (newunit=io,file="grad_dcvdx.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) (((dcvdx(n,l,m), n=1,n_solpts), l=1,nr), m=1,nq)
  close (io)
  !
  open (newunit=io,file="grad_gold_dpvdx.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) (((dpvdx(n,l,m), n=1,n_solpts), l=1,nr), m=1,nq)
  close (io)
  deallocate ( dpvdx , stat=ierr )
  !
end subroutine output_vflx_data
!
!###############################################################################
!
subroutine output_gold_vflx_data(n_solpts,vfluxes)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n_solpts
  real(wp), intent(in) :: vfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer :: io,l,m,n,ierr
  !
continue
  !
  !
  !
  open (newunit=io,file="vflx_gold_vfluxes.dat",status="replace", &
                   action="write",form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) (((vfluxes(n,l,m), n=1,n_solpts), l=1,nr), m=1,nq)
  close (io)
  !
end subroutine output_gold_vflx_data
!
!###############################################################################
!
subroutine compute_fluxes_at_qdtpts(usp,dusp,ifluxes)
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts,n_qdtpts
  use ovar,    only : governing_equations
  !
  use visc_flux_mod, only : visc_flx_cv_gpu
  use invs_flux_mod, only : invs_flx_cv_gpu
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: dusp(:,:,:)
  real(wp), intent(inout) :: ifluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer :: l,m,n
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: uqp(:,:)
  real(wp), allocatable :: gpqp(:,:)
  real(wp), allocatable :: duqp(:,:,:)
  real(wp), allocatable :: iflxqp(:,:,:)
  real(wp), allocatable :: vflxqp(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_fluxes_at_qdtpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( uqp(1:n_qdtpts,1:nq) )
  allocate ( gpqp(1:n_qdtpts,1:ng) )
  allocate ( iflxqp(1:n_qdtpts,1:nr,1:nq) )
  !$ACC ENTER DATA CREATE(uqp,gpqp,iflxqp)
  !
  ! Interpolate the conserved variables from the solution points
  ! to the quadrature points
  !
  call interp_cv_to_qdtpts( usp , uqp , gpqp )
  !
  ! Compute the inviscid fluxes at the quadrature points
  !
  call invs_flx_cv_gpu( uqp , gpqp , iflxqp )
  !
  if (governing_equations == NavierStokes_Eqns) then
    !
    ! Interpolate the derivatives of the conserved variables from
    ! the solution points to the quadrature points
    !
    allocate ( duqp(1:n_qdtpts,1:nr,1:nq) )
    allocate ( vflxqp(1:n_qdtpts,1:nr,1:nq) )
    !$ACC ENTER DATA CREATE(duqp,vflxqp)
    !
    call interp_dusp_to_qdtpts( dusp , duqp )
    !
    ! Compute the viscous fluxes at the quadrature points
    !
    call visc_flx_cv_gpu( 1 , n_qdtpts , uqp , gpqp , duqp , vflxqp )
    !
    ! Subtract the viscous fluxes from the inviscid fluxes to get
    ! the total flux at the quadrature points
    !
    !$ACC PARALLEL LOOP COLLAPSE(3) &
    !$ACC&              PRESENT(iflxqp,vflxqp)
    do m = 1,nq
      do l = 1,nr
        do n = 1,n_qdtpts
          iflxqp(n,l,m) = iflxqp(n,l,m) - vflxqp(n,l,m)
        end do
      end do
    end do
    !
    !$ACC EXIT DATA DELETE(duqp,vflxqp)
    deallocate ( vflxqp )
    deallocate ( duqp )
    !
  end if
  !
  !$ACC EXIT DATA DELETE(uqp,gpqp)
  deallocate ( gpqp )
  deallocate ( uqp )
  !
  ! Project the total flux at the quadrature points to the solution polynomial
  ! space and get the value of this polynomial at the solution points
  !
  call project_fluxes_to_solpts( iflxqp , ifluxes )
  !
  !$ACC EXIT DATA DELETE(iflxqp)
  deallocate ( iflxqp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_fluxes_at_qdtpts
!
!###############################################################################
!
subroutine interp_cv_to_qdtpts(usp,uqp,gpqp)
  !
  !.. Use Statements ..
  use geovar,        only : nr
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: uqp(:,:)
  real(wp), intent(inout) :: gpqp(:,:)
  !
  !.. Local Scalars ..
  integer :: n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_qdtpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call interp_cv_to_qdtpts_3d(usp,uqp)
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call interp_cv_to_qdtpts_2d(usp,uqp)
  else
    call interp_cv_to_qdtpts_nd(usp,uqp)
  end if
  !
  ! Update the thermodynamic and transport properties at the quadrature points
  !
  n1 = lbound(uqp,dim=1)
  n2 = ubound(uqp,dim=1)
  !
  call update_thermo_properties(n1,n2,uqp,gpqp,caller_name=pname)
  call update_transport_properties(n1,n2,uqp,gpqp)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine interp_cv_to_qdtpts
!
!###############################################################################
!
subroutine interp_cv_to_qdtpts_3d(usp,uqp)
  !
  use geovar,            only : ncell,cell
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: uqp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nc,n0,q0,nqp
  integer  :: this_geom,this_order,qdtr_order
  real(wp) :: a,cv_ec,cv_mx,cv_my,cv_mz,cv_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_qdtpts_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(uqp) &
  !$ACC&                   PRESENT(interp)
  do nc = 1,ncell
    !
    n0  = cell(nc)%beg_sp - 1
    !
    q0  = cell(nc)%beg_qp - 1
    nqp = cell(nc)%end_qp - q0
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    qdtr_order = cell(nc)%qdtr_order
    !
    ! Now interpolate the conserved variables from the solution points
    ! to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,nqp
      !
      cv_ec = zero
      cv_mx = zero
      cv_my = zero
      cv_mz = zero
      cv_ee = zero
      !
      do i = interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k)+1, &
             interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k+1)
        !
        j = interp(this_geom,this_order)%toCell(qdtr_order)%col_idx(i)
        a = interp(this_geom,this_order)%toCell(qdtr_order)%col_val(i)
        !
        cv_ec = cv_ec + a * usp(n0+j,nec)
        cv_mx = cv_mx + a * usp(n0+j,nmx)
        cv_my = cv_my + a * usp(n0+j,nmy)
        cv_mz = cv_mz + a * usp(n0+j,nmz)
        cv_ee = cv_ee + a * usp(n0+j,nee)
        !
      end do
      !
      uqp(q0+k,nec) = cv_ec
      uqp(q0+k,nmx) = cv_mx
      uqp(q0+k,nmy) = cv_my
      uqp(q0+k,nmz) = cv_mz
      uqp(q0+k,nee) = cv_ee
      !
    end do
    !
  end do
  !
end subroutine interp_cv_to_qdtpts_3d
!
!###############################################################################
!
subroutine interp_cv_to_qdtpts_2d(usp,uqp)
  !
  use geovar,            only : ncell,cell
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: uqp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nc,n0,q0,nqp
  integer  :: this_geom,this_order,qdtr_order
  real(wp) :: a,cv_ec,cv_mx,cv_my,cv_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_qdtpts_2d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(uqp) &
  !$ACC&                   PRESENT(interp)
  do nc = 1,ncell
    !
    n0  = cell(nc)%beg_sp - 1
    !
    q0  = cell(nc)%beg_qp - 1
    nqp = cell(nc)%end_qp - q0
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    qdtr_order = cell(nc)%qdtr_order
    !
    ! Now interpolate the conserved variables from the solution points
    ! to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,nqp
      !
      cv_ec = zero
      cv_mx = zero
      cv_my = zero
      cv_ee = zero
      !
      do i = interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k)+1, &
             interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k+1)
        !
        j = interp(this_geom,this_order)%toCell(qdtr_order)%col_idx(i)
        a = interp(this_geom,this_order)%toCell(qdtr_order)%col_val(i)
        !
        cv_ec = cv_ec + a * usp(n0+j,nec)
        cv_mx = cv_mx + a * usp(n0+j,nmx)
        cv_my = cv_my + a * usp(n0+j,nmy)
        cv_ee = cv_ee + a * usp(n0+j,nee)
        !
      end do
      !
      uqp(q0+k,nec) = cv_ec
      uqp(q0+k,nmx) = cv_mx
      uqp(q0+k,nmy) = cv_my
      uqp(q0+k,nee) = cv_ee
      !
    end do
    !
  end do
  !
end subroutine interp_cv_to_qdtpts_2d
!
!###############################################################################
!
subroutine interp_cv_to_qdtpts_nd(usp,uqp)
  !
  use geovar,            only : ncell,cell
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: uqp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,nc,n0,q0,nqp
  integer  :: this_geom,this_order,qdtr_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_qdtpts_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(uqp) &
  !$ACC&                   PRESENT(interp)
  do nc = 1,ncell
    !
    n0  = cell(nc)%beg_sp - 1
    !
    q0  = cell(nc)%beg_qp - 1
    nqp = cell(nc)%end_qp - q0
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    qdtr_order = cell(nc)%qdtr_order
    !
    ! Now interpolate the conserved variables from the solution points
    ! to the quadrature points
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do k = 1,nqp
      do m = 1,nq
        !
        val = zero
        !
        do i = interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k)+1, &
               interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k+1)
          !
          j = interp(this_geom,this_order)%toCell(qdtr_order)%col_idx(i)
          a = interp(this_geom,this_order)%toCell(qdtr_order)%col_val(i)
          !
          val = val + a * usp(n0+j,m)
          !
        end do
        !
        uqp(q0+k,m) = val
        !
      end do
      !
    end do
  end do
  !
end subroutine interp_cv_to_qdtpts_nd
!
!###############################################################################
!
subroutine interp_dusp_to_qdtpts(dusp,duqp)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: dusp(:,:,:)
  real(wp), intent(inout) :: duqp(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_dusp_to_qdtpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call interp_dusp_to_qdtpts_3d(dusp,duqp)
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call interp_dusp_to_qdtpts_2d(dusp,duqp)
  else
    call interp_dusp_to_qdtpts_nd(dusp,duqp)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine interp_dusp_to_qdtpts
!
!###############################################################################
!
subroutine interp_dusp_to_qdtpts_3d(dusp,duqp)
  !
  use geovar,            only : ncell,cell
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: dusp(:,:,:)
  real(wp), intent(inout) :: duqp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nc,n0,q0,nqp
  integer  :: this_geom,this_order,qdtr_order
  real(wp) :: a
  real(wp) :: decdx,dmxdx,dmydx,dmzdx,deedx
  real(wp) :: decdy,dmxdy,dmydy,dmzdy,deedy
  real(wp) :: decdz,dmxdz,dmydz,dmzdz,deedz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_dusp_to_qdtpts_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(dusp) &
  !$ACC&                   PRESENT(duqp) &
  !$ACC&                   PRESENT(interp)
  do nc = 1,ncell
    !
    n0  = cell(nc)%beg_sp - 1
    !
    q0  = cell(nc)%beg_qp - 1
    nqp = cell(nc)%end_qp - q0
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    qdtr_order = cell(nc)%qdtr_order
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,nqp
      !
      decdx = zero ; decdy = zero ; decdz = zero
      dmxdx = zero ; dmxdy = zero ; dmxdz = zero
      dmydx = zero ; dmydy = zero ; dmydz = zero
      dmzdx = zero ; dmzdy = zero ; dmzdz = zero
      deedx = zero ; deedy = zero ; deedz = zero
      !
      do i = interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k)+1, &
             interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k+1)
        !
        j = interp(this_geom,this_order)%toCell(qdtr_order)%col_idx(i)
        a = interp(this_geom,this_order)%toCell(qdtr_order)%col_val(i)
        !
        decdx = decdx + a * dusp(n0+j,1,nec)
        decdy = decdy + a * dusp(n0+j,2,nec)
        decdz = decdz + a * dusp(n0+j,3,nec)
        !
        dmxdx = dmxdx + a * dusp(n0+j,1,nmx)
        dmxdy = dmxdy + a * dusp(n0+j,2,nmx)
        dmxdz = dmxdz + a * dusp(n0+j,3,nmx)
        !
        dmydx = dmydx + a * dusp(n0+j,1,nmy)
        dmydy = dmydy + a * dusp(n0+j,2,nmy)
        dmydz = dmydz + a * dusp(n0+j,3,nmy)
        !
        dmzdx = dmzdx + a * dusp(n0+j,1,nmz)
        dmzdy = dmzdy + a * dusp(n0+j,2,nmz)
        dmzdz = dmzdz + a * dusp(n0+j,3,nmz)
        !
        deedx = deedx + a * dusp(n0+j,1,nee)
        deedy = deedy + a * dusp(n0+j,2,nee)
        deedz = deedz + a * dusp(n0+j,3,nee)
        !
      end do
      !
      duqp(q0+k,1,nec) = decdx
      duqp(q0+k,2,nec) = decdy
      duqp(q0+k,3,nec) = decdz
      !
      duqp(q0+k,1,nmx) = dmxdx
      duqp(q0+k,2,nmx) = dmxdy
      duqp(q0+k,3,nmx) = dmxdz
      !
      duqp(q0+k,1,nmy) = dmydx
      duqp(q0+k,2,nmy) = dmydy
      duqp(q0+k,3,nmy) = dmydz
      !
      duqp(q0+k,1,nmz) = dmzdx
      duqp(q0+k,2,nmz) = dmzdy
      duqp(q0+k,3,nmz) = dmzdz
      !
      duqp(q0+k,1,nee) = deedx
      duqp(q0+k,2,nee) = deedy
      duqp(q0+k,3,nee) = deedz
      !
    end do
    !
  end do
  !
end subroutine interp_dusp_to_qdtpts_3d
!
!###############################################################################
!
subroutine interp_dusp_to_qdtpts_2d(dusp,duqp)
  !
  use geovar,            only : ncell,cell
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: dusp(:,:,:)
  real(wp), intent(inout) :: duqp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nc,n0,q0,nqp
  integer  :: this_geom,this_order,qdtr_order
  real(wp) :: a
  real(wp) :: decdx,dmxdx,dmydx,deedx
  real(wp) :: decdy,dmxdy,dmydy,deedy
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_dusp_to_qdtpts_2d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(dusp) &
  !$ACC&                   PRESENT(duqp) &
  !$ACC&                   PRESENT(interp)
  do nc = 1,ncell
    !
    n0  = cell(nc)%beg_sp - 1
    !
    q0  = cell(nc)%beg_qp - 1
    nqp = cell(nc)%end_qp - q0
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    qdtr_order = cell(nc)%qdtr_order
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,nqp
      !
      decdx = zero ; decdy = zero
      dmxdx = zero ; dmxdy = zero
      dmydx = zero ; dmydy = zero
      deedx = zero ; deedy = zero
      !
      do i = interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k)+1, &
             interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k+1)
        !
        j = interp(this_geom,this_order)%toCell(qdtr_order)%col_idx(i)
        a = interp(this_geom,this_order)%toCell(qdtr_order)%col_val(i)
        !
        decdx = decdx + a * dusp(n0+j,1,nec)
        decdy = decdy + a * dusp(n0+j,2,nec)
        !
        dmxdx = dmxdx + a * dusp(n0+j,1,nmx)
        dmxdy = dmxdy + a * dusp(n0+j,2,nmx)
        !
        dmydx = dmydx + a * dusp(n0+j,1,nmy)
        dmydy = dmydy + a * dusp(n0+j,2,nmy)
        !
        deedx = deedx + a * dusp(n0+j,1,nee)
        deedy = deedy + a * dusp(n0+j,2,nee)
        !
      end do
      !
      duqp(q0+k,1,nec) = decdx
      duqp(q0+k,2,nec) = decdy
      !
      duqp(q0+k,1,nmx) = dmxdx
      duqp(q0+k,2,nmx) = dmxdy
      !
      duqp(q0+k,1,nmy) = dmydx
      duqp(q0+k,2,nmy) = dmydy
      !
      duqp(q0+k,1,nee) = deedx
      duqp(q0+k,2,nee) = deedy
      !
    end do
    !
  end do
  !
end subroutine interp_dusp_to_qdtpts_2d
!
!###############################################################################
!
subroutine interp_dusp_to_qdtpts_nd(dusp,duqp)
  !
  use geovar,            only : nr,ncell,cell
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: dusp(:,:,:)
  real(wp), intent(inout) :: duqp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,m,nc,n0,q0,nqp
  integer  :: this_geom,this_order,qdtr_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_dusp_to_qdtpts_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(dusp) &
  !$ACC&                   PRESENT(duqp) &
  !$ACC&                   PRESENT(interp)
  do nc = 1,ncell
    !
    n0  = cell(nc)%beg_sp - 1
    !
    q0  = cell(nc)%beg_qp - 1
    nqp = cell(nc)%end_qp - q0
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    qdtr_order = cell(nc)%qdtr_order
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do k = 1,nqp
      do m = 1,nq
        do l = 1,nr
          !
          val = zero
          !
          do i = interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k)+1, &
                 interp(this_geom,this_order)%toCell(qdtr_order)%col_ptr(k+1)
            !
            j = interp(this_geom,this_order)%toCell(qdtr_order)%col_idx(i)
            a = interp(this_geom,this_order)%toCell(qdtr_order)%col_val(i)
            !
            val = val + a * dusp(n0+j,l,m)
            !
          end do
          !
          duqp(q0+k,l,m) = val
          !
        end do
      end do
    end do
    !
  end do
  !
end subroutine interp_dusp_to_qdtpts_nd
!
!###############################################################################
!
subroutine project_fluxes_to_solpts(iflxqp,ifluxes)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: iflxqp(:,:,:)
  real(wp), intent(inout) :: ifluxes(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_fluxes_to_solpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call project_fluxes_to_solpts_3d(iflxqp,ifluxes)
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call project_fluxes_to_solpts_2d(iflxqp,ifluxes)
  else
    call project_fluxes_to_solpts_nd(iflxqp,ifluxes)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine project_fluxes_to_solpts
!
!###############################################################################
!
subroutine project_fluxes_to_solpts_3d(iflxqp,ifluxes)
  !
  use geovar,         only : ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: iflxqp(:,:,:)
  real(wp), intent(inout) :: ifluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nc,n0,q0,np
  integer  :: this_geom,this_order,qdtr_order
  real(wp) :: a
  real(wp) :: xflx_ec,xflx_mx,xflx_my,xflx_mz,xflx_ee
  real(wp) :: yflx_ec,yflx_mx,yflx_my,yflx_mz,yflx_ee
  real(wp) :: zflx_ec,zflx_mx,zflx_my,zflx_mz,zflx_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_fluxes_to_solpts_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(iflxqp) &
  !$ACC&                   PRESENT(ifluxes) &
  !$ACC&                   PRESENT(projct)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    q0 = cell(nc)%beg_qp - 1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    qdtr_order = cell(nc)%qdtr_order
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,np
      !
      xflx_ec = zero ; yflx_ec = zero ; zflx_ec = zero
      xflx_mx = zero ; yflx_mx = zero ; zflx_mx = zero
      xflx_my = zero ; yflx_my = zero ; zflx_my = zero
      xflx_mz = zero ; yflx_mz = zero ; zflx_mz = zero
      xflx_ee = zero ; yflx_ee = zero ; zflx_ee = zero
      !
      do i = projct(this_geom,this_order)%CellOI(qdtr_order)%row_ptr(k)+1, &
             projct(this_geom,this_order)%CellOI(qdtr_order)%row_ptr(k+1)
        !
        j = projct(this_geom,this_order)%CellOI(qdtr_order)%row_idx(i)
        a = projct(this_geom,this_order)%CellOI(qdtr_order)%row_val(i)
        !
        xflx_ec = xflx_ec + a * iflxqp(q0+j,1,nec)
        yflx_ec = yflx_ec + a * iflxqp(q0+j,2,nec)
        zflx_ec = zflx_ec + a * iflxqp(q0+j,3,nec)
        !
        xflx_mx = xflx_mx + a * iflxqp(q0+j,1,nmx)
        yflx_mx = yflx_mx + a * iflxqp(q0+j,2,nmx)
        zflx_mx = zflx_mx + a * iflxqp(q0+j,3,nmx)
        !
        xflx_my = xflx_my + a * iflxqp(q0+j,1,nmy)
        yflx_my = yflx_my + a * iflxqp(q0+j,2,nmy)
        zflx_my = zflx_my + a * iflxqp(q0+j,3,nmy)
        !
        xflx_mz = xflx_mz + a * iflxqp(q0+j,1,nmz)
        yflx_mz = yflx_mz + a * iflxqp(q0+j,2,nmz)
        zflx_mz = zflx_mz + a * iflxqp(q0+j,3,nmz)
        !
        xflx_ee = xflx_ee + a * iflxqp(q0+j,1,nee)
        yflx_ee = yflx_ee + a * iflxqp(q0+j,2,nee)
        zflx_ee = zflx_ee + a * iflxqp(q0+j,3,nee)
        !
      end do
      !
      ifluxes(n0+k,1,nec) = xflx_ec
      ifluxes(n0+k,2,nec) = yflx_ec
      ifluxes(n0+k,3,nec) = zflx_ec
      !
      ifluxes(n0+k,1,nmx) = xflx_mx
      ifluxes(n0+k,2,nmx) = yflx_mx
      ifluxes(n0+k,3,nmx) = zflx_mx
      !
      ifluxes(n0+k,1,nmy) = xflx_my
      ifluxes(n0+k,2,nmy) = yflx_my
      ifluxes(n0+k,3,nmy) = zflx_my
      !
      ifluxes(n0+k,1,nmz) = xflx_mz
      ifluxes(n0+k,2,nmz) = yflx_mz
      ifluxes(n0+k,3,nmz) = zflx_mz
      !
      ifluxes(n0+k,1,nee) = xflx_ee
      ifluxes(n0+k,2,nee) = yflx_ee
      ifluxes(n0+k,3,nee) = zflx_ee
      !
    end do
    !
  end do
  !
end subroutine project_fluxes_to_solpts_3d
!
!###############################################################################
!
subroutine project_fluxes_to_solpts_2d(iflxqp,ifluxes)
  !
  use geovar,         only : ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: iflxqp(:,:,:)
  real(wp), intent(inout) :: ifluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nc,n0,q0,np
  integer  :: this_geom,this_order,qdtr_order
  real(wp) :: a
  real(wp) :: xflx_ec,xflx_mx,xflx_my,xflx_ee
  real(wp) :: yflx_ec,yflx_mx,yflx_my,yflx_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_fluxes_to_solpts_2d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(iflxqp) &
  !$ACC&                   PRESENT(ifluxes) &
  !$ACC&                   PRESENT(projct)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    q0 = cell(nc)%beg_qp - 1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    qdtr_order = cell(nc)%qdtr_order
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,np
      !
      xflx_ec = zero ; yflx_ec = zero
      xflx_mx = zero ; yflx_mx = zero
      xflx_my = zero ; yflx_my = zero
      xflx_ee = zero ; yflx_ee = zero
      !
      do i = projct(this_geom,this_order)%CellOI(qdtr_order)%row_ptr(k)+1, &
             projct(this_geom,this_order)%CellOI(qdtr_order)%row_ptr(k+1)
        !
        j = projct(this_geom,this_order)%CellOI(qdtr_order)%row_idx(i)
        a = projct(this_geom,this_order)%CellOI(qdtr_order)%row_val(i)
        !
        xflx_ec = xflx_ec + a * iflxqp(q0+j,1,nec)
        yflx_ec = yflx_ec + a * iflxqp(q0+j,2,nec)
        !
        xflx_mx = xflx_mx + a * iflxqp(q0+j,1,nmx)
        yflx_mx = yflx_mx + a * iflxqp(q0+j,2,nmx)
        !
        xflx_my = xflx_my + a * iflxqp(q0+j,1,nmy)
        yflx_my = yflx_my + a * iflxqp(q0+j,2,nmy)
        !
        xflx_ee = xflx_ee + a * iflxqp(q0+j,1,nee)
        yflx_ee = yflx_ee + a * iflxqp(q0+j,2,nee)
        !
      end do
      !
      ifluxes(n0+k,1,nec) = xflx_ec
      ifluxes(n0+k,2,nec) = yflx_ec
      !
      ifluxes(n0+k,1,nmx) = xflx_mx
      ifluxes(n0+k,2,nmx) = yflx_mx
      !
      ifluxes(n0+k,1,nmy) = xflx_my
      ifluxes(n0+k,2,nmy) = yflx_my
      !
      ifluxes(n0+k,1,nee) = xflx_ee
      ifluxes(n0+k,2,nee) = yflx_ee
      !
    end do
    !
  end do
  !
end subroutine project_fluxes_to_solpts_2d
!
!###############################################################################
!
subroutine project_fluxes_to_solpts_nd(iflxqp,ifluxes)
  !
  use geovar,         only : nr,ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: iflxqp(:,:,:)
  real(wp), intent(inout) :: ifluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,m,nc,n0,q0,np
  integer  :: this_geom,this_order,qdtr_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_fluxes_to_solpts_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(iflxqp) &
  !$ACC&                   PRESENT(ifluxes) &
  !$ACC&                   PRESENT(projct)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    q0 = cell(nc)%beg_qp - 1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    qdtr_order = cell(nc)%qdtr_order
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do k = 1,np
      do m = 1,nq
        do l = 1,nr
          !
          val = zero
          !
          do i = projct(this_geom,this_order)%CellOI(qdtr_order)%row_ptr(k)+1, &
                 projct(this_geom,this_order)%CellOI(qdtr_order)%row_ptr(k+1)
            !
            j = projct(this_geom,this_order)%CellOI(qdtr_order)%row_idx(i)
            a = projct(this_geom,this_order)%CellOI(qdtr_order)%row_val(i)
            !
            val = val + a * iflxqp(q0+j,l,m)
            !
          end do
          !
          ifluxes(n0+k,l,m) = val
          !
        end do
      end do
    end do
    !
  end do
  !
end subroutine project_fluxes_to_solpts_nd
!
!###############################################################################
!
subroutine transform_fluxes(tfluxes,rfluxes)
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: tfluxes(:,:,:)
  real(wp), intent(inout) :: rfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "transform_fluxes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = size(tfluxes,dim=2)
  !
  if (force_nd_subroutines) then
    !
    call transform_fluxes_nd(1,nq,tfluxes,rfluxes)
    !
  else
    !
    if (nr == 3) then
      call transform_fluxes_3d(tfluxes,rfluxes)
    else if (nr == 2) then
      call transform_fluxes_2d(tfluxes,rfluxes)
    else
      call transform_fluxes_1d(tfluxes,rfluxes)
    end if
    !
    if (ncs > 0) then
      call transform_fluxes_nd(1,ncs,tfluxes,rfluxes)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine transform_fluxes
!
!###############################################################################
!
subroutine transform_fluxes_3d(tfluxes,rfluxes)
  !
  !.. Use Statements ..
  use metrics_mod, only : metrics_array,jac
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: tfluxes(:,:,:)
  real(wp), intent(inout) :: rfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,np
  real(wp) :: d1dx,d2dx,d3dx
  real(wp) :: d1dy,d2dy,d3dy
  real(wp) :: d1dz,d2dz,d3dz
  real(wp) :: xflx,yflx,zflx
  real(wp) :: jacobian
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "transform_fluxes_3d"
  !
continue
  !
  np = size(tfluxes,dim=1)
  !
  !$ACC PARALLEL LOOP PRESENT(tfluxes,rfluxes,jac,metrics_array)
  do n = 1,np
    !
    jacobian = jac(n)
    !
    d1dx = metrics_array(n,1,1)*jacobian
    d2dx = metrics_array(n,2,1)*jacobian
    d3dx = metrics_array(n,3,1)*jacobian
    !
    d1dy = metrics_array(n,1,2)*jacobian
    d2dy = metrics_array(n,2,2)*jacobian
    d3dy = metrics_array(n,3,2)*jacobian
    !
    d1dz = metrics_array(n,1,3)*jacobian
    d2dz = metrics_array(n,2,3)*jacobian
    d3dz = metrics_array(n,3,3)*jacobian
    !
    xflx = tfluxes(n,1,nec)
    yflx = tfluxes(n,2,nec)
    zflx = tfluxes(n,3,nec)
    !
    rfluxes(n,1,nec) = xflx*d1dx + yflx*d1dy + zflx*d1dz
    rfluxes(n,2,nec) = xflx*d2dx + yflx*d2dy + zflx*d2dz
    rfluxes(n,3,nec) = xflx*d3dx + yflx*d3dy + zflx*d3dz
    !
    xflx = tfluxes(n,1,nmx)
    yflx = tfluxes(n,2,nmx)
    zflx = tfluxes(n,3,nmx)
    !
    rfluxes(n,1,nmx) = xflx*d1dx + yflx*d1dy + zflx*d1dz
    rfluxes(n,2,nmx) = xflx*d2dx + yflx*d2dy + zflx*d2dz
    rfluxes(n,3,nmx) = xflx*d3dx + yflx*d3dy + zflx*d3dz
    !
    xflx = tfluxes(n,1,nmy)
    yflx = tfluxes(n,2,nmy)
    zflx = tfluxes(n,3,nmy)
    !
    rfluxes(n,1,nmy) = xflx*d1dx + yflx*d1dy + zflx*d1dz
    rfluxes(n,2,nmy) = xflx*d2dx + yflx*d2dy + zflx*d2dz
    rfluxes(n,3,nmy) = xflx*d3dx + yflx*d3dy + zflx*d3dz
    !
    xflx = tfluxes(n,1,nmz)
    yflx = tfluxes(n,2,nmz)
    zflx = tfluxes(n,3,nmz)
    !
    rfluxes(n,1,nmz) = xflx*d1dx + yflx*d1dy + zflx*d1dz
    rfluxes(n,2,nmz) = xflx*d2dx + yflx*d2dy + zflx*d2dz
    rfluxes(n,3,nmz) = xflx*d3dx + yflx*d3dy + zflx*d3dz
    !
    xflx = tfluxes(n,1,nee)
    yflx = tfluxes(n,2,nee)
    zflx = tfluxes(n,3,nee)
    !
    rfluxes(n,1,nee) = xflx*d1dx + yflx*d1dy + zflx*d1dz
    rfluxes(n,2,nee) = xflx*d2dx + yflx*d2dy + zflx*d2dz
    rfluxes(n,3,nee) = xflx*d3dx + yflx*d3dy + zflx*d3dz
    !
  end do
  !
end subroutine transform_fluxes_3d
!
!###############################################################################
!
subroutine transform_fluxes_2d(tfluxes,rfluxes)
  !
  !.. Use Statements ..
  use metrics_mod, only : metrics_array,jac
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: tfluxes(:,:,:)
  real(wp), intent(inout) :: rfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,np
  real(wp) :: d1dx,d2dx
  real(wp) :: d1dy,d2dy
  real(wp) :: xflx,yflx
  real(wp) :: jacobian
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "transform_fluxes_2d"
  !
continue
  !
  np = size(tfluxes,dim=1)
  !
  !$ACC PARALLEL LOOP PRESENT(tfluxes,rfluxes,jac,metrics_array)
  do n = 1,np
    !
    jacobian = jac(n)
    !
    d1dx = metrics_array(n,1,1)*jacobian
    d2dx = metrics_array(n,2,1)*jacobian
    !
    d1dy = metrics_array(n,1,2)*jacobian
    d2dy = metrics_array(n,2,2)*jacobian
    !
    xflx = tfluxes(n,1,nec)
    yflx = tfluxes(n,2,nec)
    !
    rfluxes(n,1,nec) = xflx*d1dx + yflx*d1dy
    rfluxes(n,2,nec) = xflx*d2dx + yflx*d2dy
    !
    xflx = tfluxes(n,1,nmx)
    yflx = tfluxes(n,2,nmx)
    !
    rfluxes(n,1,nmx) = xflx*d1dx + yflx*d1dy
    rfluxes(n,2,nmx) = xflx*d2dx + yflx*d2dy
    !
    xflx = tfluxes(n,1,nmy)
    yflx = tfluxes(n,2,nmy)
    !
    rfluxes(n,1,nmy) = xflx*d1dx + yflx*d1dy
    rfluxes(n,2,nmy) = xflx*d2dx + yflx*d2dy
    !
    xflx = tfluxes(n,1,nee)
    yflx = tfluxes(n,2,nee)
    !
    rfluxes(n,1,nee) = xflx*d1dx + yflx*d1dy
    rfluxes(n,2,nee) = xflx*d2dx + yflx*d2dy
    !
  end do
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,l,m,n
    character(len=30) :: fname
    !
    write (fname,1) "phys_fluxes",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,np
      write (io,2) n,((tfluxes(n,l,m), l=1,2), m=nec,nee)
    end do
    close (io)
    !
    write (fname,1) "comp_fluxes",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,np
      write (io,2) n,((rfluxes(n,l,m), l=1,2), m=nec,nee)
    end do
    close (io)
    !
    1 format (a,".",i0,".dat")
    2 format (i7,2x,20es14.6)
    !
  end block
#endif
  !
end subroutine transform_fluxes_2d
!
!###############################################################################
!
subroutine transform_fluxes_1d(tfluxes,rfluxes)
  !
  !.. Use Statements ..
  use metrics_mod, only : metrics_array,jac
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: tfluxes(:,:,:)
  real(wp), intent(inout) :: rfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,np
  real(wp) :: d1dx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "transform_fluxes_1d"
  !
continue
  !
  np = size(tfluxes,dim=1)
  !
  !$ACC PARALLEL LOOP PRESENT(tfluxes,rfluxes,jac,metrics_array)
  do n = 1,np
    !
    d1dx = metrics_array(n,1,1)*jac(n)
    !
    rfluxes(n,1,nec) = tfluxes(n,1,nec)*d1dx
    rfluxes(n,1,nmx) = tfluxes(n,1,nmx)*d1dx
    rfluxes(n,1,nee) = tfluxes(n,1,nee)*d1dx
    !
  end do
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,l,m,n
    character(len=30) :: fname
    !
    write (fname,1) "phys_fluxes",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,np
      write (io,2) n,((tfluxes(n,l,m), l=1,1), m=nec,nee)
    end do
    close (io)
    !
    write (fname,1) "comp_fluxes",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,np
      write (io,2) n,((rfluxes(n,l,m), l=1,1), m=nec,nee)
    end do
    close (io)
    !
    1 format (a,".",i0,".dat")
    2 format (i7,2x,20es14.6)
    !
  end block
#endif
  !
end subroutine transform_fluxes_1d
!
!###############################################################################
!
subroutine transform_fluxes_nd(mbeg,mend,tfluxes,rfluxes)
  !
  !.. Use Statements ..
  use metrics_mod, only : metrics_array,jac
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: mbeg,mend
  real(wp),    intent(in) :: tfluxes(:,:,:)
  real(wp), intent(inout) :: rfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,l,m,np,nr
  real(wp) :: val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "transform_fluxes_nd"
  !
continue
  !
  np = size(tfluxes,dim=1)
  nr = size(tfluxes,dim=2)
  !
  !$ACC LOOP VECTOR COLLAPSE(3)
  do m = mbeg,mend
    do l = 1,nr
      do k = 1,np
        val = zero
        do i = 1,nr
          val = val + tfluxes(k,i,m)*metrics_array(k,l,i)
        end do
        rfluxes(k,l,m) = val * jac(k)
      end do
    end do
  end do
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,l,m,n
    character(len=30) :: fname
    !
    write (fname,1) "phys_fluxes",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,np
      write (io,2) n,((tfluxes(n,l,m), l=1,nr), m=mbeg,mend)
    end do
    close (io)
    !
    write (fname,1) "comp_fluxes",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,np
      write (io,2) n,((rfluxes(n,l,m), l=1,nr), m=mbeg,mend)
    end do
    close (io)
    !
    1 format (a,".",i0,".dat")
    2 format (i7,2x,20es14.6)
    !
  end block
#endif
  !
end subroutine transform_fluxes_nd
!
!###############################################################################
!
subroutine compute_flux_divergence(rfluxes,iter,nst)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iter,nst
  real(wp), intent(in) :: rfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_flux_divergence"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = size(rfluxes,dim=2)
  !
  if (force_nd_subroutines) then
    !
    call compute_flux_divergence_nd(1,nq,rfluxes)
    !
  else
    !
    if (nr == 3) then
      call compute_flux_divergence_3d(rfluxes,iter,nst)
    else if (nr == 2) then
      call compute_flux_divergence_2d(rfluxes)
    else
      call compute_flux_divergence_1d(rfluxes)
    end if
    !
    if (ncs > 0) then
      call compute_flux_divergence_nd(1,ncs,rfluxes)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_flux_divergence
!
!###############################################################################
!
subroutine compute_flux_divergence_3d(rfluxes,iter,nst)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,       only : maxSP
  use geovar,          only : ncell,cell
  use derivatives_mod, only : deriv
  use flowvar,         only : residual
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iter,nst
  real(wp), intent(in) :: rfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,div_ec,div_mx,div_my,div_mz,div_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_flux_divergence_3d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&                   PRESENT(deriv) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(residual) &
  !$ACC&                   PRESENT(rfluxes)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the divergence of the discontinuous fluxes
    !
    !$ACC LOOP VECTOR
    do k = 1,np
      !
      div_ec = zero
      div_mx = zero
      div_my = zero
      div_mz = zero
      div_ee = zero
      !
      ! Add the flux divergence from the first computational
      ! direction to the running divergence sums
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(k+1)
        !
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        n = deriv(this_geom,this_order)%d(1)%row_idx(i)
        !
        div_ec = div_ec + a*rfluxes(n0+n,1,nec)
        div_mx = div_mx + a*rfluxes(n0+n,1,nmx)
        div_my = div_my + a*rfluxes(n0+n,1,nmy)
        div_mz = div_mz + a*rfluxes(n0+n,1,nmz)
        div_ee = div_ee + a*rfluxes(n0+n,1,nee)
        !
      end do
      !
      ! Add the flux divergence from the second computational
      ! direction to the running divergence sums
      !
      do i = deriv(this_geom,this_order)%d(2)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(2)%row_ptr(k+1)
        !
        a = deriv(this_geom,this_order)%d(2)%row_val(i)
        n = deriv(this_geom,this_order)%d(2)%row_idx(i)
        !
        div_ec = div_ec + a*rfluxes(n0+n,2,nec)
        div_mx = div_mx + a*rfluxes(n0+n,2,nmx)
        div_my = div_my + a*rfluxes(n0+n,2,nmy)
        div_mz = div_mz + a*rfluxes(n0+n,2,nmz)
        div_ee = div_ee + a*rfluxes(n0+n,2,nee)
        !
      end do
      !
      ! Add the flux divergence from the third computational
      ! direction to the running divergence sums
      !
      do i = deriv(this_geom,this_order)%d(3)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(3)%row_ptr(k+1)
        !
        a = deriv(this_geom,this_order)%d(3)%row_val(i)
        n = deriv(this_geom,this_order)%d(3)%row_idx(i)
        !
        div_ec = div_ec + a*rfluxes(n0+n,3,nec)
        div_mx = div_mx + a*rfluxes(n0+n,3,nmx)
        div_my = div_my + a*rfluxes(n0+n,3,nmy)
        div_mz = div_mz + a*rfluxes(n0+n,3,nmz)
        div_ee = div_ee + a*rfluxes(n0+n,3,nee)
        !
      end do
      !
      ! Save the flux divergence values to the RHS for solpt k
      !
      residual(n0+k,nec) = div_ec
      residual(n0+k,nmx) = div_mx
      residual(n0+k,nmy) = div_my
      residual(n0+k,nmz) = div_mz
      residual(n0+k,nee) = div_ee
      !
    end do
    !
  end do
  !
 !call write_array(100,iter,nst,residual)
  !
end subroutine compute_flux_divergence_3d
!
!###############################################################################
!
subroutine compute_flux_divergence_2d(rfluxes)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,       only : maxSP
  use geovar,          only : ncell,cell
  use derivatives_mod, only : deriv
  use flowvar,         only : residual
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,div_ec,div_mx,div_my,div_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_flux_divergence_2d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&                   PRESENT(deriv) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(residual) &
  !$ACC&                   PRESENT(rfluxes)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the divergence of the discontinuous fluxes
    !
    !$ACC LOOP VECTOR
    do k = 1,np
      !
      div_ec = zero
      div_mx = zero
      div_my = zero
      div_ee = zero
      !
      ! Add the flux divergence from the first computational
      ! direction to the running divergence sums
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(k+1)
        !
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        n = deriv(this_geom,this_order)%d(1)%row_idx(i)
        !
        div_ec = div_ec + a*rfluxes(n0+n,1,nec)
        div_mx = div_mx + a*rfluxes(n0+n,1,nmx)
        div_my = div_my + a*rfluxes(n0+n,1,nmy)
        div_ee = div_ee + a*rfluxes(n0+n,1,nee)
        !
      end do
      !
      ! Add the flux divergence from the second computational
      ! direction to the running divergence sums
      !
      do i = deriv(this_geom,this_order)%d(2)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(2)%row_ptr(k+1)
        !
        a = deriv(this_geom,this_order)%d(2)%row_val(i)
        n = deriv(this_geom,this_order)%d(2)%row_idx(i)
        !
        div_ec = div_ec + a*rfluxes(n0+n,2,nec)
        div_mx = div_mx + a*rfluxes(n0+n,2,nmx)
        div_my = div_my + a*rfluxes(n0+n,2,nmy)
        div_ee = div_ee + a*rfluxes(n0+n,2,nee)
        !
      end do
      !
      ! Save the flux divergence values to the RHS for solpt k
      !
      residual(n0+k,nec) = div_ec
      residual(n0+k,nmx) = div_mx
      residual(n0+k,nmy) = div_my
      residual(n0+k,nee) = div_ee
      !
    end do
    !
  end do
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,m,n
    character(len=20) :: fname
    !
    write (fname,1) itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,cell(ncell)%end_sp
      write (io,2) n,(residual(n,m), m=nec,nee)
    end do
    close (io)
    !
    1 format ("residual.",i0,".dat")
    2 format (i7,2x,10es18.10)
    !
  end block
#endif
  !
end subroutine compute_flux_divergence_2d
!
!###############################################################################
!
subroutine compute_flux_divergence_1d(rfluxes)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,       only : maxSP
  use geovar,          only : ncell,cell
  use derivatives_mod, only : deriv
  use flowvar,         only : residual
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,div_ec,div_mx,div_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_flux_divergence_1d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&                   PRESENT(deriv) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(residual) &
  !$ACC&                   PRESENT(rfluxes)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the divergence of the discontinuous fluxes
    !
    !$ACC LOOP VECTOR
    do k = 1,np
      !
      div_ec = zero
      div_mx = zero
      div_ee = zero
      !
      ! Add the flux divergence from the first computational
      ! direction to the running divergence sums
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(k+1)
        !
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        n = deriv(this_geom,this_order)%d(1)%row_idx(i)
        !
        div_ec = div_ec + a*rfluxes(n0+n,1,nec)
        div_mx = div_mx + a*rfluxes(n0+n,1,nmx)
        div_ee = div_ee + a*rfluxes(n0+n,1,nee)
        !
      end do
      !
      ! Save the flux divergence values to the RHS for solpt k
      !
      residual(n0+k,nec) = div_ec
      residual(n0+k,nmx) = div_mx
      residual(n0+k,nee) = div_ee
      !
    end do
    !
  end do
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,m,n
    character(len=20) :: fname
    !
    write (fname,1) itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,cell(ncell)%end_sp
      write (io,2) n,(residual(n,m), m=nec,nee)
    end do
    close (io)
    !
    1 format ("residual.",i0,".dat")
    2 format (i7,2x,10es18.10)
    !
  end block
#endif
  !
end subroutine compute_flux_divergence_1d
!
!###############################################################################
!
subroutine compute_flux_divergence_nd(mbeg,mend,rfluxes)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,       only : maxSP
  use geovar,          only : nr,ncell,cell
  use derivatives_mod, only : deriv
  use flowvar,         only : residual
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: mbeg,mend
  real(wp), intent(in) :: rfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,l,m,n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: rflx(1:maxSP,1:nr,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_flux_divergence_nd"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(rfluxes) &
  !$ACC&                   PRESENT(deriv) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(residual) &
  !$ACC&                   PRIVATE(rflx)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do m = mbeg,mend
      do l = 1,nr
        do k = 1,np
          rflx(k,l,m) = rfluxes(n0+k,l,m)
        end do
      end do
    end do
    !
    ! Compute the divergence of the discontinuous fluxes
    !
    ! Loop through each coordinate direction and add the divergence of
    ! the flux to the local residual array for the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do k = 1,np
      do m = mbeg,mend
        !
        ! Compute the divergence of the discontinuous fluxes for solpt k
        !
        val = zero
!       !$ACC LOOP SEQ
        do l = 1,nr
          do i = deriv(this_geom,this_order)%d(l)%row_ptr(k)+1, &
                 deriv(this_geom,this_order)%d(l)%row_ptr(k+1)
            !
            n = deriv(this_geom,this_order)%d(l)%row_idx(i)
            a = deriv(this_geom,this_order)%d(l)%row_val(i)
            !
            val = val + a * rflx(n,l,m)
            !
          end do
        end do
        !
        ! Save the divergence of the discontinuous fluxes to the RHS for solpt k
        !
        residual(n0+k,m) = val
        !
      end do
    end do
    !
  end do
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,m,n
    character(len=20) :: fname
    !
    write (fname,1) itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,cell(ncell)%end_sp
      write (io,2) n,(residual(n,m), m=mbeg,mend)
    end do
    close (io)
    !
    1 format ("residual.",i0,".dat")
    2 format (i7,2x,10es18.10)
    !
  end block
#endif
  !
end subroutine compute_flux_divergence_nd
!
!###############################################################################
!
subroutine add_flux_correction(tfluxes)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: tfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_flux_correction"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = size(tfluxes,dim=2)
  !
  if (force_nd_subroutines) then
    !
    call add_flux_correction_nd(1,nq,tfluxes)
    !
  else
    !
    if (nr == 3) then
      call add_flux_correction_3d(tfluxes)
    else if (nr == 2) then
      call add_flux_correction_2d(tfluxes)
    else
      call add_flux_correction_1d(tfluxes)
    end if
    !
    if (ncs > 0) then
      call add_flux_correction_nd(1,ncs,tfluxes)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine add_flux_correction
!
!###############################################################################
!
subroutine add_flux_correction_3d(tfluxes)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,         only : maxFP,geom_solpts
  use geovar,            only : ncell,cell,fp,face
  use interpolation_mod, only : interp
  use metrics_mod,       only : facenrm
  use flowvar,           only : residual,faceflx
  use correction_mod,    only : correct
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: tfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,nc,n0,np,kf,lf,f0
  integer  :: ncf,nfp,num_cell_faces
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,face_rot
  real(wp) :: flx_ec,flx_mx,flx_my,flx_mz,flx_ee
  real(wp) :: a,ds,rnx,rny,rnz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_flux_correction_3d"
  !
!#ifdef DEBUG_ON
!  logical(lk), save :: first_step = true
!#endif
continue
  !
 !do kf = lbound(geom_solpts,dim=2),ubound(geom_solpts,dim=2)
 !  do lf = lbound(geom_solpts,dim=1),ubound(geom_solpts,dim=1)
 !    write (*,*) lf,kf,geom_solpts(lf,kf)
 !  end do
 !end do
 !write (*,*)
 !!
 !!$ACC PARALLEL LOOP GANG NUM_GANGS(1)
 !do kf = lbound(geom_solpts,dim=2),ubound(geom_solpts,dim=2)
 !  !$ACC LOOP SEQ
 !  do lf = lbound(geom_solpts,dim=1),ubound(geom_solpts,dim=1)
 !    write (*,*) lf,kf,geom_solpts(lf,kf)
 !  end do
 !end do
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,"error")
!#ifdef DEBUG_ON
!  if (ncpu > 1) first_step = fals
!  if (ncell > 100) first_step = fals
!  if (first_step) then
!    write (100,2)
!  2 format ("nc  = cell number",/, &
!            "ncf = face of cell",/, &
!            "kf  = fp local to face of cell",/, &
!            "fp  = global fp",/, &
!            "lf  = fp of all faces of cell",/, &
!            "np  = sp local to cell",/, &
!            "sp  = global sp",/, &
!            "mv  = matrix value",/, &
!            " nc  ncf  kf    fp   lf   np    sp         mv")
!  end if
!#endif
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(correct) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(residual) &
  !$ACC&                   PRESENT(faceflx,facenrm,tfluxes) &
  !$ACC&                   PRESENT(cl)
  cell_loop: do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
!#ifdef DEBUG_ON
!  if (first_step) then
!    if (nc == 1) then
!      do n = 1,size(correct(this_geom,this_order)%mat,dim=1)
!        write (101,10) (correct(this_geom,this_order)%mat(n,i), &
!                        i=1,size(correct(this_geom,this_order)%mat,dim=2))
!      end do
!      face_order = cell(nc)%face(1)%order
!      do n = 1,size(interp(this_geom,this_order)%toFace(face_order)%mat,dim=1)
!        write (102,10) (interp(this_geom,this_order)%toFace(face_order)%mat(n,i), &
!                        i=1,size(interp(this_geom,this_order)%toFace(face_order)%mat,dim=2))
!      end do
!    end if
!  end if
!  10 format (1000es11.3)
!#endif
    !
    ! Interpolate the fluxes to the interface flux points and dot this
    ! interpolated vector with the face normal at each flux point
    ! NOTE: We dont need to exchange this information for communication and
    !       periodic boundary faces since this flux is only used within the
    !       current cell.
    !
    ! Loop over the faces of this cell
    num_cell_faces = size(cell(nc)%face)
    !
    faces_of_cell: do ncf = 1,num_cell_faces
      !
      ds = cl(this_geom) * merge( -one , one , cell(nc)%face(ncf)%side==1 )
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      nfp = geom_solpts(face_geom,face_order)
      f0 = face( cell(nc)%face(ncf)%idx )%beg_fp-1
      !
      !$ACC LOOP VECTOR
      pts_on_face: do kf = 1,nfp
        !
        rnx = facenrm(1,f0+kf)
        rny = facenrm(2,f0+kf)
        rnz = facenrm(3,f0+kf)
        !
        ! Initialize the interface flux values to the common interface flux
        ! NOTE: We initialize it to negative so that we can do a positive
        !       summation in the following loop. The sign of ds will make
        !       sure the flux values are the correct sign when using them
        !       to apply the correction to the interior solution points.
        !
        flx_ec = -faceflx(f0+kf,nec)
        flx_mx = -faceflx(f0+kf,nmx)
        flx_my = -faceflx(f0+kf,nmy)
        flx_mz = -faceflx(f0+kf,nmz)
        flx_ee = -faceflx(f0+kf,nee)
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        ! Get the normal interior discontinuous flux at this flux point
        !
       !flx_ec = rnx*a1*tx1 + rnx*a2*tx2 + rnx*a3*tx3 + rnx*a4*tx4 &
       !       + rny*a1*ty1 + rny*a2*ty2 + rny*a3*ty3 + rny*a4*ty4 &
       !       + rnz*a1*tz1 + rnz*a2*tz2 + rnz*a3*tz3 + rnz*a4*tz4
       !!
       !flx_ec = rnx*a1*tx1 + rny*a1*ty1 + rnz*a1*tz1 &
       !       + rnx*a2*tx2 + rny*a2*ty2 + rnz*a2*tz2 &
       !       + rnx*a3*tx3 + rny*a3*ty3 + rnz*a3*tz3 &
       !       + rnx*a4*tx4 + rny*a4*ty4 + rnz*a4*tz4
       !!
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          flx_ec = flx_ec + a*rnx*tfluxes(n0+n,1,nec) &
                          + a*rny*tfluxes(n0+n,2,nec) &
                          + a*rnz*tfluxes(n0+n,3,nec)
          !
          flx_mx = flx_mx + a*rnx*tfluxes(n0+n,1,nmx) &
                          + a*rny*tfluxes(n0+n,2,nmx) &
                          + a*rnz*tfluxes(n0+n,3,nmx)
          !
          flx_my = flx_my + a*rnx*tfluxes(n0+n,1,nmy) &
                          + a*rny*tfluxes(n0+n,2,nmy) &
                          + a*rnz*tfluxes(n0+n,3,nmy)
          !
          flx_mz = flx_mz + a*rnx*tfluxes(n0+n,1,nmz) &
                          + a*rny*tfluxes(n0+n,2,nmz) &
                          + a*rnz*tfluxes(n0+n,3,nmz)
          !
          flx_ee = flx_ee + a*rnx*tfluxes(n0+n,1,nee) &
                          + a*rny*tfluxes(n0+n,2,nee) &
                          + a*rnz*tfluxes(n0+n,3,nee)
          !
        end do
        !
        flx_ec = ds*flx_ec
        flx_mx = ds*flx_mx
        flx_my = ds*flx_my
        flx_mz = ds*flx_mz
        flx_ee = ds*flx_ee
        !
        ! Add the correction from the interface flux to
        ! the RHS for the points in the current cell
        !
        do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
               correct(this_geom,this_order)%row_ptr(lf+1)
          !
          n = correct(this_geom,this_order)%row_idx(i)
          a = correct(this_geom,this_order)%row_val(i)
          !
!#ifdef DEBUG_ON
!  if (first_step) then
!    write (100,1) nc,ncf,kf,f0+kf,lf,n,n0+n,a
!  end if
!  1 format (1x,i2,4x,i1,2x,i2,2x,i4,2x,i3,2x,i3,2x,i4,3x,es14.6)
!#endif
          residual(n0+n,nec) = residual(n0+n,nec) + a*flx_ec
          residual(n0+n,nmx) = residual(n0+n,nmx) + a*flx_mx
          residual(n0+n,nmy) = residual(n0+n,nmy) + a*flx_my
          residual(n0+n,nmz) = residual(n0+n,nmz) + a*flx_mz
          residual(n0+n,nee) = residual(n0+n,nee) + a*flx_ee
          !
        end do
        !
      end do pts_on_face
      !
    end do faces_of_cell
    !
  end do cell_loop
!#ifdef DEBUG_ON
!  first_step = fals
!#endif
  !
end subroutine add_flux_correction_3d
!
!###############################################################################
!
subroutine add_flux_correction_2d(tfluxes)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,         only : maxFP,geom_solpts
  use geovar,            only : ncell,cell,fp,face
  use interpolation_mod, only : interp
  use metrics_mod,       only : facenrm
  use flowvar,           only : residual,faceflx
  use correction_mod,    only : correct
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: tfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,nc,n0,np,kf,lf,f0
  integer  :: ncf,nfp,num_cell_faces
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,face_rot
  real(wp) :: flx_ec,flx_mx,flx_my,flx_ee
  real(wp) :: a,ds,rnx,rny
  !
  !.. Local Arrays ..
  integer :: fpidx(1:maxFP)
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  !.. Local Arrays ..
  real(wp) :: flux_correction(1:size(residual,dim=1),1:size(residual,dim=2))
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_flux_correction_2d"
  !
!#ifdef DEBUG_ON
!  logical(lk), save :: first_step = true
!#endif
continue
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  flux_correction(:,:) = zero
#endif
  !
!#ifdef DEBUG_ON
!  if (first_step) then
!    write (100,2)
!  2 format ("nc  = cell number",/, &
!            "ncf = face of cell",/, &
!            "kf  = fp local to face of cell",/, &
!            "fp  = global fp",/, &
!            "lf  = fp of all faces of cell",/, &
!            "np  = sp local to cell",/, &
!            "sp  = global sp",/, &
!            "mv  = matrix value",/, &
!            " nc  ncf  kf  fp  lf  np   sp         mv")
!  end if
!#endif
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(correct) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(residual) &
  !$ACC&                   PRESENT(faceflx,facenrm,tfluxes) &
  !$ACC&                   PRESENT(cl) &
  !$ACC&                   PRIVATE(fpidx)
  cell_loop: do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
!#ifdef DEBUG_ON
!  if (first_step) then
!    if (nc == 1) then
!      do n = 1,size(correct(this_geom,this_order)%mat,dim=1)
!        write (101,10) (correct(this_geom,this_order)%mat(n,i), &
!                        i=1,size(correct(this_geom,this_order)%mat,dim=2))
!      end do
!    end if
!  end if
!  10 format (20es11.3)
!#endif
    !
    ! Interpolate the fluxes to the interface flux points and dot this
    ! interpolated vector with the face normal at each flux point
    ! NOTE: We dont need to exchange this information for communication and
    !       periodic boundary faces since this flux is only used within the
    !       current cell.
    !
    ! Loop over the faces of this cell
    num_cell_faces = size(cell(nc)%face)
    !
    faces_of_cell: do ncf = 1,num_cell_faces
      !
      ds = cl(this_geom) * merge( -one , one , cell(nc)%face(ncf)%side==1 )
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      nfp = geom_solpts(face_geom,face_order)
      f0 = face( cell(nc)%face(ncf)%idx )%beg_fp-1
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
      !$ACC LOOP VECTOR
      pts_on_face: do kf = 1,nfp
        !
        rnx = facenrm(1,f0+kf)
        rny = facenrm(2,f0+kf)
        !
        ! Initialize the interface flux values to the common interface flux
        ! NOTE: We initialize it to negative so that we can do a positive
        !       summation in the following loop. The sign of ds will make
        !       sure the flux values are the correct sign when using them
        !       to apply the correction to the interior solution points.
        !
        flx_ec = -faceflx(f0+kf,nec)
        flx_mx = -faceflx(f0+kf,nmx)
        flx_my = -faceflx(f0+kf,nmy)
        flx_ee = -faceflx(f0+kf,nee)
        !
        lf = fpidx(kf)
        !
        ! Get the normal interior discontinuous flux at this flux point
        !
       !flx_ec = rnx*a1*tx1 + rnx*a2*tx2 + rnx*a3*tx3 + rnx*a4*tx4 &
       !       + rny*a1*ty1 + rny*a2*ty2 + rny*a3*ty3 + rny*a4*ty4
       !!
       !flx_ec = rnx*a1*tx1 + rny*a1*ty1 &
       !       + rnx*a2*tx2 + rny*a2*ty2 &
       !       + rnx*a3*tx3 + rny*a3*ty3 &
       !       + rnx*a4*tx4 + rny*a4*ty4
       !!
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          flx_ec = flx_ec + a*rnx*tfluxes(n0+n,1,nec) &
                          + a*rny*tfluxes(n0+n,2,nec)
          !
          flx_mx = flx_mx + a*rnx*tfluxes(n0+n,1,nmx) &
                          + a*rny*tfluxes(n0+n,2,nmx)
          !
          flx_my = flx_my + a*rnx*tfluxes(n0+n,1,nmy) &
                          + a*rny*tfluxes(n0+n,2,nmy)
          !
          flx_ee = flx_ee + a*rnx*tfluxes(n0+n,1,nee) &
                          + a*rny*tfluxes(n0+n,2,nee)
          !
        end do
        !
        flx_ec = ds*flx_ec
        flx_mx = ds*flx_mx
        flx_my = ds*flx_my
        flx_ee = ds*flx_ee
        !
        ! Add the correction from the interface flux to
        ! the RHS for the points in the current cell
        !
        do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
               correct(this_geom,this_order)%row_ptr(lf+1)
          !
          n = correct(this_geom,this_order)%row_idx(i)
          a = correct(this_geom,this_order)%row_val(i)
          !
!#ifdef DEBUG_ON
!  if (first_step) then
!    write (100,1) nc,ncf,kf,f0+kf,lf,n,n0+n,a
!  end if
!  1 format (1x,i2,4x,i1,3x,i1,2x,i2,2x,i2,2x,i2,2x,i3,3x,es14.6)
!#endif
          residual(n0+n,nec) = residual(n0+n,nec) + a*flx_ec
          residual(n0+n,nmx) = residual(n0+n,nmx) + a*flx_mx
          residual(n0+n,nmy) = residual(n0+n,nmy) + a*flx_my
          residual(n0+n,nee) = residual(n0+n,nee) + a*flx_ee
          !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
          flux_correction(n0+n,nec) = flux_correction(n0+n,nec) + a*flx_ec
          flux_correction(n0+n,nmx) = flux_correction(n0+n,nmx) + a*flx_mx
          flux_correction(n0+n,nmy) = flux_correction(n0+n,nmy) + a*flx_my
          flux_correction(n0+n,nee) = flux_correction(n0+n,nee) + a*flx_ee
#endif
          !
        end do
        !
      end do pts_on_face
      !
    end do faces_of_cell
    !
  end do cell_loop
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,n,m
    character(len=30) :: fname
    !
    write (fname,1) "residual_corr",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,cell(ncell)%end_sp
      write (io,2) n,(residual(n,m), m=nec,nee)
    end do
    close (io)
    !
    write (fname,1) "face_flux",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,size(faceflx,dim=1)
      write (io,2) n,(faceflx(n,m), m=nec,nee)
    end do
    close (io)
    !
    write (fname,1) "flux_correction",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,cell(ncell)%end_sp
      write (io,2) n,(flux_correction(n,m), m=nec,nee)
    end do
    close (io)
    !
    if (itcur == 1) then
      open (newunit=io,file="facenrm.dat",status="replace",action="write")
      do n = 1,size(facenrm,dim=2)
        write (io,2) n,(facenrm(m,n), m=1,2)
      end do
      close (io)
    end if
    !
    1 format (a,".",i0,".dat")
    2 format (i7,2x,10es18.10)
    !
  end block
#endif
!#ifdef DEBUG_ON
!  first_step = fals
!#endif
  !
end subroutine add_flux_correction_2d
!
!###############################################################################
!
subroutine add_flux_correction_1d(tfluxes)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,         only : maxFP,geom_solpts
  use geovar,            only : ncell,cell,fp,face
  use interpolation_mod, only : interp
  use metrics_mod,       only : facenrm
  use flowvar,           only : residual,faceflx
  use correction_mod,    only : correct
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: tfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,nc,n0,np,kf,lf,f0
  integer  :: ncf,nfp,num_cell_faces
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,face_rot
  real(wp) :: flx_ec,flx_mx,flx_ee
  real(wp) :: a,ds,rnx
  !
  !.. Local Arrays ..
  integer :: fpidx(1:maxFP)
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  !.. Local Arrays ..
  real(wp) :: flux_correction(1:size(residual,dim=1),1:size(residual,dim=2))
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_flux_correction_1d"
  !
continue
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  flux_correction(:,:) = zero
#endif
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(correct) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(residual) &
  !$ACC&                   PRESENT(faceflx,facenrm,tfluxes) &
  !$ACC&                   PRESENT(cl) &
  !$ACC&                   PRIVATE(fpidx)
  cell_loop: do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Interpolate the fluxes to the interface flux points and dot this
    ! interpolated vector with the face normal at each flux point
    ! NOTE: We dont need to exchange this information for communication and
    !       periodic boundary faces since this flux is only used within the
    !       current cell.
    !
    ! Loop over the faces of this cell
    num_cell_faces = size(cell(nc)%face)
    !
    faces_of_cell: do ncf = 1,num_cell_faces
      !
      ds = cl(this_geom) * merge( -one , one , cell(nc)%face(ncf)%side==1 )
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      nfp = geom_solpts(face_geom,face_order)
      f0 = face( cell(nc)%face(ncf)%idx )%beg_fp-1
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
      !$ACC LOOP VECTOR
      pts_on_face: do kf = 1,nfp
        !
        rnx = facenrm(1,f0+kf)
        !
        ! Initialize the interface flux values to the common interface flux
        ! NOTE: We initialize it to negative so that we can do a positive
        !       summation in the following loop. The sign of ds will make
        !       sure the flux values are the correct sign when using them
        !       to apply the correction to the interior solution points.
        !
        flx_ec = -faceflx(f0+kf,nec)
        flx_mx = -faceflx(f0+kf,nmx)
        flx_ee = -faceflx(f0+kf,nee)
        !
        lf = fpidx(kf)
        !
        ! Get the normal interior discontinuous flux at this flux point
        !
       !flx_ec = rnx*a1*tx1 + rnx*a2*tx2 + rnx*a3*tx3 + rnx*a4*tx4 &
       !       + rny*a1*ty1 + rny*a2*ty2 + rny*a3*ty3 + rny*a4*ty4
       !!
       !flx_ec = rnx*a1*tx1 + rny*a1*ty1 &
       !       + rnx*a2*tx2 + rny*a2*ty2 &
       !       + rnx*a3*tx3 + rny*a3*ty3 &
       !       + rnx*a4*tx4 + rny*a4*ty4
       !!
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          flx_ec = flx_ec + a*rnx*tfluxes(n0+n,1,nec)
          flx_mx = flx_mx + a*rnx*tfluxes(n0+n,1,nmx)
          flx_ee = flx_ee + a*rnx*tfluxes(n0+n,1,nee)
          !
        end do
        !
        flx_ec = ds*flx_ec
        flx_mx = ds*flx_mx
        flx_ee = ds*flx_ee
        !
        ! Add the correction from the interface flux to
        ! the RHS for the points in the current cell
        !
        do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
               correct(this_geom,this_order)%row_ptr(lf+1)
          !
          n = correct(this_geom,this_order)%row_idx(i)
          a = correct(this_geom,this_order)%row_val(i)
          !
          residual(n0+n,nec) = residual(n0+n,nec) + a*flx_ec
          residual(n0+n,nmx) = residual(n0+n,nmx) + a*flx_mx
          residual(n0+n,nee) = residual(n0+n,nee) + a*flx_ee
          !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
          flux_correction(n0+n,nec) = flux_correction(n0+n,nec) + a*flx_ec
          flux_correction(n0+n,nmx) = flux_correction(n0+n,nmx) + a*flx_mx
          flux_correction(n0+n,nee) = flux_correction(n0+n,nee) + a*flx_ee
#endif
          !
        end do
        !
      end do pts_on_face
      !
    end do faces_of_cell
    !
  end do cell_loop
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,n,m
    character(len=30) :: fname
    !
    write (fname,1) "residual_corr",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,cell(ncell)%end_sp
      write (io,2) n,(residual(n,m), m=nec,nee)
    end do
    close (io)
    !
    write (fname,1) "face_flux",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,size(faceflx,dim=1)
      write (io,2) n,(faceflx(n,m), m=nec,nee)
    end do
    close (io)
    !
    write (fname,1) "flux_correction",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,cell(ncell)%end_sp
      write (io,2) n,(flux_correction(n,m), m=nec,nee)
    end do
    close (io)
    !
    if (itcur == 1) then
      open (newunit=io,file="facenrm.dat",status="replace",action="write")
      do n = 1,size(facenrm,dim=2)
        write (io,2) n,(facenrm(m,n), m=1,1)
      end do
      close (io)
    end if
    !
    1 format (a,".",i0,".dat")
    2 format (i7,2x,10es18.10)
    !
  end block
#endif
  !
end subroutine add_flux_correction_1d
!
!###############################################################################
!
subroutine add_flux_correction_nd(mbeg,mend,tfluxes)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,         only : maxFP,geom_solpts
  use geovar,            only : nr,ncell,cell,fp,face
  use interpolation_mod, only : interp
  use metrics_mod,       only : facenrm
  use flowvar,           only : residual,faceflx
  use correction_mod,    only : correct
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: mbeg,mend
  real(wp), intent(in) :: tfluxes(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,l,m,n,nc,n0,np
  integer  :: nf,sf,kf,lf,f0
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,face_rot,nfp,ncf,num_cell_faces
  real(wp) :: a,val,fval,correct_n,ds
  !
  !.. Local Arrays ..
  integer :: fpidx(1:maxFP)
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  !.. Local Arrays ..
  real(wp) :: flux_correction(1:size(residual,dim=1),1:size(residual,dim=2))
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_flux_correction_nd"
  !
continue
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  flux_correction(:,:) = zero
#endif
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(correct) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(residual) &
  !$ACC&                   PRESENT(faceflx,facenrm,tfluxes) &
  !$ACC&                   PRESENT(cl) &
  !$ACC&                   PRIVATE(fpidx)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Interpolate the fluxes to the interface flux points and dot this
    ! interpolated vector with the face normal at each flux point
    ! NOTE: We dont need to exchange this information for communication and
    !       periodic boundary faces since this flux is only used within the
    !       current cell.
    !
    ! Loop over the faces of this cell
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      ds = cl(this_geom) * merge(one,-one,sf==1)
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      nfp = geom_solpts(face_geom,face_order)
      f0 = face(nf)%beg_fp-1
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
      !$ACC LOOP VECTOR COLLAPSE(2)
      do kf = 1,nfp
        do m = mbeg,mend
          fval = zero
          lf = fpidx(kf)
          ! Get the normal interior discontinuous flux at this flux point
          do l = 1,nr
            val = zero
            do i = interp(this_geom,this_order)% &
                             toFace(face_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)% &
                             toFace(face_order)%row_ptr(lf+1)
              !
              n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
              a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
              !
              val = val + a * tfluxes(n0+n,l,m)
              !
            end do
            fval = fval + val * facenrm(l,f0+kf)
          end do
          ! Subtract the interior discontinuous flux from the common flux
          fval = faceflx(f0+kf,m) - fval
          do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
                 correct(this_geom,this_order)%row_ptr(lf+1)
            n = n0 + correct(this_geom,this_order)%row_idx(i)
            correct_n = correct(this_geom,this_order)%row_val(i) * fval * ds
            residual(n,m) = residual(n,m) + correct_n
#ifdef DEBUG_FLUXES_AND_RESIDUAL
            flux_correction(n,m) = flux_correction(n,m) + correct_n
#endif
          end do
        end do
      end do
      !
    end do
    !
  end do
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,n,m
    character(len=30) :: fname
    !
    write (fname,1) "residual_corr",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,cell(ncell)%end_sp
      write (io,2) n,(residual(n,m), m=mbeg,mend)
    end do
    close (io)
    !
    write (fname,1) "face_flux",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,size(faceflx,dim=1)
      write (io,2) n,(faceflx(n,m), m=mbeg,mend)
    end do
    close (io)
    !
    write (fname,1) "flux_correction",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,cell(ncell)%end_sp
      write (io,2) n,(flux_correction(n,m), m=mbeg,mend)
    end do
    close (io)
    !
    if (itcur == 1) then
      open (newunit=io,file="facenrm.dat",status="replace",action="write")
      do n = 1,size(facenrm,dim=2)
        write (io,2) n,(facenrm(m,n), m=1,nr)
      end do
      close (io)
    end if
    !
    1 format (a,".",i0,".dat")
    2 format (i7,2x,10es18.10)
    !
  end block
#endif
  !
end subroutine add_flux_correction_nd
!
!###############################################################################
!
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Procedures for computing common fluxes at face points
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine get_invs_common_fluxes_interior()
  !
  !.. Use Statements ..
  use geovar,  only : using_quadrature
  use flowvar, only : faceflx
  !
  !.. Local Scalars ..
  integer :: m,n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_invs_common_fluxes_interior"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(faceflx)
  do m = 1,size(faceflx,dim=2)
    do n = 1,size(faceflx,dim=1)
      faceflx(n,m) = zero
    end do
  end do
  !
  ! Get the inviscid upwind flux at all of the interior flux points
  !
  if (using_quadrature) then
    call invs_common_fluxes_interior_at_qdtpts
  else
    call invs_common_fluxes_interior_at_flxpts
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_invs_common_fluxes_interior
!
!###############################################################################
!
subroutine invs_common_fluxes_interior_at_flxpts()
  !
  !.. Use Statements ..
  use geovar,      only : n_bnd_flxpts,n_flxpts
  use ovar,        only : invs_flux_method
  use metrics_mod, only : facenrm
  use flowvar,     only : faceusp_array
  use flowvar,     only : facegpsp
  use flowvar,     only : faceflx
  !
  use upwind_flux_mod, only : upwind_flux_gpu
  !
  !.. Local Scalars ..
  integer :: n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "invs_common_fluxes_interior_at_flxpts"
  !
continue
  !
  ! Get the inviscid upwind flux at all of the interior flux points
  !
  n1 = n_bnd_flxpts+1
  n2 = n_flxpts
  !
  call upwind_flux_gpu(n1,n2,facenrm,faceusp_array,facegpsp,faceflx)
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,n,m
    character(len=30) :: fname
    !
    write (fname,1) "faceflx",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,size(faceflx,dim=1)
      write (io,2) n,(faceflx(n,m), m=nec,nee)
    end do
    close(io)
    !
    write (fname,1) "faceusp_L",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,size(faceusp_array,dim=1)
      write (io,2) n,(faceusp_array(n,m,1), m=nec,nee)
    end do
    close(io)
    !
    write (fname,1) "faceusp_R",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,size(faceusp_array,dim=1)
      write (io,2) n,(faceusp_array(n,m,2), m=nec,nee)
    end do
    close(io)
    !
    1 format (a,"_in.",i0,".dat")
    2 format (i7,2x,10es18.10)
    !
  end block
#endif
  !
end subroutine invs_common_fluxes_interior_at_flxpts
!
!###############################################################################
!
subroutine invs_common_fluxes_interior_at_qdtpts()
  !
  !.. Use Statements ..
  use geovar,      only : nfbnd,nface,n_bnd_qdtpts,n_flx_qdtpts
  use ovar,        only : invs_flux_method
  use ovar,        only : combine_common_flux_projections
  use metrics_mod, only : facenrm_qdtpts
  use flowvar,     only : faceusp_array
  use flowvar,     only : faceflx
  !
  use upwind_flux_mod, only : upwind_flux_gpu
  !
  !.. Local Scalars ..
  integer :: n1,n2,nf1,nf2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "invs_common_fluxes_interior_at_qdtpts"
  !
continue
  !
  allocate ( faceuqp(1:n_flx_qdtpts,1:nq,1:2) )
  allocate ( facegpqp(1:n_flx_qdtpts,1:ng,1:2) )
  allocate ( flux_qdtpts(1:n_flx_qdtpts,1:nq) )
  !$ACC ENTER DATA CREATE(faceusp,facegpqp,flux_qdtpts)
  !
  ! Get the inviscid upwind flux at all of the interior flux points
  !
  nf1 = nfbnd+1
  nf2 = nface
  !
  call interp_faceusp_to_faceuqp(nf1,nf2,faceusp_array,faceuqp,facegpqp)
  !
  n1 = n_bnd_qdtpts+1
  n2 = n_flx_qdtpts
  !
  call upwind_flux_gpu(n1,n2,facenrm_qdtpts,faceuqp,facegpqp,flux_qdtpts)
  !
  if (combine_common_flux_projections == Separate_All_Projections) then
    call project_flux_qdtpts_to_faceflx(nf1,nf2,flux_qdtpts,faceflx)
  end if
  !
end subroutine invs_common_fluxes_interior_at_qdtpts
!
!###############################################################################
!
subroutine get_visc_common_fluxes_interior()
  !
  !.. Use Statements ..
  use geovar, only : using_quadrature
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_visc_common_fluxes_interior"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (using_quadrature) then
    call visc_common_fluxes_interior_at_qdtpts
  else
    call visc_common_fluxes_interior_at_flxpts
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_visc_common_fluxes_interior
!
!###############################################################################
!
subroutine visc_common_fluxes_interior_at_flxpts()
  !
  !.. Use Statements ..
  use geovar,      only : n_bnd_flxpts,n_flxpts
  use ovar,        only : visc_variables_to_ave
  use metrics_mod, only : facenrm
  use flowvar,     only : faceusp_array
  use flowvar,     only : facedusp_array
  use flowvar,     only : facegpsp
  use flowvar,     only : faceflx
  !
  use visc_flux_mod, only : common_vflx_ave_fluxes_gpu
  use visc_flux_mod, only : common_vflx_ave_sol_and_grad_gpu
  !
  !.. Local Scalars ..
  integer :: n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "visc_common_fluxes_interior_at_flxpts"
  !
continue
  !
  n1 = n_bnd_flxpts+1
  n2 = n_flxpts
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    call common_vflx_ave_fluxes_gpu(n1,n2,facenrm,faceusp_array,facegpsp, &
                                    facedusp_array,faceflx)
    !
  else ! visc_variables_to_ave == visc_solution_and_gradients
    !
    call common_vflx_ave_sol_and_grad_gpu(n1,n2,facenrm,faceusp_array, &
                                          facedusp_array,faceflx)
    !
  end if
  !
end subroutine visc_common_fluxes_interior_at_flxpts
!
!###############################################################################
!
subroutine visc_common_fluxes_interior_at_qdtpts()
  !
  !.. Use Statements ..
  use geovar,      only : nr,nfbnd,nface,n_bnd_qdtpts,n_flx_qdtpts
  use ovar,        only : visc_variables_to_ave
  use ovar,        only : combine_common_flux_projections
  use metrics_mod, only : facenrm_qdtpts
  use flowvar,     only : facedusp_array
  use flowvar,     only : faceflx
  !
  use visc_flux_mod, only : common_vflx_ave_fluxes_gpu
  use visc_flux_mod, only : common_vflx_ave_sol_and_grad_gpu
  !
  !.. Local Scalars ..
  integer :: n1,n2,nf1,nf2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "visc_common_fluxes_interior_at_qdtpts"
  !
continue
  !
  allocate ( faceduqp(1:n_flx_qdtpts,1:nr,1:nq,1:2) )
  !$ACC ENTER DATA CREATE(facedusp)
  !
  nf1 = nfbnd+1
  nf2 = nface
  !
  call interp_facedusp_to_faceduqp(nf1,nf2,facedusp_array,faceduqp)
  !
  n1 = n_bnd_qdtpts+1
  n2 = n_flx_qdtpts
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    call common_vflx_ave_fluxes_gpu(n1,n2,facenrm_qdtpts,faceuqp,facegpqp, &
                                    faceduqp,flux_qdtpts)
    !
  else ! visc_variables_to_ave == visc_solution_and_gradients
    !
    call common_vflx_ave_sol_and_grad_gpu(n1,n2,facenrm_qdtpts,faceuqp, &
                                          faceduqp,flux_qdtpts)
    !
  end if
  !
  if (combine_common_flux_projections /= Combine_All_Projections) then
    call project_flux_qdtpts_to_faceflx(nf1,nf2,flux_qdtpts,faceflx)
  end if
  !
end subroutine visc_common_fluxes_interior_at_qdtpts
!
!###############################################################################
!
!#define THX_DEBUG
subroutine get_invs_common_fluxes_bc()
  !
  !.. Use Statements ..
#ifdef THX_DEBUG
  use ovar,    only : rhoref,aref,tref,grid_scaling_factor
  use geovar,  only : nfbnd,bface,face,cell_map
  use eqn_idx, only : ntemp,npres
  use flowvar, only : facexyz_array
#endif
  use geovar,        only : using_quadrature,n_bnd_flxpts
  use flowvar,       only : faceusp_array,facegpsp
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Local Scalars ..
  integer :: n1,n2
#ifdef THX_DEBUG
  integer :: n,nf,f1,f2,nc,gc,k,ierr
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: face_for_fp(:)
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_invs_common_fluxes_bc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! First, update the thermodynamic and transport gas properties for the
  ! exterior flux points of all the boundary faces
  !
  n1 = 1
  n2 = n_bnd_flxpts
  !
  call update_thermo_properties(n1,n2,faceusp_array(:,:,2),facegpsp(:,:,2), &
                                      loc="face",caller_name=pname)
#ifdef THX_DEBUG
  !
  face_for_fp = init_array(n_bnd_flxpts,0) ! F2003 AUTO-REALLOC
  !
  do nf = 1,nfbnd
    f1 = face(nf)%beg_fp
    f2 = face(nf)%end_fp
    do n = f1,f2
      face_for_fp(n) = nf
    end do
  end do
  !
  if (any(abs(faceusp_array(n1:n2,nee,2)) > ten*ten)) then
    do n = n1,n2
      if (abs(faceusp_array(n,nee,2)) > ten*ten) then
        nf = face_for_fp(n)
        nc = bface(2,nf)
        gc = cell_map(mypnum+1)%loc_to_glb(nc)
        k = n - face(nf)%beg_fp + 1
        write (10000+mypnum,2) n,k,nf,nc,gc, &
                               facexyz_array(1,n)/grid_scaling_factor, &
                               facexyz_array(2,n)/grid_scaling_factor, &
                               facexyz_array(3,n)/grid_scaling_factor, &
                               faceusp_array(n,nmx,2)*aref, &
                               faceusp_array(n,nmy,2)*aref, &
                               faceusp_array(n,nmz,2)*aref, &
                               facegpsp(n,npres,2)*rhoref*aref*aref, &
                               facegpsp(n,ntemp,2)*tref
      end if
    end do
    write (10000+mypnum,3)
    flush (10000+mypnum)
  end if
  !
  if (allocated(face_for_fp)) deallocate ( face_for_fp , stat=ierr )
  !
  2 format (i8,i5,3i8,100(2x,es14.6,:))
  3 format (//,80("#"),//)
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  call update_transport_properties(n1,n2,faceusp_array(:,:,2),facegpsp(:,:,2))
  !
  if (using_quadrature) then
    call invs_common_fluxes_bc_at_qdtpts
  else
    call invs_common_fluxes_bc_at_flxpts
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_invs_common_fluxes_bc
!
!###############################################################################
!
subroutine invs_common_fluxes_bc_at_flxpts()
  !
  !.. Use Statements ..
  use geovar,      only : n_wall_flxpts,n_bnd_flxpts
  use ovar,        only : invs_flux_method
  use metrics_mod, only : facenrm
  use flowvar,     only : faceusp_array
  use flowvar,     only : facegpsp
  use flowvar,     only : faceflx
  !
  use invs_flux_mod, only : invs_wall_flx_gpu
  !
  use upwind_flux_mod, only : upwind_flux_gpu
  !
  !.. Local Scalars ..
  integer :: n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "invs_common_fluxes_bc_at_flxpts"
  !
continue
  !
  ! Get the inviscid common flux for all flux points on wall BCs
  !
  n1 = 1
  n2 = n_wall_flxpts
  !
  if (n2 >= n1) then
    !
    call invs_wall_flx_gpu(n1,n2,facenrm,faceusp_array,facegpsp,faceflx)
    !
  end if
  !
  ! For all other boundaries, compute the boundary
  ! fluxes in the same way as the interior faces.
  !
  n1 = n_wall_flxpts+1
  n2 = n_bnd_flxpts
  !
  call upwind_flux_gpu(n1,n2,facenrm,faceusp_array,facegpsp,faceflx)
  !
#ifdef DEBUG_FLUXES_AND_RESIDUAL
  block
    !
    use ovar, only : itcur
    integer :: io,n,m
    character(len=30) :: fname
    !
    write (fname,1) "faceflx",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,size(faceflx,dim=1)
      write (io,2) n,(faceflx(n,m), m=nec,nee)
    end do
    close(io)
    !
    write (fname,1) "faceusp_L",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,size(faceusp_array,dim=1)
      write (io,2) n,(faceusp_array(n,m,1), m=nec,nee)
    end do
    close(io)
    !
    write (fname,1) "faceusp_R",itcur
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,size(faceusp_array,dim=1)
      write (io,2) n,(faceusp_array(n,m,2), m=nec,nee)
    end do
    close(io)
    !
    1 format (a,"_bc.",i0,".dat")
    2 format (i7,2x,10es18.10)
    !
  end block
#endif
  !
end subroutine invs_common_fluxes_bc_at_flxpts
!
!###############################################################################
!
subroutine invs_common_fluxes_bc_at_qdtpts()
  !
  !.. Use Statements ..
  use geovar,      only : nfbnd,n_wall_qdtpts,n_bnd_qdtpts
  use ovar,        only : invs_flux_method
  use ovar,        only : combine_common_flux_projections
  use metrics_mod, only : facenrm_qdtpts
  use flowvar,     only : faceusp_array
  use flowvar,     only : faceflx
  !
  use invs_flux_mod, only : invs_wall_flx_gpu
  !
  use upwind_flux_mod, only : upwind_flux_gpu
  !
  !.. Local Scalars ..
  integer :: n1,n2,nf1,nf2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "invs_common_fluxes_bc_at_qdtpts"
  !
continue
  !
  ! Get the inviscid upwind flux at all of the bc flux points
  !
  nf1 = 1
  nf2 = nfbnd
  !
  call interp_faceusp_to_faceuqp(nf1,nf2,faceusp_array,faceuqp,facegpqp)
  !
  n1 = 1
  n2 = n_wall_qdtpts
  !
  if (n2 >= n1) then
    !
    call invs_wall_flx_gpu(n1,n2,facenrm_qdtpts,faceuqp,facegpqp,flux_qdtpts)
    !
  end if
  !
  n1 = n_wall_qdtpts+1
  n2 = n_bnd_qdtpts
  !
  call upwind_flux_gpu(n1,n2,facenrm_qdtpts,faceuqp,facegpqp,flux_qdtpts)
  !
  if (combine_common_flux_projections == Separate_All_Projections) then
    call project_flux_qdtpts_to_faceflx(nf1,nf2,flux_qdtpts,faceflx)
  end if
  !
end subroutine invs_common_fluxes_bc_at_qdtpts
!
!###############################################################################
!
subroutine get_visc_common_fluxes_bc()
  !
  !.. Use Statements ..
  use geovar, only : using_quadrature
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_visc_common_fluxes_bc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (using_quadrature) then
    call visc_common_fluxes_bc_at_qdtpts
  else
    call visc_common_fluxes_bc_at_flxpts
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_visc_common_fluxes_bc
!
!###############################################################################
!
subroutine visc_common_fluxes_bc_at_flxpts()
  !
  !.. Use Statements ..
  use geovar,      only : n_wall_flxpts,n_bnd_flxpts
  use ovar,        only : visc_variables_to_ave
  use metrics_mod, only : facenrm
  use flowvar,     only : faceusp_array
  use flowvar,     only : facedusp_array
  use flowvar,     only : facegpsp
  use flowvar,     only : faceflx
  !
  use visc_flux_mod, only : visc_wall_flx_gpu
  use visc_flux_mod, only : common_vflx_ave_fluxes_gpu
  use visc_flux_mod, only : common_vflx_ave_sol_and_grad_gpu
  !
  !.. Local Scalars ..
  integer :: n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "visc_common_fluxes_bc_at_flxpts"
  !
continue
  !
  ! Get the viscous common flux for all flux points on wall BCs
  !
  n1 = 1
  n2 = n_wall_flxpts
  !
  if (n2 >= n1) then
    !
    call visc_wall_flx_gpu(n1,n2,facenrm,faceusp_array,facegpsp, &
                           facedusp_array,faceflx)
    !
  end if
  !
  ! For all other boundaries, compute the boundary
  ! fluxes in the same way as the interior faces.
  !
  n1 = n_wall_flxpts+1
  n2 = n_bnd_flxpts
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    call common_vflx_ave_fluxes_gpu(n1,n2,facenrm,faceusp_array,facegpsp, &
                                    facedusp_array,faceflx)
    !
  else ! visc_variables_to_ave == visc_solution_and_gradients
    !
    call common_vflx_ave_sol_and_grad_gpu(n1,n2,facenrm,faceusp_array, &
                                          facedusp_array,faceflx)
    !
  end if
  !
end subroutine visc_common_fluxes_bc_at_flxpts
!
!###############################################################################
!
subroutine visc_common_fluxes_bc_at_qdtpts()
  !
  !.. Use Statements ..
  use geovar,      only : nfbnd,nface,n_wall_qdtpts,n_bnd_qdtpts
  use ovar,        only : visc_variables_to_ave
  use ovar,        only : combine_common_flux_projections
  use metrics_mod, only : facenrm_qdtpts
  use flowvar,     only : facedusp_array
  use flowvar,     only : faceflx
  !
  use visc_flux_mod, only : visc_wall_flx_gpu
  use visc_flux_mod, only : common_vflx_ave_fluxes_gpu
  use visc_flux_mod, only : common_vflx_ave_sol_and_grad_gpu
  !
  !.. Local Scalars ..
  integer :: n1,n2,nf1,nf2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "visc_common_fluxes_bc_at_qdtpts"
  !
continue
  !
  nf1 = 1
  nf2 = nfbnd
  !
  call interp_facedusp_to_faceduqp(nf1,nf2,facedusp_array,faceduqp)
  !
  n1 = 1
  n2 = n_wall_qdtpts
  !
  if (n2 >= n1) then
    !
    call visc_wall_flx_gpu(n1,n2,facenrm_qdtpts,faceuqp,facegpqp, &
                           faceduqp,flux_qdtpts)
    !
  end if
  !
  n1 = n_wall_qdtpts+1
  n2 = n_bnd_qdtpts
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    call common_vflx_ave_fluxes_gpu(n1,n2,facenrm_qdtpts,faceuqp,facegpqp, &
                                    faceduqp,flux_qdtpts)
    !
  else ! visc_variables_to_ave == visc_solution_and_gradients
    !
    call common_vflx_ave_sol_and_grad_gpu(n1,n2,facenrm_qdtpts,faceuqp, &
                                          faceduqp,flux_qdtpts)
    !
  end if
  !
  if (combine_common_flux_projections == Combine_All_Projections) then
    nf1 = 1
    nf2 = nface
  end if
  !
  call project_flux_qdtpts_to_faceflx(nf1,nf2,flux_qdtpts,faceflx)
  !
  !$ACC EXIT DATA DELETE(faceduqp,flux_qdtpts,facegpqp,faceuqp)
  deallocate ( faceduqp )
  deallocate ( flux_qdtpts )
  deallocate ( facegpqp )
  deallocate ( faceuqp )
  !
end subroutine visc_common_fluxes_bc_at_qdtpts
!
!###############################################################################
!
subroutine project_flux_qdtpts_to_faceflx(nf1,nf2,flxqp,flxfp)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: flxqp(:,:)
  real(wp), intent(inout) :: flxfp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_flux_qdtpts_to_faceflx"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call project_flux_qdtpts_to_faceflx_3d(nf1,nf2,flxqp,flxfp)
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call project_flux_qdtpts_to_faceflx_2d(nf1,nf2,flxqp,flxfp)
  else
    call project_flux_qdtpts_to_faceflx_nd(nf1,nf2,flxqp,flxfp)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine project_flux_qdtpts_to_faceflx
!
!###############################################################################
!
subroutine project_flux_qdtpts_to_faceflx_3d(nf1,nf2,flxqp,flxfp)
  !
  !.. Use Statements ..
  use geovar,         only : face
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: flxqp(:,:)
  real(wp), intent(inout) :: flxfp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nf,f0,q0,nfp
  integer  :: face_geom,face_order,qdtr_order
  real(wp) :: a,flx_ec,flx_mx,flx_my,flx_mz,flx_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_flux_qdtpts_to_faceflx_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(flxqp) &
  !$ACC&                   PRESENT(flxfp) &
  !$ACC&                   PRESENT(projct)
  do nf = nf1,nf2
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    qdtr_order = face(nf)%qdtr_order
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    q0 = face(nf)%beg_qp - 1
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,nfp
      !
      flx_ec = zero
      flx_mx = zero
      flx_my = zero
      flx_mz = zero
      flx_ee = zero
      !
      do i = projct(face_geom,face_order)%CellOI(qdtr_order)%row_ptr(k)+1, &
             projct(face_geom,face_order)%CellOI(qdtr_order)%row_ptr(k+1)
        !
        j = projct(face_geom,face_order)%CellOI(qdtr_order)%row_idx(i)
        a = projct(face_geom,face_order)%CellOI(qdtr_order)%row_val(i)
        !
        flx_ec = flx_ec + a * flxqp(q0+j,nec)
        flx_mx = flx_mx + a * flxqp(q0+j,nmx)
        flx_my = flx_my + a * flxqp(q0+j,nmy)
        flx_mz = flx_mz + a * flxqp(q0+j,nmz)
        flx_ee = flx_ee + a * flxqp(q0+j,nee)
        !
      end do
      !
      flxfp(f0+k,nec) = flxfp(f0+k,nec) + flx_ec
      flxfp(f0+k,nmx) = flxfp(f0+k,nmx) + flx_mx
      flxfp(f0+k,nmy) = flxfp(f0+k,nmy) + flx_my
      flxfp(f0+k,nmz) = flxfp(f0+k,nmz) + flx_mz
      flxfp(f0+k,nee) = flxfp(f0+k,nee) + flx_ee
      !
    end do
    !
  end do
  !
end subroutine project_flux_qdtpts_to_faceflx_3d
!
!###############################################################################
!
subroutine project_flux_qdtpts_to_faceflx_2d(nf1,nf2,flxqp,flxfp)
  !
  !.. Use Statements ..
  use geovar,         only : face
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: flxqp(:,:)
  real(wp), intent(inout) :: flxfp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nf,f0,q0,nfp
  integer  :: face_geom,face_order,qdtr_order
  real(wp) :: a,flx_ec,flx_mx,flx_my,flx_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_flux_qdtpts_to_faceflx_2d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(flxqp) &
  !$ACC&                   PRESENT(flxfp) &
  !$ACC&                   PRESENT(projct)
  do nf = nf1,nf2
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    qdtr_order = face(nf)%qdtr_order
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    q0 = face(nf)%beg_qp - 1
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,nfp
      !
      flx_ec = zero
      flx_mx = zero
      flx_my = zero
      flx_ee = zero
      !
      do i = projct(face_geom,face_order)%CellOI(qdtr_order)%row_ptr(k)+1, &
             projct(face_geom,face_order)%CellOI(qdtr_order)%row_ptr(k+1)
        !
        j = projct(face_geom,face_order)%CellOI(qdtr_order)%row_idx(i)
        a = projct(face_geom,face_order)%CellOI(qdtr_order)%row_val(i)
        !
        flx_ec = flx_ec + a * flxqp(q0+j,nec)
        flx_mx = flx_mx + a * flxqp(q0+j,nmx)
        flx_my = flx_my + a * flxqp(q0+j,nmy)
        flx_ee = flx_ee + a * flxqp(q0+j,nee)
        !
      end do
      !
      flxfp(f0+k,nec) = flxfp(f0+k,nec) + flx_ec
      flxfp(f0+k,nmx) = flxfp(f0+k,nmx) + flx_mx
      flxfp(f0+k,nmy) = flxfp(f0+k,nmy) + flx_my
      flxfp(f0+k,nee) = flxfp(f0+k,nee) + flx_ee
      !
    end do
    !
  end do
  !
end subroutine project_flux_qdtpts_to_faceflx_2d
!
!###############################################################################
!
subroutine project_flux_qdtpts_to_faceflx_nd(nf1,nf2,flxqp,flxfp)
  !
  !.. Use Statements ..
  use geovar,         only : face
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: flxqp(:,:)
  real(wp), intent(inout) :: flxfp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,nf,f0,q0,nfp
  integer  :: face_geom,face_order,qdtr_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_flux_qdtpts_to_faceflx_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(flxqp) &
  !$ACC&                   PRESENT(flxfp) &
  !$ACC&                   PRESENT(projct)
  do nf = nf1,nf2
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    qdtr_order = face(nf)%qdtr_order
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    q0 = face(nf)%beg_qp - 1
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do k = 1,nfp
      do m = 1,nq
        !
        val = zero
        !
        do i = projct(face_geom,face_order)%CellOI(qdtr_order)%row_ptr(k)+1, &
               projct(face_geom,face_order)%CellOI(qdtr_order)%row_ptr(k+1)
          !
          j = projct(face_geom,face_order)%CellOI(qdtr_order)%row_idx(i)
          a = projct(face_geom,face_order)%CellOI(qdtr_order)%row_val(i)
          !
          val = val + a * flxqp(q0+j,m)
          !
        end do
        !
        flxfp(f0+k,m) = flxfp(f0+k,m) + val
        !
      end do
    end do
    !
  end do
  !
end subroutine project_flux_qdtpts_to_faceflx_nd
!
!###############################################################################
!
subroutine interp_faceusp_to_faceuqp(nf1,nf2,faceusp,faceuqp,facegpqp)
  !
  !.. Use Statements ..
  use geovar,        only : nr,face
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: faceusp(:,:,:)
  real(wp), intent(inout) :: faceuqp(:,:,:)
  real(wp), intent(inout) :: facegpqp(:,:,:)
  !
  !.. Local Scalars ..
  integer :: n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_faceusp_to_faceuqp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call interp_faceusp_to_faceuqp_3d(nf1,nf2,faceusp,faceuqp)
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call interp_faceusp_to_faceuqp_2d(nf1,nf2,faceusp,faceuqp)
  else
    call interp_faceusp_to_faceuqp_nd(nf1,nf2,faceusp,faceuqp)
  end if
  !
  ! Update the thermodynamic and transport properties at the quadrature points
  !
  n1 = face(nf1)%beg_qp
  n2 = face(nf2)%end_qp
  !
  call update_thermo_properties(n1,n2,faceuqp,facegpqp,caller_name=pname)
  call update_transport_properties(n1,n2,faceuqp,facegpqp)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine interp_faceusp_to_faceuqp
!
!###############################################################################
!
subroutine interp_faceusp_to_faceuqp_3d(nf1,nf2,faceusp,faceuqp)
  !
  !.. Use Statements ..
  use geovar,            only : face
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: faceusp(:,:,:)
  real(wp), intent(inout) :: faceuqp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nf,f0,q0,nqp
  integer  :: face_geom,face_order,qdtr_order
  real(wp) :: a
  real(wp) :: l_ec,l_mx,l_my,l_mz,l_ee
  real(wp) :: r_ec,r_mx,r_my,r_mz,r_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_faceusp_to_faceuqp_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(faceuqp) &
  !$ACC&                   PRESENT(interp)
  do nf = nf1,nf2
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    qdtr_order = face(nf)%qdtr_order
    !
    f0 = face(nf)%beg_fp - 1
    !
    q0 = face(nf)%beg_qp - 1
    nqp = face(nf)%end_qp - q0
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,nqp
      !
      l_ec = zero
      l_mx = zero
      l_my = zero
      l_mz = zero
      l_ee = zero
      !
      r_ec = zero
      r_mx = zero
      r_my = zero
      r_mz = zero
      r_ee = zero
      !
      do i = interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k)+1, &
             interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k+1)
        !
        j = interp(face_geom,face_order)%toCell(qdtr_order)%col_idx(i)
        a = interp(face_geom,face_order)%toCell(qdtr_order)%col_val(i)
        !
        l_ec = l_ec + a * faceusp(f0+j,nec,1)
        l_mx = l_mx + a * faceusp(f0+j,nmx,1)
        l_my = l_my + a * faceusp(f0+j,nmy,1)
        l_mz = l_mz + a * faceusp(f0+j,nmz,1)
        l_ee = l_ee + a * faceusp(f0+j,nee,1)
        !
        r_ec = r_ec + a * faceusp(f0+j,nec,2)
        r_mx = r_mx + a * faceusp(f0+j,nmx,2)
        r_my = r_my + a * faceusp(f0+j,nmy,2)
        r_mz = r_mz + a * faceusp(f0+j,nmz,2)
        r_ee = r_ee + a * faceusp(f0+j,nee,2)
        !
      end do
      !
      faceuqp(q0+k,nec,1) = l_ec
      faceuqp(q0+k,nmx,1) = l_mx
      faceuqp(q0+k,nmy,1) = l_my
      faceuqp(q0+k,nmz,1) = l_mz
      faceuqp(q0+k,nee,1) = l_ee
      !
      faceuqp(q0+k,nec,2) = r_ec
      faceuqp(q0+k,nmx,2) = r_mx
      faceuqp(q0+k,nmy,2) = r_my
      faceuqp(q0+k,nmz,2) = r_mz
      faceuqp(q0+k,nee,2) = r_ee
      !
    end do
    !
  end do
  !
end subroutine interp_faceusp_to_faceuqp_3d
!
!###############################################################################
!
subroutine interp_faceusp_to_faceuqp_2d(nf1,nf2,faceusp,faceuqp)
  !
  !.. Use Statements ..
  use geovar,            only : face
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: faceusp(:,:,:)
  real(wp), intent(inout) :: faceuqp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nf,f0,q0,nqp
  integer  :: face_geom,face_order,qdtr_order
  real(wp) :: a
  real(wp) :: l_ec,l_mx,l_my,l_ee
  real(wp) :: r_ec,r_mx,r_my,r_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_faceusp_to_faceuqp_2d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(faceuqp) &
  !$ACC&                   PRESENT(interp)
  do nf = nf1,nf2
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    qdtr_order = face(nf)%qdtr_order
    !
    f0 = face(nf)%beg_fp - 1
    !
    q0 = face(nf)%beg_qp - 1
    nqp = face(nf)%end_qp - q0
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,nqp
      !
      l_ec = zero
      l_mx = zero
      l_my = zero
      l_ee = zero
      !
      r_ec = zero
      r_mx = zero
      r_my = zero
      r_ee = zero
      !
      do i = interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k)+1, &
             interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k+1)
        !
        j = interp(face_geom,face_order)%toCell(qdtr_order)%col_idx(i)
        a = interp(face_geom,face_order)%toCell(qdtr_order)%col_val(i)
        !
        l_ec = l_ec + a * faceusp(f0+j,nec,1)
        l_mx = l_mx + a * faceusp(f0+j,nmx,1)
        l_my = l_my + a * faceusp(f0+j,nmy,1)
        l_ee = l_ee + a * faceusp(f0+j,nee,1)
        !
        r_ec = r_ec + a * faceusp(f0+j,nec,2)
        r_mx = r_mx + a * faceusp(f0+j,nmx,2)
        r_my = r_my + a * faceusp(f0+j,nmy,2)
        r_ee = r_ee + a * faceusp(f0+j,nee,2)
        !
      end do
      !
      faceuqp(q0+k,nec,1) = l_ec
      faceuqp(q0+k,nmx,1) = l_mx
      faceuqp(q0+k,nmy,1) = l_my
      faceuqp(q0+k,nee,1) = l_ee
      !
      faceuqp(q0+k,nec,2) = r_ec
      faceuqp(q0+k,nmx,2) = r_mx
      faceuqp(q0+k,nmy,2) = r_my
      faceuqp(q0+k,nee,2) = r_ee
      !
    end do
    !
  end do
  !
end subroutine interp_faceusp_to_faceuqp_2d
!
!###############################################################################
!
subroutine interp_faceusp_to_faceuqp_nd(nf1,nf2,faceusp,faceuqp)
  !
  !.. Use Statements ..
  use geovar,            only : face
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: faceusp(:,:,:)
  real(wp), intent(inout) :: faceuqp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,n,nf,f0,q0,nqp,i1,i2
  integer  :: face_geom,face_order,qdtr_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_faceusp_to_faceuqp_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(faceuqp) &
  !$ACC&                   PRESENT(interp)
  do nf = nf1,nf2
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    qdtr_order = face(nf)%qdtr_order
    !
    f0 = face(nf)%beg_fp - 1
    !
    q0 = face(nf)%beg_qp - 1
    nqp = face(nf)%end_qp - q0
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR COLLAPSE(4)
    do k = 1,nqp
      do n = 1,2
        do m = 1,nq
          !
          val = zero
          !
          do i = interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k)+1, &
                 interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k+1)
            !
            j = interp(face_geom,face_order)%toCell(qdtr_order)%col_idx(i)
            a = interp(face_geom,face_order)%toCell(qdtr_order)%col_val(i)
            !
            val = val + a * faceusp(f0+j,m,n)
            !
          end do
          !
          faceuqp(q0+k,m,n) = val
          !
        end do
      end do
    end do
    !
  end do
  !
end subroutine interp_faceusp_to_faceuqp_nd
!
!###############################################################################
!
subroutine interp_facedusp_to_faceduqp(nf1,nf2,facedusp,faceduqp)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: facedusp(:,:,:,:)
  real(wp), intent(inout) :: faceduqp(:,:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_facedusp_to_faceduqp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call interp_facedusp_to_faceduqp_3d(nf1,nf2,facedusp,faceduqp)
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call interp_facedusp_to_faceduqp_2d(nf1,nf2,facedusp,faceduqp)
  else
    call interp_facedusp_to_faceduqp_nd(nf1,nf2,facedusp,faceduqp)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine interp_facedusp_to_faceduqp
!
!###############################################################################
!
subroutine interp_facedusp_to_faceduqp_3d(nf1,nf2,facedusp,faceduqp)
  !
  !.. Use Statements ..
  use geovar,            only : face
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: facedusp(:,:,:,:)
  real(wp), intent(inout) :: faceduqp(:,:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nf,f0,q0,nqp
  integer  :: face_geom,face_order,qdtr_order
  real(wp) :: a
  real(wp) :: decdxl,dmxdxl,dmydxl,dmzdxl,deedxl
  real(wp) :: decdyl,dmxdyl,dmydyl,dmzdyl,deedyl
  real(wp) :: decdzl,dmxdzl,dmydzl,dmzdzl,deedzl
  real(wp) :: decdxr,dmxdxr,dmydxr,dmzdxr,deedxr
  real(wp) :: decdyr,dmxdyr,dmydyr,dmzdyr,deedyr
  real(wp) :: decdzr,dmxdzr,dmydzr,dmzdzr,deedzr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_facedusp_to_faceduqp_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(facedusp) &
  !$ACC&                   PRESENT(faceduqp) &
  !$ACC&                   PRESENT(interp)
  do nf = nf1,nf2
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    qdtr_order = face(nf)%qdtr_order
    !
    f0 = face(nf)%beg_fp - 1
    !
    q0 = face(nf)%beg_qp - 1
    nqp = face(nf)%end_qp - q0
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,nqp
      !
      decdxl = zero ; decdyl = zero ; decdzl = zero
      dmxdxl = zero ; dmxdyl = zero ; dmxdzl = zero
      dmydxl = zero ; dmydyl = zero ; dmydzl = zero
      dmzdxl = zero ; dmzdyl = zero ; dmzdzl = zero
      deedxl = zero ; deedyl = zero ; deedzl = zero
      !
      decdxr = zero ; decdyr = zero ; decdzr = zero
      dmxdxr = zero ; dmxdyr = zero ; dmxdzr = zero
      dmydxr = zero ; dmydyr = zero ; dmydzr = zero
      dmzdxr = zero ; dmzdyr = zero ; dmzdzr = zero
      deedxr = zero ; deedyr = zero ; deedzr = zero
      !
      do i = interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k)+1, &
             interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k+1)
        !
        j = interp(face_geom,face_order)%toCell(qdtr_order)%col_idx(i)
        a = interp(face_geom,face_order)%toCell(qdtr_order)%col_val(i)
        !
        decdxl = decdxl + a * facedusp(f0+j,1,nec,1)
        decdyl = decdyl + a * facedusp(f0+j,2,nec,1)
        decdzl = decdzl + a * facedusp(f0+j,3,nec,1)
        !
        dmxdxl = dmxdxl + a * facedusp(f0+j,1,nmx,1)
        dmxdyl = dmxdyl + a * facedusp(f0+j,2,nmx,1)
        dmxdzl = dmxdzl + a * facedusp(f0+j,3,nmx,1)
        !
        dmydxl = dmydxl + a * facedusp(f0+j,1,nmy,1)
        dmydyl = dmydyl + a * facedusp(f0+j,2,nmy,1)
        dmydzl = dmydzl + a * facedusp(f0+j,3,nmy,1)
        !
        dmzdxl = dmzdxl + a * facedusp(f0+j,1,nmz,1)
        dmzdyl = dmzdyl + a * facedusp(f0+j,2,nmz,1)
        dmzdzl = dmzdzl + a * facedusp(f0+j,3,nmz,1)
        !
        deedxl = deedxl + a * facedusp(f0+j,1,nee,1)
        deedyl = deedyl + a * facedusp(f0+j,2,nee,1)
        deedzl = deedzl + a * facedusp(f0+j,3,nee,1)
        !
        decdxr = decdxr + a * facedusp(f0+j,1,nec,2)
        decdyr = decdyr + a * facedusp(f0+j,2,nec,2)
        decdzr = decdzr + a * facedusp(f0+j,3,nec,2)
        !
        dmxdxr = dmxdxr + a * facedusp(f0+j,1,nmx,2)
        dmxdyr = dmxdyr + a * facedusp(f0+j,2,nmx,2)
        dmxdzr = dmxdzr + a * facedusp(f0+j,3,nmx,2)
        !
        dmydxr = dmydxr + a * facedusp(f0+j,1,nmy,2)
        dmydyr = dmydyr + a * facedusp(f0+j,2,nmy,2)
        dmydzr = dmydzr + a * facedusp(f0+j,3,nmy,2)
        !
        dmzdxr = dmzdxr + a * facedusp(f0+j,1,nmz,2)
        dmzdyr = dmzdyr + a * facedusp(f0+j,2,nmz,2)
        dmzdzr = dmzdzr + a * facedusp(f0+j,3,nmz,2)
        !
        deedxr = deedxr + a * facedusp(f0+j,1,nee,2)
        deedyr = deedyr + a * facedusp(f0+j,2,nee,2)
        deedzr = deedzr + a * facedusp(f0+j,3,nee,2)
        !
      end do
      !
      faceduqp(q0+k,1,nec,1) = decdxl
      faceduqp(q0+k,2,nec,1) = decdyl
      faceduqp(q0+k,3,nec,1) = decdzl
      !
      faceduqp(q0+k,1,nmx,1) = dmxdxl
      faceduqp(q0+k,2,nmx,1) = dmxdyl
      faceduqp(q0+k,3,nmx,1) = dmxdzl
      !
      faceduqp(q0+k,1,nmy,1) = dmydxl
      faceduqp(q0+k,2,nmy,1) = dmydyl
      faceduqp(q0+k,3,nmy,1) = dmydzl
      !
      faceduqp(q0+k,1,nmz,1) = dmzdxl
      faceduqp(q0+k,2,nmz,1) = dmzdyl
      faceduqp(q0+k,3,nmz,1) = dmzdzl
      !
      faceduqp(q0+k,1,nee,1) = deedxl
      faceduqp(q0+k,2,nee,1) = deedyl
      faceduqp(q0+k,3,nee,1) = deedzl
      !
      faceduqp(q0+k,1,nec,2) = decdxr
      faceduqp(q0+k,2,nec,2) = decdyr
      faceduqp(q0+k,3,nec,2) = decdzr
      !
      faceduqp(q0+k,1,nmx,2) = dmxdxr
      faceduqp(q0+k,2,nmx,2) = dmxdyr
      faceduqp(q0+k,3,nmx,2) = dmxdzr
      !
      faceduqp(q0+k,1,nmy,2) = dmydxr
      faceduqp(q0+k,2,nmy,2) = dmydyr
      faceduqp(q0+k,3,nmy,2) = dmydzr
      !
      faceduqp(q0+k,1,nmz,2) = dmzdxr
      faceduqp(q0+k,2,nmz,2) = dmzdyr
      faceduqp(q0+k,3,nmz,2) = dmzdzr
      !
      faceduqp(q0+k,1,nee,2) = deedxr
      faceduqp(q0+k,2,nee,2) = deedyr
      faceduqp(q0+k,3,nee,2) = deedzr
      !
    end do
    !
  end do
  !
end subroutine interp_facedusp_to_faceduqp_3d
!
!###############################################################################
!
subroutine interp_facedusp_to_faceduqp_2d(nf1,nf2,facedusp,faceduqp)
  !
  !.. Use Statements ..
  use geovar,            only : face
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: facedusp(:,:,:,:)
  real(wp), intent(inout) :: faceduqp(:,:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,nf,f0,q0,nqp
  integer  :: face_geom,face_order,qdtr_order
  real(wp) :: a
  real(wp) :: decdxl,dmxdxl,dmydxl,deedxl
  real(wp) :: decdyl,dmxdyl,dmydyl,deedyl
  real(wp) :: decdxr,dmxdxr,dmydxr,deedxr
  real(wp) :: decdyr,dmxdyr,dmydyr,deedyr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_facedusp_to_faceduqp_2d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(facedusp) &
  !$ACC&                   PRESENT(faceduqp) &
  !$ACC&                   PRESENT(interp)
  do nf = nf1,nf2
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    qdtr_order = face(nf)%qdtr_order
    !
    f0 = face(nf)%beg_fp - 1
    !
    q0 = face(nf)%beg_qp - 1
    nqp = face(nf)%end_qp - q0
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR
    do k = 1,nqp
      !
      decdxl = zero ; decdyl = zero
      dmxdxl = zero ; dmxdyl = zero
      dmydxl = zero ; dmydyl = zero
      deedxl = zero ; deedyl = zero
      !
      decdxr = zero ; decdyr = zero
      dmxdxr = zero ; dmxdyr = zero
      dmydxr = zero ; dmydyr = zero
      deedxr = zero ; deedyr = zero
      !
      do i = interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k)+1, &
             interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k+1)
        !
        j = interp(face_geom,face_order)%toCell(qdtr_order)%col_idx(i)
        a = interp(face_geom,face_order)%toCell(qdtr_order)%col_val(i)
        !
        decdxl = decdxl + a * facedusp(f0+j,1,nec,1)
        decdyl = decdyl + a * facedusp(f0+j,2,nec,1)
        !
        dmxdxl = dmxdxl + a * facedusp(f0+j,1,nmx,1)
        dmxdyl = dmxdyl + a * facedusp(f0+j,2,nmx,1)
        !
        dmydxl = dmydxl + a * facedusp(f0+j,1,nmy,1)
        dmydyl = dmydyl + a * facedusp(f0+j,2,nmy,1)
        !
        deedxl = deedxl + a * facedusp(f0+j,1,nee,1)
        deedyl = deedyl + a * facedusp(f0+j,2,nee,1)
        !
        decdxr = decdxr + a * facedusp(f0+j,1,nec,2)
        decdyr = decdyr + a * facedusp(f0+j,2,nec,2)
        !
        dmxdxr = dmxdxr + a * facedusp(f0+j,1,nmx,2)
        dmxdyr = dmxdyr + a * facedusp(f0+j,2,nmx,2)
        !
        dmydxr = dmydxr + a * facedusp(f0+j,1,nmy,2)
        dmydyr = dmydyr + a * facedusp(f0+j,2,nmy,2)
        !
        deedxr = deedxr + a * facedusp(f0+j,1,nee,2)
        deedyr = deedyr + a * facedusp(f0+j,2,nee,2)
        !
      end do
      !
      faceduqp(q0+k,1,nec,1) = decdxl
      faceduqp(q0+k,2,nec,1) = decdyl
      !
      faceduqp(q0+k,1,nmx,1) = dmxdxl
      faceduqp(q0+k,2,nmx,1) = dmxdyl
      !
      faceduqp(q0+k,1,nmy,1) = dmydxl
      faceduqp(q0+k,2,nmy,1) = dmydyl
      !
      faceduqp(q0+k,1,nee,1) = deedxl
      faceduqp(q0+k,2,nee,1) = deedyl
      !
      faceduqp(q0+k,1,nec,2) = decdxr
      faceduqp(q0+k,2,nec,2) = decdyr
      !
      faceduqp(q0+k,1,nmx,2) = dmxdxr
      faceduqp(q0+k,2,nmx,2) = dmxdyr
      !
      faceduqp(q0+k,1,nmy,2) = dmydxr
      faceduqp(q0+k,2,nmy,2) = dmydyr
      !
      faceduqp(q0+k,1,nee,2) = deedxr
      faceduqp(q0+k,2,nee,2) = deedyr
      !
    end do
    !
  end do
  !
end subroutine interp_facedusp_to_faceduqp_2d
!
!###############################################################################
!
subroutine interp_facedusp_to_faceduqp_nd(nf1,nf2,facedusp,faceduqp)
  !
  !.. Use Statements ..
  use geovar,            only : nr,face
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf1
  integer,     intent(in) :: nf2
  real(wp),    intent(in) :: facedusp(:,:,:,:)
  real(wp), intent(inout) :: faceduqp(:,:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,m,n,nf,f0,q0,nqp,i1,i2
  integer  :: face_geom,face_order,qdtr_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_facedusp_to_faceduqp_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(facedusp) &
  !$ACC&                   PRESENT(faceduqp) &
  !$ACC&                   PRESENT(interp)
  do nf = nf1,nf2
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    qdtr_order = face(nf)%qdtr_order
    !
    f0 = face(nf)%beg_fp - 1
    !
    q0 = face(nf)%beg_qp - 1
    nqp = face(nf)%end_qp - q0
    !
    ! Now interpolate the derivatives of the conserved variables
    ! from the solution points to the quadrature points
    !
    !$ACC LOOP VECTOR COLLAPSE(4)
    do k = 1,nqp
      do n = 1,2
        do m = 1,nq
          do l = 1,nr
            !
            val = zero
            !
            i1 = interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k)+1
            i2 = interp(face_geom,face_order)%toCell(qdtr_order)%col_ptr(k+1)
            !
            do i = i1,i2
              !
              j = interp(face_geom,face_order)%toCell(qdtr_order)%col_idx(i)
              a = interp(face_geom,face_order)%toCell(qdtr_order)%col_val(i)
              !
              val = val + a * facedusp(f0+j,l,m,n)
              !
            end do
            !
            faceduqp(q0+k,l,m,n) = val
            !
          end do
        end do
      end do
    end do
    !
  end do
  !
end subroutine interp_facedusp_to_faceduqp_nd
!
!###############################################################################
!
!
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Procedures for getting the conservative variable gradients
!
!###############################################################################
!###############################################################################
!###############################################################################
!
!
subroutine discontinuous_gradient_at_solpts()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "discontinuous_gradient_at_solpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (force_nd_subroutines) then
    !
    call discontinuous_gradient_at_solpts_nd(1,nq)
    !
  else
    !
    if (nr == 3) then
      call discontinuous_gradient_at_solpts_3d
    else if (nr == 2) then
      call discontinuous_gradient_at_solpts_2d
    else
      call discontinuous_gradient_at_solpts_1d
    end if
    !
    if (ncs > 0) then
      call discontinuous_gradient_at_solpts_nd(1,ncs)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine discontinuous_gradient_at_solpts
!
!###############################################################################
!
subroutine discontinuous_gradient_at_solpts_3d()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use order_mod,       only : maxSP
  use geovar,          only : ncell,cell
  use metrics_mod,     only : metrics_array
  use derivatives_mod, only : deriv
  use flowvar,         only : usp,dusp
  !
  !.. Local Scalars ..
  integer :: i,j,n,nc,n0,np
  integer :: this_geom,this_order
  real(wp) :: a
  real(wp) :: decdi,dmxdi,dmydi,dmzdi,deedi,didx,didy,didz
  real(wp) :: decdj,dmxdj,dmydj,dmzdj,deedj,djdx,djdy,djdz
  real(wp) :: decdk,dmxdk,dmydk,dmzdk,deedk,dkdx,dkdy,dkdz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "discontinuous_gradient_at_solpts_3d"
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(deriv) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(metrics_array) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(dusp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1  ! beginning index for solpts in cell nc
    np = cell(nc)%end_sp-n0 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      decdi = zero
      dmxdi = zero
      dmydi = zero
      dmzdi = zero
      deedi = zero
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(1)%row_idx(i)
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        !
        decdi = decdi + a * usp(n0+j,nec)
        dmxdi = dmxdi + a * usp(n0+j,nmx)
        dmydi = dmydi + a * usp(n0+j,nmy)
        dmzdi = dmzdi + a * usp(n0+j,nmz)
        deedi = deedi + a * usp(n0+j,nee)
        !
      end do
      !
      decdj = zero
      dmxdj = zero
      dmydj = zero
      dmzdj = zero
      deedj = zero
      !
      do i = deriv(this_geom,this_order)%d(2)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(2)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(2)%row_idx(i)
        a = deriv(this_geom,this_order)%d(2)%row_val(i)
        !
        decdj = decdj + a * usp(n0+j,nec)
        dmxdj = dmxdj + a * usp(n0+j,nmx)
        dmydj = dmydj + a * usp(n0+j,nmy)
        dmzdj = dmzdj + a * usp(n0+j,nmz)
        deedj = deedj + a * usp(n0+j,nee)
        !
      end do
      !
      decdk = zero
      dmxdk = zero
      dmydk = zero
      dmzdk = zero
      deedk = zero
      !
      do i = deriv(this_geom,this_order)%d(3)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(3)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(3)%row_idx(i)
        a = deriv(this_geom,this_order)%d(3)%row_val(i)
        !
        decdk = decdk + a * usp(n0+j,nec)
        dmxdk = dmxdk + a * usp(n0+j,nmx)
        dmydk = dmydk + a * usp(n0+j,nmy)
        dmzdk = dmzdk + a * usp(n0+j,nmz)
        deedk = deedk + a * usp(n0+j,nee)
        !
      end do
      !
      didx = metrics_array(n0+n,1,1)
      djdx = metrics_array(n0+n,2,1)
      dkdx = metrics_array(n0+n,3,1)
      !
      didy = metrics_array(n0+n,1,2)
      djdy = metrics_array(n0+n,2,2)
      dkdy = metrics_array(n0+n,3,2)
      !
      didz = metrics_array(n0+n,1,3)
      djdz = metrics_array(n0+n,2,3)
      dkdz = metrics_array(n0+n,3,3)
      !
      dusp(n0+n,1,nec) = decdi*didx + decdj*djdx + decdk*dkdx
      dusp(n0+n,2,nec) = decdi*didy + decdj*djdy + decdk*dkdy
      dusp(n0+n,3,nec) = decdi*didz + decdj*djdz + decdk*dkdz
      !
      dusp(n0+n,1,nmx) = dmxdi*didx + dmxdj*djdx + dmxdk*dkdx
      dusp(n0+n,2,nmx) = dmxdi*didy + dmxdj*djdy + dmxdk*dkdy
      dusp(n0+n,3,nmx) = dmxdi*didz + dmxdj*djdz + dmxdk*dkdz
      !
      dusp(n0+n,1,nmy) = dmydi*didx + dmydj*djdx + dmydk*dkdx
      dusp(n0+n,2,nmy) = dmydi*didy + dmydj*djdy + dmydk*dkdy
      dusp(n0+n,3,nmy) = dmydi*didz + dmydj*djdz + dmydk*dkdz
      !
      dusp(n0+n,1,nmz) = dmzdi*didx + dmzdj*djdx + dmzdk*dkdx
      dusp(n0+n,2,nmz) = dmzdi*didy + dmzdj*djdy + dmzdk*dkdy
      dusp(n0+n,3,nmz) = dmzdi*didz + dmzdj*djdz + dmzdk*dkdz
      !
      dusp(n0+n,1,nee) = deedi*didx + deedj*djdx + deedk*dkdx
      dusp(n0+n,2,nee) = deedi*didy + deedj*djdy + deedk*dkdy
      dusp(n0+n,3,nee) = deedi*didz + deedj*djdz + deedk*dkdz
      !
    end do
    !
  end do
  !
end subroutine discontinuous_gradient_at_solpts_3d
!
!###############################################################################
!
subroutine discontinuous_gradient_at_solpts_2d()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use order_mod,       only : maxSP
  use geovar,          only : ncell,cell
  use metrics_mod,     only : metrics_array
  use derivatives_mod, only : deriv
  use flowvar,         only : usp,dusp
  !
  !.. Local Scalars ..
  integer :: i,j,n,nc,n0,np
  integer :: this_geom,this_order
  real(wp) :: a
  real(wp) :: decdi,dmxdi,dmydi,deedi,didx,didy
  real(wp) :: decdj,dmxdj,dmydj,deedj,djdx,djdy
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "discontinuous_gradient_at_solpts_2d"
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(deriv) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(metrics_array) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(dusp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1  ! beginning index for solpts in cell nc
    np = cell(nc)%end_sp-n0 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      decdi = zero
      dmxdi = zero
      dmydi = zero
      deedi = zero
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(1)%row_idx(i)
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        !
        decdi = decdi + a * usp(n0+j,nec)
        dmxdi = dmxdi + a * usp(n0+j,nmx)
        dmydi = dmydi + a * usp(n0+j,nmy)
        deedi = deedi + a * usp(n0+j,nee)
        !
      end do
      !
      decdj = zero
      dmxdj = zero
      dmydj = zero
      deedj = zero
      !
      do i = deriv(this_geom,this_order)%d(2)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(2)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(2)%row_idx(i)
        a = deriv(this_geom,this_order)%d(2)%row_val(i)
        !
        decdj = decdj + a * usp(n0+j,nec)
        dmxdj = dmxdj + a * usp(n0+j,nmx)
        dmydj = dmydj + a * usp(n0+j,nmy)
        deedj = deedj + a * usp(n0+j,nee)
        !
      end do
      !
      didx = metrics_array(n0+n,1,1)
      djdx = metrics_array(n0+n,2,1)
      !
      didy = metrics_array(n0+n,1,2)
      djdy = metrics_array(n0+n,2,2)
      !
      dusp(n0+n,1,nec) = decdi*didx + decdj*djdx
      dusp(n0+n,2,nec) = decdi*didy + decdj*djdy
      !
      dusp(n0+n,1,nmx) = dmxdi*didx + dmxdj*djdx
      dusp(n0+n,2,nmx) = dmxdi*didy + dmxdj*djdy
      !
      dusp(n0+n,1,nmy) = dmydi*didx + dmydj*djdx
      dusp(n0+n,2,nmy) = dmydi*didy + dmydj*djdy
      !
      dusp(n0+n,1,nee) = deedi*didx + deedj*djdx
      dusp(n0+n,2,nee) = deedi*didy + deedj*djdy
      !
    end do
    !
  end do
  !
end subroutine discontinuous_gradient_at_solpts_2d
!
!###############################################################################
!
subroutine discontinuous_gradient_at_solpts_1d()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use order_mod,       only : maxSP
  use geovar,          only : ncell,cell
  use metrics_mod,     only : metrics_array
  use derivatives_mod, only : deriv
  use flowvar,         only : usp,dusp
  !
  !.. Local Scalars ..
  integer :: i,j,n,nc,n0,np
  integer :: this_geom,this_order
  real(wp) :: a
  real(wp) :: decdi,dmxdi,deedi,didx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "discontinuous_gradient_at_solpts_1d"
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(deriv) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(metrics_array) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(dusp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1  ! beginning index for solpts in cell nc
    np = cell(nc)%end_sp-n0 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      decdi = zero
      dmxdi = zero
      deedi = zero
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(1)%row_idx(i)
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        !
        decdi = decdi + a * usp(n0+j,nec)
        dmxdi = dmxdi + a * usp(n0+j,nmx)
        deedi = deedi + a * usp(n0+j,nee)
        !
      end do
      !
      didx = metrics_array(n0+n,1,1)
      !
      dusp(n0+n,1,nec) = decdi*didx
      dusp(n0+n,1,nmx) = dmxdi*didx
      dusp(n0+n,1,nee) = deedi*didx
      !
    end do
    !
  end do
  !
end subroutine discontinuous_gradient_at_solpts_1d
!
!###############################################################################
!
subroutine discontinuous_gradient_at_solpts_nd(mbeg,mend)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use order_mod,       only : maxSP
  use geovar,          only : nr,ncell,cell
 !use metrics_mod,     only : metrics_dt
  use metrics_mod,     only : metrics_array
  use derivatives_mod, only : deriv
  use flowvar,         only : usp,dusp
  !
  !.. Formal Arguments ..
  integer, intent(in) :: mbeg,mend
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,m,n,nc,n0,np
  integer :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
 !real(wp), dimension(1:nr) :: dudr
  real(wp), dimension(1:nr,1:nq,1:maxSP) :: dudr
 !real(wp), dimension(1:maxSP,1:nr,1:2) :: dtusp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "discontinuous_gradient_at_solpts_nd"
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(deriv) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(metrics_array) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(dusp) &
  !$ACC&                   PRIVATE(dudr)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1  ! beginning index for solpts in cell nc
    np = cell(nc)%end_sp-n0 ! ending    index for solpts in cell nc
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do k = 1,np
      do m = mbeg,mend
        do l = 1,nr
          !
          val = zero
          do i = deriv(this_geom,this_order)%d(l)%row_ptr(k)+1, &
                 deriv(this_geom,this_order)%d(l)%row_ptr(k+1)
            !
            j = deriv(this_geom,this_order)%d(l)%row_idx(i)
            a = deriv(this_geom,this_order)%d(l)%row_val(i)
            !
            val = val + a * usp(n0+j,m)
            !
          end do
          dudr(l,m,k) = val
        end do
      end do
    end do
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do k = 1,np
      do m = mbeg,mend
        do l = 1,nr
          dusp(n0+k,l,m) = dot( dudr(1:nr,m,k) , metrics_array(n0+k,1:nr,l) )
        end do
      end do
    end do
    !
  end do
  !
end subroutine discontinuous_gradient_at_solpts_nd
!
!###############################################################################
!
subroutine correct_gradients_at_solpts_gpu()
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_at_solpts_gpu"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call correct_gradients_at_solpts_gpu_3d
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    call correct_gradients_at_solpts_gpu_2d
  else
    call correct_gradients_at_solpts_gpu_nd
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine correct_gradients_at_solpts_gpu
!
!###############################################################################
!
subroutine correct_gradients_at_solpts_gpu_3d()
  !
  !.. Use Statements ..
  use geovar,         only : nr,nfbnd,nface,n_flxpts
  use geovar,         only : face,cell,fp
  use metrics_mod,    only : facenrm,jac
  use correction_mod, only : correct
  use flowvar,        only : dusp,faceusp_array
  !
  !.. Local Scalars ..
  integer  :: i,k,l,n,nf,kf,l0,r0,f0,nfp
  integer  :: kl,lc_offset,lc_rot
  integer  :: kr,rc_offset,rc_rot
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: a,bias,rnx,rny,rnz
  real(wp) :: cor_ec,cor_mx,cor_my,cor_mz,cor_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_at_solpts_gpu_3d"
  !
continue
  !
  ! Precompute the biased interface jumps
  !
  !$ACC PARALLEL LOOP &
  !$ACC&              PRESENT(faceusp_array,face_jump)
  do n = 1,n_flxpts
    face_jump(n,nec) = faceusp_array(n,nec,2) - faceusp_array(n,nec,1)
    face_jump(n,nmx) = faceusp_array(n,nmx,2) - faceusp_array(n,nmx,1)
    face_jump(n,nmy) = faceusp_array(n,nmy,2) - faceusp_array(n,nmy,1)
    face_jump(n,nmz) = faceusp_array(n,nmz,2) - faceusp_array(n,nmz,1)
    face_jump(n,nee) = faceusp_array(n,nee,2) - faceusp_array(n,nee,1)
  end do
  !
  ! #######################
  ! #######################
  ! ### BOUNDARY FACES  ###
  ! #######################
  ! #######################
  !
  ! Loop through the boundary faces and correct the gradient of the
  ! continuous solution polynomial in the cells adjacent to the boundary.
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facenrm) &
  !$ACC&              PRESENT(jac) &
  !$ACC&              PRESENT(face) &
  !$ACC&              PRESENT(cell) &
  !$ACC&              PRESENT(fp) &
  !$ACC&              PRESENT(correct) &
  !$ACC&              PRESENT(dusp) &
  !$ACC&              PRESENT(cl) &
  !$ACC&              PRESENT(face_jump) &
  !$ACC&              PRESENT(cbias)
  boundary_faces: do nf = 1,nfbnd
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    host_cell = face(nf)%left%cell ! left  cell on face
    !
    host_geom = cell(host_cell)%geom
    host_order = cell(host_cell)%order
    !
    l0 = cell(host_cell)%beg_sp-1  ! offset index for sol-pts in host cell
    !
    lc_offset = face(nf)%left%fp_offset
    lc_rot = face(nf)%left%rotation
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    !$ACC LOOP VECTOR
    do kf = 1,nfp
      !
      k = f0 + kf
      !
      rnx = facenrm(1,k)
      rny = facenrm(2,k)
      rnz = facenrm(3,k)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kl = lc_offset + fp%geord(face_geom,face_order)%rot(lc_rot)%pts(kf)
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      bias = cl(host_geom) * (cbias(k) + half)
      !
      do i = correct(host_geom,host_order)%row_ptr(kl)+1, &
             correct(host_geom,host_order)%row_ptr(kl+1)
        !
        n = correct(host_geom,host_order)%row_idx(i)
        a = correct(host_geom,host_order)%row_val(i)
        !
        cor_ec = a * face_jump(k,nec) * bias / jac(l0+n)
        cor_mx = a * face_jump(k,nmx) * bias / jac(l0+n)
        cor_my = a * face_jump(k,nmy) * bias / jac(l0+n)
        cor_mz = a * face_jump(k,nmz) * bias / jac(l0+n)
        cor_ee = a * face_jump(k,nee) * bias / jac(l0+n)
        !
        dusp(l0+n,1,nec) = dusp(l0+n,1,nec) + cor_ec*rnx
        dusp(l0+n,2,nec) = dusp(l0+n,2,nec) + cor_ec*rny
        dusp(l0+n,3,nec) = dusp(l0+n,3,nec) + cor_ec*rnz
        !
        dusp(l0+n,1,nmx) = dusp(l0+n,1,nmx) + cor_mx*rnx
        dusp(l0+n,2,nmx) = dusp(l0+n,2,nmx) + cor_mx*rny
        dusp(l0+n,3,nmx) = dusp(l0+n,3,nmx) + cor_mx*rnz
        !
        dusp(l0+n,1,nmy) = dusp(l0+n,1,nmy) + cor_my*rnx
        dusp(l0+n,2,nmy) = dusp(l0+n,2,nmy) + cor_my*rny
        dusp(l0+n,3,nmy) = dusp(l0+n,3,nmy) + cor_my*rnz
        !
        dusp(l0+n,1,nmz) = dusp(l0+n,1,nmz) + cor_mz*rnx
        dusp(l0+n,2,nmz) = dusp(l0+n,2,nmz) + cor_mz*rny
        dusp(l0+n,3,nmz) = dusp(l0+n,3,nmz) + cor_mz*rnz
        !
        dusp(l0+n,1,nee) = dusp(l0+n,1,nee) + cor_ee*rnx
        dusp(l0+n,2,nee) = dusp(l0+n,2,nee) + cor_ee*rny
        dusp(l0+n,3,nee) = dusp(l0+n,3,nee) + cor_ee*rnz
        !
      end do
      !
    end do
    !
  end do boundary_faces
  !
  ! #######################
  ! #######################
  ! ### INTERIOR FACES  ###
  ! #######################
  ! #######################
  !
  ! Loop over the interior faces and correct the gradient of the
  ! continuous solution polynomial in the cells adjacent to the face
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facenrm) &
  !$ACC&              PRESENT(jac) &
  !$ACC&              PRESENT(face) &
  !$ACC&              PRESENT(cell) &
  !$ACC&              PRESENT(fp) &
  !$ACC&              PRESENT(correct) &
  !$ACC&              PRESENT(dusp) &
  !$ACC&              PRESENT(cl) &
  !$ACC&              PRESENT(face_jump) &
  !$ACC&              PRESENT(cbias)
  interior_faces: do nf = nfbnd+1,nface
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    left_cell = face(nf)%left%cell ! left  cell on face
    right_cell = face(nf)%right%cell ! right cell on face
    !
    left_geom = cell(left_cell)%geom
    left_order = cell(left_cell)%order
    !
    right_geom = cell(right_cell)%geom
    right_order = cell(right_cell)%order
    !
    l0 = cell(left_cell)%beg_sp-1  ! offset index for sol-pts in left cell
    r0 = cell(right_cell)%beg_sp-1 ! offset index for sol-pts in right cell
    !
    lc_offset = face(nf)%left%fp_offset
    lc_rot = face(nf)%left%rotation
    !
    rc_offset = face(nf)%right%fp_offset
    rc_rot = face(nf)%right%rotation
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! ##############################
    ! ###  LEFT CELL CORRECTION  ###
    ! ##############################
    !
    !$ACC LOOP VECTOR
    do kf = 1,nfp
      !
      k = f0 + kf
      !
      rnx = facenrm(1,k)
      rny = facenrm(2,k)
      rnz = facenrm(3,k)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kl = lc_offset + fp%geord(face_geom,face_order)%rot(lc_rot)%pts(kf)
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      bias = cl(left_geom) * (cbias(k) + half)
      !
      do i = correct(left_geom,left_order)%row_ptr(kl)+1, &
             correct(left_geom,left_order)%row_ptr(kl+1)
        !
        n = correct(left_geom,left_order)%row_idx(i)
        a = correct(left_geom,left_order)%row_val(i)
        !
        cor_ec = a * face_jump(k,nec) * bias / jac(l0+n)
        cor_mx = a * face_jump(k,nmx) * bias / jac(l0+n)
        cor_my = a * face_jump(k,nmy) * bias / jac(l0+n)
        cor_mz = a * face_jump(k,nmz) * bias / jac(l0+n)
        cor_ee = a * face_jump(k,nee) * bias / jac(l0+n)
        !
        dusp(l0+n,1,nec) = dusp(l0+n,1,nec) + cor_ec*rnx
        dusp(l0+n,2,nec) = dusp(l0+n,2,nec) + cor_ec*rny
        dusp(l0+n,3,nec) = dusp(l0+n,3,nec) + cor_ec*rnz
        !
        dusp(l0+n,1,nmx) = dusp(l0+n,1,nmx) + cor_mx*rnx
        dusp(l0+n,2,nmx) = dusp(l0+n,2,nmx) + cor_mx*rny
        dusp(l0+n,3,nmx) = dusp(l0+n,3,nmx) + cor_mx*rnz
        !
        dusp(l0+n,1,nmy) = dusp(l0+n,1,nmy) + cor_my*rnx
        dusp(l0+n,2,nmy) = dusp(l0+n,2,nmy) + cor_my*rny
        dusp(l0+n,3,nmy) = dusp(l0+n,3,nmy) + cor_my*rnz
        !
        dusp(l0+n,1,nmz) = dusp(l0+n,1,nmz) + cor_mz*rnx
        dusp(l0+n,2,nmz) = dusp(l0+n,2,nmz) + cor_mz*rny
        dusp(l0+n,3,nmz) = dusp(l0+n,3,nmz) + cor_mz*rnz
        !
        dusp(l0+n,1,nee) = dusp(l0+n,1,nee) + cor_ee*rnx
        dusp(l0+n,2,nee) = dusp(l0+n,2,nee) + cor_ee*rny
        dusp(l0+n,3,nee) = dusp(l0+n,3,nee) + cor_ee*rnz
        !
      end do
      !
    end do
    !
    ! ###############################
    ! ###  RIGHT CELL CORRECTION  ###
    ! ###############################
    !
    !$ACC LOOP VECTOR
    do kf = 1,nfp
      !
      k = f0 + kf
      !
      rnx = facenrm(1,k)
      rny = facenrm(2,k)
      rnz = facenrm(3,k)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kr = rc_offset + fp%geord(face_geom,face_order)%rot(rc_rot)%pts(kf)
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the right cell
      !
      bias = cl(right_geom) * (cbias(k) - half)
      !
      do i = correct(right_geom,right_order)%row_ptr(kr)+1, &
             correct(right_geom,right_order)%row_ptr(kr+1)
        !
        n = correct(right_geom,right_order)%row_idx(i)
        a = correct(right_geom,right_order)%row_val(i)
        !
        cor_ec = a * face_jump(k,nec) * bias / jac(r0+n)
        cor_mx = a * face_jump(k,nmx) * bias / jac(r0+n)
        cor_my = a * face_jump(k,nmy) * bias / jac(r0+n)
        cor_mz = a * face_jump(k,nmz) * bias / jac(r0+n)
        cor_ee = a * face_jump(k,nee) * bias / jac(r0+n)
        !
        dusp(r0+n,1,nec) = dusp(r0+n,1,nec) - cor_ec*rnx
        dusp(r0+n,2,nec) = dusp(r0+n,2,nec) - cor_ec*rny
        dusp(r0+n,3,nec) = dusp(r0+n,3,nec) - cor_ec*rnz
        !
        dusp(r0+n,1,nmx) = dusp(r0+n,1,nmx) - cor_mx*rnx
        dusp(r0+n,2,nmx) = dusp(r0+n,2,nmx) - cor_mx*rny
        dusp(r0+n,3,nmx) = dusp(r0+n,3,nmx) - cor_mx*rnz
        !
        dusp(r0+n,1,nmy) = dusp(r0+n,1,nmy) - cor_my*rnx
        dusp(r0+n,2,nmy) = dusp(r0+n,2,nmy) - cor_my*rny
        dusp(r0+n,3,nmy) = dusp(r0+n,3,nmy) - cor_my*rnz
        !
        dusp(r0+n,1,nmz) = dusp(r0+n,1,nmz) - cor_mz*rnx
        dusp(r0+n,2,nmz) = dusp(r0+n,2,nmz) - cor_mz*rny
        dusp(r0+n,3,nmz) = dusp(r0+n,3,nmz) - cor_mz*rnz
        !
        dusp(r0+n,1,nee) = dusp(r0+n,1,nee) - cor_ee*rnx
        dusp(r0+n,2,nee) = dusp(r0+n,2,nee) - cor_ee*rny
        dusp(r0+n,3,nee) = dusp(r0+n,3,nee) - cor_ee*rnz
        !
      end do
      !
    end do
    !
  end do interior_faces
  !
end subroutine correct_gradients_at_solpts_gpu_3d
!
!###############################################################################
!
subroutine correct_gradients_at_solpts_gpu_2d()
  !
  !.. Use Statements ..
  use geovar,         only : nr,nfbnd,nface,n_flxpts
  use geovar,         only : face,cell,fp
  use metrics_mod,    only : facenrm,jac
  use correction_mod, only : correct
  use flowvar,        only : dusp,faceusp_array
  !
  !.. Local Scalars ..
  integer  :: i,k,l,n,nf,kf,l0,r0,f0,nfp
  integer  :: kl,lc_offset,lc_rot
  integer  :: kr,rc_offset,rc_rot
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: a,bias,rnx,rny
  real(wp) :: cor_ec,cor_mx,cor_my,cor_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_at_solpts_gpu_2d"
  !
continue
  !
  ! Precompute the biased interface jumps
  !
  !$ACC PARALLEL LOOP &
  !$ACC&              PRESENT(faceusp_array,face_jump)
  do n = 1,n_flxpts
    face_jump(n,nec) = faceusp_array(n,nec,2) - faceusp_array(n,nec,1)
    face_jump(n,nmx) = faceusp_array(n,nmx,2) - faceusp_array(n,nmx,1)
    face_jump(n,nmy) = faceusp_array(n,nmy,2) - faceusp_array(n,nmy,1)
    face_jump(n,nee) = faceusp_array(n,nee,2) - faceusp_array(n,nee,1)
  end do
  !
  ! #######################
  ! #######################
  ! ### BOUNDARY FACES  ###
  ! #######################
  ! #######################
  !
  ! Loop through the boundary faces and correct the gradient of the
  ! continuous solution polynomial in the cells adjacent to the boundary.
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facenrm) &
  !$ACC&              PRESENT(jac) &
  !$ACC&              PRESENT(face) &
  !$ACC&              PRESENT(cell) &
  !$ACC&              PRESENT(fp) &
  !$ACC&              PRESENT(correct) &
  !$ACC&              PRESENT(dusp) &
  !$ACC&              PRESENT(cl) &
  !$ACC&              PRESENT(face_jump) &
  !$ACC&              PRESENT(cbias)
  boundary_faces: do nf = 1,nfbnd
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    host_cell = face(nf)%left%cell ! left  cell on face
    !
    host_geom = cell(host_cell)%geom
    host_order = cell(host_cell)%order
    !
    l0 = cell(host_cell)%beg_sp-1  ! offset index for sol-pts in host cell
    !
    lc_offset = face(nf)%left%fp_offset
    lc_rot = face(nf)%left%rotation
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    !$ACC LOOP VECTOR
    do kf = 1,nfp
      !
      k = f0 + kf
      !
      rnx = facenrm(1,k)
      rny = facenrm(2,k)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kl = lc_offset + fp%geord(face_geom,face_order)%rot(lc_rot)%pts(kf)
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      bias = cl(host_geom) * (cbias(k) + half)
      !
      do i = correct(host_geom,host_order)%row_ptr(kl)+1, &
             correct(host_geom,host_order)%row_ptr(kl+1)
        !
        n = correct(host_geom,host_order)%row_idx(i)
        a = correct(host_geom,host_order)%row_val(i)
        !
        cor_ec = a * face_jump(k,nec) * bias / jac(l0+n)
        cor_mx = a * face_jump(k,nmx) * bias / jac(l0+n)
        cor_my = a * face_jump(k,nmy) * bias / jac(l0+n)
        cor_ee = a * face_jump(k,nee) * bias / jac(l0+n)
        !
        dusp(l0+n,1,nec) = dusp(l0+n,1,nec) + cor_ec*rnx
        dusp(l0+n,2,nec) = dusp(l0+n,2,nec) + cor_ec*rny
        !
        dusp(l0+n,1,nmx) = dusp(l0+n,1,nmx) + cor_mx*rnx
        dusp(l0+n,2,nmx) = dusp(l0+n,2,nmx) + cor_mx*rny
        !
        dusp(l0+n,1,nmy) = dusp(l0+n,1,nmy) + cor_my*rnx
        dusp(l0+n,2,nmy) = dusp(l0+n,2,nmy) + cor_my*rny
        !
        dusp(l0+n,1,nee) = dusp(l0+n,1,nee) + cor_ee*rnx
        dusp(l0+n,2,nee) = dusp(l0+n,2,nee) + cor_ee*rny
        !
      end do
      !
    end do
    !
  end do boundary_faces
  !
  ! #######################
  ! #######################
  ! ### INTERIOR FACES  ###
  ! #######################
  ! #######################
  !
  ! Loop over the interior faces and correct the gradient of the
  ! continuous solution polynomial in the cells adjacent to the face
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facenrm) &
  !$ACC&              PRESENT(jac) &
  !$ACC&              PRESENT(face) &
  !$ACC&              PRESENT(cell) &
  !$ACC&              PRESENT(fp) &
  !$ACC&              PRESENT(correct) &
  !$ACC&              PRESENT(dusp) &
  !$ACC&              PRESENT(cl) &
  !$ACC&              PRESENT(face_jump) &
  !$ACC&              PRESENT(cbias)
  interior_faces: do nf = nfbnd+1,nface
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    left_cell = face(nf)%left%cell ! left  cell on face
    right_cell = face(nf)%right%cell ! right cell on face
    !
    left_geom = cell(left_cell)%geom
    left_order = cell(left_cell)%order
    !
    right_geom = cell(right_cell)%geom
    right_order = cell(right_cell)%order
    !
    l0 = cell(left_cell)%beg_sp-1  ! offset index for sol-pts in left cell
    r0 = cell(right_cell)%beg_sp-1 ! offset index for sol-pts in right cell
    !
    lc_offset = face(nf)%left%fp_offset
    lc_rot = face(nf)%left%rotation
    !
    rc_offset = face(nf)%right%fp_offset
    rc_rot = face(nf)%right%rotation
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! ##############################
    ! ###  LEFT CELL CORRECTION  ###
    ! ##############################
    !
    !$ACC LOOP VECTOR
    do kf = 1,nfp
      !
      k = f0 + kf
      !
      rnx = facenrm(1,k)
      rny = facenrm(2,k)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kl = lc_offset + fp%geord(face_geom,face_order)%rot(lc_rot)%pts(kf)
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      bias = cl(left_geom) * (cbias(k) + half)
      !
      do i = correct(left_geom,left_order)%row_ptr(kl)+1, &
             correct(left_geom,left_order)%row_ptr(kl+1)
        !
        n = correct(left_geom,left_order)%row_idx(i)
        a = correct(left_geom,left_order)%row_val(i)
        !
        cor_ec = a * face_jump(k,nec) * bias / jac(l0+n)
        cor_mx = a * face_jump(k,nmx) * bias / jac(l0+n)
        cor_my = a * face_jump(k,nmy) * bias / jac(l0+n)
        cor_ee = a * face_jump(k,nee) * bias / jac(l0+n)
        !
        dusp(l0+n,1,nec) = dusp(l0+n,1,nec) + cor_ec*rnx
        dusp(l0+n,2,nec) = dusp(l0+n,2,nec) + cor_ec*rny
        !
        dusp(l0+n,1,nmx) = dusp(l0+n,1,nmx) + cor_mx*rnx
        dusp(l0+n,2,nmx) = dusp(l0+n,2,nmx) + cor_mx*rny
        !
        dusp(l0+n,1,nmy) = dusp(l0+n,1,nmy) + cor_my*rnx
        dusp(l0+n,2,nmy) = dusp(l0+n,2,nmy) + cor_my*rny
        !
        dusp(l0+n,1,nee) = dusp(l0+n,1,nee) + cor_ee*rnx
        dusp(l0+n,2,nee) = dusp(l0+n,2,nee) + cor_ee*rny
        !
      end do
      !
    end do
    !
    ! ###############################
    ! ###  RIGHT CELL CORRECTION  ###
    ! ###############################
    !
    !$ACC LOOP VECTOR
    do kf = 1,nfp
      !
      k = f0 + kf
      !
      rnx = facenrm(1,k)
      rny = facenrm(2,k)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kr = rc_offset + fp%geord(face_geom,face_order)%rot(rc_rot)%pts(kf)
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the right cell
      !
      bias = cl(right_geom) * (cbias(k) - half)
      !
      do i = correct(right_geom,right_order)%row_ptr(kr)+1, &
             correct(right_geom,right_order)%row_ptr(kr+1)
        !
        n = correct(right_geom,right_order)%row_idx(i)
        a = correct(right_geom,right_order)%row_val(i)
        !
        cor_ec = a * face_jump(k,nec) * bias / jac(r0+n)
        cor_mx = a * face_jump(k,nmx) * bias / jac(r0+n)
        cor_my = a * face_jump(k,nmy) * bias / jac(r0+n)
        cor_ee = a * face_jump(k,nee) * bias / jac(r0+n)
        !
        dusp(r0+n,1,nec) = dusp(r0+n,1,nec) - cor_ec*rnx
        dusp(r0+n,2,nec) = dusp(r0+n,2,nec) - cor_ec*rny
        !
        dusp(r0+n,1,nmx) = dusp(r0+n,1,nmx) - cor_mx*rnx
        dusp(r0+n,2,nmx) = dusp(r0+n,2,nmx) - cor_mx*rny
        !
        dusp(r0+n,1,nmy) = dusp(r0+n,1,nmy) - cor_my*rnx
        dusp(r0+n,2,nmy) = dusp(r0+n,2,nmy) - cor_my*rny
        !
        dusp(r0+n,1,nee) = dusp(r0+n,1,nee) - cor_ee*rnx
        dusp(r0+n,2,nee) = dusp(r0+n,2,nee) - cor_ee*rny
        !
      end do
      !
    end do
    !
  end do interior_faces
  !
end subroutine correct_gradients_at_solpts_gpu_2d
!
!###############################################################################
!
subroutine correct_gradients_at_solpts_gpu_nd()
  !
  !.. Use Statements ..
  use geovar,         only : nr,nfbnd,nface,n_flxpts
  use geovar,         only : face,cell,fp
  use metrics_mod,    only : facenrm,jac
  use correction_mod, only : correct
  use flowvar,        only : dusp,faceusp_array
  !
  !.. Local Scalars ..
  integer  :: i,k,l,m,n,nf,l0,r0,f0,nfp
  integer  :: kl,lc_offset,lc_rot
  integer  :: kr,rc_offset,rc_rot
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: a,bias,cor_m
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_at_solpts_gpu_nd"
  !
continue
  !
  ! Precompute the biased interface jumps
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(faceusp_array,face_jump)
  do m = 1,nq
    do n = 1,n_flxpts
      face_jump(n,m) = faceusp_array(n,m,2) - faceusp_array(n,m,1)
    end do
  end do
  !
  ! #######################
  ! #######################
  ! ### BOUNDARY FACES  ###
  ! #######################
  ! #######################
  !
  ! Loop through the boundary faces and correct the gradient of the
  ! continuous solution polynomial in the cells adjacent to the boundary.
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facenrm) &
  !$ACC&              PRESENT(jac) &
  !$ACC&              PRESENT(face) &
  !$ACC&              PRESENT(cell) &
  !$ACC&              PRESENT(fp) &
  !$ACC&              PRESENT(correct) &
  !$ACC&              PRESENT(dusp) &
  !$ACC&              PRESENT(cl) &
  !$ACC&              PRESENT(face_jump) &
  !$ACC&              PRESENT(cbias)
  boundary_faces: do nf = 1,nfbnd
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    host_cell = face(nf)%left%cell ! left  cell on face
    !
    host_geom = cell(host_cell)%geom
    host_order = cell(host_cell)%order
    !
    l0 = cell(host_cell)%beg_sp-1  ! offset index for sol-pts in host cell
    !
    lc_offset = face(nf)%left%fp_offset
    lc_rot = face(nf)%left%rotation
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do k = 1,nfp
      do m = 1,nq
        do l = 1,nr
          !
          ! Get the index for the correction function
          ! corresponding to the current interface flux point
          !
          kl = lc_offset + fp%geord(face_geom,face_order)%rot(lc_rot)%pts(k)
          !
          ! Now use the correction function with the interface jump to correct
          ! the gradients of the conservative variables in the left cell
          !
          bias = cl(host_geom) * (cbias(f0+k) + half)
          !
          do i = correct(host_geom,host_order)%row_ptr(kl)+1, &
                 correct(host_geom,host_order)%row_ptr(kl+1)
            !
            n = correct(host_geom,host_order)%row_idx(i)
            a = correct(host_geom,host_order)%row_val(i)
            !
            cor_m = a * face_jump(f0+k,m) * bias / jac(l0+n)
            !
            dusp(l0+n,l,m) = dusp(l0+n,l,m) + cor_m * facenrm(l,f0+k)
            !
          end do
          !
        end do
      end do
    end do
    !
  end do boundary_faces
  !
  ! #######################
  ! #######################
  ! ### INTERIOR FACES  ###
  ! #######################
  ! #######################
  !
  ! Loop over the interior faces and correct the gradient of the
  ! continuous solution polynomial in the cells adjacent to the face
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facenrm) &
  !$ACC&              PRESENT(jac) &
  !$ACC&              PRESENT(face) &
  !$ACC&              PRESENT(cell) &
  !$ACC&              PRESENT(fp) &
  !$ACC&              PRESENT(correct) &
  !$ACC&              PRESENT(dusp) &
  !$ACC&              PRESENT(cl) &
  !$ACC&              PRESENT(face_jump) &
  !$ACC&              PRESENT(cbias)
  interior_faces: do nf = nfbnd+1,nface
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    left_cell = face(nf)%left%cell ! left  cell on face
    right_cell = face(nf)%right%cell ! right cell on face
    !
    left_geom = cell(left_cell)%geom
    left_order = cell(left_cell)%order
    !
    right_geom = cell(right_cell)%geom
    right_order = cell(right_cell)%order
    !
    l0 = cell(left_cell)%beg_sp-1  ! offset index for sol-pts in left cell
    r0 = cell(right_cell)%beg_sp-1 ! offset index for sol-pts in right cell
    !
    lc_offset = face(nf)%left%fp_offset
    lc_rot = face(nf)%left%rotation
    !
    rc_offset = face(nf)%right%fp_offset
    rc_rot = face(nf)%right%rotation
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! ##############################
    ! ###  LEFT CELL CORRECTION  ###
    ! ##############################
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do k = 1,nfp
      do m = 1,nq
        do l = 1,nr
          !
          ! Get the index for the correction function
          ! corresponding to the current interface flux point
          !
          kl = lc_offset + fp%geord(face_geom,face_order)%rot(lc_rot)%pts(k)
          !
          ! Now use the correction function with the interface jump to correct
          ! the gradients of the conservative variables in the left cell
          !
          bias = cl(left_geom) * (cbias(f0+k) + half)
          !
          do i = correct(left_geom,left_order)%row_ptr(kl)+1, &
                 correct(left_geom,left_order)%row_ptr(kl+1)
            !
            n = correct(left_geom,left_order)%row_idx(i)
            a = correct(left_geom,left_order)%row_val(i)
            !
            cor_m = a * face_jump(f0+k,m) * bias / jac(l0+n)
            !
            dusp(l0+n,l,m) = dusp(l0+n,l,m) + cor_m * facenrm(l,f0+k)
            !
          end do
          !
        end do
      end do
    end do
    !
    ! ###############################
    ! ###  RIGHT CELL CORRECTION  ###
    ! ###############################
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do k = 1,nfp
      do m = 1,nq
        do l = 1,nr
          !
          ! Get the index for the correction function
          ! corresponding to the current interface flux point
          !
          kr = rc_offset + fp%geord(face_geom,face_order)%rot(rc_rot)%pts(k)
          !
          ! Now use the correction function with the interface jump to correct
          ! the gradients of the conservative variables in the right cell
          !
          bias = cl(right_geom) * (cbias(f0+k) - half)
          !
          do i = correct(right_geom,right_order)%row_ptr(kr)+1, &
                 correct(right_geom,right_order)%row_ptr(kr+1)
            !
            n = correct(right_geom,right_order)%row_idx(i)
            a = correct(right_geom,right_order)%row_val(i)
            !
            cor_m = a * face_jump(f0+k,m) * bias / jac(r0+n)
            !
            dusp(r0+n,l,m) = dusp(r0+n,l,m) - cor_m * facenrm(l,f0+k)
            !
          end do
          !
        end do
      end do
    end do
    !
  end do interior_faces
  !
end subroutine correct_gradients_at_solpts_gpu_nd
!
!###############################################################################
!
subroutine init_cbias()
  !
  !.. Use Statements ..
  use geovar,      only : nr,n_flxpts,nfbnd,bface,face
  use ovar,        only : walls_are_exact,visc_flux_method,LDG_beta
  use metrics_mod, only : facenrm
  !
  !.. Local Scalars ..
  integer :: n,nf,n1,n2
  !
  !.. Local Arrays ..
  real(wp) :: unrm(1:nr)
  !
continue
  !
  allocate ( cbias(1:n_flxpts) )
  !$ACC UPDATE DEVICE(cbias)
  !
  if (visc_flux_method == visc_flux_BR2) then
    do n = 1,n_flxpts
      cbias(n) = half
    end do
  else
    do n = 1,n_flxpts
      unrm(1:nr) = unit_vector( facenrm(1:nr,n) )
      cbias(n) = dot_product( LDG_beta(1:nr) , unrm(1:nr) )
    end do
  end if
  !
  if (walls_are_exact) then
    do nf = 1,nfbnd
      if (any(bface(1,nf) == bc_walls)) then
        n1 = face(nf)%beg_fp
        n2 = face(nf)%end_fp
        do n = n1,n2
          cbias(n) = one
        end do
      end if
    end do
  end if
  !$ACC UPDATE DEVICE(cbias)
  !
  allocate ( face_jump(1:n_flxpts,1:nq) )
  !$ACC UPDATE DEVICE(face_jump)
  !
end subroutine init_cbias
!
!###############################################################################
!
subroutine correct_gradients_br2_gpu()
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_br2_gpu"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (force_nd_subroutines) then
    !
    call correct_gradients_br2_gpu_nd(1,nq)
    !
  else
    !
    if (nr == 3) then
      call correct_gradients_br2_gpu_3d
    else if (nr == 2) then
      call correct_gradients_br2_gpu_2d
    else
      call correct_gradients_br2_gpu_1d
    end if
    !
    if (ncs > 0) then
      call correct_gradients_br2_gpu_nd(1,ncs)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine correct_gradients_br2_gpu
!
!###############################################################################
!
subroutine correct_gradients_br2_gpu_3d()
  !
  use order_mod,         only : maxSP
  use geovar,            only : n_flxpts,ncell,cell,fp,face
  use metrics_mod,       only : jac,facenrm
  use correction_mod,    only : correct
  use flowvar,           only : dusp
  use flowvar,           only : faceusp_array
  use flowvar,           only : facedusp_array
  use interpolation_mod, only : interp
  !
  !.. Local Scalars ..
  integer :: i,n,nc,n1,n2,np,k,n0
  integer :: nf,kf,lf,sf,ncf,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,rnx,rny,rnz
  real(wp) :: decdx,dmxdx,dmydx,dmzdx,deedx
  real(wp) :: decdy,dmxdy,dmydy,dmzdy,deedy
  real(wp) :: decdz,dmxdz,dmydz,dmzdz,deedz
  real(wp) :: cor_ec,cor_mx,cor_my,cor_mz,cor_ee
  !
  !.. Local Arrays ..
  real(wp) :: face_corr(1:maxSP,1:6,1:3,1:nq)
  real(wp) :: cell_invjac(1:maxSP)
  !
  integer :: array(1:1) ! DONT REMOVE THIS
  !
  integer  :: cf_f0(1:6)
  integer  :: cf_sf(1:6)
  integer  :: cf_nfp(1:6)
  integer  :: cf_face_order(1:6)
  integer  :: cf_face_geom(1:6)
  integer  :: cf_fp_offset(1:6)
  integer  :: cf_face_rot(1:6)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_br2_gpu_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
 !!$OMP PARALLEL &
 !!$OMP&         DEFAULT(SHARED) &
 !!$OMP&         PRIVATE(n,a,nc,n1,n2,np,n0,this_geom,this_order, &
 !!$OMP&                 num_cell_faces,cell_invjac,ncf,nf,cf_sf, &
 !!$OMP&                 cf_face_geom,cf_face_order,cf_fp_offset, &
 !!$OMP&                 cf_face_rot,cf_f0,cf_nfp,face_corr, &
 !!$OMP&                 kf,k,rnx,rny,rnz,face_order,face_geom, &
 !!$OMP&                 fp_offset,face_rot,lf,i,cor_ec,cor_mx, &
 !!$OMP&                 cor_my,cor_mz,cor_ee,sf, &
 !!$OMP&                 decdx,dmxdx,dmydx,dmzdx,deedx, &
 !!$OMP&                 decdy,dmxdy,dmydy,dmzdy,deedy, &
 !!$OMP&                 decdz,dmxdz,dmydz,dmzdz,deedz)
  !
#ifdef WRITE_DATA_FOR_PROFILING
  if (num_prof_br2_outputs == 0) then
    if (mypnum == glb_root) then
      call output_br2_data
    end if
  end if
#endif
  !
#ifdef TIME_MINIAPPS_FOR_GPU_COMPARISON
  timings_counter1 = timings_counter1 + 1
  timings_points1(timings_counter1) = n_flxpts
  call gfr_timer( timings_beg1(timings_counter1) )
#endif
  !$ACC PARALLEL LOOP &
  !$ACC&              PRESENT(faceusp_array,cbias,face_jump)
 !!$OMP DO
  do n = 1,n_flxpts
    a = cbias(n)
    face_jump(n,nec) = a*( faceusp_array(n,nec,2) - faceusp_array(n,nec,1) )
    face_jump(n,nmx) = a*( faceusp_array(n,nmx,2) - faceusp_array(n,nmx,1) )
    face_jump(n,nmy) = a*( faceusp_array(n,nmy,2) - faceusp_array(n,nmy,1) )
    face_jump(n,nmz) = a*( faceusp_array(n,nmz,2) - faceusp_array(n,nmz,1) )
    face_jump(n,nee) = a*( faceusp_array(n,nee,2) - faceusp_array(n,nee,1) )
  end do
 !!$OMP END DO
  !
#ifdef TIME_MINIAPPS_FOR_GPU_COMPARISON
  call gfr_timer( timings_end1(timings_counter1) )
#endif
  !
  !
#ifdef TIME_MINIAPPS_FOR_GPU_COMPARISON
  timings_counter = timings_counter + 1
  timings_points(timings_counter) = ncell
  call gfr_timer( timings_beg(timings_counter) )
#endif
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facedusp_array,facenrm) &
  !$ACC&              PRESENT(jac) &
  !$ACC&              PRESENT(face) &
  !$ACC&              PRESENT(cell) &
  !$ACC&              PRESENT(fp) &
  !$ACC&              PRESENT(correct) &
  !$ACC&              PRESENT(interp) &
  !$ACC&              PRESENT(dusp) &
  !$ACC&              PRESENT(cl) &
  !$ACC&              PRESENT(face_jump) &
  !$ACC&              PRIVATE(face_corr,cf_f0,cf_sf,cell_invjac, &
  !$ACC&                      cf_face_geom,cf_fp_offset,cf_face_rot, &
  !$ACC&                      cf_nfp,cf_face_order,array)
 !!$OMP DO
  cell_loop: do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    n0 = n1-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now Interpolate the interior solution gradients to the flux points
    !
    num_cell_faces = size(cell(nc)%face)
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      cell_invjac(n) = cl(this_geom) / jac(n0+n)
    end do
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      !
      cf_sf(ncf) = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      cf_face_geom(ncf) = cell(nc)%face(ncf)%geom   ! geom  of current face
      cf_face_order(ncf) = cell(nc)%face(ncf)%order ! order of current face
      cf_fp_offset(ncf) = cell(nc)%face(ncf)%fp_offset
      cf_face_rot(ncf) = cell(nc)%face(ncf)%rotation
      !
      cf_f0(ncf) = face(nf)%beg_fp - 1
      cf_nfp(ncf) = face(nf)%end_fp - cf_f0(ncf)
      !
    end do
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do ncf = 1,6
      do n = 1,np
        !
        face_corr(n,ncf,1,nec) = zero
        face_corr(n,ncf,2,nec) = zero
        face_corr(n,ncf,3,nec) = zero
        !
        face_corr(n,ncf,1,nmx) = zero
        face_corr(n,ncf,2,nmx) = zero
        face_corr(n,ncf,3,nmx) = zero
        !
        face_corr(n,ncf,1,nmy) = zero
        face_corr(n,ncf,2,nmy) = zero
        face_corr(n,ncf,3,nmy) = zero
        !
        face_corr(n,ncf,1,nmz) = zero
        face_corr(n,ncf,2,nmz) = zero
        face_corr(n,ncf,3,nmz) = zero
        !
        face_corr(n,ncf,1,nee) = zero
        face_corr(n,ncf,2,nee) = zero
        face_corr(n,ncf,3,nee) = zero
        !
      end do
    end do
    !
    ! Loop over the faces of this cell
    cell_faces: do ncf = 1,num_cell_faces
      !
      ! kf: index of flx-pt local to this face
      ! lf: index of flx-pt within interp%mat
      !
      ! Apply the jump correction to the interior discontinuous gradients
      !
      !$ACC LOOP VECTOR
      face_points: do kf = 1,cf_nfp(ncf)
        !
        k = cf_f0(ncf) + kf
        !
        rnx = facenrm(1,k)
        rny = facenrm(2,k)
        rnz = facenrm(3,k)
        !
        face_order = cf_face_order(ncf)
        face_geom = cf_face_geom(ncf)
        fp_offset = cf_fp_offset(ncf)
        face_rot = cf_face_rot(ncf)
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        !$ACC LOOP SEQ
        do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
               correct(this_geom,this_order)%row_ptr(lf+1)
          !
          n = correct(this_geom,this_order)%row_idx(i)
          a = correct(this_geom,this_order)%row_val(i)
          !
          cor_ec = a * face_jump(k,nec) * cell_invjac(n)
          cor_mx = a * face_jump(k,nmx) * cell_invjac(n)
          cor_my = a * face_jump(k,nmy) * cell_invjac(n)
          cor_mz = a * face_jump(k,nmz) * cell_invjac(n)
          cor_ee = a * face_jump(k,nee) * cell_invjac(n)
          !
          face_corr(n,ncf,1,nec) = face_corr(n,ncf,1,nec) + cor_ec*rnx
          face_corr(n,ncf,2,nec) = face_corr(n,ncf,2,nec) + cor_ec*rny
          face_corr(n,ncf,3,nec) = face_corr(n,ncf,3,nec) + cor_ec*rnz
          !
          face_corr(n,ncf,1,nmx) = face_corr(n,ncf,1,nmx) + cor_mx*rnx
          face_corr(n,ncf,2,nmx) = face_corr(n,ncf,2,nmx) + cor_mx*rny
          face_corr(n,ncf,3,nmx) = face_corr(n,ncf,3,nmx) + cor_mx*rnz
          !
          face_corr(n,ncf,1,nmy) = face_corr(n,ncf,1,nmy) + cor_my*rnx
          face_corr(n,ncf,2,nmy) = face_corr(n,ncf,2,nmy) + cor_my*rny
          face_corr(n,ncf,3,nmy) = face_corr(n,ncf,3,nmy) + cor_my*rnz
          !
          face_corr(n,ncf,1,nmz) = face_corr(n,ncf,1,nmz) + cor_mz*rnx
          face_corr(n,ncf,2,nmz) = face_corr(n,ncf,2,nmz) + cor_mz*rny
          face_corr(n,ncf,3,nmz) = face_corr(n,ncf,3,nmz) + cor_mz*rnz
          !
          face_corr(n,ncf,1,nee) = face_corr(n,ncf,1,nee) + cor_ee*rnx
          face_corr(n,ncf,2,nee) = face_corr(n,ncf,2,nee) + cor_ee*rny
          face_corr(n,ncf,3,nee) = face_corr(n,ncf,3,nee) + cor_ee*rnz
          !
        end do
        !
      end do face_points
    end do cell_faces
    !
    ! Interpolate the face corrected gradient from the interior solution
    ! points to the interface flux points now that we are done with all
    ! the flux points on the current face of the current cell.
    !
    cell_faces2: do ncf = 1,num_cell_faces
      !
      ! kf: index of flx-pt local to this face
      ! lf: index of flx-pt within interp%mat
      !
      !$ACC LOOP VECTOR
      interpolate_face: do kf = 1,cf_nfp(ncf)
        !
        k = cf_f0(ncf) + kf
        sf = cf_sf(ncf)
        face_order = cf_face_order(ncf)
        face_geom = cf_face_geom(ncf)
        fp_offset = cf_fp_offset(ncf)
        face_rot = cf_face_rot(ncf)
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        decdx = zero ; decdy = zero ; decdz = zero
        dmxdx = zero ; dmxdy = zero ; dmxdz = zero
        dmydx = zero ; dmydy = zero ; dmydz = zero
        dmzdx = zero ; dmzdy = zero ; dmzdz = zero
        deedx = zero ; deedy = zero ; deedz = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          decdx = decdx + a * ( dusp(n0+n,1,nec) + face_corr(n,ncf,1,nec) )
          decdy = decdy + a * ( dusp(n0+n,2,nec) + face_corr(n,ncf,2,nec) )
          decdz = decdz + a * ( dusp(n0+n,3,nec) + face_corr(n,ncf,3,nec) )
          !
          dmxdx = dmxdx + a * ( dusp(n0+n,1,nmx) + face_corr(n,ncf,1,nmx) )
          dmxdy = dmxdy + a * ( dusp(n0+n,2,nmx) + face_corr(n,ncf,2,nmx) )
          dmxdz = dmxdz + a * ( dusp(n0+n,3,nmx) + face_corr(n,ncf,3,nmx) )
          !
          dmydx = dmydx + a * ( dusp(n0+n,1,nmy) + face_corr(n,ncf,1,nmy) )
          dmydy = dmydy + a * ( dusp(n0+n,2,nmy) + face_corr(n,ncf,2,nmy) )
          dmydz = dmydz + a * ( dusp(n0+n,3,nmy) + face_corr(n,ncf,3,nmy) )
          !
          dmzdx = dmzdx + a * ( dusp(n0+n,1,nmz) + face_corr(n,ncf,1,nmz) )
          dmzdy = dmzdy + a * ( dusp(n0+n,2,nmz) + face_corr(n,ncf,2,nmz) )
          dmzdz = dmzdz + a * ( dusp(n0+n,3,nmz) + face_corr(n,ncf,3,nmz) )
          !
          deedx = deedx + a * ( dusp(n0+n,1,nee) + face_corr(n,ncf,1,nee) )
          deedy = deedy + a * ( dusp(n0+n,2,nee) + face_corr(n,ncf,2,nee) )
          deedz = deedz + a * ( dusp(n0+n,3,nee) + face_corr(n,ncf,3,nee) )
          !
        end do
        !
        facedusp_array(k,1,nec,sf) = decdx
        facedusp_array(k,2,nec,sf) = decdy
        facedusp_array(k,3,nec,sf) = decdz
        !
        facedusp_array(k,1,nmx,sf) = dmxdx
        facedusp_array(k,2,nmx,sf) = dmxdy
        facedusp_array(k,3,nmx,sf) = dmxdz
        !
        facedusp_array(k,1,nmy,sf) = dmydx
        facedusp_array(k,2,nmy,sf) = dmydy
        facedusp_array(k,3,nmy,sf) = dmydz
        !
        facedusp_array(k,1,nmz,sf) = dmzdx
        facedusp_array(k,2,nmz,sf) = dmzdy
        facedusp_array(k,3,nmz,sf) = dmzdz
        !
        facedusp_array(k,1,nee,sf) = deedx
        facedusp_array(k,2,nee,sf) = deedy
        facedusp_array(k,3,nee,sf) = deedz
        !
      end do interpolate_face
      !
    end do cell_faces2
    !
    ! Copy the total corrected gradient accounting for all cell faces
    ! back into dusp
    !
    !$ACC LOOP VECTOR
    do n = n1,n2
      !
      dusp(n,1,nec) = dusp(n,1,nec) &
                    + face_corr(n-n0,1,1,nec) + face_corr(n-n0,2,1,nec) &
                    + face_corr(n-n0,3,1,nec) + face_corr(n-n0,4,1,nec) &
                    + face_corr(n-n0,5,1,nec) + face_corr(n-n0,6,1,nec)
      !
      dusp(n,2,nec) = dusp(n,2,nec) &
                    + face_corr(n-n0,1,2,nec) + face_corr(n-n0,2,2,nec) &
                    + face_corr(n-n0,3,2,nec) + face_corr(n-n0,4,2,nec) &
                    + face_corr(n-n0,5,2,nec) + face_corr(n-n0,6,2,nec)
      !
      dusp(n,3,nec) = dusp(n,3,nec) &
                    + face_corr(n-n0,1,3,nec) + face_corr(n-n0,2,3,nec) &
                    + face_corr(n-n0,3,3,nec) + face_corr(n-n0,4,3,nec) &
                    + face_corr(n-n0,5,3,nec) + face_corr(n-n0,6,3,nec)
      !
      dusp(n,1,nmx) = dusp(n,1,nmx) &
                    + face_corr(n-n0,1,1,nmx) + face_corr(n-n0,2,1,nmx) &
                    + face_corr(n-n0,3,1,nmx) + face_corr(n-n0,4,1,nmx) &
                    + face_corr(n-n0,5,1,nmx) + face_corr(n-n0,6,1,nmx)
      !
      dusp(n,2,nmx) = dusp(n,2,nmx) &
                    + face_corr(n-n0,1,2,nmx) + face_corr(n-n0,2,2,nmx) &
                    + face_corr(n-n0,3,2,nmx) + face_corr(n-n0,4,2,nmx) &
                    + face_corr(n-n0,5,2,nmx) + face_corr(n-n0,6,2,nmx)
      !
      dusp(n,3,nmx) = dusp(n,3,nmx) &
                    + face_corr(n-n0,1,3,nmx) + face_corr(n-n0,2,3,nmx) &
                    + face_corr(n-n0,3,3,nmx) + face_corr(n-n0,4,3,nmx) &
                    + face_corr(n-n0,5,3,nmx) + face_corr(n-n0,6,3,nmx)
      !
      dusp(n,1,nmy) = dusp(n,1,nmy) &
                    + face_corr(n-n0,1,1,nmy) + face_corr(n-n0,2,1,nmy) &
                    + face_corr(n-n0,3,1,nmy) + face_corr(n-n0,4,1,nmy) &
                    + face_corr(n-n0,5,1,nmy) + face_corr(n-n0,6,1,nmy)
      !
      dusp(n,2,nmy) = dusp(n,2,nmy) &
                    + face_corr(n-n0,1,2,nmy) + face_corr(n-n0,2,2,nmy) &
                    + face_corr(n-n0,3,2,nmy) + face_corr(n-n0,4,2,nmy) &
                    + face_corr(n-n0,5,2,nmy) + face_corr(n-n0,6,2,nmy)
      !
      dusp(n,3,nmy) = dusp(n,3,nmy) &
                    + face_corr(n-n0,1,3,nmy) + face_corr(n-n0,2,3,nmy) &
                    + face_corr(n-n0,3,3,nmy) + face_corr(n-n0,4,3,nmy) &
                    + face_corr(n-n0,5,3,nmy) + face_corr(n-n0,6,3,nmy)
      !
      dusp(n,1,nmz) = dusp(n,1,nmz) &
                    + face_corr(n-n0,1,1,nmz) + face_corr(n-n0,2,1,nmz) &
                    + face_corr(n-n0,3,1,nmz) + face_corr(n-n0,4,1,nmz) &
                    + face_corr(n-n0,5,1,nmz) + face_corr(n-n0,6,1,nmz)
      !
      dusp(n,2,nmz) = dusp(n,2,nmz) &
                    + face_corr(n-n0,1,2,nmz) + face_corr(n-n0,2,2,nmz) &
                    + face_corr(n-n0,3,2,nmz) + face_corr(n-n0,4,2,nmz) &
                    + face_corr(n-n0,5,2,nmz) + face_corr(n-n0,6,2,nmz)
      !
      dusp(n,3,nmz) = dusp(n,3,nmz) &
                    + face_corr(n-n0,1,3,nmz) + face_corr(n-n0,2,3,nmz) &
                    + face_corr(n-n0,3,3,nmz) + face_corr(n-n0,4,3,nmz) &
                    + face_corr(n-n0,5,3,nmz) + face_corr(n-n0,6,3,nmz)
      !
      dusp(n,1,nee) = dusp(n,1,nee) &
                    + face_corr(n-n0,1,1,nee) + face_corr(n-n0,2,1,nee) &
                    + face_corr(n-n0,3,1,nee) + face_corr(n-n0,4,1,nee) &
                    + face_corr(n-n0,5,1,nee) + face_corr(n-n0,6,1,nee)
      !
      dusp(n,2,nee) = dusp(n,2,nee) &
                    + face_corr(n-n0,1,2,nee) + face_corr(n-n0,2,2,nee) &
                    + face_corr(n-n0,3,2,nee) + face_corr(n-n0,4,2,nee) &
                    + face_corr(n-n0,5,2,nee) + face_corr(n-n0,6,2,nee)
      !
      dusp(n,3,nee) = dusp(n,3,nee) &
                    + face_corr(n-n0,1,3,nee) + face_corr(n-n0,2,3,nee) &
                    + face_corr(n-n0,3,3,nee) + face_corr(n-n0,4,3,nee) &
                    + face_corr(n-n0,5,3,nee) + face_corr(n-n0,6,3,nee)
      !
    end do
    !
  end do cell_loop
 !!$OMP END DO
 !!$OMP END PARALLEL
  !
#ifdef TIME_MINIAPPS_FOR_GPU_COMPARISON
  call gfr_timer( timings_end(timings_counter) )
#endif
  !
#ifdef WRITE_DATA_FOR_PROFILING
  if (num_prof_br2_outputs == 0) then
    if (mypnum == glb_root) then
      call output_gold_br2_data(dusp,facedusp_array)
    end if
    num_prof_br2_outputs = 1
  end if
#endif
  !
end subroutine correct_gradients_br2_gpu_3d
!
!###############################################################################
!
subroutine output_br2_data()
  !
  !.. Use Statements ..
  use geovar,            only : nr,n_flxpts,n_solpts
  use geovar,            only : output_cell_dt,input_cell_dt
  use geovar,            only : output_face_dt,input_face_dt
  use geovar,            only : output_fp_dt,input_fp_dt
  use metrics_mod,       only : jac,facenrm
  use correction_mod,    only : output_correct_dt,input_correct_dt
  use flowvar,           only : dusp
  use flowvar,           only : faceusp_array
  use interpolation_mod, only : output_interp_dt,input_interp_dt
  !
  !.. Local Scalars ..
  integer :: io,l,m,n,ierr
  logical(lk) :: success
  !
continue
  !
  open (newunit=io,file="br2_cbias.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_flxpts
  write (io) (cbias(n), n=1,n_flxpts)
  close (io)
  !
  open (newunit=io,file="br2_dusp.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) (((dusp(n,l,m), n=1,n_solpts), l=1,nr), m=1,nq)
  close (io)
  !
  open (newunit=io,file="br2_faceusp.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_flxpts
  write (io) (((faceusp_array(n,m,l), n=1,n_flxpts), m=1,nq), l=1,2)
  close (io)
  !
  open (newunit=io,file="br2_jac.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) (jac(n), n=1,n_solpts)
  close (io)
  !
  open (newunit=io,file="br2_facenrm.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_flxpts
  write (io) ((facenrm(l,n), l=1,nr), n=1,n_flxpts)
  close (io)
  !
  open (newunit=io,file="br2_cell.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  call output_cell_dt(io)
  close (io)
  !
  open (newunit=io,file="br2_face.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  call output_face_dt(io)
  close (io)
  !
  open (newunit=io,file="br2_fp.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  call output_fp_dt(io)
  close (io)
  !
  open (newunit=io,file="br2_correct.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  call output_correct_dt(io)
  close (io)
  !
  open (newunit=io,file="br2_interp.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  call output_interp_dt(io)
  close (io)
  !
  ! Check to make sure the derived types were output correctly
  !
  open (newunit=io,file="br2_cell.dat",status="old",action="read", &
                   form="unformatted",iostat=ierr)
  call input_cell_dt(io,success)
  close (io)
  if (success) then
    write (iout,*) "cell output correctly!"
  else
    write (iout,*) "ERROR: cell NOT output correctly"
  end if
  !
  open (newunit=io,file="br2_face.dat",status="old",action="read", &
                   form="unformatted",iostat=ierr)
  call input_face_dt(io,success)
  close (io)
  if (success) then
    write (iout,*) "face output correctly!"
  else
    write (iout,*) "ERROR: face NOT output correctly"
  end if
  !
  open (newunit=io,file="br2_fp.dat",status="old",action="read", &
                   form="unformatted",iostat=ierr)
  call input_fp_dt(io,success)
  close (io)
  if (success) then
    write (iout,*) "fp output correctly!"
  else
    write (iout,*) "ERROR: fp NOT output correctly"
  end if
  !
  open (newunit=io,file="br2_correct.dat",status="old",action="read", &
                   form="unformatted",iostat=ierr)
  call input_correct_dt(io,success)
  close (io)
  if (success) then
    write (iout,*) "correct output correctly!"
  else
    write (iout,*) "ERROR: correct NOT output correctly"
  end if
  !
  open (newunit=io,file="br2_interp.dat",status="old",action="read", &
                   form="unformatted",iostat=ierr)
  call input_interp_dt(io,success)
  close (io)
  if (success) then
    write (iout,*) "interp output correctly!"
  else
    write (iout,*) "ERROR: interp NOT output correctly"
  end if
  !
end subroutine output_br2_data
!
!###############################################################################
!
subroutine output_gold_br2_data(dusp,facedusp)
  !
  !.. Use Statements ..
  use geovar, only : nr,n_flxpts,n_solpts
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: dusp(:,:,:)
  real(wp), intent(in) :: facedusp(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: io,k,l,m,n,ierr
  !
continue
  !
  open (newunit=io,file="br2_gold_dusp.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) (((dusp(n,l,m), n=1,n_solpts), l=1,nr), m=1,nq)
  close (io)
  !
  open (newunit=io,file="br2_gold_facedusp.dat",status="replace", &
                   action="write",form="unformatted",iostat=ierr)
  write (io) n_flxpts
  write (io) ((((facedusp(n,l,m,k), n=1,n_flxpts), l=1,nr), m=1,nq), k=1,2)
  close (io)
  !
end subroutine output_gold_br2_data
!
!###############################################################################
!
subroutine correct_gradients_br2_gpu_2d()
  !
  use order_mod,         only : maxSP
  use geovar,            only : n_flxpts,ncell,cell,fp,face
  use metrics_mod,       only : jac,facenrm
  use correction_mod,    only : correct
  use flowvar,           only : dusp
  use flowvar,           only : faceusp_array
  use flowvar,           only : facedusp_array
  use interpolation_mod, only : interp
  !
  !.. Local Scalars ..
  integer :: i,n,nc,n1,n2,np,k,n0
  integer :: nf,kf,lf,sf,ncf,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,rnx,rny
  real(wp) :: decdx,dmxdx,dmydx,deedx
  real(wp) :: decdy,dmxdy,dmydy,deedy
  real(wp) :: cor_ec,cor_mx,cor_my,cor_ee
  !
  !.. Local Arrays ..
  real(wp) :: face_corr(1:maxSP,1:4,1:2,1:nq)
  !
  integer :: array(1:1) ! DONT REMOVE THIS
  !
  integer  :: cf_f0(1:4)
  integer  :: cf_sf(1:4)
  integer  :: cf_nfp(1:4)
  integer  :: cf_face_order(1:4)
  integer  :: cf_face_geom(1:4)
  integer  :: cf_fp_offset(1:4)
  integer  :: cf_face_rot(1:4)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_br2_gpu_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP &
  !$ACC&              PRESENT(faceusp_array,cbias,face_jump)
  do n = 1,n_flxpts
    a = cbias(n)
    face_jump(n,nec) = a*( faceusp_array(n,nec,2) - faceusp_array(n,nec,1) )
    face_jump(n,nmx) = a*( faceusp_array(n,nmx,2) - faceusp_array(n,nmx,1) )
    face_jump(n,nmy) = a*( faceusp_array(n,nmy,2) - faceusp_array(n,nmy,1) )
    face_jump(n,nee) = a*( faceusp_array(n,nee,2) - faceusp_array(n,nee,1) )
  end do
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              WORKER &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facedusp_array,facenrm) &
  !$ACC&              PRESENT(jac) &
  !$ACC&              PRESENT(face) &
  !$ACC&              PRESENT(cell) &
  !$ACC&              PRESENT(fp) &
  !$ACC&              PRESENT(correct) &
  !$ACC&              PRESENT(interp) &
  !$ACC&              PRESENT(dusp) &
  !$ACC&              PRESENT(cl) &
  !$ACC&              PRESENT(face_jump) &
  !$ACC&              PRIVATE(face_corr,cf_f0,cf_sf, &
  !$ACC&                      cf_face_geom,cf_fp_offset,cf_face_rot, &
  !$ACC&                      cf_nfp,cf_face_order,array)
  cell_loop: do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    n0 = n1-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now Interpolate the interior solution gradients to the flux points
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      !
      cf_sf(ncf) = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      cf_face_geom(ncf) = cell(nc)%face(ncf)%geom   ! geom  of current face
      cf_face_order(ncf) = cell(nc)%face(ncf)%order ! order of current face
      cf_fp_offset(ncf) = cell(nc)%face(ncf)%fp_offset
      cf_face_rot(ncf) = cell(nc)%face(ncf)%rotation
      !
      cf_f0(ncf) = face(nf)%beg_fp - 1
      cf_nfp(ncf) = face(nf)%end_fp - cf_f0(ncf)
      !
    end do
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do ncf = 1,4
      do n = 1,np
        !
        face_corr(n,ncf,1,nec) = zero
        face_corr(n,ncf,2,nec) = zero
        !
        face_corr(n,ncf,1,nmx) = zero
        face_corr(n,ncf,2,nmx) = zero
        !
        face_corr(n,ncf,1,nmy) = zero
        face_corr(n,ncf,2,nmy) = zero
        !
        face_corr(n,ncf,1,nee) = zero
        face_corr(n,ncf,2,nee) = zero
        !
      end do
    end do
    !
    ! Loop over the faces of this cell
    cell_faces: do ncf = 1,num_cell_faces
      !
      ! kf: index of flx-pt local to this face
      ! lf: index of flx-pt within interp%mat
      !
      ! Apply the jump correction to the interior discontinuous gradients
      !
      !$ACC LOOP VECTOR
      face_points: do kf = 1,cf_nfp(ncf)
        !
        k = cf_f0(ncf) + kf
        !
        rnx = facenrm(1,k)
        rny = facenrm(2,k)
        !
        face_order = cf_face_order(ncf)
        face_geom = cf_face_geom(ncf)
        fp_offset = cf_fp_offset(ncf)
        face_rot = cf_face_rot(ncf)
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        !$ACC LOOP SEQ
        do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
               correct(this_geom,this_order)%row_ptr(lf+1)
          !
          n = correct(this_geom,this_order)%row_idx(i)
          a = correct(this_geom,this_order)%row_val(i)
          !
          cor_ec = a * face_jump(k,nec) * cl(this_geom) / jac(n0+n)
          cor_mx = a * face_jump(k,nmx) * cl(this_geom) / jac(n0+n)
          cor_my = a * face_jump(k,nmy) * cl(this_geom) / jac(n0+n)
          cor_ee = a * face_jump(k,nee) * cl(this_geom) / jac(n0+n)
          !
          face_corr(n,ncf,1,nec) = face_corr(n,ncf,1,nec) + cor_ec*rnx
          face_corr(n,ncf,2,nec) = face_corr(n,ncf,2,nec) + cor_ec*rny
          !
          face_corr(n,ncf,1,nmx) = face_corr(n,ncf,1,nmx) + cor_mx*rnx
          face_corr(n,ncf,2,nmx) = face_corr(n,ncf,2,nmx) + cor_mx*rny
          !
          face_corr(n,ncf,1,nmy) = face_corr(n,ncf,1,nmy) + cor_my*rnx
          face_corr(n,ncf,2,nmy) = face_corr(n,ncf,2,nmy) + cor_my*rny
          !
          face_corr(n,ncf,1,nee) = face_corr(n,ncf,1,nee) + cor_ee*rnx
          face_corr(n,ncf,2,nee) = face_corr(n,ncf,2,nee) + cor_ee*rny
          !
        end do
        !
      end do face_points
    end do cell_faces
    !
    ! Interpolate the face corrected gradient from the interior solution
    ! points to the interface flux points now that we are done with all
    ! the flux points on the current face of the current cell.
    !
    cell_faces2: do ncf = 1,num_cell_faces
      !
      ! kf: index of flx-pt local to this face
      ! lf: index of flx-pt within interp%mat
      !
      !$ACC LOOP VECTOR
      interpolate_face: do kf = 1,cf_nfp(ncf)
        !
        k = cf_f0(ncf) + kf
        sf = cf_sf(ncf)
        face_order = cf_face_order(ncf)
        face_geom = cf_face_geom(ncf)
        fp_offset = cf_fp_offset(ncf)
        face_rot = cf_face_rot(ncf)
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        decdx = zero ; decdy = zero
        dmxdx = zero ; dmxdy = zero
        dmydx = zero ; dmydy = zero
        deedx = zero ; deedy = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          decdx = decdx + a * ( dusp(n0+n,1,nec) + face_corr(n,ncf,1,nec) )
          decdy = decdy + a * ( dusp(n0+n,2,nec) + face_corr(n,ncf,2,nec) )
          !
          dmxdx = dmxdx + a * ( dusp(n0+n,1,nmx) + face_corr(n,ncf,1,nmx) )
          dmxdy = dmxdy + a * ( dusp(n0+n,2,nmx) + face_corr(n,ncf,2,nmx) )
          !
          dmydx = dmydx + a * ( dusp(n0+n,1,nmy) + face_corr(n,ncf,1,nmy) )
          dmydy = dmydy + a * ( dusp(n0+n,2,nmy) + face_corr(n,ncf,2,nmy) )
          !
          deedx = deedx + a * ( dusp(n0+n,1,nee) + face_corr(n,ncf,1,nee) )
          deedy = deedy + a * ( dusp(n0+n,2,nee) + face_corr(n,ncf,2,nee) )
          !
        end do
        !
        facedusp_array(k,1,nec,sf) = decdx
        facedusp_array(k,2,nec,sf) = decdy
        !
        facedusp_array(k,1,nmx,sf) = dmxdx
        facedusp_array(k,2,nmx,sf) = dmxdy
        !
        facedusp_array(k,1,nmy,sf) = dmydx
        facedusp_array(k,2,nmy,sf) = dmydy
        !
        facedusp_array(k,1,nee,sf) = deedx
        facedusp_array(k,2,nee,sf) = deedy
        !
      end do interpolate_face
      !
    end do cell_faces2
    !
    ! Copy the total corrected gradient accounting for all cell faces
    ! back into dusp
    !
    !$ACC LOOP VECTOR
    do n = n1,n2
      !
      dusp(n,1,nec) = dusp(n,1,nec) &
                    + face_corr(n-n0,1,1,nec) + face_corr(n-n0,2,1,nec) &
                    + face_corr(n-n0,3,1,nec) + face_corr(n-n0,4,1,nec)
      !
      dusp(n,2,nec) = dusp(n,2,nec) &
                    + face_corr(n-n0,1,2,nec) + face_corr(n-n0,2,2,nec) &
                    + face_corr(n-n0,3,2,nec) + face_corr(n-n0,4,2,nec)
      !
      dusp(n,1,nmx) = dusp(n,1,nmx) &
                    + face_corr(n-n0,1,1,nmx) + face_corr(n-n0,2,1,nmx) &
                    + face_corr(n-n0,3,1,nmx) + face_corr(n-n0,4,1,nmx)
      !
      dusp(n,2,nmx) = dusp(n,2,nmx) &
                    + face_corr(n-n0,1,2,nmx) + face_corr(n-n0,2,2,nmx) &
                    + face_corr(n-n0,3,2,nmx) + face_corr(n-n0,4,2,nmx)
      !
      dusp(n,1,nmy) = dusp(n,1,nmy) &
                    + face_corr(n-n0,1,1,nmy) + face_corr(n-n0,2,1,nmy) &
                    + face_corr(n-n0,3,1,nmy) + face_corr(n-n0,4,1,nmy)
      !
      dusp(n,2,nmy) = dusp(n,2,nmy) &
                    + face_corr(n-n0,1,2,nmy) + face_corr(n-n0,2,2,nmy) &
                    + face_corr(n-n0,3,2,nmy) + face_corr(n-n0,4,2,nmy)
      !
      dusp(n,1,nee) = dusp(n,1,nee) &
                    + face_corr(n-n0,1,1,nee) + face_corr(n-n0,2,1,nee) &
                    + face_corr(n-n0,3,1,nee) + face_corr(n-n0,4,1,nee)
      !
      dusp(n,2,nee) = dusp(n,2,nee) &
                    + face_corr(n-n0,1,2,nee) + face_corr(n-n0,2,2,nee) &
                    + face_corr(n-n0,3,2,nee) + face_corr(n-n0,4,2,nee)
      !
    end do
    !
  end do cell_loop
  !
end subroutine correct_gradients_br2_gpu_2d
!
!###############################################################################
!
subroutine correct_gradients_br2_gpu_1d()
  !
  use order_mod,         only : maxSP
  use geovar,            only : n_flxpts,ncell,cell,fp,face
  use metrics_mod,       only : jac,facenrm
  use correction_mod,    only : correct
  use flowvar,           only : dusp
  use flowvar,           only : faceusp_array
  use flowvar,           only : facedusp_array
  use interpolation_mod, only : interp
  !
  !.. Local Scalars ..
  integer :: i,n,nc,n1,n2,np,k,n0
  integer :: nf,kf,lf,sf,ncf,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,rnx
  real(wp) :: decdx,dmxdx,deedx
  real(wp) :: cor_ec,cor_mx,cor_my,cor_ee
  !
  !.. Local Arrays ..
  real(wp) :: face_corr(1:maxSP,1:2,1:2,1:nq)
  !
  integer :: array(1:1) ! DONT REMOVE THIS
  !
  integer  :: cf_f0(1:2)
  integer  :: cf_sf(1:2)
  integer  :: cf_nfp(1:2)
  integer  :: cf_face_order(1:2)
  integer  :: cf_face_geom(1:2)
  integer  :: cf_fp_offset(1:2)
  integer  :: cf_face_rot(1:2)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_br2_gpu_1d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP &
  !$ACC&              PRESENT(faceusp_array,cbias,face_jump)
  do n = 1,n_flxpts
    a = cbias(n)
    face_jump(n,nec) = a*( faceusp_array(n,nec,2) - faceusp_array(n,nec,1) )
    face_jump(n,nmx) = a*( faceusp_array(n,nmx,2) - faceusp_array(n,nmx,1) )
    face_jump(n,nee) = a*( faceusp_array(n,nee,2) - faceusp_array(n,nee,1) )
  end do
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              WORKER &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facedusp_array,facenrm) &
  !$ACC&              PRESENT(jac) &
  !$ACC&              PRESENT(face) &
  !$ACC&              PRESENT(cell) &
  !$ACC&              PRESENT(fp) &
  !$ACC&              PRESENT(correct) &
  !$ACC&              PRESENT(interp) &
  !$ACC&              PRESENT(dusp) &
  !$ACC&              PRESENT(cl) &
  !$ACC&              PRESENT(face_jump) &
  !$ACC&              PRIVATE(face_corr,cf_f0,cf_sf, &
  !$ACC&                      cf_face_geom,cf_fp_offset,cf_face_rot, &
  !$ACC&                      cf_nfp,cf_face_order,array)
  cell_loop: do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    n0 = n1-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now Interpolate the interior solution gradients to the flux points
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      !
      cf_sf(ncf) = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      cf_face_geom(ncf) = cell(nc)%face(ncf)%geom   ! geom  of current face
      cf_face_order(ncf) = cell(nc)%face(ncf)%order ! order of current face
      cf_fp_offset(ncf) = cell(nc)%face(ncf)%fp_offset
      cf_face_rot(ncf) = cell(nc)%face(ncf)%rotation
      !
      cf_f0(ncf) = face(nf)%beg_fp - 1
      cf_nfp(ncf) = face(nf)%end_fp - cf_f0(ncf)
      !
    end do
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do ncf = 1,2
      do n = 1,np
        !
        face_corr(n,ncf,1,nec) = zero
        face_corr(n,ncf,1,nmx) = zero
        face_corr(n,ncf,1,nee) = zero
        !
      end do
    end do
    !
    ! Loop over the faces of this cell
    cell_faces: do ncf = 1,num_cell_faces
      !
      ! kf: index of flx-pt local to this face
      ! lf: index of flx-pt within interp%mat
      !
      ! Apply the jump correction to the interior discontinuous gradients
      !
      !$ACC LOOP VECTOR
      face_points: do kf = 1,cf_nfp(ncf)
        !
        k = cf_f0(ncf) + kf
        !
        rnx = facenrm(1,k)
        !
        face_order = cf_face_order(ncf)
        face_geom = cf_face_geom(ncf)
        fp_offset = cf_fp_offset(ncf)
        face_rot = cf_face_rot(ncf)
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        !$ACC LOOP SEQ
        do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
               correct(this_geom,this_order)%row_ptr(lf+1)
          !
          n = correct(this_geom,this_order)%row_idx(i)
          a = correct(this_geom,this_order)%row_val(i)
          !
          cor_ec = a * face_jump(k,nec) * cl(this_geom) / jac(n0+n)
          cor_mx = a * face_jump(k,nmx) * cl(this_geom) / jac(n0+n)
          cor_ee = a * face_jump(k,nee) * cl(this_geom) / jac(n0+n)
          !
          face_corr(n,ncf,1,nec) = face_corr(n,ncf,1,nec) + cor_ec*rnx
          face_corr(n,ncf,1,nmx) = face_corr(n,ncf,1,nmx) + cor_mx*rnx
          face_corr(n,ncf,1,nee) = face_corr(n,ncf,1,nee) + cor_ee*rnx
          !
        end do
        !
      end do face_points
    end do cell_faces
    !
    ! Interpolate the face corrected gradient from the interior solution
    ! points to the interface flux points now that we are done with all
    ! the flux points on the current face of the current cell.
    !
    cell_faces2: do ncf = 1,num_cell_faces
      !
      ! kf: index of flx-pt local to this face
      ! lf: index of flx-pt within interp%mat
      !
      !$ACC LOOP VECTOR
      interpolate_face: do kf = 1,cf_nfp(ncf)
        !
        k = cf_f0(ncf) + kf
        sf = cf_sf(ncf)
        face_order = cf_face_order(ncf)
        face_geom = cf_face_geom(ncf)
        fp_offset = cf_fp_offset(ncf)
        face_rot = cf_face_rot(ncf)
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        decdx = zero
        dmxdx = zero
        deedx = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          decdx = decdx + a * ( dusp(n0+n,1,nec) + face_corr(n,ncf,1,nec) )
          dmxdx = dmxdx + a * ( dusp(n0+n,1,nmx) + face_corr(n,ncf,1,nmx) )
          deedx = deedx + a * ( dusp(n0+n,1,nee) + face_corr(n,ncf,1,nee) )
          !
        end do
        !
        facedusp_array(k,1,nec,sf) = decdx
        facedusp_array(k,1,nmx,sf) = dmxdx
        facedusp_array(k,1,nee,sf) = deedx
        !
      end do interpolate_face
      !
    end do cell_faces2
    !
    ! Copy the total corrected gradient accounting for all cell faces
    ! back into dusp
    !
    !$ACC LOOP VECTOR
    do n = n1,n2
      !
      dusp(n,1,nec) = dusp(n,1,nec) &
                    + face_corr(n-n0,1,1,nec) + face_corr(n-n0,2,1,nec)
      !
      dusp(n,1,nmx) = dusp(n,1,nmx) &
                    + face_corr(n-n0,1,1,nmx) + face_corr(n-n0,2,1,nmx)
      !
      dusp(n,1,nee) = dusp(n,1,nee) &
                    + face_corr(n-n0,1,1,nee) + face_corr(n-n0,2,1,nee)
      !
    end do
    !
  end do cell_loop
  !
end subroutine correct_gradients_br2_gpu_1d
!
!###############################################################################
!
subroutine correct_gradients_br2_gpu_nd(mbeg,mend)
  !
  use order_mod,         only : maxSP
  use geovar,            only : nr,n_flxpts,ncell,cell,fp,face
  use metrics_mod,       only : jac,facenrm
  use correction_mod,    only : correct
  use flowvar,           only : dusp
  use flowvar,           only : faceusp_array
  use flowvar,           only : facedusp_array
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  integer, intent(in) :: mbeg,mend
  !
  !.. Local Scalars ..
  integer :: i,n,nc,n1,n2,np,l,m,f0,k,n0
  integer :: nf,kf,lf,sf,ncf,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,cor_m,val
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxSP,1:nr,1:nq,1:6) :: face_corr
  !
  integer :: array(1:1) ! DONT REMOVE THIS
  !
  integer  :: cf_f0(1:6)
  integer  :: cf_sf(1:6)
  integer  :: cf_nfp(1:6)
  integer  :: cf_face_order(1:6)
  integer  :: cf_face_geom(1:6)
  integer  :: cf_fp_offset(1:6)
  integer  :: cf_face_rot(1:6)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_br2_gpu_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
 !!$OMP PARALLEL &
 !!$OMP&         DEFAULT(SHARED) &
 !!$OMP&         PRIVATE(nc,n1,n2,np,n0,this_geom,this_order,num_cell_faces, &
 !!$OMP&                 ncf,nf,cf_sf,cf_face_geom,cf_face_order, &
 !!$OMP&                 cf_fp_offset,cf_face_rot,cf_f0,cf_nfp,face_corr, &
 !!$OMP&                 kf,m,l,k,face_order,face_geom,fp_offset,face_rot, &
 !!$OMP&                 lf,i,n,a,cor_m,f0,sf,val)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(faceusp_array,cbias,face_jump)
 !!$OMP DO
  do m = mbeg,mend
    do n = 1,n_flxpts
      face_jump(n,m) = cbias(n)*( faceusp_array(n,m,2) - faceusp_array(n,m,1) )
    end do
  end do
 !!$OMP END DO
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              WORKER &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facedusp_array,facenrm) &
  !$ACC&              PRESENT(jac) &
  !$ACC&              PRESENT(fp) &
  !$ACC&              PRESENT(face) &
  !$ACC&              PRESENT(correct) &
  !$ACC&              PRESENT(cell) &
  !$ACC&              PRESENT(interp) &
  !$ACC&              PRESENT(dusp) &
  !$ACC&              PRESENT(cl) &
  !$ACC&              PRESENT(face_jump) &
  !$ACC&              PRIVATE(face_corr,cf_f0,cf_sf, &
  !$ACC&                      cf_face_geom,cf_fp_offset,cf_face_rot, &
  !$ACC&                      cf_nfp,cf_face_order,array)
 !!$OMP DO
  cell_loop: do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    n0 = n1-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now Interpolate the interior solution gradients to the flux points
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      !
      cf_sf(ncf) = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      cf_face_geom(ncf) = cell(nc)%face(ncf)%geom   ! geom  of current face
      cf_face_order(ncf) = cell(nc)%face(ncf)%order ! order of current face
      cf_fp_offset(ncf) = cell(nc)%face(ncf)%fp_offset
      cf_face_rot(ncf) = cell(nc)%face(ncf)%rotation
      !
      cf_f0(ncf) = face(nf)%beg_fp - 1
      cf_nfp(ncf) = face(nf)%end_fp - cf_f0(ncf)
      !
    end do
    !
    face_corr(:,:,:,:) = zero
    !
    ! Loop over the faces of this cell
    cell_faces: do ncf = 1,num_cell_faces
      !
      ! kf: index of flx-pt local to this face
      ! lf: index of flx-pt within interp%mat
      !
      ! Apply the jump correction to the interior discontinuous gradients
      !
      !$ACC LOOP VECTOR COLLAPSE(3)
      face_points: do kf = 1,cf_nfp(ncf)
        do m = mbeg,mend
          do l = 1,nr
            !
            k = cf_f0(ncf) + kf
            !
            face_order = cf_face_order(ncf)
            face_geom = cf_face_geom(ncf)
            fp_offset = cf_fp_offset(ncf)
            face_rot = cf_face_rot(ncf)
            lf = fp_offset + fp%geord(face_geom,face_order)% &
                                              rot(face_rot)%pts(kf)
            !
            do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
                   correct(this_geom,this_order)%row_ptr(lf+1)
              !
              n = correct(this_geom,this_order)%row_idx(i)
              a = correct(this_geom,this_order)%row_val(i)
              !
              cor_m = a * face_jump(k,m) * cl(this_geom) / jac(n0+n)
              !
              face_corr(n,l,m,ncf) = face_corr(n,l,m,ncf) + cor_m * facenrm(l,k)
              !
            end do
            !
          end do
        end do
      end do face_points
      !
    end do cell_faces
    !
    ! Interpolate the face corrected gradient from the interior solution
    ! points to the interface flux points now that we are done with all
    ! the flux points on the current face of the current cell.
    !
    cell_faces2: do ncf = 1,num_cell_faces
      !
      ! kf: index of flx-pt local to this face
      ! lf: index of flx-pt within interp%mat
      !
      !$ACC LOOP VECTOR COLLAPSE(3)
      interpolate_face: do kf = 1,cf_nfp(ncf)
        do m = mbeg,mend
          do l = 1,nr
            !
            f0 = cf_f0(ncf)
            sf = cf_sf(ncf)
            face_order = cf_face_order(ncf)
            face_geom = cf_face_geom(ncf)
            fp_offset = cf_fp_offset(ncf)
            face_rot = cf_face_rot(ncf)
            lf = fp_offset + fp%geord(face_geom,face_order)% &
                                              rot(face_rot)%pts(kf)
            !
            val = zero
            !
            do i = interp(this_geom,this_order)% &
                             toFace(face_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)% &
                             toFace(face_order)%row_ptr(lf+1)
              !
              n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
              a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
              !
              val = val + a * ( dusp(n0+n,l,m) + face_corr(n,l,m,ncf) )
              !
            end do
            !
            facedusp_array(f0+kf,l,m,sf) = val
            !
          end do
        end do
      end do interpolate_face
      !
    end do cell_faces2
    !
    ! Copy the total corrected gradient accounting for all cell faces
    ! back into dusp
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do m = mbeg,mend
      do l = 1,nr
        do n = n1,n2
          dusp(n,l,m) = dusp(n,l,m) &
                      + face_corr(n-n0,l,m,1) + face_corr(n-n0,l,m,2) &
                      + face_corr(n-n0,l,m,3) + face_corr(n-n0,l,m,4) &
                      + face_corr(n-n0,l,m,5) + face_corr(n-n0,l,m,6)
        end do
      end do
    end do
    !
  end do cell_loop
 !!$OMP END DO
 !!$OMP END PARALLEL
  !
end subroutine correct_gradients_br2_gpu_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine write_2d_array(io_start,iter,nst,array)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: io_start,iter,nst
  real(wp), intent(in) :: array(:,:)
  !
  !.. Local Scalars ..
  integer :: io,n,m
  !
continue
  !
  !$ACC UPDATE SELF(array)
  io = io_start + 10*(iter-1) + nst
  do n = 1,size(array,dim=1)
    write (io,100) (array(n,m), m=1,size(array,dim=2))
  end do
  100 format (5es16.8)
  !$ACC UPDATE DEVICE(array)
  !
end subroutine write_2d_array
!
!###############################################################################
!
subroutine write_3d_array(io_start,iter,nst,array)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: io_start,iter,nst
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Local Scalars ..
  integer :: io,n,m,l
  !
continue
  !
  !$ACC UPDATE SELF(array)
  io = io_start + 10*(iter-1) + nst
  do n = 1,size(array,dim=1)
    write (io,100) ((array(n,l,m), l=1,size(array,dim=2)), &
                                   m=1,size(array,dim=3))
  end do
  100 format (15es11.3)
  !$ACC UPDATE DEVICE(array)
  !
end subroutine write_3d_array
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! START OF OLD FLUX SUBROUTINES
!
!###############################################################################
!###############################################################################
!###############################################################################
!
#ifdef ENABLE_OLD_FLUX_SUBROUTINES
!!subroutine flux_gradient_at_solpts_visc(time_step,rk_step)
!!!
!!!... this subroutine get derivative of fluxes (no interaction)
!!!... for each cell, get f, g, ftd, gtd at the solution points
!!!
!!!.. Use Statements ..
!!use order_mod,         only : maxSP,maxFP,geom_solpts
!!use geovar,            only : nr,ncell
!!use geovar,            only : cell,face,fp
!!use interpolation_mod, only : interp
!!use derivatives_mod,   only : deriv
!!use metrics_mod,       only : metrics_dt,facenrm
!!use flowvar,           only : usp,gpsp,dusp,residual
!!use ovar,              only : governing_equations
!!use correction_mod,    only : correct
!!!
!!use visc_flux_mod, only : visc_flx_cv
!!use invs_flux_mod, only : invs_flx_cv
!!!
!!!.. Formal Arguments ..
!!integer, intent(in) :: time_step
!!integer, intent(in) :: rk_step
!!!
!!!.. Local Scalars ..
!!integer :: i,j,k,l,m,n,nc,n1,n2,np
!!integer :: nf,sf,kf,lf,f0
!!integer :: this_geom,this_order
!!integer :: face_geom,face_order
!!integer :: fp_offset,face_rot,nfp,ncf
!!real(wp) :: a,ds,val,fval,correct_n
!!!
!!!.. Local Arrays ..
!!real(wp), dimension(1:nq,1:nr) :: vflx
!!real(wp), dimension(1:nq,1:nr) :: iflx
!!real(wp), dimension(1:nq,1:maxSP) :: lcl_res
!!real(wp), dimension(1:maxSP,1:nq,1:nr) :: rsflx
!!real(wp), dimension(1:maxSP,1:nq,1:nr) :: xyflx
!!integer, dimension(1:maxFP) :: fpidx
!!real(wp), dimension(1:nq) :: flx
!!!
!!!.. Local Parameters ..
!!character(len=*), parameter :: pname = "flux_gradient_at_solpts_visc"
!!!
!!continue
!!!
!!call debug_timer(entering_procedure,pname)
!!!
!!!$ACC DATA CREATE
!!vflx(:,:) = zero
!!!
!!! Loop over the cells and get the derivative of the
!!! fluxes at each solution point within each cell
!!!
! !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
! !$ACC&                   PRESENT(deriv) &
! !$ACC&                   PRESENT(correct) &
! !$ACC&                   PRESENT(fp) &
! !$ACC&                   PRESENT(cell) &
! !$ACC&                   PRESENT(interp) &
! !$ACC&                   PRESENT(residual) &
! !$ACC&                   PRESENT(dusp) &
! !$ACC&                   PRESENT(usp) &
! !$ACC&                   PRIVATE(iflx,vflx,lcl_res,rsflx,xyflx,fpidx,flx)
!!!$ACC&                   PRIVATE(lcl_res,rsflx,xyflx,fpidx,flx)
!!do nc = 1,ncell
!!  !
!!  n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
!!  n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
!!  np = n2-n1+1         ! number of sol-pts in this cell
!!  !
!!  this_geom = cell(nc)%geom   ! geometry type of this cell
!!  this_order = cell(nc)%order ! solution order of this cell
!!  !
!  !!$ACC LOOP VECTOR COLLAPSE(2)
!  !!$ACC LOOP WORKER
!   !$ACC LOOP SEQ
!!  do k = 1,np
!!    !
!!    n = k + n1-1
!!    !
!!    ! Inviscid fluxes
!!    !
!    !!$ACC LOOP SEQ
!!    do l = 1,nr
!!      ! l-directional flux
!!      call invs_flx_cv( usp(n,:) , gpsp(n,:) , l , iflx(:,l) )
!!    end do
!!    !
!!    ! Viscous fluxes
!!    !
!!    if (governing_equations == NavierStokes_Eqns) then
!!      call visc_flx_cv( usp(n,:) , gpsp(n,:) , dusp(n,:,:) , vflx )
!!    end if
!!    !
!!    ! Total fluxes
!!    !
!    !!$ACC LOOP SEQ
!!    do l = 1,nr
!      !!$ACC LOOP SEQ
!!      do m = 1,nq
!!        xyflx(k,m,l) = iflx(m,l) - vflx(m,l)
!!      end do
!!    end do
!!    !
!!  end do
!!  !
!!  ! Transform the fluxes in the x-y physical
!!  ! domain to the r-s computational domain
!!  !
!   !$ACC LOOP VECTOR COLLAPSE(3)
!!  do m = 1,nq
!!    do l = 1,nr
!!      do k = 1,np
!!        val = zero
!!        do i = 1,nr
!!          val = val + xyflx(k,m,i)*metrics_dt(nc)%met(k,l,i)
!!        end do
!!        rsflx(k,m,l) = val * metrics_dt(nc)%jac(k)
!!      end do
!!    end do
!!  end do
!!  !
!!  ! Compute the divergence of the discontinuous fluxes
!!  !
!!  ! Zero the residual array for the current cell
!!  !
!  !!$ACC LOOP VECTOR COLLAPSE(2)
!! !do n = 1,maxSP
!! !  do m = 1,nq
!! !    lcl_res(m,n) = zero
!! !  end do
!! !end do
!!  !
!!  ! Loop through each coordinate direction and add the divergence of
!!  ! the flux to the local residual array for the current cell
!!  !
!   !$ACC LOOP VECTOR COLLAPSE(2)
!!  do k = 1,np
!!    do m = 1,nq
!!      val = zero
!!      do l = 1,nr
!!        do i = deriv(this_geom,this_order)%d(l)%row_ptr(k)+1, &
!!               deriv(this_geom,this_order)%d(l)%row_ptr(k+1)
!!          !
!!          n = deriv(this_geom,this_order)%d(l)%row_idx(i)
!!          a = deriv(this_geom,this_order)%d(l)%row_val(i)
!!          !
!!          val = val + a * rsflx(n,m,l)
!!          !
!!        end do
!!      end do
!!      lcl_res(m,k) = val
!!    end do
!!  end do
!!  !
!!  ! Copy the local residual for the current cell into the residual array
!!  !
!   !$ACC LOOP VECTOR COLLAPSE(2)
!!  do m = 1,nq
!!    do n = n1,n2
!!      residual(n,m) = lcl_res(m,n-n1+1)
!!    end do
!!  end do
!!  !
!!  ! Interpolate the fluxes to the interface flux points and dot this
!!  ! interpolated vector with the face normal at each flux point
!!  ! NOTE: We dont need to exchange this information for communication and
!!  !       periodic boundary faces since this flux is only used within the
!!  !       current cell.
!!  !
!!  ! Loop over the faces of this cell
!!  !
!!  do ncf = 1,size(cell(nc)%face)
!!    !
!!    nf = cell(nc)%face(ncf)%idx  ! index of the current face
!!    sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
!!    ds = cl(this_geom) * merge(-one,one,sf==1)
!!    f0 = face(nf)%beg_fp-1
!!    !
!!    face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
!!    face_order = cell(nc)%face(ncf)%order ! order of current face
!!    fp_offset = cell(nc)%face(ncf)%fp_offset
!!    face_rot = cell(nc)%face(ncf)%rotation
!!    !
!!    nfp = geom_solpts(face_geom,face_order)
!!    fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
!!                                  rot(face_rot)%pts(1:nfp)
!!    !
!     !$ACC LOOP VECTOR COLLAPSE(2)
!!    do kf = 1,nfp
!!      !
!!      ! kf: index of the flx-pt local to this face
!!      ! lf: index of flx-pt within interp%mat
!!      !
!!      do m = 1,nq
!!        lf = fpidx(kf)
!!        fval = zero
!!        do l = 1,nr
!!          val = zero
!!          do i = interp(this_geom,this_order)% &
!!                           toFace(face_order)%row_ptr(lf)+1, &
!!                 interp(this_geom,this_order)% &
!!                           toFace(face_order)%row_ptr(lf+1)
!!            !
!!            n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
!!            a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
!!            !
!!            val = val + a * xyflx(n,m,l)
!!            !
!!          end do
!!          fval = fval + val*facenrm(l,f0+kf)
!!        end do
!!        do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
!!               correct(this_geom,this_order)%row_ptr(lf+1)
!!          n = n1 - 1 + correct(this_geom,this_order)%row_idx(i)
!!          correct_n = correct(this_geom,this_order)%row_val(i) * fval * ds
!!          residual(n,m) = residual(n,m) + correct_n
!!        end do
!!      end do
!!      !
!!    end do
!!    !
!!  end do
!!  !
!!end do
!!!
!!call debug_timer(leaving_procedure,pname)
!!!
!!end subroutine flux_gradient_at_solpts_visc
!
!###############################################################################
!
!!subroutine flux_gradient_at_solpts_visc_gpu(time_step,rk_step)
!!!
!!!... this subroutine get derivative of fluxes (no interaction)
!!!... for each cell, get f, g, ftd, gtd at the solution points
!!!
!!!.. Use Statements ..
!!use order_mod,         only : maxSP,maxFP,geom_solpts
!!use geovar,            only : nr,ncell
!!use geovar,            only : cell,face,fp,n_solpts
!!use interpolation_mod, only : interp
!!use derivatives_mod,   only : deriv
!!use metrics_mod,       only : metrics_array,jac,facenrm
!!use flowvar,           only : usp,gpsp,dusp,residual
!!use ovar,              only : governing_equations
!!use correction_mod,    only : correct
!!!
!!use visc_flux_mod, only : visc_flx_cv_gpu
!!use invs_flux_mod, only : invs_flx_cv_gpu
!!!
!!!.. Formal Arguments ..
!!integer, intent(in) :: time_step
!!integer, intent(in) :: rk_step
!!!
!!!.. Local Scalars ..
!!integer :: i,j,k,l,m,n,nc,n1,n2,np,f0
!!integer :: nf,sf,kf,lf
!!integer :: this_geom,this_order
!!integer :: face_geom,face_order
!!integer :: ierr,io
!!integer :: fp_offset,face_rot,nfp,ncf,num_cell_faces
!!real(wp) :: a,val,flxvel
!!real(wp) :: mvsc,mcdt,ave_ext_rate
!!real(wp) :: cv_inv,kflx
!!real(wp) :: fval,correct_n,ds
!!!
!!character(len=100) :: flux_fmt
!!character(len=300) :: spf_fname
!!!
!!!.. Local Arrays ..
!!real(wp), dimension(1:nr)      :: unrm
!!real(wp), dimension(1:nq,1:maxSP) :: lcl_res
!!real(wp), dimension(1:maxSP,1:nq,1:nr) :: vflx
!!real(wp), dimension(1:maxSP,1:nq,1:nr) :: rsflx
!!real(wp), dimension(1:maxSP,1:nq,1:nr) :: xyflx
!!integer, dimension(1:maxFP) :: fpidx
!!real(wp), dimension(1:nq) :: flx
!!real(wp), dimension(1:maxSP) :: pres
!!real(wp), dimension(1:maxSP,1:nq,1:nr) :: dpvdx
!!real(wp), dimension(1:nq) :: var
!!real(wp), dimension(1:nq,1:nr) :: vflux
!!real(wp), dimension(1:nq,1:nr) :: iflux
!!real(wp), dimension(1:maxSP,1:nr,1:nq) :: tflx
!!real(wp), dimension(1:maxSP,1:nr,1:nq) :: rflx
!!real(wp), dimension(1:maxFP) :: area
!!!
!!!.. Local Allocatable Arrays ..
!!real(wp), allocatable :: tfluxes(:,:,:)
!!real(wp), allocatable :: ifluxes(:,:,:)
!!real(wp), allocatable :: vfluxes(:,:,:)
!!!
!!!.. Local Parameters ..
!!character(len=*), parameter :: pname = "flux_gradient_at_solpts_visc_gpu"
!!!
!!continue
!!!
!!call debug_timer(entering_procedure,pname)
!!!
!!allocate ( tfluxes(1:n_solpts,1:nr,1:nq) )
!!allocate ( ifluxes(1:n_solpts,1:nr,1:nq) )
!!allocate ( vfluxes(1:n_solpts,1:nr,1:nq) )
!!!$ACC DATA CREATE(tfluxes,ifluxes,vfluxes)
!!!
!!call invs_flx_cv_gpu( usp , gpsp , ifluxes )
!!!
!!if (governing_equations == NavierStokes_Eqns) then
!!  !
!!  call visc_flx_cv_gpu( 1 , n_solpts , usp , gpsp , dusp , vfluxes )
!!  !
!!else
!!  !
!!  !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(vfluxes)
!!  do m = 1,nq
!!    do l = 1,nr
!!      do n = 1,n_solpts
!!        vfluxes(n,l,m) = zero
!!      end do
!!    end do
!!  end do
!!  !
!!end if
!!!
!!!$ACC PARALLEL LOOP COLLAPSE(3) &
!!!$ACC&              PRESENT(tfluxes,ifluxes,vfluxes)
!!do m = 1,nq
!!  do l = 1,nr
!!    do n = 1,n_solpts
!!      tfluxes(n,l,m) = ifluxes(n,l,m) - vfluxes(n,l,m)
!!    end do
!!  end do
!!end do
!!!
!!!#############################################################################
!!!
!!!$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
!!!$ACC&                   PRESENT(tfluxes) &
!!!$ACC&                   PRESENT(deriv) &
!!!$ACC&                   PRESENT(correct) &
!!!$ACC&                   PRESENT(fp) &
!!!$ACC&                   PRESENT(face) &
!!!$ACC&                   PRESENT(cell) &
!!!$ACC&                   PRESENT(interp) &
!!!$ACC&                   PRESENT(residual) &
!!!$ACC&                   PRESENT(metrics_array) &
!!!$ACC&                   PRESENT(jac) &
!!!$ACC&                   PRESENT(facenrm) &
!!!$ACC&                   PRESENT(cl) &
!!!$ACC& PRIVATE(rflx,lcl_res,tflx,fpidx,area)
!!do nc = 1,ncell
!!  !
!!  n1 = cell(nc)%beg_sp
!!  n2 = cell(nc)%end_sp
!!  np = n2-n1+1
!!  !
!!  this_geom = cell(nc)%geom
!!  this_order = cell(nc)%order
!!  !
!!  !$ACC LOOP VECTOR COLLAPSE(3)
!!  do k = 1,np
!!    do l = 1,nr
!!      do m = 1,nq
!!        n = k+n1-1
!!        tflx(k,l,m) = tfluxes(n,l,m)
!!      end do
!!    end do
!!  end do
!!  !
!!  !$ACC LOOP VECTOR COLLAPSE(3)
!!  do m = 1,nq
!!    do l = 1,nr
!!      do k = 1,np
!!        val = zero
!!        do i = 1,nr
!!          val = val + tflx(k,i,m)*metrics_array(n1-1+k,l,i)
!!        end do
!!        rflx(k,l,m) = val * jac(n1-1+k)
!!      end do
!!    end do
!!  end do
!!  !
!!  ! Compute the divergence of the discontinuous fluxes
!!  !
!!  ! Loop through each coordinate direction and add the divergence of
!!  ! the flux to the local residual array for the current cell
!!  !
!!  !$ACC LOOP VECTOR COLLAPSE(2)
!!  do k = 1,np
!!    do m = 1,nq
!!      val = zero
!       !$ACC LOOP SEQ
!!      do l = 1,nr
!!        do i = deriv(this_geom,this_order)%d(l)%row_ptr(k)+1, &
!!               deriv(this_geom,this_order)%d(l)%row_ptr(k+1)
!!          !
!!          n = deriv(this_geom,this_order)%d(l)%row_idx(i)
!!          a = deriv(this_geom,this_order)%d(l)%row_val(i)
!!          !
!!          val = val + a * rflx(n,l,m)
!!          !
!!        end do
!!      end do
!!      lcl_res(m,k) = val
!!    end do
!!  end do
!!  !
!!  ! Copy the local residual for the current cell into the residual array
!!  !
!!  !$ACC LOOP VECTOR COLLAPSE(2)
!!  do n = n1,n2
!!    do m = 1,nq
!!      residual(n,m) = lcl_res(m,n-n1+1)
!!    end do
!!  end do
!!  !
 !end do
 !!
 !!#############################################################################
 !!
 !!$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) PRIVATE(tflx,fpidx,area)
 !do nc = 1,ncell
 !  !
 !  n1 = cell(nc)%beg_sp
 !  n2 = cell(nc)%end_sp
 !  np = n2-n1+1
 !  !
 !  this_geom = cell(nc)%geom
 !  this_order = cell(nc)%order
!!  !
!!  ! Interpolate the fluxes to the interface flux points and dot this
!!  ! interpolated vector with the face normal at each flux point
!!  ! NOTE: We dont need to exchange this information for communication and
!!  !       periodic boundary faces since this flux is only used within the
!!  !       current cell.
!!  !
!!  ! Loop over the faces of this cell
!!  num_cell_faces = size(cell(nc)%face)
!!  !
!!  do ncf = 1,num_cell_faces
!!    !
!!    nf = cell(nc)%face(ncf)%idx  ! index of the current face
!!    sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
!!    ds = cl(this_geom) * merge(-one,one,sf==1)
!!    f0 = face(nf)%beg_fp-1
!!    !
!!    face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
!!    face_order = cell(nc)%face(ncf)%order ! order of current face
!!    fp_offset = cell(nc)%face(ncf)%fp_offset
!!    face_rot = cell(nc)%face(ncf)%rotation
!!    !
!!    nfp = geom_solpts(face_geom,face_order)
!!    fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
!!                                  rot(face_rot)%pts(1:nfp)
!!    !
!!    !$ACC LOOP VECTOR COLLAPSE(2)
!!    do kf = 1,nfp
!!      do m = 1,nq
!!        fval = zero
!!        lf = fpidx(kf)
!!        do l = 1,nr
!!          val = zero
!!          do i = interp(this_geom,this_order)% &
!!                           toFace(face_order)%row_ptr(lf)+1, &
!!                 interp(this_geom,this_order)% &
!!                           toFace(face_order)%row_ptr(lf+1)
!!            !
!!            n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
!!            a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
!!            !
!!            val = val + a * tflx(n,l,m)
!!            !
!!          end do
!!          fval = fval + val * facenrm(l,f0+kf)
!!        end do
!!        do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
!!               correct(this_geom,this_order)%row_ptr(lf+1)
!!          n = n1 - 1 + correct(this_geom,this_order)%row_idx(i)
!!          correct_n = correct(this_geom,this_order)%row_val(i) * fval * ds
!!          residual(n,m) = residual(n,m) + correct_n
!!        end do
!!      end do
!!    end do
!!    !
!!  end do
!!  !
!!end do
!!!
!!!$ACC END DATA
!!deallocate ( vfluxes )
!!deallocate ( ifluxes )
!!deallocate ( tfluxes )
!!!
!!call debug_timer(leaving_procedure,pname)
!!!
!!end subroutine flux_gradient_at_solpts_visc_gpu
!
!###############################################################################
!
!!subroutine correct_interface_flux_visc(time_step,rk_step)
!!!
!!use order_mod,         only : maxFP,geom_solpts
!!use geovar,            only : nr,nfbnd,nface
!!use geovar,            only : bface,face,cell,fp
!!use ovar,              only : invs_flux_method,governing_equations
!!use ovar,              only : output_dir,dump_fluxes
!!use metrics_mod,       only : facenrm
!!use correction_mod,    only : correct
!!use flowvar,           only : residual
!!use flowvar,           only : faceusp_array,facedusp_array
!!use flowvar,           only : facegpsp
!!!
!!use visc_flux_mod, only : common_visc_flx_cv
!!use visc_flux_mod, only : visc_wall_flx_cv
!!use invs_flux_mod, only : invs_wall_flx_cv
!!!
!!use upwind_flux_mod, only : upwind_flux_cv
!!use upwind_flux_mod, only : hllc_upwind_flux_cv
!!use upwind_flux_mod, only : ldfss_upwind_flux_cv
!!use upwind_flux_mod, only : rotated_rhll_upwind_flux_cv_2d
!!use upwind_flux_mod, only : rotated_rhll_upwind_flux_cv_3d
!!!
!!!.. Formal Arguments ..
!!integer, intent(in) :: time_step
!!integer, intent(in) :: rk_step
!!!
!!!.. Local Scalars ..
!!integer  :: i,k,l,m,n,nf,r1,r2,l1,l2
!!integer  :: kl,kr,iu,nfp,f0,f1,f2
!!integer  :: ierr,iofp
!!integer  :: face_geom,face_order
!!integer  :: host_cell,host_geom,host_order
!!integer  :: left_cell,left_geom,left_order
!!integer  :: right_cell,right_geom,right_order
!!real(wp) :: ds,dsl,dsr,cmat_val
!!!
!!character(len=300) :: fpf_fname
!!!
!!!.. Local Arrays ..
!!real(wp), dimension(1:nr) :: unrm
!!real(wp), dimension(1:nq) :: iflx,vflx
!!real(wp), dimension(1:nr,1:maxFP) :: rnrm
!!real(wp), dimension(1:nq,1:maxFP) :: lflx,rflx
!!real(wp), dimension(1:nq,1:maxFP) :: luqp,ruqp,fncq
!!real(wp), dimension(1:ng,1:maxFP) :: lgpqp,rgpqp
!!real(wp), dimension(1:nr,1:nq,1:maxFP) :: lduqp,rduqp
!!!
!!!.. Local Parameters ..
!!character(len=*), parameter :: pname = "correct_interface_flux_visc"
!!!
!!continue
!!!
!!call debug_timer(entering_procedure,pname)
!!!
!!! ##############################################################
!!! ##############################################################
!!! ##############################################################
!!! ##########              BOUNDARY FACES              ##########
!!! ##############################################################
!!! ##############################################################
!!! ##############################################################
!!!
!!! Loop through the boundary faces to apply the correction
!!! due to the boundary fluxes.
!!!
!!boundary_fluxes: do nf = 1,nfbnd
!!  !
!!  face_geom  = face(nf)%geom  ! geometry of face
!!  face_order = face(nf)%order ! order of face
!!  !
!!  f1 = face(nf)%beg_fp
!!  f2 = face(nf)%end_fp
!!  f0 = f1 - 1
!!  nfp = f2 - f0
!!  !
!!  host_cell  = face(nf)%left%cell    ! host cell on boundary face
!!  host_geom  = cell(host_cell)%geom  ! geometry of host cell
!!  host_order = cell(host_cell)%order ! order of host cell
!!  !
!!  l1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in left cell
!!  l2 = cell(host_cell)%end_sp ! ending    index for sol-pts in left cell
!!  !
!!  do k = 1,nfp
!!    do m = 1,nq
!!      ! FP solution from left  cell on face
!!      luqp(m,k) = faceusp_array(f0+k,m,1)
!!      ! FP solution from right cell on face
!!      ruqp(m,k) = faceusp_array(f0+k,m,2)
!!    end do
!!  end do
!!  !
!!  do k = 1,nfp
!!    do m = 1,ng
!!      ! FP gas properties from left  cell on face
!!      lgpqp(m,k) = facegpsp(f0+k,m,1)
!!      ! FP gas properties from right cell on face
!!      rgpqp(m,k) = facegpsp(f0+k,m,2)
!!    end do
!!  end do
!!  !
!!  if (governing_equations == NavierStokes_Eqns) then
!!    do k = 1,nfp
!!      do m = 1,nq
!!        do l = 1,nr
!!          ! Solution gradient at FP from left  cell on face
!!          lduqp(l,m,k) = facedusp_array(f0+k,l,m,1)
!!          ! Solution gradient at FP from right cell on face
!!          rduqp(l,m,k) = facedusp_array(f0+k,l,m,2)
!!        end do
!!      end do
!!    end do
!!  end if
!!  !
!!  rnrm(1:nr,1:nfp) = facenrm(1:nr,f1:f2) ! Area face normals at FP
!!  !
!!  ! Get the common flux at each flux point
!!  !
!!  do k = 1,nfp
!!    !
!!    ! Get the unit normal at this flux point
!!    !
!!    unrm(:) = unit_vector( rnrm(:,k) )
!!    !
!!    ! Get the inviscid upwind flux at the current flux point
!!    !
!!    if (any(bface(1,nf) == bc_walls)) then
!!      !
!!      ! This is a wall boundary so ensure that the inviscid
!!      ! boundary flux is only due to the pressure term.
!!      !
!!      iflx(:) = invs_wall_flx_cv( unrm(:) , luqp(:,k) )
!!      !
!!      ! Compute the viscous wall boundary flux
!!      !
!!      if (governing_equations == NavierStokes_Eqns) then
!!        vflx(:) = visc_wall_flx_cv( unrm(:) , ruqp(:,k) , &
!!                                    rgpqp(:,k) , rduqp(:,:,k) )
!!      end if
!!      !
!!    else
!!      !
!!      ! For all other boundaries, compute the boundary
!!      ! fluxes in the same way as the interior faces.
!!      !
!!      if (invs_flux_method == invs_flux_ROE) then
!!        iflx(:) = upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
!!      else if (invs_flux_method == invs_flux_HLLC) then
!!        iflx(:) = hllc_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
!!      else if (invs_flux_method == invs_flux_LDFSS) then
!!        iflx(:) = ldfss_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
!!      else if (invs_flux_method == invs_flux_Rotated_RHLL) then
!!        if (nr == 2) then
!!          iflx(:) = rotated_rhll_upwind_flux_cv_2d ( unrm(:) , luqp(:,k) , &
!!                                                               ruqp(:,k) )
!!        else
!!          iflx(:) = rotated_rhll_upwind_flux_cv_3d ( unrm(:) , luqp(:,k) , &
!!                                                               ruqp(:,k) )
!!        end if
!!      end if
!!      !
!!      ! Get the common viscous flux at the current flux point
!!      !
!!      if (governing_equations == NavierStokes_Eqns) then
!!        vflx(:) = common_visc_flx_cv( unrm(:), &
!!                                      luqp(:,k), lgpqp(:,k), lduqp(:,:,k), &
!!                                      ruqp(:,k), rgpqp(:,k), rduqp(:,:,k) )
!!      end if
!!      !
!!    end if
!!    !
!!    ! Get the total common boundary flux at this flux point
!!    !
!!    fncq(:,k) = iflx(:) - vflx(:)
!!    !
!!  end do
!!  !
!!  ! Apply the flux correction due to the flux jump
!!  ! NOTE: We only need to apply the correction to the left cell because
!!  !       this is a boundary face and the right cell is a ghost cell.
!!  !
!!  do k = 1,nfp
!!    !
!!    ! Get the index for the correction function
!!    ! corresponding to the current interface flux point
!!    !
!!    kl = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
!!                                      rot(face(nf)%left%rotation)% &
!!                                   pts(k)
!!    !
!!    ! Get the contribution to the face length/area from the current flux
!!    ! point for the left cell
!!    !
!!    ds  = norm2( rnrm(:,k) )
!!    dsl = ds * cl(host_geom)
!!    !
!!    ! Add the flux correction to the residuals of the sol-pts in cell
!!    ! host_cell
!!    !
!!    do i = correct(host_geom,host_order)%row_ptr(kl)+1, &
!!           correct(host_geom,host_order)%row_ptr(kl+1)
!!      n = l1 - 1 + correct(host_geom,host_order)%row_idx(i)
!!      cmat_val = dsl * correct(host_geom,host_order)%row_val(i)
!!      do m = 1,nq
!!        residual(n,m) = residual(n,m) + fncq(m,k) * cmat_val
!!      end do
!!    end do
!!    !
!!  end do
!!  !
!!end do boundary_fluxes
!!!
!!! ##############################################################
!!! ##############################################################
!!! ##############################################################
!!! ##########              INTERIOR FACES              ##########
!!! ##############################################################
!!! ##############################################################
!!! ##############################################################
!!!
!!! Loop over the interior faces and apply the correction from
!!! the fluxes at each face point to all the cell solution points
!!!
!!interior_fluxes: do nf = nfbnd+1,nface
!!  !
!!  face_geom = face(nf)%geom  ! geometry of face
!!  face_order = face(nf)%order ! order of face
!!  !
!!  f1 = face(nf)%beg_fp
!!  f2 = face(nf)%end_fp
!!  f0 = f1 - 1
!!  nfp = f2 - f0 ! number of flux points on this face
!!  !
!!  left_cell = face(nf)%left%cell ! left  cell on face
!!  right_cell = face(nf)%right%cell ! right cell on face
!!  !
!!  left_geom = cell(left_cell)%geom
!!  left_order = cell(left_cell)%order
!!  !
!!  right_geom = cell(right_cell)%geom
!!  right_order = cell(right_cell)%order
!!  !
!!  l1 = cell(left_cell)%beg_sp ! beginning index for sol-pts in left cell
!!  l2 = cell(left_cell)%end_sp ! ending    index for sol-pts in left cell
!!  !
!!  r1 = cell(right_cell)%beg_sp ! beginning index for sol-pts in right cell
!!  r2 = cell(right_cell)%end_sp ! ending    index for sol-pts in right cell
!!  !
!!  do k = 1,nfp
!!    do m = 1,nq
!!      ! FP solution from left  cell on face
!!      luqp(m,k) = faceusp_array(f0+k,m,1)
!!      ! FP solution from right cell on face
!!      ruqp(m,k) = faceusp_array(f0+k,m,2)
!!    end do
!!  end do
!!  !
!!  if (governing_equations == NavierStokes_Eqns) then
!!    do k = 1,nfp
!!      do m = 1,nq
!!        do l = 1,nr
!!          ! Solution gradient at FP from left  cell on face
!!          lduqp(l,m,k) = facedusp_array(f0+k,l,m,1)
!!          ! Solution gradient at FP from right cell on face
!!          rduqp(l,m,k) = facedusp_array(f0+k,l,m,2)
!!        end do
!!      end do
!!    end do
!!  end if
!!  !
!!  rnrm(1:nr,1:nfp) = facenrm(1:nr,f1:f2)
!!  !
!!  ! Get the common flux at each flux point
!!  !
!!  do k = 1,nfp
!!    !
!!    ! Get the unit normal at this flux point
!!    !
!!    unrm(:) = unit_vector( rnrm(:,k) )
!!    !
!!    ! Get the inviscid upwind flux at the current flux point
!!    !
!!    if (invs_flux_method == invs_flux_ROE) then
!!      iflx(:) = upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
!!    else if (invs_flux_method == invs_flux_HLLC) then
!!      iflx(:) = hllc_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
!!    else if (invs_flux_method == invs_flux_LDFSS) then
!!      iflx(:) = ldfss_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
!!    else if (invs_flux_method == invs_flux_Rotated_RHLL) then
!!      if (nr == 2) then
!!        iflx(:) = rotated_rhll_upwind_flux_cv_2d ( unrm(:) , luqp(:,k) , &
!!                                                             ruqp(:,k) )
!!      else
!!        iflx(:) = rotated_rhll_upwind_flux_cv_3d ( unrm(:) , luqp(:,k) , &
!!                                                             ruqp(:,k) )
!!      end if
!!    end if
!!    !
!!    ! Get the common viscous flux at the current flux point
!!    !
!!    if (governing_equations == NavierStokes_Eqns) then
!!      vflx(:) = common_visc_flx_cv( unrm(:), &
!!                                    luqp(:,k), lgpqp(:,k), lduqp(:,:,k), &
!!                                    ruqp(:,k), rgpqp(:,k), rduqp(:,:,k) )
!!    end if
!!    !
!!    ! Get the total common boundary flux at this flux point
!!    !
!!    fncq(:,k) = iflx(:) - vflx(:)
!!    !
!!  end do
!!  !
!!  ! Apply the flux correction due to the flux jump
!!  !
!!  do k = 1,nfp
!!    !
!!    ! Get the index for the correction function
!!    ! corresponding to the current interface flux point
!!    !
!!    kl = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
!!                                      rot(face(nf)%left%rotation)% &
!!                                   pts(k)
!!    kr = face(nf)%right%fp_offset + fp%geord(face_geom,face_order)% &
!!                                      rot(face(nf)%right%rotation)% &
!!                                    pts(k)
!!    !
!!    ! Get the face length for the left and right cell
!!    !
!!    ds  = norm2( rnrm(:,k) )
!!    dsl = ds * cl(left_geom)
!!    dsr = ds * cl(right_geom)
!!    !
!!    ! Add the flux correction to the residuals of the sol-pts in cell
!!    ! left_cell
!!    !
!!    do i = correct(left_geom,left_order)%row_ptr(kl)+1, &
!!           correct(left_geom,left_order)%row_ptr(kl+1)
!!      n = l1 - 1 + correct(left_geom,left_order)%row_idx(i)
!!      cmat_val = dsl * correct(left_geom,left_order)%row_val(i)
!!      do m = 1,nq
!!        residual(n,m) = residual(n,m) + fncq(m,k) * cmat_val
!!      end do
!!    end do
!!    !
!!    ! Add the flux correction to the residuals of the sol-pts in cell
!!    ! right_cell
!!    ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
!!    !
!!    do i = correct(right_geom,right_order)%row_ptr(kr)+1, &
!!           correct(right_geom,right_order)%row_ptr(kr+1)
!!      n = r1 - 1 + correct(right_geom,right_order)%row_idx(i)
!!      cmat_val = dsr * correct(right_geom,right_order)%row_val(i)
!!      do m = 1,nq
!!        residual(n,m) = residual(n,m) - fncq(m,k) * cmat_val
!!      end do
!!    end do
!!    !
!!  end do
!!  !
!!end do interior_fluxes
!!!
!!call debug_timer(leaving_procedure,pname)
!!!
!!end subroutine correct_interface_flux_visc
!
!###############################################################################
!
!!subroutine correct_gradients_br2()
!!!
!!use generic_types_mod, only : matrix
!!use order_mod,         only : maxSP,maxFP,geom_solpts
!!use geovar,            only : cell
!!use geovar,            only : nr,ncell,nfbnd,bface
!!use geovar,            only : face,fp
!!use metrics_mod,       only : facenrm,metrics_dt
!!use correction_mod,    only : correct
!!use flowvar,           only : dusp,facedusp_array,faceusp_array
!!use ovar,              only : walls_are_exact,LDG_beta
!!use interpolation_mod, only : interp
!!!
!!!.. Local Scalars ..
!!integer :: n,nc,n1,n2,np,k,l,m,f0
!!integer :: nf,kf,lf,sf,ncf
!!integer :: this_geom,this_order
!!integer :: face_geom,face_order
!!real(wp) :: dsl,bias,hfpb,correct_n,correct_m
!!real(wp) :: a,val
!!!
!!!.. Local Arrays ..
!!real(wp), dimension(1:nr) :: unrm
!!real(wp), dimension(1:nq) :: lufp,rufp,ujump
!!real(wp), dimension(1:maxSP,1:nr,1:nq) :: cell_dusp
!!real(wp), dimension(1:maxSP,1:nr,1:nq) :: face_dusp
!!real(wp), dimension(1:maxSP,1:nr,1:nq) :: dusp_copy
!!!
!!integer :: i,i1,i2
!!type(matrix) :: cmat
!!!
!!!.. Local Parameters ..
!!character(len=*), parameter :: pname = "correct_gradients_br2"
!!!
!!continue
!!!
!!call debug_timer(entering_procedure,pname)
!!!
!!! Loop through each interior cell and interpolate the conserved
!!! variables at the solution points to the flux points of each cell
!!!
!!#ifdef ENABLE_OLD_BR2
!!cell_loop: do nc = 1,ncell
!!  !
!!  n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
!!  n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
!!  np = n2-n1+1         ! number of sol-pts in this cell
!!  !
!!  this_geom = cell(nc)%geom   ! geometry type of this cell
!!  this_order = cell(nc)%order ! solution order of this cell
!!  !
!!  ! Store a local copy of the correction matrices divided by the value of
!!  ! the jacobian for each solution point. This will help reduce cache
!!  ! misses and prevent repetitive instances of division by the jacobian.
!!  !
!!  cmat = correct(this_geom,this_order)
!!  !
!!  do i = 1,size(cmat%row_ptr)-1
!!    i1 = cmat%row_ptr(i)+1
!!    i2 = cmat%row_ptr(i+1)
!!    cmat%row_val(i1:i2) = cmat%row_val(i1:i2) / &
!!                          metrics_dt(nc)%jac( cmat%row_idx(i1:i2) )
!!  end do
!!  !
!!  ! Create a local copy of dusp for this cell to which all the
!!  ! corrections will be applied. This will help reduce cache misses.
!!  !
!!  cell_dusp(1:np,1:nr,1:nq) = dusp(n1:n2,1:nr,1:nq)
!! !cell_dusp(1:nr,1:nq,1:np) = dusp(1:nr,1:nq,n1:n2)
!!  !
!!  ! Create a copy of the uncorrected dusp for this cell which
!!  ! will be used to reset the face_dusp array in the cell_faces loop. This
!!  ! will help reduce cache misses.
!!  !
!!  dusp_copy(1:np,1:nr,1:nq) = cell_dusp(1:np,1:nr,1:nq)
!! !do n = 1,np
!! !  dusp_copy(n,1:nr,1:nq) = cell_dusp(1:nr,1:nq,n)
!! !end do
!!  !
!!  ! Now Interpolate the interior solution gradients to the flux points
!!  !
!!  ! Loop over the faces of this cell
!!  !
!!  cell_faces: do ncf = 1,size(cell(nc)%face)
!!    !
!!    nf = cell(nc)%face(ncf)%idx  ! index of the current face
!!    sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
!!    !
!!    face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
!!    face_order = cell(nc)%face(ncf)%order ! order of current face
!!    !
!!    f0 = face(nf)%beg_sp - 1
!!    !
!!    ! Reset face_dusp to the original discontinuous gradient
!!    !
!!    face_dusp(1:np,1:nr,1:nq) = dusp_copy(1:np,1:nr,1:nq)
!!    !
!!    ! Loop over the points on this face
!!    !
!!    face_points: do kf = 1,geom_solpts(face_geom,face_order)
!!      !
!!      ! kf: index of flx-pt local to this face
!!      ! lf: index of flx-pt within interp%mat
!!      lf = cell(nc)%face(ncf)%fp_offset + fp%geord(face_geom,face_order)% &
!!                                        rot(cell(nc)%face(ncf)%rotation)% &
!!                                          pts(kf)
!!      !
!!      ! Solution at this flux point for the left and right cells
!!      !
!!      do m = 1,nq
!!        lufp(m) = faceusp_array(f0+kf,m,1)
!!        rufp(m) = faceusp_array(f0+kf,m,2)
!!      end do
!!      !
!!      ! Get the unit normal at this flux point
!!      !
!!      unrm(1:nr) = unit_vector( facenrm(1:nr,f0+kf) )
!!      !
!!      ! Get the face length for the current face
!!      !
!!      dsl = cl(this_geom) * norm2( facenrm(1:nr,f0+kf) )
!!      !
!!      ! Compute the BR2 jump correction for this point
!!      !    jump_left  =  half * (u_right - u_left)
!!      !    jump_right = -half * (u_right - u_left)
!!      !
!!      ! NOTE: Multiplying by the normal in the loop below that applies the
!!      !       correction
!!      !           i.e., correct_m * unrm(l)
!!      !       takes care of the minus sign if the current cell is the right
!!      !       cell on the current face.
!!      !
!!      ! NOTE: Boundary faces should never be the right cell on the interface.
!!      !       For non-wall boundaries, this should reduce to jump_left above.
!!      !       For wall boundaries, jump_left should reduce to:
!!      !         jump_left = (u_right - u_left) where u_right = u_wall
!!      !
!!      bias = zero
!!      !
!!      ! Reset bias if this is a wall boundary
!!      !
!!      if (nf <= nfbnd .and. walls_are_exact) then
!!        if (any(bface(1,nf) == bc_walls)) bias = half
!!      end if
!!      !
!!      hfpb = bias + half
!!      !
!!      ! Get the interface jumps to correct the gradients
!!      !
!!      do m = 1,nq
!!        ujump(m) = hfpb*(rufp(m) - lufp(m))
!!      end do
!!      !
!!      ! Apply the jump correction to the interior discontinuous gradients
!!      !
!!      do m = 1,nq
!!        correct_m = dsl * ujump(m)
!!        do l = 1,nr
!!          correct_n = correct_m * unrm(l)
!!          do i = cmat%row_ptr(lf)+1, cmat%row_ptr(lf+1)
!!            n = cmat%row_idx(i)
!!            a = cmat%row_val(i)
!!            face_dusp(n,l,m) = face_dusp(n,l,m) + a * correct_n
!!            cell_dusp(n,l,m) = cell_dusp(n,l,m) + a * correct_n
!!          end do
!!        end do
!!      end do
!!      !
!!    end do face_points
!!    !
!!    ! Interpolate the face corrected gradient from the interior solution
!!    ! points to the interface flux points now that we are done with all
!!    ! the flux points on the current face of the current cell.
!!    !
!!    interpolate_face: do kf = 1,geom_solpts(face_geom,face_order)
!!      !
!!      ! kf: index of flx-pt local to this face
!!      ! lf: index of flx-pt within interp%mat
!!      lf = cell(nc)%face(ncf)%fp_offset + fp%geord(face_geom,face_order)% &
!!                                        rot(cell(nc)%face(ncf)%rotation)% &
!!                                          pts(kf)
!!      !
!!     !facedusp_array(f0+kf,1:nr,1:nq,sf) = interp(this_geom,this_order)% &
!!     !                                            toFace(face_order)% &
!!     !                               dot_mat2( lf ,face_dusp(1:np,1:nr,1:nq) )
!!      do m = 1,nq
!!        do l = 1,nr
!!         !facedusp_array(f0+kf,l,m,sf) = interp(this_geom,this_order)% &
!!         !                                      toFace(face_order)% &
!!         !                            dot( lf ,face_dusp(1:np,l,m) )
!!          val = zero
!!          do i = interp(this_geom,this_order)% &
!!                           toFace(face_order)%row_ptr(lf)+1, &
!!                 interp(this_geom,this_order)% &
!!                           toFace(face_order)%row_ptr(lf+1)
!!            !
!!            n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
!!            a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
!!            !
!!            val = val + a * face_dusp(n,l,m)
!!            !
!!          end do
!!          facedusp_array(f0+kf,l,m,sf) = val
!!        end do
!!      end do
!!      !
!!    end do interpolate_face
!!    !
!!  end do cell_faces
!!  !
!!  ! Copy the total corrected gradient accounting for all cell faces
!!  ! back into dusp
!!  !
!!  dusp(n1:n2,1:nr,1:nq) = cell_dusp(1:np,1:nr,1:nq)
!!  !
!!end do cell_loop
!!#endif
!!!
!!call debug_timer(leaving_procedure,pname)
!!!
!!end subroutine correct_gradients_br2
!
!###############################################################################
!
!!subroutine correct_gradients_at_solpts()
!!!
!!!.. Use Statements ..
!!use order_mod,         only : maxpts,geom_solpts
!!use geovar,            only : nfbnd,nface,nr,n_solpts,ncell
!!use geovar,            only : face,cell,bface,fp
!!use metrics_mod,       only : facenrm,metrics_dt
!!use correction_mod,    only : correct
!!use flowvar,           only : dusp,faceusp_array
!!use ovar,              only : walls_are_exact
!!use ovar,              only : LDG_beta
!!!
!!!.. Local Scalars ..
!!integer  :: i,k,l,m,n,nk,nf,r1,r2,l1,l2,n1,n2,nc
!!integer  :: kl,kr,nfp,f0
!!integer  :: face_geom,face_order
!!integer  :: host_cell,host_geom,host_order
!!integer  :: left_cell,left_geom,left_order
!!integer  :: right_cell,right_geom,right_order
!!real(wp) :: bias,hfpb,hfmb,cmat_val
!!real(wp) :: correct_n,correct_m,ds,dsl,dsr
!!!
!!!.. Local Arrays ..
!!real(wp), dimension(1:nr)          :: unrm
!!real(wp), dimension(1:nq)          :: ucom
!!real(wp), dimension(1:nq)          :: ljfp,rjfp
!!real(wp), dimension(1:nq,1:maxpts) :: lufp,rufp
!!!
!!!.. Local Parameters ..
!!character(len=*), parameter :: pname = "correct_gradients_at_solpts"
!!!
!!continue
!!!
!!call debug_timer(entering_procedure,pname)
!!!
!!! #######################
!!! #######################
!!! ### BOUNDARY FACES  ###
!!! #######################
!!! #######################
!!!
!!! Loop through the boundary faces and correct the gradient of the
!!! continuous solution polynomial in the cells adjacent to the boundary.
!!!
!!boundary_faces: do nf = 1,nfbnd
!!  !
!!  face_geom = face(nf)%geom  ! geometry of face
!!  face_order = face(nf)%order ! order of face
!!  !
!!  f0 = face(nf)%beg_fp - 1
!!  nfp = face(nf)%end_fp - f0 ! number of flux points on this face
!!  !
!!  host_cell = face(nf)%left%cell ! left  cell on face
!!  !
!!  host_geom = cell(host_cell)%geom
!!  host_order = cell(host_cell)%order
!!  !
!!  l1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in left cell
!!  l2 = cell(host_cell)%end_sp ! ending    index for sol-pts in left cell
!!  !
!!  ! For all the flux points on the current face, extract the conservative
!!  ! conservative variables for the left and right cells from faceusp
!!  !
!!  do k = 1,nfp
!!    do m = 1,nq
!!      lufp(m,k) = faceusp_array(f0+k,m,1)
!!      rufp(m,k) = faceusp_array(f0+k,m,2)
!!    end do
!!  end do
!!  !
!!  ! Correct the gradients of the conservative variables using the
!!  ! solution jumps at the interface and the correction function
!!  !
!!  ! Use the LDG flux formulation for computing the common interface solutions
!!  !   u_com = half*(u_left + u_right) - &
!!  !                (u_left - u_right)*dot(LDG_beta,unit_normal)
!!  !
!!  do k = 1,nfp
!!    !
!!    ! Get the index for the correction function
!!    ! corresponding to the current interface flux point
!!    !
!!    kl = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
!!                                      rot(face(nf)%left%rotation)% &
!!                                   pts(k)
!!    !
!!    ! Get the unit normal at this flux point
!!    !
!!    unrm(1:nr) = unit_vector( facenrm(1:nr,f0+k) )
!!    !
!!    ! Get the face length for the current face
!!    !
!!    dsl = cl(host_geom) * norm2( facenrm(1:nr,f0+k) )
!!    !
!!    ! Compute the LDG directional parameters
!!    ! NOTE: bias is used for all boundary faces that allow for flow
!!    !       across the interface, e.g, inflow/outflow, cpu/partition
!!    !       boundary, or periodic boundary. These boundaries are
!!    !       treated the same as an interior face except that we dont
!!    !       need to worry about correcting the right cell.
!!    !
!!    bias = dot( LDG_beta(1:nr) , unrm(1:nr) )
!!    !
!!    ! Reset bias if this is a wall boundary
!!    !
!!    if (walls_are_exact) then
!!      if (any(bface(1,nf) == bc_walls)) bias = half
!!    end if
!!    !
!!    ! Use the merge function to overide bias with 1/2 if this is a
!!    ! wall boundary face meaning there is no flow across the interface,
!!    ! e.g., slip wall, symmetry boundary, no-slip wall. For these
!!    ! boundaries, the solution in rufp is the exact solution at the boundary
!!    ! and we use this exact solution for the interface common solution.
!!    ! This means the left jump is simply:
!!    !           jump_left = u_wall - u_interior
!!    !
!!    hfpb = half + bias
!!    !
!!    ! Compute the interface jumps for the left and right solutions
!!    ! NOTE: The common interface solution simplifies to
!!    !         u_com = (half - bias)*u_left + (half + bias)*u_right
!!    !       Using jump_left = u_com - u_left
!!    !         jump_left  = (half + bias)*(u_right-u_left)
!!    !
!!    do m = 1,nq
!!      ljfp(m) = hfpb*(rufp(m,k) - lufp(m,k))
!!    end do
!!    !
!!    ! Unsimplified version of the above LDG flux formulation
!!    !
!!   !if (any(bface(1,nf) == bc_walls)) then
!!   !  ucom(1:nq) = rufp(1:nq,k)
!!   !else
!!   !  do m = 1,nq
!!   !    ucom(m) = lufp(m,k)*(half - dot(LDG_beta(1:nr),+unrm(1:nr))) + &
!!   !              rufp(m,k)*(half - dot(LDG_beta(1:nr),-unrm(1:nr)))
!!   !  end do
!!   !end if
!!   !!
!!   !do m = 1,nq
!!   !  ljfp(m) = ucom(m) - lufp(m,k)
!!   !end do
!!    !
!!    ! Now use the correction function with the interface jump to correct
!!    ! the gradients of the conservative variables in the left cell
!!    !
!!    do i = correct(host_geom,host_order)%row_ptr(kl)+1, &
!!           correct(host_geom,host_order)%row_ptr(kl+1)
!!      nk = correct(host_geom,host_order)%row_idx(i)
!!      n = l1 - 1 + nk
!!      correct_n = dsl * correct(host_geom,host_order)%row_val(i) &
!!                      / metrics_dt(host_cell)%jac(nk)
!!      do m = 1,nq
!!        correct_m = correct_n * ljfp(m)
!!        do l = 1,nr
!!          dusp(n,l,m) = dusp(n,l,m) + correct_m * unrm(l)
!!        end do
!!      end do
!!    end do
!!    !
!!  end do
!!  !
!!end do boundary_faces
!!!
!!! #######################
!!! #######################
!!! ### INTERIOR FACES  ###
!!! #######################
!!! #######################
!!!
!!! Loop over the interior faces and correct the gradient of the
!!! continuous solution polynomial in the cells adjacent to the face
!!!
!!interior_faces: do nf = nfbnd+1,nface
!!  !
!!  face_geom = face(nf)%geom  ! geometry of face
!!  face_order = face(nf)%order ! order of face
!!  !
!!  f0 = face(nf)%beg_fp - 1
!!  nfp = face(nf)%end_fp - f0 ! number of flux points on this face
!!  !
!!  left_cell = face(nf)%left%cell ! left  cell on face
!!  right_cell = face(nf)%right%cell ! right cell on face
!!  !
!!  left_geom = cell(left_cell)%geom
!!  left_order = cell(left_cell)%order
!!  !
!!  right_geom = cell(right_cell)%geom
!!  right_order = cell(right_cell)%order
!!  !
!!  l1 = cell(left_cell)%beg_sp ! beginning index for sol-pts in left cell
!!  l2 = cell(left_cell)%end_sp ! ending    index for sol-pts in left cell
!!  !
!!  r1 = cell(right_cell)%beg_sp ! beginning index for sol-pts in right cell
!!  r2 = cell(right_cell)%end_sp ! ending    index for sol-pts in right cell
!!  !
!!  ! For all the flux points on the current face, extract the conservative
!!  ! conservative variables for the left and right cells from faceusp
!!  !
!!  do k = 1,nfp
!!    do m = 1,nq
!!      lufp(m,k) = faceusp_array(f0+k,m,1)
!!      rufp(m,k) = faceusp_array(f0+k,m,2)
!!    end do
!!  end do
!!  !
!!  ! Correct the gradients of the conservative variables using the
!!  ! solution jumps at the interface and the correction function
!!  !
!!  ! Use the LDG flux formulation for computing the common interface solutions
!!  !   u_com = half*(u_left + u_right) - &
!!  !                (u_left - u_right)*dot(LDG_beta,unit_normal)
!!  !
!!  do k = 1,nfp
!!    !
!!    ! Get the index for the correction function
!!    ! corresponding to the current interface flux point
!!    !
!!    kl = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
!!                                      rot(face(nf)%left%rotation)% &
!!                                   pts(k)
!!    kr = face(nf)%right%fp_offset + fp%geord(face_geom,face_order)% &
!!                                      rot(face(nf)%right%rotation)% &
!!                                   pts(k)
!!    !
!!    ! Get the unit normal at this flux point
!!    !
!!    unrm(1:nr) = unit_vector( facenrm(1:nr,f0+k) )
!!    !
!!    ! Get the face length for the current face
!!    !
!!    ds = norm2( facenrm(1:nr,f0+k) )
!!    dsl = cl(left_geom)*ds
!!    dsr = cl(right_geom)*ds
!!    !
!!    ! Compute the LDG directional parameters
!!    !
!!    bias = dot( LDG_beta(1:nr) , unrm(1:nr) )
!!    !
!!    hfpb = half + bias
!!    hfmb = half - bias
!!    !
!!    ! Compute the interface jumps for the left and right solutions
!!    ! NOTE: The common interface solution simplifies to
!!    !         u_com = (half - bias)*u_left + (half + bias)*u_right
!!    !       Using jump_left = u_com - u_left and jump_right = u_com - u_right
!!    !         jump_left  = (half + bias)*(u_right-u_left)
!!    !         jump_right = (half - bias)*(u_left-u_right)
!!    !
!!    do m = 1,nq
!!      ljfp(m) = hfpb*(rufp(m,k) - lufp(m,k))
!!      rjfp(m) = hfmb*(lufp(m,k) - rufp(m,k))
!!    end do
!!    !
!!    ! Unsimplified version of the above LDG flux formulation
!!    !
!!   !do m = 1,nq
!!   !  ucom(m) = lufp(m,k)*(half - dot(LDG_beta(1:nr),+unrm(1:nr))) + &
!!   !            rufp(m,k)*(half - dot(LDG_beta(1:nr),-unrm(1:nr)))
!!   !end do
!!   !!
!!   !do m = 1,nq
!!   !  ljfp(m) = ucom(m) - lufp(m,k)
!!   !  rjfp(m) = ucom(m) - rufp(m,k)
!!   !end do
!!    !
!!    ! Now use the correction function with the interface jump to correct
!!    ! the gradients of the conservative variables in the left cell
!!    !
!!    do i = correct(left_geom,left_order)%row_ptr(kl)+1, &
!!           correct(left_geom,left_order)%row_ptr(kl+1)
!!      nk = correct(left_geom,left_order)%row_idx(i)
!!      n = l1 - 1 + nk
!!      correct_n = dsl * correct(left_geom,left_order)%row_val(i) &
!!                      / metrics_dt(left_cell)%jac(nk)
!!      do m = 1,nq
!!        correct_m = correct_n * ljfp(m)
!!        do l = 1,nr
!!          dusp(n,l,m) = dusp(n,l,m) + correct_m * unrm(l)
!!        end do
!!      end do
!!    end do
!!    !
!!    ! Now use the correction function with the interface jump to correct
!!    ! the gradients of the conservative variables in the right cell
!!    !
!!    do i = correct(right_geom,right_order)%row_ptr(kr)+1, &
!!           correct(right_geom,right_order)%row_ptr(kr+1)
!!      nk = correct(right_geom,right_order)%row_idx(i)
!!      n = r1 - 1 + nk
!!      correct_n = dsr * correct(right_geom,right_order)%row_val(i) &
!!                      / metrics_dt(right_cell)%jac(nk)
!!      do m = 1,nq
!!        correct_m = correct_n * rjfp(m)
!!        do l = 1,nr
!!          dusp(n,l,m) = dusp(n,l,m) - correct_m * unrm(l)
!!        end do
!!      end do
!!    end do
!!    !
!!  end do
!!  !
!!end do interior_faces
!!!
!!call debug_timer(leaving_procedure,pname)
!!!
!!end subroutine correct_gradients_at_solpts
#endif
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! END OF OLD FLUX SUBROUTINES
!
!###############################################################################
!###############################################################################
!###############################################################################
!
#include "Functions/check_for_bad_number.f90"
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
!#include "Functions/inviscid_fluxes.f90"
!#include "Functions/viscous_fluxes.f90"
!
!###############################################################################
!
!#include "Functions/upwind_fluxes.f90"
!include "Functions/grad_cv_to_grad_pv.f90"
!include "Functions/viscosity.f90"
!include "Functions/pressure.f90"
!
!###############################################################################
!
end module flux_mod
