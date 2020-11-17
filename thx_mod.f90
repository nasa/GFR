module thx_mod
  !
!#define PROFILE_ON
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer => debug_timer
  use eqn_idx, only : nmx,nmy,nmz
  !
  implicit none
  !
  private
  !
  public :: compute_blowing_ratio
  !
  !.. Local Module Parameters ..
  !
  real(wp), parameter :: noz_height = 2.68_wp*0.0254_wp
  real(wp), parameter :: inj_radius = 0.25_wp*0.0254_wp
  real(wp), parameter :: lip_height = 0.05_wp*0.0254_wp + noz_height
  !
  real(wp), parameter :: noz_area = noz_height**2
  real(wp), parameter :: inj_area = PI*(inj_radius**2)
  !
  !.. Local Module Variables ..
  !
  integer, save :: br_is_not_initialized = true
  !
  !.. Local Module Arrays ..
  !
  integer, allocatable, save :: my_noz_faces(:)
  integer, allocatable, save :: my_inj_faces(:)
  !
contains
!
!###############################################################################
!
subroutine compute_blowing_ratio(blowing_ratio,mdot_nozzle,mdot_injector)
  !
  !.. Use Statements ..
  use ovar,           only : rhoref,aref
  use geovar,         only : face
  use flowvar,        only : faceusp_array
  use metrics_mod,    only : facenrm
  use quadrature_mod, only : std_elem
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: blowing_ratio
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(inout) :: mdot_nozzle
  real(wp), optional, intent(inout) :: mdot_injector
  !
  !.. Local Scalars ..
  integer  :: n,nf,sf,f0,nfp,k
  integer  :: face_geom,face_order
  real(wp) :: rsgn,rnx,rny,rnz,unx,uny,unz,ds
  real(wp) :: rhou,rhov,rhow,mdot
  real(wp) :: inj_mdot,noz_mdot
  real(wp) :: my_noz_mdot,my_inj_mdot
  !
  !.. Local Arrays ..
  real(wp) :: all_mdots(1:2,1:ncpu)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_blowing_ratio"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (br_is_not_initialized) then
    call init_blowing_ratio()
    br_is_not_initialized = fals
  end if
  !
  ! Compute the mass flow rate of the injector outflow surface
  !
  my_inj_mdot = zero
  !
  do n = 1,size(my_inj_faces)
    !
    nf = abs(my_inj_faces(n))
    sf = merge( 1 , 2 , my_inj_faces(n) > 0 )
    !
    rsgn = merge( one , -one , sf == 1 )
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    f0  = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do k = 1,nfp
      !
      rnx = facenrm(1,f0+k)
      rny = facenrm(2,f0+k)
      rnz = facenrm(3,f0+k)
     !ds  = sqrt(rnx*rnx+rny*rny+rnz*rnz) * rsgn
     !unx = rnx/ds
     !uny = rny/ds
     !unz = rnz/ds
      !
      rhou = faceusp_array(f0+k,nmx,sf)
      rhov = faceusp_array(f0+k,nmy,sf)
      rhow = faceusp_array(f0+k,nmz,sf)
      !
     !mdot = rhou*unx + rhov*uny + rhow*unz
      mdot = rhou*rnx + rhov*rny + rhow*rnz
      mdot = mdot*rsgn
      !
      my_inj_mdot = my_inj_mdot + mdot*std_elem(face_geom,face_order)%wts(k)
      !
    end do
    !
  end do
  !
  ! Compute the mass flow rate of the nozzle outflow surface
  !
  my_noz_mdot = zero
  !
  do n = 1,size(my_noz_faces)
    !
    nf = abs(my_noz_faces(n))
    sf = merge( 1 , 2 , my_noz_faces(n) > 0 )
    !
    rsgn = merge( one , -one , sf == 1 )
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    f0  = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do k = 1,nfp
      !
      rnx = facenrm(1,f0+k)
      rny = facenrm(2,f0+k)
      rnz = facenrm(3,f0+k)
     !ds  = sqrt(rnx*rnx+rny*rny+rnz*rnz) * rsgn
     !unx = rnx/ds
     !uny = rny/ds
     !unz = rnz/ds
      !
      rhou = faceusp_array(f0+k,nmx,sf)
      rhov = faceusp_array(f0+k,nmy,sf)
      rhow = faceusp_array(f0+k,nmz,sf)
      !
     !mdot = rhou*unx + rhov*uny + rhow*unz
      mdot = rhou*rnx + rhov*rny + rhow*rnz
      mdot = mdot*rsgn
      !
      my_noz_mdot = my_noz_mdot + mdot*std_elem(face_geom,face_order)%wts(k)
      !
    end do
    !
  end do
  !
  all_mdots(:,:) = zero
  !
  all_mdots(1,mypnum+1) = my_inj_mdot
  all_mdots(2,mypnum+1) = my_noz_mdot
  !
  if (ncpu > 1) then
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                       all_mdots,size(all_mdots,dim=1,kind=INT_MPI), &
                       mpi_flttyp,MPI_COMM_WORLD,mpierr)
  end if
  !
  inj_mdot = sum(all_mdots(1,:)) * rhoref*aref
  noz_mdot = sum(all_mdots(2,:)) * rhoref*aref
  !
  if (abs(noz_mdot) < eps6) then
    noz_mdot = noz_area
  end if
  !
  blowing_ratio = (inj_mdot/inj_area) / (noz_mdot/noz_area)
  !
  if (present(mdot_injector)) mdot_injector = inj_mdot*2.20462_wp
  if (present(mdot_nozzle  )) mdot_nozzle   = noz_mdot*2.20462_wp
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_blowing_ratio
!
!###############################################################################
!
subroutine init_blowing_ratio()
  !
  !.. Use Statements ..
  use geovar,         only : ncell,cell,nface,face,xyz
  use flowvar,        only : facexyz_array
  use metrics_mod,    only : facenrm
  use quadrature_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer  :: n,nc,ncf,nf,nfp,n1,n2,k,f0,f1,f2,sf,ierr
  integer  :: face_geom,face_order
  real(wp) :: my_inj_area,sum_inj_area,exp_inj_area
  real(wp) :: my_noz_area,sum_noz_area,exp_noz_area
  !
  !.. Local Arrays ..
  real(wp) :: all_areas(1:2,1:ncpu)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: inj_msk(:)
  integer, allocatable :: noz_msk(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_blowing_ratio"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  noz_msk = init_array(nface,0) ! F2003 AUTO-REALLOC
  inj_msk = init_array(nface,0) ! F2003 AUTO-REALLOC
  !
  cell_loop: &
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    ! Check to see if any of the faces of this cell
    ! are on the injector outflow surface.
    !
    inj_cell_loop: &
    do n = 1,1
      !
      if (any(xyz(2,n1:n2) > zero)) exit inj_cell_loop
      if (any(xyz(1,n1:n2) < zero)) exit inj_cell_loop
      if (any(xyz(2,n1:n2) < -0.002_wp)) exit inj_cell_loop
      !
      do ncf = 1,size(cell(nc)%face)
        !
        nf = cell(nc)%face(ncf)%idx
        sf = cell(nc)%face(ncf)%side
        !
        f1 = face(nf)%beg_fp
        f2 = face(nf)%end_fp
        !
        ! If all the y-coordinates of the points on this face are within 1.0e-6
        ! of zero, then 1) this face is part of the injector outflow surface so
        ! mark this face in the injector mask array, and 2) cycle to the next
        ! cell because this cell cant also be on the nozzle outflow surface.
        !
        if (all(abs(facexyz_array(2,f1:f2)) < eps6)) then
          inj_msk(nf) = merge(1,-1,sf==1)
          cycle cell_loop
        end if
        !
      end do
      !
    end do inj_cell_loop
    !
    ! Check to see if any of the faces of this cell
    ! are on the nozzle outflow surface.
    !
    noz_cell_loop: &
    do n = 1,1
      !
      if (any(xyz(1,n1:n2) > zero)) exit noz_cell_loop
      if (any(xyz(1,n1:n2) < -0.002_wp)) exit noz_cell_loop
      if (any(xyz(2,n1:n2) < zero)) exit noz_cell_loop
      if (any(xyz(2,n1:n2) > lip_height)) exit noz_cell_loop
      if (any(abs(xyz(3,n1:n2)) > half*lip_height)) exit noz_cell_loop
      !
      do ncf = 1,size(cell(nc)%face)
        !
        nf = cell(nc)%face(ncf)%idx
        sf = cell(nc)%face(ncf)%side
        !
        f1 = face(nf)%beg_fp
        f2 = face(nf)%end_fp
        !
        ! If all the x-coordinates of the points on this face are within
        ! 1.0e-6 of zero, then this face is part of the nozzle outflow
        ! surface so mark this face in the nozzle mask array.
        !
        if (all(abs(facexyz_array(1,f1:f2)) < eps6)) then
          noz_msk(nf) = merge(1,-1,sf==1)
          cycle cell_loop
        end if
        !
      end do
      !
    end do noz_cell_loop
    !
  end do cell_loop
  !
  ! Collect the face indices for the injector and nozzle surfaces
  ! and store them respectively in my_inj_faces and my_noz_faces.
  !
  my_inj_faces = pack( intseq(1,nface) , inj_msk /= 0 ) ! F2003 AUTO-REALLOC
  my_noz_faces = pack( intseq(1,nface) , noz_msk /= 0 ) ! F2003 AUTO-REALLOC
  !
  ! Make the face indices negative for any faces
  ! whose cell is the right cell on the face
  !
  my_inj_faces(:) = my_inj_faces(:) * inj_msk(my_inj_faces(:))
  my_noz_faces(:) = my_noz_faces(:) * noz_msk(my_noz_faces(:))
  !
  ! Deallocate inj_msk and noz_msk since we are done with them
  !
  if (allocated(inj_msk)) deallocate ( inj_msk , stat=ierr )
  if (allocated(noz_msk)) deallocate ( noz_msk , stat=ierr )
  !
  ! Compute the area of the injector face on this MPI process
  !
  my_inj_area = zero
  !
  do n = 1,size(my_inj_faces)
    !
    nf = abs(my_inj_faces(n))
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    f0  = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do k = 1,nfp
      my_inj_area = my_inj_area + norm2(facenrm(:,f0+k)) * &
                                  std_elem(face_geom,face_order)%wts(k)
    end do
    !
  end do
  !
  ! Compute the area of the nozzle face on this MPI process
  !
  my_noz_area = zero
  !
  do n = 1,size(my_noz_faces)
    !
    nf = abs(my_noz_faces(n))
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    f0  = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do k = 1,nfp
      my_noz_area = my_noz_area + norm2(facenrm(:,f0+k)) * &
                                  std_elem(face_geom,face_order)%wts(k)
    end do
    !
  end do
  !
  ! Collect the individual injector and nozzle
  ! surface areas from all MPI processes
  !
  all_areas(:,:) = zero
  !
  all_areas(1,mypnum+1) = my_inj_area
  all_areas(2,mypnum+1) = my_noz_area
  !
  if (ncpu > 1) then
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                       all_areas,size(all_areas,dim=1,kind=INT_MPI), &
                       mpi_flttyp,MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Get the total injector and nozzle surface areas
  !
  sum_inj_area = sum(all_areas(1,:))
  sum_noz_area = sum(all_areas(2,:))
  !
  if (mypnum == glb_root) then
    !
    sum_inj_area = sum_inj_area / (0.0254_wp**2)
    sum_noz_area = sum_noz_area / (0.0254_wp**2)
    exp_inj_area = PI*one4*half
    exp_noz_area = noz_area / (0.0254_wp**2)
    !
    write (iout,1) sum_inj_area,exp_inj_area, &
                   ten2*abs(sum_inj_area-exp_inj_area)/exp_inj_area, &
                   sum_noz_area,exp_noz_area, &
                   ten2*abs(sum_noz_area-exp_noz_area)/exp_noz_area
    !
    flush (iout)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/, &
            80("="),/, &
            8x,"Surface",7x,"Computed Area",3x,"Expected Area",3x,"% Diff",/, &
            80("-"),/, &
            2x,"Injector Outflow : ",es14.6,2x,es14.6,2x,f9.2,/, &
            2x,"  Nozzle Outflow : ",es14.6,2x,es14.6,2x,f9.2,/, &
            80("="),/)
  !
end subroutine init_blowing_ratio
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module thx_mod
