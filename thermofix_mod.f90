!###############################################################################
!> Module for finding and limiting cells with bad thermo properties
!> (pressure, temperature, density).
!###############################################################################
module thermofix_mod
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx, only : ncs,npres,ntemp
  !
  use ovar, only : force_nd_subroutines
  !
  implicit none
  !
  private
  !
  !
  !.. Public Procedures ..
  !
  public :: check_cell_thermo_properties
  public :: check_face_thermo_properties
  !
  !
  !.. Private Procedures ..
  !
  private :: reinterpolate_bad_cells
  private :: reinterpolate_bad_cells_3d
  private :: reinterpolate_bad_cells_2d
  private :: reinterpolate_bad_cells_1d
  private :: reinterpolate_bad_cells_nd
  private :: emergency_limit_bad_cells
  private :: emergency_limit_bad_cells_3d
  private :: emergency_limit_bad_cells_2d
  private :: emergency_limit_bad_cells_1d
  private :: emergency_limit_bad_cells_nd
  !
contains
!
!###############################################################################
!> Check volume points for bad thermo properties
!> (pressure, temperature, density) and limit the cell. Do not update faces.
!###############################################################################
subroutine check_cell_thermo_properties()
  !
  !.. Use Statements ..
  use geovar,  only : ncell,cell
  use flowvar, only : usp,gpsp
  use ovar,    only : rhoref,tref,aref
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_is_negative
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,n
  logical(lk) :: bad_cell_found
  real(wp) :: dens_limit,temp_limit,pres_limit
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable :: bad_cell(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_cell_thermo_properties"
  !
continue
  !
  if (ncs > 0) return
  !
  dens_limit = 0.1_wp/rhoref
  temp_limit = 100.0_wp/tref
  pres_limit = 40000.0_wp/(rhoref*aref*aref)
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Find any cells with solution points that have
  ! a negative density, temperature, or pressure
  !
  bad_cell_found = fals
  !
  bad_cell = init_array(ncell,fals) ! F2003 AUTO-REALLOC
  !$ACC DATA CREATE(bad_cell)
  !
  cell_loop: &
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    do n = n1,n2
     !if ( ieee_is_negative(  usp(n,nec)   ) .or. &
     !     ieee_is_negative( gpsp(n,npres) ) .or. &
     !     ieee_is_negative( gpsp(n,ntemp) ) ) then
      if (  usp(n,nec)   < dens_limit .or. &
           gpsp(n,ntemp) < temp_limit .or. &
           gpsp(n,npres) < pres_limit ) then
        bad_cell(nc)   = true
        bad_cell_found = true
        cycle cell_loop
      end if
    end do
    !
  end do cell_loop
  !$ACC UPDATE DEVICE(bad_cell)
  !
  if (bad_cell_found) then
    !
    ! Warning message
    !
    write (iout,1) mypnum
    flush (iout)
    !
    ! Use the emergency filter to limit any cells with
    ! a negative density, temperature, or pressure
    !
    call emergency_limit_bad_cells(bad_cell,usp,gpsp)
  end if
  !
  !$ACC END DATA
  if (allocated(bad_cell)) deallocate ( bad_cell )
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/,2x,"WARNING!!! CPU # ",i6," FOUND A CELL WITH ", &
                 "A NEGATIVE DENSITY, TEMPERATURE, OR PRESSURE,",/, &
              2x,"WARNING!!!",14x,"AND IS ATTEMPTING TO USE THE ", &
                 "EMERGENCY LIMITER TO CONTINUE RUNNING!",/)
  !
end subroutine check_cell_thermo_properties
!
!###############################################################################
!> Check face points for bad thermo properties
!> (pressure, temperature, density), limit the cell, and update the face.
!###############################################################################
subroutine check_face_thermo_properties()
  !
  !.. Use Statements ..
  use geovar,  only : ncell,nfbnd,nface,face
  use flowvar, only : usp,gpsp,faceusp,facegpsp
  use ovar,    only : rhoref,tref,aref
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_is_negative
  !
  !.. Local Scalars ..
  integer :: nf,f1,f2,n,left_cell,right_cell
  logical(lk) :: bad_cell_found
  real(wp) :: dens_limit,temp_limit,pres_limit
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable :: bad_cell(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_face_thermo_properties"
  !
continue
  !
  if (ncs > 0) return
  !
  dens_limit = 0.1_wp/rhoref
  temp_limit = 100.0_wp/tref
  pres_limit = 40000.0_wp/(rhoref*aref*aref)
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Find any cells with solution points that have
  ! a negative density, temperature, or pressure
  !
  bad_cell_found = fals
  !
  bad_cell = init_array(ncell,fals) ! F2003 AUTO-REALLOC
  !$ACC DATA CREATE(bad_cell)
  !
  ! Search the boundary faces for any negative density, temperature, or pressure
  !
  bnd_face_loop: &
  do nf = 1,nfbnd
    !
    f1 = face(nf)%beg_fp
    f2 = face(nf)%end_fp
    !
    left_cell = face(nf)%left%cell
    !
    if (.not. bad_cell(left_cell)) then
      !
      do n = f1,f2
       !if ( ieee_is_negative(  faceusp(n,nec  ,1) ) .or. &
       !     ieee_is_negative( facegpsp(n,ntemp,1) ) .or. &
       !     ieee_is_negative( facegpsp(n,npres,1) ) ) then
        if (  faceusp(n,nec  ,1) < dens_limit .or. &
             facegpsp(n,ntemp,1) < temp_limit .or. &
             facegpsp(n,npres,1) < pres_limit ) then
          bad_cell(left_cell) = true
          bad_cell_found      = true
          cycle bnd_face_loop
        end if
      end do
      !
    end if
    !
  end do bnd_face_loop
  !
  ! Search the interior faces for any negative density, temperature, or pressure
  !
  int_face_loop: &
  do nf = nfbnd+1,nface
    !
    f1 = face(nf)%beg_fp
    f2 = face(nf)%end_fp
    !
    left_cell = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    !
    ! Check the left side of this face
    !
    left_cell_check: &
    if (.not. bad_cell(left_cell)) then
      !
      do n = f1,f2
       !if ( ieee_is_negative(  faceusp(n,nec  ,1) ) .or. &
       !     ieee_is_negative( facegpsp(n,ntemp,1) ) .or. &
       !     ieee_is_negative( facegpsp(n,npres,1) ) ) then
        if (  faceusp(n,nec  ,1) < dens_limit .or. &
             facegpsp(n,ntemp,1) < temp_limit .or. &
             facegpsp(n,npres,1) < pres_limit ) then
          bad_cell(left_cell) = true
          bad_cell_found      = true
          exit left_cell_check
        end if
      end do
      !
    end if left_cell_check
    !
    ! Check the right side of this face
    !
    right_cell_check: &
    if (.not. bad_cell(right_cell)) then
      !
      do n = f1,f2
       !if ( ieee_is_negative(  faceusp(n,nec  ,2) ) .or. &
       !     ieee_is_negative( facegpsp(n,ntemp,2) ) .or. &
       !     ieee_is_negative( facegpsp(n,npres,2) ) ) then
        if (  faceusp(n,nec  ,2) < dens_limit .or. &
             facegpsp(n,ntemp,2) < temp_limit .or. &
             facegpsp(n,npres,2) < pres_limit ) then
          bad_cell(right_cell) = true
          bad_cell_found       = true
          exit right_cell_check
        end if
      end do
      !
    end if right_cell_check
    !
  end do int_face_loop
  !$ACC UPDATE DEVICE(bad_cell)
  !
  if (bad_cell_found) then
    !
    ! Warning message
    !
    write (iout,1) mypnum
    flush (iout)
    !
    ! Use the emergency filter to limit any cells with
    ! a negative density, temperature, or pressure
    !
    call emergency_limit_bad_cells(bad_cell,usp,gpsp)
    !
    ! Now go back and for all bad cells, re-interpolate
    ! the conserved variables to the face points
    !
    call reinterpolate_bad_cells(bad_cell,usp,faceusp,facegpsp)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/,2x,"WARNING!!! CPU # ",i6," FOUND A FACE WITH ", &
                 "A NEGATIVE DENSITY, TEMPERATURE, OR PRESSURE,",/, &
              2x,"WARNING!!!",14x,"AND IS ATTEMPTING TO USE THE ", &
                 "EMERGENCY LIMITER TO CONTINUE RUNNING!",/)
  !
end subroutine check_face_thermo_properties
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!> For cells marked as bad, reinterpolate the cell volume solution data to the
!> face points and update the face gas property values.
!###############################################################################
subroutine reinterpolate_bad_cells(bad_cell,usp,faceusp,facegpsp)
  !
  use geovar,     only : nfbnd,nface,face
  use thermo_mod, only : update_thermo_properties
  !
  !.. Optional Arguments ..
  logical(lk),    intent(in) :: bad_cell(:)
  real(wp),       intent(in) :: usp(:,:)
  real(wp),    intent(inout) :: faceusp(:,:,:)
  real(wp),    intent(inout) :: facegpsp(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr,nf,sf,f1,f2
  !
  !.. Local Allocatable Arrays
  logical(lk), allocatable :: faces_to_update(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reinterpolate_bad_cells"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  faces_to_update = init_array(nface,2,fals) ! F2003 AUTO-REALLOC
  !$ACC DATA CREATE(faces_to_update)
  !
  if (force_nd_subroutines) then
    !
    call reinterpolate_bad_cells_nd(1,nq,bad_cell,usp, &
                                    faceusp,faces_to_update)
    !
  else
    !
    if (nr == 3) then
      call reinterpolate_bad_cells_3d(bad_cell,usp, &
                                      faceusp,faces_to_update)
    else if (nr == 2) then
      call reinterpolate_bad_cells_2d(bad_cell,usp, &
                                      faceusp,faces_to_update)
    else
      call reinterpolate_bad_cells_1d(bad_cell,usp, &
                                      faceusp,faces_to_update)
    end if
    !
  end if
  !$ACC UPDATE HOST(faces_to_update)
  !
  ! Update boundary faces, side-of-face=1 for interior solution.
  ! Do not update......... side-of-face=2 for boundary values.
  do nf = 1,nfbnd  ! loop boundary faces
    if (faces_to_update(nf,1)) then
      f1 = face(nf)%beg_fp
      f2 = face(nf)%end_fp
      call update_thermo_properties(f1,f2, &
                                    faceusp(:,:,1),facegpsp(:,:,1), &
                                    loc="face",caller_name=pname)
    end if
  end do
  !
  ! Update interior faces, side-of-face=1:2 for left:right solutions.
  do sf = 1,2  ! loop side of face
    do nf = nfbnd+1,nface  ! loop interior faces
      if (faces_to_update(nf,sf)) then
        f1 = face(nf)%beg_fp
        f2 = face(nf)%end_fp
        call update_thermo_properties(f1,f2, &
                                      faceusp(:,:,sf),facegpsp(:,:,sf), &
                                      loc="face",caller_name=pname)
      end if
    end do
  end do
  !
  !$ACC END DATA
  if (allocated(faces_to_update)) deallocate ( faces_to_update )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine reinterpolate_bad_cells
!
!###############################################################################
!> For cells marked as bad, reinterpolate the cell volume solution data to the
!> face points and update the face gas property values.
!>
!> Optimized routine for 3d.
!###############################################################################
pure &
subroutine reinterpolate_bad_cells_3d(bad_cell,usp, &
                                      faceusp,faces_to_update)
  !
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use geovar,            only : face
  !
  !.. Optional Arguments ..
  logical(lk),    intent(in) :: bad_cell(:)
  real(wp),       intent(in) :: usp(:,:)
  real(wp),    intent(inout) :: faceusp(:,:,:)
  logical(lk), intent(inout) :: faces_to_update(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,face_ec,face_mx,face_my,face_mz,face_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reinterpolate_bad_cells_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(bad_cell,faces_to_update)
  cell_loop: &
  do nc = 1,ncell
    !
    if (.not. bad_cell(nc)) cycle cell_loop
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
      ! Mark this side of this face as needing
      ! an update of the thermo properties
      !
      faces_to_update(nf,sf) = true
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
        faceusp(f0+kf,nec,sf) = face_ec
        faceusp(f0+kf,nmx,sf) = face_mx
        faceusp(f0+kf,nmy,sf) = face_my
        faceusp(f0+kf,nmz,sf) = face_mz
        faceusp(f0+kf,nee,sf) = face_ee
        !
      end do
      !
    end do
    !
  end do cell_loop
  !
end subroutine reinterpolate_bad_cells_3d
!
!###############################################################################
!> For cells marked as bad, reinterpolate the cell volume solution data to the
!> face points and update the face gas property values.
!>
!> Optimized routine for 2d.
!###############################################################################
pure &
subroutine reinterpolate_bad_cells_2d(bad_cell,usp, &
                                      faceusp,faces_to_update)
  !
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use geovar,            only : face
  !
  !.. Optional Arguments ..
  logical(lk),    intent(in) :: bad_cell(:)
  real(wp),       intent(in) :: usp(:,:)
  real(wp),    intent(inout) :: faceusp(:,:,:)
  logical(lk), intent(inout) :: faces_to_update(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,face_ec,face_mx,face_my,face_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reinterpolate_bad_cells_2d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(bad_cell,faces_to_update)
  cell_loop: &
  do nc = 1,ncell
    !
    if (.not. bad_cell(nc)) cycle cell_loop
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
      ! Mark this side of this face as needing
      ! an update of the thermo properties
      !
      faces_to_update(nf,sf) = true
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
        faceusp(f0+kf,nec,sf) = face_ec
        faceusp(f0+kf,nmx,sf) = face_mx
        faceusp(f0+kf,nmy,sf) = face_my
        faceusp(f0+kf,nee,sf) = face_ee
        !
      end do
      !
    end do
    !
  end do cell_loop
  !
end subroutine reinterpolate_bad_cells_2d
!
!###############################################################################
!> For cells marked as bad, reinterpolate the cell volume solution data to the
!> face points and update the face gas property values.
!>
!> Optimized routine for 1d.
!###############################################################################
pure &
subroutine reinterpolate_bad_cells_1d(bad_cell,usp, &
                                      faceusp,faces_to_update)
  !
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use geovar,            only : face
  !
  !.. Optional Arguments ..
  logical(lk),    intent(in) :: bad_cell(:)
  real(wp),       intent(in) :: usp(:,:)
  real(wp),    intent(inout) :: faceusp(:,:,:)
  logical(lk), intent(inout) :: faces_to_update(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,face_ec,face_mx,face_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reinterpolate_bad_cells_1d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(bad_cell,faces_to_update)
  cell_loop: &
  do nc = 1,ncell
    !
    if (.not. bad_cell(nc)) cycle cell_loop
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
      ! Mark this side of this face as needing
      ! an update of the thermo properties
      !
      faces_to_update(nf,sf) = true
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
        faceusp(f0+kf,nec,sf) = face_ec
        faceusp(f0+kf,nmx,sf) = face_mx
        faceusp(f0+kf,nee,sf) = face_ee
        !
      end do
      !
    end do
    !
  end do cell_loop
  !
end subroutine reinterpolate_bad_cells_1d
!
!###############################################################################
!> For cells marked as bad, reinterpolate the cell volume solution data to the
!> face points and update the face gas property values.
!>
!> General routine for nd.
!###############################################################################
pure &
subroutine reinterpolate_bad_cells_nd(mbeg,mend,bad_cell,usp, &
                                      faceusp,faces_to_update)
  !
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use geovar,            only : face
  !
  !.. Optional Arguments ..
  integer,        intent(in) :: mbeg,mend
  logical(lk),    intent(in) :: bad_cell(:)
  real(wp),       intent(in) :: usp(:,:)
  real(wp),    intent(inout) :: faceusp(:,:,:)
  logical(lk), intent(inout) :: faces_to_update(:,:)
  !
  !.. Local Scalars ..
  integer :: m,nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reinterpolate_bad_cells_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(bad_cell,faces_to_update)
  cell_loop: &
  do nc = 1,ncell
    !
    if (.not. bad_cell(nc)) cycle cell_loop
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
      ! Mark this side of this face as needing
      ! an update of the thermo properties
      !
      faces_to_update(nf,sf) = true
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
      !$ACC LOOP VECTOR COLLAPSE(2)
      do kf = 1,nfp
        do m = mbeg,mend
          !
          ! kf: index of the flx-pt local to this face
          ! lf: index of flx-pt within interp%mat
          !
          lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
          !
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
          faceusp(f0+kf,m,sf) = val
          !
        end do
      end do
      !
    end do
    !
  end do cell_loop
  !
end subroutine reinterpolate_bad_cells_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!> Apply the sponge filter to bad cells and
!> update the volume solution and thermo in usp and gpsp.
!###############################################################################
subroutine emergency_limit_bad_cells(bad_cell,usp,gpsp)
  !
  !.. Use Statements ..
  use geovar,     only : ncell,cell
  use thermo_mod, only : update_thermo_properties
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: gpsp(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,nr,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit_bad_cells"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (force_nd_subroutines) then
    !
    call emergency_limit_bad_cells_nd(1,nq,bad_cell,usp)
    !
  else
    !
    if (nr == 3) then
      call emergency_limit_bad_cells_3d(bad_cell,usp)
    else if (nr == 2) then
      call emergency_limit_bad_cells_2d(bad_cell,usp)
    else
      call emergency_limit_bad_cells_1d(bad_cell,usp)
    end if
    !
  end if
  !
  ! Update the thermodynamic properties for any cells that had to be limited
  !
  do nc = 1,ncell
    if (bad_cell(nc)) then
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%beg_sp
      call update_thermo_properties(n1,n2,usp,gpsp, &
                                    loc="cell",caller_name=pname)
    end if
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine emergency_limit_bad_cells
!
!###############################################################################
!> Apply the sponge filter to bad cells and
!> update the volume solution and thermo in usp and gpsp.
!>
!> Optimized routine for 3d.
!###############################################################################
pure &
subroutine emergency_limit_bad_cells_3d(bad_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,     only : ncell,cell
  use filter_mod, only : filter
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_mz,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit_bad_cells_3d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nmz) = usp(n0+n,nmz)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_my = zero
      filt_mz = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
        !
        filt_ec = filt_ec + a * old_usp(k,nec)
        filt_mx = filt_mx + a * old_usp(k,nmx)
        filt_my = filt_my + a * old_usp(k,nmy)
        filt_mz = filt_mz + a * old_usp(k,nmz)
        filt_ee = filt_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nmy) = filt_my
      usp(n0+n,nmz) = filt_mz
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_limit_bad_cells_3d
!
!###############################################################################
!> Apply the sponge filter to bad cells and
!> update the volume solution and thermo in usp and gpsp.
!>
!> Optimized routine for 2d.
!###############################################################################
pure &
subroutine emergency_limit_bad_cells_2d(bad_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,     only : ncell,cell
  use filter_mod, only : filter
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit_bad_cells_2d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_my = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
        !
        filt_ec = filt_ec + a * old_usp(k,nec)
        filt_mx = filt_mx + a * old_usp(k,nmx)
        filt_my = filt_my + a * old_usp(k,nmy)
        filt_ee = filt_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nmy) = filt_my
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_limit_bad_cells_2d
!
!###############################################################################
!> Apply the sponge filter to bad cells and
!> update the volume solution and thermo in usp and gpsp.
!>
!> Optimized routine for 1d.
!###############################################################################
pure &
subroutine emergency_limit_bad_cells_1d(bad_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,     only : ncell,cell
  use filter_mod, only : filter
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit_bad_cells_1d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
        !
        filt_ec = filt_ec + a * old_usp(k,nec)
        filt_mx = filt_mx + a * old_usp(k,nmx)
        filt_ee = filt_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_limit_bad_cells_1d
!
!###############################################################################
!> Apply the sponge filter to bad cells and
!> update the volume solution and thermo in usp and gpsp.
!>
!> General routine for nd.
!###############################################################################
pure &
subroutine emergency_limit_bad_cells_nd(mbeg,mend,bad_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,     only : ncell,cell
  use filter_mod, only : filter
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: mbeg,mend
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit_bad_cells_nd"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        old_usp(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        !
        val = zero
        !
        do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
               filter(this_geom,this_order)%sponge%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%sponge%col_idx(i)
          a = filter(this_geom,this_order)%sponge%col_val(i)
          !
          val = val + a * old_usp(k,m)
          !
        end do
        !
        usp(n0+n,m) = val
        !
      end do
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_limit_bad_cells_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module thermofix_mod
