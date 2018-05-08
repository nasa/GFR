module flux_mod
  !
!#define PROFILE_ON
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer => debug_timer
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,ntk,ntl,nq
  !
  implicit none
  !
  private
  !
  public :: interp_cv_to_flxpts
  public :: interp_dusp_to_flxpts
  public :: flux_gradient_at_solpts_visc
  public :: flux_gradient_at_solpts_visc_gpu
  public :: flux_gradient_at_solpts
  public :: correct_interface_flux_visc
  public :: discontinuous_gradient_at_solpts
  public :: correct_gradients_at_solpts
  public :: correct_gradients_at_solpts_gpu
  public :: correct_gradients_br2
  public :: correct_gradients_br2_gpu
  public :: get_invs_common_fluxes_interior
  public :: get_visc_common_fluxes_interior
  public :: get_invs_common_fluxes_bc
  public :: get_visc_common_fluxes_bc
  public :: init_cbias
  !
  real(wp), save, allocatable :: cbias(:)
  !$ACC DECLARE CREATE(cbias)
  real(wp), save, allocatable :: face_jump(:,:)
  !$ACC DECLARE CREATE(face_jump)
  !
  interface write_array
    module procedure write_2d_array, write_3d_array
  end interface write_array
  !
contains
!
!###############################################################################
!
subroutine interp_cv_to_flxpts(iter,nst)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iter,nst
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_flxpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3) then
    call interp_cv_to_flxpts_3d(iter,nst)
  else if (nr == 2) then
    call interp_cv_to_flxpts_2d
  else
    call interp_cv_to_flxpts_nd
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine interp_cv_to_flxpts
!
!###############################################################################
!
subroutine interp_cv_to_flxpts_3d(iter,nst)
  !
  use order_mod,         only : maxSP,maxFP
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use flowvar,           only : usp
  use geovar,            only : face
  use flowvar,           only : faceusp_array
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iter,nst
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
 !if (iter > 0 .and. nst > 0) then
 !  call write_array(100,iter,nst,faceusp_array(:,:,1))
 !end if
  !
end subroutine interp_cv_to_flxpts_3d
!
!###############################################################################
!
subroutine interp_cv_to_flxpts_2d()
  !
  use order_mod,         only : maxSP,maxFP
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
subroutine interp_cv_to_flxpts_nd()
  !
  use order_mod,         only : maxSP,maxFP
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use flowvar,           only : usp
  use geovar,            only : face
  use flowvar,           only : faceusp_array
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
        do m = 1,nq
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
  if (nr == 3) then
    call interp_dusp_to_flxpts_3d
  else if (nr == 2) then
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
subroutine flux_gradient_at_solpts_visc(time_step,rk_step)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,maxFP,geom_solpts
  use geovar,            only : nr,ncell
  use geovar,            only : cell,fp
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  use metrics_mod,       only : metrics_dt,geofa
  use flowvar,           only : usp,dusp,residual
  use ovar,              only : governing_equations
  use correction_mod,    only : correct
  !
  use visc_flux_mod, only : visc_flx_cv
  use invs_flux_mod, only : invs_flx_cv
  !
  !.. Formal Arguments ..
  integer, intent(in) :: time_step
  integer, intent(in) :: rk_step
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,m,n,nc,n1,n2,np
  integer :: nf,sf,kf,lf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot,nfp,ncf
  real(wp) :: a,ds,val,fval,correct_n
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq,1:nr) :: vflx
  real(wp), dimension(1:nq,1:nr) :: iflx
  real(wp), dimension(1:nq,1:maxSP) :: lcl_res
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: rsflx
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: xyflx
  integer, dimension(1:maxFP) :: fpidx
  real(wp), dimension(1:nq) :: flx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "flux_gradient_at_solpts_visc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
!!!$ACC DATA CREATE
  vflx(:,:) = zero
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
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
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
!  !!$ACC LOOP VECTOR COLLAPSE(2)
!  !!$ACC LOOP WORKER
!   !$ACC LOOP SEQ
    do k = 1,np
      !
      n = k + n1-1
      !
      ! Inviscid fluxes
      !
!    !!$ACC LOOP SEQ
      do l = 1,nr
        call invs_flx_cv( usp(n,:) , l , iflx(:,l) ) ! l-directional flux
      end do
      !
      ! Viscous fluxes
      !
      if (governing_equations == NavierStokes_Eqns) then
        call visc_flx_cv( usp(n,:) , dusp(n,:,:) , vflx )
      end if
      !
      ! Total fluxes
      !
!    !!$ACC LOOP SEQ
      do l = 1,nr
!      !!$ACC LOOP SEQ
        do m = 1,nq
          xyflx(k,m,l) = iflx(m,l) - vflx(m,l)
        end do
      end do
      !
    end do
    !
    ! Transform the fluxes in the x-y physical
    ! domain to the r-s computational domain
    !
!   !$ACC LOOP VECTOR COLLAPSE(3)
    do m = 1,nq
      do l = 1,nr
        do k = 1,np
          val = zero
          do i = 1,nr
            val = val + xyflx(k,m,i)*metrics_dt(nc)%met(k,l,i)
          end do
          rsflx(k,m,l) = val * metrics_dt(nc)%jac(k)
        end do
      end do
    end do
    !
    ! Compute the divergence of the discontinuous fluxes
    !
    ! Zero the residual array for the current cell
    !
!  !!$ACC LOOP VECTOR COLLAPSE(2)
   !do n = 1,maxSP
   !  do m = 1,nq
   !    lcl_res(m,n) = zero
   !  end do
   !end do
    !
    ! Loop through each coordinate direction and add the divergence of
    ! the flux to the local residual array for the current cell
    !
!   !$ACC LOOP VECTOR COLLAPSE(2)
    do k = 1,np
      do m = 1,nq
        val = zero
        do l = 1,nr
          do i = deriv(this_geom,this_order)%d(l)%row_ptr(k)+1, &
                 deriv(this_geom,this_order)%d(l)%row_ptr(k+1)
            !
            n = deriv(this_geom,this_order)%d(l)%row_idx(i)
            a = deriv(this_geom,this_order)%d(l)%row_val(i)
            !
            val = val + a * rsflx(n,m,l)
            !
          end do
        end do
        lcl_res(m,k) = val
      end do
    end do
    !
    ! Copy the local residual for the current cell into the residual array
    !
!   !$ACC LOOP VECTOR COLLAPSE(2)
    do m = 1,nq
      do n = n1,n2
        residual(n,m) = lcl_res(m,n-n1+1)
      end do
    end do
    !
    ! Interpolate the fluxes to the interface flux points and dot this
    ! interpolated vector with the face normal at each flux point
    ! NOTE: We dont need to exchange this information for communication and
    !       periodic boundary faces since this flux is only used within the
    !       current cell.
    !
    ! Loop over the faces of this cell
    !
    do ncf = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      ds = cl(this_geom) * merge(-one,one,sf==1)
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      nfp = geom_solpts(face_geom,face_order)
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
!     !$ACC LOOP VECTOR COLLAPSE(2)
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        do m = 1,nq
          lf = fpidx(kf)
          fval = zero
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
              val = val + a * xyflx(n,m,l)
              !
            end do
            fval = fval + val*geofa(nf)%v(l,kf)
          end do
          do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
                 correct(this_geom,this_order)%row_ptr(lf+1)
            n = n1 - 1 + correct(this_geom,this_order)%row_idx(i)
            correct_n = correct(this_geom,this_order)%row_val(i) * fval * ds
            residual(n,m) = residual(n,m) + correct_n
          end do
        end do
        !
      end do
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine flux_gradient_at_solpts_visc
!
!###############################################################################
!
subroutine flux_gradient_at_solpts_visc_gpu(time_step,rk_step)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,maxFP,geom_solpts
  use geovar,            only : nr,ncell
  use geovar,            only : cell,face,fp,n_solpts
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  use metrics_mod,       only : metrics_array,jac,facenrm
  use flowvar,           only : usp,dusp,residual
  use ovar,              only : governing_equations
  use correction_mod,    only : correct
  !
  use visc_flux_mod, only : visc_flx_cv_gpu
  use invs_flux_mod, only : invs_flx_cv_gpu
  !
  !.. Formal Arguments ..
  integer, intent(in) :: time_step
  integer, intent(in) :: rk_step
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,m,n,nc,n1,n2,np,f0
  integer :: nf,sf,kf,lf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: ierr,io
  integer :: fp_offset,face_rot,nfp,ncf,num_cell_faces
  real(wp) :: a,val,flxvel
  real(wp) :: mvsc,mcdt,ave_ext_rate
  real(wp) :: cv_inv,kflx
  real(wp) :: fval,correct_n,ds
  !
  character(len=100) :: flux_fmt
  character(len=300) :: spf_fname
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:nq,1:maxSP) :: lcl_res
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: vflx
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: rsflx
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: xyflx
  integer, dimension(1:maxFP) :: fpidx
  real(wp), dimension(1:nq) :: flx
  real(wp), dimension(1:maxSP) :: pres
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: dpvdx
  real(wp), dimension(1:nq) :: var
  real(wp), dimension(1:nq,1:nr) :: vflux
  real(wp), dimension(1:nq,1:nr) :: iflux
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tflx
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: rflx
  real(wp), dimension(1:maxFP) :: area
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tfluxes(:,:,:)
  real(wp), allocatable :: ifluxes(:,:,:)
  real(wp), allocatable :: vfluxes(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "flux_gradient_at_solpts_visc_gpu"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( tfluxes(1:n_solpts,1:nr,1:nq) )
  allocate ( ifluxes(1:n_solpts,1:nr,1:nq) )
  allocate ( vfluxes(1:n_solpts,1:nr,1:nq) )
  !$ACC DATA CREATE(tfluxes,ifluxes,vfluxes)
  !
  call invs_flx_cv_gpu( usp , ifluxes )
  !
  if (governing_equations == NavierStokes_Eqns) then
    !
    call visc_flx_cv_gpu( 1 , n_solpts , usp , dusp , vfluxes )
    !
  else
    !
    !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(vfluxes)
    do m = 1,nq
      do l = 1,nr
        do n = 1,n_solpts
          vfluxes(n,l,m) = zero
        end do
      end do
    end do
    !
  end if
  !
  !$ACC PARALLEL LOOP COLLAPSE(3) &
  !$ACC&              PRESENT(tfluxes,ifluxes,vfluxes)
  do m = 1,nq
    do l = 1,nr
      do n = 1,n_solpts
        tfluxes(n,l,m) = ifluxes(n,l,m) - vfluxes(n,l,m)
      end do
    end do
  end do
  !
  !#############################################################################
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(tfluxes) &
  !$ACC&                   PRESENT(deriv) &
  !$ACC&                   PRESENT(correct) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(residual) &
  !$ACC&                   PRESENT(metrics_array) &
  !$ACC&                   PRESENT(jac) &
  !$ACC&                   PRESENT(facenrm) &
  !$ACC&                   PRESENT(cl) &
  !$ACC& PRIVATE(rflx,lcl_res,tflx,fpidx,area)
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    np = n2-n1+1
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do k = 1,np
      do l = 1,nr
        do m = 1,nq
          n = k+n1-1
          tflx(k,l,m) = tfluxes(n,l,m)
        end do
      end do
    end do
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do m = 1,nq
      do l = 1,nr
        do k = 1,np
          val = zero
          do i = 1,nr
            val = val + tflx(k,i,m)*metrics_array(n1-1+k,l,i)
          end do
          rflx(k,l,m) = val * jac(n1-1+k)
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
      do m = 1,nq
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
        lcl_res(m,k) = val
      end do
    end do
    !
    ! Copy the local residual for the current cell into the residual array
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do n = n1,n2
      do m = 1,nq
        residual(n,m) = lcl_res(m,n-n1+1)
      end do
    end do
    !
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
      ds = cl(this_geom) * merge(-one,one,sf==1)
      f0 = face(nf)%beg_fp-1
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      nfp = geom_solpts(face_geom,face_order)
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
      !$ACC LOOP VECTOR COLLAPSE(2)
      do kf = 1,nfp
        do m = 1,nq
          fval = zero
          lf = fpidx(kf)
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
              val = val + a * tflx(n,l,m)
              !
            end do
            fval = fval + val * facenrm(l,f0+kf)
          end do
          do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
                 correct(this_geom,this_order)%row_ptr(lf+1)
            n = n1 - 1 + correct(this_geom,this_order)%row_idx(i)
            correct_n = correct(this_geom,this_order)%row_val(i) * fval * ds
            residual(n,m) = residual(n,m) + correct_n
          end do
        end do
      end do
      !
    end do
    !
  end do
  !
  !$ACC END DATA
  deallocate ( vfluxes )
  deallocate ( ifluxes )
  deallocate ( tfluxes )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine flux_gradient_at_solpts_visc_gpu
!
!###############################################################################
!
subroutine flux_gradient_at_solpts(iter,nst)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts
  use flowvar, only : usp,dusp
  use ovar,    only : governing_equations
  !
  use visc_flux_mod, only : visc_flx_cv_gpu
  use invs_flux_mod, only : invs_flx_cv_gpu
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
  call invs_flx_cv_gpu( usp , ifluxes )
  !
  if (governing_equations == NavierStokes_Eqns) then
    !
    call visc_flx_cv_gpu( 1 , n_solpts , usp , dusp , vfluxes )
    !
  else
    !
    !$ACC PARALLEL LOOP COLLAPSE(3)
    do m = 1,nq
      do l = 1,nr
        do n = 1,n_solpts
          vfluxes(n,l,m) = zero
        end do
      end do
    end do
    !
  end if
  !
  !$ACC PARALLEL LOOP COLLAPSE(3)
  do m = 1,nq
    do l = 1,nr
      do n = 1,n_solpts
        ifluxes(n,l,m) = ifluxes(n,l,m) - vfluxes(n,l,m)
      end do
    end do
  end do
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
  if (nr == 3) then
    call transform_fluxes_3d(tfluxes,rfluxes)
  else if (nr == 2) then
    call transform_fluxes_2d(tfluxes,rfluxes)
  else
    call transform_fluxes_nd(tfluxes,rfluxes)
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
  real(wp) :: d1dx,d2dx,d3dx
  real(wp) :: d1dy,d2dy,d3dy
  real(wp) :: d1dz,d2dz,d3dz
  real(wp) :: xflx,yflx,zflx
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
    xflx = tfluxes(n,1,nmz)
    yflx = tfluxes(n,2,nmz)
    !
    rfluxes(n,1,nmz) = xflx*d1dx + yflx*d1dy
    rfluxes(n,2,nmz) = xflx*d2dx + yflx*d2dy
    !
    xflx = tfluxes(n,1,nee)
    yflx = tfluxes(n,2,nee)
    !
    rfluxes(n,1,nee) = xflx*d1dx + yflx*d1dy
    rfluxes(n,2,nee) = xflx*d2dx + yflx*d2dy
    !
  end do
  !
end subroutine transform_fluxes_2d
!
!###############################################################################
!
subroutine transform_fluxes_nd(tfluxes,rfluxes)
  !
  !.. Use Statements ..
  use metrics_mod, only : metrics_array,jac
  !
  !.. Formal Arguments ..
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
  do m = 1,nq
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
  if (nr == 3) then
    call compute_flux_divergence_3d(rfluxes,iter,nst)
  else if (nr == 2) then
    call compute_flux_divergence_2d(rfluxes)
  else
    call compute_flux_divergence_nd(rfluxes)
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
end subroutine compute_flux_divergence_2d
!
!###############################################################################
!
subroutine compute_flux_divergence_nd(rfluxes)
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
    do m = 1,nq
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
      do m = 1,nq
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
  if (nr == 3) then
    call add_flux_correction_3d(tfluxes)
  else if (nr == 2) then
    call add_flux_correction_2d(tfluxes)
  else
    call add_flux_correction_nd(tfluxes)
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
 !call stop_recon(stop_mpi,pname,__LINE__,__FILE__,"error")
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_flux_correction_2d"
  !
continue
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
          residual(n0+n,nec) = residual(n0+n,nec) + a*flx_ec
          residual(n0+n,nmx) = residual(n0+n,nmx) + a*flx_mx
          residual(n0+n,nmy) = residual(n0+n,nmy) + a*flx_my
          residual(n0+n,nee) = residual(n0+n,nee) + a*flx_ee
          !
        end do
        !
      end do pts_on_face
      !
    end do faces_of_cell
    !
  end do cell_loop
  !
end subroutine add_flux_correction_2d
!
!###############################################################################
!
subroutine add_flux_correction_nd(tfluxes)
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_flux_correction_nd"
  !
continue
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
        do m = 1,nq
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
          end do
        end do
      end do
      !
    end do
    !
  end do
  !
end subroutine add_flux_correction_nd
!
!###############################################################################
!
subroutine correct_interface_flux_visc(time_step,rk_step)
  !
  use order_mod,         only : maxFP,geom_solpts
  use geovar,            only : nr,nfbnd,nface
  use geovar,            only : bface,face,cell,fp
  use ovar,              only : invs_flux_method,governing_equations
  use ovar,              only : output_dir,dump_fluxes
  use metrics_mod,       only : geofa
  use correction_mod,    only : correct
  use flowvar,           only : residual
  use flowvar,           only : faceusp_array,facedusp_array
  !
  use visc_flux_mod, only : common_visc_flx_cv
  use visc_flux_mod, only : visc_wall_flx_cv
  use invs_flux_mod, only : invs_wall_flx_cv
  !
  use upwind_flux_mod, only : upwind_flux_cv
  use upwind_flux_mod, only : hllc_upwind_flux_cv
  use upwind_flux_mod, only : ldfss_upwind_flux_cv
  use upwind_flux_mod, only : rotated_rhll_upwind_flux_cv_2d
  use upwind_flux_mod, only : rotated_rhll_upwind_flux_cv_3d
  !
  !.. Formal Arguments ..
  integer, intent(in) :: time_step
  integer, intent(in) :: rk_step
  !
  !.. Local Scalars ..
  integer  :: i,k,l,m,n,nf,r1,r2,l1,l2
  integer  :: kl,kr,iu,nfp,f0
  integer  :: ierr,iofp
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: ds,dsl,dsr,cmat_val
  !
  character(len=300) :: fpf_fname
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:nq) :: iflx,vflx
  real(wp), dimension(1:nr,1:maxFP) :: rnrm
  real(wp), dimension(1:nq,1:maxFP) :: lflx,rflx
  real(wp), dimension(1:nq,1:maxFP) :: luqp,ruqp,fncq
  real(wp), dimension(1:nr,1:nq,1:maxFP) :: lduqp,rduqp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_interface_flux_visc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! ##############################################################
  ! ##############################################################
  ! ##############################################################
  ! ##########              BOUNDARY FACES              ##########
  ! ##############################################################
  ! ##############################################################
  ! ##############################################################
  !
  ! Loop through the boundary faces to apply the correction
  ! due to the boundary fluxes.
  !
  boundary_fluxes: do nf = 1,nfbnd
    !
    face_geom  = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    host_cell  = face(nf)%left%cell    ! host cell on boundary face
    host_geom  = cell(host_cell)%geom  ! geometry of host cell
    host_order = cell(host_cell)%order ! order of host cell
    !
    l1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in left cell
    l2 = cell(host_cell)%end_sp ! ending    index for sol-pts in left cell
    !
    do k = 1,nfp
      do m = 1,nq
        ! FP solution from left  cell on face
        luqp(m,k) = faceusp_array(f0+k,m,1)
        ! FP solution from right cell on face
        ruqp(m,k) = faceusp_array(f0+k,m,2)
      end do
    end do
    !
    if (governing_equations == NavierStokes_Eqns) then
      do k = 1,nfp
        do m = 1,nq
          do l = 1,nr
            ! Solution gradient at FP from left  cell on face
            lduqp(l,m,k) = facedusp_array(f0+k,l,m,1)
            ! Solution gradient at FP from right cell on face
            rduqp(l,m,k) = facedusp_array(f0+k,l,m,2)
          end do
        end do
      end do
    end if
    !
    rnrm(1:nr,1:nfp) = geofa(nf)%v(1:nr,1:nfp) ! Area face normals at FP
    !
    ! Get the common flux at each flux point
    !
    do k = 1,nfp
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( rnrm(:,k) )
      !
      ! Get the inviscid upwind flux at the current flux point
      !
      if (any(bface(1,nf) == bc_walls)) then
        !
        ! This is a wall boundary so ensure that the inviscid
        ! boundary flux is only due to the pressure term.
        !
        iflx(:) = invs_wall_flx_cv( unrm(:) , luqp(:,k) )
        !
        ! Compute the viscous wall boundary flux
        !
        if (governing_equations == NavierStokes_Eqns) then
          vflx(:) = visc_wall_flx_cv( unrm(:) , ruqp(:,k) , rduqp(:,:,k) )
        end if
        !
      else
        !
        ! For all other boundaries, compute the boundary
        ! fluxes in the same way as the interior faces.
        !
        if (invs_flux_method == invs_flux_ROE) then
          iflx(:) = upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
        else if (invs_flux_method == invs_flux_HLLC) then
          iflx(:) = hllc_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
        else if (invs_flux_method == invs_flux_LDFSS) then
          iflx(:) = ldfss_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
        else if (invs_flux_method == invs_flux_Rotated_RHLL) then
          if (nr == 2) then
            iflx(:) = rotated_rhll_upwind_flux_cv_2d ( unrm(:) , luqp(:,k) , &
                                                                 ruqp(:,k) )
          else
            iflx(:) = rotated_rhll_upwind_flux_cv_3d ( unrm(:) , luqp(:,k) , &
                                                                 ruqp(:,k) )
          end if
        end if
        !
        ! Get the common viscous flux at the current flux point
        !
        if (governing_equations == NavierStokes_Eqns) then
          vflx(:) = common_visc_flx_cv( unrm(:), luqp(:,k), lduqp(:,:,k), &
                                                 ruqp(:,k), rduqp(:,:,k) )
        end if
        !
      end if
      !
      ! Get the total common boundary flux at this flux point
      !
      fncq(:,k) = iflx(:) - vflx(:)
      !
    end do
    !
    ! Apply the flux correction due to the flux jump
    ! NOTE: We only need to apply the correction to the left cell because
    !       this is a boundary face and the right cell is a ghost cell.
    !
    do k = 1,nfp
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kl = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(nf)%left%rotation)% &
                                     pts(k)
      !
      ! Get the contribution to the face length/area from the current flux
      ! point for the left cell
      !
      ds  = norm2( rnrm(:,k) )
      dsl = ds * cl(host_geom)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! host_cell
      !
      do i = correct(host_geom,host_order)%row_ptr(kl)+1, &
             correct(host_geom,host_order)%row_ptr(kl+1)
        n = l1 - 1 + correct(host_geom,host_order)%row_idx(i)
        cmat_val = dsl * correct(host_geom,host_order)%row_val(i)
        do m = 1,nq
          residual(n,m) = residual(n,m) + fncq(m,k) * cmat_val
        end do
      end do
      !
    end do
    !
  end do boundary_fluxes
  !
  ! ##############################################################
  ! ##############################################################
  ! ##############################################################
  ! ##########              INTERIOR FACES              ##########
  ! ##############################################################
  ! ##############################################################
  ! ##############################################################
  !
  ! Loop over the interior faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  interior_fluxes: do nf = nfbnd+1,nface
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0 ! number of flux points on this face
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
    l1 = cell(left_cell)%beg_sp ! beginning index for sol-pts in left cell
    l2 = cell(left_cell)%end_sp ! ending    index for sol-pts in left cell
    !
    r1 = cell(right_cell)%beg_sp ! beginning index for sol-pts in right cell
    r2 = cell(right_cell)%end_sp ! ending    index for sol-pts in right cell
    !
    do k = 1,nfp
      do m = 1,nq
        ! FP solution from left  cell on face
        luqp(m,k) = faceusp_array(f0+k,m,1)
        ! FP solution from right cell on face
        ruqp(m,k) = faceusp_array(f0+k,m,2)
      end do
    end do
    !
    if (governing_equations == NavierStokes_Eqns) then
      do k = 1,nfp
        do m = 1,nq
          do l = 1,nr
            ! Solution gradient at FP from left  cell on face
            lduqp(l,m,k) = facedusp_array(f0+k,l,m,1)
            ! Solution gradient at FP from right cell on face
            rduqp(l,m,k) = facedusp_array(f0+k,l,m,2)
          end do
        end do
      end do
    end if
    !
    rnrm(1:nr,1:nfp) = geofa(nf)%v(1:nr,1:nfp)
    !
    ! Get the common flux at each flux point
    !
    do k = 1,nfp
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( rnrm(:,k) )
      !
      ! Get the inviscid upwind flux at the current flux point
      !
      if (invs_flux_method == invs_flux_ROE) then
        iflx(:) = upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (invs_flux_method == invs_flux_HLLC) then
        iflx(:) = hllc_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (invs_flux_method == invs_flux_LDFSS) then
        iflx(:) = ldfss_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (invs_flux_method == invs_flux_Rotated_RHLL) then
        if (nr == 2) then
          iflx(:) = rotated_rhll_upwind_flux_cv_2d ( unrm(:) , luqp(:,k) , &
                                                               ruqp(:,k) )
        else
          iflx(:) = rotated_rhll_upwind_flux_cv_3d ( unrm(:) , luqp(:,k) , &
                                                               ruqp(:,k) )
        end if
      end if
      !
      ! Get the common viscous flux at the current flux point
      !
      if (governing_equations == NavierStokes_Eqns) then
        vflx(:) = common_visc_flx_cv( unrm(:), luqp(:,k), lduqp(:,:,k), &
                                               ruqp(:,k), rduqp(:,:,k) )
      end if
      !
      ! Get the total common boundary flux at this flux point
      !
      fncq(:,k) = iflx(:) - vflx(:)
      !
    end do
    !
    ! Apply the flux correction due to the flux jump
    !
    do k = 1,nfp
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kl = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(nf)%left%rotation)% &
                                     pts(k)
      kr = face(nf)%right%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(nf)%right%rotation)% &
                                      pts(k)
      !
      ! Get the face length for the left and right cell
      !
      ds  = norm2( rnrm(:,k) )
      dsl = ds * cl(left_geom)
      dsr = ds * cl(right_geom)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! left_cell
      !
      do i = correct(left_geom,left_order)%row_ptr(kl)+1, &
             correct(left_geom,left_order)%row_ptr(kl+1)
        n = l1 - 1 + correct(left_geom,left_order)%row_idx(i)
        cmat_val = dsl * correct(left_geom,left_order)%row_val(i)
        do m = 1,nq
          residual(n,m) = residual(n,m) + fncq(m,k) * cmat_val
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! right_cell
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      do i = correct(right_geom,right_order)%row_ptr(kr)+1, &
             correct(right_geom,right_order)%row_ptr(kr+1)
        n = r1 - 1 + correct(right_geom,right_order)%row_idx(i)
        cmat_val = dsr * correct(right_geom,right_order)%row_val(i)
        do m = 1,nq
          residual(n,m) = residual(n,m) - fncq(m,k) * cmat_val
        end do
      end do
      !
    end do
    !
  end do interior_fluxes
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine correct_interface_flux_visc
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
  use geovar,      only : nr,n_bnd_flxpts,n_flxpts
  use ovar,        only : invs_flux_method
  use metrics_mod, only : facenrm
  use flowvar,     only : faceusp_array
  use flowvar,     only : faceflx
  !
  use upwind_flux_mod, only : roe_upwind_flux_gpu
  use upwind_flux_mod, only : hllc_upwind_flux_gpu
  use upwind_flux_mod, only : ldfss_upwind_flux_gpu
  use upwind_flux_mod, only : rotated_rhll_upwind_flux_gpu
  !
  !.. Local Scalars ..
  integer :: n,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_invs_common_fluxes_interior"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the inviscid upwind flux at all of the interior flux points
  !
  n1 = n_bnd_flxpts+1
  n2 = n_flxpts
#ifdef DEBUG_FLUXES
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  do n = 0,ncpu-1
    if (mypnum == n) then
      write (iout,100) mypnum,n1,n2
      100 format ("[",i4,"] : (n_bnd_flxpts+1,n_flxpts) = (",i0,",",i0,")")
      flush (iout)
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  if (invs_flux_method == invs_flux_ROE) then
    !
    call roe_upwind_flux_gpu(n1,n2,facenrm,faceusp_array,faceflx)
    !
  else if (invs_flux_method == invs_flux_HLLC) then
    !
    call hllc_upwind_flux_gpu(n1,n2,facenrm,faceusp_array,faceflx)
    !
  else if (invs_flux_method == invs_flux_LDFSS) then
    !
    call ldfss_upwind_flux_gpu(n1,n2,facenrm,faceusp_array,faceflx)
    !
  else if (invs_flux_method == invs_flux_Rotated_RHLL) then
    !
    call rotated_rhll_upwind_flux_gpu(n1,n2,facenrm,faceusp_array,faceflx)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_invs_common_fluxes_interior
!
!###############################################################################
!
subroutine get_visc_common_fluxes_interior()
  !
  use geovar,      only : nr,n_bnd_flxpts,n_flxpts
  use metrics_mod, only : facenrm
  use flowvar,     only : faceusp_array
  use flowvar,     only : facedusp_array
  use flowvar,     only : faceflx
  use ovar,        only : visc_variables_to_ave
  !
  use visc_flux_mod, only : common_vflx_ave_fluxes_gpu
  use visc_flux_mod, only : common_vflx_ave_sol_and_grad_gpu
  !
  !.. Local Scalars ..
  integer :: n,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_visc_common_fluxes_interior"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  n1 = n_bnd_flxpts+1
  n2 = n_flxpts
#ifdef DEBUG_FLUXES
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  do n = 0,ncpu-1
    if (mypnum == n) then
      write (iout,100) mypnum,n1,n2
      100 format ("[",i4,"] : (n_bnd_flxpts+1,n_flxpts) = (",i0,",",i0,")")
      flush (iout)
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    call common_vflx_ave_fluxes_gpu(n1,n2,facenrm,faceusp_array, &
                                    facedusp_array,faceflx)
    !
  else ! visc_variables_to_ave == visc_solution_and_gradients
    !
    call common_vflx_ave_sol_and_grad_gpu(n1,n2,facenrm,faceusp_array, &
                                          facedusp_array,faceflx)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_visc_common_fluxes_interior
!
!###############################################################################
!
subroutine get_invs_common_fluxes_bc()
  !
  use geovar,      only : nr,n_wall_flxpts,n_bnd_flxpts
  use ovar,        only : invs_flux_method
  use metrics_mod, only : facenrm
  use flowvar,     only : faceusp_array
  use flowvar,     only : faceflx
  !
  use invs_flux_mod, only : invs_wall_flx_gpu
  !
  use upwind_flux_mod, only : roe_upwind_flux_gpu
  use upwind_flux_mod, only : hllc_upwind_flux_gpu
  use upwind_flux_mod, only : ldfss_upwind_flux_gpu
  use upwind_flux_mod, only : rotated_rhll_upwind_flux_gpu
  !
  !.. Local Scalars ..
  integer :: i,j,n,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_invs_common_fluxes_bc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
#ifdef DEBUG_FLUXES
  if (mypnum == 998) then
    do n = 1,n_bnd_flxpts
      write (mypnum+5000,102) ((faceusp_array(n,i,j),i=1,nq),j=1,2)
    end do
  end if
  102 format (5es10.2,6x,5es10.2)
#endif
  !
  ! Get the inviscid common flux for all flux points on wall BCs
  !
  n1 = 1
  n2 = n_wall_flxpts
#ifdef DEBUG_FLUXES
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  do n = 0,ncpu-1
    if (mypnum == n) then
      if (n2 >= n1) write (iout,100) mypnum,n1,n2
      100 format ("[",i4,"] : (1,n_wall_flxpts) = (",i0,",",i0,")")
      flush (iout)
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  if (n2 >= n1) then
    !
    call invs_wall_flx_gpu(n1,n2,facenrm,faceusp_array,faceflx)
    !
  end if
  !
  ! For all other boundaries, compute the boundary
  ! fluxes in the same way as the interior faces.
  !
  n1 = n_wall_flxpts+1
  n2 = n_bnd_flxpts
#ifdef DEBUG_FLUXES
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  do n = 0,ncpu-1
    if (mypnum == n) then
      if (n2 >= n1) write (iout,101) mypnum,n1,n2
      101 format ("[",i4,"] : (n_wall_flxpts+1,n_bnd_flxpts) = (",i0,",",i0,")")
      flush (iout)
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  if (invs_flux_method == invs_flux_ROE) then
    !
    call roe_upwind_flux_gpu(n1,n2,facenrm,faceusp_array,faceflx)
    !
  else if (invs_flux_method == invs_flux_HLLC) then
    !
    call hllc_upwind_flux_gpu(n1,n2,facenrm,faceusp_array,faceflx)
    !
  else if (invs_flux_method == invs_flux_LDFSS) then
    !
    call ldfss_upwind_flux_gpu(n1,n2,facenrm,faceusp_array,faceflx)
    !
  else if (invs_flux_method == invs_flux_Rotated_RHLL) then
    !
    call rotated_rhll_upwind_flux_gpu(n1,n2,facenrm,faceusp_array,faceflx)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_invs_common_fluxes_bc
!
!###############################################################################
!
subroutine get_visc_common_fluxes_bc()
  !
  use geovar,      only : nr,n_wall_flxpts,n_bnd_flxpts
  use metrics_mod, only : facenrm
  use flowvar,     only : faceusp_array
  use flowvar,     only : facedusp_array
  use flowvar,     only : faceflx
  use ovar,        only : visc_variables_to_ave
  !
  use visc_flux_mod, only : visc_wall_flx_gpu
  use visc_flux_mod, only : common_vflx_ave_fluxes_gpu
  use visc_flux_mod, only : common_vflx_ave_sol_and_grad_gpu
  !
  !.. Local Scalars ..
  integer :: n,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_visc_common_fluxes_bc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the viscous common flux for all flux points on wall BCs
  !
  n1 = 1
  n2 = n_wall_flxpts
#ifdef DEBUG_FLUXES
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  do n = 0,ncpu-1
    if (mypnum == n) then
      if (n2 >= n1) write (iout,100) mypnum,n1,n2
      100 format ("[",i4,"] : (1,n_wall_flxpts) = (",i0,",",i0,")")
      flush (iout)
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  if (n2 >= n1) then
    !
    call visc_wall_flx_gpu(n1,n2,facenrm,faceusp_array, &
                           facedusp_array,faceflx)
    !
  end if
  !
  ! For all other boundaries, compute the boundary
  ! fluxes in the same way as the interior faces.
  !
  n1 = n_wall_flxpts+1
  n2 = n_bnd_flxpts
#ifdef DEBUG_FLUXES
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  do n = 0,ncpu-1
    if (mypnum == n) then
      if (n2 >= n1) write (iout,101) mypnum,n1,n2
      101 format ("[",i4,"] : (n_wall_flxpts+1,n_bnd_flxpts) = (",i0,",",i0,")")
      flush (iout)
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    call common_vflx_ave_fluxes_gpu(n1,n2,facenrm,faceusp_array, &
                                    facedusp_array,faceflx)
    !
  else ! visc_variables_to_ave == visc_solution_and_gradients
    !
    call common_vflx_ave_sol_and_grad_gpu(n1,n2,facenrm,faceusp_array, &
                                          facedusp_array,faceflx)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_visc_common_fluxes_bc
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
  use geovar,          only : nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "discontinuous_gradient_at_solpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3) then
    call discontinuous_gradient_at_solpts_3d
  else if (nr == 2) then
    call discontinuous_gradient_at_solpts_2d
  else
    call discontinuous_gradient_at_solpts_nd
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
subroutine discontinuous_gradient_at_solpts_nd()
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
      do m = 1,nq
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
      do m = 1,nq
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
  if (nr == 3) then
    call correct_gradients_at_solpts_gpu_3d
  else if (nr == 2) then
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
#define ENABLE_BR2_3D
#define ENABLE_BR2_2D
#define ENABLE_BR2_ND
!#define ENABLE_OLD_BR2
  if (nr == 3) then
#define USE_NEW_GPU_BR2
#ifdef USE_NEW_GPU_BR2
    call correct_gradients_br2_gpu_3d_new
#else
    call correct_gradients_br2_gpu_3d
#endif
  else if (nr == 2) then
    call correct_gradients_br2_gpu_2d
  else
    call correct_gradients_br2_gpu_nd
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine correct_gradients_br2_gpu
!
!###############################################################################
!
subroutine correct_gradients_br2_gpu_3d_new()
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
#ifndef DISABLE_FOR_GPU
  integer :: cell_info(1:6,1:ncell)
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_br2_gpu_3d_new"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
#ifdef ENABLE_BR2_3D
  !
  !$ACC PARALLEL LOOP &
  !$ACC&              PRESENT(faceusp_array,cbias,face_jump)
  do n = 1,n_flxpts
    a = cbias(n)
    face_jump(n,nec) = a*( faceusp_array(n,nec,2) - faceusp_array(n,nec,1) )
    face_jump(n,nmx) = a*( faceusp_array(n,nmx,2) - faceusp_array(n,nmx,1) )
    face_jump(n,nmy) = a*( faceusp_array(n,nmy,2) - faceusp_array(n,nmy,1) )
    face_jump(n,nmz) = a*( faceusp_array(n,nmz,2) - faceusp_array(n,nmz,1) )
    face_jump(n,nee) = a*( faceusp_array(n,nee,2) - faceusp_array(n,nee,1) )
  end do
  !
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
  cell_loop: do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    n0 = n1-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
#ifndef DISABLE_FOR_GPU
   !write (150,101) nc
   !write (160,101) nc
#endif
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
#ifndef DISABLE_FOR_GPU
         !write (150,100) ncf,kf,lf,i,n
#endif
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
#ifndef DISABLE_FOR_GPU
         !write (160,100) ncf,kf,lf,i,n
#endif
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
#ifndef DISABLE_FOR_GPU
 !call stop_recon(stop_mpi,pname,__LINE__,__FILE__,"end of br2_gpu_3d_new")
#endif
  !
#endif
  !
  100 format(i4,4i5)
  101 format(/,8x,"Cell # ",i0,/, &
               " ncf   kf   lf    i    n")
  !
end subroutine correct_gradients_br2_gpu_3d_new
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
#ifdef ENABLE_BR2_3D
  !
  !$ACC PARALLEL LOOP &
  !$ACC&              PRESENT(faceusp_array,cbias,face_jump)
  do n = 1,n_flxpts
    a = cbias(n)
    face_jump(n,nec) = a*( faceusp_array(n,nec,2) - faceusp_array(n,nec,1) )
    face_jump(n,nmx) = a*( faceusp_array(n,nmx,2) - faceusp_array(n,nmx,1) )
    face_jump(n,nmy) = a*( faceusp_array(n,nmy,2) - faceusp_array(n,nmy,1) )
    face_jump(n,nmz) = a*( faceusp_array(n,nmz,2) - faceusp_array(n,nmz,1) )
    face_jump(n,nee) = a*( faceusp_array(n,nee,2) - faceusp_array(n,nee,1) )
  end do
  !
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
  !
#endif
  !
end subroutine correct_gradients_br2_gpu_3d
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
#ifdef ENABLE_BR2_2D
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
#endif
  !
end subroutine correct_gradients_br2_gpu_2d
!
!###############################################################################
!
subroutine correct_gradients_br2_gpu_nd()
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
#ifdef ENABLE_BR2_ND
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(faceusp_array,cbias,face_jump)
  do m = 1,nq
    do n = 1,n_flxpts
      face_jump(n,m) = cbias(n)*( faceusp_array(n,m,2) - faceusp_array(n,m,1) )
    end do
  end do
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
        do m = 1,nq
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
        do m = 1,nq
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
    do m = 1,nq
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
  !
#endif
  !
end subroutine correct_gradients_br2_gpu_nd
!
!###############################################################################
!
subroutine correct_gradients_br2()
  !
  use generic_types_mod, only : matrix
  use order_mod,         only : maxSP,maxFP,geom_solpts
  use geovar,            only : cell
  use geovar,            only : nr,ncell,nfbnd,bface
  use geovar,            only : face,fp
  use metrics_mod,       only : geofa,metrics_dt
  use correction_mod,    only : correct
  use flowvar,           only : dusp,facedusp_array,faceusp_array
  use ovar,              only : walls_are_exact,LDG_beta
  use interpolation_mod, only : interp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,l,m,f0
  integer :: nf,kf,lf,sf,ncf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  real(wp) :: dsl,bias,hfpb,correct_n,correct_m
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:nq) :: lufp,rufp,ujump
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: cell_dusp
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: face_dusp
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: dusp_copy
  !
  integer :: i,i1,i2
  type(matrix) :: cmat
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_br2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
#ifdef ENABLE_OLD_BR2
  cell_loop: do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Store a local copy of the correction matrices divided by the value of
    ! the jacobian for each solution point. This will help reduce cache
    ! misses and prevent repetitive instances of division by the jacobian.
    !
    cmat = correct(this_geom,this_order)
    !
    do i = 1,size(cmat%row_ptr)-1
      i1 = cmat%row_ptr(i)+1
      i2 = cmat%row_ptr(i+1)
      cmat%row_val(i1:i2) = cmat%row_val(i1:i2) / &
                            metrics_dt(nc)%jac( cmat%row_idx(i1:i2) )
    end do
    !
    ! Create a local copy of dusp for this cell to which all the
    ! corrections will be applied. This will help reduce cache misses.
    !
    cell_dusp(1:np,1:nr,1:nq) = dusp(n1:n2,1:nr,1:nq)
   !cell_dusp(1:nr,1:nq,1:np) = dusp(1:nr,1:nq,n1:n2)
    !
    ! Create a copy of the uncorrected dusp for this cell which
    ! will be used to reset the face_dusp array in the cell_faces loop. This
    ! will help reduce cache misses.
    !
    dusp_copy(1:np,1:nr,1:nq) = cell_dusp(1:np,1:nr,1:nq)
   !do n = 1,np
   !  dusp_copy(n,1:nr,1:nq) = cell_dusp(1:nr,1:nq,n)
   !end do
    !
    ! Now Interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    cell_faces: do ncf = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      !
      f0 = face(nf)%beg_sp - 1
      !
      ! Reset face_dusp to the original discontinuous gradient
      !
      face_dusp(1:np,1:nr,1:nq) = dusp_copy(1:np,1:nr,1:nq)
      !
      ! Loop over the points on this face
      !
      face_points: do kf = 1,geom_solpts(face_geom,face_order)
        !
        ! kf: index of flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        lf = cell(nc)%face(ncf)%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(cell(nc)%face(ncf)%rotation)% &
                                            pts(kf)
        !
        ! Solution at this flux point for the left and right cells
        !
        do m = 1,nq
          lufp(m) = faceusp_array(f0+kf,m,1)
          rufp(m) = faceusp_array(f0+kf,m,2)
        end do
        !
        ! Get the unit normal at this flux point
        !
        unrm(1:nr) = unit_vector( geofa(nf)%v(1:nr,kf) )
        !
        ! Get the face length for the current face
        !
        dsl = cl(this_geom) * norm2( geofa(nf)%v(1:nr,kf) )
        !
        ! Compute the BR2 jump correction for this point
        !    jump_left  =  half * (u_right - u_left)
        !    jump_right = -half * (u_right - u_left)
        !
        ! NOTE: Multiplying by the normal in the loop below that applies the
        !       correction
        !           i.e., correct_m * unrm(l)
        !       takes care of the minus sign if the current cell is the right
        !       cell on the current face.
        !
        ! NOTE: Boundary faces should never be the right cell on the interface.
        !       For non-wall boundaries, this should reduce to jump_left above.
        !       For wall boundaries, jump_left should reduce to:
        !         jump_left = (u_right - u_left) where u_right = u_wall
        !
        bias = zero
        !
        ! Reset bias if this is a wall boundary
        !
        if (nf <= nfbnd .and. walls_are_exact) then
          if (any(bface(1,nf) == bc_walls)) bias = half
        end if
        !
        hfpb = bias + half
        !
        ! Get the interface jumps to correct the gradients
        !
        do m = 1,nq
          ujump(m) = hfpb*(rufp(m) - lufp(m))
        end do
        !
        ! Apply the jump correction to the interior discontinuous gradients
        !
        do m = 1,nq
          correct_m = dsl * ujump(m)
          do l = 1,nr
            correct_n = correct_m * unrm(l)
            do i = cmat%row_ptr(lf)+1, cmat%row_ptr(lf+1)
              n = cmat%row_idx(i)
              a = cmat%row_val(i)
              face_dusp(n,l,m) = face_dusp(n,l,m) + a * correct_n
              cell_dusp(n,l,m) = cell_dusp(n,l,m) + a * correct_n
            end do
          end do
        end do
        !
      end do face_points
      !
      ! Interpolate the face corrected gradient from the interior solution
      ! points to the interface flux points now that we are done with all
      ! the flux points on the current face of the current cell.
      !
      interpolate_face: do kf = 1,geom_solpts(face_geom,face_order)
        !
        ! kf: index of flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        lf = cell(nc)%face(ncf)%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(cell(nc)%face(ncf)%rotation)% &
                                            pts(kf)
        !
       !facedusp_array(f0+kf,1:nr,1:nq,sf) = interp(this_geom,this_order)% &
       !                                            toFace(face_order)% &
       !                               dot_mat2( lf ,face_dusp(1:np,1:nr,1:nq) )
        do m = 1,nq
          do l = 1,nr
           !facedusp_array(f0+kf,l,m,sf) = interp(this_geom,this_order)% &
           !                                      toFace(face_order)% &
           !                            dot( lf ,face_dusp(1:np,l,m) )
            val = zero
            do i = interp(this_geom,this_order)% &
                             toFace(face_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)% &
                             toFace(face_order)%row_ptr(lf+1)
              !
              n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
              a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
              !
              val = val + a * face_dusp(n,l,m)
              !
            end do
            facedusp_array(f0+kf,l,m,sf) = val
          end do
        end do
        !
      end do interpolate_face
      !
    end do cell_faces
    !
    ! Copy the total corrected gradient accounting for all cell faces
    ! back into dusp
    !
    dusp(n1:n2,1:nr,1:nq) = cell_dusp(1:np,1:nr,1:nq)
    !
  end do cell_loop
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine correct_gradients_br2
!
!###############################################################################
!
subroutine correct_gradients_at_solpts()
  !
  !.. Use Statements ..
  use order_mod,         only : maxpts,geom_solpts
  use geovar,            only : nfbnd,nface,nr,n_solpts,ncell
  use geovar,            only : face,cell,bface,fp
  use metrics_mod,       only : geofa,metrics_dt
  use correction_mod,    only : correct
  use flowvar,           only : dusp,faceusp_array
  use ovar,              only : walls_are_exact
  use ovar,              only : LDG_beta
  !
  !.. Local Scalars ..
  integer  :: i,k,l,m,n,nk,nf,r1,r2,l1,l2,n1,n2,nc
  integer  :: kl,kr,nfp,f0
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: bias,hfpb,hfmb,cmat_val
  real(wp) :: correct_n,correct_m,ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)          :: unrm
  real(wp), dimension(1:nq)          :: ucom
  real(wp), dimension(1:nq)          :: ljfp,rjfp
  real(wp), dimension(1:nq,1:maxpts) :: lufp,rufp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "correct_gradients_at_solpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
  boundary_faces: do nf = 1,nfbnd
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0 ! number of flux points on this face
    !
    host_cell = face(nf)%left%cell ! left  cell on face
    !
    host_geom = cell(host_cell)%geom
    host_order = cell(host_cell)%order
    !
    l1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in left cell
    l2 = cell(host_cell)%end_sp ! ending    index for sol-pts in left cell
    !
    ! For all the flux points on the current face, extract the conservative
    ! conservative variables for the left and right cells from faceusp
    !
    do k = 1,nfp
      do m = 1,nq
        lufp(m,k) = faceusp_array(f0+k,m,1)
        rufp(m,k) = faceusp_array(f0+k,m,2)
      end do
    end do
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! Use the LDG flux formulation for computing the common interface solutions
    !   u_com = half*(u_left + u_right) - &
    !                (u_left - u_right)*dot(LDG_beta,unit_normal)
    !
    do k = 1,nfp
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kl = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(nf)%left%rotation)% &
                                     pts(k)
      !
      ! Get the unit normal at this flux point
      !
      unrm(1:nr) = unit_vector( geofa(nf)%v(1:nr,k) )
      !
      ! Get the face length for the current face
      !
      dsl = cl(host_geom) * norm2( geofa(nf)%v(1:nr,k) )
      !
      ! Compute the LDG directional parameters
      ! NOTE: bias is used for all boundary faces that allow for flow
      !       across the interface, e.g, inflow/outflow, cpu/partition
      !       boundary, or periodic boundary. These boundaries are
      !       treated the same as an interior face except that we dont
      !       need to worry about correcting the right cell.
      !
      bias = dot( LDG_beta(1:nr) , unrm(1:nr) )
      !
      ! Reset bias if this is a wall boundary
      !
      if (walls_are_exact) then
        if (any(bface(1,nf) == bc_walls)) bias = half
      end if
      !
      ! Use the merge function to overide bias with 1/2 if this is a
      ! wall boundary face meaning there is no flow across the interface,
      ! e.g., slip wall, symmetry boundary, no-slip wall. For these
      ! boundaries, the solution in rufp is the exact solution at the boundary
      ! and we use this exact solution for the interface common solution.
      ! This means the left jump is simply:
      !           jump_left = u_wall - u_interior
      !
      hfpb = half + bias
      !
      ! Compute the interface jumps for the left and right solutions
      ! NOTE: The common interface solution simplifies to
      !         u_com = (half - bias)*u_left + (half + bias)*u_right
      !       Using jump_left = u_com - u_left
      !         jump_left  = (half + bias)*(u_right-u_left)
      !
      do m = 1,nq
        ljfp(m) = hfpb*(rufp(m,k) - lufp(m,k))
      end do
      !
      ! Unsimplified version of the above LDG flux formulation
      !
     !if (any(bface(1,nf) == bc_walls)) then
     !  ucom(1:nq) = rufp(1:nq,k)
     !else
     !  do m = 1,nq
     !    ucom(m) = lufp(m,k)*(half - dot(LDG_beta(1:nr),+unrm(1:nr))) + &
     !              rufp(m,k)*(half - dot(LDG_beta(1:nr),-unrm(1:nr)))
     !  end do
     !end if
     !!
     !do m = 1,nq
     !  ljfp(m) = ucom(m) - lufp(m,k)
     !end do
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      do i = correct(host_geom,host_order)%row_ptr(kl)+1, &
             correct(host_geom,host_order)%row_ptr(kl+1)
        nk = correct(host_geom,host_order)%row_idx(i)
        n = l1 - 1 + nk
        correct_n = dsl * correct(host_geom,host_order)%row_val(i) &
                        / metrics_dt(host_cell)%jac(nk)
        do m = 1,nq
          correct_m = correct_n * ljfp(m)
          do l = 1,nr
            dusp(n,l,m) = dusp(n,l,m) + correct_m * unrm(l)
          end do
        end do
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
  interior_faces: do nf = nfbnd+1,nface
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0 ! number of flux points on this face
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
    l1 = cell(left_cell)%beg_sp ! beginning index for sol-pts in left cell
    l2 = cell(left_cell)%end_sp ! ending    index for sol-pts in left cell
    !
    r1 = cell(right_cell)%beg_sp ! beginning index for sol-pts in right cell
    r2 = cell(right_cell)%end_sp ! ending    index for sol-pts in right cell
    !
    ! For all the flux points on the current face, extract the conservative
    ! conservative variables for the left and right cells from faceusp
    !
    do k = 1,nfp
      do m = 1,nq
        lufp(m,k) = faceusp_array(f0+k,m,1)
        rufp(m,k) = faceusp_array(f0+k,m,2)
      end do
    end do
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! Use the LDG flux formulation for computing the common interface solutions
    !   u_com = half*(u_left + u_right) - &
    !                (u_left - u_right)*dot(LDG_beta,unit_normal)
    !
    do k = 1,nfp
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kl = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(nf)%left%rotation)% &
                                     pts(k)
      kr = face(nf)%right%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(nf)%right%rotation)% &
                                     pts(k)
      !
      ! Get the unit normal at this flux point
      !
      unrm(1:nr) = unit_vector( geofa(nf)%v(1:nr,k) )
      !
      ! Get the face length for the current face
      !
      ds = norm2( geofa(nf)%v(1:nr,k) )
      dsl = cl(left_geom)*ds
      dsr = cl(right_geom)*ds
      !
      ! Compute the LDG directional parameters
      !
      bias = dot( LDG_beta(1:nr) , unrm(1:nr) )
      !
      hfpb = half + bias
      hfmb = half - bias
      !
      ! Compute the interface jumps for the left and right solutions
      ! NOTE: The common interface solution simplifies to
      !         u_com = (half - bias)*u_left + (half + bias)*u_right
      !       Using jump_left = u_com - u_left and jump_right = u_com - u_right
      !         jump_left  = (half + bias)*(u_right-u_left)
      !         jump_right = (half - bias)*(u_left-u_right)
      !
      do m = 1,nq
        ljfp(m) = hfpb*(rufp(m,k) - lufp(m,k))
        rjfp(m) = hfmb*(lufp(m,k) - rufp(m,k))
      end do
      !
      ! Unsimplified version of the above LDG flux formulation
      !
     !do m = 1,nq
     !  ucom(m) = lufp(m,k)*(half - dot(LDG_beta(1:nr),+unrm(1:nr))) + &
     !            rufp(m,k)*(half - dot(LDG_beta(1:nr),-unrm(1:nr)))
     !end do
     !!
     !do m = 1,nq
     !  ljfp(m) = ucom(m) - lufp(m,k)
     !  rjfp(m) = ucom(m) - rufp(m,k)
     !end do
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      do i = correct(left_geom,left_order)%row_ptr(kl)+1, &
             correct(left_geom,left_order)%row_ptr(kl+1)
        nk = correct(left_geom,left_order)%row_idx(i)
        n = l1 - 1 + nk
        correct_n = dsl * correct(left_geom,left_order)%row_val(i) &
                        / metrics_dt(left_cell)%jac(nk)
        do m = 1,nq
          correct_m = correct_n * ljfp(m)
          do l = 1,nr
            dusp(n,l,m) = dusp(n,l,m) + correct_m * unrm(l)
          end do
        end do
      end do
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the right cell
      !
      do i = correct(right_geom,right_order)%row_ptr(kr)+1, &
             correct(right_geom,right_order)%row_ptr(kr+1)
        nk = correct(right_geom,right_order)%row_idx(i)
        n = r1 - 1 + nk
        correct_n = dsr * correct(right_geom,right_order)%row_val(i) &
                        / metrics_dt(right_cell)%jac(nk)
        do m = 1,nq
          correct_m = correct_n * rjfp(m)
          do l = 1,nr
            dusp(n,l,m) = dusp(n,l,m) - correct_m * unrm(l)
          end do
        end do
      end do
      !
    end do
    !
  end do interior_faces
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine correct_gradients_at_solpts
!
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
! Procedures for inviscid fluxes
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
