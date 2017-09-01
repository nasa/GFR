module flux_mod
  !
  !.. Use Statements ..
  use module_kind_types
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
  public :: flux_gradient_at_solpts_visc_gpu1
  public :: correct_interface_flux_visc
  public :: discontinuous_gradient_at_solpts
  public :: correct_gradients_at_solpts
  public :: correct_gradients_br2
  !
contains
!
!###############################################################################
!
subroutine interp_cv_to_flxpts()
  !
  use order_mod,         only : maxSP,maxFP,geom_solpts
  use geovar,            only : ncell,cell,fp
  use interpolation_mod, only : interp
  use flowvar,           only : usp,faceusp
  !
 !use parallel_mod,      only : exchange_faceusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,m
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  real(wp) :: val
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxSP,1:nq) :: tusp
  integer, dimension(1:maxFP) :: fpidx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_flxpts"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL DO GANG VECTOR_LENGTH(32) PRIVATE(tusp,fpidx)
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Create transpose of usp array to reduce cache misses
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do k = 1,np
      do m = 1,nq
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
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
      nfp = geom_solpts(face_geom,face_order)
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
      !$ACC LOOP VECTOR COLLAPSE(2)
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        do m = 1,nq
          lf = fpidx(kf)
          val = zero
          do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
                 interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
            j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
            val = val + tusp(j,m) * interp(this_geom,this_order)% &
                                              toFace(face_order)% &
                                    row_val(i)
          end do
          faceusp(nf)%v(m,kf,sf) = val
        end do
        !
      end do
      !
    end do
    !
  end do
  !
  ! Exchange the left and right states for faces on communication boundaries
  !
 !if (ncpu > 1) call exchange_faceusp(faceusp)
  !
end subroutine interp_cv_to_flxpts
!
!###############################################################################
!
subroutine interp_dusp_to_flxpts()
  !
  use order_mod,         only : maxSP,maxFP,geom_solpts
  use geovar,            only : nr,ncell,cell,fp
  use interpolation_mod, only : interp
  use flowvar,           only : dusp,facedusp
  !
 !use parallel_mod,      only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,l,m
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  real(wp) :: val
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tdusp
  integer, dimension(1:maxFP) :: fpidx
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL DO GANG VECTOR_LENGTH(32) PRIVATE(tdusp,fpidx)
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1               ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Create transpose of dusp array to reduce cache misses
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do k = 1,np
      do m = 1,nq
        do l = 1,nr
          tdusp(k,l,m) = dusp(l,m,k+n1-1)
        end do
      end do
    end do
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
      nfp = geom_solpts(face_geom,face_order)
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
      !$ACC LOOP VECTOR COLLAPSE(3)
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
       !lf = cell(nc)%face(n)%get_fp_idx(kf)
        !
        do m = 1,nq
          do l = 1,nr
            lf = fpidx(kf)
            val = zero
            do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
              j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
              val = val + tdusp(j,l,m) * interp(this_geom,this_order)% &
                                                   toFace(face_order)% &
                                         row_val(i)
            end do
            facedusp(nf)%v(l,m,kf,sf) = val
          end do
        end do
       !facedusp(nf)%v(1:nr,1:nq,kf,sf) = interp(this_geom,this_order)% &
       !                                            toFace(face_order)% &
       !                                  dot_mat2( lf , tdusp(1:np,1:nr,1:nq) )
      end do
      !
    end do
    !
  end do
  !
  ! Exchange the left and right states for faces on communication boundaries
  !
 !if (ncpu > 1) call exchange_facedusp(facedusp)
  !
end subroutine interp_dusp_to_flxpts
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
  use ovar, only : gam
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
  real(wp) :: ds,val,fval,correct_n
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
 !!$ACC DATA CREATE
  vflx(:,:) = zero
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  !$ACC PARALLEL DO GANG VECTOR_LENGTH(32) &
  !$ACC&                 PRIVATE(iflx,vflx,lcl_res,rsflx,xyflx,fpidx,flx)
 !!$ACC&                 PRIVATE(lcl_res,rsflx,xyflx,fpidx,flx)
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
   !!$ACC LOOP VECTOR COLLAPSE(2)
   !!$ACC LOOP WORKER
    !$ACC LOOP SEQ
    do k = 1,np
      !
      n = k + n1-1
      !
      ! Inviscid fluxes
      !
     !!$ACC LOOP SEQ
      do l = 1,nr
        call invs_flx_cv( usp(:,n) , l , iflx(:,l) ) ! l-directional flux
      end do
      !
      ! Viscous fluxes
      !
      if (governing_equations == NavierStokes_Eqns) then
        call visc_flx_cv( usp(:,n) , dusp(:,:,n) , vflx )
      end if
      !
      ! Total fluxes
      !
     !!$ACC LOOP SEQ
      do l = 1,nr
       !!$ACC LOOP SEQ
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
    !$ACC LOOP VECTOR COLLAPSE(3)
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
   !!$ACC LOOP VECTOR COLLAPSE(2)
   !do n = 1,maxSP
   !  do m = 1,nq
   !    lcl_res(m,n) = zero
   !  end do
   !end do
    !
    ! Loop through each coordinate direction and add the divergence of
    ! the flux to the local residual array for the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do k = 1,np
      do m = 1,nq
        val = zero
        do l = 1,nr
          do i = deriv(this_geom,this_order)%d(l)%row_ptr(k)+1, &
                 deriv(this_geom,this_order)%d(l)%row_ptr(k+1)
            n = deriv(this_geom,this_order)%d(l)%row_idx(i)
            val = val + rsflx(n,m,l) * deriv(this_geom,this_order)% &
                                                              d(l)% &
                                       row_val(i)
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
        residual(m,n) = lcl_res(m,n-n1+1)
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
      !$ACC LOOP VECTOR COLLAPSE(2)
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
            do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
              n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
              val = val + xyflx(n,m,l) * interp(this_geom,this_order)% &
                                                   toFace(face_order)% &
                                         row_val(i)
            end do
            fval = fval + val*geofa(nf)%v(l,kf)
          end do
          do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
                 correct(this_geom,this_order)%row_ptr(lf+1)
            n = n1 - 1 + correct(this_geom,this_order)%row_idx(i)
            correct_n = correct(this_geom,this_order)%row_val(i) * fval * ds
            residual(m,n) = residual(m,n) + correct_n
          end do
        end do
        !
      end do
      !
    end do
    !
  end do
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
  use geovar,            only : nr,ncell,using_quadrature
  use geovar,            only : cell,fp,n_solpts
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  use metrics_mod,       only : metrics_dt,geofa
  use flowvar,           only : usp,dusp,residual,faceusp
  use ovar,              only : governing_equations
  use ovar,              only : output_dir,dump_fluxes
  use ovar,              only : gam,rgasref_nd,cpref_nd,Pr
  use correction_mod,    only : correct
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
  integer :: ierr,io,meq
  integer :: fp_offset,face_rot,nfp,ncf,num_cell_faces
  real(wp) :: val,pres_sp,flxvel,q2
  real(wp) :: mvsc,mcdt,ave_ext_rate
  real(wp) :: cv_inv,rho_inv,kflx
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
  real(wp), dimension(1:maxSP,1:nq) :: tusp
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: tdusp
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: dpvdx
  real(wp), dimension(1:nq) :: var
  real(wp), dimension(1:nq,1:nr,1:n_solpts) :: fluxes
  real(wp), dimension(1:nq,1:nr) :: vflux
  real(wp), dimension(1:nq,1:nr) :: iflux
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tflx
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: rflx
  real(wp), dimension(1:maxFP) :: area
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "flux_gradient_at_solpts_visc_gpu"
  !
continue
  !
  vflux(1:nq,1:nr) = zero
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) PRIVATE(iflux,vflux)
  do n = 1,n_solpts
    !
   !!$ACC LOOP SEQ
    do l = 1,nr
      call invs_flx_cv( usp(:,n) , l , iflux(:,l) )
    end do
    !
    if (governing_equations == NavierStokes_Eqns) then
      call visc_flx_cv( usp(:,n) , dusp(:,:,n) , vflux )
    end if
    !
   !!$ACC LOOP SEQ
    do l = 1,nr
     !!$ACC LOOP SEQ
      do m = 1,nq
        fluxes(m,l,n) = iflux(m,l) - vflux(m,l)
      end do
    end do
    !
  end do
  !
  !#############################################################################
  !
  !$ACC PARALLEL DO GANG VECTOR_LENGTH(32) &
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
          tflx(k,l,m) = fluxes(m,l,n)
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
            val = val + tflx(k,i,m)*metrics_dt(nc)%met(k,l,i)
          end do
          rflx(k,l,m) = val * metrics_dt(nc)%jac(k)
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
            n = deriv(this_geom,this_order)%d(l)%row_idx(i)
            val = val + rflx(n,l,m) * deriv(this_geom,this_order)% &
                                                             d(l)% &
                                      row_val(i)
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
        residual(m,n) = lcl_res(m,n-n1+1)
      end do
    end do
    !
 !end do
 !!
 !!#############################################################################
 !!
 !!$ACC PARALLEL DO GANG VECTOR_LENGTH(32) PRIVATE(tflx,fpidx,area)
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
            do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
              n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
              val = val + tflx(n,l,m) * interp(this_geom,this_order)% &
                                                  toFace(face_order)% &
                                        row_val(i)
            end do
            fval = fval + val * geofa(nf)%v(l,kf)
          end do
          do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
                 correct(this_geom,this_order)%row_ptr(lf+1)
            n = n1 - 1 + correct(this_geom,this_order)%row_idx(i)
            correct_n = correct(this_geom,this_order)%row_val(i) * fval * ds
            residual(m,n) = residual(m,n) + correct_n
          end do
        end do
      end do
      !
    end do
    !
  end do
 !write (iout,*) "finished step #3"
 !!$ACC END DATA
  !
!  do nc = 1,ncell
!    !
!    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
!    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
!    np = n2-n1+1         ! number of sol-pts in this cell
!    !
!    this_geom = cell(nc)%geom   ! geometry type of this cell
!    this_order = cell(nc)%order ! solution order of this cell
!    !
!#ifdef USE_NEW_FLUX_METHOD
!    do n = n1,n2
!      pres(n-n1+1) = pressure_cv_sp( usp(1:nq,n) )
!    end do
!    !
!   !!$ACC LOOP VECTOR COLLAPSE(2)
!    do k = 1,np
!      do m = 1,nq
!        tusp(k,m) = usp(m,k+n1+1)
!      end do
!    end do
!    !
!   !!$ACC LOOP VECTOR COLLAPSE(3)
!    do k = 1,np
!      do m = 1,nq
!        do l = 1,nr
!          tdusp(k,m,l) = dusp(l,m,k+n1+1)
!        end do
!      end do
!    end do
!    !
!    ! Get the flux at each quadrature point in the current cell
!    !
!    ! Inviscid fluxes
!    !
!   !!$ACC LOOP VECTOR COLLAPSE(3)
!    do l = 1,nr
!      do m = 1,nq
!        do k = 1,np
!          meq = nec + l
!          xyflx(k,m,l) = tusp(k,m) * tusp(k,meq) / tusp(k,nec)
!        end do
!      end do
!    end do
!    !
!   !!$ACC LOOP VECTOR COLLAPSE(2)
!    do l = 1,nr
!      do k = 1,np
!        meq = nec + l
!        xyflx(k,meq,l) = xyflx(k,meq,l) + pres(k)
!      end do
!    end do
!    !
!   !!$ACC LOOP VECTOR COLLAPSE(2)
!    do l = 1,nr
!      do k = 1,np
!        meq = nec + l
!        xyflx(k,nee,l) = xyflx(k,nee,l) + pres(k)*tusp(k,meq)/tusp(k,nec)
!      end do
!    end do
!    !
!    ! Viscous fluxes
!    !
!    if (governing_equations == NavierStokes_Eqns) then
!      !
!        do l = 1,nr
!          do m = 1,nq
!            do k = 1,np
!              vflx(k,m,l) = zero
!            end do
!          end do
!        end do
!        !
!        ! Get the gradients of the primitive variables using the 
!        ! gradients of the conservative variables and the product rule
!        !
!       !call grad_cv_to_grad_pv_sp_sub(usp(1:nq,n),dcvdx,dpvdx)
!        rho_inv = one/usp(nec,n)
!        !
!        do m = 1,nq
!          var(m) = rho_inv*usp(m,n)
!        end do
!        !
!        ! Copy the density gradient
!        !
!        do l = 1,nr
!          dpvdx(l,nec) = dcvdx(l,nec)
!        end do
!        !
!        ! Using the chain rule:
!        !
!        ! 1. Compute the velocity gradients (m = 2,nq-1)
!        ! 2. Compute the total energy gradient and store in the
!        !    location for the temperature gradient (m = nq)
!        !
!        do m = nmb,nee
!          do l = 1,nr
!            dpvdx(l,m) = rho_inv*dcvdx(l,m) - rho_inv*var(m)*dcvdx(l,nec)
!          end do
!        end do
!        !
!        ! The location for the temperature gradient currently contains
!        ! the total energy gradient from the previous loop.
!        ! Now use the velocity gradients to subtract the gradient of
!        ! kinetic-energy-per-unit-mass from the total energy gradient.
!        ! 
!        do m = nmb,nme
!          do l = 1,nr
!            dpvdx(l,nee) = dpvdx(l,nee) - var(m)*dpvdx(l,m)
!          end do
!        end do
!        !
!        ! Finish computing the temperature gradient by dividing through by
!        ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
!        !
!        do l = 1,nr
!          dpvdx(l,nee) = dpvdx(l,nee) * cv_inv
!        end do
!        !### FINISHED COMPUTING GRADIENTS OF PRIMITIVE VARIABLES ###
!        !
!        ! Compute the molecular viscosity and molecular conductivity
!        !
!        mvsc = viscosity_cv_sp( usp(1:nq,n) )
!        mcdt = cpref_nd * mvsc / Pr
!        !
!        ! Compute the average extension rate
!        !
!        ave_ext_rate = zero
!        do l = 1,nr
!          ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
!        end do
!        ave_ext_rate = two3 * ave_ext_rate
!        !
!        ! Compute the viscous fluxes for the momentum equations which are just
!        ! the components of the viscous stress tensor.
!        ! NOTE: We are initially excluding the molecular viscosity term which
!        !       will be multiplied to the whole tau array after these loops.
!        !
!        do j = 1,nr
!          do i = 1,nr
!              vflx(k,nmb-1+i,j) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
!          end do
!        end do
!        !
!        ! Subtract the average extension rate from the normal stresses
!        !
!        do l = 1,nr
!          vflx(k,nmb-1+l,l) = vflx(k,nmb-1+l,l) - ave_ext_rate
!        end do
!        !
!        ! Multiply the viscous stresses by the molecular viscosity
!        !
!        do m = nmb,nme
!          do l = 1,nr
!            vflx(k,m,l) = mvsc * vflx(k,m,l)
!          end do
!        end do
!        !
!        ! Compute the viscous fluxes for the energy equation
!        !
!        do l = 1,nr
!          vflx(k,nq,l) = mcdt*dpvdx(l,nq) + dot( var(nmb:nme) , vflx(k,nmb:nme,l) )
!        end do
!      !
!      ! Subtract the viscous flux from the inviscid flux
!      !
!      do m = 1,nq
!        do l = 1,nr
!          do k = 1,np
!            xyflx(k,m,l) = xyflx(k,m,l) - vflx(k,m,l)
!          end do
!        end do
!      end do
!      !
!    end if
!    !
!    ! Transform the fluxes in the x-y physical
!    ! domain to the r-s computational domain
!    !
!   !!$ACC LOOP VECTOR COLLAPSE(3)
!    do l = 1,nr
!      do m = 1,nq
!        do k = 1,np
!          val = zero
!          do i = 1,nr
!            val = val + xyflx(k,m,i)*metrics_dt(nc)%met(k,l,i)
!          end do
!          rsflx(k,m,l) = val * metrics_dt(nc)%jack(k)
!        end do
!      end do
!    end do
!#else
!    !
!#define USE_INLINED_INVISCID_FLUXES
!    !
!    ! Inviscid fluxes
!    !
!   !!$ACC LOOP SEQ
!    do k = 1,np
!      do l = 1,nr
!        !
!        n = k + n1-1
!        !
!#ifdef USE_INLINED_FLUXES
!        meq = nec + l
!        !
!        flxvel = usp(meq,n)/usp(nec,n)
!        !
!        pres_sp = pressure_cv_sp( usp(1:nq,n) )
!        !
!        xyflx(k,nec,l) = usp(meq,n)
!        xyflx(k,nmb:nme,l) = flxvel * usp(nmb:nme,n)
!        xyflx(k,meq,l) = xyflx(k,meq,l) + pres_sp
!        xyflx(k,nee,l) = flxvel * (usp(nee,n) + pres_sp)
!#else
!        call invs_flx_cv( usp(1,n) , l , flx ) ! x-directional flux
!        xyflx(k,1:nq,l) = flx(1:nq)
!#endif
!        !
!      end do
!    end do
!    !
!    ! Viscous fluxes
!    !
!    if (governing_equations == NavierStokes_Eqns) then
!      do k = 1,np
!        !
!        n = k + n1-1
!        !
!#ifdef USE_INLINED_FLUXES
!        do l = 1,nr
!          do m = 1,nq
!            vflx(m,l) = zero
!          end do
!        end do
!        !
!        do m = 1,nq
!          do l = 1,nr
!            dcvdx(l,m) = dusp(l,m,n)
!          end do
!        end do
!        !
!        ! Get the gradients of the primitive variables using the 
!        ! gradients of the conservative variables and the product rule
!        !
!       !call grad_cv_to_grad_pv_sp_sub(usp(1:nq,n),dcvdx,dpvdx)
!        rho_inv = one/usp(nec,n)
!        !
!        do m = 1,nq
!          var(m) = rho_inv*usp(m,n)
!        end do
!        !
!        ! Copy the density gradient
!        !
!        do l = 1,nr
!          dpvdx(l,nec) = dcvdx(l,nec)
!        end do
!        !
!        ! Using the chain rule:
!        !
!        ! 1. Compute the velocity gradients (m = 2,nq-1)
!        ! 2. Compute the total energy gradient and store in the
!        !    location for the temperature gradient (m = nq)
!        !
!        do m = nmb,nee
!          do l = 1,nr
!            dpvdx(l,m) = rho_inv*dcvdx(l,m) - rho_inv*var(m)*dcvdx(l,nec)
!          end do
!        end do
!        !
!        ! The location for the temperature gradient currently contains
!        ! the total energy gradient from the previous loop.
!        ! Now use the velocity gradients to subtract the gradient of
!        ! kinetic-energy-per-unit-mass from the total energy gradient.
!        ! 
!        do m = nmb,nme
!          do l = 1,nr
!            dpvdx(l,nee) = dpvdx(l,nee) - var(m)*dpvdx(l,m)
!          end do
!        end do
!        !
!        ! Finish computing the temperature gradient by dividing through by
!        ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
!        !
!        do l = 1,nr
!          dpvdx(l,nee) = dpvdx(l,nee) * cv_inv
!        end do
!        !### FINISHED COMPUTING GRADIENTS OF PRIMITIVE VARIABLES ###
!        !
!        ! Compute the molecular viscosity and molecular conductivity
!        !
!        mvsc = viscosity_cv_sp( usp(1:nq,n) )
!        mcdt = cpref_nd * mvsc / Pr
!        !
!        ! Compute the average extension rate
!        !
!        ave_ext_rate = zero
!        do l = 1,nr
!          ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
!        end do
!        ave_ext_rate = two3 * ave_ext_rate
!        !
!        ! Compute the viscous fluxes for the momentum equations which are just
!        ! the components of the viscous stress tensor.
!        ! NOTE: We are initially excluding the molecular viscosity term which
!        !       will be multiplied to the whole tau array after these loops.
!        !
!        do j = 1,nr
!          do i = 1,nr
!              vflx(nmb-1+i,j) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
!          end do
!        end do
!        !
!        ! Subtract the average extension rate from the normal stresses
!        !
!        do l = 1,nr
!          vflx(nmb-1+l,l) = vflx(nmb-1+l,l) - ave_ext_rate
!        end do
!        !
!        ! Multiply the viscous stresses by the molecular viscosity
!        !
!        do m = nmb,nme
!          do l = 1,nr
!            vflx(m,l) = mvsc * vflx(m,l)
!          end do
!        end do
!        !
!        ! Compute the viscous fluxes for the energy equation
!        !
!        do l = 1,nr
!          vflx(nq,l) = mcdt*dpvdx(l,nq) + dot( var(nmb:nme) , vflx(nmb:nme,l) )
!        end do
!#else
!        call visc_flx_cv( usp(1,n) , dusp(1,1,n) , vflx )
!#endif
!        !
!        ! Subtract the viscous flux from the inviscid flux
!        !
!        do m = 1,nq
!          do l = 1,nr
!            xyflx(k,m,l) = xyflx(k,m,l) - vflx(m,l)
!          end do
!        end do
!        !
!      end do
!    end if
!    !
!    ! Transform the fluxes in the x-y physical
!    ! domain to the r-s computational domain
!    !
!   !!$ACC LOOP VECTOR COLLAPSE(3)
!    do m = 1,nq
!      do l = 1,nr
!        do k = 1,np
!          val = zero
!         !!$ACC LOOP SEQ
!          do i = 1,nr
!            val = val + xyflx(k,m,i)*metrics_dt(nc)%met(k,l,i)
!          end do
!          rsflx(k,m,l) = val * metrics_dt(nc)%jac(k)
!        end do
!      end do
!    end do
!#endif
!    !
!    ! Compute the divergence of the discontinuous fluxes 
!    !
!    ! Zero the residual array for the current cell
!    !
!    do n = 1,maxSP
!      do m = 1,nq
!        lcl_res(m,n) = zero
!      end do
!    end do
!    !
!    ! Loop through each coordinate direction and add the divergence of
!    ! the flux to the local residual array for the current cell
!    !
!   !!$ACC LOOP VECTOR COLLAPSE(2)
!    do k = 1,np
!      do m = 1,nq
!        val = zero
!        do l = 1,nr
!          do i = deriv(this_geom,this_order)%d(l)%row_ptr(k)+1, &
!                 deriv(this_geom,this_order)%d(l)%row_ptr(k+1)
!            n = deriv(this_geom,this_order)%d(l)%row_idx(i)
!            val = val + rsflx(n,m,l) * deriv(this_geom,this_order)% &
!                                                              d(l)% &
!                                       row_val(i)
!          end do
!        end do
!        lcl_res(m,k) = val
!      end do
!    end do
!    !
!    ! Copy the local residual for the current cell into the residual array
!    !
!   !!$ACC LOOP VECTOR COLLAPSE(2)
!    do n = n1,n2
!      do m = 1,nq
!        residual(m,n) = lcl_res(m,n-n1+1)
!      end do
!    end do
!    !
!    ! Interpolate the fluxes to the interface flux points and dot this
!    ! interpolated vector with the face normal at each flux point
!    ! NOTE: We dont need to exchange this information for communication and
!    !       periodic boundary faces since this flux is only used within the
!    !       current cell.
!    !
!    ! Loop over the faces of this cell
!    !
!    do ncf = 1,size(cell(nc)%face)
!      !
!      nf = cell(nc)%face(ncf)%idx  ! index of the current face
!      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
!      !
!      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
!      face_order = cell(nc)%face(ncf)%order ! order of current face
!      fp_offset = cell(nc)%face(ncf)%fp_offset
!      face_rot = cell(nc)%face(ncf)%rotation
!      !
!      nfp = geom_solpts(face_geom,face_order)
!      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
!                                    rot(face_rot)%pts(1:nfp)
!      !
!      do kf = 1,nfp
!        !
!        ! kf: index of the flx-pt local to this face
!        ! lf: index of flx-pt within interp%mat
!       !lf = cell(nc)%face(ncf)%get_fp_idx(kf)
!        !
!        ! Get the unit normal for this face flux/solution point
!        !
!        val = zero
!       !!$ACC LOOP SEQ
!        do l = 1,nr
!          val = val + geofa(nf)%v(l,kf) * geofa(nf)%v(l,kf)
!        end do
!        unrm(1:nr) = geofa(nf)%v(1:nr,kf) / sqrt(val)
!        !
!        do l = 1,nr
!          do m = 1,nq
!            lf = fpidx(kf)
!            val = zero
!           !!$ACC LOOP SEQ
!            do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
!                   interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
!              n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
!              val = val + xyflx(n,m,l) * interp(this_geom,this_order)% &
!                                                   toFace(face_order)% &
!                                         row_val(i)
!            end do
!            vflx(m,l) = val
!          end do
!        end do
!        !
!        do m = 1,nq
!          val = zero
!         !!$ACC LOOP SEQ
!          do l = 1,nr
!            val = val + vflx(m,l) * unrm(l)
!          end do
!        end do
!        !
!      end do
!      !
!    end do
!    !
!  end do
  !
end subroutine flux_gradient_at_solpts_visc_gpu
!
!###############################################################################
!
subroutine flux_gradient_at_solpts_visc_gpu1(time_step,rk_step)
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,maxFP,geom_solpts
  use geovar,            only : nr,ncell,using_quadrature
  use geovar,            only : cell,fp,n_solpts
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  use metrics_mod,       only : metrics_dt,geofa
  use flowvar,           only : usp,dusp,residual
  use ovar,              only : governing_equations
  use ovar,              only : output_dir,dump_fluxes
  use ovar,              only : gam,rgasref_nd,cpref_nd,Pr
  use correction_mod,    only : correct
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
  integer :: ierr,io
  integer :: fp_offset,face_rot,nfp,ncf,num_cell_faces
  real(wp) :: val,flxvel
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
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: tdusp
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: dpvdx
  real(wp), dimension(1:nq) :: var
  real(wp), dimension(1:nq,1:nr) :: vflux
  real(wp), dimension(1:nq,1:nr) :: iflux
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tflx
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: rflx
  real(wp), dimension(1:maxFP) :: area
  !
  real(wp), dimension(1:n_solpts,1:nr,1:nq) :: tfluxes
  real(wp), dimension(1:n_solpts,1:nr,1:nq) :: ifluxes
  real(wp), dimension(1:n_solpts,1:nr,1:nq) :: vfluxes
  real(wp), dimension(1:n_solpts,1:nq) :: tusp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "flux_gradient_at_solpts_visc_gpu1"
  !
continue
  !
 !vflux(1:nq,1:nr) = zero
 !!
 !!$ACC PARALLEL LOOP VECTOR_LENGTH(32) PRIVATE(iflux,vflux)
 !do n = 1,n_solpts
 !  !
 ! !!$ACC LOOP SEQ
 !  do l = 1,nr
 !    call invs_flx_cv( usp(:,n) , l , iflux(:,l) )
 !  end do
 !  !
 !  if (governing_equations == NavierStokes_Eqns) then
 !    call visc_flx_cv( usp(:,n) , dusp(:,:,n) , vflux )
 !  end if
 !  !
 ! !!$ACC LOOP SEQ
 !  do l = 1,nr
 !   !!$ACC LOOP SEQ
 !    do m = 1,nq
 !      fluxes(m,l,n) = iflux(m,l) - vflux(m,l)
 !    end do
 !  end do
 !  !
 !end do
  !
  call invs_flx_cv_gpu( usp , ifluxes )
  !
  if (governing_equations == NavierStokes_Eqns) then
    !
#define USE_GPU_FLUXES
#ifndef USE_GPU_FLUXES
   !!$ACC PARALLEL LOOP
    !$ACC PARALLEL LOOP PRIVATE(vflux)
    do n = 1,n_solpts
      call visc_flx_cv( usp(:,n) , dusp(:,:,n) , vflux )
      do l = 1,nr
        do m = 1,nq
          vfluxes(n,l,m) = vflux(m,l)
        end do
      end do
    end do
#else
    call visc_flx_cv_gpu( usp , dusp , vfluxes )
#endif
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
        tfluxes(n,l,m) = ifluxes(n,l,m) - vfluxes(n,l,m)
      end do
    end do
  end do
  !
  !#############################################################################
  !
  !$ACC PARALLEL DO GANG VECTOR_LENGTH(32) &
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
            val = val + tflx(k,i,m)*metrics_dt(nc)%met(k,l,i)
          end do
          rflx(k,l,m) = val * metrics_dt(nc)%jac(k)
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
            n = deriv(this_geom,this_order)%d(l)%row_idx(i)
            val = val + rflx(n,l,m) * deriv(this_geom,this_order)% &
                                                             d(l)% &
                                      row_val(i)
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
        residual(m,n) = lcl_res(m,n-n1+1)
      end do
    end do
    !
 !end do
 !!
 !!#############################################################################
 !!
 !!$ACC PARALLEL DO GANG VECTOR_LENGTH(32) PRIVATE(tflx,fpidx,area)
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
            do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
              n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
              val = val + tflx(n,l,m) * interp(this_geom,this_order)% &
                                                  toFace(face_order)% &
                                        row_val(i)
            end do
            fval = fval + val * geofa(nf)%v(l,kf)
          end do
          do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
                 correct(this_geom,this_order)%row_ptr(lf+1)
            n = n1 - 1 + correct(this_geom,this_order)%row_idx(i)
            correct_n = correct(this_geom,this_order)%row_val(i) * fval * ds
            residual(m,n) = residual(m,n) + correct_n
          end do
        end do
      end do
      !
    end do
    !
  end do
  !
end subroutine flux_gradient_at_solpts_visc_gpu1
!
!###############################################################################
!
subroutine correct_interface_flux_visc(time_step,rk_step)
  !
  use generic_types_mod, only : matrix
  use order_mod,         only : maxFP,geom_solpts
  use geovar,            only : nr,nfbnd,nface
  use geovar,            only : bface,face,cell
  use ovar,              only : invs_flux_method,governing_equations
  use ovar,              only : output_dir,dump_fluxes
  use metrics_mod,       only : geofa
  use correction_mod,    only : correct
  use flowvar,           only : usp,residual,faceusp,facedusp
  !
  !.. Formal Arguments ..
  integer, intent(in) :: time_step
  integer, intent(in) :: rk_step
  !
  !.. Local Scalars ..
  integer  :: i,k,l,m,n,nf,r1,r2,l1,l2
  integer  :: kl,kr,iu,nfp
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
    nfp = geom_solpts(face_geom,face_order) ! number of flux points on this face
    !
    host_cell  = face(nf)%left%cell    ! host cell on boundary face
    host_geom  = cell(host_cell)%geom  ! geometry of host cell
    host_order = cell(host_cell)%order ! order of host cell
    !
    l1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in left cell
    l2 = cell(host_cell)%end_sp ! ending    index for sol-pts in left cell
    !
    ! FP solution from left  cell on face
    luqp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,1)
    ! FP solution from right cell on face
    ruqp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,2)
    !
    if (governing_equations == NavierStokes_Eqns) then
      ! Solution gradient at FP from left  cell on face
      lduqp(1:nr,1:nq,1:nfp) = facedusp(nf)%v(1:nr,1:nq,1:nfp,1)
      ! Solution gradient at FP from right cell on face
      rduqp(1:nr,1:nq,1:nfp) = facedusp(nf)%v(1:nr,1:nq,1:nfp,2)
    end if
    !
    rnrm(1:nr,1:nfp) = geofa(nf)%v(1:nr,1:nfp) ! Area face normals at FP
    !
    ! Get the common flux at each quadrature flux point
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
      kl = face(nf)%left%get_fp_idx(face_geom,face_order,k)
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
          residual(m,n) = residual(m,n) + fncq(m,k) * cmat_val
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
    nfp = geom_solpts(face_geom,face_order) ! number of flux points on this face
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
    luqp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,1)
    ruqp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,2)
    !
    if (governing_equations == NavierStokes_Eqns) then
      lduqp(1:nr,1:nq,1:nfp) = facedusp(nf)%v(1:nr,1:nq,1:nfp,1)
      rduqp(1:nr,1:nq,1:nfp) = facedusp(nf)%v(1:nr,1:nq,1:nfp,2)
    end if
    !
    rnrm(1:nr,1:nfp) = geofa(nf)%v(1:nr,1:nfp)
    !
    ! Get the common flux at each quadrature flux point
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
      kl = face(nf)%left%get_fp_idx(face_geom,face_order,k)
      kr = face(nf)%right%get_fp_idx(face_geom,face_order,k)
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
          residual(m,n) = residual(m,n) + fncq(m,k) * cmat_val
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
          residual(m,n) = residual(m,n) - fncq(m,k) * cmat_val
        end do
      end do
      !
    end do
    !
  end do interior_fluxes
  !
end subroutine correct_interface_flux_visc
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
subroutine discontinuous_gradient_at_solpts()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use order_mod,       only : maxSP
  use geovar,          only : nr,ncell,cell
  use metrics_mod,     only : metrics_dt
  use derivatives_mod, only : deriv
  use flowvar,         only : usp,dusp
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,m,n,nc,n1,n2,np
  integer :: this_geom,this_order
  real(wp) :: val
  !
  !.. Local Arrays ..
 !real(wp), dimension(1:nr) :: dudr
  real(wp), dimension(1:nr,1:nq,1:maxSP) :: dudr
  real(wp), dimension(1:maxSP,1:nq) :: tusp
 !real(wp), dimension(1:maxSP,1:nr,1:2) :: dtusp
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  !$ACC PARALLEL DO GANG VECTOR_LENGTH(32) PRIVATE(dudr,tusp)
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index for solpts in cell nc
    n2 = cell(nc)%end_sp ! ending    index for solpts in cell nc
    np = n2-n1+1 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Store the transpose of usp for this cell to reduce cache misses
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do k = 1,np
      do m = 1,nq
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
   !do l = 1,nr
   !  dtusp(1:np,l,1) = deriv(this_geom,this_order)%d(l)%vm_mult(tusp(1:np,nec))
   !end do
   !!
   !do l = 1,nr
   !  dtusp(1:np,l,2) = matmul(transpose(deriv(this_geom,this_order)%d(l)%mat), &
   !                           tusp(1:np,nec))
   !end do
   !!
   !write (170,100) ((dtusp(k,l,1),k=1,np),l=1,nr)
   !write (171,100) ((dtusp(k,l,2),k=1,np),l=1,nr)
   !100 format (/,40("-"),/,3(16(4(es14.6),/),/))
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
            j = deriv(this_geom,this_order)%d(l)%row_idx(i)
            val = val + tusp(j,m) * deriv(this_geom,this_order)% &
                                                           d(l)% &
                                    row_val(i)
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
          n = k+n1-1 ! total index of solution point on current partition
          dusp(l,m,n) = dot( dudr(1:nr,m,k) , metrics_dt(nc)%met(k,1:nr,l) )
        end do
      end do
    end do
    !
  end do
  !
end subroutine discontinuous_gradient_at_solpts
!
!###############################################################################
!
subroutine correct_gradients_at_solpts()
  !
  !.. Use Statements ..
  use generic_types_mod, only : matrix
  use order_mod,         only : maxpts,geom_solpts
  use geovar,            only : nfbnd,nface,nr,n_solpts,ncell
  use geovar,            only : face,cell,bface
  use metrics_mod,       only : geofa,metrics_dt
  use correction_mod,    only : correct
  use flowvar,           only : dusp,faceusp
  use ovar,              only : walls_are_exact
  use ovar,              only : LDG_beta
  !
  !.. Local Scalars ..
  integer  :: i,k,l,m,n,nk,nf,r1,r2,l1,l2,n1,n2,nc
  integer  :: kl,kr,nfp
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: cbias,hfpb,hfmb,cmat_val
  real(wp) :: correct_n,correct_m,ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)          :: unrm
  real(wp), dimension(1:nq)          :: ucom
  real(wp), dimension(1:nq)          :: ljfp,rjfp
  real(wp), dimension(1:nq,1:maxpts) :: lufp,rufp
  !
continue
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
    nfp = geom_solpts(face_geom,face_order) ! number of flux points on this face
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
    lufp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,1)
    rufp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,2)
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
      kl = face(nf)%left%get_fp_idx(face_geom,face_order,k)
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
      ! NOTE: cbias is used for all boundary faces that allow for flow
      !       across the interface, e.g, inflow/outflow, cpu/partition
      !       boundary, or periodic boundary. These boundaries are
      !       treated the same as an interior face except that we dont
      !       need to worry about correcting the right cell.
      !
      cbias = dot( LDG_beta(1:nr) , unrm(1:nr) )
      !
      ! Reset cbias if this is a wall boundary
      !
      if (walls_are_exact) then
        if (any(bface(1,nf) == bc_walls)) cbias = half
      end if
      !
      ! Use the merge function to overide cbias with 1/2 if this is a
      ! wall boundary face meaning there is no flow across the interface,
      ! e.g., slip wall, symmetry boundary, no-slip wall. For these
      ! boundaries, the solution in rufp is the exact solution at the boundary
      ! and we use this exact solution for the interface common solution.
      ! This means the left jump is simply:
      !           jump_left = u_wall - u_interior
      !
      hfpb = half + cbias
      !
      ! Compute the interface jumps for the left and right solutions
      ! NOTE: The common interface solution simplifies to
      !         u_com = (half - cbias)*u_left + (half + cbias)*u_right
      !       Using jump_left = u_com - u_left
      !         jump_left  = (half + cbias)*(u_right-u_left)
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
            dusp(l,m,n) = dusp(l,m,n) + correct_m * unrm(l)
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
    nfp = geom_solpts(face_geom,face_order) ! number of flux points on this face
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
    lufp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,1)
    rufp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,2)
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
      kl = face(nf)%left%get_fp_idx(face_geom,face_order,k)
      kr = face(nf)%right%get_fp_idx(face_geom,face_order,k)
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
      cbias = dot( LDG_beta(1:nr) , unrm(1:nr) )
      !
      hfpb = half + cbias
      hfmb = half - cbias
      !
      ! Compute the interface jumps for the left and right solutions
      ! NOTE: The common interface solution simplifies to
      !         u_com = (half - cbias)*u_left + (half + cbias)*u_right
      !       Using jump_left = u_com - u_left and jump_right = u_com - u_right
      !         jump_left  = (half + cbias)*(u_right-u_left)
      !         jump_right = (half - cbias)*(u_left-u_right)
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
            dusp(l,m,n) = dusp(l,m,n) + correct_m * unrm(l)
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
            dusp(l,m,n) = dusp(l,m,n) - correct_m * unrm(l)
          end do
        end do
      end do
      !
    end do
    !
  end do interior_faces
  !
end subroutine correct_gradients_at_solpts
!
!###############################################################################
!
subroutine correct_gradients_br2()
  !
  use order_mod,         only : maxSP,maxFP,geom_solpts
  use geovar,            only : cell
  use geovar,            only : nr,ncell,nfbnd,bface
  use geovar,            only : fp
  use metrics_mod,       only : geofa,metrics_dt
  use correction_mod,    only : correct
  use flowvar,           only : dusp,facedusp,faceusp
  use ovar,              only : walls_are_exact
  use interpolation_mod, only : interp
  !
 !use parallel_mod,      only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,l,m
  integer :: nf,kf,lf,sf,ncf,j
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: dsl,cbias,hfpb,correct_n,correct_m,val,face_sum
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxSP,1:nr,1:nq,1:6) :: face_dusp_corr
  !
  integer :: i,nfp
  integer :: array(1:1) ! DONT REMOVE THIS
  !
  integer, dimension(1:6) :: cf_nf,cf_sf,cf_face_order,cf_nfp
  real(wp), dimension(1:6) :: cf_hfpb
  integer, dimension(1:maxFP,1:6) :: cf_fpidx
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL DO GANG VECTOR_LENGTH(32) PRIVATE(face_dusp_corr,cf_fpidx, &
  !$ACC&                                           cf_nf,cf_sf,cf_nfp,cf_hfpb, &
  !$ACC&                                           cf_face_order,array)
  cell_loop: do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
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
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      cf_nf(ncf) = nf
      cf_sf(ncf) = sf
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      cf_face_order(ncf) = face_order
      nfp = geom_solpts(face_geom,face_order)
      cf_nfp(ncf) = nfp
      cf_fpidx(1:nfp,ncf) = fp_offset + fp%geord(face_geom,face_order)% &
                                           rot(face_rot)%pts(1:nfp)
      !
      ! Reset cbias if this is a wall boundary
      !
      cbias = zero
      if (nf <= nfbnd .and. walls_are_exact) then
        !$ACC LOOP SEQ
        do m = 1,size(bc_walls)
          if (bface(1,nf) == bc_walls(m)) then
            cbias = half
            exit
          end if
        end do
      end if
      !
      cf_hfpb(ncf) = ( cbias + half ) * cl(this_geom)
      !
    end do 
    !
    face_dusp_corr(:,:,:,:) = zero
   !!$ACC LOOP VECTOR COLLAPSE(4)
   !do ncf = 1,num_cell_faces
   !  do m = 1,nq
   !    do l = 1,nr
   !      do n = 1,np
   !        face_dusp_corr(n,l,m,ncf) = zero
   !      end do
   !    end do
   !  end do
   !end do
    !
    ! Loop over the faces of this cell
#ifndef ENABLE_LOOP
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
            nf = cf_nf(ncf)
           !nf = cell(nc)%face(ncf)%idx
            lf = cf_fpidx(kf,ncf)
            do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
                   correct(this_geom,this_order)%row_ptr(lf+1)
              n = correct(this_geom,this_order)%row_idx(i)
              correct_n = faceusp(nf)%v(m,kf,2) - faceusp(nf)%v(m,kf,1)
              correct_m = correct(this_geom,this_order)%row_val(i) * &
                          correct_n * cf_hfpb(ncf) / metrics_dt(nc)%jac(n)
              face_dusp_corr(n,l,m,ncf) = face_dusp_corr(n,l,m,ncf) + &
                                          correct_m * geofa(nf)%v(l,kf)
            end do
          end do
        end do
        !
      end do face_points
    end do cell_faces
#endif
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
            nf = cf_nf(ncf)
           !nf = cell(nc)%face(ncf)%idx
            sf = cf_sf(ncf)
           !sf = cell(nc)%face(ncf)%side
            face_order = cf_face_order(ncf)
            lf = cf_fpidx(kf,ncf)
            val = zero
            do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
              n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
              val = val + ( dusp(l,m,n-1+n1) + face_dusp_corr(n,l,m,ncf) ) &
                                      * interp(this_geom,this_order)% &
                                                   toFace(face_order)% &
                                        row_val(i)
            end do
            facedusp(nf)%v(l,m,kf,sf) = val
          end do
        end do
        !
      end do interpolate_face
      !
    end do cell_faces2
    !
    ! Copy the total corrected gradient accounting for all cell faces
    ! back into dusp
    !
    !$ACC LOOP VECTOR COLLAPSE(3)
    do n = n1,n2
      do m = 1,nq
        do l = 1,nr
          do ncf = 1,num_cell_faces
            dusp(l,m,n) = dusp(l,m,n) + face_dusp_corr(n-n1+1,l,m,ncf)
          end do
        end do
      end do
    end do
    !
   !!$ACC LOOP VECTOR COLLAPSE(3)
   !do n = n1,n2
   !  do i = 1,nq
   !    do j = 1,nr
   !      face_sum = zero
   !      do ncf = 1,num_cell_faces
   !        face_sum = face_sum + face_dusp_corr(n-n1+1,j,i,ncf)
   !      end do
   !      dusp(j,i,n) = dusp(j,i,n) + face_sum
   !    end do
   !  end do
   !end do
    !
  end do cell_loop
  !
  ! Exchange the left and right states for faces on communication boundaries
  !
 !if (ncpu > 1) call exchange_facedusp(facedusp)
  !
end subroutine correct_gradients_br2
!
!###############################################################################
!
!subroutine correct_gradients_br2()
! !
! use generic_types_mod, only : matrix
! use order_mod,         only : maxSP,maxFP,geom_solpts
! use geovar,            only : cell
! use geovar,            only : nr,ncell,nfbnd,bface
! use geovar,            only : face,fp
! use metrics_mod,       only : geofa,metrics_dt
! use correction_mod,    only : correct
! use flowvar,           only : dusp,facedusp,faceusp
! use ovar,              only : walls_are_exact,LDG_beta
! use interpolation_mod, only : interp
! !
!!use parallel_mod,      only : exchange_facedusp
! !
! !.. Local Scalars ..
! integer :: n,nc,n1,n2,np,k,l,m
! integer :: nf,kf,lf,sf,ncf
! integer :: this_geom,this_order
! integer :: face_geom,face_order
! integer :: fp_offset,face_rot
! real(wp) :: dsl,cbias,hfpb,correct_n,correct_m,val
! !
! !.. Local Arrays ..
! real(wp), dimension(1:nr) :: unrm
! real(wp), dimension(1:nq) :: lufp,rufp,ujump
! real(wp), dimension(1:nr,1:nq,1:maxSP) :: cell_dusp
! real(wp), dimension(1:maxSP,1:nr,1:nq) :: face_dusp
! real(wp), dimension(1:maxSP,1:nr,1:nq) :: dusp_copy
! !
! integer :: i,i1,i2,nfp,j,j1,j2,r1,r2,r3
! type(matrix) :: cmat
! integer :: array(1:1) ! DONT REMOVE THIS
! integer :: fpidx(1:maxFP)
! real(wp) :: fp_unrm(1:nr,1:maxFP)
! real(wp) :: fp_area(1:maxFP)
! !
! type :: matrix2
!   real(wp), allocatable :: row_val(:)
!   integer, allocatable :: row_idx(:)
!   integer, allocatable :: row_ptr(:)
! end type matrix2
! !
! type(matrix2), allocatable :: correct2(:,:)
! type :: interp2_t
!   type(matrix2), allocatable :: toFace(:)
! end type interp2_t
! type(interp2_t), allocatable :: interp2(:,:)
!continue
! !
! ! Loop through each interior cell and interpolate the conserved
! ! variables at the solution points to the flux points of each cell
! !
!!n1 = lbound(correct,dim=2)
!!n2 = ubound(correct,dim=2)
!!i1 = lbound(correct,dim=1)
!!i2 = ubound(correct,dim=1)
!!allocate (correct2(i1:i2,n1:n2))
!!do n = lbound(correct,dim=2),ubound(correct,dim=2)
!!  do i = lbound(correct,dim=1),ubound(correct,dim=1)
!!    if (allocated(correct(i,n)%row_val)) &
!!      correct2(i,n)%row_val = correct(i,n)%row_val
!!    if (allocated(correct(i,n)%row_idx)) &
!!      correct2(i,n)%row_idx = correct(i,n)%row_idx
!!    if (allocated(correct(i,n)%row_ptr)) &
!!      correct2(i,n)%row_ptr = correct(i,n)%row_ptr
!!  end do
!!end do
!!n1 = lbound(interp,dim=2)
!!n2 = ubound(interp,dim=2)
!!i1 = lbound(interp,dim=1)
!!i2 = ubound(interp,dim=1)
!!allocate (interp2(i1:i2,n1:n2))
!!do n = n1,n2
!!  do i = i1,i2
!!    if (allocated(interp(i,n)%toFace)) then
!!      j1 = lbound(interp(i,n)%toFace,dim=1)
!!      j2 = ubound(interp(i,n)%toFace,dim=1)
!!      allocate ( interp2(i,n)%toFace(j1:j2) )
!!      do j = j1,j2
!!        if (allocated(interp(i,n)%toFace(j)%row_val)) &
!!          interp2(i,n)%toFace(j)%row_val = interp(i,n)%toFace(j)%row_val
!!        if (allocated(interp(i,n)%toFace(j)%row_idx)) &
!!          interp2(i,n)%toFace(j)%row_idx = interp(i,n)%toFace(j)%row_idx
!!        if (allocated(interp(i,n)%toFace(j)%row_ptr)) &
!!          interp2(i,n)%toFace(j)%row_ptr = interp(i,n)%toFace(j)%row_ptr
!!      end do
!!    end if
!!  end do
!!end do
! !
!!!$ACC PARALLEL DO GANG VECTOR_LENGTH(64) PRIVATE(unrm,lufp,rufp, &
!!!$ACC PARALLEL DO GANG PRIVATE(unrm,lufp,rufp, &
! !$ACC PARALLEL DO GANG VECTOR_LENGTH(32) PRIVATE(unrm,lufp,rufp, &
! !$ACC&                    fpidx,fp_unrm,fp_area, &
! !$ACC&                    ujump,cell_dusp,face_dusp,dusp_copy,array)
! cell_loop: do nc = 1,ncell
!   !
!   n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
!   n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
!   np = n2-n1+1         ! number of sol-pts in this cell
!   !
!   this_geom = cell(nc)%geom   ! geometry type of this cell
!   this_order = cell(nc)%order ! solution order of this cell
!   !
!   ! Create a local copy of dusp for this cell to which all the
!   ! corrections will be applied. This will help reduce cache misses.
!   !
!   ! Create a transposed copy of the uncorrected dusp for this cell which
!   ! will be used to reset the face_dusp array in the cell_faces loop. This
!   ! will help reduce cache misses.
!   !
!   !$ACC LOOP VECTOR COLLAPSE(3)
!   do n = 1,np
!     do i = 1,nq
!       do j = 1,nr
!         val = dusp(j,i,n-1+n1)
!         cell_dusp(j,i,n) = val
!         dusp_copy(n,j,i) = val
!       end do
!     end do
!   end do
!   !
!   ! Now Interpolate the interior solution gradients to the flux points
!   !
!   ! Loop over the faces of this cell
!   !
!   cell_faces: do ncf = 1,size(cell(nc)%face)
!     !
!     nf = cell(nc)%face(ncf)%idx  ! index of the current face
!     sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
!     !
!     face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
!     face_order = cell(nc)%face(ncf)%order ! order of current face
!     fp_offset = cell(nc)%face(ncf)%fp_offset
!     face_rot = cell(nc)%face(ncf)%rotation
!     !
!     ! Reset face_dusp to the original discontinuous gradient
!     !
!     face_dusp(1:np,1:nr,1:nq) = dusp_copy(1:np,1:nr,1:nq)
!    !!$ACC LOOP VECTOR COLLAPSE(3)
!    !do i = 1,nq
!    !  do j = 1,nr
!    !    do n = 1,np
!    !      face_dusp(n,j,i) = dusp_copy(n,j,i)
!    !    end do
!    !  end do
!    !end do
!     !
!     ! Loop over the points on this face
!     !
!     nfp = geom_solpts(face_geom,face_order)
!     fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
!                                   rot(face_rot)%pts(1:nfp)
!     !
!     ! Compute the BR2 jump correction for this point
!     !    jump_left  =  half * (u_right - u_left)
!     !    jump_right = -half * (u_right - u_left)
!     !
!     ! NOTE: Multiplying by the normal in the loop below that applies the
!     !       correction
!     !           i.e., correct_m * unrm(l)
!     !       takes care of the minus sign if the current cell is the right
!     !       cell on the current face.
!     !
!     ! NOTE: Boundary faces should never be the right cell on the interface.
!     !       For non-wall boundaries, this should reduce to jump_left above.
!     !       For wall boundaries, jump_left should reduce to:
!     !         jump_left = (u_right - u_left) where u_right = u_wall
!     !
!     cbias = zero
!     !
!     ! Reset cbias if this is a wall boundary
!     !
!     if (nf <= nfbnd .and. walls_are_exact) then
!       if (any(bface(1,nf) == bc_walls)) cbias = half
!     end if
!     !
!     hfpb = cbias + half
!     !
!     do kf = 1,nfp
!       dsl = norm2( geofa(nf)%v(1:nr,kf) )
!       fp_unrm(1:nr,kf) = geofa(nf)%v(1:nr,kf)/dsl
!       fp_area(kf) = cl(this_geom) * dsl
!     end do
!     !
!     !$ACC LOOP VECTOR COLLAPSE(3)
!     face_points: do kf = 1,nfp
!       !
!       ! kf: index of flx-pt local to this face
!       ! lf: index of flx-pt within interp%mat
!       !
!       ! Apply the jump correction to the interior discontinuous gradients
!       !
!       do m = 1,nq
!         do l = 1,nr
!           lf = fpidx(kf)
!           do i = correct(this_geom,this_order)%row_ptr(lf)+1, &
!                  correct(this_geom,this_order)%row_ptr(lf+1)
!             n = correct(this_geom,this_order)%row_idx(i)
!             correct_n = correct(this_geom,this_order)%row_val(i) * &
!                         fp_area(kf) / metrics_dt(nc)%jac(n)
!             correct_m = correct_n * hfpb * (faceusp(nf)%v(m,kf,2) - &
!                                             faceusp(nf)%v(m,kf,1))
!            !correct_m = correct(this_geom,this_order)%row_val(i) * &
!            !            (faceusp(nf)%v(m,kf,2) - faceusp(nf)%v(m,kf,1)) * &
!            !            hfpb * fp_area(kf) / metrics_dt(nc)%jac(n)
!             face_dusp(n,l,m) = face_dusp(n,l,m) + correct_m * fp_unrm(l,kf)
!             cell_dusp(l,m,n) = cell_dusp(l,m,n) + correct_m * fp_unrm(l,kf)
!           end do
!         end do
!       end do
!       !
!     end do face_points
!     !
!     ! Interpolate the face corrected gradient from the interior solution
!     ! points to the interface flux points now that we are done with all
!     ! the flux points on the current face of the current cell.
!     !
!     !$ACC LOOP VECTOR COLLAPSE(3)
!     interpolate_face: do k = 1,nfp
!       !
!       ! k: index of flx-pt local to this face
!       ! lf: index of flx-pt within interp%mat
!       !
!       do m = 1,nq
!         do l = 1,nr
!           lf = fpidx(k)
!           val = zero
!           do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
!                  interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
!             n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
!             val = val + face_dusp(n,l,m) * interp(this_geom,this_order)% &
!                                                       toFace(face_order)% &
!                                            row_val(i)
!           end do
!           facedusp(nf)%v(l,m,k,sf) = val
!         end do
!       end do
!       !
!     end do interpolate_face
!     !
!   end do cell_faces
!   !
!   ! Copy the total corrected gradient accounting for all cell faces
!   ! back into dusp
!   !
!  !dusp(1:nr,1:nq,n1:n2) = cell_dusp(1:nr,1:nq,1:np)
!  !$ACC LOOP VECTOR COLLAPSE(3)
!   do n = 1,np
!     do i = 1,nq
!       do j = 1,nr
!         dusp(j,i,n-1+n1) = cell_dusp(j,i,n)
!       end do
!     end do
!   end do
!   !
! end do cell_loop
! !
! ! Exchange the left and right states for faces on communication boundaries
! !
!!if (ncpu > 1) call exchange_facedusp(facedusp)
! !
!end subroutine correct_gradients_br2
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
subroutine check_for_bad_number(line_num)
  !
#ifndef SPECIAL_FOR_GCC
  use, intrinsic :: ieee_exceptions, only : ieee_usual
  use, intrinsic :: ieee_exceptions, only : ieee_get_flag
#endif
  use ovar, only : itcur
  !
  !.. Formal Arguments ..
  integer, intent(in) :: line_num
  !
#ifndef SPECIAL_FOR_GCC
  !
  !.. Local Scalars ..
  character(len=300) :: overflow_error = ""
  character(len=300) :: divide_by_zero_error = ""
  character(len=300) :: invalid_error = ""
  !
  !.. Local Arrays ..
  integer      :: bad_number_found(1:size(ieee_usual))
  logical(ldk) :: signaling_exceptions(1:size(ieee_usual))
  !
  character(len=*), parameter :: pname = "check_for_bad_number"
  !
continue
  !
  ! First check to make sure that no floating-point exceptions are signaling
  !
  call ieee_get_flag(ieee_usual,signaling_exceptions)
  !
  bad_number_found = merge(1,0,signaling_exceptions)
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,bad_number_found, &
                       size(bad_number_found,kind=int_mpi), &
                       mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
  end if
  !
  if (sum(bad_number_found) > 0) then
    !
    if (bad_number_found(1) > 0) write (overflow_error,101)
    if (bad_number_found(2) > 0) write (divide_by_zero_error,102)
    if (bad_number_found(3) > 0) write (invalid_error,103)
    write (error_message,100) itcur, &
                              trim(overflow_error), &
                              trim(divide_by_zero_error), &
                              trim(invalid_error)
    !
    call stop_recon(stop_mpi,pname,line_num,__FILE__,error_message)
    !
  end if
  !
  100 format ("A SIGNALING FLOATING POINT EXCEPTION WAS ", &
              "DETECTED AFTER UPDATING THE SOLUTION AFTER ", &
              "TIME STEP # ",i0," !!! DUMPING SOLUTION AND ", &
              "RESTART FILES FOR THE PREVIOUS TIME STEP!!!", &
              3(2x,a,:))
  !
  101 format ("OVERFLOW EXCEPTION DETECTED!!!")
  102 format ("DIVIDE-BY-ZERO EXCEPTION DETECTED!!!")
  103 format ("INVALID NUMBER (NaN) EXCEPTION DETECTED!!!")
#endif
  !
end subroutine check_for_bad_number
!
!###############################################################################
!
!
!###############################################################################
!
#include "Functions/inviscid_fluxes.f90"
#include "Functions/viscous_fluxes.f90"
!
!###############################################################################
!
include "Functions/upwind_fluxes.f90"
include "Functions/grad_cv_to_grad_pv.f90"
include "Functions/viscosity.f90"
include "Functions/pressure.f90"
!
!###############################################################################
!
end module flux_mod
