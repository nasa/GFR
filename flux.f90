module flux_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use ovar, only : invs_flux_method
  !
  implicit none
  !
  private
  !
 !include "kind_defs.h"
 !include "constants.h"
  !
  public :: correct_flux_project_cv
  public :: correct_flux_project_pv
  public :: correct_flux_project_flx
  public :: flux_divergence_lp
  public :: flux_divergence_cr
  public :: flux_divergence_jac_cr
  public :: dg_weak_form
  public :: dg_strong_form
  public :: interp_cv_to_flxpts
  public :: interp_dusp_to_flxpts
  public :: derflux1_qdtr
  public :: derflux1_qdtr_with_visc
  public :: correct_interface_flux_qdtr
  public :: correct_interface_flux_qdtr_with_visc
  public :: discontinuous_gradient_at_solpts
  public :: correct_gradients_at_solpts
  public :: correct_gradients_br2
  !
  real(wp), parameter :: amco2m = 1.0e-4_wp
  real(wp), parameter :: alphax = 0.1875_wp
  real(wp), parameter :: betax  = 0.125_wp
  !
  ! Generalization parameters for the diffusion scheme
  !
  ! LDG_BETA : Constant that determines the directional bias of
  !            the diffusion scheme for each coordinate direction
  !        = [0,0,0] -> Results in BR2, a centered biased scheme
  !
  ! In 1D, LDG_BETA is a scalar and the resulting schemes are
  !    LDG_BETA = 0    -> BR2, a centered biasing
  !    LDG_BETA = +1/2 -> LDG
  !    LDG_BETA = -1/2 -> CDG
  !   
  real(wp), parameter :: LDG_beta(3) = (/ zero, zero, zero /)
 !real(wp), parameter :: LDG_beta(3) = (/ half, half, half /)
  !
  ! LDG_TAU : Constant that controls the penalty term
  !           term of the common diffusive flux
  !        = 
  real(wp), parameter :: LDG_tau = zero
 !real(wp), parameter :: LDG_tau = one
  !
  ! Parameters giving what set of variables to
  ! average when computing the common viscous flux.
  !
  ! VISC_INTERFACE_FLUXES : Compute the common viscous flux by averaging the
  !                         left and right viscous interface fluxes which are
  !                         computed using the left and right interface
  !                         solutions and gradients.
  ! VISC_SOLUTION_AND_GRADIENTS : Compute common solutions and gradients by
  !                               averaging the left and right interface
  !                               solutions and gradients and use these
  !                               values to compute the common viscous flux.
  !
  integer, parameter :: visc_interface_fluxes = 1
  integer, parameter :: visc_solution_and_gradients = 2
 !integer, parameter :: visc_variables_to_ave = visc_interface_fluxes
 !integer, parameter :: visc_variables_to_ave = visc_solution_and_gradients
  !
  integer, parameter :: ROEFLUX = 1
  integer, parameter :: HLLC = 2
  integer, parameter :: LDFSS = 3
  !
  logical, parameter :: output_boundary_fluxes = fals
 !logical, parameter :: output_boundary_fluxes = true
  logical, save :: first_output_done = fals
  !
  real(wp), parameter :: debug_disc_visc_flux = one
  real(wp), parameter :: debug_common_visc_bnd_flux = one
  real(wp), parameter :: debug_common_visc_flux = one
  !
  integer, save :: flux_cell
  integer, save :: cell_flux_point
  !
contains
!
!###############################################################################
!
include "Functions/upwind_fluxes.f90"
!
!###############################################################################
!
subroutine interp_cv_to_flxpts()
  !
  use geovar,  only : neqns,ncell,n_sp2dt,n_sp2dq,pts_in_cell_ptr
  use geovar,  only : flxpts_of_cell_ptr,flxpts_of_cell
  use geovar,  only : interp_s2f,interp_s2g,n_sp1de,n_fp1de
  use flowvar, only : usp,faceusp
  !
  use parallel_mod, only : exchange_faceusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,nl,k,m
  integer :: nf,f1,f2,kf,lf,sf
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_sp2dq,1:neqns) :: tusp
  real(wp), dimension(1:n_sp2dq,1:4*n_fp1de,1:2) :: interp_mat
  !
continue
  !
  !######################################################################
  !######################################################################
  ! WE NEED TO DO SOMETHING FOR THE COMMUNICATION AND PERIODIC BOUNDARY
  ! CONDITIONS TO INTERPOLATE THE GHOST CELL SOLUTION POINTS TO THE FLUX
  ! POINTS ON THE BOUNDARY FACE.
  !
  ! AFTER CALLING THIS ROUTINE, WE NEED TO EXCHANGE THE FACE DATA
  ! FOR COMMUNICATION AND PERIODIC BOUNDARY CONDITIONS USING MPI
  ! NON-BLOCKING ROUTINES. THIS WAY WE CAN INITIALIZE THIS COMMUNICATION
  ! AND THEN CONTINUE ON TO COMPUTING THE DISCONTINUOUS FLUX DERIVATIVES
  ! AT THE INTERNAL SOLUTION POINTS WHILE THE COMMUNICATION FINISHES.
  ! THIS CAN BE DONE BECAUSE WE DONT NEED THE INTERFACE INFORMATION UNTIL
  ! AFTER THAT STEP.
  !######################################################################
  !######################################################################
  !
  if (n_fp1de /= n_sp1de) then
    interp_mat(:,:,:) = interp_s2g(:,:,:)
  else
    interp_mat(:,:,:) = interp_s2f(:,:,:)
  end if
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    np = n2-n1+1
    nl = merge(1,2,np == n_sp2dt)
    !
    f1 = flxpts_of_cell_ptr(nc)+1
    f2 = flxpts_of_cell_ptr(nc+1)
    !
    ! Create transpose of usp array to reduce cache misses
    !
    do k = 1,np
      do m = 1,neqns
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
    ! Now Interpolate the interior solution gradients to the flux points
    !
    !
    do n = f1,f2
      !
      nf = flxpts_of_cell(1,n)
      kf = flxpts_of_cell(2,n)
      lf = flxpts_of_cell(3,n) ! ideally this should evaluate to 1:f2-f1+1 in this loop
      sf = flxpts_of_cell(4,n)
      !
      !######################################################################
      !######################################################################
      ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
      ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
      ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
      !######################################################################
      !######################################################################
      !
     !write (mypnum+110,'(*(i6))') nc,n,f1,f2,nf,kf,lf,sf
      !
      lmat(1:np) = interp_mat(1:np,lf,nl)
      !
      do m = 1,neqns
        faceusp(m,kf,sf,nf) = dot_product(lmat(1:np),tusp(1:np,m))
      end do
      !
    end do
    !
  end do
  !
  ! Exchange the left and right states for faces on communication boundaries
  !
  if (ncpu > 1) call exchange_faceusp(faceusp)
  !
end subroutine interp_cv_to_flxpts
!
!###############################################################################
!
subroutine interp_dusp_to_flxpts()
  !
  use geovar,  only : nr,neqns,ncell,n_sp2dt,n_sp2dq,pts_in_cell_ptr
  use geovar,  only : flxpts_of_cell_ptr,flxpts_of_cell
  use geovar,  only : interp_s2f,interp_s2g,n_sp1de,n_fp1de
  use flowvar, only : dusp,facedusp
  !
  use parallel_mod, only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,nl,k,l,m
  integer :: nf,f1,f2,kf,lf,sf
  integer :: fp1,fp2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_sp2dq,1:nr,1:neqns) :: tdusp
  real(wp), dimension(1:n_sp2dq,1:4*n_fp1de,1:2) :: interp_mat
  logical, dimension(1:n_sp1de) :: face_check
  !
continue
  !
  !######################################################################
  !######################################################################
  ! WE NEED TO DO SOMETHING FOR THE COMMUNICATION AND PERIODIC BOUNDARY
  ! CONDITIONS TO INTERPOLATE THE GHOST CELL SOLUTION POINTS TO THE FLUX
  ! POINTS ON THE BOUNDARY FACE.
  !
  ! AFTER CALLING THIS ROUTINE, WE NEED TO EXCHANGE THE FACE DATA
  ! FOR COMMUNICATION AND PERIODIC BOUNDARY CONDITIONS USING MPI
  ! NON-BLOCKING ROUTINES. THIS WAY WE CAN INITIALIZE THIS COMMUNICATION
  ! AND THEN CONTINUE ON TO COMPUTING THE DISCONTINUOUS FLUX DERIVATIVES
  ! AT THE INTERNAL SOLUTION POINTS WHILE THE COMMUNICATION FINISHES.
  ! THIS CAN BE DONE BECAUSE WE DONT NEED THE INTERFACE INFORMATION UNTIL
  ! AFTER THAT STEP.
  !######################################################################
  !######################################################################
  !
  if (n_fp1de /= n_sp1de) then
    interp_mat(:,:,:) = interp_s2g(:,:,:)
  else
    interp_mat(:,:,:) = interp_s2f(:,:,:)
  end if
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    np = n2-n1+1
    nl = merge(1,2,np == n_sp2dt)
    !
    f1 = flxpts_of_cell_ptr(nc)+1
    f2 = flxpts_of_cell_ptr(nc+1)
    !
    ! Create transpose of dusp array to reduce cache misses
    !
    do k = 1,np
      do m = 1,neqns
        do l = 1,nr
          tdusp(k,l,m) = dusp(l,m,k+n1-1)
        end do
      end do
    end do
    !
    ! Now Interpolate the interior solution gradients to the flux points
    !
   !do n = 1,(f2-f1+1)/n_sp1de
   !  fp1 = f1 + (n-1)*n_sp1de
   !  fp2 = f1 + n*n_sp1de - 1
   !  write (106,2) (flxpts_of_cell(1,k),k=fp1,fp2)
   !  2 format (*(2x,i6))
   !  do k = fp1,fp2
   !    face_check(k-fp1+1) = all(flxpts_of_cell(1,k)==flxpts_of_cell(1,fp1:fp2))
   !  end do
   !  if (.not. all(face_check)) then
   !    write (iout,3) nc,n,(flxpts_of_cell(1,k),k=fp1,fp2)
   !    3 format ("Error: Cell=",i0,"; Face=",i0,"; FLXPTS =",*(1x,i0))
   !  end if
   !end do
    do n = f1,f2
      !
      nf = flxpts_of_cell(1,n)
      kf = flxpts_of_cell(2,n)
      lf = flxpts_of_cell(3,n) ! ideally this should evaluate to 1:f2-f1+1 in this loop
      sf = flxpts_of_cell(4,n)
     !write (105,1) nc,nf,n
     !1 format ("Cell = ",i4,"; Face = ",i4,"; Cell-FP = ",i0)
      !
      !######################################################################
      !######################################################################
      ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
      ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
      ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
      !######################################################################
      !######################################################################
      !
     !write (mypnum+110,'(*(i6))') nc,n,f1,f2,nf,kf,lf,sf
      !
      lmat(1:np) = interp_mat(1:np,lf,nl)
      !
      do m = 1,neqns
        do l = 1,nr
          facedusp(l,m,kf,sf,nf) = dot_product(lmat(1:np),tdusp(1:np,l,m))
        end do
      end do
      !
    end do
    !
  end do
  !
  ! Exchange the left and right states for faces on communication boundaries
  !
  if (ncpu > 1) call exchange_facedusp(facedusp)
  !
end subroutine interp_dusp_to_flxpts
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Procedures for using quadrature to compute
! the divergence of only inviscid fluxes
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine derflux1_qdtr
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,nr
  use geovar,  only : pts_in_cell_ptr,interp_s2q,n_qp2dq,n_qp2dt
  use geovar,  only : flxpts_of_cell,flxpts_of_cell_ptr,projct_q2s
  use geovar,  only : interp_s2f,interp_s2g,using_quadrature
  use metvar,  only : dx,dy,dr_mat,ds_mat,geofa
  use metvar,  only : dqr_LPmat,dqs_LPmat,jacobian
  use flowvar, only : usp,residual,faceflx,faceusp
  use ovar,    only : dflux_opt
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nw,n1,n2,nsp,nqp
  integer :: f1,f2,nf,sf,kf,lf
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:neqns)   :: ff,gg
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_qp2dq) :: dxdr,dxds
  real(wp), dimension(1:n_qp2dq) :: dydr,dyds
  real(wp), dimension(1:n_qp2dq) :: qjac
  !
  real(wp), dimension(1:n_sp2dq,1:neqns) :: tusp
  real(wp), dimension(1:neqns,1:n_qp2dq) :: uqp
  !
  real(wp), dimension(1:n_sp2dq,1:neqns,1:nr) :: rsflx_sp
  real(wp), dimension(1:n_sp2dq,1:neqns,1:nr) :: xyflx_sp
  real(wp), dimension(1:n_sp2dq,1:neqns,1:nr) :: dflx_sp
  real(wp), dimension(1:n_qp2dq,1:neqns,1:nr) :: rsflx_qp
  real(wp), dimension(1:n_qp2dq,1:neqns,1:nr) :: xyflx_qp
  !
continue
  !
  ! Initialize faceflx to zero before we start
  !
  faceflx(:,:,:,:) = zero
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    nsp = n2-n1+1 ! number of solution points in this cell
    nqp = merge(n_qp2dt,n_qp2dq,nsp==n_sp2dt) ! no. of QP in this cell
    !
    nw = merge( 1 , 2 , nsp == n_sp2dt ) ! check if cell n is tria or quad
    !
    f1 = flxpts_of_cell_ptr(n)+1
    f2 = flxpts_of_cell_ptr(n+1)
    !
    ! Store the transpose of usp to reduce cache misses
    !
    tusp(1:nsp,1:neqns) = transpose(usp(1:neqns,n1:n2))
    !
    ! Get the information local to this cell
    !
    if (using_quadrature) then
      !
      do k = 1,nqp
        !
        ! Get the conservative solution at the quadrature points
        !
        do m = 1,neqns
          uqp(m,k) = dot_product( interp_s2q(1:nsp,k,nw) , tusp(1:nsp,m) )
        end do
        !
        ! Get the metrics for this cell
        !
        dxdr(k) = dot_product( interp_s2q(1:nsp,k,nw) , dx(n1:n2)%dr )
        dxds(k) = dot_product( interp_s2q(1:nsp,k,nw) , dx(n1:n2)%ds )
        dydr(k) = dot_product( interp_s2q(1:nsp,k,nw) , dy(n1:n2)%dr )
        dyds(k) = dot_product( interp_s2q(1:nsp,k,nw) , dy(n1:n2)%ds )
        qjac(k) = dot_product( interp_s2q(1:nsp,k,nw) , jacobian(n1:n2) )
        !
      end do
      !
    else
      !
      uqp(1:neqns,1:nsp) = usp(1:neqns,n1:n2)
      !
      dxdr(1:nsp) = dx(n1:n2)%dr
      dxds(1:nsp) = dx(n1:n2)%ds
      dydr(1:nsp) = dy(n1:n2)%dr
      dyds(1:nsp) = dy(n1:n2)%ds
      qjac(1:nsp) = jacobian(n1:n2)
      !
    end if
    !
    rsflx_qp(:,:,:) = zero ! zero for dot product later with hybrid tria/quad
    xyflx_qp(:,:,:) = zero ! zero for dot product later with hybrid tria/quad
    !
    ! Get the flux at each quadrature point in the current cell
    !
    do k = 1,nqp
      !
      xyflx_qp(k,:,1) = f_invs_flx_cv( uqp(:,k) ) ! x-directional flux
      xyflx_qp(k,:,2) = g_invs_flx_cv( uqp(:,k) ) ! y-directional flux
      !
      rsflx_qp(k,:,1) =  xyflx_qp(k,:,1)*dyds(k) - xyflx_qp(k,:,2)*dxds(k)
      rsflx_qp(k,:,2) = -xyflx_qp(k,:,1)*dydr(k) + xyflx_qp(k,:,2)*dxdr(k)
      !
    end do
    !
    ! Two different methods for getting the derivative at the solution points
    !
    if (.not.using_quadrature .or. dflux_opt == derivative_of_proj_fluxes) then
      !
      ! Project the fluxes back to the solution polynomial order
      !
      if (using_quadrature) then
        !
        do l = 1,nr
          do m = 1,neqns
            do k = 1,nsp
              xyflx_sp(k,m,l) = dot_product( projct_q2s(1:nqp,k,nw) , xyflx_qp(1:nqp,m,l) )
              rsflx_sp(k,m,l) = dot_product( projct_q2s(1:nqp,k,nw) , rsflx_qp(1:nqp,m,l) )
            end do
          end do
        end do
        !
      else
        !
        xyflx_sp(1:nsp,1:neqns,1:nr) = xyflx_qp(1:nsp,1:neqns,1:nr)
        rsflx_sp(1:nsp,1:neqns,1:nr) = rsflx_qp(1:nsp,1:neqns,1:nr)
        !
      end if
      !
      ! Take the derivative of the projected fluxes
      !
      do m = 1,neqns
        do k = 1,nsp
          dflx_sp(k,m,1) = dot_product( dr_mat(1:nsp,k,nw) , rsflx_sp(1:nsp,m,1) )
          dflx_sp(k,m,2) = dot_product( ds_mat(1:nsp,k,nw) , rsflx_sp(1:nsp,m,2) )
        end do
      end do
      !
      ! Interpolate the fluxes to the interface flux points and dot this
      ! interpolated vector with the face normal at each flux point
      ! NOTE: We dont need to exchange this information for communication and
      !       periodic boundary faces since this flux is only used within the
      !       current cell.
      !
      do k = f1,f2
        !
        nf = flxpts_of_cell(1,k)
        kf = flxpts_of_cell(2,k)
        lf = flxpts_of_cell(3,k)
        sf = flxpts_of_cell(4,k)
        !
        ! Get the unit normal for this face flux/solution point
        !
        unrm(:) = unit_vector( geofa(:,kf,nf) )
        !
        ! WE NEED TO ADD SOMETHING HERE IN CASE LOBATTO NODES ARE USED. LOOK AT 
        ! ROUTINE interp_cv_to_flxpts FOR MORE DETAILS.
        !
        lmat(1:nsp) = interp_s2f(1:nsp,lf,nw)
        !
        do l = 1,nr
          do m = 1,neqns
            faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
                                  dot_product(lmat(1:nsp),xyflx_sp(1:nsp,m,l))*unrm(l)
          end do
        end do
        !
      end do
      !
    else if (dflux_opt == project_derivatives) then
      !
      ! Compute the derivative of the fluxes at the quadrature points
      ! and project these derivatives to the solution points
      !
      do k = 1,nsp
        do m = 1,neqns
          dflx_sp(k,m,1) = dot_product( dqr_LPmat(1:nqp,k,nw) , rsflx_qp(1:nqp,m,1) )
          dflx_sp(k,m,2) = dot_product( dqs_LPmat(1:nqp,k,nw) , rsflx_qp(1:nqp,m,2) )
        end do
      end do
      !
      ! Interpolate the fluxes to the interface flux points and dot this
      ! interpolated vector with the face normal at each flux point
      ! NOTE: We dont need to exchange this information for communication and
      !       periodic boundary faces since this flux is only used within the
      !       current cell.
      !
     !do k = f1,f2
     !  !
     !  nf = flxpts_of_cell(1,k)
     !  kf = flxpts_of_cell(2,k)
     !  lf = flxpts_of_cell(3,k)
     !  sf = flxpts_of_cell(4,k)
     !  !
     !  ! #################################################################
     !  ! #################################################################
     !  ! #################################################################
     !  ! ##  ERROR! THIS WILL NOT WORK BECAUSE WE NEED TO INTERPOLATE  ##
     !  ! ##         THE FACE NORMALS TO THE FACE QUADRATURE POINTS!    ##
     !  ! #################################################################
     !  ! #################################################################
     !  ! #################################################################
     !  !
     !  ! Get the unit normal for this face flux/solution point
     !  !
     !  unrm(:) = unit_vector( geofa(:,kf,nf) )
     !  !
     !  ! WE NEED TO ADD SOMETHING HERE IN CASE LOBATTO NODES ARE USED. LOOK AT 
     !  ! ROUTINE interp_cv_to_flxpts FOR MORE DETAILS.
     !  !
     !  lmat(1:nsp) = interp_s2g(1:nsp,lf,nw)
     !  !
     !  do l = 1,nr
     !    do m = 1,neqns
     !      faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
     !                            dot_product(lmat(1:nsp),xyflx_qp(1:nsp,m,l))*unrm(l)
     !    end do
     !  end do
     !  !
     !end do
      !
    end if
    !
    ! Add the derivatives of the flux to the residual
    !
    do m = 1,neqns
      do k = 1,nsp
        residual(m,k+n1-1) = dflx_sp(k,m,1) + dflx_sp(k,m,2)
      end do
    end do
    !
   !do k = 1,nsp
   !  do m = 1,neqns
   !   !residual(m,k+n1-1) = dot_product(dr_mat(1:nsp,k,nw),ftdc(1:nsp,m)) - &
   !   !                     dot_product(ds_mat(1:nsp,k,nw),gtdc(1:nsp,m))
   !   !residual(m,k+n1-1) = dot_product(dqr_mat(:,k,nw),ftdc(:,m)) - &
   !   !                     dot_product(dqs_mat(:,k,nw),gtdc(:,m))
   !    residual(m,k+n1-1) = dot_product(dr_mat(1:nsp,k,nw),dfsp(1:nsp,m)) - &
   !                         dot_product(ds_mat(1:nsp,k,nw),dgsp(1:nsp,m))
   !  end do
   !end do
    !
   !do m = 1,neqns
   !  residual(m,n1:n2) = matmul(invM_mat(:,:,nw),residual(m,n1:n2))!/jacobian(n1:n2)
   !end do
    !
  end do
  !
end subroutine derflux1_qdtr
!
!###############################################################################
!
subroutine correct_interface_flux_qdtr()
  !
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,nr
  use geovar,  only : bface,cells_on_face,pts_in_cell_ptr,xyz,flxwts_1d
  use geovar,  only : interp_s2f,n_qp1de,interp_f2g,projct_g2f,using_quadrature
  use ovar,    only : flux_divergence,dflux_opt,itcur,nst
  use metvar,  only : cormat,geofa,volj
  use flowvar, only : usp,residual,faceflx,faceusp
  !
  use generic_mod,  only : pressure_cv_sp
 !use parallel_mod, only : wait_faceusp
  !
  !.. Local Scalars ..
  integer   :: k,l,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,wl,wr,nsp,nqp
  integer   :: kbl,kbr,kl,kr,iu
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)    :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr,fnc
  real(wp), dimension(1:neqns) :: srcl,srcr
  real(wp), dimension(1:nr,1:n_qp1de)    :: rnrm
  real(wp), dimension(1:neqns,1:n_sp1de) :: fjl,fjr
  real(wp), dimension(1:neqns,1:n_fp1de) :: lflx,rflx
  real(wp), dimension(1:neqns,1:n_qp1de) :: luqp,ruqp,fncq
  !
continue
  !
  nsp = n_sp1de
  nqp = n_qp1de
  !
  ! Make sure that the non-blocking exchange of faceusp has completed
  !
 !if (ncpu > 1) call wait_faceusp(faceusp)
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  if (output_boundary_fluxes) write (102,110) itcur,nst
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! ending    flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    wl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , r2-r1+1 == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Interpolate the conserved variables and face normals at
    ! the solution flux points to the quadrature flux points
    !
    if (using_quadrature) then
      !
      do k = 1,nqp
        !
        do m = 1,neqns
          luqp(m,k) = dot_product( interp_f2g(1:nsp,k) , faceusp(m,1:nsp,1,nf) )
          ruqp(m,k) = dot_product( interp_f2g(1:nsp,k) , faceusp(m,1:nsp,2,nf) )
        end do
        !
        do l = 1,nr
          rnrm(l,k) = dot_product( interp_f2g(1:nsp,k) , geofa(l,1:nsp,nf) )
        end do
        !
      end do
      !
    else
      !
      lflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,1,nf)
      rflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,2,nf)
      luqp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,1,nf)
      ruqp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,2,nf)
      rnrm(1:nr,1:nsp) = geofa(1:nr,1:nsp,nf)
      !
    end if
    !
    if (output_boundary_fluxes) then
      if (nf <= nfbnd) then
        write (102,111) nf,trim(bc_integer_to_string(bface(1,nf)))
        do k = 1,nsp
          write (102,112) k,luqp(1:neqns,k),pressure_cv_sp(luqp(1:neqns,k))
          write (102,113) k,ruqp(1:neqns,k),pressure_cv_sp(ruqp(1:neqns,k))
        end do
      end if
    end if
    !
    ! Get the discontinuous fluxes at the interface
    !
    if (dflux_opt == project_derivatives) then
      !
      ! Get the discontinuous fluxes at the quadrature flux points
      !
      do k = 1,nqp
        lflx(1:neqns,k) = normal_invs_flx_cv( luqp(:,k) , rnrm(:,k) )
        rflx(1:neqns,k) = normal_invs_flx_cv( ruqp(:,k) , rnrm(:,k) )
      end do
      !
    else
      !
      ! Get the discontinuous fluxes at the solution flux points
      !
      lflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,1,nf)
      rflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,2,nf)
      !
    end if
    !
    ! Get the common flux at each quadrature flux point
    !
    do k = 1,nqp
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( rnrm(:,k) )
      !
      ! Get the numerical flux at the current flux point
      !
      if (invs_flux_method == ROEFLUX) then
        fncq(:,k) = upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (invs_flux_method == HLLC) then
        fncq(:,k) = hllc_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (invs_flux_method == LDFSS) then
        fncq(:,k) = ldfss_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      end if
      !
    end do
    !
    ! Project the common flux to the solution flux points
    !
    if (using_quadrature) then
      !
      if (dflux_opt == project_derivatives) then
        !
        lflx(1:neqns,1:nqp) = fncq(1:neqns,1:nqp) - lflx(1:neqns,1:nqp)
        rflx(1:neqns,1:nqp) = fncq(1:neqns,1:nqp) - rflx(1:neqns,1:nqp)
        !
        do k = 1,nsp
          do m = 1,neqns
            fjl(m,k) = dot_product( projct_g2f(1:nqp,k) , lflx(m,1:nqp) )
            fjr(m,k) = dot_product( projct_g2f(1:nqp,k) , rflx(m,1:nqp) )
          end do
        end do
        !
      else
        !
        do k = 1,nsp
          !
          do m = 1,neqns
            fnc(m) = dot_product( projct_g2f(1:nqp,k) , fncq(m,1:nqp) )
          end do
          !
          fjl(:,k) = fnc(:) - lflx(:,k)
          fjr(:,k) = fnc(:) - rflx(:,k)
          !
        end do
        !
      end if
      !
    else
      !
      fjl(1:neqns,1:nsp) = fncq(1:neqns,1:nsp) - lflx(1:neqns,1:nsp)
      fjr(1:neqns,1:nsp) = fncq(1:neqns,1:nsp) - rflx(1:neqns,1:nsp)
      ! 
    end if
    !
    ! Compute the flux jump at each solution flux point
    !
   !do k = 1,nsp
   !  fjl(:,k) = fnc(:,k) - lflx(:,k)
   !  fjr(:,k) = fnc(:,k) - rflx(:,k)
   !end do
   !!
   !! Initialize the conservation source terms
   !!
   !srcl(:) = zero
   !srcr(:) = zero
   !!
   !1 format ("     lflx        rflx        fjl         fjr         fnc",/,*(i8))
   !4 format ("     luqp        ruqp",/,*(i8))
   !3 format (*(i6))
   !2 format (*(es12.4))
   !iu = 103 
   !write (iu,1) nf,jl,jr
   !write (iu+1,4) nf,jl,jr
    do k = 1,nsp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = kbl-1 + k ! index for left  cell
      kr = kbr+1 - k ! index for right cell
      !
     !write (iu,3) k,kl,kr
     !write (iu+1,3) k,kl,kr
     !do m = 1,neqns
     !  write (iu,2) lflx(m,k),rflx(m,k),fjl(m,k),fjr(m,k),fnc(m,k)
     !  write (iu+1,2) luqp(m,k),ruqp(m,k)
     !end do
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( rnrm(:,k) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m,k)*cormat(n-l1+1,kl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m,k)*cormat(n-r1+1,kr,wr)
          end do
        end do
      end if
      !
      ! From ZJ, to ensure conservation with chain-rule approach
      !
     !if (flux_divergence == ChainRuleCV) then
     !  do m = 1,neqns   ! ARE DSL AND DSR CORRECT HERE?
     !    srcl(m) = srcl(m) + dsl*lflx(m,k)*flxwts_1d(k)
     !    srcr(m) = srcr(m) + dsr*rflx(m,k)*flxwts_1d(k)
     !  end do
     !end if
      !
    end do
    !
    ! Now add this faces contribution to the surface conservation
    ! source term to each solution point in the left and right cells
    !
   !if (flux_divergence == ChainRuleCV) then
   !  !
   !  do n = l1,l2
   !    do m = 1,neqns
   !      residual(m,n) = residual(m,n) + srcl(m)/volj(n)
   !    end do
   !  end do
   !  !
   !  if (nf > nfbnd) then
   !    do n = r1,r2
   !      do m = 1,neqns
   !        residual(m,n) = residual(m,n) + srcr(m)/volj(n)
   !      end do
   !    end do
   !  end if
   !  !
   !end if
    !
  end do
  !
  110 format (//"itcur = ",i0,"; rkstep = ",i0)
  111 format ("Boundary face ",i0," with BC = ",a)
  112 format ("FP-",i0,"; Left : ",*(es14.6))
  113 format ("FP-",i0,"; Right: ",*(es14.6))
  !
end subroutine correct_interface_flux_qdtr
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
subroutine derflux1_qdtr_with_visc
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,nr
  use geovar,  only : pts_in_cell_ptr,interp_s2q,n_qp2dq,n_qp2dt
  use geovar,  only : flxpts_of_cell,flxpts_of_cell_ptr,projct_q2s
  use geovar,  only : interp_s2f,interp_s2g,using_quadrature
  use metvar,  only : dx,dy,dr_mat,ds_mat,geofa
  use metvar,  only : dqr_LPmat,dqs_LPmat,jacobian
  use flowvar, only : usp,dusp,residual,faceflx,faceusp
  use ovar,    only : dflux_opt
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nw,n1,n2,np,nq
  integer :: f1,f2,nf,sf,kf,lf
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:neqns,1:nr)   :: vflx
  real(wp), dimension(1:n_sp2dq) :: lmat
  !
  real(wp), dimension(1:n_sp2dq,1:neqns,1:nr) :: rsflx
  real(wp), dimension(1:n_sp2dq,1:neqns,1:nr) :: xyflx
  !
continue
  !
  nq = neqns
  !
  ! Initialize faceflx to zero before we start
  !
  faceflx(:,:,:,:) = zero
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    f1 = flxpts_of_cell_ptr(n)+1 ! beginning index of the flx-pts for cell n
    f2 = flxpts_of_cell_ptr(n+1) ! ending    index of the flx-pts for cell n
    !
    if (.not. first_output_done .and. output_boundary_fluxes) then
      write (120,110) n
      do k = n1,n2
        write (120,111) k,dx(k)%dr,dx(k)%ds,dy(k)%dr,dy(k)%ds,jacobian(k) 
      end do
    end if
    110 format ("Cell ",i0)
    111 format ("SP = ",i0,"; dxdr=",es14.6,", dxds=",es14.6, &
                ", dydr=",es14.6,", dyds=",es14.6,", jacobian=",es14.6)
    !
    ! Get the flux at each quadrature point in the current cell
    !
    do k = 1,np
      !
      l = k + n1-1
      !
      ! Inviscid fluxes
      !
      xyflx(k,:,1) = f_invs_flx_cv( usp(:,l) ) ! x-directional flux
      xyflx(k,:,2) = g_invs_flx_cv( usp(:,l) ) ! y-directional flux
      !
      ! Viscous fluxes
      !
      flux_cell = n
      cell_flux_point = k
      vflx(1:nq,1:nr) = visc_flx_cv( usp(1:nq,l) , dusp(1:nr,1:nq,l) )
      !
      xyflx(k,1:nq,1:nr) = xyflx(k,1:nq,1:nr) - debug_disc_visc_flux*vflx(1:nq,1:nr)
      !
      ! Transform the fluxes in the x-y physical
      ! domain to the r-s computational domain
      !
      rsflx(k,:,1) =  xyflx(k,:,1)*dy(l)%ds - xyflx(k,:,2)*dx(l)%ds
      rsflx(k,:,2) = -xyflx(k,:,1)*dy(l)%dr + xyflx(k,:,2)*dx(l)%dr
      !
    end do
    !
    ! Compute the divergence of the discontinuous fluxes 
    !
    do m = 1,nq
      do k = 1,np
        residual(m,k+n1-1) = dot_product( dr_mat(1:np,k,nw) , rsflx(1:np,m,1) ) + &
                             dot_product( ds_mat(1:np,k,nw) , rsflx(1:np,m,2) )
      end do
    end do
    !
    ! Interpolate the fluxes to the interface flux points and dot this
    ! interpolated vector with the face normal at each flux point
    ! NOTE: We dont need to exchange this information for communication and
    !       periodic boundary faces since this flux is only used within the
    !       current cell.
    !
    do k = f1,f2
      !
      nf = flxpts_of_cell(1,k)
      kf = flxpts_of_cell(2,k)
      lf = flxpts_of_cell(3,k)
      sf = flxpts_of_cell(4,k)
      !
      ! Get the unit normal for this face flux/solution point
      !
      unrm(:) = unit_vector( geofa(:,kf,nf) )
      !
      ! WE NEED TO ADD SOMETHING HERE IN CASE LOBATTO NODES ARE USED. LOOK AT 
      ! ROUTINE interp_cv_to_flxpts FOR MORE DETAILS.
      !
      lmat(1:np) = interp_s2f(1:np,lf,nw)
      !
      do l = 1,nr
        do m = 1,nq
          faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
                                dot_product(lmat(1:np),xyflx(1:np,m,l))*unrm(l)
        end do
      end do
      !
    end do
    !
  end do
  !
  first_output_done = true
  !
end subroutine derflux1_qdtr_with_visc
!
!###############################################################################
!
subroutine correct_interface_flux_qdtr_with_visc()
  !
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,nr
  use geovar,  only : bface,cells_on_face,pts_in_cell_ptr,xyz,flxwts_1d
  use geovar,  only : interp_s2f,n_qp1de,interp_f2g,projct_g2f,using_quadrature
  use ovar,    only : flux_divergence,dflux_opt,itcur,nst
  use metvar,  only : cormat,geofa,volj
  use flowvar, only : usp,residual,faceflx,faceusp,facedusp
  !
  use generic_mod, only : pressure_cv_sp
 !use parallel_mod, only : wait_facedusp
  !
  !.. Local Scalars ..
  integer   :: k,l,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,wl,wr,nsp
  integer   :: kbl,kbr,kl,kr,iu
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)    :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr,fnc
  real(wp), dimension(1:neqns) :: srcl,srcr
  real(wp), dimension(1:neqns) :: iflx,vflx
  real(wp), dimension(1:nr,1:n_sp1de)    :: rnrm
  real(wp), dimension(1:neqns,1:n_sp1de) :: fjl,fjr
  real(wp), dimension(1:neqns,1:n_sp1de) :: lflx,rflx
  real(wp), dimension(1:neqns,1:n_sp1de) :: luqp,ruqp,fncq
  real(wp), dimension(1:nr,1:neqns,1:n_sp1de) :: lduqp,rduqp
  !
continue
  !
  nsp = n_sp1de
  !
  ! Make sure that the non-blocking exchange of faceusp has completed
  !
 !if (ncpu > 1) call wait_facedusp(facedusp)
  !
  ! ##############################################################
  ! ##############################################################
  ! ##############################################################
  ! ##########              BOUNDARY FACES              ##########
  ! ##############################################################
  ! ##############################################################
  ! ##############################################################
  !
  if (output_boundary_fluxes) write (101,110) itcur,nst
  !
  ! Loop through the boundary faces to apply the correction
  ! due to the boundary fluxes.
  !
  boundary_fluxes: do nf = 1,nfbnd
    !
   !write (102,'("Boundary face ",i0)') nf
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    wl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if jl is a tria or quad
    !
    lflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,1,nf)
    !
    luqp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,1,nf)
    ruqp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,2,nf)
    !
    if (output_boundary_fluxes) then
      if (nf <= nfbnd) then
        write (101,111) nf,trim(bc_integer_to_string(bface(1,nf)))
        do k = 1,nsp
          write (101,112) k,luqp(1:neqns,k),pressure_cv_sp(luqp(1:neqns,k))
          write (101,113) k,ruqp(1:neqns,k),pressure_cv_sp(ruqp(1:neqns,k))
        end do
      end if
    end if
    !
    lduqp(1:nr,1:neqns,1:nsp) = facedusp(1:nr,1:neqns,1:nsp,1,nf)
    rduqp(1:nr,1:neqns,1:nsp) = facedusp(1:nr,1:neqns,1:nsp,2,nf)
    !
    rnrm(1:nr,1:nsp) = geofa(1:nr,1:nsp,nf)
    !
    ! Get the common flux at each quadrature flux point
    !
    do k = 1,nsp
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
        vflx(:) = visc_wall_flx_cv( unrm(:) , ruqp(:,k) , rduqp(:,:,k) )
        !
      else
        !
        ! For all other boundaries, compute the boundary
        ! fluxes in the same way as the interior faces.
        !
        if (invs_flux_method == ROEFLUX) then
          iflx(:) = upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
        else if (invs_flux_method == HLLC) then
          iflx(:) = hllc_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
        else if (invs_flux_method == LDFSS) then
          iflx(:) = ldfss_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
        end if
        !
        ! Get the common viscous flux at the current flux point
        !
        flux_cell = jl
        cell_flux_point = k
        vflx(:) = common_visc_flx_cv( unrm(:), luqp(:,k), lduqp(:,:,k), &
                                               ruqp(:,k), rduqp(:,:,k) )
        !
      end if
      if (output_boundary_fluxes) then
        write (103,201) nf,k,trim(bc_integer_to_string(bface(1,nf)))
        write (103,203) iflx(1:neqns)
        write (103,204) vflx(1:neqns)
        write (103,205) lflx(1:neqns,k)
      end if
      !
      ! Get the total common boundary flux at this flux point
      !
      fncq(:,k) = iflx(:) - debug_common_visc_bnd_flux*vflx(:)
      !
    end do
    !
    ! Compute the flux jump at the interface
    !
    fjl(1:neqns,1:nsp) = fncq(1:neqns,1:nsp) - lflx(1:neqns,1:nsp)
    !
    ! Apply the flux correction due to the flux jump
    ! NOTE: We only need to apply the correction to the left cell because
    !       this is a boundary face and the right cell is a ghost cell.
    !
    do k = 1,nsp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = kbl-1 + k ! index for left  cell
      !
      ! Get the face length for the left cell
      !
      ds  = rnorm2( rnrm(:,k) )
      dsl = ds * cl(wl)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m,k)*cormat(n-l1+1,kl,wl)
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
   !write (102,'("Interior face ",i0)') nf
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! ending    flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    wl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , r2-r1+1 == n_sp2dt ) ! check if jr is a tria or quad
    !
    lflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,1,nf)
    rflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,2,nf)
    !
    luqp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,1,nf)
    ruqp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,2,nf)
    !
    lduqp(1:nr,1:neqns,1:nsp) = facedusp(1:nr,1:neqns,1:nsp,1,nf)
    rduqp(1:nr,1:neqns,1:nsp) = facedusp(1:nr,1:neqns,1:nsp,2,nf)
    !
    rnrm(1:nr,1:nsp) = geofa(1:nr,1:nsp,nf)
    !
    ! Get the common flux at each quadrature flux point
    !
    do k = 1,nsp
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( rnrm(:,k) )
      !
      ! Get the inviscid upwind flux at the current flux point
      !
      if (invs_flux_method == ROEFLUX) then
        iflx(:) = upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (invs_flux_method == HLLC) then
        iflx(:) = hllc_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (invs_flux_method == LDFSS) then
        iflx(:) = ldfss_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      end if
      !
      ! Get the common viscous flux at the current flux point
      !
      flux_cell = jl
      cell_flux_point = jr
      vflx(:) = common_visc_flx_cv( unrm(:), luqp(:,k), lduqp(:,:,k), &
                                             ruqp(:,k), rduqp(:,:,k) )
      if (output_boundary_fluxes) then
        write (104,202) nf,k
        write (104,203) iflx(1:neqns)
        write (104,204) vflx(1:neqns)
        write (104,205) lflx(1:neqns,k)
        write (104,206) rflx(1:neqns,k)
      end if
      !
      ! Get the total common boundary flux at this flux point
      !
      fncq(:,k) = iflx(:) - debug_common_visc_flux*vflx(:)
      !
    end do
    !
    ! Compute the flux jump at the interface
    !
    fjl(1:neqns,1:nsp) = fncq(1:neqns,1:nsp) - lflx(1:neqns,1:nsp)
    fjr(1:neqns,1:nsp) = fncq(1:neqns,1:nsp) - rflx(1:neqns,1:nsp)
    !
    ! Apply the flux correction due to the flux jump
    !
    do k = 1,nsp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = kbl-1 + k ! index for left  cell
      kr = kbr+1 - k ! index for right cell
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( rnrm(:,k) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m,k)*cormat(n-l1+1,kl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      do n = r1,r2
        do m = 1,neqns
          residual(m,n) = residual(m,n) - dsr*fjr(m,k)*cormat(n-r1+1,kr,wr)
        end do
      end do
      !
    end do
    !
  end do interior_fluxes
  !
  ! Format Statements
  !
  110 format (//"itcur = ",i0,"; rkstep = ",i0)
  111 format ("Boundary face ",i0," with BC = ",a)
  112 format ("FP-",i0,"; Left : ",*(es14.6))
  113 format ("FP-",i0,"; Right: ",*(es14.6))
  201 format ("Boundary Face = ",i0,"; FP-",i0,"; with BC = ",a)
  202 format ("Interior Face = ",i0,"; FP-",i0)
  203 format ("Invs  Fluxes: ",4(1x,es23.15))
  204 format ("Visc  Fluxes: ",4(1x,es23.15))
  205 format ("Left  Fluxes: ",4(1x,es23.15))
  206 format ("Right Fluxes: ",4(1x,es23.15))
  !
end subroutine correct_interface_flux_qdtr_with_visc
!
!###############################################################################
!
subroutine correct_flux_project_cv()
  !
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_sp2dt,n_sp2dq
  use geovar,  only : bface,cells_on_face,pts_in_cell_ptr,xyz,flxwts_1d
  use ovar,    only : flux_divergence
  use metvar,  only : cormat,geofa,volj
  use geovar,  only : interp_s2f
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  logical   :: is_1d_bc
  integer   :: k,l,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,nl,nr,wl,wr
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr
  real(wp), dimension(1:neqns) :: fjl,fjr
  real(wp), dimension(1:neqns) :: fnc
  real(wp), dimension(1:neqns) :: srcl,srcr
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rusp
  !
continue
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    ! Check if this is a boundary face needing a 1D treatment of the flux points
    !
    is_1d_bc = fals
    if (nf <= nfbnd) is_1d_bc = .not. any(bface(1,nf) == bc_comm)
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Store the transpose of usp for the left
    ! and right cells to reduce cache misses
    !
    do k = l1,l2
      lusp (k-l1+1,:) = usp(:,k)
    end do
    !
    do k = r1,r2
      rusp(k-r1+1,:) = usp(:,k)
    end do
   !!
   !! Initialize the conservation source terms
   !!
   !srcl(:) = zero
   !srcr(:) = zero
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the conservative variables
      ! from the left cells to the current flux point on the interface
      !
      do m = 1,neqns
        ul(m) = dot_product(interp_s2f(1:nl,kcorl,wl),usp(m,l1:l2)) ! left  state
      end do
      !
      ! Now get the information at flux point on the right side
      !
      if (is_1d_bc) then
        !
        ! Copy the specified boundary condition in usp into ur and use
        ! this value to get the discontinuous fluxes at this flux point
        !
        ur(:) = rusp(kcorr,:)
        !
      else
        !
        ! Use the Lagrange polynomial to project the conserved solution
        ! variables and the discontinuous fluxes from the right (non-ghost)
        ! cell to the current flux point on the interface.
        !
        do m = 1,neqns
          ur(m)  = dot_product(interp_s2f(1:nr,kcorr,wr), rusp(1:nr,m)) ! right solution
        end do
        !
      end if
      !
      ! Get the numerical flux at the current flux point
      !
      if (invs_flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , ul , ur )
      else if (invs_flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , ul , ur )
      else if (invs_flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , ul , ur )
      end if
      !
      ! Get left and right discontinuous fluxes at the current flux point
      !
      dfl(:) = f_invs_flx_cv(ul)
      dgl(:) = g_invs_flx_cv(ul)
      dfr(:) = f_invs_flx_cv(ur)
      dgr(:) = g_invs_flx_cv(ur)
      !
      ! Get the discontinuous flux function at the current
      ! interface point in the direction of the face normal
      !
      fnl(:) = dfl(:)*unrm(1) + dgl(:)*unrm(2)
      fnr(:) = dfr(:)*unrm(1) + dgr(:)*unrm(2)
      !
      ! The jump in the flux at the current interface point
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m)*cormat(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m)*cormat(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
      ! From ZJ, to ensure conservation with chain-rule approach
      !
     !if (flux_divergence == ChainRuleCV) then
     !  do m = 1,neqns   ! ARE DSL AND DSR CORRECT HERE?
     !    srcl(m) = srcl(m) + dsl*fnl(m)*flxwts_1d(k)
     !    srcr(m) = srcr(m) + dsr*fnr(m)*flxwts_1d(k)
     !  end do
     !end if
      !
    end do
    !
    ! Now add this faces contribution to the surface conservation
    ! source term to each solution point in the left and right cells
    !
   !if (flux_divergence == ChainRuleCV) then
   !  !
   !  do n = l1,l2
   !    do m = 1,neqns
   !      residual(m,n) = residual(m,n) + srcl(m)/volj(n)
   !    end do
   !  end do
   !  !
   !  if (nf > nfbnd) then
   !    do n = r1,r2
   !      do m = 1,neqns
   !        residual(m,n) = residual(m,n) + srcr(m)/volj(n)
   !      end do
   !    end do
   !  end if
   !  !
   !end if
    !
  end do
  !
end subroutine correct_flux_project_cv
!
!###############################################################################
!
subroutine correct_flux_project_pv()
  !
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_sp2dt, &
                      cells_on_face,pts_in_cell_ptr
  use geovar,  only : interp_s2f
  use metvar,  only : cormat,geofa
  use flowvar, only : vsp,residual
  use generic_mod, only : vsp2u_sp
  !
  !.. Local Scalars ..
  integer   :: k,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,nl,nr,wl,wr
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: vl,vr
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr
  real(wp), dimension(1:neqns) :: fjl,fjr
  real(wp), dimension(1:neqns) :: fnc
  !
continue
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the primitive variables from
      ! the left and right cells to the current flux point on the interface
      !
      do m = 1,neqns
        vl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),vsp(m,l1:l2)) ! left  state
        vr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),vsp(m,r1:r2)) ! right state
      end do
      !
      ! Get the numerical flux at the current flux point
      !
      if (invs_flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , vsp2u_sp(vl) , vsp2u_sp(vr) )
      else if (invs_flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , vsp2u_sp(vl) , vsp2u_sp(vr) )
      else if (invs_flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , vsp2u_sp(vl) , vsp2u_sp(vr) )
      end if
      !
      ! Get left and right discontinuous fluxes at the current flux point
      !
      dfl(:) = f_invs_flx_pv(vl)
      dgl(:) = g_invs_flx_pv(vl)
      dfr(:) = f_invs_flx_pv(vr)
      dgr(:) = g_invs_flx_pv(vr)
      !
      ! Get the discontinuous flux function at the current
      ! interface point in the direction of the face normal
      !
      fnl(:) = dfl(:)*unrm(1) + dgl(:)*unrm(2)
      fnr(:) = dfr(:)*unrm(1) + dgr(:)*unrm(2)
      !
      ! The jump in the flux at the current interface point
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m)*cormat(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m)*cormat(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine correct_flux_project_pv
!
!###############################################################################
!
subroutine correct_flux_project_flx()
  !
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_sp2dt,n_sp2dq, &
                      cells_on_face,pts_in_cell_ptr,bface,xyz
  use geovar,  only : interp_s2f
  use metvar,  only : cormat,geofa,dx,dy
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  logical   :: is_1d_bc
  integer   :: k,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,nl,nr,wl,wr,iu
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr
  real(wp), dimension(1:neqns) :: fjl,fjr
  real(wp), dimension(1:neqns) :: fnc
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxg
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxg
  !
continue
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    ! Check if this is a boundary face needing a 1D treatment of the flux points
    !
    is_1d_bc = fals
    if (nf <= nfbnd) is_1d_bc = .not. any(bface(1,nf) == bc_comm)
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Compute the fluxes at each solution point and store the
    ! transpose of usp for the left cell to reduce cache misses
    !
    do k = l1,l2
      lflxf(k-l1+1,:) = f_invs_flx_cv( usp(:,k) )
      lflxg(k-l1+1,:) = g_invs_flx_cv( usp(:,k) )
      lusp (k-l1+1,:) = usp(:,k)
    end do
    !
    ! Store the transpose of usp for the right cell to reduce cache misses
    !
    do k = r1,r2
      rusp(k-r1+1,:) = usp(:,k)
    end do
    !
    ! Compute the right cell fluxes at each solution point for only the
    ! interior cells and communication/periodic boundary conditions
    !
   !1 format ("     lflx        rflx        fjl         fjr         fnc",/,*(i8))
   !4 format ("     luqp        ruqp",/,*(i8))
   !3 format (*(i6))
   !2 format (i6,*(es12.4))
   !iu = 105+10*mypnum
   !write (iu,1) nf,jl,jr
   !write (iu+1,3) nf,jl,jr
   !if (nf <= nfbnd) write (iu+1,3) (bface(m,nf),m=1,size(bface,1))
   !do k = r1,r2
   !  write (iu+1,2) k,(usp(m,k),m=1,neqns),(xyz(m,k),m=1,2)
   !end do
    if (.not. is_1d_bc) then
      do k = r1,r2
        rflxf(k-r1+1,:) = f_invs_flx_cv( usp(:,k) )
        rflxg(k-r1+1,:) = g_invs_flx_cv( usp(:,k) )
      end do
    end if
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the conserved
      ! solution variables and the discontinuous fluxes from the
      ! left cell to the current flux point on the interface
      !
      do m = 1,neqns
        ul(m)  = dot_product(interp_s2f(1:nl,kcorl,wl), lusp(1:nl,m)) ! left solution
        dfl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxf(1:nl,m)) ! left f flux
        dgl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxg(1:nl,m)) ! left g flux
      end do
      !
      ! Now get the information at flux point on the right side
      !
      if (is_1d_bc) then
        !
        ! Copy the specified boundary condition in usp into ur and use
        ! this value to get the discontinuous fluxes at this flux point
        !
        !########################################################
        !########################################################
        ! NEED TO CHECK IF kcorr IS THE RIGHT INDEX THAT NEEDS TO
        ! BE USED HERE. I THINK kcorr VARIES FROM 1:3*n_sp1de
        ! AND NOT FROM 1:n_sp2dt OR 1:n_sp2dq. 
        !########################################################
        !########################################################
        !
       !write (*,*) "kcorr = ",kcorr
        ur(:) = rusp(kcorr,:)
        dfr(:) = f_invs_flx_cv( ur )
        dgr(:) = g_invs_flx_cv( ur )
        !
      else
        !
        ! Use the Lagrange polynomial to project the conserved solution
        ! variables and the discontinuous fluxes from the right (non-ghost)
        ! cell to the current flux point on the interface.
        !
        do m = 1,neqns
          ur(m)  = dot_product(interp_s2f(1:nr,kcorr,wr), rusp(1:nr,m)) ! right solution
          dfr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxf(1:nr,m)) ! right f flux
          dgr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxg(1:nr,m)) ! right g flux
        end do
        !
      end if
      !
      ! Get the numerical flux at the current flux point
      !
      if (invs_flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , ul , ur )
      else if (invs_flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , ul , ur )
      else if (invs_flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , ul , ur )
      end if
      !
      ! Get the discontinuous flux function at the current
      ! interface point in the direction of the face normal
      !
      fnl(:) = dfl(:)*unrm(1) + dgl(:)*unrm(2)
      fnr(:) = dfr(:)*unrm(1) + dgr(:)*unrm(2)
      !
      ! The jump in the flux at the current interface point
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
     !write (iu,3) k,kcorl,kcorr
     !write (iu+1,3) k,kcorl,kcorr
     !do m = 1,neqns
     !  write (iu,2) fnl(m),fnr(m),fjl(m),fjr(m),fnc(m)
     !  write (iu+1,2) ul(m),ur(m)
     !end do
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m)*cormat(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m)*cormat(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine correct_flux_project_flx
!
!###############################################################################
!
subroutine flux_divergence_lp
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp2dt,n_sp2dq, &
                      pts_in_cell_ptr,write_variable
  use metvar,  only : dx,dy,dr_mat,ds_mat
  use flowvar, only : usp,vsp,residual
  !
  !.. Local Scalars ..
  integer :: k,m,n,nw,n1,n2,np
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns)           :: ff,gg
  real(wp), dimension(1:n_sp2dq,1:neqns) :: ftdc,gtdc
  !
continue
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    ftdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    gtdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    !
    ! Get the flux at each solution point in the current cell
    !
    do k = n1,n2
      !
      ff(:) = f_invs_flx_cv( usp(:,k) ) ! x-directional flux
      gg(:) = g_invs_flx_cv( usp(:,k) ) ! y-directional flux
      !
      ftdc(k-n1+1,:) =  ff(:)*dy(k)%ds - gg(:)*dx(k)%ds
      gtdc(k-n1+1,:) =  ff(:)*dy(k)%dr - gg(:)*dx(k)%dr
      !
    end do
    !
    ! Compute the derivatives of the fluxes at each
    ! solution point and add to the residual
    !
    do k = 1,np
      do m = 1,neqns
        residual(m,k+n1-1) = dot_product(dr_mat(1:np,k,nw),ftdc(1:np,m)) - &
                             dot_product(ds_mat(1:np,k,nw),gtdc(1:np,m))
      end do
    end do
    !
  end do
  !
end subroutine flux_divergence_lp
!
!###############################################################################
!
subroutine flux_divergence_cr
  !
  !... this subroutine get derivative of fluxes (no interaction) by
  !... the chain rule (non-conservative formula,
  !...  eq. 95 p. 124 Thompson et al. Grid Generation book)
  !... for each cell, get uxi and uet at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp2dt,n_sp2dq, &
                      pts_in_cell_ptr
  use metvar,  only : dx,dy,dr_mat,ds_mat
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  integer :: k,m,n,n1,n2,nw,np,ip
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns)           :: src
  real(wp), dimension(1:neqns,1:neqns)   :: dfdu,dgdu
  real(wp), dimension(1:neqns,1:neqns)   :: c_f,c_g
  real(wp), dimension(1:neqns,1:n_sp2dq) :: ucr,ucs
  real(wp), dimension(1:n_sp2dq,1:neqns) :: tusp
  !
continue
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge(1,2,np == n_sp2dt) ! check if cell n is tria or quad
    !
    ! Store the transpose of usp for this cell to reduce cache misses
    !
    do k = 1,np
      do m = 1,neqns
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
    ! Compute the derivative of the conserved variables at each solution point
    !
    do k = 1,np
      do m = 1,neqns
        ucr(m,k) = dot_product( dr_mat(1:np,k,nw) , tusp(1:np,m) )
        ucs(m,k) = dot_product( ds_mat(1:np,k,nw) , tusp(1:np,m) )
      end do
    end do
   !!
   !! Initialize the conservation source term
   !!
   !src(:) = zero
    !
    ! Compute the derivatives of the fluxes at each solution point
    !
    do k = n1,n2
      !
      ip = k-n1+1 ! local index of current solution point
      !
      ! Compute the Jacobian of the fluxes at the current sol-pt
      !
      dfdu(:,:) = f_invs_flx_jacobian_cv( usp(:,k) )
      dgdu(:,:) = g_invs_flx_jacobian_cv( usp(:,k) )
      !
      ! Compute the Jacobian in the direction of the face normal
      !
      c_f(:,:) = dy(k)%ds*dfdu(:,:) - dx(k)%ds*dgdu(:,:)
      c_g(:,:) = dy(k)%dr*dfdu(:,:) - dx(k)%dr*dgdu(:,:)
      !
      ! Use the chain rule and add to the residual
      !
      do m = 1,neqns
        residual(m,k) = dot_product(ucr(:,ip),c_f(:,m)) - &
                        dot_product(ucs(:,ip),c_g(:,m))
      end do
   !  !
   !  ! From ZJ, to ensure conservation with chain-rule approach
   !  !
   !  do m = 1,neqns
   !    src(m) = src(m) + residual(m,k)*jacobian(k)*weights2d(ip,nw)
   !  end do
      !
    end do
   !!
   !! Now at the volume conservation source term
   !! to each solution point in this cell
   !!
   !do k = n1,n2
   !  do m = 1,neqns
   !    residual(m,k) = residual(m,k) - src(m)/volj(k)
   !  end do
   !end do
    !
  end do
  !
end subroutine flux_divergence_cr
!
!###############################################################################
!
subroutine flux_divergence_jac_cr
  !
  !... this subroutine get derivative of fluxes (no interaction) by
  !... applying the chain rule on UTD
  !... for each cell, get uxi and uet at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp2dt,n_sp2dq, &
                      pts_in_cell_ptr
  use metvar,  only : dx,dy,dr_mat,ds_mat, &
                      jacobian,dJ
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  integer :: k,m,n,nw,n1,n2,np,ip
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns,1:neqns)      :: dfdu,dgdu
  real(wp), dimension(1:neqns,1:neqns)      :: c_f,c_g
  real(wp), dimension(1:neqns,1:n_sp2dq) :: ucr,ucs
  real(wp), dimension(1:n_sp2dq,1:neqns) :: uja
  real(wp), dimension(1:neqns,1:n_sp2dq) :: ujar,ujas
  !
continue
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1
    !
    nw = merge(1,2,np == n_sp2dt) ! check if cell n is tria or quad
    !
    ! Compute u_tilde (= conserved variables * jacobian) at each solution point
    !
    do k = n1,n2
      uja(k-n1+1,:) = usp(:,k)*jacobian(k)
    end do
    !
    ! Compute the derivative of u_tilde at each solution point
    !
    do k = 1,np
      do m = 1,neqns
        ujar(m,k) = dot_product( dr_mat(1:np,k,nw) , uja(1:np,m) )
        ujas(m,k) = dot_product( ds_mat(1:np,k,nw) , uja(1:np,m) )
      end do
    end do
    !
    ! Compute the derivatives of the fluxes at each solution point
    !
    do k = n1,n2
      !
      ip = k-n1+1 ! local index of current solution point
      !
      ! Compute the derivative of the conserved
      ! variables for the current sol-pt
      !
      ! = 1/J * { [d(u*J)/dr] - [u * d(J)/dr] }
      ! = 1/J * { [J * d(u)/dr] + [u * d(J)/dr] - [u * d(J)/dr] }
      ! = 1/J * { [J * d(u)/dr] }
      ! = d(u)/dr
      !
      ucr(:,ip) = (one/jacobian(k)) * (ujar(:,ip) - usp(:,k)*dJ(k)%dr)
      ucs(:,ip) = (one/jacobian(k)) * (ujas(:,ip) - usp(:,k)*dJ(k)%ds)
      !
      ! Compute the Jacobian of the fluxes at the current sol-pt
      !
      dfdu(:,:) = f_invs_flx_jacobian_cv( usp(:,k) )
      dgdu(:,:) = g_invs_flx_jacobian_cv( usp(:,k) )
      !
      ! Compute the Jacobian in the direction of the face normal
      !
      c_f(:,:) = dy(k)%ds*dfdu(:,:) - dx(k)%ds*dgdu(:,:)
      c_g(:,:) = dy(k)%dr*dfdu(:,:) - dx(k)%dr*dgdu(:,:)
      !
      ! Use the chain rule and add to the residual
      !
      do m = 1,neqns
        residual(m,k) = dot_product(ucr(:,ip), c_f(:,m)) - &
                        dot_product(ucs(:,ip), c_g(:,m))
      end do
      !
    end do
    !
  end do
  !
end subroutine flux_divergence_jac_cr
!
!###############################################################################
!
subroutine dg_strong_form_old()
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : n_sp1de,n_sp2dt,n_sp2dq, &
                      neqns,ncell,nface,nfbnd, &
                      pts_in_cell_ptr,cells_on_face
  use geovar,  only : interp_s2f
  use metvar,  only : cormat_dg,geofa, &
                      dx,dy,dr_mat,ds_mat
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  integer   :: k,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,nl,nr,nw,n1,n2,np,wl,wr
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr
  real(wp), dimension(1:neqns) :: fjl,fjr
  real(wp), dimension(1:neqns) :: fnc
  real(wp), dimension(1:neqns) :: ff,gg
  real(wp), dimension(1:n_sp2dq,1:neqns) :: ftdc
  real(wp), dimension(1:n_sp2dq,1:neqns) :: gtdc
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxg
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxg
  !
continue
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    ftdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    gtdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    !
    ! Get the flux at each solution point in the current cell
    !
    do k = n1,n2
      !
      ff(:) = f_invs_flx_cv( usp(:,k) ) ! x-directional flux
      gg(:) = g_invs_flx_cv( usp(:,k) ) ! y-directional flux
      !
      ftdc(k-n1+1,:) =  ff(:)*dy(k)%ds - gg(:)*dx(k)%ds
      gtdc(k-n1+1,:) =  ff(:)*dy(k)%dr - gg(:)*dx(k)%dr
      !
    end do
    !
    ! Compute the derivatives of the fluxes at each
    ! solution point and add to the residual
    !
    do k = 1,np
      do m = 1,neqns
        residual(m,k+n1-1) = dot_product(dr_mat(1:np,k,nw),ftdc(1:np,m)) - &
                             dot_product(ds_mat(1:np,k,nw),gtdc(1:np,m))
      end do
    end do
    !
  end do
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Compute the fluxes at each solution point and store the
    ! transpose of usp for the left cell to reduce cache misses
    !
    do k = l1,l2
      lflxf(k-l1+1,:) = f_invs_flx_cv( usp(:,k) )
      lflxg(k-l1+1,:) = g_invs_flx_cv( usp(:,k) )
      lusp (k-l1+1,:) = usp(:,k)
    end do
    !
    ! Compute the fluxes at each solution point and store the
    ! transpose of usp for the right cell to reduce cache misses
    ! 
    do k = r1,r2
      rflxf(k-r1+1,:) = f_invs_flx_cv( usp(:,k) )
      rflxg(k-r1+1,:) = g_invs_flx_cv( usp(:,k) )
      rusp (k-r1+1,:) = usp(:,k)
    end do
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the conservative variables from
      ! the left and right cells to the current flux point on the interface
      ! in order to compute the numerical flux at the current flux point
      !
      do m = 1,neqns
        ul(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lusp(1:nl,m)) ! left  state
        ur(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rusp(1:nr,m)) ! right state
      end do
      !
      ! Get the numerical flux at the current flux point
      !
      if (invs_flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , ul , ur )
      else if (invs_flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , ul , ur )
      else if (invs_flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , ul , ur )
      end if
      !
      ! Get left and right discontinuous fluxes at the current flux point
      !
      do m = 1,neqns
        dfl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxf(1:nl,m))
        dgl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxg(1:nl,m))
        dfr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxf(1:nr,m))
        dgr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxg(1:nr,m))
      end do
      !
      ! Get the discontinuous flux function at the current
      ! interface point in the direction of the face normal
      !
      fnl(:) = dfl(:)*unrm(1) + dgl(:)*unrm(2)
      fnr(:) = dfr(:)*unrm(1) + dgr(:)*unrm(2)
      !
      ! The jump in the flux at the current interface point
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m)*cormat_dg(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m)*cormat_dg(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine dg_strong_form_old
!
!###############################################################################
!
subroutine dg_strong_form()
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : n_sp1de,n_sp2dt,n_sp2dq, &
                      neqns,ncell,nface,nfbnd, &
                      pts_in_cell_ptr,cells_on_face
  use geovar,  only : n_qp1de,n_qp2dq,n_qp2dt
  use geovar,  only : interp_f2g,interp_s2q,qdtwts_2d
  use geovar,  only : interp_s2f,interp_s2g
  use metvar,  only : cormat_dg,geofa,dqr_mat,dqs_mat, &
                      dx,dy,dr_mat,ds_mat,jacobian
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  integer   :: k,m,n,nf,r1,r2,l1,l2,nsp,nqp
  integer   :: jl,jr,nl,nr,nw,n1,n2,np,wl,wr
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr
  real(wp), dimension(1:neqns) :: fjl,fjr
  real(wp), dimension(1:neqns) :: fnc
  real(wp), dimension(1:neqns) :: ff,gg
  real(wp), dimension(1:n_qp2dq) :: dxdr
  real(wp), dimension(1:n_qp2dq) :: dxds
  real(wp), dimension(1:n_qp2dq) :: dydr
  real(wp), dimension(1:n_qp2dq) :: dyds
  real(wp), dimension(1:neqns,1:n_qp2dq) :: uqp
  real(wp), dimension(1:n_qp2dq,1:neqns) :: ftdc
  real(wp), dimension(1:n_qp2dq,1:neqns) :: gtdc
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxg
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxg
  real(wp), dimension(1:n_qp2dq) :: qjac
  real(wp), dimension(1:2,1:n_qp1de) :: facenorms
  !
continue
  !
  nsp = n_sp1de
  nqp = n_qp1de
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    do k = 1,n_qp2dq
      !
      ! Get the conservative solution at the quadrature points
      !
      do m = 1,neqns
        uqp(m,k) = dot_product( interp_s2q(1:np,k,nw) , usp(m,n1:n2) )
      end do
      !
      ! Get the metrics for this cell
      !
      dxdr(k) = dot_product( interp_s2q(1:np,k,nw) , dx(n1:n2)%dr )
      dxds(k) = dot_product( interp_s2q(1:np,k,nw) , dx(n1:n2)%ds )
      dydr(k) = dot_product( interp_s2q(1:np,k,nw) , dy(n1:n2)%dr )
      dyds(k) = dot_product( interp_s2q(1:np,k,nw) , dy(n1:n2)%ds )
      qjac(k) = dot_product( interp_s2q(1:np,k,nw) , jacobian(n1:n2) )
      !
    end do
    !
    !
    ftdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    gtdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    !
    ! Get the flux at each solution point in the current cell
    !
    do k = 1,n_qp2dq
      !
      ff(:) = f_invs_flx_cv( uqp(:,k) ) ! x-directional flux
      gg(:) = g_invs_flx_cv( uqp(:,k) ) ! y-directional flux
      !
      ftdc(k,:) =  ff(:)*dyds(k) - gg(:)*dxds(k)
      gtdc(k,:) =  ff(:)*dydr(k) - gg(:)*dxdr(k)
      !
      ftdc(k,:) = ftdc(k,:) * qdtwts_2d(k,nw)
      gtdc(k,:) = gtdc(k,:) * qdtwts_2d(k,nw)
      !
    end do
    !
    ! Compute the derivatives of the fluxes at each
    ! solution point and add to the residual
    !
    do k = 1,np
      do m = 1,neqns
        residual(m,k+n1-1) = dot_product(dqr_mat(:,k,nw),ftdc(:,m)) - &
                             dot_product(dqs_mat(:,k,nw),gtdc(:,m))
      end do
    end do
    !
  end do
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Compute the fluxes at each solution point and store the
    ! transpose of usp for the left cell to reduce cache misses
    !
    do k = l1,l2
      lflxf(k-l1+1,:) = f_invs_flx_cv( usp(:,k) )
      lflxg(k-l1+1,:) = g_invs_flx_cv( usp(:,k) )
      lusp (k-l1+1,:) = usp(:,k)
    end do
    !
    ! Compute the fluxes at each solution point and store the
    ! transpose of usp for the right cell to reduce cache misses
    ! 
    do k = r1,r2
      rflxf(k-r1+1,:) = f_invs_flx_cv( usp(:,k) )
      rflxg(k-r1+1,:) = g_invs_flx_cv( usp(:,k) )
      rusp (k-r1+1,:) = usp(:,k)
    end do
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_qp1de
      do n = 1,2
        facenorms(n,k) = dot_product( interp_f2g(1:nsp,k) , geofa(n,1:nsp,nf) )
      end do
    end do
    !
   !do k = 1,n_sp1de
    do k = 1,n_qp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( facenorms(:,k) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( facenorms(:,k) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the conservative variables from
      ! the left and right cells to the current flux point on the interface
      ! in order to compute the numerical flux at the current flux point
      !
      do m = 1,neqns
       !ul(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lusp(1:nl,m)) ! left  state
        ul(m) = dot_product(interp_s2g(1:nl,kcorl,wl),lusp(1:nl,m)) ! left  state
       !ur(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rusp(1:nr,m)) ! right state
        ur(m) = dot_product(interp_s2g(1:nr,kcorr,wr),rusp(1:nr,m)) ! right state
      end do
      !
      ! Get the numerical flux at the current flux point
      !
      if (invs_flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , ul , ur )
      else if (invs_flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , ul , ur )
      else if (invs_flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , ul , ur )
      end if
      !
      ! Get left and right discontinuous fluxes at the current flux point
      !
      do m = 1,neqns
       !dfl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxf(1:nl,m))
        dfl(m) = dot_product(interp_s2g(1:nl,kcorl,wl),lflxf(1:nl,m))
       !dgl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxg(1:nl,m))
        dgl(m) = dot_product(interp_s2g(1:nl,kcorl,wl),lflxg(1:nl,m))
       !dfr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxf(1:nr,m))
        dfr(m) = dot_product(interp_s2g(1:nr,kcorr,wr),rflxf(1:nr,m))
       !dgr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxg(1:nr,m))
        dgr(m) = dot_product(interp_s2g(1:nr,kcorr,wr),rflxg(1:nr,m))
      end do
      !
      ! Get the discontinuous flux function at the current
      ! interface point in the direction of the face normal
      !
      fnl(:) = dfl(:)*unrm(1) + dgl(:)*unrm(2)
      fnr(:) = dfr(:)*unrm(1) + dgr(:)*unrm(2)
      !
      ! The jump in the flux at the current interface point
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m)*cormat_dg(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m)*cormat_dg(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine dg_strong_form
!
!###############################################################################
!
subroutine dg_weak_form()
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : n_sp1de,n_sp2dt,n_sp2dq, &
                      neqns,ncell,nface,nfbnd,bface, &
                      cells_on_face,pts_in_cell_ptr
  use geovar,  only : interp_s2f
  use metvar,  only : cormat_dg,geofa, &
                      dx,dy,dr_hat,ds_hat
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  logical   :: is_1d_bc
  integer   :: k,l,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,nl,nr,nw,n1,n2,np,wl,wr
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)            :: unrm
  real(wp), dimension(1:neqns)           :: fnc
  real(wp), dimension(1:neqns)           :: ul,ur
  real(wp), dimension(1:neqns)           :: ff,gg
  real(wp), dimension(1:neqns)           :: fjl,fjr
  real(wp), dimension(1:n_sp2dq,1:neqns) :: ftdc,gtdc
  !
continue
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    ftdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    gtdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    !
    ! Get the flux at each solution point in the current cell
    !
    do k = n1,n2
      !
      ff(:) = f_invs_flx_cv( usp(:,k) ) ! x-directional flux
      gg(:) = g_invs_flx_cv( usp(:,k) ) ! y-directional flux
      !
      ftdc(k-n1+1,:) =  ff(:)*dy(k)%ds - gg(:)*dx(k)%ds
      gtdc(k-n1+1,:) =  ff(:)*dy(k)%dr - gg(:)*dx(k)%dr
      !
    end do
    !
    ! Compute the derivatives of the fluxes at each
    ! solution point and add to the residual
    !
    do k = 1,np
      do m = 1,neqns
        residual(m,k+n1-1) = dot_product(dr_hat(1:np,k,nw),ftdc(1:np,m)) - &
                             dot_product(ds_hat(1:np,k,nw),gtdc(1:np,m))
      end do
    end do
    !
  end do
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    ! Check if this is a boundary face needing a 1D treatment of the flux points
    !
    is_1d_bc = fals
    if (nf <= nfbnd) is_1d_bc = .not. any(bface(1,nf) == bc_comm)
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the conservative variables from
      ! the left and right cells to the current flux point on the interface
      !
      do m = 1,neqns
        ul(m) = dot_product(interp_s2f(1:nl,kcorl,wl),usp(m,l1:l2)) ! left  state
        ur(m) = dot_product(interp_s2f(1:nr,kcorr,wr),usp(m,r1:r2)) ! right state
      end do
      !
      ! Get the numerical flux at the current flux point
      !
      if (invs_flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , ul , ur )
      else if (invs_flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , ul , ur )
      else if (invs_flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , ul , ur )
      end if
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fnc(m)*cormat_dg(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fnc(m)*cormat_dg(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine dg_weak_form
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
include "Functions/inviscid_fluxes.f90"
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Procedures for viscous fluxes
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function visc_flx_cv(uu,dcvdx) result(return_value) 
  !
  !.. Use Statements ..
  use geovar, only : nr,neqns
  use ovar,   only : cpref_nd,Pr
  !
  use generic_mod, only : viscosity_cv_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:neqns),      intent(in) :: uu
  real(wp), dimension(1:nr,1:neqns), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:neqns,1:nr) :: return_value
  !
  !.. Local Scalars ..
  integer  :: nq,vb,ve
  real(wp) :: dudx,dudy,dvdx,dvdy
  real(wp) :: mvsc,mcdt,dTdx,dTdy
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vel
  real(wp), dimension(1:nr,1:neqns) :: dpvdx
  !
continue
  !
  nq = neqns ! alias for neqns
  vb = 2     ! first momentum equation
  ve = nq-1  ! last  momentum equation 
  !
  ! Get the gradients of the primitive variables using the 
  ! gradients of the conservative variables and the product rule
  !
  dpvdx(1:nr,1:nq) = grad_cv_to_grad_pv_sp(uu,dcvdx)
  !
  ! Get the velocity vector
  !
  vel(:) = uu(vb:ve)/uu(1)
  !
  mvsc = viscosity_cv_sp(uu)
  mcdt = cpref_nd * mvsc / Pr
  !
  dudx = dpvdx(1,2)
  dudy = dpvdx(2,2)
  dvdx = dpvdx(1,3)
  dvdy = dpvdx(2,3)
  dTdx = dpvdx(1,nq)
  dTdy = dpvdx(2,nq)
  !
  ! x-directional viscous fluxes
  !
  return_value( 1,1) = zero
  return_value( 2,1) = mvsc * two3 * (two*dudx - dvdy)
  return_value( 3,1) = mvsc * (dudy + dvdx)
  return_value(nq,1) = mcdt*dTdx + dot_product( vel , return_value(vb:ve,1) )
  !
  ! y-directional viscous fluxes
  !
  return_value( 1,2) = zero
  return_value( 2,2) = mvsc * (dudy + dvdx)
  return_value( 3,2) = mvsc * two3 * (two*dvdy - dudx)
  return_value(nq,2) = mcdt*dTdy + dot_product( vel , return_value(vb:ve,2) )
  !
end function visc_flx_cv
!
!###############################################################################
!
pure function normal_visc_flx_cv(rnrm,uu,dcvdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : neqns,nr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:neqns),      intent(in) :: uu
  real(wp), dimension(1:nr),         intent(in) :: rnrm
  real(wp), dimension(1:nr,1:neqns), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:neqns) :: return_value
  !
  !.. Local Scalars ..
  integer :: l
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)         :: unrm
  real(wp), dimension(1:neqns,1:nr) :: vflx
  !
continue
  !
  unrm(:) = unit_vector( rnrm(:) )
  !
  vflx(1:neqns,1:nr) = visc_flx_cv(uu,dcvdx)
  !
  return_value(1:neqns) = unrm(1)*vflx(1:neqns,1)
  do l = 2,nr
    return_value(1:neqns) = return_value(1:neqns) + unrm(l)*vflx(1:neqns,l)
  end do
  !
end function normal_visc_flx_cv
!
!###############################################################################
!
pure function visc_wall_flx_cv(unrm,uwall,dpvdx) result(return_value) 
  !
  ! UNRM  : Outward pointing unit normal to the wall
  !         boundary face for current flux point
  ! UWALL : Conservative variables exact solution at the wall
  ! DPVDX : Prescribed wall gradients of the primitive variables at the
  !         current flux point on the wall boundary face. Most gradients
  !         are extrapolated from the interior cell flux point on the
  !         boundary except for certain special cases, e.g., an adiabatic
  !         wall where dT/dx = dT/dy = 0.
  !
  !.. Use Statements ..
  use geovar, only : nr,neqns
  use ovar,   only : cpref_nd,Pr
  !
  use generic_mod, only : viscosity_cv_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),         intent(in) :: unrm
  real(wp), dimension(1:neqns),      intent(in) :: uwall
  real(wp), dimension(1:nr,1:neqns), intent(in) :: dpvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:neqns) :: return_value
  !
  !.. Local Scalars ..
  integer  :: nq,vb,ve
  real(wp) :: dudx,dudy,dvdx,dvdy
  real(wp) :: mvsc,mcdt,dTdx,dTdy
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)    :: vel
  real(wp), dimension(1:neqns) :: xflx,yflx
  !
continue
  !
  ! For slip wall: all gradients /= zero
  !                uwall /= zero
  !
  ! For adiabatic wall: Temperature gradient == zero
  !                     Velocity gradients /= zero
  !                     uwall(2:neqns-1) [i.e., wall velocities] == zero
  !
  ! For isothermal wall: all gradients /= zero
  !                      uwall(2:neqns-1) [i.e., wall velocities] == zero
  !                
  !
  nq = neqns ! alias for neqns
  vb = 2     ! first momentum equation
  ve = nq-1  ! last  momentum equation 
  !
  ! Get the velocity vector from the wall boundary conditions
  !
  vel(:) = uwall(vb:ve)/uwall(1)
  !
  ! Compute the molecular viscosity and
  ! conductivity using the wall temperature
  !
  mvsc = viscosity_cv_sp(uwall)
  mcdt = cpref_nd * mvsc / Pr
  !
  ! Load the necessary gradient variables
  !
  dudx = dpvdx(1,2)
  dudy = dpvdx(2,2)
  dvdx = dpvdx(1,3)
  dvdy = dpvdx(2,3)
  dTdx = dpvdx(1,nq)
  dTdy = dpvdx(2,nq)
  !
  ! x-directional viscous fluxes
  !
  xflx( 1) = zero
  xflx( 2) = mvsc * two3 * (two*dudx - dvdy)
  xflx( 3) = mvsc * (dudy + dvdx)
  xflx(nq) = mcdt*dTdx + dot_product( vel , xflx(vb:ve) )
  !
  ! y-directional viscous fluxes
  !
  yflx( 1) = zero
  yflx( 2) = mvsc * (dudy + dvdx)
  yflx( 3) = mvsc * two3 * (two*dvdy - dudx)
  yflx(nq) = mcdt*dTdy + dot_product( vel , yflx(vb:ve) )
  !
  return_value(1:nq) = xflx(1:nq)*unrm(1) + yflx(1:nq)*unrm(2)
  !
end function visc_wall_flx_cv
!
!###############################################################################
!
pure function common_visc_flx_cv(unrm,ul,duldx,ur,durdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr,neqns
  use ovar,   only : visc_variables_to_ave
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),         intent(in) :: unrm
  real(wp), dimension(1:neqns),      intent(in) :: ul,ur
  real(wp), dimension(1:nr,1:neqns), intent(in) :: duldx,durdx
  !
  !.. Function Result ..
  real(wp), dimension(1:neqns) :: return_value
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: vjump,qjump,penalty
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns)      :: ucom
  real(wp), dimension(1:nr,1:neqns) :: qcom
  real(wp), dimension(1:neqns,1:nr) :: lvflx
  real(wp), dimension(1:neqns,1:nr) :: rvflx
  real(wp), dimension(1:nr,1:neqns) :: vflux
  real(wp), dimension(1:neqns,1:nr) :: fvcom
  !
continue
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    ! Get the left and right viscous fluxes
    !
    lvflx(1:neqns,1:nr) = visc_flx_cv(ul,duldx)
    rvflx(1:neqns,1:nr) = visc_flx_cv(ur,durdx)
    !
    ! Use the LDG method to compute the common viscous flux vector
    !
    ! Start by applying the averaging operator {{.}} to the viscous flux vectors
    !
    do m = 1,neqns
      do l = 1,nr
        vflux(l,m) = half * (lvflx(m,l) + rvflx(m,l))
      end do
    end do
    !
    ! Apply the differencing (jump) operator [[.]] to the interface viscous
    ! flux vectors and then add it to the common viscous flux vector
    !
    do m = 1,neqns
      !
      vjump = dot_product( unrm(1:nr) , lvflx(m,1:nr)-rvflx(m,1:nr) )
      !
      do l = 1,nr
        vflux(l,m) = vflux(l,m) + LDG_beta(l)*vjump
      end do
      !
    end do
    !
    ! Finally, add the penalty term, (simply the differencing (jump) operator
    ! [[.]] applied to the interface solutions) to the common viscous flux vector
    !
    do m = 1,neqns
      penalty = LDG_tau * ( ul(m) - ur(m) )
      do l = 1,nr
        vflux(l,m) = vflux(l,m) + penalty * unrm(l)
      end do
    end do
    !
    ! Finally, dot the common viscous flux vector with the face unit normal
    !
    do m = 1,neqns
      return_value(m) = dot_product( vflux(1:nr,m) , unrm(1:nr) )
    end do
    !
  else if (visc_variables_to_ave == visc_solution_and_gradients) then
    !
    ! Get the common solution at the interface using a simple average
    !
    do m = 1,neqns
      ucom(m) = half * ( ul(m) + ur(m) )
    end do
    !
    ! Use the LDG method to compute the common gradient
    !
    ! Start by applying the averaging operator {{.}} to the gradients
    !
    do m = 1,neqns
      do l = 1,nr
        qcom(l,m) = half * (duldx(l,m) + durdx(l,m))
      end do
    end do
    !
    ! Apply the differencing (jump) operator [[.]] to the
    ! gradients and then add it to the common gradient
    !
    do m = 1,neqns
      !
      qjump = dot_product( unrm(1:nr) , duldx(1:nr,m)-durdx(1:nr,m) )
      !
      do l = 1,nr
        qcom(l,m) = qcom(l,m) + LDG_beta(l)*qjump
      end do
      !
    end do
    !
    ! Finally, add the penalty term, (simply the differencing (jump) operator
    ! [[.]] applied to the interface solutions) to the common gradient
    !
    do m = 1,neqns
      penalty = LDG_tau * ( ul(m) - ur(m) )
      do l = 1,nr
        qcom(l,m) = qcom(l,m) + penalty * unrm(l)
      end do
    end do
    !
    ! Compute the common viscous flux vector using
    ! the common solution and common gradient
    !
    fvcom(1:neqns,1:nr) = visc_flx_cv(ucom,qcom)
    !
    ! Finally, dot the common viscous flux vector with the face unit normal
    !
    do m = 1,neqns
      return_value(m) = dot_product( fvcom(m,1:nr) , unrm(1:nr) )
    end do
    !
  end if
  !
end function common_visc_flx_cv
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
  use geovar,  only : neqns,ncell,n_sp2dt,n_sp2dq
  use geovar,  only : pts_in_cell_ptr
  use metvar,  only : dr_mat,ds_mat,dx,dy,jacobian
  use flowvar, only : usp,dusp
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,n1,n2,np,nw
  real(wp) :: drdx,drdy,dsdx,dsdy
  real(wp) :: dudr,duds,dudx,dudy
  real(wp) :: jac_inv
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq,1:neqns) :: tusp
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for solpts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for solpts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge(1,2,np == n_sp2dt) ! check if cell n is tria or quad
    !
    ! Store the transpose of usp for this cell to reduce cache misses
    !
    do k = 1,np
      do m = 1,neqns
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    do k = 1,np
      !
      l = k+n1-1 ! total index of solution point on current partition
      !
      jac_inv = one/jacobian(l) ! inverse of det(J)
      !
      ! Compute the contravariant metrics
      !
      drdx =  dy(l)%ds
      drdy = -dx(l)%ds
      dsdx = -dy(l)%dr
      dsdy =  dx(l)%dr
     !drdx =  jac_inv * dy(l)%ds
     !drdy = -jac_inv * dx(l)%ds
     !dsdx = -jac_inv * dy(l)%dr
     !dsdy =  jac_inv * dx(l)%dr
      !
      do m = 1,neqns
        !
        ! Compute the derivatives in r-s space
        !
        dudr = dot_product( dr_mat(1:np,k,nw) , tusp(1:np,m) )
        duds = dot_product( ds_mat(1:np,k,nw) , tusp(1:np,m) )
        !
        ! Transform the derivatives back to x-y space
        !
        dudx = drdx*dudr + dsdx*duds
        dudy = drdy*dudr + dsdy*duds
        !
        dusp(1,m,l) = dudx*jac_inv
        dusp(2,m,l) = dudy*jac_inv
        !
      end do
      !
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
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_sp2dt,nr,n_solpts,ncell
  use geovar,  only : cells_on_face,pts_in_cell_ptr,bface
  use metvar,  only : cormat,geofa,jacobian,volj
  use flowvar, only : dusp,faceusp
  use ovar,    only : walls_are_exact
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nf,r1,r2,l1,l2,n1,n2,nc
  integer  :: jl,jr,wl,wr,nsp
  integer  :: kbl,kbr,kl,kr
  real(wp) :: cbias,hfpb,hfmb
  real(wp) :: correct_n,correct_m,ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)              :: unrm
  real(wp), dimension(1:neqns)           :: ljfp,rjfp
  real(wp), dimension(1:neqns,1:n_sp1de) :: lufp,rufp
  !
continue
  !
  nsp = n_sp1de
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
    jl = cells_on_face(1,nf) ! left  cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    wl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if jl is a tria or quad
    !
    ! For all the flux points on the current face, extract the conservative
    ! conservative variables for the left and right cells from faceusp
    !
    lufp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,1,nf)
    rufp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,2,nf)
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! Use the LDG flux formulation for computing the common interface solutions
    !   u_com = half*(u_left + u_right) - &
    !                (u_left - u_right)*dot_product(LDG_beta,unit_normal)
    !
    do k = 1,nsp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = kbl-1 + k ! index for left cell
      !
      ! Get the unit normal at this flux point
      !
      unrm(1:nr) = unit_vector( geofa(1:nr,k,nf) )
      !
      ! Get the face length for the current face
      !
      dsl = cl(wl) * rnorm2( geofa(1:nr,k,nf) )
      !
      ! Compute the LDG directional parameters
      ! NOTE: cbias is used for all boundary faces that allow for flow
      !       across the interface, e.g, inflow/outflow, cpu/partition
      !       boundary, or periodic boundary. These boundaries are
      !       treated the same as an interior face except that we dont
      !       need to worry about correcting the right cell.
      !
      cbias = dot_product( LDG_beta(1:nr) , unrm(1:nr) )
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
      do m = 1,neqns
        ljfp(m) = hfpb*(rufp(m,k) - lufp(m,k))
      end do
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      do n = l1,l2
        correct_n = dsl * cormat(n-l1+1,kl,wl) / jacobian(n)
        do m = 1,neqns
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
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! ending    flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    wl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , r2-r1+1 == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! For all the flux points on the current face, extract the conservative
    ! conservative variables for the left and right cells from faceusp
    !
    lufp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,1,nf)
    rufp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,2,nf)
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! Use the LDG flux formulation for computing the common interface solutions
    !   u_com = half*(u_left + u_right) - &
    !                (u_left - u_right)*dot_product(LDG_beta,unit_normal)
    !
    do k = 1,nsp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = kbl-1 + k ! index for left  cell
      kr = kbr+1 - k ! index for right cell
      !
      ! Get the unit normal at this flux point
      !
      unrm(1:nr) = unit_vector( geofa(1:nr,k,nf) )
      !
      ! Get the face length for the current face
      !
      ds = rnorm2( geofa(1:nr,k,nf) )
      dsl = cl(wl)*ds
      dsr = cl(wr)*ds
      !
      ! Compute the LDG directional parameters
      !
      cbias = dot_product( LDG_beta(1:nr) , unrm(1:nr) )
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
      do m = 1,neqns
        ljfp(m) = hfpb*(rufp(m,k) - lufp(m,k))
        rjfp(m) = hfmb*(lufp(m,k) - rufp(m,k))
      end do
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      do n = l1,l2
        correct_n = dsl * cormat(n-l1+1,kl,wl) / jacobian(n)
        do m = 1,neqns
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
      do n = r1,r2
        correct_n = dsr * cormat(n-r1+1,kr,wr) / jacobian(n)
        do m = 1,neqns
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
  use geovar,  only : nr,neqns,ncell,n_sp2dt,n_sp2dq,pts_in_cell_ptr
  use geovar,  only : flxpts_of_cell_ptr,flxpts_of_cell,n_order
  use geovar,  only : interp_s2f,n_sp1de,nfbnd,bface
  use metvar,  only : geofa,cormat,jacobian,volj
  use flowvar, only : dusp,facedusp,faceusp
  use geovar,  only : cells_on_face
  use ovar,    only : walls_are_exact
  !
  use parallel_mod, only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,wl,k,l,m
  integer :: nf,f1,f2,kf,lf,sf
  integer :: ncf,fp1,fp2,kfp,kbl,kl
  real(wp) :: dsl,cbias,hfpb,correct_n,correct_m
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:neqns) :: lufp,rufp,ujump
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_sp2dq,1:nr,1:neqns) :: face_dusp
  real(wp), dimension(1:nr,1:neqns,1:n_sp2dq) :: cell_dusp
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  cell_loop: do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    np = n2-n1+1
    wl = merge(1,2,np == n_sp2dt)
    !
    f1 = flxpts_of_cell_ptr(nc)+1
    f2 = flxpts_of_cell_ptr(nc+1)
    !
    ! Create transpose of dusp array to reduce cache misses
    !
    do n = 1,np
      do m = 1,neqns
        do l = 1,nr
          cell_dusp(l,m,n) = dusp(l,m,n+n1-1)
        end do
      end do
    end do
    !
    ! Now Interpolate the interior solution gradients to the flux points
    ! NOTE: I believe that the ordering of the flux points from f1 to f2 is
    !       such that all flux points on a given face are grouped together.
    !       Basically:
    !   FPs for face A are (f1+0*n_sp1de):(f1-1+1*n_sp1de)
    !   FPs for face B are (f1+1*n_sp1de):(f1-1+2*n_sp1de)
    !   FPs for face C are (f1+2*n_sp1de):(f1-1+3*n_sp1de)
    !   FPs for face D are (f1+3*n_sp1de):(f1-1+4*n_sp1de) (quad cells only)
    !
    cell_faces: do ncf = 1,merge(3,4,np==n_sp2dt) ! number of faces for this cell
      !
      fp1 = f1 + (ncf-1)*n_sp1de
      fp2 = f1 + ncf*n_sp1de - 1
      !
      ! Reset face_dusp to the original discontinuous gradient
      !
      do n = 1,np
        do m = 1,neqns
          do l = 1,nr
            face_dusp(n,l,m) = dusp(l,m,n+n1-1)
          end do
        end do
      end do
      !
      face_points: do kfp = fp1,fp2
        !
        nf = flxpts_of_cell(1,kfp) ! current face number for the current flux point
        kf = flxpts_of_cell(2,kfp) ! FP index local to the face
        lf = flxpts_of_cell(3,kfp) ! ideally this should evaluate to 1:f2-f1+1 in this loop
        sf = flxpts_of_cell(4,kfp) ! side of the face for this cell (1=left,2=right)
        !
        ! Solution at this flux point for the left and right cells
        !
        lufp(:) = faceusp(:,kf,1,nf)
        rufp(:) = faceusp(:,kf,2,nf)
        !
        ! Get the unit normal at this flux point
        !
        unrm(1:nr) = unit_vector( geofa(1:nr,kf,nf) )
        !
        ! Get the face length for the current face
        !
        dsl = cl(wl) * rnorm2( geofa(1:nr,kf,nf) )
        !
        ! Compute the directional bias coefficient. The merge term accounts for
        ! when the current cell is the right cell. In both cases, the jump
        ! equations should reduce to:
        !    jump_left  = (cbias + half)*(u_right - u_left)
        !    jump_right = (cbias - half)*(u_right - u_left)
        ! NOTE: Boundary faces should never be the right cell on the interface.
        !       For non-wall boundaries, this should reduce to jump_left above.
        !       For wall boundaries, jump_left should reduce to:
        !         jump_left = (u_right - u_left) where u_right = u_wall
        !
        cbias = dot_product( LDG_beta(1:nr) , unrm(1:nr) )
        !
        ! Reset cbias if this is a wall boundary
        !
        if (nf <= nfbnd .and. walls_are_exact) then
          if (any(bface(1,nf) == bc_walls)) cbias = half
        end if
        !
        hfpb = cbias + half
        !
        ! Get the interface jumps to correct the gradients
        !
        do m = 1,neqns
          ujump(m) = hfpb*(rufp(m) - lufp(m))
        end do
        !
        ! Apply the jump correction to the interior discontinuous gradients
        !
        do n = 1,np
          correct_n = dsl * cormat(n,lf,wl) / jacobian(n+n1-1)
          do m = 1,neqns
            correct_m = correct_n * ujump(m)
            do l = 1,nr
              face_dusp(n,l,m) = face_dusp(n,l,m) + correct_m * unrm(l)
              cell_dusp(l,m,n) = cell_dusp(l,m,n) + correct_m * unrm(l)
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
      interpolate_face: do kfp = fp1,fp2
        !
        nf = flxpts_of_cell(1,kfp)
        kf = flxpts_of_cell(2,kfp)
        lf = flxpts_of_cell(3,kfp)
        sf = flxpts_of_cell(4,kfp)
        !
        lmat(1:np) = interp_s2f(1:np,lf,wl)
        !
        do m = 1,neqns
          do l = 1,nr
            facedusp(l,m,kf,sf,nf) = dot_product(lmat(1:np),face_dusp(1:np,l,m))
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
    do n = 1,np
      do m = 1,neqns
        do l = 1,nr
          dusp(l,m,n+n1-1) = cell_dusp(l,m,n)
        end do
      end do
    end do
    !
  end do cell_loop
  !
  ! Exchange the left and right states for faces on communication boundaries
  !
  if (ncpu > 1) call exchange_facedusp(facedusp)
  !
end subroutine correct_gradients_br2
!
!###############################################################################
!
subroutine test_face_connectivity
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,ncell,nface,n_sp1de,bface,n_sp2dt
  use geovar, only : flxpts_of_cell_ptr,flxpts_of_cell
  use geovar, only : cells_on_face,nodes_of_cell_ptr
  use geovar, only : pts_in_cell_ptr
  !
  !.. Local Scalars ..
  integer :: nf,nc,ncf,jl,jr,kbl,kbr,k,kl,kr,n1,n2,np,f1,f2,fp1,fp2,kf,lf,sf,ierr,i
  integer :: kfp
  !
  !.. Local Arrays ..
  integer, dimension(1:2,1:n_sp1de,1:nface) :: face_loop
  integer, dimension(1:2,1:n_sp1de,1:nface) :: cell_loop
  !
  !.. Local Parameters ..
  integer, parameter :: iboth = 1011
  integer, parameter :: iface = 1111
  integer, parameter :: icell = 1211
  integer, parameter, dimension(1:2) :: side_switch = [1,-1]
  !
continue
  !
  face_loop(:,:,:) = 0
  cell_loop(:,:,:) = 0
  !
  ! Old style loop over nface using cells_on_face
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf)
    jr = cells_on_face(2,nf)
    !
    kbl = cells_on_face(5,nf)
    kbr = cells_on_face(6,nf)
    !
    do k = 1,n_sp1de
      !
      kl = kbl-1 + k
      kr = kbr+1 - k
      !
      face_loop(1,k,nf) = kl
      face_loop(2,k,nf) = kr
      !
    end do
    !
  end do
  !
  ! New style loop over ncell using flxpts_of_cell
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    np = n2-n1+1
    !
    f1 = flxpts_of_cell_ptr(nc)+1
    f2 = flxpts_of_cell_ptr(nc+1)
    !
    do ncf = 1,merge(3,4,np==n_sp2dt)
      !
      fp1 = f1 + (ncf-1)*n_sp1de
      fp2 = f1 + ncf*n_sp1de-1
      !
      do kfp = fp1,fp2
        !
        nf = flxpts_of_cell(1,kfp)
        kf = flxpts_of_cell(2,kfp)
        lf = flxpts_of_cell(3,kfp)
        sf = flxpts_of_cell(4,kfp)
        !
        kbl = cells_on_face(4+sf,nf)
        k = kfp-fp1+1
        kl = kbl - side_switch(sf) + side_switch(sf)*k
        !
        cell_loop(sf,k,nf) = kl
        !
        write (102,5) sf,ncf,k,kf,k==kf,kl,lf,kl==lf,kfp
        5 format ("sf=",i1,", ncf=",i1,";     k=",i2,", kf=",i2,", ",l1, &
                                       "     kl=",i2,", lf=",i2,", ",l1, &
                                       "     kfp=",i0)
      end do
      !
    end do
    !
  end do
  !
  open (iboth,file="face_vs_cell_loop.dat",status="replace",action="write", &
             position="rewind",form="formatted",iostat=ierr)
  call io_error("test_face_connectivity","face_vs_cell_loop.dat",1,ierr)
  !
  open (iface,file="face_loop.dat",status="replace",action="write", &
             position="rewind",form="formatted",iostat=ierr)
  call io_error("test_face_connectivity","face_loop.dat",1,ierr)
  !
  open (icell,file="cell_loop.dat",status="replace",action="write", &
             position="rewind",form="formatted",iostat=ierr)
  call io_error("test_face_connectivity","cell_loop.dat",1,ierr)
  !
  do nf = 1,nface
    !
    write (iboth,1) nf
    if (nf > nfbnd) then
      write (iface,2) nf
      write (icell,2) nf
    end if
    !
    do k = 1,n_sp1de
      write (iboth,3) k,(face_loop(i,k,nf),cell_loop(i,k,nf),i=1,2)
      if (nf > nfbnd) then
        write (iface,4) k,(face_loop(i,k,nf),i=1,2)
        write (icell,4) k,(cell_loop(i,k,nf),i=1,2)
      end if
    end do
    !
  end do
  !
  close (iboth,iostat=ierr)
  call io_error("test_face_connectivity","face_vs_cell_loop.dat",2,ierr)
  !
  close (iface,iostat=ierr)
  call io_error("test_face_connectivity","face_loop.dat",2,ierr)
  !
  close (icell,iostat=ierr)
  call io_error("test_face_connectivity","cell_loop.dat",2,ierr)
  !
  ! Kill execution after outputing the connectivities
  !
  call stop_recon(stop_mpi,"test_face_connectivity",&
                           "Finished test of face connectivity")
  !
  1 format ("----------------------------------------------",/, &
            "Face: ",i5,";    left_side         right_side",/, &
            "----------------------------------------------")
  2 format ("Face: ",i5)
  3 format ("  flxpt: ",i1,";  face=",i2,",cell=",i2,"   face=",i2,",cell=",i2)
  4 format ("  flxpt: ",i1,";  left=",i2,"  right=",i2)
  !
end subroutine test_face_connectivity
!
!###############################################################################
!
include "Functions/grad_cv_to_grad_pv.f90"
!
!###############################################################################
!
end module flux_mod
