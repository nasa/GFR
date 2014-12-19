module flux_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  !
  implicit none
  !
  private
  !
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
subroutine interp_cv_to_flxpts()
  !
  use order_mod,  only : n_sp2dt,n_sp2dq
  use order_mod,  only : n_sp1de,n_fp1de
  use geovar,  only : ncell,pts_in_cell_ptr
  use geovar,  only : flxpts_of_cell_ptr,flxpts_of_cell
  use geovar,  only : interp_s2f,interp_s2g
  use flowvar, only : usp,faceusp
  !
  use parallel_mod, only : exchange_faceusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,wlocl,k,m
  integer :: nf,f1,f2,kf,lf,sf
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_sp2dq,1:nq) :: tusp
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
    wlocl = merge(1,2,np == n_sp2dt)
    !
    f1 = flxpts_of_cell_ptr(nc)+1
    f2 = flxpts_of_cell_ptr(nc+1)
    !
    ! Create transpose of usp array to reduce cache misses
    !
    do k = 1,np
      do m = 1,nq
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
      lmat(1:np) = interp_mat(1:np,lf,wlocl)
      !
      do m = 1,nq
        faceusp(m,kf,sf,nf) = dot(lmat(1:np),tusp(1:np,m))
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
  use order_mod,  only : n_sp2dt,n_sp2dq
  use order_mod,  only : n_sp1de,n_fp1de
  use geovar,  only : nr,ncell,pts_in_cell_ptr
  use geovar,  only : flxpts_of_cell_ptr,flxpts_of_cell
  use geovar,  only : interp_s2f,interp_s2g
  use flowvar, only : dusp,facedusp
  !
  use parallel_mod, only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,wlocl,k,l,m
  integer :: nf,f1,f2,kf,lf,sf
  integer :: fp1,fp2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_sp2dq,1:nr,1:nq) :: tdusp
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
    wlocl = merge(1,2,np == n_sp2dt)
    !
    f1 = flxpts_of_cell_ptr(nc)+1
    f2 = flxpts_of_cell_ptr(nc+1)
    !
    ! Create transpose of dusp array to reduce cache misses
    !
    do k = 1,np
      do m = 1,nq
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
      lmat(1:np) = interp_mat(1:np,lf,wlocl)
      !
      do m = 1,nq
        do l = 1,nr
          facedusp(l,m,kf,sf,nf) = dot(lmat(1:np),tdusp(1:np,l,m))
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
  use order_mod,  only : n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,n_qp2dq,n_qp2dt
  use geovar,  only : nr,ncell,pts_in_cell_ptr,interp_s2q
  use geovar,  only : flxpts_of_cell,flxpts_of_cell_ptr,projct_q2s
  use geovar,  only : interp_s2f,interp_s2g,using_quadrature
  use metvar,  only : metrics
  use metvar,  only : dr_mat,ds_mat
  use metvar,  only : metrics,geofa
  use metvar,  only : dqr_LPmat,dqs_LPmat,jacobian
  use flowvar, only : usp,residual,faceflx,faceusp
  use ovar,    only : dflux_opt
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,wlocl,n1,n2,nsp,nqp
  integer :: f1,f2,nf,sf,kf,lf
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:nq)   :: ff,gg
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_qp2dq) :: dxdr,dxds
  real(wp), dimension(1:n_qp2dq) :: dydr,dyds
  real(wp), dimension(1:n_qp2dq) :: qjac
  !
  real(wp), dimension(1:n_sp2dq,1:nq) :: tusp
  real(wp), dimension(1:nq,1:n_qp2dq) :: uqp
  !
  real(wp), dimension(1:n_sp2dq,1:nq,1:nr) :: rsflx_sp
  real(wp), dimension(1:n_sp2dq,1:nq,1:nr) :: xyflx_sp
  real(wp), dimension(1:n_sp2dq,1:nq,1:nr) :: dflx_sp
  real(wp), dimension(1:n_qp2dq,1:nq,1:nr) :: rsflx_qp
  real(wp), dimension(1:n_qp2dq,1:nq,1:nr) :: xyflx_qp
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
    wlocl = merge( 1 , 2 , nsp == n_sp2dt ) ! check if cell n is tria or quad
    !
    f1 = flxpts_of_cell_ptr(n)+1
    f2 = flxpts_of_cell_ptr(n+1)
    !
    ! Store the transpose of usp to reduce cache misses
    !
    tusp(1:nsp,1:nq) = transpose(usp(1:nq,n1:n2))
    !
    ! Get the information local to this cell
    !
    if (using_quadrature) then
      !
      do k = 1,nqp
        !
        ! Get the conservative solution at the quadrature points
        !
        do m = 1,nq
          uqp(m,k) = dot( interp_s2q(1:nsp,k,wlocl) , tusp(1:nsp,m) )
        end do
        !
        ! Get the metrics for this cell
        !
        dxdr(k) = dot( interp_s2q(1:nsp,k,wlocl) , metrics(1,1,n1:n2) )
        dxds(k) = dot( interp_s2q(1:nsp,k,wlocl) , metrics(1,2,n1:n2) )
        dydr(k) = dot( interp_s2q(1:nsp,k,wlocl) , metrics(2,1,n1:n2) )
        dyds(k) = dot( interp_s2q(1:nsp,k,wlocl) , metrics(2,2,n1:n2) )
        qjac(k) = dot( interp_s2q(1:nsp,k,wlocl) , jacobian(n1:n2) )
        !
      end do
      !
    else
      !
      uqp(1:nq,1:nsp) = usp(1:nq,n1:n2)
      !
      dxdr(1:nsp) = metrics(1,1,n1:n2)
      dxds(1:nsp) = metrics(1,2,n1:n2)
      dydr(1:nsp) = metrics(2,1,n1:n2)
      dyds(1:nsp) = metrics(2,2,n1:n2)
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
      do l = 1,nr
        xyflx_qp(k,:,l) = invs_flx_cv( uqp(:,k) , l )
      end do
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
          do m = 1,nq
            do k = 1,nsp
              xyflx_sp(k,m,l) = dot( projct_q2s(1:nqp,k,wlocl) , xyflx_qp(1:nqp,m,l) )
              rsflx_sp(k,m,l) = dot( projct_q2s(1:nqp,k,wlocl) , rsflx_qp(1:nqp,m,l) )
            end do
          end do
        end do
        !
      else
        !
        xyflx_sp(1:nsp,1:nq,1:nr) = xyflx_qp(1:nsp,1:nq,1:nr)
        rsflx_sp(1:nsp,1:nq,1:nr) = rsflx_qp(1:nsp,1:nq,1:nr)
        !
      end if
      !
      ! Take the derivative of the projected fluxes
      !
      do m = 1,nq
        do k = 1,nsp
          dflx_sp(k,m,1) = dot( dr_mat(1:nsp,k,wlocl) , rsflx_sp(1:nsp,m,1) )
          dflx_sp(k,m,2) = dot( ds_mat(1:nsp,k,wlocl) , rsflx_sp(1:nsp,m,2) )
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
        lmat(1:nsp) = interp_s2f(1:nsp,lf,wlocl)
        !
        do l = 1,nr
          do m = 1,nq
            faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
                                  dot(lmat(1:nsp),xyflx_sp(1:nsp,m,l))*unrm(l)
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
        do m = 1,nq
          dflx_sp(k,m,1) = dot( dqr_LPmat(1:nqp,k,wlocl) , rsflx_qp(1:nqp,m,1) )
          dflx_sp(k,m,2) = dot( dqs_LPmat(1:nqp,k,wlocl) , rsflx_qp(1:nqp,m,2) )
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
     !  lmat(1:nsp) = interp_s2g(1:nsp,lf,wlocl)
     !  !
     !  do l = 1,nr
     !    do m = 1,nq
     !      faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
     !                            dot(lmat(1:nsp),xyflx_qp(1:nsp,m,l))*unrm(l)
     !    end do
     !  end do
     !  !
     !end do
      !
    end if
    !
    ! Add the derivatives of the flux to the residual
    !
    do m = 1,nq
      do k = 1,nsp
        residual(m,k+n1-1) = dflx_sp(k,m,1) + dflx_sp(k,m,2)
      end do
    end do
    !
   !do k = 1,nsp
   !  do m = 1,nq
   !   !residual(m,k+n1-1) = dot(dr_mat(1:nsp,k,wlocl),ftdc(1:nsp,m)) - &
   !   !                     dot(ds_mat(1:nsp,k,wlocl),gtdc(1:nsp,m))
   !   !residual(m,k+n1-1) = dot(dqr_mat(:,k,wlocl),ftdc(:,m)) - &
   !   !                     dot(dqs_mat(:,k,wlocl),gtdc(:,m))
   !    residual(m,k+n1-1) = dot(dr_mat(1:nsp,k,wlocl),dfsp(1:nsp,m)) - &
   !                         dot(ds_mat(1:nsp,k,wlocl),dgsp(1:nsp,m))
   !  end do
   !end do
    !
   !do m = 1,nq
   !  residual(m,n1:n2) = matmul(invM_mat(:,:,wlocl),residual(m,n1:n2))!/jacobian(n1:n2)
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
  use order_mod,  only : n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,n_qp1de
  use geovar,  only : nfbnd,nface,nr
  use geovar,  only : bface,cells_on_face,pts_in_cell_ptr,xyz,flxwts_1d
  use geovar,  only : interp_f2g,projct_g2f,using_quadrature
  use ovar,    only : flux_divergence,dflux_opt,itcur,nst
  use ovar,    only : invs_flux_method
  use metvar,  only : cormat,geofa,volj
  use flowvar, only : usp,residual,faceflx,faceusp
  use geovar,  only : cell_geom
  !
 !use parallel_mod, only : wait_faceusp
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nf,r1,r2,l1,l2
  integer  :: wlocl,wlocr,nsp,nqp
  integer  :: kbl,kbr,kl,kr,iu
  integer  :: left_cell,right_cell
  integer  :: left_geom,right_geom
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)    :: unrm
  real(wp), dimension(1:nq) :: ul,ur
  real(wp), dimension(1:nq) :: dfl,dfr
  real(wp), dimension(1:nq) :: dgl,dgr
  real(wp), dimension(1:nq) :: fnl,fnr,fnc
  real(wp), dimension(1:nq) :: srcl,srcr
  real(wp), dimension(1:nr,1:n_qp1de)    :: rnrm
  real(wp), dimension(1:nq,1:n_sp1de) :: fjl,fjr
  real(wp), dimension(1:nq,1:n_fp1de) :: lflx,rflx
  real(wp), dimension(1:nq,1:n_qp1de) :: luqp,ruqp,fncq
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
    left_cell = cells_on_face(1,nf) ! left  cell on face
    right_cell = cells_on_face(2,nf) ! right cell on face
    !
    left_geom = cell_geom(left_cell)
    right_geom = cell_geom(right_cell)
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! ending    flux point for right cell
    !
    l1 = pts_in_cell_ptr(left_cell)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(left_cell+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(right_cell)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(right_cell+1) ! ending    index for sol-pts in right cell
    !
    wlocl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if left_cell is a tria or quad
    wlocr = merge( 1 , 2 , r2-r1+1 == n_sp2dt ) ! check if right_cell is a tria or quad
    !
    ! Interpolate the conserved variables and face normals at
    ! the solution flux points to the quadrature flux points
    !
    if (using_quadrature) then
      !
      do k = 1,nqp
        !
        do m = 1,nq
          luqp(m,k) = dot( interp_f2g(1:nsp,k) , faceusp(m,1:nsp,1,nf) )
          ruqp(m,k) = dot( interp_f2g(1:nsp,k) , faceusp(m,1:nsp,2,nf) )
        end do
        !
        do l = 1,nr
          rnrm(l,k) = dot( interp_f2g(1:nsp,k) , geofa(l,1:nsp,nf) )
        end do
        !
      end do
      !
    else
      !
      lflx(1:nq,1:nsp) = faceflx(1:nq,1:nsp,1,nf)
      rflx(1:nq,1:nsp) = faceflx(1:nq,1:nsp,2,nf)
      luqp(1:nq,1:nsp) = faceusp(1:nq,1:nsp,1,nf)
      ruqp(1:nq,1:nsp) = faceusp(1:nq,1:nsp,2,nf)
      rnrm(1:nr,1:nsp) = geofa(1:nr,1:nsp,nf)
      !
    end if
    !
    if (output_boundary_fluxes) then
      if (nf <= nfbnd) then
        write (102,111) nf,trim(bc_integer_to_string(bface(1,nf)))
        do k = 1,nsp
          write (102,112) k,luqp(1:nq,k),pressure_cv_sp(luqp(1:nq,k))
          write (102,113) k,ruqp(1:nq,k),pressure_cv_sp(ruqp(1:nq,k))
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
        lflx(1:nq,k) = normal_invs_flx_cv( luqp(:,k) , rnrm(:,k) )
        rflx(1:nq,k) = normal_invs_flx_cv( ruqp(:,k) , rnrm(:,k) )
      end do
      !
    else
      !
      ! Get the discontinuous fluxes at the solution flux points
      !
      lflx(1:nq,1:nsp) = faceflx(1:nq,1:nsp,1,nf)
      rflx(1:nq,1:nsp) = faceflx(1:nq,1:nsp,2,nf)
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
        lflx(1:nq,1:nqp) = fncq(1:nq,1:nqp) - lflx(1:nq,1:nqp)
        rflx(1:nq,1:nqp) = fncq(1:nq,1:nqp) - rflx(1:nq,1:nqp)
        !
        do k = 1,nsp
          do m = 1,nq
            fjl(m,k) = dot( projct_g2f(1:nqp,k) , lflx(m,1:nqp) )
            fjr(m,k) = dot( projct_g2f(1:nqp,k) , rflx(m,1:nqp) )
          end do
        end do
        !
      else
        !
        do k = 1,nsp
          !
          do m = 1,nq
            fnc(m) = dot( projct_g2f(1:nqp,k) , fncq(m,1:nqp) )
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
      fjl(1:nq,1:nsp) = fncq(1:nq,1:nsp) - lflx(1:nq,1:nsp)
      fjr(1:nq,1:nsp) = fncq(1:nq,1:nsp) - rflx(1:nq,1:nsp)
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
   !write (iu,1) nf,left_cell,right_cell
   !write (iu+1,4) nf,left_cell,right_cell
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
     !do m = 1,nq
     !  write (iu,2) lflx(m,k),rflx(m,k),fjl(m,k),fjr(m,k),fnc(m,k)
     !  write (iu+1,2) luqp(m,k),ruqp(m,k)
     !end do
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( rnrm(:,k) )
      dsl = ds * cl(left_geom)
      dsr = ds * cl(right_geom)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! left_cell
      !
      do n = l1,l2
        do m = 1,nq
          residual(m,n) = residual(m,n) + dsl*fjl(m,k)*cormat(n-l1+1,kl,wlocl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! right_cell
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,nq
            residual(m,n) = residual(m,n) - dsr*fjr(m,k)*cormat(n-r1+1,kr,wlocr)
          end do
        end do
      end if
      !
      ! From ZJ, to ensure conservation with chain-rule approach
      !
     !if (flux_divergence == ChainRuleCV) then
     !  do m = 1,nq   ! ARE DSL AND DSR CORRECT HERE?
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
   !    do m = 1,nq
   !      residual(m,n) = residual(m,n) + srcl(m)/volj(n)
   !    end do
   !  end do
   !  !
   !  if (nf > nfbnd) then
   !    do n = r1,r2
   !      do m = 1,nq
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
  use order_mod,  only : n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,n_qp2dq,n_qp2dt
  use geovar,  only : nr,ncell,pts_in_cell_ptr
  use geovar,  only : flxpts_of_cell,flxpts_of_cell_ptr
  use geovar,  only : interp_s2f,using_quadrature
  use metvar,  only : metrics
  use metvar,  only : dr_mat,ds_mat
  use metvar,  only : metrics,geofa,jacobian
  use flowvar, only : usp,dusp,residual,faceflx,faceusp
  use ovar,    only : dflux_opt
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,wlocl,n1,n2,np
  integer :: f1,f2,nf,sf,kf,lf
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:nq,1:nr)   :: vflx
  real(wp), dimension(1:n_sp2dq) :: lmat
  !
  real(wp), dimension(1:n_sp2dq,1:nq,1:nr) :: rsflx
  real(wp), dimension(1:n_sp2dq,1:nq,1:nr) :: xyflx
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
    np = n2-n1+1 ! number of solution points in this cell
    !
    wlocl = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    f1 = flxpts_of_cell_ptr(n)+1 ! beginning index of the flx-pts for cell n
    f2 = flxpts_of_cell_ptr(n+1) ! ending    index of the flx-pts for cell n
    !
    if (.not. first_output_done .and. output_boundary_fluxes) then
      write (120,110) n
      do k = n1,n2
        write (120,111) k,((metrics(l,m,k),m=1,2),l=1,2),jacobian(k) 
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
      do m = 1,nr
        xyflx(k,:,m) = invs_flx_cv( usp(:,l) , m ) ! x-directional flux
      end do
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
      rsflx(k,:,1) =  xyflx(k,:,1)*metrics(2,2,l) - xyflx(k,:,2)*metrics(1,2,l)
      rsflx(k,:,2) = -xyflx(k,:,1)*metrics(2,1,l) + xyflx(k,:,2)*metrics(1,1,l)
      !
    end do
    !
    ! Compute the divergence of the discontinuous fluxes 
    !
    do m = 1,nq
      do k = 1,np
        residual(m,k+n1-1) = dot( dr_mat(1:np,k,wlocl) , rsflx(1:np,m,1) ) + &
                             dot( ds_mat(1:np,k,wlocl) , rsflx(1:np,m,2) )
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
      lmat(1:np) = interp_s2f(1:np,lf,wlocl)
      !
      do l = 1,nr
        do m = 1,nq
          faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
                                dot(lmat(1:np),xyflx(1:np,m,l))*unrm(l)
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
  use order_mod,  only : n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,n_qp1de
  use geovar,  only : nr,nfbnd,nface
  use geovar,  only : bface,cells_on_face,pts_in_cell_ptr,xyz,flxwts_1d
  use geovar,  only : using_quadrature
  use geovar,  only : cell_geom
  use ovar,    only : flux_divergence,dflux_opt,itcur,nst
  use ovar,    only : invs_flux_method
  use metvar,  only : cormat,geofa,volj
  use flowvar, only : usp,residual,faceflx,faceusp,facedusp
  !
 !use parallel_mod, only : wait_facedusp
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nf,r1,r2,l1,l2
  integer  :: wlocl,wlocr,nsp
  integer  :: kbl,kbr,kl,kr,iu
  integer  :: host_cell,left_cell,right_cell
  integer  :: host_geom,left_geom,right_geom
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)    :: unrm
  real(wp), dimension(1:nq) :: ul,ur
  real(wp), dimension(1:nq) :: dfl,dfr
  real(wp), dimension(1:nq) :: dgl,dgr
  real(wp), dimension(1:nq) :: fnl,fnr,fnc
  real(wp), dimension(1:nq) :: srcl,srcr
  real(wp), dimension(1:nq) :: iflx,vflx
  real(wp), dimension(1:nr,1:n_sp1de)    :: rnrm
  real(wp), dimension(1:nq,1:n_sp1de) :: fjl,fjr
  real(wp), dimension(1:nq,1:n_sp1de) :: lflx,rflx
  real(wp), dimension(1:nq,1:n_sp1de) :: luqp,ruqp,fncq
  real(wp), dimension(1:nr,1:nq,1:n_sp1de) :: lduqp,rduqp
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
    host_cell = cells_on_face(1,nf) ! left  cell on face
    !
    host_geom = cell_geom(host_cell)
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    !
    l1 = pts_in_cell_ptr(host_cell)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(host_cell+1) ! ending    index for sol-pts in left cell
    !
    wlocl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if host_cell is a tria or quad
    !
    lflx(1:nq,1:nsp) = faceflx(1:nq,1:nsp,1,nf)
    !
    luqp(1:nq,1:nsp) = faceusp(1:nq,1:nsp,1,nf)
    ruqp(1:nq,1:nsp) = faceusp(1:nq,1:nsp,2,nf)
    !
    if (output_boundary_fluxes) then
      if (nf <= nfbnd) then
        write (101,111) nf,trim(bc_integer_to_string(bface(1,nf)))
        do k = 1,nsp
          write (101,112) k,luqp(1:nq,k),pressure_cv_sp(luqp(1:nq,k))
          write (101,113) k,ruqp(1:nq,k),pressure_cv_sp(ruqp(1:nq,k))
        end do
      end if
    end if
    !
    lduqp(1:nr,1:nq,1:nsp) = facedusp(1:nr,1:nq,1:nsp,1,nf)
    rduqp(1:nr,1:nq,1:nsp) = facedusp(1:nr,1:nq,1:nsp,2,nf)
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
        flux_cell = host_cell
        cell_flux_point = k
        vflx(:) = common_visc_flx_cv( unrm(:), luqp(:,k), lduqp(:,:,k), &
                                               ruqp(:,k), rduqp(:,:,k) )
        !
      end if
      if (output_boundary_fluxes) then
        write (103,201) nf,k,trim(bc_integer_to_string(bface(1,nf)))
        write (103,203) iflx(1:nq)
        write (103,204) vflx(1:nq)
        write (103,205) lflx(1:nq,k)
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
    fjl(1:nq,1:nsp) = fncq(1:nq,1:nsp) - lflx(1:nq,1:nsp)
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
      dsl = ds * cl(host_geom)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! host_cell
      !
      do n = l1,l2
        do m = 1,nq
          residual(m,n) = residual(m,n) + dsl*fjl(m,k)*cormat(n-l1+1,kl,wlocl)
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
    left_cell = cells_on_face(1,nf) ! left  cell on face
    right_cell = cells_on_face(2,nf) ! right cell on face
    !
    left_geom = cell_geom(left_cell)
    right_geom = cell_geom(right_cell)
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! ending    flux point for right cell
    !
    l1 = pts_in_cell_ptr(left_cell)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(left_cell+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(right_cell)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(right_cell+1) ! ending    index for sol-pts in right cell
    !
    wlocl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if left_cell is a tria or quad
    wlocr = merge( 1 , 2 , r2-r1+1 == n_sp2dt ) ! check if right_cell is a tria or quad
    !
    lflx(1:nq,1:nsp) = faceflx(1:nq,1:nsp,1,nf)
    rflx(1:nq,1:nsp) = faceflx(1:nq,1:nsp,2,nf)
    !
    luqp(1:nq,1:nsp) = faceusp(1:nq,1:nsp,1,nf)
    ruqp(1:nq,1:nsp) = faceusp(1:nq,1:nsp,2,nf)
    !
    lduqp(1:nr,1:nq,1:nsp) = facedusp(1:nr,1:nq,1:nsp,1,nf)
    rduqp(1:nr,1:nq,1:nsp) = facedusp(1:nr,1:nq,1:nsp,2,nf)
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
      flux_cell = left_cell
      cell_flux_point = right_cell
      vflx(:) = common_visc_flx_cv( unrm(:), luqp(:,k), lduqp(:,:,k), &
                                             ruqp(:,k), rduqp(:,:,k) )
      if (output_boundary_fluxes) then
        write (104,202) nf,k
        write (104,203) iflx(1:nq)
        write (104,204) vflx(1:nq)
        write (104,205) lflx(1:nq,k)
        write (104,206) rflx(1:nq,k)
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
    fjl(1:nq,1:nsp) = fncq(1:nq,1:nsp) - lflx(1:nq,1:nsp)
    fjr(1:nq,1:nsp) = fncq(1:nq,1:nsp) - rflx(1:nq,1:nsp)
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
      dsl = ds * cl(left_geom)
      dsr = ds * cl(right_geom)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! left_cell
      !
      do n = l1,l2
        do m = 1,nq
          residual(m,n) = residual(m,n) + dsl*fjl(m,k)*cormat(n-l1+1,kl,wlocl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! right_cell
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      do n = r1,r2
        do m = 1,nq
          residual(m,n) = residual(m,n) - dsr*fjr(m,k)*cormat(n-r1+1,kr,wlocr)
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
  use geovar, only : nr
  use ovar,   only : cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
  real(wp), dimension(1:nr,1:nq), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq,1:nr) :: return_value
  !
  !.. Local Scalars ..
  integer  :: vb,ve
  real(wp) :: dudx,dudy,dvdx,dvdy
  real(wp) :: mvsc,mcdt,dTdx,dTdy
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vel
  real(wp), dimension(1:nr,1:nq) :: dpvdx
  !
continue
  !
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
  return_value(nq,1) = mcdt*dTdx + dot( vel , return_value(vb:ve,1) )
  !
  ! y-directional viscous fluxes
  !
  return_value( 1,2) = zero
  return_value( 2,2) = mvsc * (dudy + dvdx)
  return_value( 3,2) = mvsc * two3 * (two*dvdy - dudx)
  return_value(nq,2) = mcdt*dTdy + dot( vel , return_value(vb:ve,2) )
  !
end function visc_flx_cv
!
!###############################################################################
!
pure function normal_visc_flx_cv(rnrm,uu,dcvdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
  real(wp), dimension(1:nr),         intent(in) :: rnrm
  real(wp), dimension(1:nr,1:nq), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: l
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)         :: unrm
  real(wp), dimension(1:nq,1:nr) :: vflx
  !
continue
  !
  unrm(:) = unit_vector( rnrm(:) )
  !
  vflx(1:nq,1:nr) = visc_flx_cv(uu,dcvdx)
  !
  return_value(1:nq) = unrm(1)*vflx(1:nq,1)
  do l = 2,nr
    return_value(1:nq) = return_value(1:nq) + unrm(l)*vflx(1:nq,l)
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
  use geovar, only : nr
  use ovar,   only : cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),         intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: uwall
  real(wp), dimension(1:nr,1:nq), intent(in) :: dpvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: vb,ve
  real(wp) :: dudx,dudy,dvdx,dvdy
  real(wp) :: mvsc,mcdt,dTdx,dTdy
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)    :: vel
  real(wp), dimension(1:nq) :: xflx,yflx
  !
continue
  !
  ! For slip wall: all gradients /= zero
  !                uwall /= zero
  !
  ! For adiabatic wall: Temperature gradient == zero
  !                     Velocity gradients /= zero
  !                     uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! For isothermal wall: all gradients /= zero
  !                      uwall(nmb:nme) [i.e., wall velocities] == zero
  !                
  !
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
  xflx(nq) = mcdt*dTdx + dot( vel , xflx(vb:ve) )
  !
  ! y-directional viscous fluxes
  !
  yflx( 1) = zero
  yflx( 2) = mvsc * (dudy + dvdx)
  yflx( 3) = mvsc * two3 * (two*dvdy - dudx)
  yflx(nq) = mcdt*dTdy + dot( vel , yflx(vb:ve) )
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
  use geovar, only : nr
  use ovar,   only : visc_variables_to_ave
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),         intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: ul,ur
  real(wp), dimension(1:nr,1:nq), intent(in) :: duldx,durdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: vjump,qjump,penalty
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq)      :: ucom
  real(wp), dimension(1:nr,1:nq) :: qcom
  real(wp), dimension(1:nq,1:nr) :: lvflx
  real(wp), dimension(1:nq,1:nr) :: rvflx
  real(wp), dimension(1:nr,1:nq) :: vflux
  real(wp), dimension(1:nq,1:nr) :: fvcom
  !
continue
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    ! Get the left and right viscous fluxes
    !
    lvflx(1:nq,1:nr) = visc_flx_cv(ul,duldx)
    rvflx(1:nq,1:nr) = visc_flx_cv(ur,durdx)
    !
    ! Use the LDG method to compute the common viscous flux vector
    !
    ! Start by applying the averaging operator {{.}} to the viscous flux vectors
    !
    do m = 1,nq
      do l = 1,nr
        vflux(l,m) = half * (lvflx(m,l) + rvflx(m,l))
      end do
    end do
    !
    ! Apply the differencing (jump) operator [[.]] to the interface viscous
    ! flux vectors and then add it to the common viscous flux vector
    !
    do m = 1,nq
      !
      vjump = dot( unrm(1:nr) , lvflx(m,1:nr)-rvflx(m,1:nr) )
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
    do m = 1,nq
      penalty = LDG_tau * ( ul(m) - ur(m) )
      do l = 1,nr
        vflux(l,m) = vflux(l,m) + penalty * unrm(l)
      end do
    end do
    !
    ! Finally, dot the common viscous flux vector with the face unit normal
    !
    do m = 1,nq
      return_value(m) = dot( vflux(1:nr,m) , unrm(1:nr) )
    end do
    !
  else if (visc_variables_to_ave == visc_solution_and_gradients) then
    !
    ! Get the common solution at the interface using a simple average
    !
    do m = 1,nq
      ucom(m) = half * ( ul(m) + ur(m) )
    end do
    !
    ! Use the LDG method to compute the common gradient
    !
    ! Start by applying the averaging operator {{.}} to the gradients
    !
    do m = 1,nq
      do l = 1,nr
        qcom(l,m) = half * (duldx(l,m) + durdx(l,m))
      end do
    end do
    !
    ! Apply the differencing (jump) operator [[.]] to the
    ! gradients and then add it to the common gradient
    !
    do m = 1,nq
      !
      qjump = dot( unrm(1:nr) , duldx(1:nr,m)-durdx(1:nr,m) )
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
    do m = 1,nq
      penalty = LDG_tau * ( ul(m) - ur(m) )
      do l = 1,nr
        qcom(l,m) = qcom(l,m) + penalty * unrm(l)
      end do
    end do
    !
    ! Compute the common viscous flux vector using
    ! the common solution and common gradient
    !
    fvcom(1:nq,1:nr) = visc_flx_cv(ucom,qcom)
    !
    ! Finally, dot the common viscous flux vector with the face unit normal
    !
    do m = 1,nq
      return_value(m) = dot( fvcom(m,1:nr) , unrm(1:nr) )
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
  use order_mod,  only : n_sp2dt,n_sp2dq
  use geovar,  only : ncell
  use geovar,  only : pts_in_cell_ptr
  use metvar,  only : jacobian
  use metvar,  only : metrics
  use metvar,  only : dr_mat,ds_mat
  use flowvar, only : usp,dusp
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,n1,n2,np,wlocl
  real(wp) :: drdx,drdy,dsdx,dsdy
  real(wp) :: dudr,duds,dudx,dudy
  real(wp) :: jac_inv
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq,1:nq) :: tusp
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
    wlocl = merge(1,2,np == n_sp2dt) ! check if cell n is tria or quad
    !
    ! Store the transpose of usp for this cell to reduce cache misses
    !
    do k = 1,np
      do m = 1,nq
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
      drdx =  metrics(2,2,l)
      drdy = -metrics(1,2,l)
      dsdx = -metrics(2,1,l)
      dsdy =  metrics(1,1,l)
     !drdx =  jac_inv * metrics(2,2,l)
     !drdy = -jac_inv * metrics(1,2,l)
     !dsdx = -jac_inv * metrics(2,1,l)
     !dsdy =  jac_inv * metrics(1,1,l)
      !
      do m = 1,nq
        !
        ! Compute the derivatives in r-s space
        !
        dudr = dot( dr_mat(1:np,k,wlocl) , tusp(1:np,m) )
        duds = dot( ds_mat(1:np,k,wlocl) , tusp(1:np,m) )
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
  use order_mod,  only : n_sp1de,n_sp2dt
  use geovar,  only : nfbnd,nface,nr,n_solpts,ncell
  use geovar,  only : cells_on_face,pts_in_cell_ptr,bface
  use geovar,  only : cell_geom
  use metvar,  only : cormat,geofa,jacobian,volj
  use flowvar, only : dusp,faceusp
  use ovar,    only : walls_are_exact
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nf,r1,r2,l1,l2,n1,n2,nc
  integer  :: wlocl,wlocr,nsp
  integer  :: kbl,kbr,kl,kr
  integer  :: host_cell,left_cell,right_cell
  integer  :: host_geom,left_geom,right_geom
  real(wp) :: cbias,hfpb,hfmb
  real(wp) :: correct_n,correct_m,ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)              :: unrm
  real(wp), dimension(1:nq)           :: ljfp,rjfp
  real(wp), dimension(1:nq,1:n_sp1de) :: lufp,rufp
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
    host_cell = cells_on_face(1,nf) ! left  cell on face
    !
    host_geom = cell_geom(host_cell)
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    !
    l1 = pts_in_cell_ptr(host_cell)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(host_cell+1) ! ending    index for sol-pts in left cell
    !
    wlocl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if host_cell is a tria or quad
    !
    ! For all the flux points on the current face, extract the conservative
    ! conservative variables for the left and right cells from faceusp
    !
    lufp(1:nq,1:nsp) = faceusp(1:nq,1:nsp,1,nf)
    rufp(1:nq,1:nsp) = faceusp(1:nq,1:nsp,2,nf)
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! Use the LDG flux formulation for computing the common interface solutions
    !   u_com = half*(u_left + u_right) - &
    !                (u_left - u_right)*dot(LDG_beta,unit_normal)
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
      dsl = cl(host_geom) * rnorm2( geofa(1:nr,k,nf) )
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
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      do n = l1,l2
        correct_n = dsl * cormat(n-l1+1,kl,wlocl) / jacobian(n)
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
    left_cell = cells_on_face(1,nf) ! left  cell on face
    right_cell = cells_on_face(2,nf) ! right cell on face
    !
    left_geom = cell_geom(left_cell)
    right_geom = cell_geom(right_cell)
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! ending    flux point for right cell
    !
    l1 = pts_in_cell_ptr(left_cell)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(left_cell+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(right_cell)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(right_cell+1) ! ending    index for sol-pts in right cell
    !
    wlocl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if left_cell is a tria or quad
    wlocr = merge( 1 , 2 , r2-r1+1 == n_sp2dt ) ! check if right_cell is a tria or quad
    !
    ! For all the flux points on the current face, extract the conservative
    ! conservative variables for the left and right cells from faceusp
    !
    lufp(1:nq,1:nsp) = faceusp(1:nq,1:nsp,1,nf)
    rufp(1:nq,1:nsp) = faceusp(1:nq,1:nsp,2,nf)
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! Use the LDG flux formulation for computing the common interface solutions
    !   u_com = half*(u_left + u_right) - &
    !                (u_left - u_right)*dot(LDG_beta,unit_normal)
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
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      do n = l1,l2
        correct_n = dsl * cormat(n-l1+1,kl,wlocl) / jacobian(n)
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
      do n = r1,r2
        correct_n = dsr * cormat(n-r1+1,kr,wlocr) / jacobian(n)
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
  use order_mod,  only : n_sp2dt,n_sp2dq,n_sp1de,n_order
  use geovar,  only : flxpts_of_cell_ptr,flxpts_of_cell
  use geovar,  only : nr,ncell,nfbnd,bface
  use geovar,  only : interp_s2f,pts_in_cell_ptr
  use geovar,  only : cells_on_face
  use geovar,  only : cell_geom
  use metvar,  only : geofa,cormat,jacobian,volj
  use flowvar, only : dusp,facedusp,faceusp
  use ovar,    only : walls_are_exact
  !
  use parallel_mod, only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,wlocl,k,l,m
  integer :: nf,f1,f2,kf,lf,sf
  integer :: ncf,fp1,fp2,kfp,kbl,kl
  integer :: this_geom
  real(wp) :: dsl,cbias,hfpb,correct_n,correct_m
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:nq) :: lufp,rufp,ujump
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_sp2dq,1:nr,1:nq) :: face_dusp
  real(wp), dimension(1:nr,1:nq,1:n_sp2dq) :: cell_dusp
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
    wlocl = merge(1,2,np == n_sp2dt)
    !
    this_geom = cell_geom(nc)
    !
    f1 = flxpts_of_cell_ptr(nc)+1
    f2 = flxpts_of_cell_ptr(nc+1)
    !
    ! Create transpose of dusp array to reduce cache misses
    !
    do n = 1,np
      do m = 1,nq
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
        do m = 1,nq
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
        dsl = cl(this_geom) * rnorm2( geofa(1:nr,kf,nf) )
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
        cbias = dot( LDG_beta(1:nr) , unrm(1:nr) )
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
        do m = 1,nq
          ujump(m) = hfpb*(rufp(m) - lufp(m))
        end do
        !
        ! Apply the jump correction to the interior discontinuous gradients
        !
        do n = 1,np
          correct_n = dsl * cormat(n,lf,wlocl) / jacobian(n+n1-1)
          do m = 1,nq
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
        lmat(1:np) = interp_s2f(1:np,lf,wlocl)
        !
        do m = 1,nq
          do l = 1,nr
            facedusp(l,m,kf,sf,nf) = dot(lmat(1:np),face_dusp(1:np,l,m))
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
      do m = 1,nq
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
!
include "Functions/upwind_fluxes.f90"
include "Functions/grad_cv_to_grad_pv.f90"
include "Functions/vsp2u.f90"
include "Functions/usp2v.f90"
include "Functions/pressure.f90"
include "Functions/viscosity.f90"
!
!###############################################################################
!
end module flux_mod
