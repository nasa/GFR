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
  public :: flux_gradient_at_solpts_with_qdtr
  public :: flux_gradient_at_solpts_visc
  public :: correct_interface_flux_with_qdtr
  public :: correct_interface_flux_visc
  public :: discontinuous_gradient_at_solpts
  public :: correct_gradients_at_solpts
  public :: correct_gradients_br2
  public :: correct_solution_for_output
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
 !real(wp), parameter :: LDG_beta(3) = (/ zero, zero, zero /)
 !real(wp), parameter :: LDG_beta(3) = (/ half, half, half /)
  !
  ! LDG_TAU : Constant that controls the penalty term
  !           term of the common diffusive flux
  !        = 
 !real(wp), parameter :: LDG_tau = zero
 !real(wp), parameter :: LDG_tau = one
  !
contains
!
!###############################################################################
!
subroutine interp_cv_to_flxpts()
  !
  use order_mod,  only : maxSP
  use geovar,  only : ncell,cell
  use interpolation_mod, only : interp
  use flowvar, only : usp,faceusp
  !
  use parallel_mod, only : exchange_faceusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,m
  integer :: nf,kf,lf,sf
  integer :: this_geom,this_order,face_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxSP,1:nq) :: tusp
#ifndef USE_MATRIX_PROCEDURES
  real(wp), dimension(1:maxSP) :: lmat
#endif
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
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
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      do kf = 1,size(cell(nc)%face(n)%cell_idx) ! index of the flx-pt
        !                                         local to this face
        lf = cell(nc)%face(n)%cell_idx(kf) ! index of flx-pt within interp%mat
        !
        !######################################################################
        !######################################################################
        ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
        ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
        ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
        !######################################################################
        !######################################################################
        !
#ifdef USE_MATRIX_PROCEDURES
        do m = 1,nq
          faceusp(nf)%v(m,kf,sf) = interp(this_geom,this_order)% &
                                             toFace(face_order)% &
                                   dot( lf , tusp(1:np,m) )
        end do
#else
        lmat(1:np) = interp(this_geom,this_order)% &
                               toFace(face_order)% &
                     mat(1:np,lf)
        !
        do m = 1,nq
          faceusp(nf)%v(m,kf,sf) = dot( lmat(1:np) , tusp(1:np,m) )
        end do
#endif
        !
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
  use order_mod,  only : maxSP
  use geovar,  only : nr,ncell,cell
  use interpolation_mod, only : interp
  use flowvar, only : dusp,facedusp
  !
  use parallel_mod, only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,l,m
  integer :: nf,kf,lf,sf
  integer :: this_geom,this_order,face_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tdusp
#ifndef USE_MATRIX_PROCEDURES
  real(wp), dimension(1:maxSP) :: lmat
#endif
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
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
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      do kf = 1,size(cell(nc)%face(n)%cell_idx) ! index of the flx-pt
        !                                         local to this face
        lf = cell(nc)%face(n)%cell_idx(kf) ! index of flx-pt within interp%mat
        !
        !######################################################################
        !######################################################################
        ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
        ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
        ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
        !######################################################################
        !######################################################################
        !
#ifdef USE_MATRIX_PROCEDURES
        do m = 1,nq
          do l = 1,nr
            facedusp(nf)%v(l,m,kf,sf) = interp(this_geom,this_order)% &
                                                  toFace(face_order)% &
                                        dot( lf , tdusp(1:np,l,m) )
          end do
        end do
#else
        lmat(1:np) = interp(this_geom,this_order)% &
                               toFace(face_order)% &
                     mat(1:np,lf)
        !
        do m = 1,nq
          do l = 1,nr
            facedusp(nf)%v(l,m,kf,sf) = dot( lmat(1:np) , tdusp(1:np,l,m) )
          end do
        end do
#endif
        !
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
! Procedures for computing the divergence of inviscid and viscous fluxes
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine flux_gradient_at_solpts_visc
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,  only : nr,ncell,using_quadrature
  use geovar,  only : cell
  use interpolation_mod, only : interp
  use derivatives_mod, only : deriv
  use metrics_mod,  only : metrics,geofa,jacobian
  use flowvar, only : usp,dusp,residual,faceflx,faceusp
  use ovar,    only : governing_equations
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nc,n1,n2,np
  integer :: nf,sf,kf,lf
  integer :: this_geom,this_order,face_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:nq,1:nr) :: vflx
  real(wp), dimension(1:nq,1:maxSP) :: lcl_res
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: rsflx
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: xyflx
#ifndef USE_MATRIX_PROCEDURES
  real(wp), dimension(1:maxSP)   :: lmat
  real(wp), dimension(1:maxSP,1:maxSP)   :: dmat
#endif
  !
continue
  !
  ! Initialize faceflx to zero before we start
  !
  call faceflx%zero
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! #########
    ! #########
    ! ### 1 ###
    ! #########
    ! #########
    !
    ! Get the flux at each quadrature point in the current cell
    !
    do k = 1,np
      !
      n = k + n1-1
      !
      ! Inviscid fluxes
      !
      do m = 1,nr
        xyflx(k,:,m) = invs_flx_cv( usp(:,n) , m ) ! x-directional flux
      end do
      !
      if (governing_equations == NavierStokes_Eqns) then
        !
        ! Viscous fluxes
        !
        vflx(1:nq,1:nr) = visc_flx_cv( usp(1:nq,n) , dusp(1:nr,1:nq,n) )
        !
        ! Subtract the viscous flux from the inviscid flux
        !
        xyflx(k,1:nq,1:nr) = xyflx(k,1:nq,1:nr) - vflx(1:nq,1:nr)
        !
      end if
      !
      ! Transform the fluxes in the x-y physical
      ! domain to the r-s computational domain
      !
      do m = 1,nq
        do l = 1,nr
          rsflx(k,m,l) = dot( xyflx(k,m,1:nr) , metrics(1:nr,l,n,2) )
        end do
      end do
      !
    end do
    !
    ! #########
    ! #########
    ! ### 2 ###
    ! #########
    ! #########
    !
    ! Compute the divergence of the discontinuous fluxes 
    !
    ! Zero the residual array for the current cell
    !
    lcl_res = zero
    !
    ! Loop through each coordinate direction and add the divergence of
    ! the flux to the local residual array for the current cell
    !
#ifdef USE_MATRIX_PROCEDURES
    do k = 1,np
      do l = 1,nr
        do m = 1,nq
          vflx(m,l) = deriv(this_geom,this_order)% &
                                             d(l)% &
                      dot( k , rsflx(1:np,m,l) )
        end do
      end do
      lcl_res(1:nq,k) = sum( vflx , dim=2 )
    end do
#else
    do l = 1,nr
      !
      dmat(1:np,1:np) = deriv(this_geom,this_order)%d(l)%mat(1:np,1:np)
      !
      do k = 1,np
        do m = 1,nq
          lcl_res(m,k) = lcl_res(m,k) + dot( dmat(1:np,k) , rsflx(1:np,m,l) )
        end do
      end do
      !
    end do
#endif
    !
    ! Copy the local residual for the current cell into the residual array
    !
    residual(1:nq,n1:n2) = lcl_res(1:nq,1:np)
    !
    !
    ! #########
    ! #########
    ! ### 3 ###
    ! #########
    ! #########
    !
    ! Interpolate the fluxes to the interface flux points and dot this
    ! interpolated vector with the face normal at each flux point
    ! NOTE: We dont need to exchange this information for communication and
    !       periodic boundary faces since this flux is only used within the
    !       current cell.
    !
    ! Loop over the faces of this cell
    !
    do n = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(n)%idx  ! index of the current face
      sf = cell(nc)%face(n)%side ! side of face that cell nc is on
      !
      face_order = cell(nc)%face(n)%order ! order of current face
      !
     !nf = size(interp(this_geom,this_order)%toFace(face_order)%mat,dim=2)
     !write (161,15) ((interp(this_geom,this_order)% &
     !                           toFace(face_order)% &
     !                 mat(k,l), k=1,np), l=1,nf)
     !15 format (25(5es16.8,/),/)
     !call stop_recon(stop_mpi,"flux_gradient_at_solpts_visc")
      !
      do kf = 1,size(cell(nc)%face(n)%cell_idx) ! index of the flx-pt
        !                                         local to this face
        lf = cell(nc)%face(n)%cell_idx(kf) ! index of flx-pt within interp%mat
        !
        ! Get the unit normal for this face flux/solution point
        !
        unrm(:) = unit_vector( geofa(nf)%v(:,kf) )
        !
        !######################################################################
        !######################################################################
        ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
        ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
        ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
        !######################################################################
        !######################################################################
        !
#ifdef USE_MATRIX_PROCEDURES
        !
        do l = 1,nr
          do m = 1,nq
            vflx(m,l) = interp(this_geom,this_order)% &
                                  toFace(face_order)% &
                        dot( lf , xyflx(1:np,m,l) )
          end do
        end do
        !
        do m = 1,nq
          faceflx(nf)%v(m,kf,sf) = dot( vflx(m,1:nr) , unrm(1:nr) )
        end do
    !   do l = 1,nr
    !     do m = 1,nq
    !       faceflx(nf)%v(m,kf,sf) = faceflx(nf)%v(m,kf,sf) + unrm(l) * &
    !interp(this_geom,this_order)%toFace(face_order)%dot( lf , xyflx(1:np,m,l) )
    !     end do
    !   end do
#else
        lmat(1:np) = interp(this_geom,this_order)% &
                               toFace(face_order)% &
                     mat(1:np,lf)
        !
        do l = 1,nr
          do m = 1,nq
            faceflx(nf)%v(m,kf,sf) = faceflx(nf)%v(m,kf,sf) + &
                            dot( lmat(1:np) , xyflx(1:np,m,l) ) * unrm(l)
          end do
        end do
#endif
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
subroutine correct_interface_flux_visc()
  !
  use generic_types_mod, only : matrix
  use order_mod,  only : maxFP,geom_solpts
  use geovar,  only : nr,nfbnd,nface
  use geovar,  only : bface,face,cell
  use ovar,    only : invs_flux_method,governing_equations
  use metrics_mod,  only : geofa,volj
  use correction_mod, only : correct
  use flowvar, only : usp,residual,faceflx,faceusp,facedusp
  !
  !.. Local Scalars ..
  integer  :: i,k,l,m,n,nf,r1,r2,l1,l2
  integer  :: kl,kr,iu,nfp
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: ds,dsl,dsr,cmat_val
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:nq) :: iflx,vflx
  real(wp), dimension(1:nr,1:maxFP) :: rnrm
  real(wp), dimension(1:nq,1:maxFP) :: fjl,fjr
  real(wp), dimension(1:nq,1:maxFP) :: lflx,rflx
  real(wp), dimension(1:nq,1:maxFP) :: luqp,ruqp,fncq
  real(wp), dimension(1:nr,1:nq,1:maxFP) :: lduqp,rduqp
  !
continue
  !
  vflx = zero
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
    ! value of interior flux polynomial at FP
    lflx(1:nq,1:nfp) = faceflx(nf)%v(1:nq,1:nfp,1)
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
    ! Compute the flux jump at the interface
    !
    fjl(1:nq,1:nfp) = fncq(1:nq,1:nfp) - lflx(1:nq,1:nfp)
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
      kl = face(nf)%left%flxpts(k)
      !
      ! Get the contribution to the face length/area from the current flux
      ! point for the left cell
      !
      ds  = rnorm2( rnrm(:,k) )
      dsl = ds * cl(host_geom)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! host_cell
      !
#ifdef USE_MATRIX_PROCEDURES
      do i = correct(host_geom,host_order)%row_ptr(kl)+1, &
             correct(host_geom,host_order)%row_ptr(kl+1)
        n = l1 - 1 + correct(host_geom,host_order)%row_idx(i)
        cmat_val = dsl * correct(host_geom,host_order)%row_val(i)
        do m = 1,nq
          residual(m,n) = residual(m,n) + fjl(m,k) * cmat_val
        end do
      end do
#else
      do n = l1,l2
        cmat_val = correct(host_geom,host_order)%mat(n-l1+1,kl)
        if (cmat_val /= zero) then
          do m = 1,nq
            residual(m,n) = residual(m,n) + dsl * fjl(m,k) * cmat_val
          end do
        end if
      end do
#endif
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
    lflx(1:nq,1:nfp) = faceflx(nf)%v(1:nq,1:nfp,1)
    rflx(1:nq,1:nfp) = faceflx(nf)%v(1:nq,1:nfp,2)
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
    ! Compute the flux jump at the interface
    !
    fjl(1:nq,1:nfp) = fncq(1:nq,1:nfp) - lflx(1:nq,1:nfp)
    fjr(1:nq,1:nfp) = fncq(1:nq,1:nfp) - rflx(1:nq,1:nfp)
    !
    ! Apply the flux correction due to the flux jump
    !
    do k = 1,nfp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = face(nf)%left%flxpts(k)
      kr = face(nf)%right%flxpts(k)
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
#ifdef USE_MATRIX_PROCEDURES
      do i = correct(left_geom,left_order)%row_ptr(kl)+1, &
             correct(left_geom,left_order)%row_ptr(kl+1)
        n = l1 - 1 + correct(left_geom,left_order)%row_idx(i)
        cmat_val = dsl * correct(left_geom,left_order)%row_val(i)
        do m = 1,nq
          residual(m,n) = residual(m,n) + fjl(m,k) * cmat_val
        end do
      end do
#else
      do n = l1,l2
        cmat_val = correct(left_geom,left_order)%mat(n-l1+1,kl)
        if (cmat_val /= zero) then
          do m = 1,nq
            residual(m,n) = residual(m,n) + dsl * fjl(m,k) * cmat_val
          end do
        end if
      end do
#endif
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! right_cell
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
#ifdef USE_MATRIX_PROCEDURES
      do i = correct(right_geom,right_order)%row_ptr(kr)+1, &
             correct(right_geom,right_order)%row_ptr(kr+1)
        n = r1 - 1 + correct(right_geom,right_order)%row_idx(i)
        cmat_val = dsr * correct(right_geom,right_order)%row_val(i)
        do m = 1,nq
          residual(m,n) = residual(m,n) - fjr(m,k) * cmat_val
        end do
      end do
#else
      do n = r1,r2
        cmat_val = correct(right_geom,right_order)%mat(n-r1+1,kr)
        if (cmat_val /= zero) then
          do m = 1,nq
            residual(m,n) = residual(m,n) - dsr * fjr(m,k) * cmat_val
          end do
        endif
      end do
#endif
      !
    end do
    !
  end do interior_fluxes
  !
  1 format (a,i0)
  !
end subroutine correct_interface_flux_visc
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
  integer :: i,j,l
  real(wp) :: mvsc,mcdt
  real(wp) :: ave_ext_rate
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vel
  real(wp), dimension(1:nr,1:nq) :: dpvdx
  !
continue
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  return_value = zero
  !
  ! Get the gradients of the primitive variables using the 
  ! gradients of the conservative variables and the product rule
  !
  dpvdx(1:nr,1:nq) = grad_cv_to_grad_pv_sp(uu,dcvdx)
  !
  ! Get the velocity vector
  !
  vel(:) = uu(nmb:nme)/uu(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = viscosity_cv_sp(uu)
  mcdt = cpref_nd * mvsc / Pr
  !
  ! Compute the average extension rate
  !
  ave_ext_rate = zero
  do l = 1,nr
    ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
  end do
  ave_ext_rate = two3 * ave_ext_rate
  !
  ! Compute the viscous fluxes for the momentum equations which are just
  ! the components of the viscous stress tensor.
  ! NOTE: We are initially excluding the molecular viscosity term which
  !       will be multiplied to the whole tau array after these loops.
  !
  do j = 1,nr
    do i = 1,nr
        return_value(nmb-1+i,j) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
    end do
  end do
  !
  ! Subtract the average extension rate from the normal stresses
  !
  do l = 1,nr
    return_value(nmb-1+l,l) = return_value(nmb-1+l,l) - ave_ext_rate
  end do
  !
  ! Multiply the viscous stresses by the molecular viscosity
  !
  return_value(nmb:nme,1:nr) = mvsc * return_value(nmb:nme,1:nr)
  !
  ! Compute the viscous fluxes for the energy equation
  !
  do l = 1,nr
    return_value(nq,l) = mcdt*dpvdx(l,nq) + dot( vel , return_value(nmb:nme,l) )
  end do
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
  real(wp), dimension(1:nr),      intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: uwall
  real(wp), dimension(1:nr,1:nq), intent(in) :: dpvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,l,m
  real(wp) :: mvsc,mcdt
  real(wp) :: ave_ext_rate
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: vel
  real(wp), dimension(1:nr,1:nq) :: vflx
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
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  vflx = zero
  !
  ! Get the velocity vector
  !
  vel(:) = uwall(nmb:nme)/uwall(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = viscosity_cv_sp(uwall)
  mcdt = cpref_nd * mvsc / Pr
  !
  ! Compute the average extension rate
  !
  ave_ext_rate = zero
  do l = 1,nr
    ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
  end do
  ave_ext_rate = two3 * ave_ext_rate
  !
  ! Compute the viscous fluxes for the momentum equations which are just
  ! the components of the viscous stress tensor.
  ! NOTE: We are initially excluding the molecular viscosity term which
  !       will be multiplied to the whole tau array after these loops.
  !
  do j = 1,nr
    do i = 1,nr
        vflx(i,nmb-1+j) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
    end do
  end do
  !
  ! Subtract the average extension rate from the normal stresses
  !
  do l = 1,nr
    vflx(l,nmb-1+l) = vflx(l,nmb-1+l) - ave_ext_rate
  end do
  !
  ! Multiply the viscous stresses by the molecular viscosity
  !
  vflx(1:nr,nmb:nme) = mvsc * vflx(1:nr,nmb:nme)
  !
  ! Compute the viscous fluxes for the energy equation
  !
  do l = 1,nr
    vflx(l,nq) = mcdt*dpvdx(l,nq) + dot( vel , vflx(l,nmb:nme) )
  end do
  !
  ! Compute the normal viscous flux
  !
  do m = 1,nq
    return_value(m) = dot( vflx(1:nr,m) , unrm(1:nr) )
  end do
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
  use ovar,   only : LDG_tau, LDG_beta
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
    ! Finally, add the penalty term, (simply the differencing (jump)
    ! operator [[.]] applied to the interface solutions) to the
    ! common viscous flux vector
    !
    do m = 1,nq
      penalty = LDG_tau * ( ur(m) - ul(m) )
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
      penalty = LDG_tau * ( ur(m) - ul(m) )
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
  use order_mod,  only : maxpts
  use geovar,  only : nr,ncell,cell
  use metrics_mod,  only : jacobian
  use metrics_mod,  only : metrics
  use derivatives_mod, only : deriv
  use flowvar, only : usp,dusp
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nc,n1,n2,np
  integer  :: this_geom,this_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: dudr
  real(wp), dimension(1:nr) :: dudx
  real(wp), dimension(1:maxpts,1:nq) :: tusp
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
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
      n = k+n1-1 ! total index of solution point on current partition
      !
      do m = 1,nq
        !
        ! Compute the derivatives in r-s space
        !
        do l = 1,nr
#ifdef USE_MATRIX_PROCEDURES
          dudr(l) = deriv(this_geom,this_order)% &
                                           d(l)% &
                    dot( k , tusp(1:np,m) )
#else
          dudr(l) = dot( deriv(this_geom,this_order)% &
                                                d(l)% &
                         mat(1:np,k) , tusp(1:np,m) )
#endif
        end do
        !
        ! Transform the derivatives back to x-y space
        !
        do l = 1,nr
          dusp(l,m,n) = dot( dudr(1:nr) , metrics(1:nr,l,n,2) ) / jacobian(n)
        end do
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
  use generic_types_mod, only : matrix
  use order_mod,  only : maxpts,geom_solpts
  use geovar,  only : nfbnd,nface,nr,n_solpts,ncell
  use geovar,  only : face,cell,bface
  use metrics_mod,  only : geofa,jacobian,volj
  use correction_mod, only : correct
  use flowvar, only : dusp,faceusp
  use ovar,    only : walls_are_exact
  use ovar,   only : LDG_beta
  !
  !.. Local Scalars ..
  integer  :: i,k,l,m,n,nf,r1,r2,l1,l2,n1,n2,nc
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
      kl = face(nf)%left%flxpts(k)
      !
      ! Get the unit normal at this flux point
      !
      unrm(1:nr) = unit_vector( geofa(nf)%v(1:nr,k) )
      !
      ! Get the face length for the current face
      !
      dsl = cl(host_geom) * rnorm2( geofa(nf)%v(1:nr,k) )
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
      if (any(bface(1,nf) == bc_walls)) then
        ucom(1:nq) = rufp(1:nq,k)
      else
        do m = 1,nq
          ucom(m) = lufp(m,k)*(half - dot(LDG_beta(1:nr),+unrm(1:nr))) + &
                    rufp(m,k)*(half - dot(LDG_beta(1:nr),-unrm(1:nr)))
        end do
      end if
      !
      do m = 1,nq
        ljfp(m) = ucom(m) - lufp(m,k)
      end do
      !
     !do m = 1,nq
     !  ljfp(m) = hfpb*(rufp(m,k) - lufp(m,k))
     !end do
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
#ifdef USE_MATRIX_PROCEDURES
      do i = correct(host_geom,host_order)%row_ptr(kl)+1, &
             correct(host_geom,host_order)%row_ptr(kl+1)
        n = l1 - 1 + correct(host_geom,host_order)%row_idx(i)
        correct_n = dsl * correct(host_geom,host_order)%row_val(i) &
                        / jacobian(n)
        do m = 1,nq
          correct_m = correct_n * ljfp(m)
          do l = 1,nr
            dusp(l,m,n) = dusp(l,m,n) + correct_m * unrm(l)
          end do
        end do
      end do
#else
      do n = l1,l2
        cmat_val = correct(host_geom,host_order)%mat(n-l1+1,kl)
        if (cmat_val /= zero) then
          correct_n = dsl * cmat_val / jacobian(n)
          do m = 1,nq
            correct_m = correct_n * ljfp(m)
            do l = 1,nr
              dusp(l,m,n) = dusp(l,m,n) + correct_m * unrm(l)
            end do
          end do
        end if
      end do
#endif
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
      kl = face(nf)%left%flxpts(k)
      kr = face(nf)%right%flxpts(k)
      !
      ! Get the unit normal at this flux point
      !
      unrm(1:nr) = unit_vector( geofa(nf)%v(1:nr,k) )
      !
      ! Get the face length for the current face
      !
      ds = rnorm2( geofa(nf)%v(1:nr,k) )
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
        ucom(m) = lufp(m,k)*(half - dot(LDG_beta(1:nr),+unrm(1:nr))) + &
                  rufp(m,k)*(half - dot(LDG_beta(1:nr),-unrm(1:nr)))
      end do
      !
      do m = 1,nq
        ljfp(m) = ucom(m) - lufp(m,k)
        rjfp(m) = ucom(m) - rufp(m,k)
      end do
      !
     !do m = 1,nq
     !  ljfp(m) = hfpb*(rufp(m,k) - lufp(m,k))
     !  rjfp(m) = hfmb*(lufp(m,k) - rufp(m,k))
     !end do
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
#ifdef USE_MATRIX_PROCEDURES
      do i = correct(left_geom,left_order)%row_ptr(kl)+1, &
             correct(left_geom,left_order)%row_ptr(kl+1)
        n = l1 - 1 + correct(left_geom,left_order)%row_idx(i)
        correct_n = dsl * correct(left_geom,left_order)%row_val(i) &
                        / jacobian(n)
        do m = 1,nq
          correct_m = correct_n * ljfp(m)
          do l = 1,nr
            dusp(l,m,n) = dusp(l,m,n) + correct_m * unrm(l)
          end do
        end do
      end do
#else
      do n = l1,l2
        cmat_val = correct(left_geom,left_order)%mat(n-l1+1,kl)
        if (cmat_val /= zero) then
          correct_n = dsl * cmat_val / jacobian(n)
          do m = 1,nq
            correct_m = correct_n * ljfp(m)
            do l = 1,nr
              dusp(l,m,n) = dusp(l,m,n) + correct_m * unrm(l)
            end do
          end do
        end if
      end do
#endif
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the right cell
      !
#ifdef USE_MATRIX_PROCEDURES
      do i = correct(right_geom,right_order)%row_ptr(kr)+1, &
             correct(right_geom,right_order)%row_ptr(kr+1)
        n = r1 - 1 + correct(right_geom,right_order)%row_idx(i)
        correct_n = dsr * correct(right_geom,right_order)%row_val(i) &
                        / jacobian(n)
        do m = 1,nq
          correct_m = correct_n * rjfp(m)
          do l = 1,nr
            dusp(l,m,n) = dusp(l,m,n) - correct_m * unrm(l)
          end do
        end do
      end do
#else
      do n = r1,r2
        cmat_val = correct(right_geom,right_order)%mat(n-r1+1,kr)
        if (cmat_val /= zero) then
          correct_n = dsr * cmat_val / jacobian(n)
          do m = 1,nq
            correct_m = correct_n * rjfp(m)
            do l = 1,nr
              dusp(l,m,n) = dusp(l,m,n) - correct_m * unrm(l)
            end do
          end do
        end if
      end do
#endif
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
  use generic_types_mod, only : matrix
  use order_mod,  only : maxSP,maxFP
  use geovar,  only : cell
  use geovar,  only : nr,ncell,nfbnd,bface
  use geovar,  only : face
  use metrics_mod,  only : geofa,jacobian
  use correction_mod, only : correct
  use flowvar, only : dusp,facedusp,faceusp
  use ovar,    only : walls_are_exact,LDG_beta
  use interpolation_mod, only : interp
  !
  use parallel_mod, only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,l,m
  integer :: nf,kf,lf,sf,ncf
  integer :: this_geom,this_order,face_order
  real(wp) :: dsl,cbias,hfpb,correct_n,correct_m
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:nq) :: lufp,rufp,ujump
  real(wp), dimension(1:nr,1:nq,1:maxSP) :: cell_dusp
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: face_dusp
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: dusp_copy
  !
#ifdef USE_MATRIX_PROCEDURES
  integer :: i,i1,i2
  type(matrix) :: cmat
#else
  integer :: nfp
  real(wp), dimension(1:maxSP) :: lmat
  real(wp), dimension(1:maxSP,1:maxFP) :: cmat_jac
#endif
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
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
#ifdef USE_MATRIX_PROCEDURES
    cmat = correct(this_geom,this_order)
    !
   !do i = 1,size(cmat%mat,dim=2)
    do i = 1,size(cmat%row_ptr)-1
      i1 = cmat%row_ptr(i)+1
      i2 = cmat%row_ptr(i+1)
      cmat%row_val(i1:i2) = cmat%row_val(i1:i2) / &
                            jacobian(n1 + cmat%row_idx(i1:i2) - 1)
    end do
#else
    nfp = size(correct(this_geom,this_order)%mat,dim=2)
    cmat_jac(1:np,1:nfp) = correct(this_geom,this_order)%mat
    !
    do kf = 1,nfp
      cmat_jac(1:np,kf) = cmat_jac(1:np,kf) / jacobian(n1:n2)
    end do
#endif
    !
    ! Create a local copy of dusp for this cell to which all the
    ! corrections will be applied. This will help reduce cache misses.
    !
    cell_dusp(1:nr,1:nq,1:np) = dusp(1:nr,1:nq,n1:n2)
    !
    ! Create a transposed copy of the uncorrected dusp for this cell which
    ! will be used to reset the face_dusp array in the cell_faces loop. This
    ! will help reduce cache misses.
    !
    do n = 1,np
      dusp_copy(n,1:nr,1:nq) = cell_dusp(1:nr,1:nq,n)
    end do
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
      face_order = cell(nc)%face(ncf)%order ! order of current face
      !
      ! Reset face_dusp to the original discontinuous gradient
      !
!###############################################################################
!###############################################################################
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
!-------------------------------------------------------------------------------
     !face_dusp = dusp_copy
!-------------------------------------------------------------------------------
      face_dusp(1:np,1:nr,1:nq) = dusp_copy(1:np,1:nr,1:nq)
!-------------------------------------------------------------------------------
!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!###############################################################################
!###############################################################################
      !
      face_points: do kf = 1,size(cell(nc)%face(ncf)%cell_idx)
        !
        !kf: index of flx-pt local to this face
        lf = cell(nc)%face(ncf)%cell_idx(kf) ! index of flx-pt within interp%mat
        !
        ! Solution at this flux point for the left and right cells
        !
        lufp(:) = faceusp(nf)%v(:,kf,1)
        rufp(:) = faceusp(nf)%v(:,kf,2)
        !
        ! Get the unit normal at this flux point
        !
        unrm(1:nr) = unit_vector( geofa(nf)%v(1:nr,kf) )
        !
        ! Get the face length for the current face
        !
        dsl = cl(this_geom) * rnorm2( geofa(nf)%v(1:nr,kf) )
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
        cbias = zero
       !cbias = dot( LDG_beta(1:nr) , unrm(1:nr) )
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
#ifdef USE_MATRIX_PROCEDURES
       !i1 = cmat%row_ptr(lf)+1
       !i2 = cmat%row_ptr(lf+1)
       !do i = i1,i2
        do i = cmat%row_ptr(lf)+1, cmat%row_ptr(lf+1)
          n = cmat%row_idx(i)
          correct_n = dsl * cmat%row_val(i)
          do m = 1,nq
            correct_m = correct_n * ujump(m)
            do l = 1,nr
              face_dusp(n,l,m) = face_dusp(n,l,m) + correct_m * unrm(l)
              cell_dusp(l,m,n) = cell_dusp(l,m,n) + correct_m * unrm(l)
            end do
          end do
        end do
#else
        do n = 1,np
         !correct_n = dsl * correct(this_geom,this_order)%mat(n,lf) / &
         !                  jacobian(n+n1-1)
          if (cmat_jac(n,lf) /= zero) then
            correct_n = dsl * cmat_jac(n,lf)
            do m = 1,nq
              correct_m = correct_n * ujump(m)
              do l = 1,nr
                face_dusp(n,l,m) = face_dusp(n,l,m) + correct_m * unrm(l)
                cell_dusp(l,m,n) = cell_dusp(l,m,n) + correct_m * unrm(l)
              end do
            end do
          end if
        end do
#endif
        !
      end do face_points
      !
      ! Interpolate the face corrected gradient from the interior solution
      ! points to the interface flux points now that we are done with all
      ! the flux points on the current face of the current cell.
      !
      interpolate_face: do kf = 1,size(cell(nc)%face(ncf)%cell_idx)
        !
        !kf: index of flx-pt local to this face
        lf = cell(nc)%face(ncf)%cell_idx(kf) ! index of flx-pt within interp%mat
        !
#ifdef USE_MATRIX_PROCEDURES
        do m = 1,nq
          do l = 1,nr
            facedusp(nf)%v(l,m,kf,sf) = interp(this_geom,this_order)% &
                                                  toFace(face_order)% &
                                        dot( lf ,face_dusp(1:np,l,m) )
          end do
        end do
#else
        lmat(1:np) = interp(this_geom,this_order)% &
                               toFace(face_order)% &
                     mat(1:np,lf)
        !
        do m = 1,nq
          do l = 1,nr
            facedusp(nf)%v(l,m,kf,sf) = dot( lmat(1:np) , face_dusp(1:np,l,m) )
          end do
        end do
#endif
        !
      end do interpolate_face
      !
    end do cell_faces
    !
    ! Copy the total corrected gradient accounting for all cell faces
    ! back into dusp
    !
    dusp(1:nr,1:nq,n1:n2) = cell_dusp(1:nr,1:nq,1:np)
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
! Procedure for correcting the solution and gradient before output
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine correct_solution_for_output(interp_usp,interp_dusp)
  !
  use generic_types_mod, only : matrix
  use order_mod,  only : maxFP,geom_solpts
  use geovar,  only : nr,nfbnd,nface,n_solpts
  use geovar,  only : bface,face,cell
  use ovar,    only : invs_flux_method,visc_flux_method
  use metrics_mod,  only : geofa,jacobian
  use correction_mod, only : correct_sol
  use flowvar, only : usp,dusp,uoldsp,faceusp
  !
  use module_bc, only : get_faceusp_bc
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),   intent(inout) :: interp_usp
  real(wp), dimension(:,:,:), intent(inout) :: interp_dusp
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nf,r1,r2,l1,l2
  integer  :: kl,kr,nfp
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: ds,dsl,dsr,correct_n
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq) :: ucom,ljfp,rjfp
  real(wp), dimension(1:nr,1:maxFP) :: rnrm
  real(wp), dimension(1:nq,1:maxFP) :: lufp,rufp
  !
  real(wp), parameter :: addsub_cor =  one
 !real(wp), parameter :: addsub_cor = -one
  !
 !logical, parameter :: use_jacobian = true
  logical, parameter :: use_jacobian = fals
  !
  real(wp), parameter :: dscoef = one
 !real(wp), parameter :: dscoef = half
  !
  integer, parameter :: dsexp = 1
 !integer, parameter :: dsexp = 0
  !
  real(wp), parameter :: left_coef =  one
 !real(wp), parameter :: left_coef = -one
  !
  real(wp), parameter :: right_coef =  one
 !real(wp), parameter :: right_coef = -one
  !
  integer :: iloops
  integer, parameter :: npasses = 1
  !
continue
  !
  ! Save the updated solution into uoldsp since this has yet to be done
  !
  uoldsp(:,1:n_solpts) = usp(:,1:n_solpts)
  !
  npass_loop: do iloops = 1,npasses
    !
  call interp_cv_to_flxpts
  call get_faceusp_bc
  call discontinuous_gradient_at_solpts
 !if (visc_flux_method == visc_flux_BR2) then
 !  call correct_gradients_br2
 !else
    call correct_gradients_at_solpts
 !end if
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
    lufp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,1)
    ! FP solution from right cell on face
    rufp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,2)
    !
    rnrm(1:nr,1:nfp) = geofa(nf)%v(1:nr,1:nfp) ! Area face normals at FP
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
      kl = face(nf)%left%flxpts(k)
      !
      ! Get the contribution to the face length/area from the current flux
      ! point for the left cell
      !
      ds  = ( dscoef * rnorm2(rnrm(:,k)) )**dsexp
      !
      dsl = ds * cl(host_geom)
      !
      ! Compute the common interface solutions
      !
      if (any(bface(1,nf) == bc_walls)) then
        ucom(1:nq) = rufp(1:nq,k)
      else
        ucom(1:nq) = half * (lufp(1:nq,k) + rufp(1:nq,k))
      end if
      !
      do m = 1,nq
        ljfp(m) = addsub_cor * dsl * (ucom(m) - lufp(m,k))
      end do
      !
      ! Now use the correction function with the interface jump
      ! to correct the conservative variables in the left cell
      !
      do i = correct_sol(host_geom,host_order)%row_ptr(kl)+1, &
             correct_sol(host_geom,host_order)%row_ptr(kl+1)
        n = l1 - 1 + correct_sol(host_geom,host_order)%row_idx(i)
        correct_n = correct_sol(host_geom,host_order)%row_val(i)
        if (use_jacobian) correct_n = correct_n / jacobian(n)
        do m = 1,nq
          usp(m,n) = usp(m,n) + ljfp(m) * correct_n * left_coef
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
 !interior_fluxes: do nf = nfbnd+1,nface
  interior_fluxes: do nf = nface,nfbnd+1,-1
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
    lufp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,1)
    rufp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,2)
    !
    rnrm(1:nr,1:nfp) = geofa(nf)%v(1:nr,1:nfp)
    !
    ! Apply the flux correction due to the flux jump
    !
    do k = 1,nfp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = face(nf)%left%flxpts(k)
      kr = face(nf)%right%flxpts(k)
      !
      ! Get the face length for the left and right cell
      !
      ds  = ( dscoef * rnorm2(rnrm(:,k)) )**dsexp
      !
      dsl = ds * cl(left_geom)
      dsr = ds * cl(right_geom)
      !
      ucom(1:nq) = half * (lufp(1:nq,k) + rufp(1:nq,k))
      !
      do m = 1,nq
        ljfp(m) = addsub_cor * dsl * (ucom(m) - lufp(m,k))
        rjfp(m) = addsub_cor * dsr * (ucom(m) - rufp(m,k))
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! left_cell
      !
      do i = correct_sol(left_geom,left_order)%row_ptr(kl)+1, &
             correct_sol(left_geom,left_order)%row_ptr(kl+1)
        n = l1 - 1 + correct_sol(left_geom,left_order)%row_idx(i)
        correct_n = correct_sol(left_geom,left_order)%row_val(i)
        if (use_jacobian) correct_n = correct_n / jacobian(n)
        do m = 1,nq
          usp(m,n) = usp(m,n) + ljfp(m) * correct_n * left_coef
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! right_cell
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      do i = correct_sol(right_geom,right_order)%row_ptr(kr)+1, &
             correct_sol(right_geom,right_order)%row_ptr(kr+1)
        n = r1 - 1 + correct_sol(right_geom,right_order)%row_idx(i)
        correct_n = correct_sol(right_geom,right_order)%row_val(i)
        if (use_jacobian) correct_n = correct_n / jacobian(n)
        do m = 1,nq
          usp(m,n) = usp(m,n) + rjfp(m) * correct_n * right_coef
        end do
      end do
      !
    end do
    !
  end do interior_fluxes
  !
  end do npass_loop
  !
  ! Copy the corrected solution and derivatives into interp_usp and interp_dusp
  !
  interp_usp(:,1:n_solpts) = usp(:,1:n_solpts)
  interp_dusp(:,:,1:n_solpts) = dusp(:,:,1:n_solpts)
  !
  ! Copy the saved updated solution back into usp
  !
  usp(:,1:n_solpts) = uoldsp(:,1:n_solpts)
  !
end subroutine correct_solution_for_output
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
subroutine flux_gradient_at_solpts_with_qdtr
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use order_mod, only : maxSP,maxQP,q_order,geom_solpts
  use geovar,  only : nr,ncell,cell,using_quadrature
  use projection_mod, only : projct
  use interpolation_mod, only : interp
  use derivatives_mod, only : deriv,dqmat
  use metrics_mod,  only : metrics,geofa,jacobian
  use flowvar, only : usp,residual,faceflx,faceusp
  use ovar,    only : dflux_opt
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nc,n1,n2,nsp,nqp
  integer  :: nf,sf,kf,lf
  integer  :: this_geom,this_order,face_order
  real(wp) :: drdxJ,drdyJ,dsdxJ,dsdyJ
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)    :: unrm
  real(wp), dimension(1:maxQP) :: qjac
  real(wp), dimension(1:maxSP) :: lmat
  !
  real(wp), dimension(1:maxQP,1:maxSP) :: pmat
  real(wp), dimension(1:maxQP,1:maxSP) :: dmat
  !
  real(wp), dimension(1:maxSP,1:nq) :: tusp
  real(wp), dimension(1:nq,1:maxQP) :: uqp
  !
  real(wp), dimension(1:nr,1:nr,1:maxQP) :: dxdr
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: rsflx_sp
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: xyflx_sp
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: dflx_sp
  real(wp), dimension(1:maxQP,1:nq,1:nr) :: rsflx_qp
  real(wp), dimension(1:maxQP,1:nq,1:nr) :: xyflx_qp
  !
continue
  !
  ! Initialize faceflx to zero before we start
  !
  call faceflx%zero
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index for sol-pts in cell nc
    n2 = cell(nc)%end_sp ! ending    index for sol-pts in cell nc
    nsp = n2-n1+1 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! no. of quadrature points for this cell
    nqp = geom_solpts(this_geom,q_order)
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
        lmat(1:nsp) = interp(this_geom,this_order)%toCell(q_order)%mat(1:nsp,k)
        !
        do m = 1,nq
          uqp(m,k) = dot( lmat(1:nsp) , tusp(1:nsp,m) )
        end do
        !
        ! Get the metrics for this cell
        !
        do m = 1,nr
          do l = 1,nr
            dxdr(l,m,k) = dot( lmat(1:nsp) , metrics(l,m,n1:n2,2) )
          end do
        end do
        qjac(k) = dot( lmat(1:nsp) , jacobian(n1:n2) )
        !
      end do
      !
    else
      !
      uqp(1:nq,1:nsp) = usp(1:nq,n1:n2)
      !
      dxdr(1:nr,1:nr,1:nsp) = metrics(1:nr,1:nr,n1:n2,2)
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
      ! Transform the fluxes in the x-y physical
      ! domain to the r-s computational domain
      ! NOTE : The 'J' at the end of each these variables means that the
      !        quantity is implicitly multiplied by the Jacobian, e.g.,
      !        dsdyJ = dsdy*Jacobian = dxdr
      !
      do m = 1,nq
        do l = 1,nr
          rsflx_qp(k,m,l) = dot( xyflx_qp(k,m,1:nr) , dxdr(1:nr,l,k) )
        end do
      end do
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
        pmat(1:nqp,1:nsp) = projct(this_geom,q_order)% &
                                   toCell(this_order)%mat(1:nqp,1:nsp)
        !
        do l = 1,nr
          do m = 1,nq
            do k = 1,nsp
              xyflx_sp(k,m,l) = dot( pmat(1:nqp,k) , xyflx_qp(1:nqp,m,l) )
              rsflx_sp(k,m,l) = dot( pmat(1:nqp,k) , rsflx_qp(1:nqp,m,l) )
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
          do l = 1,nr
            dflx_sp(k,m,l) = dot( deriv(this_geom,this_order)%d(l)% &
                                  mat(1:nsp,k) , rsflx_sp(1:nsp,m,l) )
          end do
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
      do n = 1,size(cell(nc)%face)
        !
        nf = cell(nc)%face(n)%idx  ! index of the current face
        sf = cell(nc)%face(n)%side ! side of face that cell nc is on
        !
        face_order = cell(nc)%face(n)%order ! order of current face
        !
        do kf = 1,size(cell(nc)%face(n)%cell_idx) ! index of the flx-pt
          !                                         local to this face
          lf = cell(nc)%face(n)%cell_idx(kf) ! index of flx-pt within interp%mat
          !
          ! Get the unit normal for this face flux/solution point
          !
          unrm(:) = unit_vector( geofa(nf)%v(:,kf) )
          !
          !#####################################################################
          !#####################################################################
          ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
          ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
          ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
          !#####################################################################
          !#####################################################################
          !
          lmat(1:nsp) = interp(this_geom,this_order)% &
                                  toFace(face_order)% &
                        mat(1:nsp,lf)
          !
          do l = 1,nr
            do m = 1,nq
              faceflx(nf)%v(m,kf,sf) = faceflx(nf)%v(m,kf,sf) + &
                              dot( lmat(1:nsp) , xyflx_sp(1:nsp,m,l) ) * unrm(l)
            end do
          end do
          !
        end do
        !
      end do
      !
    else if (dflux_opt == project_derivatives) then
      !
      ! Compute the derivative of the fluxes at the quadrature points
      ! and project these derivatives to the solution points
      !
      do l = 1,nr
        dmat(1:nqp,1:nsp) = dqmat(this_geom,q_order)%d(l)%mat(1:nqp,1:nsp)
        do m = 1,nq
          do k = 1,nsp
            dflx_sp(k,m,l) = dot( dmat(1:nqp,k) , rsflx_qp(1:nqp,m,l) )
          end do
        end do
      end do
      !
     !! Interpolate the fluxes to the interface flux points and dot this
     !! interpolated vector with the face normal at each flux point
     !! NOTE: We dont need to exchange this information for communication and
     !!       periodic boundary faces since this flux is only used within the
     !!       current cell.
     !!
     !! Loop over the faces of this cell
     !!
     !do n = 1,size(cell(nc)%face)
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
     !  nf = cell(nc)%face(n)%idx  ! index of the current face
     !  sf = cell(nc)%face(n)%side ! side of face that cell nc is on
     !  !
     !  face_order = cell(nc)%face(n)%order ! order of current face
     !  !
     !  do kf = 1,size(cell(nc)%face(n)%cell_idx) ! index of the flx-pt
     !    !                                         local to this face
     !    lf = cell(nc)%face(n)%cell_idx(kf) ! index of flx-pt within interp%mat
     !    !
     !    ! Get the unit normal for this face flux/solution point
     !    !
     !    unrm(:) = unit_vector( geofa(nf)%v(:,kf) )
     !    !
     !    !#####################################################################
     !    !#####################################################################
     !    ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
     !    ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
     !    ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
     !    !#####################################################################
     !    !#####################################################################
     !    !
     !    lmat(1:nsp) = interp(this_geom,this_order)% &
     !                               toFace(q_order)% &
     !                  mat(1:nsp,lf)
     !    !
     !    do l = 1,nr
     !      do m = 1,nq
     !        faceflx(nf)%v(m,kf,sf) = faceflx(nf)%v(m,kf,sf) + &
     !                        dot( lmat(1:nsp) , xyflx_qp(1:nsp,m,l) ) * unrm(l)
     !      end do
     !    end do
     !    !
     !  end do
     !  !
     !end do
     !!
    end if
    !
    ! Add the derivatives of the flux to the residual
    !
    do m = 1,nq
      do k = 1,nsp
        residual(m,k+n1-1) = sum( dflx_sp(k,m,1:nr) ) 
      end do
    end do
    !
   !do k = 1,nsp
   !  do m = 1,nq
   !   !
   !   !residual(m,k+n1-1) = dot( deriv(this_geom,this_order)% &
   !   !                                                 d(1)% &
   !   !                          mat(1:nsp,k) , ftdc(1:nsp,m) )
   !    residual(m,k+n1-1) = dot( deriv(this_geom,this_order)% &
   !                                                     d(1)% &
   !                              mat(1:nsp,k) , dfsp(1:nsp,m) )
   !    do l = 2,nr
   !   !  residual(m,k+n1-1) = residual(m,k+n1-1) + &
   !   !                       dot( deriv(this_geom,this_order)% &
   !   !                                                   d(l)% &
   !   !                            mat(1:nsp,k) , gtdc(1:nsp,m) )
   !      residual(m,k+n1-1) = residual(m,k+n1-1) + &
   !                           dot( deriv(this_geom,this_order)% &
   !                                                       d(l)% &
   !                                mat(1:nsp,k) , dgsp(1:nsp,m) )
   !    end do
   !  end do
   !end do
    !
  end do
  !
end subroutine flux_gradient_at_solpts_with_qdtr
!
!###############################################################################
!
subroutine correct_interface_flux_with_qdtr()
  !
  use order_mod,  only : maxFP,maxGP,q_order,geom_solpts
  use geovar,  only : nfbnd,nface,nr
  use geovar,  only : bface,face,cell,xyz
  use geovar,  only : using_quadrature
  use ovar,    only : flux_divergence,dflux_opt,itcur,nst
  use ovar,    only : invs_flux_method
  use flowvar, only : usp,residual,faceflx,faceusp
  use metrics_mod,  only : geofa,volj
  use correction_mod, only : correct
  use interpolation_mod, only : interp
  use projection_mod, only : projct
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nf,r1,r2,l1,l2
  integer  :: nfp,ngp,kl,kr,iu
  integer  :: face_geom,face_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:nq) :: fnl,fnr,fnc
  real(wp), dimension(1:nq) :: srcl,srcr
  !
  real(wp), dimension(1:maxFP) :: lmat
  !
  real(wp), dimension(1:nr,1:maxGP) :: rnrm
  real(wp), dimension(1:nq,1:maxFP) :: fjl,fjr
  real(wp), dimension(1:nq,1:maxGP) :: lflx,rflx
  real(wp), dimension(1:nq,1:maxGP) :: luqp,ruqp,fncq
  !
  real(wp), dimension(1:maxGP,1:maxFP) :: pmat
  !
continue
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    nfp = geom_solpts(face_geom,face_order) ! number of solution   flux points
    ngp = geom_solpts(face_geom,q_order)    ! number of quadrature flux points
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
    ! Interpolate the conserved variables and face normals at
    ! the solution flux points to the quadrature flux points
    !
    if (using_quadrature) then
      !
      do k = 1,ngp
        !
        lmat(1:nfp) = interp(face_geom,face_order)%toCell(q_order)%mat(1:nfp,k)
        !
        do m = 1,nq
          luqp(m,k) = dot( lmat(1:nfp) , faceusp(nf)%v(m,1:nfp,1) )
          ruqp(m,k) = dot( lmat(1:nfp) , faceusp(nf)%v(m,1:nfp,2) )
        end do
        !
        do l = 1,nr
          rnrm(l,k) = dot( lmat(1:nfp) , geofa(nf)%v(l,1:nfp) )
        end do
        !
      end do
      !
    else
      !
      lflx(1:nq,1:nfp) = faceflx(nf)%v(1:nq,1:nfp,1)
      rflx(1:nq,1:nfp) = faceflx(nf)%v(1:nq,1:nfp,2)
      luqp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,1)
      ruqp(1:nq,1:nfp) = faceusp(nf)%v(1:nq,1:nfp,2)
      rnrm(1:nr,1:nfp) = geofa(nf)%v(1:nr,1:nfp)
      !
    end if
    !
    ! Get the discontinuous fluxes at the interface
    !
    if (dflux_opt == project_derivatives) then
      !
      ! Get the discontinuous fluxes at the quadrature flux points
      !
      do k = 1,ngp
        lflx(1:nq,k) = normal_invs_flx_cv( luqp(:,k) , rnrm(:,k) )
        rflx(1:nq,k) = normal_invs_flx_cv( ruqp(:,k) , rnrm(:,k) )
      end do
      !
    else
      !
      ! Get the discontinuous fluxes at the solution flux points
      !
      lflx(1:nq,1:nfp) = faceflx(nf)%v(1:nq,1:nfp,1)
      rflx(1:nq,1:nfp) = faceflx(nf)%v(1:nq,1:nfp,2)
      !
    end if
    !
    ! Get the common flux at each quadrature flux point
    !
    do k = 1,ngp
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( rnrm(:,k) )
      !
      ! Get the numerical flux at the current flux point
      !
      if (invs_flux_method == invs_flux_ROE) then
        fncq(:,k) = upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (invs_flux_method == invs_flux_HLLC) then
        fncq(:,k) = hllc_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (invs_flux_method == invs_flux_LDFSS) then
        fncq(:,k) = ldfss_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      end if
      !
    end do
    !
    ! Project the common flux to the solution flux points
    !
    if (using_quadrature) then
      !
      pmat(1:ngp,1:nfp) = projct(face_geom,q_order)% &
                                 toCell(face_order)% &
                          mat(1:ngp,1:nfp)
      !
      if (dflux_opt == project_derivatives) then
        !
        lflx(1:nq,1:ngp) = fncq(1:nq,1:ngp) - lflx(1:nq,1:ngp)
        rflx(1:nq,1:ngp) = fncq(1:nq,1:ngp) - rflx(1:nq,1:ngp)
        !
        do k = 1,nfp
          !
          do m = 1,nq
            fjl(m,k) = dot( pmat(1:ngp,k) , lflx(m,1:ngp) )
            fjr(m,k) = dot( pmat(1:ngp,k) , rflx(m,1:ngp) )
          end do
          !
        end do
        !
      else
        !
        do k = 1,nfp
          !
          do m = 1,nq
            fnc(m) = dot( pmat(1:ngp,k) , fncq(m,1:ngp) )
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
      fjl(1:nq,1:nfp) = fncq(1:nq,1:nfp) - lflx(1:nq,1:nfp)
      fjr(1:nq,1:nfp) = fncq(1:nq,1:nfp) - rflx(1:nq,1:nfp)
      ! 
    end if
    !
    ! Compute the flux jump at each solution flux point
    !
   !do k = 1,nfp
   !  fjl(:,k) = fnc(:,k) - lflx(:,k)
   !  fjr(:,k) = fnc(:,k) - rflx(:,k)
   !end do
   !!
   !! Initialize the conservation source terms
   !!
   !srcl(:) = zero
   !srcr(:) = zero
   !!
    do k = 1,nfp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = face(nf)%left%flxpts(k)
      kr = face(nf)%right%flxpts(k)
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
          residual(m,n) = residual(m,n) + &
               dsl * fjl(m,k) * correct(left_geom,left_order)%mat(n-l1+1,kl)
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
            residual(m,n) = residual(m,n) - &
                 dsr * fjr(m,k) * correct(right_geom,right_order)%mat(n-r1+1,kr)
          end do
        end do
      end if
      !
      ! From ZJ, to ensure conservation with chain-rule approach
      !
     !if (flux_divergence == ChainRuleCV) then
     !  do m = 1,nq   ! ARE DSL AND DSR CORRECT HERE?
     !    srcl(m) = srcl(m) + dsl*lflx(m,k)*flxwts(k)
     !    srcr(m) = srcr(m) + dsr*rflx(m,k)*flxwts(k)
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
end subroutine correct_interface_flux_with_qdtr
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
include "Functions/viscosity.f90"
!
!###############################################################################
!
end module flux_mod
