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
contains
!
!###############################################################################
!
subroutine interp_cv_to_flxpts()
  !
  use order_mod,  only : maxpts
  use geovar,  only : ncell
  use geovar,  only : pts_in_cell_ptr
  use geovar,  only : flxpts_of_cell_ptr
  use geovar,  only : flxpts_of_cell
  use geovar,  only : cell_geom,cell_order
  use interpolation_mod, only : interp
  use flowvar, only : usp,faceusp
  !
  use parallel_mod, only : exchange_faceusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,m
  integer :: nf,f1,f2,kf,lf,sf
  integer :: this_geom,this_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxpts) :: lmat
  real(wp), dimension(1:maxpts,1:nq) :: tusp
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1 ! beginning index of sol-pts for this cell
    n2 = pts_in_cell_ptr(nc+1) ! ending    index of sol-pts for this cell
    np = n2-n1+1               ! number of sol-pts in this cell
    !
    this_geom = cell_geom(nc)   ! geometry type of this cell
    this_order = cell_order(nc) ! solution order of this cell
    !
    f1 = flxpts_of_cell_ptr(nc)+1 ! beginning index of flx-pts for this cell
    f2 = flxpts_of_cell_ptr(nc+1) ! ending    index of flx-pts for this cell 
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
    do n = f1,f2
      !
      nf = flxpts_of_cell(1,n) ! index of the face that this flux point is on
      kf = flxpts_of_cell(2,n) ! index of flux point local to this face
      lf = flxpts_of_cell(3,n) ! index of flux point local to this cell
      sf = flxpts_of_cell(4,n) ! side of face cell nc is on
      !
      !######################################################################
      !######################################################################
      ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
      ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
      ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
      !######################################################################
      !######################################################################
      !
      lmat(1:np) = interp(this_geom,this_order)%toFace(this_order)%mat(1:np,lf)
      !
      do m = 1,nq
        faceusp(m,kf,sf,nf) = dot( lmat(1:np) , tusp(1:np,m) )
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
  use order_mod,  only : maxpts
  use geovar,  only : nr,ncell
  use geovar,  only : pts_in_cell_ptr
  use geovar,  only : flxpts_of_cell_ptr
  use geovar,  only : flxpts_of_cell
  use geovar,  only : cell_geom,cell_order
  use interpolation_mod, only : interp
  use flowvar, only : dusp,facedusp
  !
  use parallel_mod, only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,l,m
  integer :: nf,f1,f2,kf,lf,sf
  integer :: fp1,fp2
  integer :: this_geom,this_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxpts) :: lmat
  real(wp), dimension(1:maxpts,1:nr,1:nq) :: tdusp
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1 ! beginning index of sol-pts for this cell
    n2 = pts_in_cell_ptr(nc+1) ! ending    index of sol-pts for this cell
    np = n2-n1+1               ! number of sol-pts in this cell
    !
    this_geom = cell_geom(nc)   ! geometry type of this cell
    this_order = cell_order(nc) ! solution order of this cell
    !
    f1 = flxpts_of_cell_ptr(nc)+1 ! beginning index of flx-pts for this cell
    f2 = flxpts_of_cell_ptr(nc+1) ! ending    index of flx-pts for this cell 
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
    do n = f1,f2
      !
      nf = flxpts_of_cell(1,n) ! index of the face that this flux point is on
      kf = flxpts_of_cell(2,n) ! index of flux point local to this face
      lf = flxpts_of_cell(3,n) ! index of flux point local to this cell
      sf = flxpts_of_cell(4,n) ! side of face cell nc is on
      !
      !######################################################################
      !######################################################################
      ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
      ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
      ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
      !######################################################################
      !######################################################################
      !
      lmat(1:np) = interp(this_geom,this_order)%toFace(this_order)%mat(1:np,lf)
      !
      do m = 1,nq
        do l = 1,nr
          facedusp(l,m,kf,sf,nf) = dot( lmat(1:np) , tdusp(1:np,l,m) )
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
  use order_mod, only : maxpts,q_order,geom_solpts
  use geovar,  only : nr,ncell,pts_in_cell_ptr,using_quadrature
  use geovar,  only : flxpts_of_cell,flxpts_of_cell_ptr
  use geovar,  only : cell_geom,cell_order
  use projection_mod, only : projct
  use interpolation_mod, only : interp
  use derivatives_mod, only : deriv,dqmat
  use metrics_mod,  only : metrics,geofa,jacobian
  use flowvar, only : usp,residual,faceflx,faceusp
  use ovar,    only : dflux_opt
  !
  !.. Local Scalars ..
  integer :: k,l,m,nc,n1,n2,nsp,nqp
  integer :: f1,f2,nf,sf,kf,lf
  integer :: this_geom,this_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)     :: unrm
  real(wp), dimension(1:maxpts) :: dxdr,dxds
  real(wp), dimension(1:maxpts) :: dydr,dyds
  real(wp), dimension(1:maxpts) :: qjac
  real(wp), dimension(1:maxpts) :: lmat
  !
  real(wp), dimension(1:maxpts,1:maxpts) :: pmat
  real(wp), dimension(1:maxpts,1:maxpts) :: dmat
  !
  real(wp), dimension(1:maxpts,1:nq) :: tusp
  real(wp), dimension(1:nq,1:maxpts) :: uqp
  !
  real(wp), dimension(1:maxpts,1:nq,1:nr) :: rsflx_sp
  real(wp), dimension(1:maxpts,1:nq,1:nr) :: xyflx_sp
  real(wp), dimension(1:maxpts,1:nq,1:nr) :: dflx_sp
  real(wp), dimension(1:maxpts,1:nq,1:nr) :: rsflx_qp
  real(wp), dimension(1:maxpts,1:nq,1:nr) :: xyflx_qp
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
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1 ! beginning index for sol-pts in cell nc
    n2 = pts_in_cell_ptr(nc+1) ! ending    index for sol-pts in cell nc
    nsp = n2-n1+1 ! number of solution points in this cell
    !
    this_geom = cell_geom(nc)
    this_order = cell_order(nc)
    !
    nqp = geom_solpts(this_geom,q_order) ! no. of quadrature points for this cell
    !
    f1 = flxpts_of_cell_ptr(nc)+1
    f2 = flxpts_of_cell_ptr(nc+1)
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
        dxdr(k) = dot( lmat(1:nsp) , metrics(1,1,n1:n2) )
        dxds(k) = dot( lmat(1:nsp) , metrics(1,2,n1:n2) )
        dydr(k) = dot( lmat(1:nsp) , metrics(2,1,n1:n2) )
        dyds(k) = dot( lmat(1:nsp) , metrics(2,2,n1:n2) )
        qjac(k) = dot( lmat(1:nsp) , jacobian(n1:n2) )
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
        pmat(1:nqp,1:nsp) = projct(this_geom,q_order)%toCell(this_order)%mat(1:nqp,1:nsp)
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
            dflx_sp(k,m,l) = dot( deriv(this_geom,this_order)%d(l)%mat(1:nsp,k) , &
                                  rsflx_sp(1:nsp,m,l) )
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
        lmat(1:nsp) = interp(this_geom,this_order)%toFace(this_order)%mat(1:nsp,lf)
        !
        do l = 1,nr
          do m = 1,nq
            faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
                                  dot( lmat(1:nsp) , xyflx_sp(1:nsp,m,l) ) * unrm(l)
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
      do l = 1,nr
        dmat(1:nqp,1:nsp) = dqmat(this_geom,q_order)%d(l)%mat(1:nqp,1:nsp)
        do m = 1,nq
          do k = 1,nsp
            dflx_sp(k,m,l) = dot( dmat(1:nqp,k) , rsflx_qp(1:nqp,m,l) )
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
     !  lmat(1:nsp) = interp(this_geom,this_order)%toFace(q_order)%mat(1:nsp,lf)
     !  !
     !  do l = 1,nr
     !    do m = 1,nq
     !      faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
     !                            dot( lmat(1:nsp) , xyflx_qp(1:nsp,m,l) ) * unrm(l)
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
   !   !
   !   !residual(m,k+n1-1) = dot( deriv(this_geom,this_order)%d(1)%mat(1:nsp,k) , &
   !   !                          ftdc(1:nsp,m) )
   !    residual(m,k+n1-1) = dot( deriv(this_geom,this_order)%d(1)%mat(1:nsp,k) , &
   !                              dfsp(1:nsp,m) )
   !    do l = 2,nr
   !   !  residual(m,k+n1-1) = residual(m,k+n1-1) + &
   !   !                       dot( deriv(this_geom,this_order)%d(l)%mat(1:nsp,k) , &
   !   !                            gtdc(1:nsp,m) )
   !      residual(m,k+n1-1) = residual(m,k+n1-1) + &
   !                           dot( deriv(this_geom,this_order)%d(l)%mat(1:nsp,k) , &
   !                                dgsp(1:nsp,m) )
   !    end do
   !  end do
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
  use order_mod,  only : maxpts,q_order,geom_solpts
  use geovar,  only : nfbnd,nface,nr
  use geovar,  only : bface,cells_on_face,pts_in_cell_ptr,xyz
  use geovar,  only : cell_geom,cell_order
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
  integer  :: nsp,nqp,kbl,kbr,kl,kr,iu
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
  real(wp), dimension(1:maxpts) :: lmat
  !
  real(wp), dimension(1:nr,1:maxpts) :: rnrm
  real(wp), dimension(1:nq,1:maxpts) :: fjl,fjr
  real(wp), dimension(1:nq,1:maxpts) :: lflx,rflx
  real(wp), dimension(1:nq,1:maxpts) :: luqp,ruqp,fncq
  !
  real(wp), dimension(1:maxpts,1:maxpts) :: pmat
  !
continue
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    face_geom = cells_on_face(7,nf)  ! geometry of face
    face_order = cells_on_face(8,nf) ! order of face
    !
    nsp = geom_solpts(face_geom,face_order) ! number of solution   flux points
    nqp = geom_solpts(face_geom,q_order)    ! number of quadrature flux points
    !
    left_cell = cells_on_face(1,nf) ! left  cell on face
    right_cell = cells_on_face(2,nf) ! right cell on face
    !
    left_geom = cell_geom(left_cell)
    left_order = cell_order(left_cell)
    !
    right_geom = cell_geom(right_cell)
    right_order = cell_order(right_cell)
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
    ! Interpolate the conserved variables and face normals at
    ! the solution flux points to the quadrature flux points
    !
    if (using_quadrature) then
      !
      do k = 1,nqp
        !
        lmat(1:nsp) = interp(face_geom,face_order)%toCell(q_order)%mat(1:nsp,k)
        !
        do m = 1,nq
          luqp(m,k) = dot( lmat(1:nsp) , faceusp(m,1:nsp,1,nf) )
          ruqp(m,k) = dot( lmat(1:nsp) , faceusp(m,1:nsp,2,nf) )
        end do
        !
        do l = 1,nr
          rnrm(l,k) = dot( lmat(1:nsp) , geofa(l,1:nsp,nf) )
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
      pmat(1:nqp,1:nsp) = projct(face_geom,q_order)%toCell(face_order)%mat(1:nqp,1:nsp)
      !
      if (dflux_opt == project_derivatives) then
        !
        lflx(1:nq,1:nqp) = fncq(1:nq,1:nqp) - lflx(1:nq,1:nqp)
        rflx(1:nq,1:nqp) = fncq(1:nq,1:nqp) - rflx(1:nq,1:nqp)
        !
        do k = 1,nsp
          !
          do m = 1,nq
            fjl(m,k) = dot( pmat(1:nqp,k) , lflx(m,1:nqp) )
            fjr(m,k) = dot( pmat(1:nqp,k) , rflx(m,1:nqp) )
          end do
          !
        end do
        !
      else
        !
        do k = 1,nsp
          !
          do m = 1,nq
            fnc(m) = dot( pmat(1:nqp,k) , fncq(m,1:nqp) )
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
  use order_mod,  only : maxpts
  use geovar,  only : nr,ncell,pts_in_cell_ptr,using_quadrature
  use geovar,  only : flxpts_of_cell,flxpts_of_cell_ptr
  use geovar,  only : cell_geom,cell_order
  use interpolation_mod, only : interp
  use derivatives_mod, only : deriv
  use metrics_mod,  only : metrics,geofa,jacobian
  use flowvar, only : usp,dusp,residual,faceflx,faceusp
  use ovar,    only : dflux_opt
  !
  !.. Local Scalars ..
  integer :: k,l,m,nc,n1,n2,np
  integer :: f1,f2,nf,sf,kf,lf
  integer :: this_geom,this_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:nq,1:nr) :: vflx
  real(wp), dimension(1:maxpts)  :: lmat
  !
  real(wp), dimension(1:maxpts,1:nq,1:nr) :: rsflx
  real(wp), dimension(1:maxpts,1:nq,1:nr) :: xyflx
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
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1 ! beginning index of sol-pts for this cell
    n2 = pts_in_cell_ptr(nc+1) ! ending    index of sol-pts for this cell
    np = n2-n1+1               ! number of sol-pts in this cell
    !
    this_geom = cell_geom(nc)   ! geometry type of this cell
    this_order = cell_order(nc) ! solution order of this cell
    !
    f1 = flxpts_of_cell_ptr(nc)+1 ! beginning index of flx-pts for this cell
    f2 = flxpts_of_cell_ptr(nc+1) ! ending    index of flx-pts for this cell
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
      vflx(1:nq,1:nr) = visc_flx_cv( usp(1:nq,l) , dusp(1:nr,1:nq,l) )
      !
      ! Subtract the viscous flux from the inviscid flux
      !
      xyflx(k,1:nq,1:nr) = xyflx(k,1:nq,1:nr) - vflx(1:nq,1:nr)
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
        !
        ! Initialize residual at this point by calculating the
        ! divergence of the flux in the first (i.e. x/xi/r) direction
        !
        residual(m,k+n1-1) = dot( deriv(this_geom,this_order)%d(1)%mat(1:np,k) , &
                                  rsflx(1:np,m,1) )
        !
        ! Loop through the remaining coordinate directions and add the
        ! divergence of the flux from each direction to the previous result
        !
        do l = 2,nr
          residual(m,k+n1-1) = residual(m,k+n1-1) + &
                             dot( deriv(this_geom,this_order)%d(l)%mat(1:np,k) , &
                                  rsflx(1:np,m,l) )
        end do
        !
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
      nf = flxpts_of_cell(1,k) ! index of the face that this flux point is on
      kf = flxpts_of_cell(2,k) ! index of flux point local to this face
      lf = flxpts_of_cell(3,k) ! index of flux point local to this cell
      sf = flxpts_of_cell(4,k) ! side of face cell nc is on
      !
      ! Get the unit normal for this face flux/solution point
      !
      unrm(:) = unit_vector( geofa(:,kf,nf) )
      !
      !######################################################################
      !######################################################################
      ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
      ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
      ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
      !######################################################################
      !######################################################################
      !
      lmat(1:np) = interp(this_geom,this_order)%toFace(this_order)%mat(1:np,lf)
      !
      do l = 1,nr
        do m = 1,nq
          faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
                                dot( lmat(1:np) , xyflx(1:np,m,l) ) * unrm(l)
        end do
      end do
      !
    end do
    !
  end do
  !
end subroutine derflux1_qdtr_with_visc
!
!###############################################################################
!
subroutine correct_interface_flux_qdtr_with_visc()
  !
  use order_mod,  only : maxpts,geom_solpts
  use geovar,  only : nr,nfbnd,nface
  use geovar,  only : bface,cells_on_face,pts_in_cell_ptr
  use geovar,  only : cell_geom,cell_order
  use ovar,    only : invs_flux_method
  use metrics_mod,  only : geofa,volj
  use correction_mod, only : correct
  use flowvar, only : usp,residual,faceflx,faceusp,facedusp
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nf,r1,r2,l1,l2
  integer  :: kbl,kbr,kl,kr,iu,nfp
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:nq) :: iflx,vflx
  real(wp), dimension(1:nr,1:maxpts) :: rnrm
  real(wp), dimension(1:nq,1:maxpts) :: fjl,fjr
  real(wp), dimension(1:nq,1:maxpts) :: lflx,rflx
  real(wp), dimension(1:nq,1:maxpts) :: luqp,ruqp,fncq
  real(wp), dimension(1:nr,1:nq,1:maxpts) :: lduqp,rduqp
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
    face_geom = cells_on_face(7,nf)  ! geometry of face
    face_order = cells_on_face(8,nf) ! order of face
    !
    nfp = geom_solpts(face_geom,face_order) ! number of flux points on this face
    !
    host_cell = cells_on_face(1,nf) ! left  cell on face
    !
    host_geom = cell_geom(host_cell)
    host_order = cell_order(host_cell)
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    !
    l1 = pts_in_cell_ptr(host_cell)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(host_cell+1) ! ending    index for sol-pts in left cell
    !
    lflx(1:nq,1:nfp) = faceflx(1:nq,1:nfp,1,nf)
    !
    luqp(1:nq,1:nfp) = faceusp(1:nq,1:nfp,1,nf)
    ruqp(1:nq,1:nfp) = faceusp(1:nq,1:nfp,2,nf)
    !
    lduqp(1:nr,1:nq,1:nfp) = facedusp(1:nr,1:nq,1:nfp,1,nf)
    rduqp(1:nr,1:nq,1:nfp) = facedusp(1:nr,1:nq,1:nfp,2,nf)
    !
    rnrm(1:nr,1:nfp) = geofa(1:nr,1:nfp,nf)
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
        vflx(:) = common_visc_flx_cv( unrm(:), luqp(:,k), lduqp(:,:,k), &
                                               ruqp(:,k), rduqp(:,:,k) )
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
          residual(m,n) = residual(m,n) + &
            dsl * fjl(m,k) * correct(host_geom,host_order)%mat(n-l1+1,kl)
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
    face_geom = cells_on_face(7,nf)  ! geometry of face
    face_order = cells_on_face(8,nf) ! order of face
    !
    nfp = geom_solpts(face_geom,face_order) ! number of flux points on this face
    !
    left_cell = cells_on_face(1,nf) ! left  cell on face
    right_cell = cells_on_face(2,nf) ! right cell on face
    !
    left_geom = cell_geom(left_cell)
    left_order = cell_order(left_cell)
    !
    right_geom = cell_geom(right_cell)
    right_order = cell_order(right_cell)
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
    lflx(1:nq,1:nfp) = faceflx(1:nq,1:nfp,1,nf)
    rflx(1:nq,1:nfp) = faceflx(1:nq,1:nfp,2,nf)
    !
    luqp(1:nq,1:nfp) = faceusp(1:nq,1:nfp,1,nf)
    ruqp(1:nq,1:nfp) = faceusp(1:nq,1:nfp,2,nf)
    !
    lduqp(1:nr,1:nq,1:nfp) = facedusp(1:nr,1:nq,1:nfp,1,nf)
    rduqp(1:nr,1:nq,1:nfp) = facedusp(1:nr,1:nq,1:nfp,2,nf)
    !
    rnrm(1:nr,1:nfp) = geofa(1:nr,1:nfp,nf)
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
      vflx(:) = common_visc_flx_cv( unrm(:), luqp(:,k), lduqp(:,:,k), &
                                             ruqp(:,k), rduqp(:,:,k) )
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
          residual(m,n) = residual(m,n) + &
            dsl * fjl(m,k) * correct(left_geom,left_order)%mat(n-l1+1,kl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell
      ! right_cell
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      do n = r1,r2
        do m = 1,nq
          residual(m,n) = residual(m,n) - &
            dsr * fjr(m,k) * correct(right_geom,right_order)%mat(n-r1+1,kr)
        end do
      end do
      !
    end do
    !
  end do interior_fluxes
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
  use order_mod,  only : maxpts
  use geovar,  only : nr,ncell
  use geovar,  only : pts_in_cell_ptr
  use geovar,  only : cell_geom,cell_order
  use metrics_mod,  only : jacobian
  use metrics_mod,  only : metrics
  use derivatives_mod, only : deriv
  use flowvar, only : usp,dusp
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nc,n1,n2,np
  integer :: this_geom,this_order
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
    n1 = pts_in_cell_ptr(nc)+1 ! beginning index for solpts in cell nc
    n2 = pts_in_cell_ptr(nc+1) ! ending    index for solpts in cell nc
    np = n2-n1+1 ! number of solution points in this cell
    !
    this_geom = cell_geom(nc)
    this_order = cell_order(nc)
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
          dudr(l) = dot( deriv(this_geom,this_order)%d(l)%mat(1:np,k) , tusp(1:np,m) )
        end do
        !
        ! Transform the derivatives back to x-y space
        !
        dudx(1) =  dudr(1)*metrics(2,2,n) - dudr(2)*metrics(2,1,n)
        dudx(2) = -dudr(1)*metrics(1,2,n) + dudr(2)*metrics(1,1,n)
        !
        dusp(1,m,n) = dudx(1)/jacobian(n)
        dusp(2,m,n) = dudx(2)/jacobian(n)
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
  use order_mod,  only : maxpts,geom_solpts
  use geovar,  only : nfbnd,nface,nr,n_solpts,ncell
  use geovar,  only : cells_on_face,pts_in_cell_ptr,bface
  use geovar,  only : cell_geom,cell_order
  use metrics_mod,  only : geofa,jacobian,volj
  use correction_mod, only : correct
  use flowvar, only : dusp,faceusp
  use ovar,    only : walls_are_exact
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nf,r1,r2,l1,l2,n1,n2,nc
  integer  :: kbl,kbr,kl,kr,nfp
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: cbias,hfpb,hfmb
  real(wp) :: correct_n,correct_m,ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)          :: unrm
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
    face_geom = cells_on_face(7,nf)  ! geometry of face
    face_order = cells_on_face(8,nf) ! order of face
    !
    nfp = geom_solpts(face_geom,face_order) ! number of flux points on this face
    !
    host_cell = cells_on_face(1,nf) ! left  cell on face
    !
    host_geom = cell_geom(host_cell)
    host_order = cell_order(host_cell)
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    !
    l1 = pts_in_cell_ptr(host_cell)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(host_cell+1) ! ending    index for sol-pts in left cell
    !
    ! For all the flux points on the current face, extract the conservative
    ! conservative variables for the left and right cells from faceusp
    !
    lufp(1:nq,1:nfp) = faceusp(1:nq,1:nfp,1,nf)
    rufp(1:nq,1:nfp) = faceusp(1:nq,1:nfp,2,nf)
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
        correct_n = dsl * correct(host_geom,host_order)%mat(n-l1+1,kl) / jacobian(n)
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
    face_geom = cells_on_face(7,nf)  ! geometry of face
    face_order = cells_on_face(8,nf) ! order of face
    !
    nfp = geom_solpts(face_geom,face_order) ! number of flux points on this face
    !
    left_cell = cells_on_face(1,nf) ! left  cell on face
    right_cell = cells_on_face(2,nf) ! right cell on face
    !
    left_geom = cell_geom(left_cell)
    left_order = cell_order(left_cell)
    !
    right_geom = cell_geom(right_cell)
    right_order = cell_order(right_cell)
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
    ! For all the flux points on the current face, extract the conservative
    ! conservative variables for the left and right cells from faceusp
    !
    lufp(1:nq,1:nfp) = faceusp(1:nq,1:nfp,1,nf)
    rufp(1:nq,1:nfp) = faceusp(1:nq,1:nfp,2,nf)
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
        correct_n = dsl * correct(left_geom,left_order)%mat(n-l1+1,kl) / jacobian(n)
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
        correct_n = dsr * correct(right_geom,right_order)%mat(n-r1+1,kr) / jacobian(n)
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
  use order_mod,  only : maxpts
  use geovar,  only : flxpts_of_cell_ptr,flxpts_of_cell
  use geovar,  only : nr,ncell,nfbnd,bface
  use geovar,  only : cells_on_face,pts_in_cell_ptr
  use geovar,  only : cell_geom,cell_order
  use metrics_mod,  only : geofa,jacobian,volj
  use correction_mod, only : correct
  use flowvar, only : dusp,facedusp,faceusp
  use ovar,    only : walls_are_exact
  use interpolation_mod, only : interp
  !
  use parallel_mod, only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,l,m
  integer :: nf,f1,f2,kf,lf,sf
  integer :: ncf,fp1,fp2,kfp,kbl,kl
  integer :: this_geom,this_order
  real(wp) :: dsl,cbias,hfpb,correct_n,correct_m
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:nq) :: lufp,rufp,ujump
  real(wp), dimension(1:maxpts) :: lmat
  real(wp), dimension(1:maxpts,1:nr,1:nq) :: face_dusp
  real(wp), dimension(1:nr,1:nq,1:maxpts) :: cell_dusp
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
    !
    this_geom = cell_geom(nc)
    this_order = cell_order(nc)
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
          correct_n = dsl * correct(this_geom,this_order)%mat(n,lf) / jacobian(n+n1-1)
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
        lmat(1:np) = interp(this_geom,this_order)%toFace(this_order)%mat(1:np,lf)
        !
        do m = 1,nq
          do l = 1,nr
            facedusp(l,m,kf,sf,nf) = dot( lmat(1:np) , face_dusp(1:np,l,m) )
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
include "Functions/viscosity.f90"
!
!###############################################################################
!
end module flux_mod
