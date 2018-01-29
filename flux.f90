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
  use order_mod,         only : maxSP,geom_solpts
  use geovar,            only : ncell,cell
  use interpolation_mod, only : interp
  use flowvar,           only : usp,faceusp
  !
  use parallel_mod,      only : exchange_faceusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,m
  integer :: nf,kf,lf,sf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxSP,1:nq) :: tusp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_flxpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
      face_geom = cell(nc)%face(n)%geom   ! geom  of current face
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      do kf = 1,geom_solpts(face_geom,face_order)
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        lf = cell(nc)%face(n)%get_fp_idx(kf)
        !
        faceusp(nf)%v(1:nq,kf,sf) = interp(this_geom,this_order)% &
                                              toFace(face_order)% &
                                    dot_mat( lf , tusp(1:np,1:nq) )
       !do m = 1,nq
       !  faceusp(nf)%v(m,kf,sf) = interp(this_geom,this_order)% &
       !                                     toFace(face_order)% &
       !                           dot( lf , tusp(1:np,m) )
       !end do
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine interp_cv_to_flxpts
!
!###############################################################################
!
subroutine interp_dusp_to_flxpts()
  !
  use order_mod,         only : maxSP,geom_solpts
  use geovar,            only : nr,ncell,cell
  use interpolation_mod, only : interp
  use flowvar,           only : dusp,facedusp
  !
  use parallel_mod,      only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,l,m
  integer :: nf,kf,lf,sf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tdusp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_dusp_to_flxpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
      face_geom = cell(nc)%face(n)%geom   ! geom  of current face
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      do kf = 1,geom_solpts(face_geom,face_order)
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        lf = cell(nc)%face(n)%get_fp_idx(kf)
        !
        facedusp(nf)%v(1:nr,1:nq,kf,sf) = interp(this_geom,this_order)% &
                                                    toFace(face_order)% &
                                          dot_mat2( lf , tdusp(1:np,1:nr,1:nq) )
       !do m = 1,nq
       !  do l = 1,nr
       !    facedusp(nf)%v(l,m,kf,sf) = interp(this_geom,this_order)% &
       !                                          toFace(face_order)% &
       !                                dot( lf , tdusp(1:np,l,m) )
       !  end do
       !end do
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
  call debug_timer(leaving_procedure,pname)
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
  use order_mod,         only : maxSP,geom_solpts
  use geovar,            only : nr,ncell,using_quadrature
  use geovar,            only : cell
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  use metrics_mod,       only : metrics_dt,geofa
  use flowvar,           only : usp,dusp,residual,faceflx,faceusp
  use ovar,              only : governing_equations
  use ovar,              only : output_dir,dump_fluxes
  !
  !.. Formal Arguments ..
  integer, intent(in) :: time_step
  integer, intent(in) :: rk_step
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nc,n1,n2,np
  integer :: nf,sf,kf,lf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: ierr,iosp
  !
  character(len=100) :: flux_fmt
  character(len=300) :: spf_fname
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:nq,1:nr) :: vflx
  real(wp), dimension(1:nq,1:maxSP) :: lcl_res
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: rsflx
  real(wp), dimension(1:maxSP,1:nq,1:nr) :: xyflx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "flux_gradient_at_solpts_visc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (dump_fluxes) then
    write (flux_fmt,4) nr,nr
    write (spf_fname,1) trim(output_dir),time_step,rk_step
    spf_fname = trim(adjustl(spf_fname))
    open (newunit=iosp,file=trim(spf_fname),status="replace",action="write", &
                       position="rewind",form="formatted", &
                       iostat=ierr,iomsg=error_message)
    call io_error(pname,trim(spf_fname),1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Initialize faceflx to zero before we start
  !
#ifdef SPECIAL_FOR_PGI
  call faceflx%zero_face_var
#else
  call faceflx%zero
#endif
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
      do l = 1,nr
        xyflx(k,:,l) = invs_flx_cv( usp(:,n) , l ) ! x-directional flux
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
          rsflx(k,m,l) = dot(xyflx(k,m,1:nr),metrics_dt(nc)%met(k,l,1:nr)) * &
                         metrics_dt(nc)%jac(k)
        end do
      end do
      !
    end do
    !
    if (dump_fluxes) then
      write (iosp,2) nc
      do k = 1,np
        write (iosp,3) k,k+n1-1
        do m = 1,nq
          write (iosp,flux_fmt) (xyflx(k,m,l),l=1,nr), (rsflx(k,m,l),l=1,nr)
        end do
      end do
    end if
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
    do m = 1,nq
      lcl_res(m,1:np) = deriv(this_geom,this_order)%div( rsflx(1:np,m,1:nr) )
    end do
    !
    ! Copy the local residual for the current cell into the residual array
    !
    residual(1:nq,n1:n2) = lcl_res(1:nq,1:np)
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
      face_geom = cell(nc)%face(n)%geom   ! geom  of current face
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      do kf = 1,geom_solpts(face_geom,face_order)
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        lf = cell(nc)%face(n)%get_fp_idx(kf)
        !
        ! Get the unit normal for this face flux/solution point
        !
        unrm(:) = unit_vector( geofa(nf)%v(:,kf) )
        !
        vflx(1:nq,1:nr) = interp(this_geom,this_order)% &
                                    toFace(face_order)% &
                          dot_mat2( lf , xyflx(1:np,1:nq,1:nr) )
       !do l = 1,nr
       !  do m = 1,nq
       !    vflx(m,l) = interp(this_geom,this_order)% &
       !                          toFace(face_order)% &
       !                dot( lf , xyflx(1:np,m,l) )
       !  end do
       !end do
        !
        do m = 1,nq
          faceflx(nf)%v(m,kf,sf) = dot( vflx(m,1:nr) , unrm(1:nr) )
        end do
        !
      end do
      !
    end do
    !
  end do
  !
  if (dump_fluxes) then
    close (iosp,iostat=ierr,iomsg=error_message)
    call io_error(pname,trim(spf_fname),2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a,"flux_files/sp_",i0,".",i0,"_fluxes.dat")
  2 format (100("#"),/,15x,"CELL ",i0,/,100("#"))
  3 format (100("-"),/,5x,"SP ",i0," (",i0,")",/,100("-"))
  4 format ("(",i0,"es15.6,5x,",i0,"es15.6)")
  !
end subroutine flux_gradient_at_solpts_visc
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
  use flowvar,           only : usp,residual,faceflx,faceusp,facedusp
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
  real(wp), dimension(1:nq,1:maxFP) :: fjl,fjr
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
  if (dump_fluxes) then
    write (fpf_fname,2) trim(output_dir),time_step,rk_step
    fpf_fname = trim(adjustl(fpf_fname))
    open (newunit=iofp,file=trim(fpf_fname),status="replace",action="write", &
                       position="rewind",form="formatted", &
                       iostat=ierr,iomsg=error_message)
    call io_error(pname,trim(fpf_fname),1,__LINE__,__FILE__,ierr,error_message)
  end if
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
    if (dump_fluxes) then
      write (iofp,3) nf,host_cell
      do k = 1,nfp
        write (iofp,4) k,face(nf)%left%get_fp_idx(face_geom,face_order,k)
        do m = 1,nq
          write (iofp,7) fncq(m,k), lflx(m,k)
        end do
      end do
    end if
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
          residual(m,n) = residual(m,n) + fjl(m,k) * cmat_val
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
  ! Loop over the interior faces and apply the correction from
  ! the fluxes at each face point to all the cell solution points
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
    if (dump_fluxes) then
      write (iofp,5) nf,left_cell,right_cell
      do k = 1,nfp
        write (iofp,6) k,face(nf)%left%get_fp_idx(face_geom,face_order,k), &
                         face(nf)%right%get_fp_idx(face_geom,face_order,k)
        do m = 1,nq
          write (iofp,7) fncq(m,k), lflx(m,k), rflx(m,k)
        end do
      end do
    end if
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
          residual(m,n) = residual(m,n) + fjl(m,k) * cmat_val
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
          residual(m,n) = residual(m,n) - fjr(m,k) * cmat_val
        end do
      end do
      !
    end do
    !
  end do interior_fluxes
  !
  if (dump_fluxes) then
    close (iofp,iostat=ierr,iomsg=error_message)
    call io_error(pname,trim(fpf_fname),2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a,i0)
  2 format (a,"flux_files/fp_",i0,".",i0,"_fluxes.dat")
  3 format (100("#"),/, &
            15x,"BOUNDARY FACE ",i0," :   HC=",i0,/, &
            100("#"))
  4 format (100("-"),/, &
            5x,"FP ",i0," :   HP=",i0,/, &
            100("-"))
  5 format (100("#"),/, &
            15x,"INTERIOR FACE ",i0," :   LC=",i0,", RC=",i0,/, &
            100("#"))
  6 format (100("-"),/, &
            5x,"FP ",i0," :   LP=",i0,", RP=",i0,/, &
            100("-"))
  7 format (es16.6,5x,es15.6,:,es15.6)
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
  use order_mod,       only : maxpts
  use geovar,          only : nr,ncell,cell
  use metrics_mod,     only : metrics_dt
  use derivatives_mod, only : deriv
  use flowvar,         only : usp,dusp
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nc,n1,n2,np
  integer :: this_geom,this_order
  !
  !.. Local Arrays ..
 !real(wp), dimension(1:nr) :: dudr
  real(wp), dimension(1:nr,1:nq) :: dudr
  real(wp), dimension(1:maxpts,1:nq) :: tusp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "discontinuous_gradient_at_solpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
      dudr(1:nr,1:nq) = deriv(this_geom,this_order)%grad( k , tusp(1:np,1:nq) )
      !
      do m = 1,nq
        do l = 1,nr
          dusp(l,m,n) = dot( dudr(1:nr,m) , metrics_dt(nc)%met(k,1:nr,l) )
        end do
      end do
      !
     !do m = 1,nq
     !  !
     !  ! Compute the derivatives in r-s space
     !  !
     !  do l = 1,nr
     !    dudr(l) = deriv(this_geom,this_order)% &
     !                                     d(l)% &
     !              dot( k , tusp(1:np,m) )
     !  end do
     !  !
     !  ! Transform the derivatives back to x-y space
     !  !
     !  do l = 1,nr
     !    dusp(l,m,n) = dot( dudr(1:nr) , metrics_dt(nc)%met(k,1:nr,l) )
     !  end do
     !  !
     !end do
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine correct_gradients_at_solpts
!
!###############################################################################
!
subroutine correct_gradients_br2()
  !
  use generic_types_mod, only : matrix
  use order_mod,         only : maxSP,maxFP,geom_solpts
  use geovar,            only : cell
  use geovar,            only : nr,ncell,nfbnd,bface
  use geovar,            only : face
  use metrics_mod,       only : geofa,metrics_dt
  use correction_mod,    only : correct
  use flowvar,           only : dusp,facedusp,faceusp
  use ovar,              only : walls_are_exact,LDG_beta
  use interpolation_mod, only : interp
  !
  use parallel_mod,      only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,k,l,m
  integer :: nf,kf,lf,sf,ncf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  real(wp) :: dsl,cbias,hfpb,correct_n,correct_m
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:nq) :: lufp,rufp,ujump
  real(wp), dimension(1:nr,1:nq,1:maxSP) :: cell_dusp
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
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
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
        lf = cell(nc)%face(ncf)%get_fp_idx(kf)
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
        cbias = zero
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
        lf = cell(nc)%face(ncf)%get_fp_idx(kf)
        !
        facedusp(nf)%v(1:nr,1:nq,kf,sf) = interp(this_geom,this_order)% &
                                                    toFace(face_order)% &
                                       dot_mat2( lf ,face_dusp(1:np,1:nr,1:nq) )
       !do m = 1,nq
       !  do l = 1,nr
       !    facedusp(nf)%v(l,m,kf,sf) = interp(this_geom,this_order)% &
       !                                          toFace(face_order)% &
       !                                dot( lf ,face_dusp(1:np,l,m) )
       !  end do
       !end do
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
  call debug_timer(leaving_procedure,pname)
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
#include "Functions/inviscid_fluxes.f90"
#include "Functions/viscous_fluxes.f90"
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
