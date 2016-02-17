module module_bc
  !
  ! Module to contain the procedures and variables related
  ! to the evaluating the boundary conditions for the flow.
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ntk,ntl
  !
  implicit none
  !
  private
  !
  !
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: get_faceusp_bc
  public :: get_facedusp_bc
  !
  !
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  !
  !
  ! ###################################
  ! ###  Module Allocatable Arrays  ###
  ! ###################################
  !
  !
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  !
contains
!
!###############################################################################
!
subroutine get_faceusp_bc
  !
  ! Routine to compute the conservative variables at the boundary flux points
  ! and store the results in the faceusp array. This should be called after
  ! the interp_cv_to_flxpts routine so the internal states at the boundary
  ! flux points should already be defined.  We use these internal states
  ! to compute the boundary states based on the face boundary condition.
  !
  !.. Use Statements ..
  use geovar,  only : bface,nfbnd
  use ovar,    only : output_bface_array
  use metrics_mod,  only : geofa
  use flowvar, only : faceusp,facexyz
  !
  use parallel_mod, only : wait_faceusp
  !
  !.. Local Scalars ..
  integer  :: nf,pf
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_faceusp_bc"
  !
continue
  !
  if (output_bface_array) call output_bface
  !
  ! Loop over the boundary faces and compute the conservative variables
  ! at the boundary flux points based on the face boundary condition.
  !
  do nf = 1,nfbnd
    !
    ! Compute the solution at the boundary flux points
    ! based on the specified boundary condition
    !
    select case (bface(1,nf))
      !
      case (bc_freestream)
        !
        ! Free stream boundary condition
        !
        call get_bc_freestream_solution(faceusp(nf)%v(:,:,2))
        !
      case (bc_characteristic)
        !
        ! Characteristic boundary condition
        !
        call get_bc_characteristic_solution(faceusp(nf)%v(:,:,1), &
                                            faceusp(nf)%v(:,:,2), &
                                            geofa(nf)%v(:,:))
        !
      case (bc_generic_inflow)
        !
        ! Generic inflow boundary condition
        !
        call get_bc_generic_inflow_solution(faceusp(nf)%v(:,:,1), &
                                            faceusp(nf)%v(:,:,2), &
                                            geofa(nf)%v(:,:))
        !
      case (bc_generic_outflow)
        !
        ! Generic outflow boundary condition
        !
        call get_bc_generic_outflow_solution(faceusp(nf)%v(:,:,1), &
                                             faceusp(nf)%v(:,:,2), &
                                             geofa(nf)%v(:,:))
        !
      case (bc_generic_freeflow)
        !
        ! Generic freeflow boundary condition
        !
        call get_bc_generic_freeflow_solution(faceusp(nf)%v(:,:,1), &
                                              faceusp(nf)%v(:,:,2), &
                                              geofa(nf)%v(:,:))
        !
      case (bc_sub_inflow)
        !
        ! Subsonic inflow boundary condition
        !
        call get_bc_sub_inflow_solution(faceusp(nf)%v(:,:,1), &
                                        faceusp(nf)%v(:,:,2))
        !
      case (bc_sub_outflow)
        !
        ! Subsonic outflow boundary condition
        !
        call get_bc_sub_outflow_solution(faceusp(nf)%v(:,:,1), &
                                         faceusp(nf)%v(:,:,2))
        !
      case (bc_sup_inflow)
        !
        ! Supersonic inflow boundary condition
        !
        call get_bc_sup_inflow_solution(faceusp(nf)%v(:,:,2))
        !
      case (bc_sup_outflow)
        !
        ! Supersonic outflow boundary condition
        !
        call get_bc_sup_outflow_solution(faceusp(nf)%v(:,:,1), &
                                         faceusp(nf)%v(:,:,2))
        !
      case (bc_slip_wall,bc_symmetry)
        !
        ! Slip wall / symmetry boundary condition
        !
        call get_bc_slip_wall_solution(faceusp(nf)%v(:,:,1), &
                                       faceusp(nf)%v(:,:,2), &
                                       geofa(nf)%v(:,:))
        !
      case (bc_adiabatic_wall)
        !
        ! No-slip adiabatic wall
        !
        call get_bc_adiabatic_wall_solution(faceusp(nf)%v(:,:,1), &
                                            faceusp(nf)%v(:,:,2))
        !
      case (bc_isothermal_wall)
        !
        ! No-slip isothermal wall
        !
        call get_bc_isothermal_wall_solution(faceusp(nf)%v(:,:,1), &
                                             faceusp(nf)%v(:,:,2))
        !
      case (bc_periodic)
        !
        ! Periodic boundary condition
        !
        pf = bface(5,nf) ! index of connecting periodic boundary face
        !
        call get_bc_periodic_solution(faceusp(pf)%v(:,:,1), &
                                      faceusp(nf)%v(:,:,2))
        !
      case (bc_mms_dirichlet)
        !
        ! Dirichlet boundary condition for MMS
        !
        call get_bc_mms_dirichlet_solution(faceusp(nf)%v(:,:,1), &
                                           faceusp(nf)%v(:,:,2), &
                                           facexyz(nf)%v(:,:), &
                                           geofa(nf)%v(:,:))
        !
      case (bc_cpu_bnd)
        !
        ! We dont need to do anything here for communication boundaries
        !
      case default
        !
        ! This shouldn't be here and is an unknown boundary condition
        !
        write (iout,1) bface(1,nf)
        call stop_recon(abort,pname)
        !
    end select
    !
  end do
  !
  ! Now that we have finished getting the boundary conditions local to
  ! this cell, make sure that the communication of faceusp between
  ! processors has completed before continuing onwards.
  !
  if (ncpu > 1) call wait_faceusp(faceusp)
  !
  call exchange_bc_conflicts_usp(faceusp)
  !
  ! Format Statements
  !
  1 format (/," ********************* ERROR! **********************",/, &
              "  UNKNOWN BOUNDARY CONDITION IN SUB GET_FACEUSP_BC!",/, &
              "  BOUNDARY CONDITION VALUE IS: ",i0,/, &
              " ********************* ERROR! **********************",/)
  !
end subroutine get_faceusp_bc
!
!###############################################################################
!
subroutine get_facedusp_bc
  !
  ! Routine to compute the derivatives of the conservative variables at the
  ! boundary flux points and store the results in the facedusp array. This
  ! should be called after the interp_dusp_to_flxpts routine so the internal
  ! gradients at the boundary flux points should already be defined.  We use
  ! these internal gradients to compute the boundary gradients based on the
  ! face boundary condition.
  !
  !.. Use Statements ..
  use geovar,  only : bface,nfbnd
  use flowvar, only : faceusp,facedusp
  !
  use parallel_mod, only : wait_facedusp
  !
  !.. Local Scalars ..
  integer  :: nf,pf
  !
  !.. Local Parameters ..
  integer, parameter :: local_bc_ignore(1:3) = [ bc_cpu_bnd, &
                                                 bc_unknown, &
                                                 not_a_bc ]
  character(len=*), parameter :: pname = "get_facedusp_bc"
  !
continue
  !
  ! Loop over the boundary faces and compute the conservative variables
  ! at the boundary flux points based on the face boundary condition.
  !
  do nf = 1,nfbnd
    !
    ! Compute the solution at the boundary flux points
    ! based on the specified boundary condition
    !
    if (any(bface(1,nf) == bc_walls)) then
      !
      ! Wall boundary condition
      !
      if (bface(1,nf) == bc_adiabatic_wall) then
        !
        ! No-slip adiabatic wall boundary conditions
        !
        call get_bc_adiabatic_wall_gradient(faceusp(nf)%v(:,:,1), &
                                            facedusp(nf)%v(:,:,:,1), &
                                            facedusp(nf)%v(:,:,:,2))
        !
      else
        !
        ! Non-adiabatic wall boundary condition
        !
        call get_bc_wall_gradient(faceusp(nf)%v(:,:,1), &
                                  facedusp(nf)%v(:,:,:,1), &
                                  facedusp(nf)%v(:,:,:,2))
        !
      end if
      !
    else if (bface(1,nf) == bc_periodic) then
      !
      ! Periodic boundary condition
      !
      pf = bface(5,nf) ! index of connecting periodic boundary face
      !
      ! We just want to copy the gradients from face pf to face nf which is
      ! exactly what the subroutine get_bc_flowthrough_gradient does, so
      ! just use that instead of making a completely separate subroutine.
      !
      call get_bc_flowthrough_gradient(facedusp(pf)%v(:,:,:,1), &
                                       facedusp(nf)%v(:,:,:,2))
      !
    else if (all(bface(1,nf) /= local_bc_ignore)) then
      !
      ! General flow-through boundary conditions
      !
      call get_bc_flowthrough_gradient(facedusp(nf)%v(:,:,:,1), &
                                       facedusp(nf)%v(:,:,:,2))
      !
   !else
   !  !
   !  ! This shouldn't be here and is an unknown boundary condition
   !  !
   !  write (iout,1) bface(1,nf)
   !  call stop_recon(abort,pname)
   !  !
    end if
    !
  end do
  !
  ! Now that we have finished getting the boundary conditions local to
  ! this cell, make sure that the communication of facedusp between
  ! processors has completed before continuing onwards.
  !
  if (ncpu > 1) call wait_facedusp(facedusp)
  !
  call exchange_bc_conflicts_dusp(facedusp)
  !
  ! Format Statements
  !
  1 format (/," ********************** ERROR! **********************",/, &
              "  UNKNOWN BOUNDARY CONDITION IN SUB GET_FACEDUSP_BC!",/, &
              "  BOUNDARY CONDITION VALUE IS: ",i0,/, &
              " ********************** ERROR! **********************",/)
  !
end subroutine get_facedusp_bc
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! Functions that compute the boundary condition
! solution flow variables in conservative form
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine get_bc_freestream_solution(gufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Free stream boundary condition
  ! #####################################################
  ! #####################################################
  ! 
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar, only : cv_in
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: k
  !
continue
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Set all variables to free stream
    !
    gufp(:,k) = cv_in(:)
    !
  end do
  !
end subroutine get_bc_freestream_solution
!
!###############################################################################
!
subroutine get_bc_characteristic_solution_old(hufp,gufp,fnrm)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FNRM : Face outward normal vector at the flux points
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar,        only : pv_in,cv_in,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: k
  real(wp) :: rp,rm,utc,unc,ac,sc,rho,pres
  real(wp) :: mach_int,gm1,gm1_inv
  real(wp) :: aspd_int,vt_int,vn_int
  real(wp) :: aspd_ext,vt_ext,vn_ext
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(fnrm,dim=1)) :: rtan
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  ! Precompute gamma-1 and its inverse
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Inward pointing unit normal at this flux point
    !
    rnrm(:) = -unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound in the host cell
    !
    aspd_int = sqrt( max( gam*hvfp(nq)/hvfp(1) , rndoff ) )
    !
    ! Compute the interior normal velocity at this flux point
    !
    vn_int = dot( rnrm , hvfp(nmb:nme) )
    !
    ! Apply boundary conditions based on the
    ! normal Mach number to the boundary face
    !
    mach_int = vn_int/aspd_int
    !
    if (mach_int >= one) then ! Supersonic inflow
      !
      ! Set ghost cell flux point to free stream conditions
      !
      gufp(:,k) = cv_in(:)
      !
    else if (mach_int <= -one) then ! Supersonic outflow
      !
      ! Extrapolate interior to ghost cell flux point
      !
      gufp(:,k) = hufp(:,k)
      !
    else ! Subsonic inflow/outflow using Riemann invariants
      !
      ! Unit tangent at this flux point computed from unit normal
      ! NOTE: Usually t = (tx,ty) = (ny,-nx) but we computed the inward
      !       pointing unit normal, so here (tx,ty) = (-ny,nx).
      !       
      rtan(:) = [-rnrm(2),rnrm(1)]
      !
      ! Compute the exterior speed of sound
      !
      aspd_ext = sqrt( max( gam*pv_in(nq)/pv_in(1) , rndoff ) )
      !
      ! Compute the exterior normal velocity at this flux point
      !
      vn_ext = dot( rnrm , pv_in(nmb:nme) )
      !
     !if (mach_int > zero) then ! this is an inflow boundary
     !  !
     !  rp = vn_ext + two*aspd_ext*gm1_inv
     !  rm = vn_int - two*aspd_int*gm1_inv
     !  unc = half*(rp+rm)
     !  utc = dot( rtan , pv_in(nmb:nme) )
     !  ac = one4*gm1*(rp-rm)
     !  sc = pv_in(nq)/(pv_in(1)**gam)
     !  !
     !else ! this is an outflow boundary
     !  !
     !  rp = -abs(vn_ext) + two*aspd_ext*gm1_inv
     !  rm =      vn_int  - two*aspd_int*gm1_inv
     !  unc = half*(rp+rm)
     !  utc = dot( rtan , hvfp(nmb:nme) )
     !  ac = one4*gm1*(rp-rm)
     !  sc = hvfp(nq)/(hvfp(1)**gam)
     !  !
     !end if
      !
      ! Compute the variables dependent on the interior flow direction
      !
      if (mach_int > zero) then
        rp   = vn_ext
        utc  = dot( rtan , pv_in(nmb:nme) )
        rho  = pv_in(1)
        pres = pv_in(nq)
      else
        rp   = -abs(vn_ext)
        utc  = dot( rtan , hvfp(nmb:nme) )
        rho  = hvfp(1)
        pres = hvfp(nq)
      end if
      !
      ! Compute the Riemann invariants
      !
      rp  = rp     + two*aspd_ext*gm1_inv
      rm  = vn_int - two*aspd_int*gm1_inv
      unc = half*(rp+rm)
      ac  = one4*gm1*(rp-rm)
      sc  = pres/(rho**gam)
      !
      ! Compute the density and pressure from the Riemann invariants
      !
      rho = ( (ac*ac) / (gam*sc) )**(gm1_inv)
      pres = sc * (rho**gam)
      !
      ! Compute the conservative variables at the boundary flux points
      !
      gufp(1,k) = rho
      gufp(nmb:nme,k) = rho * (unc*rnrm(:) + utc*rtan(:))
      gufp(nq,k) = pres*gm1_inv + half * selfdot( gufp(nmb:nme,k) ) / rho
      !
    end if
    !
  end do
  !
end subroutine get_bc_characteristic_solution_old
!
!###############################################################################
!
subroutine get_bc_characteristic_solution(hufp,gufp,fnrm)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FNRM : Face outward normal vector at the flux points
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar,        only : pv_in,cv_in,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: k
  real(wp) :: rp,rm,utc,unc,ac,sc,rho,pres
  real(wp) :: mach_int,mach_ext,gm1,gm1_inv
  real(wp) :: aspd_int,vt_int,vn_int
  real(wp) :: aspd_ext,vt_ext,vn_ext
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(fnrm,dim=1)) :: rtan
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  ! Precompute gamma-1 and its inverse
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Inward pointing unit normal at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound in the host cell
    !
    ! interior speed of sound 
    aspd_int = sqrt( max( gam*hvfp(nee)/hvfp(nec) , rndoff ) )
    ! exterior speed of sound 
    aspd_ext = sqrt( max( gam*pv_in(nee)/pv_in(nec) , rndoff ) )
    !
    ! Compute the interior normal velocity at this flux point
    !
    vn_int = dot( rnrm , hvfp(nmb:nme) )  ! interior normal velocity
    vn_ext = dot( rnrm , pv_in(nmb:nme) ) ! exterior normal velocity
    !
    ! Apply boundary conditions based on the
    ! normal Mach number to the boundary face
    !
    mach_int = vn_int/aspd_int ! interior normal Mach number
    mach_ext = vn_ext/aspd_ext ! exterior normal Mach number
    !
    if (abs(vn_ext) >= aspd_ext .and. vn_int >= zero) then
      rp = vn_int - two*aspd_int*gm1_inv
    else
      rp = vn_ext - two*aspd_ext*gm1_inv
    end if
    !
    if (abs(vn_ext) >= aspd_ext .and. vn_int < zero) then
      rm = vn_ext + two*aspd_ext*gm1_inv
    else
      rm = vn_int + two*aspd_int*gm1_inv
    end if
    !
    unc = half*(rp+rm)
    ac  = one4*gm1*(rm-rp)
    sc  = pv_in(nee) / (pv_in(nec)**gam)
    !
    if (vn_int < zero) then
      rho = ( ac * ac / (gam*sc) )**gm1_inv
    else
      rho = hvfp(nec) * (hvfp(nec)*ac*ac/(gam*hvfp(nee)))**gm1_inv
    end if
    !
    pres = rho * ac * ac / gam
    !
    ! Compute the conservative variables at the boundary flux points
    !
    gufp(nec,k) = rho
    if (vn_int >= zero) then
      gufp(nmb:nme,k) = rho * (hvfp(nmb:nme) + (unc - vn_int)*rnrm(:))
    else
      gufp(nmb:nme,k) = rho * (pv_in(nmb:nme) + (unc - vn_ext)*rnrm(:))
    end if
    gufp(nee,k) = pres*gm1_inv + half * selfdot( gufp(nmb:nme,k) ) / rho
    !
  end do
  !
end subroutine get_bc_characteristic_solution
!
!###############################################################################
!
subroutine get_bc_generic_inflow_solution(hufp,gufp,fnrm)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic inflow boundary condition
  ! ##    - Use subsonic or supersonic inflow boundary
  ! ##      conditions based on the Mach number of the
  ! ##      normal inflow velocity
  ! #####################################################
  ! #####################################################
  ! 
  ! FNRM : Face outward normal vector at the flux points
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar,        only : pv_in,cv_in,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nr,k
  real(wp) :: gm1_inv,ske_inf,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  nr = size(fnrm,dim=1)
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  ske_inf = half * pv_in(1) * selfdot( pv_in(nmb:nme) )
  !
  gm1_inv = one / (gam - one)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Unit normal and tangent at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*hvfp(nq)/hvfp(1) , rndoff ) )
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    mach_int = dot( hvfp(nmb:nme) , rnrm(1:nr) ) / aspd_int
    !
    if (mach_int > -one .and. mach_int <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables but pressure to free stream
      !
      gufp(1:nq-1,k) = cv_in(1:nq-1)
      !
      ! Compute the total energy using interior pressure
      ! and free stream specific kinetic energy
      !
      gufp(nq,k) = hvfp(nq) * gm1_inv + ske_inf
      !
    else if (mach_int <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      gufp(:,k) = cv_in(:)
      !
    else
      !
      ! The normal velocity at this flux point seems to
      ! indicate this should be an outflow boundary.
      !
      write (iout,*) "Error with bc_generic_inflow boundary condition!"
      !
    end if
    !
  end do
  !
end subroutine get_bc_generic_inflow_solution
!
!###############################################################################
!
subroutine get_bc_generic_outflow_solution(hufp,gufp,fnrm)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic outflow boundary condition
  ! ##    - Use subsonic or supersonic outflow boundary
  ! ##      conditions based on the Mach number of the
  ! ##      normal outflow velocity
  ! #####################################################
  ! #####################################################
  ! 
  ! FNRM : Face outward normal vector at the flux points
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar,        only : pv_in,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nr,k
  real(wp) :: gm1_inv,ske_int,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  nr = size(fnrm,dim=1)
  !
  ! Precompute the inverse of gamma-1
  !
  gm1_inv = one / (gam - one)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Unit normal and tangent at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Now convert the current host cell flux point
    ! conservative variables to primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*hvfp(nq)/hvfp(1) , rndoff ) )
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    mach_int = dot( hvfp(nmb:nme) , rnrm(1:nr) ) / aspd_int
    !
    if (mach_int < one .and. mach_int >= zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      !
      gufp(1:nq-1,k) = hufp(1:nq-1,k)
      !
      ! Compute the interior specific kinetic energy
      !
      ske_int = half * selfdot( hufp(nmb:nme,k) ) / hufp(1,k)
      !
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      gufp(nq,k) = pv_in(nq) * gm1_inv + ske_int
      !
    else if (mach_int >= one) then
      !
      ! Supersonic outflow
      !
      ! Extrapolate all interior variables to the boundary
      !
      gufp(1:nq,k) = hufp(1:nq,k)
      !
    else
      !
      ! The normal velocity at this flux point seems to
      ! indicate this should be an inflow boundary.
      !
      write (iout,*) "Error with bc_generic_outflow boundary condition!"
      !
    end if
    !
  end do
  !
end subroutine get_bc_generic_outflow_solution
!
!###############################################################################
!
subroutine get_bc_generic_freeflow_solution(hufp,gufp,fnrm)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic freeflow boundary condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #####################################################
  ! #####################################################
  ! 
  ! FNRM : Face outward normal vector at the flux points
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar,        only : pv_in,cv_in,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nr,k
  real(wp) :: gm1_inv,ske_inf,ske_int,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  nr = size(fnrm,dim=1)
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  ske_inf = half * pv_in(1) * selfdot( pv_in(nmb:nme) )
  !
  gm1_inv = one / (gam - one)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Unit normal and tangent at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*hvfp(nq)/hvfp(1) , rndoff ) )
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    mach_int = dot( hvfp(nmb:nme) , rnrm(1:nr) ) / aspd_int
    !
    if (mach_int >= one) then
      !
      ! Supersonic outflow
      !
      ! Extrapolate all interior variables to the boundary
      !
      gufp(1:nq,k) = hufp(1:nq,k)
      !
    else if (mach_int <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      gufp(:,k) = cv_in(:)
      !
    else if (mach_int > -one .and. mach_int <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables but pressure to free stream
      !
      gufp(1:nq-1,k) = cv_in(1:nq-1)
      !
      ! Compute the total energy using interior pressure
      ! and free stream specific kinetic energy
      !
      gufp(nq,k) = hvfp(nq) * gm1_inv + ske_inf
      !
    else if (mach_int < one .and. mach_int > zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      !
      gufp(1:nq-1,k) = hufp(1:nq-1,k)
      !
      ! Compute the interior specific kinetic energy
      !
      ske_int = half * selfdot( hufp(nmb:nme,k) ) / hufp(1,k)
      !
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      gufp(nq,k) = pv_in(nq) * gm1_inv + ske_int
      !
    end if
    !
  end do
  !
end subroutine get_bc_generic_freeflow_solution
!
!###############################################################################
!
subroutine get_bc_sub_inflow_solution(hufp,gufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Subsonic inflow boundary condition
  ! ##    - Fix density to free stream
  ! ##    - Fix velocities to free stream
  ! ##    - Extrapolate pressure
  ! #####################################################
  ! #####################################################
  !
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar, only : pv_in,cv_in,gam
  use ovar, only : ptotref_nd,ttotref_nd
  use ovar, only : rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer  :: k
  real(wp) :: pres_int,ske_inf,gm1_inv
  real(wp) :: gm1,rt,rhot,a,b,c,velx2
  !
continue
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  rt = rgasref_nd*ttotref_nd
  !
  rhot = ptotref_nd/rt
  a = rt*gm1_inv
  b = two*gam*a
  c = half*gm1/gam
  !
 !do k = 1,size(gufp,dim=2)
 !  !
 !  velx2 = hufp(nmx,k)**2
 !  !
 !  gufp(nec,k) = rhot * ( (one + velx2/(b-velx2))**-gm1_inv )
 !  gufp(nmb:nme,k) = zero
 !  gufp(nmx,k) = gufp(nec,k)*hufp(nmx,k)
 !  gufp(nee,k) = gufp(nec,k) * (a + c*velx2)
 !  !
 !end do
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  ske_inf = half * pv_in(1) * selfdot( pv_in(nmb:nme) )
  !
  gm1_inv = one / (gam - one)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Set all variables but pressure to free stream
    !
    gufp(1:nq-1,k) = cv_in(1:nq-1)
    !
    ! Get the interior pressure
    !
    pres_int = pressure_cv_sp( hufp(:,k) )
    !
    ! Compute the total energy using interior pressure
    ! and free stream specific kinetic energy
    !
    gufp(nq,k) = pres_int * gm1_inv + ske_inf
    !
  end do
  !
end subroutine get_bc_sub_inflow_solution
!
!###############################################################################
!
subroutine get_bc_sub_outflow_solution(hufp,gufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Subsonic outflow boundary condition
  ! ##    - Extrapolate density
  ! ##    - Extrapolate velocities
  ! ##    - Fix pressure to free stream
  ! #####################################################
  ! #####################################################
  ! 
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar, only : pv_in,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: k
  real(wp) :: gm1_inv,ske_int
  !
continue
  !
  ! Precompute the inverse of gamma-1
  !
  gm1_inv = one / (gam - one)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Extrapolate interior density and velocities to boundary
    !
    gufp(1:nq-1,k) = hufp(1:nq-1,k)
    !
    ! Compute the interior specific kinetic energy
    !
    ske_int = half * selfdot( hufp(nmb:nme,k) ) / hufp(1,k)
    !
    ! Compute the total energy using free stream pressure
    ! and interior specific kinetic energy
    !
    gufp(nq,k) = pv_in(nq) * gm1_inv + ske_int
    !
  end do
  !
end subroutine get_bc_sub_outflow_solution
!
!###############################################################################
!
subroutine get_bc_sup_inflow_solution(gufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic inflow boundary condition
  ! ##    - Fix all variables to free stream
  ! #####################################################
  ! #####################################################
  ! 
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar, only : cv_in
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: k
  !
continue
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Set all variables to free stream
    !
    gufp(:,k) = cv_in(:)
    !
  end do
  !
end subroutine get_bc_sup_inflow_solution
!
!###############################################################################
!
subroutine get_bc_sup_outflow_solution(hufp,gufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic outflow boundary condition
  ! ##    - Extrapolate all variables
  ! #####################################################
  ! #####################################################
  ! 
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: k
  !
continue
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Extrapolate all interior variables to the boundary
    !
    gufp(1:nq,k) = hufp(1:nq,k)
    !
  end do
  !
end subroutine get_bc_sup_outflow_solution
!
!###############################################################################
!
subroutine get_bc_slip_wall_solution(hufp,gufp,fnrm)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FNRM : Face outward normal vector at the flux points
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar,        only : walls_are_exact
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: k
  real(wp) :: vn,wall_solution_location
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  wall_solution_location = merge(zero,one,walls_are_exact)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Unit normal at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Get the interior primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the normal velocity at this boundary flux point
    !
    vn = dot( hvfp(nmb:nme) , rnrm )
    !
    ! For velocities, subtract (wall_solution_location+one) * (the wall
    ! normal velocity) from the interior velocity to get the wall
    ! boundary velocity.
    !
    hvfp(nmb:nme) = hvfp(nmb:nme) - (wall_solution_location+one)*vn*rnrm(:)
    !
    ! Convert the primitive variables to conservative.
    ! This ensures that the interior density and pressure are
    ! extrapolated to the boundary and the total energy at the
    ! boundary is consistent with these two states and the
    ! computed boundary velocity.
    !
    gufp(:,k) = vsp2u_sp( hvfp(:) )
    !
  end do
  !
end subroutine get_bc_slip_wall_solution
!
!###############################################################################
!
subroutine get_bc_adiabatic_wall_solution(hufp,gufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar,        only : walls_are_exact
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: k
  real(wp) :: wall_solution_location
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  wall_solution_location = merge(zero,one,walls_are_exact)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Get the interior primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Apply the no-slip condition to the boundary velocities
    !
    hvfp(nmb:nme) = -hvfp(nmb:nme) * wall_solution_location
    !
    ! Convert the primitive variables to conservative.
    ! This ensures that the interior density and pressure are
    ! extrapolated to the boundary and the total energy at the
    ! boundary is consistent with these two states and the
    ! computed boundary velocity.
    !
    gufp(:,k) = vsp2u_sp( hvfp(:) )
    !
  end do
  !
end subroutine get_bc_adiabatic_wall_solution
!
!###############################################################################
!
subroutine get_bc_isothermal_wall_solution(hufp,gufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar,        only : iso_wall_temp,tref,rgasref_nd
  use ovar,        only : walls_are_exact
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: k
  real(wp) :: iso_wall_temp_nd,wall_solution_location
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  iso_wall_temp_nd = iso_wall_temp / tref
  !
  wall_solution_location = merge(zero,one,walls_are_exact)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Get the interior primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Apply the no-slip condition to the boundary velocities
    !
    hvfp(nmb:nme) = -hvfp(nmb:nme) * wall_solution_location
    !
    ! Keeping the wall pressure constant, compute the wall
    ! density using the isothermal wall temperature.
    !
    hvfp(1) = hvfp(nq) / rgasref_nd / iso_wall_temp_nd
    !
    ! Convert the primitive variables to conservative.
    ! This ensures that the total energy at the boundary is
    ! consistent with the extrapolated internal pressure, the
    ! isothermal wall temperature, and the computed boundary
    ! velocity.
    !
    gufp(:,k) = vsp2u_sp( hvfp(:) )
    !
  end do
  !
end subroutine get_bc_isothermal_wall_solution
!
!###############################################################################
!
subroutine get_bc_periodic_solution(pufp,gufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Periodic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! PUFP : Conservative variables at the periodic cell flux points
  ! GUFP : Conservative variables at the ghost    cell flux points
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(inout) :: pufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: k
  !
continue
  !
  do k = 1,size(gufp,dim=2)
    !
    ! The conservative variables at this boundary flux point correspond
    ! to the conservative variables at the interior boundary flux points
    ! within the matching periodic cell. Therefore just copy these
    ! values over.
    ! NOTE: We may need to reverse the copy order, i.e., 
    !          gufp(:,k) = pufp(:,size(gufp,dim=2)+1-k)
    !
    gufp(:,k) = pufp(:,k)
    !
  end do
  !
end subroutine get_bc_periodic_solution
!
!###############################################################################
!
subroutine get_bc_mms_dirichlet_solution(hufp,gufp,fxyz,fnrm)
  !
  ! #########################################################
  ! #########################################################
  ! ##  MMS Dirichlet (generic freeflow) Boundary Condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #########################################################
  ! #########################################################
  !
  ! GUFP : Conservative variables at the ghost cell flux points
  ! HUFP : Conservative variables at the host  cell flux points
  ! FXYZ : Coordinates of the flux points
  ! FNRM : Face outward normal vector at the flux points
  !
  !.. Use Statements ..
  use ovar, only : gam
  use mms_mod, only : mms_solution_at_xyz
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:),    intent(in) :: fxyz
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nr,k
  real(wp) :: gm1_inv,ske_inf,ske_int,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  real(wp), dimension(1:size(gufp,dim=1)) :: mms_sol
  !
continue
  !
  nr = size(fnrm,dim=1)
  !
  gm1_inv = one / (gam - one)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Get the 'freestream' conditions from the MMS solution
    !
    mms_sol(1:nq) = mms_solution_at_xyz( fxyz(:,k) , &
                                         return_nondim=true , &
                                         return_cv=true )
    !
    ! Unit normal and tangent at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*hvfp(nq)/hvfp(1) , rndoff ) )
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    mach_int = dot( hvfp(nmb:nme) , rnrm(1:nr) ) / aspd_int
    !
    if (mach_int >= one) then
      !
      ! Supersonic outflow
      !
      ! Extrapolate all interior variables to the boundary
      !
      gufp(1:nq,k) = hufp(1:nq,k)
      !
    else if (mach_int <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      gufp(:,k) = mms_sol(1:nq)
      !
    else if (mach_int > -one .and. mach_int <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables but pressure to free stream
      !
      gufp(1:nq-1,k) = mms_sol(1:nq-1)
      !
      ! Compute the total energy using interior pressure
      ! and free stream specific kinetic energy
      !
      ske_inf = half * selfdot( mms_sol(nmb:nme) ) / mms_sol(nec)
      !
      gufp(nq,k) = hvfp(nq) * gm1_inv + ske_inf
      !
    else if (mach_int < one .and. mach_int > zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      !
      gufp(1:nq-1,k) = hufp(1:nq-1,k)
      !
      ! Compute the interior specific kinetic energy
      !
      ske_int = half * selfdot( hufp(nmb:nme,k) ) / hufp(1,k)
      !
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      gufp(nq,k) = pressure_cv_sp( mms_sol(1:nq) ) * gm1_inv + ske_int
      !
    end if
    !
  end do
  !
end subroutine get_bc_mms_dirichlet_solution
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! Functions that compute the gradients of the
! conservative variables at the boundary flux points.
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine get_bc_flowthrough_gradient(hdufp,gdufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  General flow-through boundary conditions
  ! #####################################################
  ! #####################################################
  ! 
  ! HDUFP : Gradients of the conservative variables at the host  cell flux pts
  ! GDUFP : Gradients of the conservative variables at the ghost cell flux pts
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:,:),    intent(in) :: hdufp
  real(wp), dimension(:,:,:), intent(inout) :: gdufp
  !
  !.. Local Scalars ..
  integer :: i,j,k
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  do k = 1,size(gdufp,dim=3)
    do j = 1,size(gdufp,dim=2)
      do i = 1,size(gdufp,dim=1)
        gdufp(i,j,k) = hdufp(i,j,k)
      end do
    end do
  end do
  !
end subroutine get_bc_flowthrough_gradient
!
!###############################################################################
!
subroutine get_bc_adiabatic_wall_gradient(hufp,hdufp,gdufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! HUFP  : Conservative variables at the host cell flux points
  ! HDUFP : Gradients of the conservative variables at the host  cell flux pts
  ! GDUFP : Gradients of the conservative variables at the ghost cell flux pts
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),      intent(in) :: hufp
  real(wp), dimension(:,:,:),    intent(in) :: hdufp
  real(wp), dimension(:,:,:), intent(inout) :: gdufp
  !
  !.. Local Scalars ..
  integer :: nr,k,m,l
  real(wp) :: rho_inv,vel
 !real(wp) :: rho_inv2,vel_over_rho
  !
continue
  !
  nr = size(gdufp,dim=1)
  !
  do k = 1,size(gdufp,dim=3)
    !
    ! Copy over the density gradients
    !
    gdufp(1:nr,1,k) = hdufp(1:nr,1,k)
    !
    ! Convert the momentum gradients to velocity gradients
    !
    rho_inv = one/hufp(1,k)
    !
    do m = 2,nq-1
      vel = rho_inv * hufp(m,k)
      do l = 1,nr
        gdufp(l,m,k) = hdufp(l,m,k) - vel*hdufp(l,1,k)
        gdufp(l,m,k) = rho_inv * gdufp(l,m,k)
      end do
    end do
    !
    ! Another algorithm to compute the velocity gradients. Unsure if
    ! performing the above subtraction step first before multiplying
    ! through by rho_inv could cause any precision errors.
    !
   !rho_inv = one/hufp(1,k)
   !rho_inv2 = rho_inv*rho_inv
   !!
   !do m = 2,nq-1
   !  vel_over_rho = rho_inv2 * hufp(m,k)
   !  do l = 1,nr
   !    gdufp(l,m,k) = rho_inv*hdufp(l,m,k) - vel_over_rho*hdufp(l,1,k)
   !  end do
   !end do
    !
    ! Finally, since this is an adiabatic wall boundary,
    ! set the temperature gradients to zero.
    !
    gdufp(1:nr,nq,k) = zero
    !
  end do
  !
end subroutine get_bc_adiabatic_wall_gradient
!
!###############################################################################
!
subroutine get_bc_wall_gradient(hufp,hdufp,gdufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Non-adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! HUFP  : Conservative variables at the host cell flux points
  ! HDUFP : Gradients of the conservative variables at the host  cell flux pts
  ! GDUFP : Gradients of the conservative variables at the ghost cell flux pts
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),      intent(in) :: hufp
  real(wp), dimension(:,:,:),    intent(in) :: hdufp
  real(wp), dimension(:,:,:), intent(inout) :: gdufp
  !
  !.. Local Scalars ..
  integer :: np,nr,k,m,l
  real(wp) :: rho_inv,vel,pvar,cv_inv
 !real(wp) :: rho_inv2,pvar_over_rho
  !
continue
  !
  nr = size(gdufp,dim=1)
  np = size(gdufp,dim=3)
  !
  ! Precompute the inverse of specific heat at constant volume
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Copy over the interior density gradient to the boundary flux points
  !
  do k = 1,np
    do l = 1,nr
      gdufp(l,1,k) = hdufp(l,1,k)
    end do
  end do
  !
  ! We will also copy the other conservative variable gradients from
  ! the interior cell to the boundary flux points. However, since
  ! we treat the viscous fluxes for wall boundaries differently from
  ! the flow-through boundaries, we need to ensure that the wall
  ! gradients are in the right form for the wall viscous flux function.
  ! This means pre-converting the conservative variable gradients to
  ! primitive variable gradients.
  !
  do k = 1,np
    !
    ! Compute the inverse of density
    !
    rho_inv = one/hufp(1,k)
    !
    ! Using the chain rule:
    !
    ! 1. Compute the velocity gradients (m = 2,nq-1)
    ! 2. Compute the total energy gradient and store in the
    !    location for the temperature gradient (m = nq)
    !
    do m = 2,nq
      pvar = rho_inv * hufp(m,k)
      do l = 1,nr
        gdufp(l,m,k) = hdufp(l,m,k) - pvar*hdufp(l,1,k)
        gdufp(l,m,k) = rho_inv * gdufp(l,m,k)
      end do
    end do
   !!
   !! Another algorithm to compute the velocity gradients. Unsure if
   !! performing the above subtraction step first before multiplying
   !! through by rho_inv could cause any precision errors.
   !!
   !rho_inv = one/hufp(1,k)
   !rho_inv2 = rho_inv*rho_inv
   !!
   !do m = 2,nq
   !  pvar_over_rho = rho_inv2 * hufp(m,k)
   !  do l = 1,nr
   !    gdufp(l,m,k) = rho_inv*hdufp(l,m,k) - pvar_over_rho*hdufp(l,1,k)
   !  end do
   !end do
    !
    ! The location for the temperature gradient currently contains
    ! the total energy gradient from the previous loop.
    ! Now use the velocity gradients to subtract the gradient of
    ! kinetic-energy-per-unit-mass from the total energy gradient.
    ! 
    do m = 2,nq-1
      vel = rho_inv * hufp(m,k)
      do l = 1,nr
        gdufp(l,nq,k) = gdufp(l,nq,k) - vel * gdufp(l,m,k)
      end do
    end do
    !
    ! Finish computing the temperature gradient by dividing through by
    ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
    !
    gdufp(1:nr,nq,k) = gdufp(1:nr,nq,k) * cv_inv
    !
  end do
  !
end subroutine get_bc_wall_gradient
!
!###############################################################################
!
subroutine output_bface
  !
  !.. Use Statements ..
  use geovar, only : bface
  !
  !.. Local Scalars ..
  integer :: i,n,maxint,myio,ierr
  character(len=100) :: write_format
  character(len=100) :: file_format
  character(len=100) :: my_file
  !
  !.. Local Arrays ..
  integer, dimension(1:size(bface,dim=1)+1) :: intlen
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_bface"
  !
continue
  !
  myio = 1456 + mypnum
  !
  ! Get the maximum length of the string containing the name of the BC types
  !
  intlen(1) = maxval( len_trim( bc_integer_to_string(bface(1,:)) ) )
  !
  do i = 1,size(intlen)-1
    n = 1
    maxint = maxval(bface(i,:))
    check_bface_size: do
      if (maxint/(10**n) == 0) then
        intlen(i+1) = n+1
        exit check_bface_size
      else
        n = n + 1
        cycle check_bface_size
      end if
    end do check_bface_size
  end do
  !
  ! Primary write of the format statement which has an unlimited 
  ! '*' repetition multiplier on the last quantity to allow for
  ! an arbitrarily large number of values contained in intlen.
  ! The ":" within the part that is unlimitedly repeated serves
  ! to terminate the format statement at that point preventing 
  ! it from going back and adding on an extra ',3x,i' at the end.
  !
  write (write_format,1) (intlen(i),i=1,size(intlen))
 !write (iout,'(a,a,a)') "write_format = '",trim(write_format),"'"
  !
  ! Finally, add on the parentheses to the format statement
  !
  write (write_format,2) trim(write_format)
 !write (iout,'(a,a,a)') "write_format = '",trim(write_format),"'"
  !
  ! Get the number of digits of the integer for the number of cpus
  !
  n = 1
  check_ncpu_size: do
    if (ncpu/(10**n) == 0) then
      write (file_format,3) n
      exit check_ncpu_size
    else
      n = n + 1
      cycle check_ncpu_size
    end if
  end do check_ncpu_size
  !
  ! Create the name of the output file
  !
  write (my_file,fmt=trim(file_format)) mypnum
  !
  ! Open the output file
  !
  open (myio,file=trim(my_file),status="replace",action="write", &
             form="formatted",position="rewind",iostat=ierr)
  call io_error(pname,trim(my_file),1,ierr)
  !
  ! Output the bface array
  !
  do n = 1,size(bface,dim=2)
    write (myio,fmt=trim(write_format)) &
                               trim(bc_integer_to_string(bface(1,n))), &
                                   (bface(i,n),i=1,size(bface,dim=1))
  end do
  !
  ! Close the output file
  !
  close (myio,iostat=ierr)
  call io_error(pname,trim(my_file),2,ierr)
  !
  ! Format Statements
  !
  1 format ('a',i0,',i',i0,',": ",1x,i',i0,',1x,i',i0,',1x,i',i0, &
            ',5x,1x,i',i0,',1x,i',i0,',1x,i',i0,',1x,i',i0,',1x,i',i0, &
            *(',3x,i',i0,:))
  2 format ('(',a,')')
  3 format ('("bface.",i0.',i0,',".dat")')
  !
end subroutine output_bface
!
!###############################################################################
!
subroutine exchange_bc_conflicts_usp(faceusp,override)
  !
  !.. Use Statements ..
  use geovar, only : bc_conflict
  use ovar,   only : iop_flxpts
  use ovar,   only : use_bc_conflicts
  use flowvar, only : face_var
  !
  !.. Formal Arguments ..
  type(face_var), dimension(:), intent(inout) :: faceusp
  logical,        optional,        intent(in) :: override
  !
  !.. Local Scalars ..
  integer :: i,n,nf,km,ks,nfm,nfs
  integer :: itag,send_cpu,recv_cpu
  !
  !.. Local Arrays ..
  integer, dimension(1:mpi_status_size) :: istat
  real(r8), dimension(1:nq) :: exch_send,exch_recv
  !
continue
  !
  ! We have finished computing the gradients on the non-communication
  ! boundary faces and the communication of all other gradients between
  ! processors of adjacent partitions has completed. Now we have to
  ! correct any flux points on boundary faces that have conflicting
  ! boundary conditions.
  !
  ! This only needs to be done if we are using flux points that include
  ! the end points of the face (e.g., Lobatto) so skip all of this and
  ! return if we are not using a quadrature with flux points at the
  ! end points UNLESS override is present. If override is present, only
  ! exit return_loop if it is true. If it is not true, we dont care
  ! that it is present and we only care if we are using Lobatto points.
  !
  return_loop: do i = 1,1
    !
    if (use_bc_conflicts) then
      !
      if (present(override)) then
        if (override) exit return_loop
      end if
      !
      if (any(iop_flxpts == Gauss_Lobatto)) exit return_loop
      !
    end if
    !
    ! Return if we have reached this point
    !
    return
    !
  end do return_loop
  !
  ! Now go through and exchange faceusp at bc conflicts
  !
  do n = 1,size(bc_conflict)
    !
    ! Continue to the next bc_conflict if this processor
    ! has zero involvment with this conflict.
    ! 
    if (bc_conflict(n)%not_my_conflict) cycle
    !
    ! Overwrite the boundary gradient depending on if this conflict
    ! is local or needs to be sent to another processor via MPI.
    !
    if (bc_conflict(n)%is_local) then
      !
      ! This conflict is local to this processor so just
      ! overwrite the boundary gradient in the slave flux points
      ! with the boundary gradient in the master flux point.
      !
      km = bc_conflict(n)%master%flux_point
      nfm = bc_conflict(n)%master%face
      !
      do i = 1,size(bc_conflict(n)%slave)
        !
        ks = bc_conflict(n)%slave(i)%flux_point
        nfs = bc_conflict(n)%slave(i)%face
        !
        faceusp(nfs)%v(:,ks,2) = faceusp(nfm)%v(:,km,2)
        !
      end do
      !
    else
      !
      ! This conflict is on a partition boundary so we need to
      ! use MPI to exchange the boundary gradients.
      !
      send_cpu = bc_conflict(n)%master%cpu
      !
      if (mypnum == send_cpu) then
        !
        km = bc_conflict(n)%master%flux_point
        nfm = bc_conflict(n)%master%face
        !
        exch_send = real(faceusp(nfm)%v(:,km,2),kind=r8)
        !
        do i = 1,size(bc_conflict(n)%slave)
          !
          recv_cpu = bc_conflict(n)%slave(i)%cpu
          !
          ! Off-chance that we have multiple slaves on a conflict
          ! but some are local and some need MPI communication
          !
          if (mypnum == recv_cpu) then
            !
            ks = bc_conflict(n)%slave(i)%flux_point
            nfs = bc_conflict(n)%slave(i)%face
            !
            faceusp(nfs)%v(:,ks,2) = real(exch_send,kind=wp)
            !
          else
            !
            itag = (send_cpu + n + 1)*i + (recv_cpu + i + 1)*n
            call mpi_send(exch_send,size(exch_send),mpi_flttyp,recv_cpu, &
                          itag,mpi_comm_world,mpierr)
            !
          end if
          !
        end do
        !
      else
        !
        do i = 1,size(bc_conflict(n)%slave)
          !
          recv_cpu = bc_conflict(n)%slave(i)%cpu
          !
          if (mypnum == recv_cpu) then
            !
            itag = (send_cpu + n + 1)*i + (recv_cpu + i + 1)*n
            call mpi_recv(exch_recv,size(exch_recv),mpi_flttyp,send_cpu, &
                          itag,mpi_comm_world,istat,mpierr)
            !
            ks = bc_conflict(n)%slave(i)%flux_point
            nfs = bc_conflict(n)%slave(i)%face
            !
            faceusp(nfs)%v(:,ks,2) = real(exch_recv,kind=wp)
            !
          end if
          !
        end do
        !
      end if
      !
    end if
    !
  end do
  !
end subroutine exchange_bc_conflicts_usp
!
!###############################################################################
!
subroutine exchange_bc_conflicts_dusp(facedusp,override)
  !
  !.. Use Statements ..
  use geovar, only : nr,bc_conflict
  use ovar,   only : iop_flxpts
  use ovar,   only : use_bc_conflicts
  use flowvar, only : face_vec
  !
  !.. Formal Arguments ..
  type(face_vec), dimension(:), intent(inout) :: facedusp
  logical,        optional,        intent(in) :: override
  !
  !.. Local Scalars ..
  integer :: i,n,nf,km,ks,nfm,nfs
  integer :: itag,send_cpu,recv_cpu
  !
  !.. Local Arrays ..
  integer, dimension(1:mpi_status_size) :: istat
  real(r8), dimension(1:nr,1:nq) :: exch_send,exch_recv
  !
continue
  !
  ! We have finished computing the gradients on the non-communication
  ! boundary faces and the communication of all other gradients between
  ! processors of adjacent partitions has completed. Now we have to
  ! correct any flux points on boundary faces that have conflicting
  ! boundary conditions.
  !
  ! This only needs to be done if we are using flux points that include
  ! the end points of the face (e.g., Lobatto) so skip all of this and
  ! return if we are not using a quadrature with flux points at the
  ! end points UNLESS override is present. If override is present, only
  ! exit return_loop if it is true. If it is not true, we dont care
  ! that it is present and we only care if we are using Lobatto points.
  !
  return_loop: do i = 1,1
    !
    if (use_bc_conflicts) then
      !
      if (present(override)) then
        if (override) exit return_loop
      end if
      !
      if (any(iop_flxpts == Gauss_Lobatto)) exit return_loop
      !
    end if
    !
    ! Return if we have reached this point
    !
    return
    !
  end do return_loop
  !
  ! Now go through and exchange facedusp at bc conflicts
  !
  do n = 1,size(bc_conflict)
    !
    ! Continue to the next bc_conflict if this processor
    ! has zero involvment with this conflict.
    ! 
    if (bc_conflict(n)%not_my_conflict) cycle
    !
    ! Overwrite the boundary gradient depending on if this conflict
    ! is local or needs to be sent to another processor via MPI.
    !
    if (bc_conflict(n)%is_local) then
      !
      ! This conflict is local to this processor so just
      ! overwrite the boundary gradient in the slave flux points
      ! with the boundary gradient in the master flux point.
      !
      km = bc_conflict(n)%master%flux_point
      nfm = bc_conflict(n)%master%face
      !
      do i = 1,size(bc_conflict(n)%slave)
        !
        ks = bc_conflict(n)%slave(i)%flux_point
        nfs = bc_conflict(n)%slave(i)%face
        !
        facedusp(nfs)%v(:,:,ks,2) = facedusp(nfm)%v(:,:,km,2)
        !
      end do
      !
    else
      !
      ! This conflict is on a partition boundary so we need to
      ! use MPI to exchange the boundary gradients.
      !
      send_cpu = bc_conflict(n)%master%cpu
      !
      if (mypnum == send_cpu) then
        !
        km = bc_conflict(n)%master%flux_point
        nfm = bc_conflict(n)%master%face
        !
        exch_send = real(facedusp(nfm)%v(:,:,km,2),kind=r8)
        !
        do i = 1,size(bc_conflict(n)%slave)
          !
          recv_cpu = bc_conflict(n)%slave(i)%cpu
          !
          ! Off-chance that we have multiple slaves on a conflict
          ! but some are local and some need MPI communication
          !
          if (mypnum == recv_cpu) then
            !
            ks = bc_conflict(n)%slave(i)%flux_point
            nfs = bc_conflict(n)%slave(i)%face
            !
            facedusp(nfs)%v(:,:,ks,2) = real(exch_send,kind=wp)
            !
          else
            !
            itag = (send_cpu + n + 1)*i + (recv_cpu + i + 1)*n
            call mpi_send(exch_send,size(exch_send),mpi_flttyp,recv_cpu, &
                          itag,mpi_comm_world,mpierr)
            !
          end if
          !
        end do
        !
      else
        !
        do i = 1,size(bc_conflict(n)%slave)
          !
          recv_cpu = bc_conflict(n)%slave(i)%cpu
          !
          if (mypnum == recv_cpu) then
            !
            itag = (send_cpu + n + 1)*i + (recv_cpu + i + 1)*n
            call mpi_recv(exch_recv,size(exch_recv),mpi_flttyp,send_cpu, &
                          itag,mpi_comm_world,istat,mpierr)
            !
            ks = bc_conflict(n)%slave(i)%flux_point
            nfs = bc_conflict(n)%slave(i)%face
            !
            facedusp(nfs)%v(:,:,ks,2) = real(exch_recv,kind=wp)
            !
          end if
          !
        end do
        !
      end if
      !
    end if
    !
  end do
  !
end subroutine exchange_bc_conflicts_dusp
!
!###############################################################################
!
include "Functions/pressure.f90"
include "Functions/usp2v.f90"
include "Functions/vsp2u.f90"
!
!###############################################################################
!
end module module_bc
