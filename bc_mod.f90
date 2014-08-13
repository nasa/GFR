module module_bc
  !
  ! Module to contain the procedures and variables related
  ! to the evaluating the boundary conditions for the flow.
  !
  !.. Use Statements ..
  use module_kind_types
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
 !public :: init_bc_module
  !
  !
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  ! WALL_SOLUTION_LOCATION : Factor to determine the location of the
  !                          solution for wall boundaries.
  !
  ! with u_int   being the interior velocity at the wall flux point
  !  and u_wall  being the boundary velocity at the wall flux point
  !  and u_exact being the exact    velocity at the wall flux point,
  !
  ! wall_solution_location = 0 : u_wall = u_exact
  ! wall_solution_location = 1 : u_exact = 1/2 * (u_int + u_wall)
  !
  real(wp), parameter :: wall_solution_location = zero
 !real(wp), parameter :: wall_solution_location = one
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
!subroutine init_bc_module
!  !
!  !.. Use Statements ..
!  use geovar, only : n_sp1de,n_fp1de,neqns,n_sp2dq
!  use geovar, only : interp_s2g,interp_s2f
!  !
!  !.. Local Scalars ..
!  integer :: ierr
!  !
!continue
!  !
!end subroutine init_bc_module
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
  use metvar,  only : geofa
  use flowvar, only : faceusp
  !
  !.. Local Scalars ..
  integer  :: nf,pf
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
    select case (bface(1,nf))
      !
      case (bc_freestream)
        !
        ! Free stream boundary condition
        !
        call get_bc_freestream_solution(faceusp(:,:,2,nf))
        !
      case (bc_characteristic)
        !
        ! Characteristic boundary condition
        !
        call get_bc_characteristic_solution(faceusp(:,:,1,nf), &
                                            faceusp(:,:,2,nf), &
                                            geofa(:,:,nf))
        !
      case (bc_generic_inflow)
        !
        ! Generic inflow boundary condition
        !
        call get_bc_generic_inflow_solution(faceusp(:,:,1,nf), &
                                            faceusp(:,:,2,nf), &
                                            geofa(:,:,nf))
        !
      case (bc_generic_outflow)
        !
        ! Generic outflow boundary condition
        !
        call get_bc_generic_outflow_solution(faceusp(:,:,1,nf), &
                                             faceusp(:,:,2,nf), &
                                             geofa(:,:,nf))
        !
      case (bc_generic_freeflow)
        !
        ! Generic freeflow boundary condition
        !
        call get_bc_generic_freeflow_solution(faceusp(:,:,1,nf), &
                                              faceusp(:,:,2,nf), &
                                              geofa(:,:,nf))
        !
      case (bc_sub_inflow)
        !
        ! Subsonic inflow boundary condition
        !
        call get_bc_sub_inflow_solution(faceusp(:,:,1,nf), &
                                        faceusp(:,:,2,nf))
        !
      case (bc_sub_outflow)
        !
        ! Subsonic outflow boundary condition
        !
        call get_bc_sub_outflow_solution(faceusp(:,:,1,nf), &
                                         faceusp(:,:,2,nf))
        !
      case (bc_sup_inflow)
        !
        ! Supersonic inflow boundary condition
        !
        call get_bc_sup_inflow_solution(faceusp(:,:,2,nf))
        !
      case (bc_sup_outflow)
        !
        ! Supersonic outflow boundary condition
        !
        call get_bc_sup_outflow_solution(faceusp(:,:,1,nf), &
                                         faceusp(:,:,2,nf))
        !
      case (bc_slip_wall,bc_symmetry)
        !
        ! Slip wall / symmetry boundary condition
        !
        call get_bc_slip_wall_solution(faceusp(:,:,1,nf), &
                                       faceusp(:,:,2,nf), &
                                       geofa(:,:,nf))
        !
      case (bc_adiabatic_wall)
        !
        ! No-slip adiabatic wall
        !
        call get_bc_adiabatic_wall_solution(faceusp(:,:,1,nf), &
                                            faceusp(:,:,2,nf))
        !
      case (bc_isothermal_wall)
        !
        ! No-slip isothermal wall
        !
        call get_bc_isothermal_wall_solution(faceusp(:,:,1,nf), &
                                             faceusp(:,:,2,nf))
        !
      case (bc_periodic)
        !
        ! Periodic boundary condition
        !
        pf = bface(5,nf) ! index of connecting periodic boundary face
        !
        call get_bc_periodic_solution(faceusp(:,:,1,pf), &
                                      faceusp(:,:,2,nf))
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
        call stop_recon(abort,"get_faceusp_bc")
        !
    end select
    !
  end do
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
  use metvar,  only : geofa
  use flowvar, only : facedusp
  !
  !.. Local Scalars ..
  integer  :: nf,pf
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
    select case (bface(1,nf))
      !
     !case (bc_flowthrough)
        !
        ! General flow-through boundary conditions
        !
       !call get_bc_flowthrough_gradient(facedusp(:,:,:,1,nf), &
       !                                 facedusp(:,:,:,2,nf))
        !
      case (bc_slip_wall,bc_symmetry)
        !
        ! Slip wall / symmetry boundary condition
        !
       !call get_bc_slip_wall_gradient(facedusp(:,:,:,1,nf), &
       !                               facedusp(:,:,:,2,nf))
        !
      case (bc_adiabatic_wall)
        !
        ! No-slip adiabatic wall
        !
       !call get_bc_adiabatic_wall_gradient(facedusp(:,:,:,1,nf), &
       !                                    facedusp(:,:,:,2,nf))
        !
      case (bc_isothermal_wall)
        !
        ! No-slip isothermal wall
        !
       !call get_bc_isothermal_wall_gradient(facedusp(:,:,:,1,nf), &
       !                                     facedusp(:,:,:,2,nf))
        !
      case (bc_periodic)
        !
        ! Periodic boundary condition
        !
        pf = bface(5,nf) ! index of connecting periodic boundary face
        !
       !call get_bc_periodic_gradient(facedusp(:,:,:,1,pf), &
       !                              facedusp(:,:,:,2,nf))
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
        call stop_recon(abort,"get_facedusp_bc")
        !
    end select
    !
  end do
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
  use generic_mod, only : usp2v_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nq,k
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
  nq = size(gufp,dim=1)
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
    vn_int = dot_product( rnrm , hvfp(2:nq-1) )
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
      vn_ext = dot_product( rnrm , pv_in(2:nq-1) )
      !
     !if (mach_int > zero) then ! this is an inflow boundary
     !  !
     !  rp = vn_ext + two*aspd_ext*gm1_inv
     !  rm = vn_int - two*aspd_int*gm1_inv
     !  unc = half*(rp+rm)
     !  utc = dot_product( rtan , pv_in(2:nq-1) )
     !  ac = one4*gm1*(rp-rm)
     !  sc = pv_in(nq)/(pv_in(1)**gam)
     !  !
     !else ! this is an outflow boundary
     !  !
     !  rp = -abs(vn_ext) + two*aspd_ext*gm1_inv
     !  rm =      vn_int  - two*aspd_int*gm1_inv
     !  unc = half*(rp+rm)
     !  utc = dot_product( rtan , hvfp(2:nq-1) )
     !  ac = one4*gm1*(rp-rm)
     !  sc = hvfp(nq)/(hvfp(1)**gam)
     !  !
     !end if
      !
      ! Compute the variables dependent on the interior flow direction
      !
      if (mach_int > zero) then
        rp   = vn_ext
        utc  = dot_product( rtan , pv_in(2:nq-1) )
        rho  = pv_in(1)
        pres = pv_in(nq)
      else
        rp   = -abs(vn_ext)
        utc  = dot_product( rtan , hvfp(2:nq-1) )
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
      gufp(2:nq-1,k) = rho * (unc*rnrm(:) + utc*rtan(:))
      gufp(nq,k) = pres*gm1_inv + half * selfdot( gufp(2:nq-1,k) ) / rho
      !
    end if
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
  use generic_mod, only : usp2v_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nr,nq,k
  real(wp) :: gm1_inv,ske_inf,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  nr = size(fnrm,dim=1)
  nq = size(gufp,dim=1)
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  ske_inf = half * pv_in(1) * selfdot( pv_in(2:nq-1) )
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
    mach_int = dot_product( hvfp(2:nq-1) , rnrm(1:nr) ) / aspd_int
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
  use generic_mod, only : usp2v_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nr,nq,k
  real(wp) :: gm1_inv,ske_int,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  nr = size(fnrm,dim=1)
  nq = size(gufp,dim=1)
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
    mach_int = dot_product( hvfp(2:nq-1) , rnrm(1:nr) ) / aspd_int
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
      ske_int = half * selfdot( hufp(2:nq-1,k) ) / hufp(1,k)
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
  use generic_mod, only : usp2v_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nr,nq,k
  real(wp) :: gm1_inv,ske_inf,ske_int,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  nr = size(fnrm,dim=1)
  nq = size(gufp,dim=1)
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  ske_inf = half * pv_in(1) * selfdot( pv_in(2:nq-1) )
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
    mach_int = dot_product( hvfp(2:nq-1) , rnrm(1:nr) ) / aspd_int
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
      ske_int = half * selfdot( hufp(2:nq-1,k) ) / hufp(1,k)
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
  use ovar,        only : pv_in,cv_in,gam
  use generic_mod, only : pressure_cv_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer  :: nq,k
  real(wp) :: pres_int,ske_inf,gm1_inv
  !
continue
  !
  nq = size(gufp,dim=1)
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  ske_inf = half * pv_in(1) * selfdot( pv_in(2:nq-1) )
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
  integer :: nq,k
  real(wp) :: gm1_inv,ske_int
  !
continue
  !
  nq = size(gufp,dim=1)
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
    ske_int = half * selfdot( hufp(2:nq-1,k) ) / hufp(1,k)
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
  integer :: nq,k
  !
continue
  !
  nq = size(gufp,dim=1)
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
  use generic_mod, only : usp2v_sp,vsp2u_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nq,k
  real(wp) :: vn
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  nq = size(gufp,dim=1)
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
    vn = dot_product( hvfp(2:nq-1) , rnrm )
    !
    ! For velocities, subtract (wall_solution_location+one) * (the wall
    ! normal velocity) from the interior velocity to get the wall
    ! boundary velocity.
    !
    hvfp(2:nq-1) = hvfp(2:nq-1) - (wall_solution_location+one)*vn*rnrm(:)
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
  use generic_mod, only : usp2v_sp,vsp2u_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nq,k
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  nq = size(gufp,dim=1)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Get the interior primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Apply the no-slip condition to the boundary velocities
    !
    hvfp(2:nq-1) = -hvfp(2:nq-1) * wall_solution_location
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
  use generic_mod, only : usp2v_sp,vsp2u_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nq,k
  real(wp) :: iso_wall_temp_nd
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
continue
  !
  nq = size(gufp,dim=1)
  !
  iso_wall_temp_nd = iso_wall_temp / tref
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Get the interior primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Apply the no-slip condition to the boundary velocities
    !
    hvfp(2:nq-1) = -hvfp(2:nq-1) * wall_solution_location
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
subroutine get_bc_slip_wall_gradient(hufp,gufp,fnrm)
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
  use generic_mod, only : usp2v_sp,vsp2u_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: fnrm
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nq,k
  real(wp) :: vn
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(fnrm,dim=1)) :: rnrm
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
  ! WALL_VELOCITY_FACTOR : Factor to determine what velocity to
  !                        compute for the wall boundary condition.
  !
  ! with u_int   being the interior velocity at the wall flux point
  !  and u_wall  being the boundary velocity at the wall flux point
  !  and u_exact being the exact, tangent-only velocity at the wall flux point,
  !
  ! wall_velocity_factor = 1 : u_wall = u_exact
  ! wall_velocity_factor = 2 : u_exact = 1/2 * (u_int + u_wall)
  !
 !real(wp), parameter :: wall_velocity_factor = one
  real(wp), parameter :: wall_velocity_factor = two
  !
continue
  !
  nq = size(gufp,dim=1)
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
    vn = dot_product( hvfp(2:nq-1) , rnrm )
    !
    ! For velocities, subtract (wall_velocity_factor) * (the wall
    ! normal velocity) from the interior velocity to get the wall
    ! boundary velocity.
    !
    hvfp(2:nq-1) = hvfp(2:nq-1) - wall_velocity_factor*vn*rnrm(:)
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
end subroutine get_bc_slip_wall_gradient
!
!###############################################################################
!
subroutine get_bc_adiabatic_wall_gradient(hdufp,gdufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! HDUFP : Gradients of the conservative variables at the host  cell flux pts
  ! GDUFP : Gradients of the conservative variables at the ghost cell flux pts
  !
  !.. Use Statements ..
  use generic_mod, only : usp2v_sp,vsp2u_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:,:),    intent(in) :: hdufp
  real(wp), dimension(:,:,:), intent(inout) :: gdufp
  !
  !.. Local Scalars ..
  integer :: nq,k
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(gdufp,dim=1)) :: hvfp
  !
continue
  !
 !nr = size(gdufp,dim=1)
  nq = size(gdufp,dim=2)
  !
  do k = 1,size(gdufp,dim=3)
   !!
   !! Get the interior primitive variables
   !!
   !hvfp(:) = usp2v_sp( hufp(:,k) )
   !!
   !! Apply the no-slip condition to the boundary velocities
   !!
   !hvfp(2:nq-1) = -hvfp(2:nq-1) * wall_velocity_factor
   !!
   !! Convert the primitive variables to conservative.
   !! This ensures that the interior density and pressure are
   !! extrapolated to the boundary and the total energy at the
   !! boundary is consistent with these two states and the
   !! computed boundary velocity.
   !!
   !gufp(:,k) = vsp2u_sp( hvfp(:) )
   !!
  end do
  !
end subroutine get_bc_adiabatic_wall_gradient
!
!###############################################################################
!
subroutine get_bc_isothermal_wall_gradient(hufp,gufp)
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
  use generic_mod, only : usp2v_sp,vsp2u_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: hufp
  real(wp), dimension(:,:), intent(inout) :: gufp
  !
  !.. Local Scalars ..
  integer :: nq,k
  real(wp) :: iso_wall_temp_nd
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(gufp,dim=1)) :: hvfp
  !
  !
  ! WALL_VELOCITY_FACTOR : Factor to determine what velocity to
  !                        compute for the wall boundary condition.
  !
  ! with u_int   being the interior velocity at the wall flux point
  !  and u_wall  being the boundary velocity at the wall flux point
  !  and u_exact being the exact, no-slip velocity at the wall flux point,
  !
  ! wall_velocity_factor = 0 : u_wall = u_exact
  ! wall_velocity_factor = 1 : u_exact = 1/2 * (u_int + u_wall)
  !
 !real(wp), parameter :: wall_velocity_factor = zero
  real(wp), parameter :: wall_velocity_factor = one
  !
continue
  !
  nq = size(gufp,dim=1)
  !
  iso_wall_temp_nd = iso_wall_temp / tref
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Get the interior primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Apply the no-slip condition to the boundary velocities
    !
    hvfp(2:nq-1) = -hvfp(2:nq-1) * wall_velocity_factor
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
end subroutine get_bc_isothermal_wall_gradient
!
!###############################################################################
!
subroutine get_bc_periodic_gradient(pufp,gufp)
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
end subroutine get_bc_periodic_gradient
!
!###############################################################################
!
end module module_bc
