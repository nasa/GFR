!!!##############################################################################
!
! > A collection of routines for enforcing specific boundedness or minimum value
! > constratins on the solution. This is motivated by the 
! > Positivity-Preservation, entropy limiting, and Invariant Domain Preserving
! > Bodies of work. (see Zhang & Shu for literature)
!
!
! > Density pressure positivity 
!       [Zhang Shu]
!
! > Internal energy positivty
!       [Jenre: Johnson, Ching]
!
! SETH: THIS AUTOMATICALLY CHECKS FACES, WHICH IS NEEDED TO REMAIN ROBUST.
!       IT MIGHT BE A GOOD IDEA TO ALSO CHECK CELL P0 NODES.
!!!##############################################################################
!
module positivity_mod 
  !
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  implicit none
  !
  private
  !
  public :: positivity_preserving_limiter
  !
  private :: limit_density
  private :: limit_species_densities_by_mixture
  private :: limit_species_densities_to_zero
  private :: limit_cv_by_pressure
  private :: limit_internal_energy
  !
contains
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
! //////////////////////////////////////////////////////////////////////////// !
! //////////////////////////////////////////////////////////////////////////// !
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!
subroutine positivity_preserving_limiter(nst)
  ! Apply limiters on thermodynamic quantities
  !
  !.. Use Statements ..
  !
  use cell_ave_mod, only : allocate_cell_ave, deallocate_cell_ave
  use cell_ave_mod, only : compute_cell_ave_sol, cell_ave_sol, cell_ave_gp
  use eqn_idx_mod,  only : nq,ncs
  use flowvar_mod,  only : usp,gpsp,faceusp,facegpsp
  use flowvar_mod,  only : get_marker_index,nmarker
  use flowvar_mod,  only : positivity_limiter_theta
  use flux_mod,     only : interp_cv_to_flxpts
  use geomvar_mod,  only : ncell
  use ovar,         only : positivity,rhoref,aref
  use thermo_mod,   only : update_thermo_properties
  !
  !.. Input ..
  integer, intent(in)   :: nst
  !
  !.. Local Scalars ..
  real(wp)    :: dens_tolerance, U_tolerance, limpres
  integer     :: imarker,ierr
  logical(lk) :: do_lim_flxpts
  !
  character(len=*), parameter :: pname = "positivity_preserving_limiter"
  !
  continue
  !
  ! Allocate the array storing the amount of limitng occuring
  if (.not. allocated(positivity_limiter_theta)) then
    allocate ( positivity_limiter_theta(1:ncell,1:nmarker) , source=one , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"positivity_limiter_theta",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !$ACC UPDATE DEVICE(positivity_limiter_theta)
  end if
  !
  imarker = get_marker_index(nst)
  do_lim_flxpts = true ! Checking face points is needed for robustness
  !
  ! ___________________________________
  ! >>> Get Cell Average quantities <<<
  !
  call allocate_cell_ave()  !!! NOTE: Make sure cell ave is deallocated
  ! NOTE: This also updates cell_ave_gp
  call compute_cell_ave_sol()
  !
  ! _______________________________________
  ! >>> Get Limiting Pressure / Density <<<
  !
  dens_tolerance  = positivity%Density_Limit / rhoref
  limpres = positivity%Pressure_Limit / (rhoref*aref*aref) 
  U_tolerance  = limpres 
  ! >>> Pressure is used as a stand in for internal energy
  ! >>>> this works because 0 < P < U
  !
  call interp_cv_to_flxpts(0,0,true) ! only interp CV, do not calculate GP
  !
  ! Density Limiter
  if (ncs>0) then  ! Special order of operations / procedure for multispecies
    !
    call limit_species_densities_by_mixture(usp,cell_ave_sol,dens_tolerance,&
                    do_lim_flxpts,faceusp,positivity_limiter_theta(:,imarker))
    call limit_species_densities_to_zero(usp,cell_ave_sol,true,faceusp,&
                                          positivity_limiter_theta(:,imarker))
    !
    !
    !
  else ! single species flow
    !
    call limit_density(usp,cell_ave_sol,dens_tolerance,do_lim_flxpts,faceusp,&
                       positivity_limiter_theta(:,imarker))
    !
  end if
  !
  !
  ! NOTE: Technically internal energy limiter is the only one needed for 
  !       a robust solution. The pressure limiter enables extra
  !       dissipation via the user input limit. This lets you run at higher
  !       CFL numbers than strict positivity.
  !
  ! Internal Energy Limiter 
  call limit_internal_energy(usp,cell_ave_sol,U_tolerance,do_lim_flxpts,&
                             faceusp,positivity_limiter_theta(:,imarker))
  !
  !
  ! Pressure Limiter 
  ! SETH: Add update for gpsp and facegpsp
  ! SETH: Make optional, probably enabled by default, more useful for higher speed flows with
  ! shocks.
  call limit_cv_by_pressure(usp,gpsp,cell_ave_sol,cell_ave_gp,&
    limpres,do_lim_flxpts,facegpsp,faceusp,positivity_limiter_theta(:,imarker))
  !
  ! Deallocate the cell_ave_sol and cell_ave_gp arrays
  !#############################################################################
  !#############################################################################
  !#############################################################################
  ! NOTE: Not sure if this is needed
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  call deallocate_cell_ave()
  !
end subroutine positivity_preserving_limiter
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
! //////////////////////////////////////////////////////////////////////////// !
! //////////////////////////////////////////////////////////////////////////// !
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine limit_density(usp,cell_ave_sol,limrho,limit_flxpts,faceusp,plt)
  !
  !-----------------------------------------------------------------------------
  ! PURPOSE: 
  !   Ensure the positivity of density by scaling about the cell average
  !-----------------------------------------------------------------------------
  ! ASSUMPTIONS: 
  !   INPUT:
  !   -> This is a single-species flow (ncs > 0)
  !   -> CV has been interpolated to flux points if limit_flxpts=T
  !   OUTPUT:
  !   <- faceusp  is     valid after this routine
  !-----------------------------------------------------------------------------
  ! Notes:
  !-----------------------------------------------------------------------------
  ! .. Use Statements  ..
  use geomvar_mod,   only : ncell,cell,face,n_solpts
  use eqn_idx_mod,   only : nec
  !
  ! .. Arguments       ..
  logical(lk), intent(in)                 :: limit_flxpts ! do we treat faces
  real(wp), intent(in)                    :: limrho      ! limiting density
  real(wp), dimension(:,:), intent(inout) :: usp
  real(wp), dimension(:,:), intent(in)    :: cell_ave_sol
  real(wp), dimension(:,:,:),intent(inout),optional   :: faceusp
  real(wp), dimension(:), intent(inout)   :: plt
  !
  ! .. Local Variables ..
  character(len=*), parameter :: pname = "limit_density"
  integer  :: ic,n,n1,n2,m,ncf,nf,sf,f1,f2
  real(wp) :: averho,minrho,theta
  !
continue
  !
  ! Loop Over Elements
  do ic = 1,ncell
    n1 = cell(ic)%beg_sp ! beginning solution point index for this cell
    n2 = cell(ic)%end_sp ! ending    solution point index for this cell
    !
    ! Find scaling factor based on mixture density
    averho = cell_ave_sol(ic,nec)
    minrho = minval(usp(n1:n2,nec))
    !
    if (limit_flxpts .and. present(faceusp)) then
    do ncf = 1,size(cell(ic)%face)
      nf  = cell(ic)%face(ncf)%idx  ! index of the current face
      sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
      f1  = face(nf)%beg_fp
      f2  = face(nf)%end_fp
      ! over the face points
      minrho  = min(minrho ,minval(faceusp(f1:f2,nec,sf)))
    end do
    end if
    !
    if (minrho >= limrho) cycle
    if (abs(averho-minrho) <= ten*rndoff) cycle
    !
    !SETH: MAKE SURE DENOMINATOR IS NONZERO
    theta  = (averho - limrho) / (averho - minrho)  ! Z-S scaling factor
    theta = max(min(theta,one),zero)
    !
    ! >> Perform the Limiting <<
    !
    ! Scale the solution density
    usp(n1:n2,nec) = cell_ave_sol(ic,nec) + &
                             theta*(usp(n1:n2,nec)-cell_ave_sol(ic,nec))
    !
    ! Update the output marker
    plt(ic) = min(plt(ic),theta)
    !
    !
    !
    if (limit_flxpts .and. present(faceusp)) then
      do ncf = 1,size(cell(ic)%face)
        nf  = cell(ic)%face(ncf)%idx  ! index of the current face
        sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
        f1  = face(nf)%beg_fp
        f2  = face(nf)%end_fp
        faceusp(f1:f2,nec,sf) = cell_ave_sol(ic,nec) + &
                           theta*(faceusp(f1:f2,nec,sf)-cell_ave_sol(ic,nec))
      end do
    end if
    !
  end do
  !
end subroutine limit_density
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
! //////////////////////////////////////////////////////////////////////////// !
! //////////////////////////////////////////////////////////////////////////// !
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine limit_species_densities_by_mixture(usp,cell_ave_sol,limrho, &
                                              limit_flxpts,faceusp,plt)
  !
  !-----------------------------------------------------------------------------
  ! PURPOSE: 
  !   Ensure the positivity of mixture density by modifying the species'.
  !-----------------------------------------------------------------------------
  ! ASSUMPTIONS: 
  !   INPUT:
  !   -> Mixture density has already been calculated -> u(:,nec) is valid
  !   -> This is a multispecies flow (ncs > 0)
  !   -> CV has been interpolated to flux points if limit_flxpts=T
  !   OUTPUT:
  !   <- u(:,nec) is     valid after this routine
  !   <- faceusp  is     valid after this routine
  !-----------------------------------------------------------------------------
  ! Notes:
  !   It's okay if an individual species density is negative. It's important
  !    to NOT treat / fix them before this.
  !-----------------------------------------------------------------------------
  ! .. Use Statements  ..
  use geomvar_mod,   only : ncell,cell,face,n_solpts
  use eqn_idx_mod,   only : nec
  !
  ! .. Arguments       ..
  logical(lk), intent(in)                 :: limit_flxpts ! do we treat faces
  real(wp), intent(in)                    :: limrho      ! limiting density
  real(wp), dimension(:,:), intent(inout) :: usp
  real(wp), dimension(:,:), intent(in)    :: cell_ave_sol
  real(wp), dimension(:,:,:),intent(inout),optional   :: faceusp
  real(wp), dimension(:), intent(inout)   :: plt
  !
  ! .. Local Variables ..
  integer  :: ic,n,n1,n2,m,ncf,nf,sf,f1,f2
  real(wp) :: averho,minrho,theta
  !
continue
  !
  ! Loop Over Elements
  do ic = 1,ncell
    n1 = cell(ic)%beg_sp ! beginning solution point index for this cell
    n2 = cell(ic)%end_sp ! ending    solution point index for this cell
    !
    !
    ! >> Find minimum density on set of points of interest
    averho = cell_ave_sol(ic,nec)
    minrho = minval(usp(n1:n2,nec))
    !
    if (limit_flxpts .and. present(faceusp)) then
    do ncf = 1,size(cell(ic)%face)
      nf  = cell(ic)%face(ncf)%idx  ! index of the current face
      sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
      f1  = face(nf)%beg_fp
      f2  = face(nf)%end_fp
      ! over the face points
      minrho  = min(minrho ,minval(faceusp(f1:f2,nec,sf)))
    end do
    end if
    !
    ! >> Find scaling factor needed to limit density
    !
    if (minrho > limrho) cycle
    if ((averho-minrho) < ten*rndoff) then
      !
      theta = zero
      !
      else
      !
      theta  = 0.99_wp * (averho-limrho) / (averho - minrho)  ! Z-S scaling factor
      theta = max(min(theta,one),zero)
      !
    end if
    !
    ! Update the output marker
    plt(ic) = min(plt(ic),theta)
    !
    !
    ! >> Perform the Limiting <<
    !
    ! Scale the individual species
    do m=1,nec
      usp(n1:n2,m) = cell_ave_sol(ic,m) + &
                                     theta*(usp(n1:n2,m)-cell_ave_sol(ic,m))
    end do
    !
    if (limit_flxpts .and. present(faceusp)) then
      do m=1,nec
      do ncf = 1,size(cell(ic)%face)
        nf  = cell(ic)%face(ncf)%idx  ! index of the current face
        sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
        f1  = face(nf)%beg_fp
        f2  = face(nf)%end_fp
        faceusp(f1:f2,m,sf) = cell_ave_sol(ic,m) + &
                           theta*(faceusp(f1:f2,m,sf)-cell_ave_sol(ic,m))
      end do
      end do
    end if
    !
  end do
  !
end subroutine limit_species_densities_by_mixture
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
! //////////////////////////////////////////////////////////////////////////// !
! //////////////////////////////////////////////////////////////////////////// !
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine limit_species_densities_to_zero(usp,cell_ave_sol,lim_flxpts,&
                                           faceusp,plt)
  !
  !-----------------------------------------------------------------------------
  ! PURPOSE: 
  !   Ensure the non-negativity of species density.
  !-----------------------------------------------------------------------------
  ! ASSUMPTIONS: 
  !   INPUT:
  !   -> This is a multispecies flow (ncs > 0)
  !   -> CV has been interpolated to flux points if limit_flxpts=T
  !   -> Mixture density has already been limited. This is for species positivit
  !   OUTPUT:
  !   <- usp     is     valid after call
  !   <- faceusp is     valid after call
  !-----------------------------------------------------------------------------
  ! Notes:
  !   This is where we fix the negative species density AFTER limiting mixture 
  !     density. I.E. this help ensure a valid state for chemistry.
  !-----------------------------------------------------------------------------
  !
  ! .. Use Statements  ..
  use geomvar_mod,   only : ncell,cell,face,n_solpts
  use thermo_mod,    only : update_mixture_density
  use eqn_idx_mod,   only : ncs
  !
  ! .. Arguments       ..
  logical(lk), intent(in)                 :: lim_flxpts ! do we treat faces
  real(wp), dimension(:,:), intent(inout) :: usp
  real(wp), dimension(:,:), intent(in)    :: cell_ave_sol
  real(wp), dimension(:,:,:),intent(inout),optional   :: faceusp
  real(wp), dimension(:), intent(inout)   :: plt
  !
  ! .. Local Variables ..
  integer  :: ic,n,n1,n2,m,ncf,nf,sf,f1,f2
  logical(lk), allocatable :: bad_cel_loc(:)
  real(wp) :: averho,minrho,theta
  !
continue
  !
  if (.not. allocated(bad_cel_loc)) then
    allocate(bad_cel_loc(1:ncell), source = fals)
  end if
  !
  ! Loop Over Species and Elements
  do m=1,ncs
  do ic = 1,ncell
    n1 = cell(ic)%beg_sp ! beginning solution point index for this cell
    n2 = cell(ic)%end_sp ! ending    solution point index for this cell
    !
    ! >> Find minimum species density on set of points of interest
    !
    averho = cell_ave_sol(ic,m)
    minrho = minval(usp(n1:n2,m))
    !
    if (lim_flxpts .and. present(faceusp)) then
    do ncf = 1,size(cell(ic)%face)
      nf  = cell(ic)%face(ncf)%idx  ! index of the current face
      sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
      f1  = face(nf)%beg_fp
      f2  = face(nf)%end_fp
      ! over the face points
      minrho  = min(minrho ,minval(faceusp(f1:f2,m,sf)))
    end do
    end if
    !
    !
    ! >> Find scaling factor needed to limit species density
    !
    if (minrho > zero) cycle
    if ((averho-minrho) < ten*rndoff) then
      !
      theta = zero
      !
      else
      !
      theta  = 0.99_wp*(averho-rndoff) / (averho - minrho)  ! Z-S scaling factor
      theta = max(min(theta,one),zero)
      !
    end if
    !
    if (averho < zero) then
      averho = zero
      theta = zero
    end if
    !
    ! Update the output marker
    plt(ic) = min(plt(ic),theta)
    bad_cel_loc(ic) = true
    !
    !
    ! >> Perform the Limiting <<
    !
    ! Scale the individual species
    usp(n1:n2,m) = cell_ave_sol(ic,m) + &
                              theta*(usp(n1:n2,m)-cell_ave_sol(ic,m))
    usp(n1:n2,m) = max(usp(n1:n2,m), zero)
    !
    if (lim_flxpts .and. present(faceusp)) then
      do ncf = 1,size(cell(ic)%face)
        nf  = cell(ic)%face(ncf)%idx  ! index of the current face
        sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
        f1  = face(nf)%beg_fp
        f2  = face(nf)%end_fp
        faceusp(f1:f2,m,sf) = cell_ave_sol(ic,m) + &
                           theta*(faceusp(f1:f2,m,sf)-cell_ave_sol(ic,m))
        faceusp(f1:f2,m,sf) = max(faceusp(f1:f2,m,sf), zero)
      end do
    end if
    !
  end do
  end do
  !
  ! >> Update the mixture density
  do ic = 1,ncell
    if (bad_cel_loc(ic)) then
      n1 = cell(ic)%beg_sp ! beginning solution point index for this cell
      n2 = cell(ic)%end_sp ! ending    solution point index for this cell
      call update_mixture_density(n1,n2,usp)
      !
      if (lim_flxpts .and. present(faceusp)) then
        do ncf = 1,size(cell(ic)%face)
          nf  = cell(ic)%face(ncf)%idx  ! index of the current face
          sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
          f1  = face(nf)%beg_fp
          f2  = face(nf)%end_fp
          call update_mixture_density(1,f2-f1+1,faceusp(f1:f2,:,sf))
        end do
      end if
    end if
  end do
  !
  if ( allocated(bad_cel_loc)) then
    deallocate(bad_cel_loc)
  end if
  !
end subroutine limit_species_densities_to_zero
!
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
! //////////////////////////////////////////////////////////////////////////// !
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
! //////////////////////////////////////////////////////////////////////////// !
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine limit_cv_by_pressure(usp,gpsp,cell_ave_sol,cell_ave_gp,&
                                limpres,lim_flxpts,facegpsp,faceusp,plt)
  !
  !-----------------------------------------------------------------------------
  ! PURPOSE: 
  !   Ensure the non-negativity of pressure. (Step 3 of Zhang-huS PP Limiter) 
  !-----------------------------------------------------------------------------
  ! ASSUMPTIONS: 
  !   INPUT:
  !   -> This can be a multispecies flow (ncs > 0)
  !   -> Pressure, gpsp(:,npres), is valid
  !   -> Has been interpolated to flux points if limit_flxpts=T
  !   -> facegpsp is valid (dependant on bool in interp_to_cv call)
  !   OUTPUT:
  !   <- faceusp  IS     valid after routine
  !   <- facegsp  IS     valid after routine
  !   <- usp      IS     valid after routine
  !   <- gpsp     IS     valid after routine
  !-----------------------------------------------------------------------------
  ! Notes:
  !  -In theory this preceeded by the density limiting should ensure robustness
  !     of the solution.
  !-----------------------------------------------------------------------------
  !
  ! .. Use Statements  ..
  use eqn_idx_mod,   only : npres,nq
  use geomvar_mod,   only : ncell,cell,face,n_solpts
  use thermo_mod,    only : update_thermo_properties
  !
  ! .. Arguments       ..
  logical(lk), intent(in)                           :: lim_flxpts   ! faces?
  real(wp), intent(in)                              :: limpres     ! limiting P
  real(wp), dimension(:,:), intent(inout)           :: usp
  real(wp), dimension(:,:), intent(inout)           :: gpsp
  real(wp), dimension(:,:), intent(in)              :: cell_ave_sol
  real(wp), dimension(:,:), intent(in)              :: cell_ave_gp
  real(wp), dimension(:,:,:),intent(inout),optional :: facegpsp
  real(wp), dimension(:,:,:),intent(inout),optional :: faceusp
  real(wp), dimension(:), intent(inout)             :: plt
  !
  ! .. Local Variables ..
  integer  :: ic,n,n1,n2,ncf,nf,sf,f1,f2,m
  integer  :: numbad, irank, ierr
  real(wp) :: avepres,minpres,theta,mincellpres,minfacepres
  !
  character(len=*), parameter :: pname = "limit_cv_by_pressure"
continue
  !
  numbad = 0
  !
  ! Loop all Cells
  do ic = 1,ncell
    n1 = cell(ic)%beg_sp ! beginning solution point index for this cell
    n2 = cell(ic)%end_sp ! ending    solution point index for this cell
    !
    ! >> Find the Minimum pressure onthe point set of interest <<
    !
    avepres = cell_ave_gp(ic,npres)
    mincellpres = minval(gpsp(n1:n2,npres))
    minfacepres = 999.0_wp + 999_wp*mincellpres
    !
    if (lim_flxpts .and. present(facegpsp)) then
    do ncf = 1,size(cell(ic)%face)
      nf  = cell(ic)%face(ncf)%idx  ! index of the current face
      sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
      f1  = face(nf)%beg_fp
      f2  = face(nf)%end_fp
      ! over the face points
      minfacepres  = min(minfacepres ,minval(facegpsp(f1:f2,npres,sf)))
    end do
    end if
    !
    minpres = min(minfacepres,mincellpres)
    !
    !
    ! >> Find the damping ratio needed to limit the cell <<
    !
    if (minpres > limpres) cycle
    numbad = numbad + 1
    if ((avepres-minpres) < ten*rndoff) then
        theta = zero
    else
        !
        theta  = 0.99_wp*(avepres-limpres) / (avepres-minpres)  ! Z-S scaling factor
        theta = max(min(theta,one),zero)
        !
    end if
    !
    ! Update the output marker
    plt(ic) = min(plt(ic),theta)
    !
    !
    ! >> Perform the Pressure Limiting <<
    !
    ! Scale all CV
    do m=1,nq
    usp(n1:n2,m) = cell_ave_sol(ic,m) + &
                              theta*(usp(n1:n2,m)-cell_ave_sol(ic,m))
    end do
    !
    call update_thermo_properties(n1,n2,usp,gpsp)
    !
    if (lim_flxpts .and. present(faceusp)) then
      do ncf = 1,size(cell(ic)%face)
        nf  = cell(ic)%face(ncf)%idx  ! index of the current face
        sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
        f1  = face(nf)%beg_fp
        f2  = face(nf)%end_fp
        !
        do m=1,nq
          faceusp(f1:f2,m,sf) = cell_ave_sol(ic,m) + &
                             theta*(faceusp(f1:f2,m,sf)-cell_ave_sol(ic,m))
        end do
        !
        call update_thermo_properties(f1,f2,faceusp(:,:,sf),&
                                facegpsp(:,:,sf),"pressure_positivity")
        !
      end do
    end if
    !
  end do
  !
end subroutine limit_cv_by_pressure
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
! //////////////////////////////////////////////////////////////////////////// !
! //////////////////////////////////////////////////////////////////////////// !
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine limit_internal_energy(usp,cell_ave_sol,Ulim,lim_flxpts,faceusp,&
                                 plt)
  !
  !-----------------------------------------------------------------------------
  ! PURPOSE: 
  !   Ensure the non-negativity of internal energy.
  !-----------------------------------------------------------------------------
  ! ASSUMPTIONS: 
  !   INPUT:
  !   -> Either this is a single species flow, or the mixture density is good
  !   -> CV has been interpolated to flux points if limit_flxpts=T
  !   OUTPUT:
  !   <- usp     is     valid after call
  !   <- faceusp is     valid after call
  !   <- any other gp at points are now invalid
  !   <- cell ave gp is invalid since it might have been calculated using 
  !           clipped temperature
  !-----------------------------------------------------------------------------
  ! Notes:
  !  >  This is done because evaluating a bad pressure for limiting is invalid
  !       given the invalidity of the curve fits for temperature.
  !  >  Idea is for this to guarantee positivity and pressure limiting to be an 
  !       optional mechanism for increased dissipation.
  !  >  A modified form of internal energy is used for TPG, this version is
  !       shifted by the T=0 internal energy such that U_in positivity implies T
  !       positivity
  !-----------------------------------------------------------------------------
  !
  ! .. Use Statements  ..
  use geomvar_mod,   only : ncell,cell,face,n_solpts
  use thermo_mod,    only : update_mixture_density, compute_internal_energy
  use eqn_idx_mod,   only : ncs,nee,nmx,nmz,nec,nq
  use thermo_mod,    only : compute_minimum_internal_energy
  !
  ! .. Arguments       ..
  logical(lk), intent(in)                 :: lim_flxpts ! do we treat faces
  real(wp), intent(in)                    :: Ulim
  real(wp), dimension(:,:), intent(inout) :: usp
  real(wp), dimension(:,:), intent(in)    :: cell_ave_sol
  real(wp), dimension(:,:,:),intent(inout),optional   :: faceusp
  real(wp), dimension(:), intent(inout)   :: plt
  !
  ! .. Local Variables ..
  character(len=*), parameter :: pname = "limit_internal_energy"
  integer  :: ic,n,n1,n2,ncf,nf,sf,f1,f2,m,nface
  integer  :: numbad, irank, ierr
  logical(lk), allocatable :: bad_cel_loc(:)
  real(wp), dimension(:), allocatable :: e,eface,Uthermominimum,Utm_ave
  real(wp) :: aveU,minU,theta,Uminface_approx
  !
continue
  !
  e           = init_array(n_solpts,zero) ! F2003 AUTO-REALLOC
  eface       = init_array(10000,zero) ! F2003 AUTO-REALLOC
  Uthermominimum = init_array(n_solpts,zero) ! F2003 AR
  Utm_ave = init_array(ncell,zero) ! F2003 AR
  !
  !
  !
  ! Find minimum internal energy to help TPG limit properly
  if (ncs > 0) then
    call compute_minimum_internal_energy(1,n_solpts,usp,Uthermominimum)
    call compute_minimum_internal_energy(1,ncell,cell_ave_sol,Utm_ave) 
  else
    Uthermominimum(:) = zero
    Utm_ave(:) = zero
  end if
  !
  ! Calculate internal energy at all solution points
  call compute_internal_energy(1,n_solpts,usp,e)
  e(:) = (e(:) - Uthermominimum(:))*usp(:,nec)
  !
  ! Loop Over Elements
  do ic = 1,ncell
    n1 = cell(ic)%beg_sp ! beginning solution point index for this cell
    n2 = cell(ic)%end_sp ! ending    solution point index for this cell
    !
    !
    ! >> Find the minimum internal energy on points of interest <<
    !
    ! Find scaling factor based on mixture density
    aveU = cell_ave_sol(ic,nee) - &
            half*DOT_PRODUCT(cell_ave_sol(ic,nmx:nmz),cell_ave_sol(ic,nmx:nmz))/cell_ave_sol(ic,nec)
    aveU = aveU - Utm_ave(ic)*cell_ave_sol(ic,nec)
    minU = minval(e(n1:n2))
    !
    Uminface_approx = minval(Uthermominimum(n1:n2))
    !
    if (lim_flxpts .and. present(faceusp)) then
    do ncf = 1,size(cell(ic)%face)
      nf  = cell(ic)%face(ncf)%idx  ! index of the current face
      sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
      f1  = face(nf)%beg_fp
      f2  = face(nf)%end_fp
      nface = f2-f1+1
      ! over the face points
      call compute_internal_energy(1,nface,faceusp(f1:f2,:,sf),eface(1:nface))
      eface(1:nface) = (eface(1:nface) - Uminface_approx) * faceusp(f1:f2,nec,sf)
      minU  = min(minU ,minval(eface(1:nface)))
      !
    end do
    end if
    !
    !
    ! >> Limit Internal Energy <<
    !
    if (minU > Ulim) cycle
    if ((aveU-minU) < ten*rndoff) then
        theta = zero
    else
        !
        theta  = 0.99_wp*(aveU - Ulim) / (aveU - minU)  ! Z-S scaling factor
        theta = max(min(theta-eps2,one),zero)
        !
    end if
    !
    ! Update the output marker
    plt(ic) = min(plt(ic),theta)
    !
    !
    ! Scale the individual species
    do m=1,nq
      usp(n1:n2,m) = cell_ave_sol(ic,m) + &
                              theta*(usp(n1:n2,m)-cell_ave_sol(ic,m))
    end do
    !
    !
    if (lim_flxpts .and. present(faceusp)) then
      do ncf = 1,size(cell(ic)%face)
        nf  = cell(ic)%face(ncf)%idx  ! index of the current face
        sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
        f1  = face(nf)%beg_fp
        f2  = face(nf)%end_fp
        nface = f2-f1+1
        do m=1,nq
        faceusp(f1:f2,m,sf) = cell_ave_sol(ic,m) + &
                           theta*(faceusp(f1:f2,m,sf)-cell_ave_sol(ic,m))
        end do
        !
      end do
    end if
    !
  end do
  !end do outer_loop
  !
  deallocate(e)
  deallocate(Uthermominimum)
  deallocate(Utm_ave)
  deallocate(eface)
  !
end subroutine limit_internal_energy
!
!
!  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
! //////////////////////////////////////////////////////////////////////////// !
! //////////////////////////////////////////////////////////////////////////// !
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!
end module positivity_mod
