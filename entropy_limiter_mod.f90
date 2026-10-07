! =============================================================================
! >>> Limiting Method Using a Local Minimum Entropy Principal <<<
!
!  >> Procedure:
!     1 - Calculate minimum entropy bounds by taking the min on an element &
!         face neighbors
!     2 - Update the solution (No entropy-related source terms[reactions])
!     3 - Limit to enforce thermodynamic positivity (need valid rho, P)
!     4 - Limit to enforce the local entropy principal
!     5 - Add contribution from entropy-related source terms (chem. reactions)
!     6 - Additional positivity limiter for good luck
!
!  >> Input Options:
!     1 - use_entropy_limiter = T/F
!
!     2 - entropy_limiter%Option = 0 (scale to average) / 1 (Dynamic filter)
!         >> Need to have limiter filter set (i.e. LIMITER_FILTER_OPTION...)
!
!     3 - entropy_limiter%Interpolant = 0-4
!         >> How damping ratio is calculated, probably doesn't matter, use 0
!
!     4 - entropy_limiter%Tolerance = [real(wp)]
!         >> Relax the entropy limit to lower dissipation
!
!     5 - entropy_limiter%Check_Faces = T/F
!         >> Whether to include face points in damping ratio calculation
!
!     6 - entropy_limiter%Check_Nodes  = T/F  !!! NOT IMPLEMENTED > sim. to faces
!     7 - entropy_limiter%Check_Global = T/F  !!! NOT IMPLEMENTED > global min
!
!  << Output:  output_marker_variables = T
!     - outputs entropy limiter filter strength per stage
!     - outputs local minimum entropy limit after the last stage
!
! =============================================================================
module entropy_limiter_mod
  !.. Use Statements ..
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx_mod, only : ncs,npres,ntemp,ng
  !
  use ovar, only : force_nd_subroutines
  use ovar, only : entropy_limiter
  !
  implicit none
  !
  private
  !
  real(wp), save, public, allocatable, dimension(:) :: entropy_minimum
  !.. Module Arrays ..
  real(wp), save, allocatable, dimension(:) :: entropy_sp
  real(wp), save, allocatable, dimension(:) :: entropy_fp
  real(wp), save, allocatable, dimension(:) :: cell_ave_entropy
  !
  !.. Public Procedures ..
  public :: perform_entropy_limiting
  public :: update_parallel_entropy_bounds
  public :: calculate_entropy_bounds
  !
contains
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
! //////////////////////////////////////////////////////////////////////////// !
! //////////////////////////////////////////////////////////////////////////// !
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !
!
subroutine perform_entropy_limiting(nst)
  !> Carries out a limiting procedure based on a local minimum entropy
  !
  !> Zhang-Shu entropy limiting
  !> Dzanic, Trojak, Witherden stuff n things
  !
  !.. Use Statements ..
  use geomvar_mod,  only : ncell,cell,face,n_solpts
  use ovar,         only : rhoref, aref
  use flux_mod,     only : interp_cv_to_flxpts
  use cell_ave_mod, only : allocate_cell_ave, deallocate_cell_ave, &
                           compute_cell_ave_sol, cell_ave_sol, cell_ave_gp
  use limiters_mod,   only : apply_filter_limiter_to_usp
  use flowvar_mod,  only : scaled_marker, get_marker_index, nmarker
  use flowvar_mod,  only : usp,gpsp,faceusp,facegpsp
  use flowvar_mod,  only : entropy_limiter_strength
  use thermo_mod,   only : update_thermo_properties 
  use transport_mod,only : update_transport_properties
  use eqn_idx_mod,  only : ngamm,nmx,nmy,nec,nee,ncr
  use ovar,         only : gam, positivity
  !
  !.. Arguments ..
  integer,     intent(in)    :: nst
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "perform_entropy_limiting"
  !
  !.. Local Scalars ..
  integer       :: n,n1,n2,eq,ic,midx,ics,nn,imarker
  integer       :: nf,sf,f1,f2,ncf,nfp
  real(wp)      :: fil_ave_jump, sol_ave_jump
  real(wp)      :: s_ave, s_delta, denom
  real(wp)      :: rrat = one
  real(wp)      :: badness, point_entropy, filter_strength
  real(wp)      :: var_jump, bad_variable, badsp, sbeef
  real(wp)      :: entropy_filtered, fssp
  logical(lk)   :: bad_density, bad_pressure, bad_entropy, filter_type
  logical(lk)   :: bad_densisp
  !
  integer :: entropy_filter_type,ierr
  !
  !.. Local Arrays ..
  logical(ldk), allocatable :: bad_cell(:), troubled_cell(:)
  real(wp), allocatable     :: usp_filtered(:,:), gp_filtered(:) &
                              ,ufp_filtered(:,:,:), ufp_temp(:,:,:)
  !
  !
continue!
  !
  entropy_filter_type = entropy_limiter%Option
  !
  ! Make sure we immediately return without doing
  ! anything if entropy limiting is not turned on.
  if (entropy_filter_type == ENTLIM_OPT__OFF) return
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate the array storing the amount of limitng occuring
  !
  if (.not. allocated(entropy_limiter_strength)) then
    allocate ( entropy_limiter_strength(1:ncell,1:nmarker) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"entropy_limiter_strength",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !$ACC UPDATE DEVICE(entropy_limiter_strength)
  end if
  imarker = get_marker_index(nst)
  !
  midx = get_marker_index(nst)
  !
  bad_cell = init_array(ncell,fals) ! F2003 AUTO-REALLOC
  !
  !
  !
  !>>  Make sure needed arrays are allocated
  !
  call allocate_cell_ave()
  !
  if (.not. allocated(cell_ave_entropy)) then
    allocate(cell_ave_entropy(ncell), source=one)
    !$ACC UPDATE DEVICE(cell_ave_entropy)
  end if
  !
  if (.not. allocated(scaled_marker)) then
    allocate ( scaled_marker(1:ncell,1:nmarker) , source=one )
    !$ACC UPDATE DEVICE(scaled_marker)
  end if
  !
  if (positivity%Option == POSLIM_OPT__OFF) then
    ! If we are not using the positivity preservation limiter,
    ! we need to compute the cell average solution
    ! NOTE: This also updates cell_ave_gp
    call compute_cell_ave_sol()
  end if
  !
  ! > Get entropy at every solpt and facepoint
  ! entropycalculation
  ! SETH : MAKE A SUBROUTINE TO COMPUTE ENTROPY OVER A LOOP OF POINTS
  entropy_sp(1:n_solpts)= (gpsp(1:n_solpts,npres)*&
                          (usp(1:n_solpts,nec)**(-gpsp(1:n_solpts,ngamm))))
  !
  if (entropy_limiter%Check_Faces) then
    !SETH: CHANGE THIS TO JUST DO ONE BIG LOOP OVER ALL FACE POINTS
    !      INSTEAD OF OVER CELLS AND THEN OVER FACE POINTS OF EACH CELL
    do ic = 1,ncell
    do ncf = 1,size(cell(ic)%face)
      nf  = cell(ic)%face(ncf)%idx  ! index of the current face
      sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
      f1  = face(nf)%beg_fp
      f2  = face(nf)%end_fp
      ! over the face points
      entropy_fp(f1:f2)= (facegpsp(f1:f2,npres,sf)*&
                            (faceusp(f1:f2,nec,sf)**(-facegpsp(f1:f2,ngamm,sf))))
    end do
    end do
  end if
  !
  !
  ! -----------------------------------
  ! >>> LOOP TO CHECK FOR BAD CELLS <<<
  ! -----------------------------------
  cellloop : do ic = 1,ncell
    n1 = cell(ic)%beg_sp ! beginning solution point index for this cell
    n2 = cell(ic)%end_sp ! ending    solution point index for this cell
    !
    do n = n1,n2
      !
      bad_entropy = entropy_sp(n)    < entropy_minimum(ic)
      !
      ! Skip rest of cell if we know it needs limiting
      if (bad_entropy) then
        bad_cell(ic) = true
        cycle cellloop
      end if
    end do
    !
    if (entropy_limiter%Check_Faces) then
    do ncf = 1,size(cell(ic)%face)
      nf  = cell(ic)%face(ncf)%idx  ! index of the current face
      sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
      f1  = face(nf)%beg_fp
      f2  = face(nf)%end_fp
      ! over the face points
      ! SETH: PROBABLY DONT NEED THE "bad_entropy .or." PART SINCE WE ARE
      !       CYCLING TO THE NEW CELL AS SOON AS BAD_ENTROPY IS TRUE.
      bad_entropy = bad_entropy .or. ANY(entropy_fp(f1:f2) < entropy_minimum(ic))
      !
      ! Skip rest of cell if we know it needs limiting
      if (bad_entropy) then
        bad_cell(ic) = true
        cycle cellloop
      end if
    end do
    end if
    !
  end do cellloop
  !
  !
  ! ---------------------------------
  ! >>> GET CELL AVERAGED ENTROPY <<<
  ! ---------------------------------
  do ic=1,ncell
    if (.not. bad_cell(ic)) cycle
      ! entropycalculation
      cell_ave_entropy(ic) = (cell_ave_gp(ic,npres)*(cell_ave_sol(ic,nec)**(-cell_ave_gp(ic,ngamm))))
      !
  end do
  !
  ! -------------------------------
  ! >>> GET FILTERED QUANTITIES <<<
  ! -------------------------------
  if (entropy_filter_type == ENTLIM_OPT__ScaledFilter) then
    !
    ! Compute the reference filtered solution if using it
    usp_filtered = usp            ! F2003 Auto-Reallocation
    troubled_cell = bad_cell      ! F2003 Auto
    !
    allocate(gp_filtered(1:nq), source=zero)
    !
    ! get filtered solution in bad cells
    scaled_marker(:,midx) = one
    call apply_filter_limiter_to_usp(troubled_cell,usp_filtered)
    scaled_marker(:,midx) = zero
    !
    if (entropy_limiter%Check_Faces) then
      !
      !SETH: SHOULD PROBABLY ADD USP AND FACEUSP AS INPUTS TO
      !      INTERP_CV_TO_FLXPTS TO PREVENT HAVING TO SWAP ARRAYS LIKE THIS.
      ufp_filtered = faceusp        ! F2003 Auto-Reallocation
      ufp_temp     = faceusp        ! F2003 Auto-Reallocation
      ! Get facepoints from filtered solution
      call interp_cv_to_flxpts(999,nst,cv_only=true)
      !
      ufp_filtered = faceusp
      faceusp      = ufp_temp
      !
    end if
    !
  end if
  !
  !
  ! -------------------------------
  ! >>> LOOP TO LIMIT BAD CELLS <<<
  ! -------------------------------
  cells_of_grid: &
  do ic = 1,ncell
    if (.not.  bad_cell(ic)) cycle
    bad_cell(ic) = fals
    !
    n1 = cell(ic)%beg_sp ! beginning solution point index for this cell
    n2 = cell(ic)%end_sp ! ending    solution point index for this cell
    !
    ! >>> badness = point local damping ratio needed to meet limit
    !       closer to zero => more bad...
    !
    if (entropy_filter_type == ENTLIM_OPT__ScaledFilter) then
      filter_strength = zero   ! weakest = no filter applied
    else ! entropy_filter_type = ENTLIM_OPT__Flat
      filter_strength = one    ! weakest = no damping of H-O modes
    end if
    !
    !
    ! -------------------------------------------------------------------------
    ! >>>           Find filter strength for solution points                <<<
    ! -------------------------------------------------------------------------
    !
    points_of_cell_loop: &
    do n = n1,n2
      !
      sol_ave_jump = eps8
      bad_variable = 1
      badness = one
      !
      bad_entropy = entropy_sp(n) < entropy_minimum(ic)
      !
      if (.not. bad_entropy) then
        cycle  ! cell is bad but not this particular point
      end if
      !
      ! ------------------------------------------------
      ! >>> Calculate Damping Ratio (BADNESS) NEEDED <<<
      ! ------------------------------------------------
      !
      ! >>> Entropy <<<
      if (bad_entropy) then
        !
        s_ave    = cell_ave_entropy(ic)
        var_jump = s_ave - entropy_sp(n)
        !
        !
        select case (entropy_limiter%Interpolant)
          case (ENTLIM_FUN__Linear)
            ! standard linear interpolant
            s_delta  = (s_ave - entropy_minimum(ic))
            badsp = s_delta/(entropy_sp(n)-entropy_minimum(ic))
            !
          case (ENTLIM_FUN__DensityWeighted)
            ! Density weighted interpolant
            s_delta  = cell_ave_sol(ic,nec)*(s_ave - entropy_minimum(ic))
            badsp = s_delta/(s_delta-usp(n,nec)*(entropy_sp(n)-entropy_minimum(ic)))
            !
          case (ENTLIM_FUN__LogScaled)
            ! Totally obvious formula which adds some extra concavity to the interpolant fxn 
            !     NOTE: these equations probably merit some kind of explanation either here or in docs
            s_delta  = log(s_ave/entropy_minimum(ic))
            badsp = s_delta/(log(s_ave/entropy_sp(n)))
            !
          case (ENTLIM_FUN__LogScaledPade)
            ! Pade approx
            s_delta =          (s_ave - entropy_minimum(ic)) / max(eps12, s_ave - entropy_sp(n))
            badsp   = s_delta *(s_ave + entropy_sp(n))        / (s_ave + entropy_minimum(ic)) 
            badsp = min(max(badsp,zero),one)
            !
          case default ! (ENTLIM_FUN__LinearLvIhme)
            ! entropy bounded dg scheme for euler eqns - Lv, Ihme - plus some algebra
            s_delta =          (entropy_sp(n) - entropy_minimum(ic))
            rrat    = (cell_ave_sol(ic,nec) / usp(n,nec))**gpsp(n,ngamm)
            denom   = (s_delta - rrat*(s_ave-entropy_minimum(ic)))
            denom   = denom + sign(eps12, denom)
            badsp = s_delta / denom 
            !
        end select
        !
        badness = min(max(badsp,zero),one)
        !
      end if
      !
      !
      ! ----------------------------------------------------
      ! >>> CALCULATE FILTER STRENGTH FROM DAMPING RATIO <<<
      ! ----------------------------------------------------
      ! > Now we have the badness of the point and it's delta 
      !     from the solution average. Need to translate that from the
      !     Zhang & Shu theta (badness), to a filter strength
      if (entropy_filter_type == ENTLIM_OPT__Flat) then
        ! If we're doing scaling that is equiv to zhang and shu
        filter_strength = min(filter_strength,badness)
        !
      else if (entropy_filter_type == ENTLIM_OPT__ScaledFilter) then
        !
        sol_ave_jump = entropy_sp(n) - cell_ave_entropy(ic)
        !
        ! Need to get the jumps between the filtered solution and cell average
        call update_thermo_properties(usp_filtered(n,:), gp_filtered, pname) 
        call calculate_entropy_1p(usp_filtered(n,:), gp_filtered, entropy_filtered)
        fil_ave_jump = entropy_filtered - cell_ave_entropy(ic)
        
        ! Now to calculate filter strength needed to filter the solution
        denom = (sol_ave_jump - fil_ave_jump)
        denom = denom + sign(eps12, denom)
        fssp = (one-badness)*sol_ave_jump / denom 
        !
        !
        ! If the filter strength is outside of bounds, use 2nd interpolation
        !   betweena filter and cell average
        if (fssp > one) then
          denom = fil_ave_jump
          denom = (sol_ave_jump - fil_ave_jump)
          denom = denom + sign(eps12, denom)
          fssp = (two*fil_ave_jump - badness*sol_ave_jump ) / denom
        end if
        !
        ! filter strength should be between zero and two
        ! fssp = max(fssp, zero)
        filter_strength = max(fssp,  filter_strength)
        !
      end if ! scaled filter
      !
      !!
    end do points_of_cell_loop
    !
    ! -------------------------------------------------------------------------
    ! >>>           Find filter strength for face points                    <<<
    ! -------------------------------------------------------------------------
    !
    if (entropy_limiter%Check_Faces) then
    points_of_cell_face_loop: &
    do ncf = 1,size(cell(ic)%face)
      nf  = cell(ic)%face(ncf)%idx  ! index of the current face
      sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
      f1  = face(nf)%beg_fp
      f2  = face(nf)%end_fp
      do nfp = f1,f2
        !
        sol_ave_jump = eps8
        bad_variable = 1
        badness = one
        !
        bad_entropy = entropy_fp(nfp) < entropy_minimum(ic)
        !
        if (.not. bad_entropy) then
          cycle  ! cell is bad but not this particular point
        end if
        !
        ! ------------------------------------------------
        ! >>> Calculate Damping Ratio (BADNESS) NEEDED <<<
        ! ------------------------------------------------
        !
        ! >>> Entropy <<<
        !
        s_ave    = cell_ave_entropy(ic)
        var_jump = s_ave - entropy_fp(nfp)
        !
        !
        select case (entropy_limiter%Interpolant)
          case (ENTLIM_FUN__Linear)
            ! standard linear interpolant
            s_delta  = (s_ave - entropy_minimum(ic))
            badsp = s_delta/(entropy_fp(nfp)-entropy_minimum(ic))
            !
          case (ENTLIM_FUN__DensityWeighted)
            ! Density weighted interpolant
            s_delta  = cell_ave_sol(ic,nec)*(s_ave - entropy_minimum(ic))
            badsp = s_delta/(s_delta-faceusp(nfp,nec,sf)*(entropy_fp(nfp)-entropy_minimum(ic)))
            !
          case (ENTLIM_FUN__LogScaled)
            ! Totally obvious formula which adds some extra concavity to the interpolant fxn 
            !     NOTE: these equations probably merit some kind of explanation either here or in docs
            s_delta  = log(s_ave/entropy_minimum(ic))
            ! SETH: ADD PROTECTION ON LOG IN DENOMINATOR
            badsp = s_delta/(log(s_ave/entropy_fp(nfp)))
            !
          case (ENTLIM_FUN__LogScaledPade)
            ! Pade approx (log scaled using an approximation for the log function)
            s_delta =          (s_ave - entropy_minimum(ic)) / max(eps12, s_ave - entropy_fp(nfp))
            badsp   = s_delta *(s_ave + entropy_fp(nfp))        / (s_ave + entropy_minimum(ic)) 
            badsp = min(max(badsp,zero),one)
            !
          case default ! (ENTLIM_FUN__LinearLvIhme)
            ! entropy bounded dg scheme for euler eqns - Lv, Ihme - plus some algebra
            ! if rrat == 1, it should reproduce ENTLIM_FUN__Linear
            ! if rrat <  1, it will add more dissipation
            s_delta =          (entropy_fp(nfp) - entropy_minimum(ic))
            rrat    = (cell_ave_sol(ic,nec) / faceusp(nfp,nec,sf))**facegpsp(nfp,ngamm,sf)
            denom   = (s_delta - rrat*(s_ave-entropy_minimum(ic)))
            denom   = denom + sign(eps12, denom)
            badsp = s_delta / denom 
            !
        end select
        !
        badness = min(max(badsp,zero),one)
        !
        !
        !
        ! ----------------------------------------------------
        ! >>> CALCULATE FILTER STRENGTH FROM DAMPING RATIO <<<
        ! ----------------------------------------------------
        ! > Now we have the badness of the point and it's delta 
        !     from the solution average. Need to translate that from the
        !     Zhang & Shu theta (badness), to a filter strength
        if (entropy_filter_type == ENTLIM_OPT__Flat) then
          ! If we're doing scaling that is equiv to zhang and shu
          filter_strength = min(filter_strength,badness)
          !
        else if (entropy_filter_type == ENTLIM_OPT__ScaledFilter) then
          !
          sol_ave_jump = entropy_fp(nfp) - cell_ave_entropy(ic)
          !
          ! Need to get the jumps between the filtered solution and cell average
          call update_thermo_properties(ufp_filtered(nfp,:,sf), gp_filtered, pname) 
          call calculate_entropy_1p(ufp_filtered(nfp,:,sf), gp_filtered, entropy_filtered)
          fil_ave_jump = entropy_filtered - cell_ave_entropy(ic)
          
          ! Now to calculate filter strength needed to filter the solution
          denom = (sol_ave_jump - fil_ave_jump)
          denom = denom + sign(eps12, denom)
          fssp = (one-badness)*sol_ave_jump / denom 
          !
          !
          ! If the filter strength is outside of bounds, use 2nd interpolation
          !   betweena filter and cell average
          if (fssp > one) then
            denom = fil_ave_jump
            denom = (sol_ave_jump - fil_ave_jump)
            denom = denom + sign(eps12, denom)
            fssp = (two*fil_ave_jump - badness*sol_ave_jump ) / denom
          end if
          !
          ! filter strength should be between zero and two
          ! fssp = max(fssp, zero)
          filter_strength = max(fssp,  filter_strength)
          !
        end if ! scaled filter
        !
        !
      end do
    end do points_of_cell_face_loop
    end if
    !
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    !
    !
    ! -------------------------------------------
    ! >>> APPLY FILTER WITH COMPUTED STRENGTH <<<
    ! -------------------------------------------
    !
    if (entropy_filter_type == ENTLIM_OPT__Flat) then
      filter_strength = min(max(zero,filter_strength),one)
      !
      !write(iout,*) "fs: ", filter_strength
      do n=n1,n2
        usp(n,:) = filter_strength * (usp(n,:) - cell_ave_sol(ic,:))&
                                               + cell_ave_sol(ic,:)
      end do
      !
      entropy_limiter_strength(ic,imarker) = one - filter_strength
      !
    else ! filter == ENTLIM_OPT__ScaledFilter
      filter_strength = max(min(two,filter_strength),zero)
      !
      if (filter_strength <= one) then
        ! 0 < fs < 1
        do n=n1,n2
          usp(n,:) = (one-filter_strength)*usp(n,:) + filter_strength*usp_filtered(n,:)
        end do
        !
      else   ! 1 < fs < 2
        do n=n1,n2
          usp(n,:) = (two-filter_strength)* & 
                (usp_filtered(n,:) - cell_ave_sol(ic,:)) + cell_ave_sol(ic,:)
        end do
      end if
      !
      entropy_limiter_strength(ic,imarker) = filter_strength
      !
    end if ! Filter IF
    !
    !
    ! Update solution point entropy for bad cells
    call update_thermo_properties(n1,n2,usp,gpsp)
    !
    ! SETH : MAKE A SUBROUTINE TO COMPUTE ENTROPY OVER A LOOP OF POINTS
    ! entropycalculation
    entropy_sp(n1:n2)= (gpsp(n1:n2,npres)*(usp(n1:n2,nec)**(-gpsp(n1:n2,ngamm))))
    !
    if (entropy_limiter%Check_Faces) then
    do ncf = 1,size(cell(ic)%face)
      nf  = cell(ic)%face(ncf)%idx  ! index of the current face
      sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
      f1  = face(nf)%beg_fp
      f2  = face(nf)%end_fp
      ! over the face points
      entropy_fp(f1:f2)= (facegpsp(f1:f2,npres,sf)*&
                            (faceusp(f1:f2,nec,sf)**(-facegpsp(f1:f2,ngamm,sf))))
    end do
    end if
    !
    
  end do cells_of_grid
  !
  if (entropy_limiter%Check_Faces) then
    call interp_cv_to_flxpts(999,nst,cv_only=true)
    do ic = 1,ncell
      if (.not.  bad_cell(ic)) cycle
      do ncf = 1,size(cell(ic)%face)
        nf  = cell(ic)%face(ncf)%idx  ! index of the current face
        sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
        f1  = face(nf)%beg_fp
        f2  = face(nf)%end_fp
        ! over the face points
        entropy_fp(f1:f2)= (facegpsp(f1:f2,npres,sf)*&
                              (faceusp(f1:f2,nec,sf)**(-facegpsp(f1:f2,ngamm,sf))))
      end do
    end do
  end if
  ! SETH: NEED TO RECOMPUTE ENTROPY_SP AND ENTROPY_FP AT THE END OF STABILIZE_USP
  !       JUST ADD THE IS_CURRENT CRAP
  !
  call deallocate_cell_ave()
  deallocate(cell_ave_entropy)
  if (allocated(gp_filtered)) deallocate(gp_filtered)
  if (allocated(usp_filtered)) deallocate(usp_filtered)
  if (allocated(bad_cell)) deallocate ( bad_cell )
  if (allocated(troubled_cell)) deallocate ( troubled_cell )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine perform_entropy_limiting
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine calculate_entropy_bounds(is_first_iter)
  ! Find the minimum entropy on an element and it's face neighbors
  !
  ! -Comunicate the cell minimum to MPI partition neighbors
  ! -Find the neighborhood minimum amonng MPI Rank partners
  !
  !.. Use Statements ..
  use geomvar_mod,  only : ncell,cell,face,n_solpts,bface,n_flxpts
  use flowvar_mod,  only : gpsp, usp, faceusp, facegpsp
  use eqn_idx_mod,  only : ngamm
  use ovar,         only : gam
  use parallel_mod, only : entropy_comm
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "calculate_entropy_bounds"
  logical, intent(in)     :: is_first_iter
  !
  !.. Local Scalars ..
  real(wp)      :: point_entropy, stol
  logical(lk)   :: update_min = true
  integer       :: ic,lif,gif,icfn,side,n1,n2,n,nf,f1,f2
  integer       :: ncf,sf
  !.. Local Arrays ..
  real(wp), allocatable      :: smin_tmp(:)
  !
  continue
  !
  call debug_timer(entering_procedure,pname)
  !
  stol = entropy_limiter%Tolerance
  !
  !
  if (.not. allocated(entropy_minimum)) then
    allocate(entropy_minimum(ncell), source=1.0E99_wp)
  end if
  !
  if (.not. allocated(entropy_sp)) then
    allocate(entropy_sp(n_solpts), source=1.0E99_wp)
  end if
  !
  if (.not. allocated(entropy_fp)) then
    allocate(entropy_fp(n_flxpts), source=1.0E99_wp)
  end if
  !
  if (.not. allocated(entropy_comm)) then
    allocate(entropy_comm(n_flxpts,2), source=1.0E99_wp)
  end if
  !
  ! ===============================================
  ! >>> Calculate the Minimum Entropy on a Cell <<<
  ! ===============================================
  ! NOTE: Can optimize by saving entropy calculated during filter application
  !
  ! > Reset at each evaulation of fxn
  entropy_minimum(:) = 1.0E99_wp
  !
  if (ncs == 0) then
    ! regular plain old air, no funny business
    gpsp(1:n_solpts,ngamm) = gam
  end if
    !
  do ic =1,ncell
    n1 = cell(ic)%beg_sp ! beginning solution point index for this cell
    n2 = cell(ic)%end_sp ! ending    solution point index for this cell
    !
    do n=n1,n2
      !
      if (is_first_iter) then 
        ! entropycalculation
        entropy_sp(n)= (gpsp(n,npres)*(usp(n,nec)**(-gpsp(n,ngamm))))
      end if
      !
      entropy_minimum(ic) = min(entropy_sp(n), entropy_minimum(ic))
      !
    end do
    !
    !
    ! >> get minimum from face points too
    ! ..... disabled because this kills the dissipation from the method
!!!    if (entropy_limiter%Check_Faces) then
!!!    do ncf = 1,size(cell(ic)%face)
!!!      nf  = cell(ic)%face(ncf)%idx  ! index of the current face
!!!      sf  = cell(ic)%face(ncf)%side ! side of face that cell nc is on
!!!      f1  = face(nf)%beg_fp
!!!      f2  = face(nf)%end_fp
!!!      ! over the face points
!!!      if (is_first_iter) then 
!!!        entropy_fp(f1:f2)= (facegpsp(f1:f2,npres,sf)*&
!!!                         (faceusp(f1:f2,nec,sf)**(-facegpsp(f1:f2,ngamm,sf))))
!!!      end if
!!!      !
!!!!!!      entropy_minimum(ic) = min(minval(entropy_fp(f1:f2)), entropy_minimum(ic))
!!!      !
!!!    end do
!!!    end if
    !
  end do
  !
  ! > Step 1.5: Set up the entropy comms to be shared with face neighbors
  ! entropy_comm(1:nfp*nf,1:2) = entropy minimum of the cell on sides of a face
  !
  ! > Loop parallel boundary faces
  do nf = 1,size(bface,dim=2)
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    ! This is a communication / boundary element
    ! Set the appropriate entry in the comm array
    !
    f1  = face(nf)%beg_fp
    f2  = face(nf)%end_fp
    !
    entropy_comm(f1:f2,1) = entropy_minimum(bface(2,nf))
    !
    !
  end do
  !
  !
  ! > Step 2: Compare face neighbors and update rank local minimum bounds
  smin_tmp = entropy_minimum ! Auto
  !
  do ic = 1,ncell
    ! > Get face neighbors
    !smin_tmp(ic) = entropy_minimum(ic)
    !
    !
    do lif = 1, size(cell(ic)%face(:))
      ! lif = local index of face
      ! gif = global index of face
      !
      side = cell(ic)%face(lif)%side
      gif  = cell(ic)%face(lif)%idx
      !
      ! icfn = index of cell face neighbor
      if (side == 1) then ! we are left,  look at right
        icfn = face(gif)%right%cell
      else                ! we are right, look at left
        icfn = face(gif)%left%cell
      end if
      !
      if (icfn > ncell) then
        cycle
      end if
      !
      smin_tmp(ic) = min(smin_tmp(ic), entropy_minimum(icfn))
      !
    end do
    !
  end do
  !
  ! First try exp(Entropy) - stol    (good tolerancing for low entropy)
  entropy_minimum(:) = max(eps12,smin_tmp(:) - stol)   
  ! Then try exp(Entropy - stop)     (tolerancing for high entropy)
  entropy_minimum(:) = min(entropy_minimum(:),(smin_tmp(:)*EXP(-one+stol)))
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine calculate_entropy_bounds
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine update_parallel_entropy_bounds
  !
  ! communication happens during exchange faceusp
  !
  ! After communication of entropy is done, update minimums for face neighbors
  ! Also check for physical boundary conditions which will lower the entropy 
  !                 (isothermal cold walls)
  !
  use geomvar_mod,  only : bface,face,n_phys_bnd_flxpts
  use parallel_mod, only : entropy_comm
  use eqn_idx_mod,  only : nec,ngamm,npres
  use flowvar_mod,  only : facegpsp,faceusp
  !
  ! .. Variable Declarations ..
  !
  integer  :: nf,f1,f2,ie 
  real(wp), dimension(:), allocatable :: p,rho,g
  real(wp) :: emin0
  !
  ! bc_comm = + periodic
  !
continue
  ! > Loop parallel boundary faces
  do nf = 1,size(bface,dim=2)
    if (bface(1,nf) == bc_cpu_bnd) then
    ! This is a communication / boundary element
    ! Set the appropriate entry in the comm array
    !
    f1  = face(nf)%beg_fp
    !
    ie = bface(2,nf)
    emin0 = entropy_minimum(ie)
    entropy_minimum(ie) = min( entropy_comm(f1,2), emin0)
    !
    else if (face(nf)%beg_fp < n_phys_bnd_flxpts) then
      ! This face is a phycal bondary and we need to make it see the entropy there
      f1 = face(nf)%beg_fp
      f2 = face(nf)%end_fp
      !
      p   = facegpsp(f1:f2,npres,2)
      g   = facegpsp(f1:f2,ngamm,2)
      rho = faceusp(f1:f2,nec,2)
      !
      ie = bface(2,nf)
      emin0 = entropy_minimum(ie)
      entropy_minimum(ie) = min(emin0, minval(p*rho**(-g)))
    end if
    !
  end do
  !
end subroutine update_parallel_entropy_bounds
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine calculate_entropy_1p(cv,gp,entropy)
  ! Calculates entropy at a single point with given cv
  !
  use eqn_idx_mod,  only : ncs,nec,ncr,ngamm
  use thermo_mod,   only : rspec,wtmol
  use thermo_mod,   only : compute_mixture_rgas
  use ovar,         only : gam
  !
  ! .. Arguments .. 
  real(wp)  :: cv(:), gp(:), entropy
  !
  ! .. Locals ..
  integer   :: ics 
  real(wp)  :: wtmol_mix_inv, rgas_mix
  real(wp), allocatable  :: ysp(:), xsp(:)
  real(wp), allocatable  :: s_sp(:,:)
  !
  continue
  !
  if (ncs == 0) then
    ! regular plain old air, no funny business
    gp(ngamm) = gam
  end if
  !
  ! entropycalculation
  entropy = (gp(npres)*(cv(nec)**(-gp(ngamm))))
  !
end subroutine calculate_entropy_1p
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
end module entropy_limiter_mod
