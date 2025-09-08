module ignition_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use eqn_idx_mod, only : nec,nee,npres,ntemp,ncp
  !
  implicit none
  !
  private
  !
  public :: compute_ignition_source_term
  !
  ! Ignition variables
  !
  real(wp), save, allocatable :: ign_src_loc(:,:)
  real(wp), save              :: ign_src_radius
  !
  logical(lk), save :: ignition_is_on = fals
  logical(lk), save :: ignition_reached_goal = fals
  logical(lk), save :: ignition_reached_temp_goal = fals
  logical(lk), save :: ignition_reached_L2_goal = fals
  logical(lk), save :: ignition_reached_max_temp = fals
  integer,     save :: ignition_iter_goal = -1
  !
  integer, save :: ignition_first_timestep = huge(0)
  integer, save :: ignition_final_timestep = huge(0)
  !
  logical(lk), save :: ignition_needs_initializing = true
  !
  real(wp), save :: ign_max_temp  = zero
  real(wp), save :: ign_src_mag   = zero
  real(wp), save :: ign_spark_cap = zero
  real(wp), save :: ign_p2r_L2_goal = -one
  !
contains
!
!###############################################################################
!
subroutine compute_ignition_source_term(n1,n2,cv,gp,residual,react_temp, &
                                        reactant_species_list, &
                                        product_species_list)
  !
  !.. Use Statements ..
  use ovar, only : ignition_method
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: residual(:,:)
  real(wp), intent(inout) :: react_temp(:)
  integer,     intent(in) :: reactant_species_list(:)
  integer,     intent(in) :: product_species_list(:)
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_ignition_source_term"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call check_for_ignition()
  !
  ! Copy the solution point temperature into the reaction temperature
  !
  do n = n1,n2
    react_temp(n) = gp(n,ntemp)
  end do
  !
  ! Compute the ignition source term or reaction temperature
  !
  if (ignition_is_on) then
    !
    if (ignition_method == IGN_DELTA_TEMP) then
      call compute_ignition_source_ncc(n1,n2,cv,gp,residual)
    else if (ignition_method == IGN_EREF_SPARK) then
      call compute_ignition_source_eref(n1,n2,gp,residual)
    else if (ignition_method == IGN_REACT_TEMP) then
      call compute_ignition_source_reactmb(n1,n2,gp,react_temp)
    end if
    !
    ! Check to see if any of the ignition success criteria have been met
    !
    call check_ignition_status(reactant_species_list,product_species_list)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_ignition_source_term
!
!###############################################################################
!
subroutine check_for_ignition()
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : itcur,grid_scaling_factor
  !
  !.. Local Scalars ..
  integer :: i
  !
continue
  !
  !
  if (ignition_needs_initializing) then
    ignition_needs_initializing = fals
    call initialize_ignition
  end if
  !
  if (itcur <= ignition_final_timestep) then
    if (itcur >= ignition_first_timestep) then
      !
      if (itcur == ignition_first_timestep .and. .not.ignition_is_on) then
        if (mypnum == glb_root) then
          write (iout,*)
          write (iout,2) 1,(ign_src_loc(i,1)/grid_scaling_factor, i=1,nr)
          write (iout,2) 2,(ign_src_loc(i,2)/grid_scaling_factor, i=1,nr)
          write (iout,3) ign_src_radius/grid_scaling_factor
          write (iout,*)
          write (iout,1) itcur
          flush (iout)
        end if
      end if
      !
      ignition_is_on = true
      !
    else
      !
      ignition_is_on = fals
      !
    end if
  end if
  !
  ! Format Statements
  !
  1 format (//,2x,80("#"),//, &
               6x,"STARTING IGNITION SOURCE TERM AT TIMESTEP # ",i0,//, &
               2x,80("#"),//)
  2 format ("Ignition source location ",i1," = ",3(f8.4,:,", "))
  3 format ("Ignition source radius          = ",f8.4)
  !
end subroutine check_for_ignition
!
!###############################################################################
!
subroutine initialize_ignition()
  !
  !.. Use Statements ..
  use ovar,        only : itcur,cv_ref,rhoref,aref,tref
  use ovar,        only : num_rk_stages
  use ovar,        only : ignition_method
  use ovar,        only : ignition_start_timestep
  use ovar,        only : ignition_prod2reac_L2_ratio_goal
  use ovar,        only : ignition_max_temperature
  use ovar,        only : ignition_spark_cap
  use ovar,        only : ignition_source_magnitude
  use ovar,        only : ignition_source_radius
  use ovar,        only : ignition_source_location
  use ovar,        only : ignition_duration
  use ovar,        only : itestcase,grid_scaling_factor
  use eqn_idx_mod, only : NCD_Axial,NCD_Vert,NCD_Span
  use geomvar_mod, only : nr
  !
  !.. Local Scalars ..
  real(wp) :: r,x,theta
  !
  logical(lk) :: disable_ignition
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_ignition"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Set the first ignition time step to the max between the current
  ! time step and the input ignition_start_timestep
  !
  ignition_first_timestep = max(itcur,ignition_start_timestep)
  !
  ! Initialize the final ignition time step to a large number
  !
  ignition_final_timestep = huge(0)
  !
  ! Check to see if ignition is disabled
  !
  disable_ignition = fals
  !
  if (ignition_method == IGN_OFF) then
    !
    ! If the ignition method is simply off, dont write anything
    ! about ignition being disabled to output
    !
    disable_ignition = true
    !
  else
    !
    if (ignition_source_radius == zero) then
      if (mypnum == glb_root) write (iout,1)
      if (mypnum == glb_root) write (iout,2)
      disable_ignition = true
    end if
    !
    if (ignition_duration <= 0) then
      if (.not. disable_ignition) then
        if (mypnum == glb_root) write (iout,1)
      end if
      if (mypnum == glb_root) write (iout,3)
      disable_ignition = true
    end if
    !
    if (ignition_source_magnitude == zero) then
      if (.not. disable_ignition) then
        if (mypnum == glb_root) write (iout,1)
      end if
      if (mypnum == glb_root) write (iout,4)
      disable_ignition = true
    end if
    !
    if (ignition_max_temperature <= zero .and. &
        ignition_prod2reac_L2_ratio_goal <= zero) then
      if (.not. disable_ignition) then
        if (mypnum == glb_root) write (iout,1)
      end if
      if (mypnum == glb_root) write (iout,5)
      disable_ignition = true
    end if
    !
    if (mypnum == glb_root) then
      write (iout,6)
      flush (iout)
    end if
    !
  end if
  !
  ! If ignition is disabled, set the ignition final time step to 0
  ! so that the subroutine check_for_ignition never turns ignition on.
  !
  if (disable_ignition) then
    ignition_final_timestep = 0
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  ! Initialize the ignition input values
  !
  ignition_reached_max_temp = fals
  !
  if (ignition_method == IGN_DELTA_TEMP) then
    ign_src_mag = ignition_source_magnitude / tref / real(num_rk_stages,kind=wp)
  else if (ignition_method == IGN_EREF_SPARK) then
    ign_src_mag = abs(ignition_source_magnitude)*cv_ref(nee)
  else if (ignition_method == IGN_REACT_TEMP) then
    ign_src_mag = ignition_source_magnitude / tref
  else if (ignition_method == IGN_REACT_PRES) then
    ign_src_mag = ignition_source_magnitude / (rhoref*aref*aref)
  end if
  !
  ign_max_temp = ignition_max_temperature / tref
  ign_spark_cap = ignition_spark_cap / (rhoref*aref*aref)
  !
  ign_p2r_L2_goal = ignition_prod2reac_L2_ratio_goal
  !
  ! Scale the ignition source location and radius from grid units 
  !
  if (abs(itestcase) == UIUC_Scramjet) then
    !
    ign_src_loc = init_array(nr,2,zero) ! F2003 AUTO-REALLOC
    !
    x = ignition_source_location(NCD_Axial)*grid_scaling_factor
    r = ignition_source_location(NCD_Vert)*grid_scaling_factor
    !
    if (x < zero) x = -one/eps8
    !
    theta = 22.5_wp*degrees2radians
    ign_src_loc(NCD_Axial,1) = x
    ign_src_loc(NCD_Vert ,1) = r*cos(theta)
    ign_src_loc(NCD_Span ,1) = r*sin(theta)
    !
    theta = 67.5_wp*degrees2radians
    ign_src_loc(NCD_Axial,2) = x
    ign_src_loc(NCD_Vert ,2) = r*cos(theta)
    ign_src_loc(NCD_Span ,2) = r*sin(theta)
    !
  else
    !
    ign_src_loc = init_array(nr,1,zero) ! F2003 AUTO-REALLOC
    !
    ign_src_loc(1:nr,1) = ignition_source_location(1:nr)*grid_scaling_factor
    !
  end if
  !
  ign_src_radius = ignition_source_radius*grid_scaling_factor
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (//,2x,80("#"),//, &
               6x,"DISABLING THE IGNITION SOURCE TERM BECAUSE:")
  2 format (8x,"THE IGNITION SOURCE TERM RADIUS IS ZERO!!")
  3 format (8x,"THE IGNITION SOURCE TERM DURATION IS ZERO!!")
  4 format (8x,"THE IGNITION SOURCE TERM MAGNITUDE IS ZERO!!")
  5 format (8x,"THE IGNITION SOURCE TERM TARGET GOALS ARE BOTH <= ZERO!!",/, &
            12x,"1) MAX TARGET TEMPERATURE, AND",/, &
            12x,"2) TARGET RATIO OF PRODUCTS-TO-REACTANTS L2 NORMS",/, &
            8x,"AN IGNITION TARGET IS NEEDED FOR THE IGNITION SOURCE TERMS!!")
  6 format (/,2x,80("#"),//)
  !
end subroutine initialize_ignition
!
!###############################################################################
!
subroutine check_ignition_status(reactant_species_list,product_species_list)
  !
  !.. Use Statements ..
  use ovar, only : uverr,itcur
  use ovar, only : ignition_duration
  !
  !.. Formal Arguments ..
  integer, intent(in) :: reactant_species_list(:)
  integer, intent(in) :: product_species_list(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: p2r,products_L2_sum,reactants_L2_sum
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_ignition_status"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Check if the ratio of [the sum of the product species L2 norms] to
  ! [the sum of the reactant species L2 norms] has met the input goal.
  !
  if (.not. ignition_reached_L2_goal) then
    if (ign_p2r_L2_goal > zero) then
      if (allocated(uverr)) then
        !
        ! Compute the ratio of [the sum of the product species L2 norms] to
        ! [the sum of the reactant species L2 norms]
        !
        products_L2_sum  = sum( uverr( product_species_list(:),2,2) )
        reactants_L2_sum = sum( uverr(reactant_species_list(:),2,2) )
        !
        if (reactants_L2_sum /= zero) then
          p2r = products_L2_sum / reactants_L2_sum
        else
          p2r = zero
        end if
        !
        ! If the above ratio is greater than the value of the input goal, we
        ! will assume that ignition was successful and reactions are occurring.
        ! Set the value of the input goal to a negative number so that ignition
        ! source terms will no longer be computed.
        !
        if (p2r > ign_p2r_L2_goal) then
          ign_p2r_L2_goal = -one
          ignition_reached_goal = true
          ignition_reached_L2_goal = true
          if (mypnum == glb_root) then
            write (iout,1) itcur
          end if
        end if
        !
      end if
    end if
  end if
  !
  ! We only need to perform an allreduce on ignition_reached_max_temp
  ! if the max temperature hasnt been reached yet.  The saved variable
  ! ignition_iter_goal will retain its initial value of -1 until the
  ! maximum temperature is reached, after which it will be set to the
  ! current iteration number which will be > 0.
  !
  if (.not. ignition_reached_temp_goal) then
    !
    if (ncpu > 1) then
      n = merge(1,0,ignition_reached_max_temp)
      call mpi_allreduce(MPI_IN_PLACE,n,1_int_mpi,mpi_inttyp, &
                         MPI_MAX,MPI_COMM_WORLD,mpierr)
      if (n == 1) then
        ignition_reached_max_temp = true
        ignition_reached_goal = true
        ignition_reached_temp_goal = true
        if (mypnum == glb_root) then
          write (iout,2) itcur
        end if
      end if
    else
      if (ignition_reached_max_temp) then
        ignition_reached_goal = true
        ignition_reached_temp_goal = true
      end if
    end if
    !
  end if
  !
  !
  !
  if (ignition_reached_goal) then
    !
    if (ignition_iter_goal < 0) then
      ignition_iter_goal = itcur
      ignition_final_timestep = itcur + ignition_duration
      if (mypnum == glb_root) then
        write (iout,3) itcur,ignition_final_timestep
      end if
    end if
    !
    if (itcur == ignition_final_timestep) then
      if (mypnum == glb_root) then
        write (iout,4) itcur
      end if
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (//,2x,80("#"),//, &
               6x,"REACHED IGNITION GOAL FOR RATIO OF L2 NORMS OF",/, &
               6x,"PRODUCTS-TO-REACTANTS AT TIMESTEP # ",i0,//, &
               2x,80("#"),//)
  2 format (//,2x,80("#"),//, &
               6x,"REACHED IGNITION MAX TEMPERATURE AT TIMESTEP # ",i0,//, &
               2x,80("#"),//)
  3 format (//,2x,80("#"),//, &
               6x,"REACHED IGNITION GOAL AT TIMESTEP # ",i0,/, &
               6x,"THE FINAL IGNITION TIME STEP WILL BE # ",i0,/, &
               2x,80("#"),//)
  4 format (//,2x,80("#"),//, &
               6x,"TURNING OFF IGNITION SOURCE TERM AT TIMESTEP # ",i0,//, &
               2x,80("#"),//)
  !
end subroutine check_ignition_status
!
!###############################################################################
!
subroutine compute_ignition_source_ncc(n1,n2,cv,gp,residual)
  !
  !======================================================================
  !
  !   Ignition method based on adding a source term to the energy
  !   equation computed as rho*cp*(delta_Temperature).
  !
  !======================================================================
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr,xyz
  use flowvar_mod, only : dtsp
  use metrics_mod, only : jac
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: residual(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: r,spark,energy_source
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_ignition_source_ncc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  !
  !
  do n = n1,n2
    !
    r = norm2(xyz(1:nr,n)-ign_src_loc(1:nr,1))
    do i = 2,size(ign_src_loc,dim=2)
      r = min( r , norm2(xyz(1:nr,n)-ign_src_loc(1:nr,i)) )
    end do
    !
    energy_source = zero
    !
    if (gp(n,ntemp) >= ign_max_temp) then
      !
      ! Set the energy source term to zero if the current temperature
      ! is greater than or equal to the ignition max temperature. Also,
      ! mark that the ignition max temperature has been reached to
      ! begin the countdown of the ignition duration.
      !
      ignition_reached_max_temp = true
      !
      energy_source = zero
      !
    else
      !
      ! Compute the value of the energy source term, with the method used
      ! depending on the value of the input ignition_source_radius.
      !
      if (ign_src_radius < zero) then
        !
        ! If ignition_source_radius is less than zero, use a Gaussian function
        ! based on radial distance from the input ignition location as the
        ! magnitude of the energy source term. The input value of
        ! ignition_source_radius should be approximately the radial distance
        ! from the ignition location at which the Gaussian energy source term
        ! reaches machine precision.
        !
        if (r <= abs(ign_src_radius)) then
          !
          r = r/abs(ign_src_radius)
          !
          spark = cv(n,nec) * gp(n,ncp) * ign_src_mag
          !
          energy_source = min(spark,ign_spark_cap)
          !
          energy_source = energy_source*exp(-36.0_wp*r*r)
          !
        end if
        !
        ! If the ignition_source_radius is greater than zero, every point that
        ! is within the ignition_source_radius from the input ignition location
        ! will get the full value of the energy source term based on the input
        ! value of ignition_source_magnitude. The energy source term for all
        ! points outside the ignition_source_radius from the input ignition
        ! location will be an exponential decay based on radius.
        !
      else if (ign_src_radius > zero) then
        !
        spark = cv(n,nec) * gp(n,ncp) * ign_src_mag
        !
        energy_source = min(spark,ign_spark_cap)
        !
        if (r >= ign_src_radius) then
          !
          r = (r-ign_src_radius) / ign_src_radius
          !
          energy_source = energy_source*exp(-36.0_wp*r*r*r)
          !
        end if
        !
      end if
      !
    end if
    !
    ! Add the ignition source term to the RHS of the energy equation
    ! NOTE: We actually need to subtract it from the RHS because the
    !       residual array will be multiplied by -1 at the end of the
    !       current time stage/step before using it to update the
    !       solution variables.
    !
    residual(n,nee) = residual(n,nee) - energy_source*jac(n)/dtsp(n)
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_ignition_source_ncc
!
!###############################################################################
!
subroutine compute_ignition_source_eref(n1,n2,gp,residual)
  !
  !======================================================================
  !
  !   Ignition method based on adding a percentage of the reference
  !   total energy as a source term to the energy equation.
  !
  !======================================================================
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr,xyz
  use flowvar_mod, only : dtsp
  use metrics_mod, only : jac
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: residual(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: r,spark,energy_source
  !
  !.. Local Saved Scalars ..
  real(wp), save :: ign_mag = zero
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_ignition_source_eref"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  !
  !
  ign_mag = ign_mag + ign_src_mag
  !
  spark = min(ign_mag,ign_spark_cap)
  !
  !
  !
  do n = n1,n2
    !
    r = norm2(xyz(1:nr,n)-ign_src_loc(1:nr,1))
    do i = 2,size(ign_src_loc,dim=2)
      r = min( r , norm2(xyz(1:nr,n)-ign_src_loc(1:nr,i)) )
    end do
    !
    energy_source = zero
    !
    if (ign_src_radius < zero) then
      !
      ! If ignition_source_radius is less than zero, use a Gaussian function
      ! based on radial distance from the input ignition location as the
      ! magnitude of the energy source term. The input value of
      ! ignition_source_radius should be approximately the radial distance from
      ! the ignition location at which the Gaussian energy source term reaches
      ! machine precision.
      !
      if (r <= abs(ign_src_radius)) then
        !
        r = r/abs(ign_src_radius)
        !
        energy_source = spark*exp(-36.0_wp*r*r)
        !
        ! If the temperature at this point is greater than the input maximum
        ! ignition temperature, change ignition_prod2reac_L2_ratio_goal to
        ! a negative value to turn off ignition.
        !
        if (gp(n,ntemp) >= ign_max_temp) then
          ignition_reached_max_temp = true
        end if
        !
      end if
      !
      ! If the ignition_source_radius is greater than zero, every point that is
      ! within the ignition_source_radius from the input ignition location will
      ! get the full value of the energy source term based on the input value of
      ! ignition_source_magnitude. The energy source term for all points outside
      ! the ignition_source_radius from the input ignition location will be an
      ! exponential decay based on radius.
      !
    else if (ign_src_radius > zero) then
      !
      energy_source = spark
      !
      if (r >= ign_src_radius) then
        !
        r = (r-ign_src_radius)/ign_src_radius
        !
        energy_source = energy_source*exp(-36.0_wp*r*r*r)
        !
      end if
      !
      ! If the temperature at this point is greater than the input maximum
      ! ignition temperature, change ignition_prod2reac_L2_ratio_goal to
      ! a negative value to turn off ignition.
      !
      if (gp(n,ntemp) >= ign_max_temp) then
        ignition_reached_max_temp = true
      end if
      !
    end if
    !
    ! Add the ignition source term to the RHS of the energy equation
    ! NOTE: We actually need to subtract it from the RHS because the
    !       residual array will be multiplied by -1 at the end of the
    !       current time stage/step before using it to update the
    !       solution variables.
    !
    residual(n,nee) = residual(n,nee) - energy_source*jac(n)/dtsp(n)
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_ignition_source_eref
!
!###############################################################################
!
subroutine compute_ignition_source_reactmb(n1,n2,gp,react_temp)
  !
  !======================================================================
  !
  !   Ignition method based on modifying the temperature value used
  !   to compute the reaction source terms. This reaction temperature
  !   is saved to the array react_temp.
  !
  !======================================================================
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : NCD_Axial
  use geomvar_mod, only : nr,xyz
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: react_temp(:)
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: r
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_ignition_source_reactmb"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (ign_src_loc(NCD_Axial,1) < -ten3) then
    !
    ! If the axial location of the source term is VERY negative,
    ! apply the ignition to all points
    !
    do n = n1,n2
      !
      ! Set the reaction temperatrure to the maximum between the current
      ! temperature at this point and the reaction temperature specified
      ! by the input ignition_source_magnitude.
      !
      react_temp(n) = max(gp(n,ntemp),ign_src_mag)
      !
      ! If the current temperature is greater than or equal to the ignition max
      ! temperature, mark that the ignition max temperature has been reached to
      ! begin the countdown of the ignition duration.
      !
      if (gp(n,ntemp) >= ign_max_temp) then
        ignition_reached_max_temp = true
      end if
      !
    end do
    !
  else
    !
    ! If the axial location is not very negative, only apply the ignition
    ! to points within the ignition source radius
    !
    do n = n1,n2
      !
      r = norm2(xyz(1:nr,n)-ign_src_loc(1:nr,1))
      do i = 2,size(ign_src_loc,dim=2)
        r = min( r , norm2(xyz(1:nr,n)-ign_src_loc(1:nr,i)) )
      end do
      !
      if (gp(n,ntemp) >= ign_max_temp) then
        !
        ! Set the reaction temperature to the current temperature if the current
        ! temperature is greater than or equal to the ignition max temperature.
        ! Also, mark that the ignition max temperature has been reached to
        ! begin the countdown of the ignition duration.
        !
        ignition_reached_max_temp = true
        !
        react_temp(n) = gp(n,ntemp)
        !
      else
        !
        if (r <= ign_src_radius) then
          !
          ! If this point is within ignition_source_radius from the input
          ! ignition location, set the reaction temperatrure to the maximum
          ! between the current temperature at this point and the reaction
          ! temperature specified by the input ignition_source_magnitude.
          !
          react_temp(n) = max(gp(n,ntemp),ign_src_mag)
          !
        else
          !
          react_temp(n) = gp(n,ntemp)
          !
        end if
        !
      end if
      !
    end do
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_ignition_source_reactmb
!
!###############################################################################
!
subroutine compute_ignition_pressure_scaling(n1,n2,gp,sc_pscale)
  !
  !======================================================================
  !
  !   Ignition method based on scaling the species concentrations used
  !   to compute the reaction source terms to achieve a specified reaction
  !   pressure. The argument sc_pscale will returned the scaling value
  !   needed for multiplying the species concentrations to implicitly
  !   set the reaction pressure.
  !
  !======================================================================
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : NCD_Axial
  use geomvar_mod, only : nr,xyz
  use ovar,        only : pref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: sc_pscale(:)
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: r,min_val_pscale,max_val_pscale
  real(wp) :: min_lim_pscale,max_lim_pscale
  !
  logical(lk), save :: this_is_first_output = true
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_ignition_pressure_scaling"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Set limits on the scaling of the species concentrations
  !
  min_lim_pscale = one
  max_lim_pscale = ign_src_mag/pref_nd
  !
  if (ign_src_loc(NCD_Axial,1) < -ten3) then
    !
    ! If the axial location of the source term is VERY negative,
    ! apply the ignition to all points
    !
    do n = n1,n2
      !
      if (gp(n,ntemp) < ign_max_temp) then
        !
        ! If the current temperature is less than the ignition max temperature,
        ! compute the species concentration scaling needed to achieve the
        ! specified reaction pressure.
        !
        if (gp(n,npres) > zero) then
          sc_pscale(n) = max( min_lim_pscale, &
                              min(max_lim_pscale,ign_src_mag/gp(n,npres)) )
        end if
        !
      else
        !
        ! If the current temperature is greater than or equal to the ignition
        ! max temperature, mark that the ignition max temperature has been
        ! reached to begin the countdown of the ignition duration.
        !
        ignition_reached_max_temp = true
        !
        sc_pscale(n) = one
        !
      end if
      !
    end do
    !
  else
    !
    ! If the axial location is not very negative, only apply the ignition
    ! to points within the ignition source radius
    !
    do n = n1,n2
      !
      r = norm2(xyz(1:nr,n)-ign_src_loc(1:nr,1))
      do i = 2,size(ign_src_loc,dim=2)
        r = min( r , norm2(xyz(1:nr,n)-ign_src_loc(1:nr,i)) )
      end do
      !
      if (r <= ign_src_radius) then
        if (gp(n,ntemp) < ign_max_temp) then
          if (gp(n,npres) > zero) then
            sc_pscale(n) = max( min_lim_pscale, &
                                min(max_lim_pscale,ign_src_mag/gp(n,npres)) )
          end if
        end if
      end if
      !
      ! If the current temperature is greater than or equal to the ignition max
      ! temperature, mark that the ignition max temperature has been reached to
      ! begin the countdown of the ignition duration.
      !
      if (gp(n,ntemp) >= ign_max_temp) then
        ignition_reached_max_temp = true
      end if
      !
    end do
    !
  end if
  !
  if (this_is_first_output) then
    !
    this_is_first_output = fals
    !
    min_val_pscale = huge(zero)
    max_val_pscale = -huge(zero)
    !
    do n = n1,n2
      min_val_pscale = min( min_val_pscale , sc_pscale(n) )
      max_val_pscale = max( max_val_pscale , sc_pscale(n) )
    end do
    !
    call mpi_allreduce(MPI_IN_PLACE,min_val_pscale,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,max_val_pscale,1_int_mpi,mpi_flttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
    !
    if (mypnum == glb_root) then
      write (iout,1) "min_lim_pscale",min_lim_pscale, &
                     "max_lim_pscale",max_lim_pscale, &
                     "min_val_pscale",min_val_pscale, &
                     "max_val_pscale",max_val_pscale
      flush (iout)
    end if
    flush (iout)
    !
  end if
  1 format (//,4(2x,a," = ",f16.8,/),/)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_ignition_pressure_scaling
!
!###############################################################################
!
end module ignition_mod
