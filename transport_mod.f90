module transport_mod
  !
  use module_kind_types
  use eqn_idx, only : nec,ncs,nq
  use eqn_idx, only : ntemp,ncp,nvisc,ncond,ng
  !
  implicit none
  !
  private
  !
  public :: init_transport_mod
  public :: update_transport_properties
  public :: compute_mole_fractions
  public :: compute_viscosity_cpg
  !
  integer, protected, save :: nmvc = 4
  integer, protected, save :: nmcc = 4
  !
  real(wp), allocatable, save :: visc_coef(:,:)
  real(wp), allocatable, save :: cond_coef(:,:)
  !
  real(wp), save :: visc_temp_min = 10.0_wp
  real(wp), save :: cond_temp_min = 10.0_wp
  !
  integer, parameter :: CHEMKIN_Curve_Fit = 1
  integer, parameter :: Sutherlands_Law = 2
  !
  integer, parameter :: Constant_Prandtl_Model = 0
  integer, parameter :: NCC_Method = 1
  integer, parameter :: Wassiljewas_Formula = 2
  !
  !
  !
  integer, save :: cond_mix_model  = NCC_Method
  integer, save :: Transport_Model = CHEMKIN_Curve_Fit
  !
  !
  ! #########################################
  ! #####   MODULE GENERIC INTERFACES   #####
  ! #########################################
  !
  interface update_transport_properties
    module procedure update_transport_properties_sp, &
                     update_transport_properties_fp, &
                     update_transport_properties_1p
  end interface update_transport_properties
  !
  interface compute_viscosity_cpg
    module procedure compute_viscosity_cpg_sp, &
                     compute_viscosity_cpg_1p
  end interface compute_viscosity_cpg
  !
contains
!
!###############################################################################
!
subroutine init_transport_mod(cv,gp)
  !
  !.. Use Statements ..
  use ovar, only : itestcase
  use ovar, only : governing_equations
  use ovar, only : suth_Tref,suth_Sref
  use ovar, only : suth_muref,suth_c1
  use ovar, only : muref,muref_nd,reyref
  use ovar, only : rhoref,aref,machref,tref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: cv(:)
  real(wp), intent(inout) :: gp(:)
  !
  !.. Local Scalars ..
  integer  :: ierr
  real(wp) :: muref_scale,tref_scale,reyn
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_transport_mod"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (governing_equations == NavierStokes_Eqns) then
    !
    if (ncs > 0) then
      !
      if (abs(itestcase) == CH4_3D_Diffusion_Flame .or. &
          abs(itestcase) == CH4_Flame_Speed_Test) then
        cond_mix_model = NCC_Method
        Transport_Model = CHEMKIN_Curve_Fit
      else if (itestcase == Reacting_Backstep) then
        cond_mix_model = NCC_Method
        Transport_Model = CHEMKIN_Curve_Fit
      else if (itestcase == H2_2D_Diffusion_Flame) then
        cond_mix_model = Wassiljewas_Formula
        Transport_Model = Sutherlands_Law
      else if (itestcase == Compression_Corner) then
        cond_mix_model = Wassiljewas_Formula
        Transport_Model = Sutherlands_Law
      else if (itestcase == N2_O2_Diffusion) then
        cond_mix_model = Wassiljewas_Formula
        Transport_Model = Sutherlands_Law
      end if
      !
      if (Transport_Model == CHEMKIN_Curve_Fit) then
        nmvc = 4
        nmcc = 4
      else if (Transport_Model == Sutherlands_Law) then
        nmvc = 2
        nmcc = 2
      end if
      !
      if (.not. allocated(visc_coef)) then
        allocate ( visc_coef(1:nmvc,1:ncs) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"visc_coef",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      if (.not. allocated(cond_coef)) then
        allocate ( cond_coef(1:nmcc,1:ncs) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"cond_coef",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      if (Transport_Model == CHEMKIN_Curve_Fit) then
        !
        visc_temp_min = 10.0_wp
        cond_temp_min = 10.0_wp
        !
        call hardwired_CHEMKIN_coefficients
        !
      else if (Transport_Model == Sutherlands_Law) then
        !
        call hardwired_Sutherlands_coefficients
        !
      end if
      !
    else
      !
      muref_scale = (one + suth_Sref/suth_Tref) / sqrt(suth_Tref)
      !
      if (reyref <= zero) then
        !
        suth_c1 = suth_muref * muref_scale
        !
        call compute_viscosity_cpg(gp(ntemp),muref_nd)
        !
        muref = muref_nd * rhoref * aref
        !
        ! NOTE: reciprocal returns 0 if muref is 0, and 1/muref otherwise.
        reyref = rhoref * machref * aref * reciprocal( muref )
        !
      else
        !
        muref = rhoref * machref * aref / reyref
        !
        muref_nd = muref / (rhoref * aref)
        !
        tref_scale = tref * sqrt(tref) / (tref + suth_Sref)
        !
        suth_muref = muref / muref_scale / tref_scale
        !
        suth_c1 = suth_muref * muref_scale
        !
      end if
      !
    end if
    !
    call update_transport_properties(cv,gp)
    !
    if (ncs > 0) then
      !
      muref_nd = gp(nvisc)
      muref = muref_nd * rhoref * aref
      !
      reyn = rhoref * machref * aref / muref
      !
      if (reyref <= zero) then
        reyref = reyn
      else
        if (abs(reyref-reyn) > eps6) then
          write (error_message,100) reyref,reyn
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
        end if
      end if
      !
    end if
    !
  else
    !
    muref = zero
    muref_nd = zero
    suth_c1 = zero
    reyref = zero
    !
  end if
  !
  !$ACC UPDATE DEVICE(suth_Tref,suth_Sref,suth_muref,suth_c1)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  100 format ( &
     " Free-stream conditions in the input file are not consistent!",'/&', &
     " Reference Reynolds number given    : ",es14.6,'/&', &
     " Reference Reynolds number computed : ",es14.6,'/&', &
     " Please correct the input file or set the reference Reynolds number", &
     " in the input file to a negative value so it is explicitly computed.")
  !
end subroutine init_transport_mod
!
!###############################################################################
!
!pure &
subroutine update_transport_properties_1p(cv_1p,gp_1p)
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: cv_1p(:)
  real(wp), intent(inout) :: gp_1p(:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: cv(:,:)
  real(wp), allocatable :: gp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_transport_properties_1p"
  !
continue
  !
  allocate ( cv(1:1,1:nq) , source=zero )
  allocate ( gp(1:1,1:ng) , source=zero )
  !
  cv(1,1:nq) = cv_1p(1:nq)
  gp(1,1:ng) = gp_1p(1:ng)
  !$ACC DATA CREATE(cv,gp)
  !
  call update_transport_properties_sp(1,1,cv,gp)
  !
  gp_1p(1:ng) = gp(1,1:ng)
  !
  !$ACC END DATA
  deallocate ( gp )
  deallocate ( cv )
  !
end subroutine update_transport_properties_1p
!
!###############################################################################
!
!pure &
subroutine update_transport_properties_fp(n1,n2,cv,gp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: cv(:,:,:)
  real(wp), intent(inout) :: gp(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_transport_properties_fp"
  !
continue
  !
  call update_transport_properties_sp(n1,n2,cv(:,:,1),gp(:,:,1))
  if (all([size(cv,dim=3),size(gp,dim=3)] > 1)) then
    call update_transport_properties_sp(n1,n2,cv(:,:,2),gp(:,:,2))
  end if
  !
end subroutine update_transport_properties_fp
!
!###############################################################################
!
!pure &
subroutine update_transport_properties_sp(n1,n2,cv,gp)
  !
  !.. Use Statements ..
  use ovar,    only : Pr,cpref_nd
  use ovar,    only : governing_equations
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp), intent(inout) :: gp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: Pr_inv,cp_divided_by_Pr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xi(:,:)
  real(wp), allocatable :: rmui(:,:)
  real(wp), allocatable :: rmci(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_transport_properties_sp"
  !
continue
  !
  if (governing_equations == Euler_Eqns) return
  !
  if (ncs > 0) then
    !
    allocate ( xi(n1:n2,1:ncs) , source=zero )
    allocate ( rmui(n1:n2,1:ncs) , source=zero )
    !$ACC DATA CREATE(xi,rmui)
    !
    ! 1. Compute the mole fraction of each species
    !
    call compute_mole_fractions(n1,n2,cv,xi)
    !
    ! 2. Compute the species molecular viscosities
    !
    call compute_species_viscosities(n1,n2,gp(:,ntemp),rmui, &
                                     imodel=Transport_Model)
    !
    ! 3. Compute the mixture molecular viscosity
    !
    call compute_mixture_viscosity(n1,n2,rmui,xi,gp(:,nvisc))
    !
    !$ACC EXIT DATA DELETE(rmui)
    deallocate ( rmui )
    !
    ! 4. Compute the species molecular thermal conductivities if the
    !    mixture molecular thermal conductivity model requires them
    !
    if (cond_mix_model /= Constant_Prandtl_Model) then
      !
      ! Allocate memory for the species molecular conductivity
      !
      allocate ( rmci(n1:n2,1:ncs) , source=zero )
      !$ACC ENTER DATA CREATE(rmci)
      !
      call compute_species_conductivities(n1,n2,gp(:,ntemp),rmci, &
                                          imodel=Transport_Model)
      !
      call compute_mixture_conductivity(n1,n2,rmci,xi,gp(:,ncond), &
                                        imodel=cond_mix_model)
      !
      !$ACC EXIT DATA DELETE(rmci)
      deallocate ( rmci )
      !
    else
      !
      Pr_inv = one / Pr
      !
      !$ACC PARALLEL LOOP PRESENT(gp)
      do n = n1,n2
        gp(n,ncond) = gp(n,nvisc) * gp(n,ncp) * Pr_inv
      end do
      !
    end if
    !
    !$ACC EXIT DATA DELETE(xi)
    deallocate ( xi )
    !
  else
    !
    ! Compute the viscosity for the calorically perfect gas
    !
    call compute_viscosity_cpg(n1,n2,gp(:,ntemp),gp(:,nvisc))
    !
    cp_divided_by_Pr = cpref_nd / Pr
    !
    !$ACC PARALLEL LOOP PRESENT(gp)
    do n = n1,n2
      gp(n,ncond) = gp(n,nvisc) * cp_divided_by_Pr
    end do
    !
  end if
  !
end subroutine update_transport_properties_sp
!
!###############################################################################
!
!pure &
subroutine compute_viscosity_cpg_sp(n1,n2,t,visc)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref
  use ovar, only : suth_c1,suth_sref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: visc(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: tdim,mudim,rhoaref_inv
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_viscosity_cpg_sp"
  !
continue
  !
  rhoaref_inv = one / rhoref / aref
  !
  !$ACC PARALLEL LOOP PRESENT(t,visc)
  do n = n1,n2
    tdim = t(n)*tref
    mudim = suth_c1 * tdim * sqrt(tdim) / (tdim + suth_sref)
    visc(n) = mudim * rhoaref_inv
  end do
  !
end subroutine compute_viscosity_cpg_sp
!
!###############################################################################
!
!pure &
subroutine compute_viscosity_cpg_1p(t,visc)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref
  use ovar, only : suth_c1,suth_sref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: t
  real(wp), intent(inout) :: visc
  !
  !.. Local Scalars ..
  real(wp) :: tdim,mudim
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_viscosity_cpg_1p"
  !
continue
  !
  tdim = t*tref
  mudim = suth_c1 * tdim * sqrt(tdim) / (tdim + suth_sref)
  visc = mudim / rhoref / aref
  !
end subroutine compute_viscosity_cpg_1p
!
!###############################################################################
!
pure subroutine compute_mixture_viscosity(n1,n2,rmui,mole_frac,visc)
  !
  !.. Use Statements ..
  use thermo_mod, only : wtmol
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: rmui(n1:n2,1:ncs)
  real(wp),    intent(in) :: mole_frac(n1:n2,1:ncs)
  real(wp), intent(inout) :: visc(:)
  !
  !.. Local Scalars ..
  integer  :: i,j,n
  real(wp) :: val,rmw_ji,rrmw_ji,rmu_ij,phi_ij,denomi
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: phi_sum(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_viscosity"
  !
continue
  !
  allocate ( phi_sum(n1:n2,1:ncs) , source=zero )
  !$ACC DATA CREATE(phi_sum)
  !
  ! Compute the mixture molecular viscosity
  !
  !$ACC PARALLEL LOOP PRESENT(visc)
  do n = n1,n2
    visc(n) = zero
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(phi_sum,rmui,mole_frac,wtmol)
  do i = 1,ncs
    do n = n1,n2
      !
      ! Compute SUM( mole_frac(j)*Phi(ij) , j=1,ncs , j/=i )
      !
      val = zero
      !$ACC LOOP SEQ REDUCTION(+:val)
      do j = 1,ncs
        !
        ! Exclude species diffusing into itself from loop
        ! since this is already built into the mixture formula
        !
        if (j /= i) then
          !
          rmw_ji = wtmol(j) / wtmol(i)
          rrmw_ji = sqrt(rmw_ji)
          !
          denomi = one / sqrt(eight+eight/rmw_ji)
          !
          rmu_ij = rmui(n,i) / rmui(n,j)
          phi_ij = ((one+sqrt(rmu_ij*rrmw_ji))**2) * denomi
          val = val + mole_frac(n,j)*phi_ij
          !
        end if
        !
      end do
      !
      phi_sum(n,i) = val
      !
    end do
  end do
  !
  !$ACC PARALLEL LOOP PRESENT(phi_sum,rmui,mole_frac)
  do n = n1,n2
    val = zero
    !$ACC LOOP SEQ REDUCTION(+:val)
    do i = 1,ncs
      val = val + rmui(n,i)*mole_frac(n,i)/(mole_frac(n,i)+phi_sum(n,i))
    end do
    visc(n) = val
  end do
  !
  !$ACC END DATA
  deallocate ( phi_sum )
  !
end subroutine compute_mixture_viscosity
!
!###############################################################################
!
pure subroutine compute_mixture_conductivity(n1,n2,rmci,mole_frac,cond,imodel)
  !
  !.. Use Statements ..
  use thermo_mod, only : wtmol
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: rmci(n1:n2,1:ncs)
  real(wp),    intent(in) :: mole_frac(n1:n2,1:ncs)
  real(wp), intent(inout) :: cond(:)
  integer,     intent(in) :: imodel
  !
  !.. Local Scalars ..
  integer  :: i,j,n
  real(wp) :: val,sum_1,sum_2,denomi
  real(wp) :: rmw_ji,rrmw_ji,rmc_ij,phi_ij
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: phi_sum(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_conductivity"
  !
continue
  !
  ! Compute the mixture thermal conductivity
  !
  !$ACC PARALLEL LOOP PRESENT(cond)
  do n = n1,n2
    cond(n) = zero
  end do
  !
  if (imodel == Wassiljewas_Formula) then
    !
    allocate ( phi_sum(n1:n2,1:ncs) , source=zero )
    !$ACC DATA CREATE(phi_sum)
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(phi_sum,rmci,mole_frac,wtmol)
    do i = 1,ncs
      do n = n1,n2
        !
        ! Compute SUM( mole_frac(j)*Phi(ij) , j=1,ncs , j/=i )
        !
        val = zero
        !$ACC LOOP SEQ REDUCTION(+:val)
        do j = 1,ncs
          !
          ! Exclude species diffusing into itself from loop
          ! since this is already built into the mixture formula
          !
          if (j /= i) then
            !
            rmw_ji = wtmol(j) / wtmol(i)
            rrmw_ji = sqrt(rmw_ji)
            !
            denomi = one / sqrt(eight+eight/rmw_ji)
            !
            rmc_ij = rmci(n,i) / rmci(n,j)
            phi_ij = ((one+sqrt(rmc_ij*rrmw_ji))**2) * denomi
            val = val + mole_frac(n,j)*phi_ij
            !
          end if
          !
        end do
        !
        phi_sum(n,i) = 1.065_wp*val
        !
      end do
    end do
    !
    !$ACC PARALLEL LOOP PRESENT(phi_sum,rmci,mole_frac,cond)
    do n = n1,n2
      val = zero
      !$ACC LOOP SEQ REDUCTION(+:val)
      do i = 1,ncs
        val = val + rmci(n,i)*mole_frac(n,i)/(mole_frac(n,i)+phi_sum(n,i))
      end do
      cond(n) = val
    end do
    !
    !$ACC END DATA
    deallocate ( phi_sum )
    !
  else
    !
    !$ACC PARALLEL LOOP PRESENT(mole_frac,rmci,cond)
    do n = n1,n2
      sum_1 = zero
      sum_2 = zero
      !$ACC LOOP SEQ REDUCTION(+:sum_1,sum_2)
      do i = 1,ncs
        sum_1 = sum_1 + mole_frac(n,i)*rmci(n,i)
        sum_2 = sum_2 + mole_frac(n,i)/rmci(n,i)
      end do
      cond(n) = half*(sum_1 + one/sum_2)
    end do
    !
  end if
  !
end subroutine compute_mixture_conductivity
!
!###############################################################################
!
pure subroutine compute_mole_fractions(n1,n2,cv,mole_frac)
  !
  !.. Use Statements ..
  use thermo_mod, only : wtmol
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp), intent(inout) :: mole_frac(n1:n2,1:ncs)
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: val
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: inv_mole_sum(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mole_fractions"
  !
continue
  !
  allocate ( inv_mole_sum(n1:n2) , source=zero )
  !$ACC DATA CREATE(inv_mole_sum)
  !
  ! NOTE: IT DOES NOT MATTER WHETHER THE SPECIES VALUES CONTAINED
  !       IN cv(:,1:ncs) ARE THE SPECIES DENSITIES OR THE SPECIES
  !       MASS FRACTIONS. EITHER WAY, THIS ALGORITHM WILL CORRECTLY
  !       COMPUTE THE MOLE FRACTIONS FOR THE GAS MIXTURE.
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cv,wtmol,mole_frac)
  do i = 1,ncs
    do n = n1,n2
      mole_frac(n,i) = (cv(n,i)/cv(n,nec)) / wtmol(i)
    end do
  end do
  !
  !$ACC PARALLEL LOOP PRESENT(mole_frac,inv_mole_sum)
  do n = n1,n2
    val = zero
    !$ACC LOOP SEQ REDUCTION(+:val)
    do i = 1,ncs
      val = val + mole_frac(n,i)
    end do
    inv_mole_sum(n) = one/val
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(mole_frac,inv_mole_sum)
  do i = 1,ncs
    do n = n1,n2
      mole_frac(n,i) = mole_frac(n,i) * inv_mole_sum(n)
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( inv_mole_sum )
  !
end subroutine compute_mole_fractions
!
!###############################################################################
!
pure subroutine compute_species_viscosities(n1,n2,t,rmui,imodel)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: rmui(n1:n2,1:ncs)
  integer,     intent(in) :: imodel
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: logt,a1,a2,a3,a4
  real(wp) :: tnd,c1_suth,s_suth
  real(wp) :: tdim,viscref_inv
  !
  !.. Local Parameters ..
  real(wp), parameter :: chemkin_cgs_to_mks = 1.0e-1_wp
  !
continue
  !
  if (imodel == Sutherlands_Law) then
    !
    ! Sutherlands law fit
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(visc_coef,t,rmui)
    do i = 1,ncs
      do n = n1,n2
        !
        tnd = t(n) ! non-dimensional temperature
        !
        c1_suth = visc_coef(1,i) ! non-dimensional Sutherland constants
        s_suth  = visc_coef(2,i) ! non-dimensional Sutherland constants
        !
        rmui(n,i) = c1_suth * tnd * sqrt(tnd) / (tnd + s_suth)
        !
      end do
    end do
    !
  else if (imodel == CHEMKIN_Curve_Fit) then
    !
    ! ChemKin fit
    !
    ! Constant to convert the curve fit results from
    ! CGS (centimeters-grams-seconds) units to
    ! SI/MKS (meters-kilograms-seconds) units,
    ! and non-dimensionalize the MKS results.
    !
    viscref_inv = chemkin_cgs_to_mks / rhoref / aref
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(visc_coef,t,rmui)
    do i = 1,ncs
      do n = n1,n2
        !
        tdim = max( t(n)*tref , visc_temp_min ) ! dimensional temperature
        !
        logt = log(tdim)
        !
        a1 = visc_coef(1,i) ! ChemKin coefficients.
        a2 = visc_coef(2,i) ! I have no idea what the dimensions of
        a3 = visc_coef(3,i) ! these are, but the curve fit requires
        a4 = visc_coef(4,i) ! dimensional temperature as input.
        !
       !rmui(n,i) = chemkin_cgs_to_mks * &
       !            exp(a1 + logt*(a2 + logt*(a3 + logt*a4)))
        rmui(n,i) = viscref_inv * exp(a1 + logt*(a2 + logt*(a3 + logt*a4)))
        !
      end do
    end do
    !
  end if
  !
end subroutine compute_species_viscosities
!
!###############################################################################
!
pure subroutine compute_species_conductivities(n1,n2,t,rmci,imodel)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: rmci(n1:n2,1:ncs)
  integer,     intent(in) :: imodel
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: logt,a1,a2,a3,a4
  real(wp) :: tnd,c1_suth,s_suth
  real(wp) :: tdim,condref_inv
  !
  !.. Local Parameters ..
  real(wp), parameter :: chemkin_cgs_to_mks = 1.0e-5_wp
  !
continue
  !
  if (imodel == Sutherlands_Law) then
    !
    ! Sutherlands law fit
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cond_coef,t,rmci)
    do i = 1,ncs
      do n = n1,n2
        !
        tnd = t(n) ! non-dimensional temperature
        !
        c1_suth = cond_coef(1,i) ! non-dimensional Sutherland constants
        s_suth  = cond_coef(2,i) ! non-dimensional Sutherland constants
        !
        rmci(n,i) = c1_suth * tnd * sqrt(tnd) / (tnd + s_suth)
        !
      end do
    end do
    !
  else if (imodel == CHEMKIN_Curve_Fit) then
    !
    ! ChemKin fit
    !
    ! Constant to convert the curve fit results from
    ! CGS (centimeters-grams-seconds) units to
    ! SI/MKS (meters-kilograms-seconds) units,
    ! and non-dimensionalize the MKS results.
    !
    condref_inv = chemkin_cgs_to_mks * tref/(rhoref*aref*aref*aref)
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cond_coef,t,rmci)
    do i = 1,ncs
      do n = n1,n2
        !
        tdim = max( t(n)*tref , cond_temp_min ) ! dimensional temperature
        !
        logt = log(tdim)
        !
        a1 = cond_coef(1,i) ! ChemKin coefficients.
        a2 = cond_coef(2,i) ! I have no idea what the dimensions of
        a3 = cond_coef(3,i) ! these are, but the curve fit requires
        a4 = cond_coef(4,i) ! dimensional temperature as input.
        !
       !rmci(n,i) = chemkin_cgs_to_mks * &
       !            exp(a1 + logt*(a2 + logt*(a3 + logt*a4)))
        rmci(n,i) = condref_inv * exp(a1 + logt*(a2 + logt*(a3 + logt*a4)))
        !
      end do
    end do
    !
  end if
  !
end subroutine compute_species_conductivities
!
!###############################################################################
!
pure subroutine mixture_transport_properties(imodel,t,rhoi,rmui,rmci,rmu,rmc, &
                                             mixture_method)
  !
  !.. Use Statements ..
  use thermo_mod, only : wtmol
  !
  ! Routine to calculate the molecular viscosity and thermal conductivity
  ! for a multi-component gas
  !
  ! T     : static temperature
  ! RHOI  : species mass fractions
  ! RMU   : molecular viscosity
  ! RMC   : molecular thermal conductivity
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: imodel
  real(wp),    intent(in) :: t
  real(wp),    intent(in) :: rhoi(:)
  real(wp), intent(inout) :: rmui(:)
  real(wp), intent(inout) :: rmci(:)
  real(wp), intent(inout) :: rmc
  real(wp), intent(inout) :: rmu
  !
  integer, optional, intent(in) :: mixture_method
  !
  !.. Local Scalars ..
  integer :: i,ics,jcs,mix_method
  !
  real(wp) :: logt
  real(wp) :: phi_ij,rmw_ji,rrmw_ji
  real(wp) :: rmc0,rmu0,s0,t0,tot0
  real(wp) :: denomi
  real(wp) :: rmu_ij
  real(wp) :: work,sum_1,sum_2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:4) :: a
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xi(:)
  !
  !.. Local Parameters ..
  integer, parameter :: NCC_Method = 1
  integer, parameter :: VULCAN_Method = 2
  !
continue
  !
  mix_method = NCC_Method
  if (present(mixture_method)) then
    mix_method = min( max(mixture_method,NCC_Method) , VULCAN_Method )
  end if
  !
  ! Allocate memory for the species mole fractions
  !
  allocate ( xi(1:ncs) )
  !
  ! Compute the mole fraction of each species
  !
  work = zero
  !
  do ics = 1,ncs
    xi(ics) = rhoi(ics) / wtmol(ics)
    work = work + xi(ics)
  end do
  !
  do ics = 1,ncs
    xi(ics) = xi(ics) / work
  end do
  !
  ! Compute the species molecular viscosities
  !
  do ics = 1,ncs
    !
    if (imodel == Sutherlands_Law) then
      !
      ! Sutherlands law fit, 1 temperature interval
      !
      rmu0 = visc_coef(1,ics)
      t0   = visc_coef(2,ics)
      s0   = visc_coef(3,ics)
      !
      rmu0 = rmu0 * (t0+s0)
      !
      tot0 = t / t0
      rmui(ics) = rmu0 * tot0 * sqrt(tot0) / (t+s0)
      !
    else if (imodel == CHEMKIN_Curve_Fit) then
      !
      ! ChemKin fit, 1 temperature interval
      !
      logt = log(t)
      !
      a(1:4) = visc_coef(1:4,ics)
      !
      rmui(ics) = exp(a(1) + logt*(a(2) + logt*(a(3) + logt*a(4))))
      rmui(ics) = 1.0e-01_wp * rmui(ics)
      !
    end if
    !
  end do
  !
  ! Compute the species molecular thermal conductivities
  ! if the mixture molecular thermal conductivity model requires them
  !
  do ics = 1,ncs
    !
    if (imodel == Sutherlands_Law) then
      !
      ! Sutherlands law fit, 1 temperature interval
      !
      rmc0 = cond_coef(1,ics)
      t0   = cond_coef(2,ics)
      s0   = cond_coef(3,ics)
      !
      rmc0 = rmc0 * (t0+s0)
      !
      tot0 = t / t0
      rmci(ics) = rmc0 * tot0 * sqrt(tot0) / (t+s0)
      !
    else if (imodel == CHEMKIN_Curve_Fit) then
      !
      ! ChemKin fit, 1 temperature interval
      !
      logt = log(t)
      !
      a(1:4) = cond_coef(1:4,ics)
      !
      rmci(ics) = exp(a(1) + logt*(a(2) + logt*(a(3) + logt*a(4))))
      rmci(ics) = 1.0e-05_wp * rmci(ics)
      !
    end if
    !
  end do
  !
  ! Compute the mixture molecular viscosity and thermal conductivity
  !
  rmu = zero
  rmc = zero
  !
  do ics = 1,ncs
    !
    ! Compute Phi(ij)
    !
    work = zero
    !
    do jcs = 1,ncs
      !
      ! Exclude species diffusing into itself from loop
      ! since this is already built into the mixture formula
      !
      if (jcs /= ics) then
        !
        rmw_ji = wtmol(jcs) / wtmol(ics)
        rrmw_ji = sqrt(rmw_ji)
        !
        denomi = one / sqrt(eight+eight/rmw_ji)
        !
        rmu_ij = rmui(ics) / rmui(jcs)
        phi_ij = ((one+sqrt(rmu_ij*rrmw_ji))**2) * denomi
        work = work + xi(jcs)*phi_ij
        !
      end if
      !
    end do
    !
    ! Add each species contribution to the mixture molecular viscosity
    !
    rmu = rmu + rmui(ics)*xi(ics)/(xi(ics)+work)
    !
    ! Add each species contribution to the mixture molecular conductivity
    !
    if (mix_method == VULCAN_Method) then
      rmc = rmc + rmci(ics)*xi(ics)/(xi(ics)+1.065_wp*work)
    end if
    !
  end do
  !
  if (mix_method == NCC_Method) then
    !
    sum_1 = zero
    sum_2 = zero
    do ics = 1,ncs
      sum_1 = sum_1 + xi(ics)*rmci(ics)
      sum_2 = sum_2 + xi(ics)/rmci(ics)
    end do
    !
    rmc = half*(sum_1 + one/sum_2)
    !
  end if
  !
  deallocate ( xi )
  !
end subroutine mixture_transport_properties
!
!###############################################################################
!
subroutine hardwired_CHEMKIN_coefficients
  !
  !.. Use Statements ..
  use thermo_mod, only : species_name
  !
  !======================================================================
  !
  ! COEFFICIENTS FOR SPECIES VISCOSITIES
  !
  !  MAXIMUM FITTING ERROR =    3.645e-03
  !
  !  CO2         -22.7   4.580  -0.492   0.02050
  !  CO          -15.6   1.950  -0.169   0.00730
  !  H2O         -14.9   0.561   0.137  -0.01040
  !  CH4         -18.6   2.940  -0.295   0.01260
  !  N2          -15.5   1.930  -0.167   0.00723
  !  O2          -16.0   2.170  -0.198   0.00854
  !
  ! CH4         -0.186e+02   0.294e+01  -0.295e+00   0.126e-01
  ! C11H21      -0.232e+02   0.396e+01  -0.352e+00   0.123e-01
  ! CH          -0.141e+02   0.144e+01  -0.103e+00   0.440e-02
  ! O2          -0.160e+02   0.217e+01  -0.198e+00   0.854e-02
  ! CO          -0.156e+02   0.195e+01  -0.169e+00   0.730e-02
  ! H2O         -0.149e+02   0.561e+00   0.137e+00  -0.104e-01
  ! N2          -0.155e+02   0.193e+01  -0.167e+00   0.723e-02
  ! O           -0.140e+02   0.144e+01  -0.103e+00   0.440e-02
  ! N           -0.140e+02   0.128e+01  -0.820e-01   0.350e-02
  ! NO          -0.155e+02   0.193e+01  -0.167e+00   0.723e-02
  ! OH          -0.140e+02   0.144e+01  -0.103e+00   0.440e-02
  ! H           -0.189e+02   0.300e+01  -0.303e+00   0.129e-01
  ! H2          -0.136e+02   0.884e+00  -0.325e-01   0.147e-02
  ! CO2         -0.227e+02   0.458e+01  -0.492e+00   0.205e-01
  ! N2O         -0.224e+02   0.446e+01  -0.479e+00   0.200e-01
  !
  !======================================================================
  !
  real(wp), parameter :: CO2_visc_temp_min = 10.0_wp
  !
  real(wp), parameter :: CO2_visc(1:4) = [-0.227e+02_wp, &
                                           0.458e+01_wp, &
                                          -0.492e+00_wp, &
                                           0.205e-01_wp]
  !
  real(wp), parameter :: CO_visc_temp_min = 10.0_wp
  !
  real(wp), parameter :: CO_visc(1:4) = [-0.156e+02_wp, &
                                          0.195e+01_wp, &
                                         -0.169e+00_wp, &
                                          0.730e-02_wp]
  !
  real(wp), parameter :: H2O_visc_temp_min = 10.0_wp
  !
  real(wp), parameter :: H2O_visc(1:4) = [-0.149e+02_wp, &
                                           0.561e+00_wp, &
                                           0.137e+00_wp, &
                                          -0.104e-01_wp]
  !
  real(wp), parameter :: CH4_visc_temp_min = 10.0_wp
  !
  real(wp), parameter :: CH4_visc(1:4) = [-0.186e+02_wp, &
                                           0.294e+01_wp, &
                                          -0.295e+00_wp, &
                                           0.126e-01_wp]
  !
  real(wp), parameter :: N2_visc_temp_min = 10.0_wp
  !
  real(wp), parameter :: N2_visc(1:4) = [-0.155e+02_wp, &
                                          0.193e+01_wp, &
                                         -0.167e+00_wp, &
                                          0.723e-02_wp]
  !
  real(wp), parameter :: O2_visc_temp_min = 10.0_wp
  !
  real(wp), parameter :: O2_visc(1:4) = [-0.160e+02_wp, &
                                          0.217e+01_wp, &
                                         -0.198e+00_wp, &
                                          0.854e-02_wp]
  !
  real(wp), parameter :: H_visc_temp_min = 10.0_wp
  !
  real(wp), parameter :: H_visc(1:4) = [-0.189e+02_wp, &
                                         0.300e+01_wp, &
                                        -0.303e+00_wp, &
                                         0.129e-01_wp]
  !
  real(wp), parameter :: H2_visc_temp_min = 10.0_wp
  !
  real(wp), parameter :: H2_visc(1:4) = [-0.136e+02_wp, &
                                          0.884e+00_wp, &
                                         -0.325e-01_wp, &
                                          0.147e-02_wp]
  !
  real(wp), parameter :: O_visc_temp_min = 10.0_wp
  !
  real(wp), parameter :: O_visc(1:4) = [-0.140e+02_wp, &
                                         0.144e+01_wp, &
                                        -0.103e+00_wp, &
                                         0.440e-02_wp]
  !
  real(wp), parameter :: OH_visc_temp_min = 10.0_wp
  !
  real(wp), parameter :: OH_visc(1:4) = [-0.140e+02_wp, &
                                          0.144e+01_wp, &
                                         -0.103e+00_wp, &
                                          0.440e-02_wp]
  !
  !======================================================================
  !
  ! COEFFICIENTS FOR SPECIES CONDUCTIVITIES
  !
  !  MAXIMUM FITTING ERROR =    1.098e-02
  !
  ! CH4         -0.469e+00   0.108e+01   0.159e+00  -0.148e-01
  ! C11H21      -0.244e+02   0.101e+02  -0.994e+00   0.346e-01
  ! CH           0.245e+02  -0.799e+01   0.125e+01  -0.590e-01
  ! O2          -0.213e+01   0.299e+01  -0.287e+00   0.124e-01
  ! CO           0.582e+01  -0.481e+00   0.211e+00  -0.114e-01
  ! H2O          0.163e+02  -0.589e+01   0.109e+01  -0.551e-01
  ! N2           0.760e+01  -0.118e+01   0.303e+00  -0.154e-01
  ! O            0.278e+01   0.144e+01  -0.103e+00   0.440e-02
  ! N            0.292e+01   0.128e+01  -0.820e-01   0.350e-02
  ! NO           0.423e+01   0.217e+00   0.113e+00  -0.688e-02
  ! OH           0.165e+02  -0.432e+01   0.696e+00  -0.314e-01
  ! H            0.604e+00   0.300e+01  -0.303e+00   0.129e-01
  ! H2           0.132e+02  -0.225e+01   0.381e+00  -0.157e-01
  ! CO2         -0.140e+02   0.706e+01  -0.744e+00   0.289e-01
  ! N2O         -0.130e+02   0.667e+01  -0.692e+00   0.266e-01
  !
  !======================================================================
  !
  real(wp), parameter :: CO2_cond_temp_min = 10.0_wp
  !
  real(wp), parameter :: CO2_cond(1:4) = [-0.140e+02_wp, &
                                           0.706e+01_wp, &
                                          -0.744e+00_wp, &
                                           0.289e-01_wp]
  !
  real(wp), parameter :: CO_cond_temp_min = 10.0_wp
  !
  real(wp), parameter :: CO_cond(1:4) = [ 0.582e+01_wp, &
                                         -0.481e+00_wp, &
                                          0.211e+00_wp, &
                                         -0.114e-01_wp]
  !
  real(wp), parameter :: H2O_cond_temp_min = 10.0_wp
  !
  real(wp), parameter :: H2O_cond(1:4) = [ 0.163e+02_wp, &
                                          -0.589e+01_wp, &
                                           0.109e+01_wp, &
                                          -0.551e-01_wp]
  !
  real(wp), parameter :: CH4_cond_temp_min = 10.0_wp
  !
  real(wp), parameter :: CH4_cond(1:4) = [-0.469e+00_wp, &
                                           0.108e+01_wp, &
                                           0.159e+00_wp, &
                                          -0.148e-01_wp]
  !
  real(wp), parameter :: N2_cond_temp_min = 10.0_wp
  !
  real(wp), parameter :: N2_cond(1:4) = [ 0.760e+01_wp, &
                                         -0.118e+01_wp, &
                                          0.303e+00_wp, &
                                         -0.154e-01_wp]
  !
  real(wp), parameter :: O2_cond_temp_min = 10.0_wp
  !
  real(wp), parameter :: O2_cond(1:4) = [-0.213e+01_wp, &
                                          0.299e+01_wp, &
                                         -0.287e+00_wp, &
                                          0.124e-01_wp]
  !
  real(wp), parameter :: H_cond_temp_min = 10.0_wp
  !
  real(wp), parameter :: H_cond(1:4) = [ 0.604e+00_wp, &
                                         0.300e+01_wp, &
                                        -0.303e+00_wp, &
                                         0.129e-01_wp]
  !
  real(wp), parameter :: H2_cond_temp_min = 10.0_wp
  !
  real(wp), parameter :: H2_cond(1:4) = [ 0.132e+02_wp, &
                                         -0.225e+01_wp, &
                                          0.381e+00_wp, &
                                         -0.157e-01_wp]
  !
  real(wp), parameter :: O_cond_temp_min = 10.0_wp
  !
  real(wp), parameter :: O_cond(1:4) = [ 0.278e+01_wp, &
                                         0.144e+01_wp, &
                                        -0.103e+00_wp, &
                                         0.440e-02_wp]
  !
  real(wp), parameter :: OH_cond_temp_min = 10.0_wp
  !
  real(wp), parameter :: OH_cond(1:4) = [ 0.165e+02_wp, &
                                         -0.432e+01_wp, &
                                          0.696e+00_wp, &
                                         -0.314e-01_wp]
  !
  integer :: i
  !
  character(len=len(species_name)) :: species
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "hardwired_CHEMKIN_coefficients"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  do i = 1,ncs
    !
    species = uppercase( trim(adjustl(species_name(i))) )
    !
    select case (trim(adjustl(species)))
      !
      case ("CO2")
        !
        visc_coef(1:4,i) = CO2_visc(1:4)
        cond_coef(1:4,i) = CO2_cond(1:4)
        !
      case ("CO")
        !
        visc_coef(1:4,i) = CO_visc(1:4)
        cond_coef(1:4,i) = CO_cond(1:4)
        !
      case ("H2O")
        !
        visc_coef(1:4,i) = H2O_visc(1:4)
        cond_coef(1:4,i) = H2O_cond(1:4)
        !
      case ("CH4")
        !
        visc_coef(1:4,i) = CH4_visc(1:4)
        cond_coef(1:4,i) = CH4_cond(1:4)
        !
      case ("N2")
        !
        visc_coef(1:4,i) = N2_visc(1:4)
        cond_coef(1:4,i) = N2_cond(1:4)
        !
      case ("O2")
        !
        visc_coef(1:4,i) = O2_visc(1:4)
        cond_coef(1:4,i) = O2_cond(1:4)
        !
      case ("H")
        !
        visc_coef(1:4,i) = H_visc(1:4)
        cond_coef(1:4,i) = H_cond(1:4)
        !
      case ("H2")
        !
        visc_coef(1:4,i) = H2_visc(1:4)
        cond_coef(1:4,i) = H2_cond(1:4)
        !
      case ("O")
        !
        visc_coef(1:4,i) = O_visc(1:4)
        cond_coef(1:4,i) = O_cond(1:4)
        !
      case ("OH")
        !
        visc_coef(1:4,i) = OH_visc(1:4)
        cond_coef(1:4,i) = OH_cond(1:4)
        !
    end select
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine hardwired_CHEMKIN_coefficients
!
!###############################################################################
!
subroutine hardwired_Sutherlands_coefficients
  !
  !.. Use Statements ..
  use ovar,       only : rhoref,aref,tref
  use thermo_mod, only : species_name
  !
  !======================================================================
  !
  ! COEFFICIENTS FOR SPECIES CONDUCTIVITIES
  !
  !======================================================================
  !
  real(wp), parameter :: CO2_cond(1:3) = [1.66135964E-02_wp, &
                                          3.00000000E+02_wp, &
                                          1.06147831E+03_wp]
  !
  real(wp), parameter :: CO_cond(1:3) = [2.61504985E-02_wp, &
                                         3.00000000E+02_wp, &
                                         1.97055369E+02_wp]
  !
  real(wp), parameter :: H2O_cond(1:3) = [2.50391238E-02_wp, &
                                          3.73200000E+02_wp, &
                                          5.51060270E+03_wp]
  !
  real(wp), parameter :: CH4_cond(1:3) = [3.46243818E-02_wp, &
                                          3.00000000E+02_wp, &
                                          3.06537750E+03_wp]
  !
  real(wp), parameter :: N2_cond(1:3) = [2.61504985E-02_wp, &
                                         3.00000000E+02_wp, &
                                         1.97055369E+02_wp]
  !
  real(wp), parameter :: O2_cond(1:3) = [2.65963673E-02_wp, &
                                         3.00000000E+02_wp, &
                                         2.53134992E+02_wp]
  !
  real(wp), parameter :: H_cond(1:3) = [4.40462538E-01_wp, &
                                        1.00000000E+03_wp, &
                                        8.34029004E+02_wp]
  !
  real(wp), parameter :: H2_cond(1:3) = [1.87808749E-01_wp, &
                                         3.00000000E+02_wp, &
                                         1.66395528E+02_wp]
  !
  real(wp), parameter :: O_cond(1:3) = [9.73731485E-02_wp, &
                                        1.00000000E+03_wp, &
                                        6.09170765E+02_wp]
  !
  real(wp), parameter :: OH_cond(1:3) = [2.65963673E-02_wp, &
                                         3.00000000E+02_wp, &
                                         2.53134992E+02_wp]
  !
  !======================================================================
  !
  ! COEFFICIENTS FOR SPECIES VISCOSITIES
  !
  !======================================================================
  !
  real(wp), parameter :: CO2_visc(1:3) = [1.50146066E-05_wp, &
                                          3.00000000E+02_wp, &
                                          2.79396589E+02_wp]
  !
  real(wp), parameter :: CO_visc(1:3) = [1.79059092E-05_wp, &
                                         3.00000000E+02_wp, &
                                         1.36773324E+02_wp]
  !
  real(wp), parameter :: H2O_visc(1:3) = [1.22841825E-05_wp, &
                                          3.73200000E+02_wp, &
                                          1.06739885E+03_wp]
  !
  real(wp), parameter :: CH4_visc(1:3) = [1.12002411E-05_wp, &
                                          3.00000000E+02_wp, &
                                          1.81765271E+02_wp]
  !
  real(wp), parameter :: N2_visc(1:3) = [1.79059092E-05_wp, &
                                         3.00000000E+02_wp, &
                                         1.36773324E+02_wp]
  !
  real(wp), parameter :: O2_visc(1:3) = [2.07429854E-05_wp, &
                                         3.00000000E+02_wp, &
                                         1.49774296E+02_wp]
  !
  real(wp), parameter :: H_visc(1:3) = [1.42353809E-05_wp, &
                                        1.00000000E+03_wp, &
                                        8.35086051E+02_wp]
  !
  real(wp), parameter :: H2_visc(1:3) = [8.97383478E-06_wp, &
                                         3.00000000E+02_wp, &
                                         1.30667324E+02_wp]
  !
  real(wp), parameter :: O_visc(1:3) = [4.99655681E-05_wp, &
                                        1.00000000E+03_wp, &
                                        6.09169837E+02_wp]
  !
  real(wp), parameter :: OH_visc(1:3) = [2.07429854E-05_wp, &
                                         3.00000000E+02_wp, &
                                         1.49774296E+02_wp]
  !
  integer  :: i
  real(wp) :: mu_suth,k_suth,t_suth,s_suth
  real(wp) :: tref_inv,viscref_inv,condref_inv
  !
  character(len=len(species_name)) :: species
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "hardwired_Sutherlands_coefficients"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  tref_inv = one / tref
  viscref_inv = one / rhoref / aref
  condref_inv = tref / (rhoref*aref*aref*aref)
  !
  do i = 1,ncs
    !
    species = uppercase( trim(adjustl(species_name(i))) )
    !
    select case (trim(adjustl(species)))
      !
      case ("CO2")
        !
        mu_suth = CO2_visc(1) * viscref_inv
        t_suth  = CO2_visc(2) * tref_inv
        s_suth  = CO2_visc(3) * tref_inv
        !
        visc_coef(1,i) = mu_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        visc_coef(2,i) = s_suth
        !
        k_suth = CO2_cond(1) * condref_inv
        t_suth = CO2_cond(2) * tref_inv
        s_suth = CO2_cond(3) * tref_inv
        !
        cond_coef(1,i) = k_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        cond_coef(2,i) = s_suth
        !
      case ("CO")
        !
        mu_suth = CO_visc(1) * viscref_inv
        t_suth  = CO_visc(2) * tref_inv
        s_suth  = CO_visc(3) * tref_inv
        !
        visc_coef(1,i) = mu_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        visc_coef(2,i) = s_suth
        !
        k_suth = CO_cond(1) * condref_inv
        t_suth = CO_cond(2) * tref_inv
        s_suth = CO_cond(3) * tref_inv
        !
        cond_coef(1,i) = k_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        cond_coef(2,i) = s_suth
        !
      case ("H2O")
        !
        mu_suth = H2O_visc(1) * viscref_inv
        t_suth  = H2O_visc(2) * tref_inv
        s_suth  = H2O_visc(3) * tref_inv
        !
        visc_coef(1,i) = mu_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        visc_coef(2,i) = s_suth
        !
        k_suth = H2O_cond(1) * condref_inv
        t_suth = H2O_cond(2) * tref_inv
        s_suth = H2O_cond(3) * tref_inv
        !
        cond_coef(1,i) = k_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        cond_coef(2,i) = s_suth
        !
      case ("CH4")
        !
        mu_suth = CH4_visc(1) * viscref_inv
        t_suth  = CH4_visc(2) * tref_inv
        s_suth  = CH4_visc(3) * tref_inv
        !
        visc_coef(1,i) = mu_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        visc_coef(2,i) = s_suth
        !
        k_suth = CH4_cond(1) * condref_inv
        t_suth = CH4_cond(2) * tref_inv
        s_suth = CH4_cond(3) * tref_inv
        !
        cond_coef(1,i) = k_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        cond_coef(2,i) = s_suth
        !
      case ("N2")
        !
        mu_suth = N2_visc(1) * viscref_inv
        t_suth  = N2_visc(2) * tref_inv
        s_suth  = N2_visc(3) * tref_inv
        !
        visc_coef(1,i) = mu_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        visc_coef(2,i) = s_suth
        !
        k_suth = N2_cond(1) * condref_inv
        t_suth = N2_cond(2) * tref_inv
        s_suth = N2_cond(3) * tref_inv
        !
        cond_coef(1,i) = k_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        cond_coef(2,i) = s_suth
        !
      case ("O2")
        !
        mu_suth = O2_visc(1) * viscref_inv
        t_suth  = O2_visc(2) * tref_inv
        s_suth  = O2_visc(3) * tref_inv
        !
        visc_coef(1,i) = mu_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        visc_coef(2,i) = s_suth
        !
        k_suth = O2_cond(1) * condref_inv
        t_suth = O2_cond(2) * tref_inv
        s_suth = O2_cond(3) * tref_inv
        !
        cond_coef(1,i) = k_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        cond_coef(2,i) = s_suth
        !
      case ("H")
        !
        mu_suth = H_visc(1) * viscref_inv
        t_suth  = H_visc(2) * tref_inv
        s_suth  = H_visc(3) * tref_inv
        !
        visc_coef(1,i) = mu_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        visc_coef(2,i) = s_suth
        !
        k_suth = H_cond(1) * condref_inv
        t_suth = H_cond(2) * tref_inv
        s_suth = H_cond(3) * tref_inv
        !
        cond_coef(1,i) = k_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        cond_coef(2,i) = s_suth
        !
      case ("H2")
        !
        mu_suth = H2_visc(1) * viscref_inv
        t_suth  = H2_visc(2) * tref_inv
        s_suth  = H2_visc(3) * tref_inv
        !
        visc_coef(1,i) = mu_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        visc_coef(2,i) = s_suth
        !
        k_suth = H2_cond(1) * condref_inv
        t_suth = H2_cond(2) * tref_inv
        s_suth = H2_cond(3) * tref_inv
        !
        cond_coef(1,i) = k_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        cond_coef(2,i) = s_suth
        !
      case ("O")
        !
        mu_suth = O_visc(1) * viscref_inv
        t_suth  = O_visc(2) * tref_inv
        s_suth  = O_visc(3) * tref_inv
        !
        visc_coef(1,i) = mu_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        visc_coef(2,i) = s_suth
        !
        k_suth = O_cond(1) * condref_inv
        t_suth = O_cond(2) * tref_inv
        s_suth = O_cond(3) * tref_inv
        !
        cond_coef(1,i) = k_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        cond_coef(2,i) = s_suth
        !
      case ("OH")
        !
        mu_suth = OH_visc(1) * viscref_inv
        t_suth  = OH_visc(2) * tref_inv
        s_suth  = OH_visc(3) * tref_inv
        !
        visc_coef(1,i) = mu_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        visc_coef(2,i) = s_suth
        !
        k_suth = OH_cond(1) * condref_inv
        t_suth = OH_cond(2) * tref_inv
        s_suth = OH_cond(3) * tref_inv
        !
        cond_coef(1,i) = k_suth * (one + s_suth/t_suth) / sqrt(t_suth)
        cond_coef(2,i) = s_suth
        !
    end select
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine hardwired_Sutherlands_coefficients
!
!###############################################################################
!
end module transport_mod
