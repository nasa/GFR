module derivatives_mod
  !
!!!#define PROFILE_ON
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  use matrix_types_mod, only : old_matrix_t
  use matrix_types_mod, only : matrix
  use matrix_types_mod, only : create_matrix_sparsity
  use matrix_types_mod, only : output_sparsity
  use matrix_types_mod, only : output_matrix
  !
  implicit none
  !
  private
  !
  !
  !.. Public Procedures ..
  !
  public :: init_derivative_matrices
  public :: derivatives_memory_usage
  public :: move_deriv_to_device
  public :: move_dqmat_to_device
  !
  !
  !
  ! derivative_matrices : derived type that stores the derivative matrices of
  !                       each local coordinate direction for a given
  !                       combination of cell geometry and cell order
  !
  type, public :: derivative_matrices_t
    type(old_matrix_t), allocatable :: d(:)
  contains
    private
    procedure :: grad_p1 => gradient_for_point_rank1
    procedure :: grad_p2 => gradient_for_point_rank2
    procedure :: grad_c1 => gradient_for_cell_rank1
    procedure :: grad_c2 => gradient_for_cell_rank2
    procedure, public  :: div => divergence_derivative_matrices
    generic,   public  :: grad => grad_p1,grad_p2,grad_c1,grad_c2
  end type derivative_matrices_t
  !
  type, public :: derivative_matrices
    type(matrix), allocatable :: d(:)
  end type derivative_matrices
  !
  ! deriv : derivative matrix for all combinations of cell geometry and
  !        and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The derivative matrix (giving the contribution of each solution
  !         point within cell 'n' to the solution point 'k' that is currently
  !         of interest) is given by the expression
  !
  !         deriv(this_geom,this_order)%d(l)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell(n)%geom
  !        this_order = cell(n)%order
  !        l = the local coordinate direction within the reference cell,
  !            i.e., (1,2,3) for (xi,eta,zeta) or (r,s,t)
  !        k = the solution point within cell 'n' that is currently of interest
  !
  type(derivative_matrices), public, save, allocatable :: deriv(:,:)
  !
  ! dqmat : same as deriv but gives the derivative matrices for the quadrature
  !         points
  !
  type(derivative_matrices), public, save, allocatable :: dqmat(:,:)
  !
  !.. Module Parameter Constants ..
  !
  ! Coefficient to multiply with the unstructured element derivative
  ! matrices. I believe the coefficient should be 1 if using bi-unit
  ! elements, or 2 if using unit elements.
  real(mp), parameter :: coef_dmat = real(one,kind=kind(coef_dmat))
 !real(mp), parameter :: coef_dmat = real(two,kind=kind(coef_dmat))
  !
contains
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function divergence_derivative_matrices(this,f) result(return_value)
  !
  !.. Formal Arguments ..
  class(derivative_matrices_t), intent(in) :: this
  real(wp),                     intent(in) :: f(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(f,dim=1))
  !
  !.. Local Scalars ..
  integer :: n,l
  !
continue
  !
  return_value = zero
  !
  do l = 1,size(this%d)
    do n = 1,size(f,dim=1)
      ! Inner loop over n means:
      !    we reuse f(:,l) and d(l) as much as possible
      !    once they are loaded into cache
      ! instead of:
      !    for n, loading f(:,l) and d(l) into cache and using it, loading
      !    f(:,l+1) and d(l+1) into cache and using it, and repeating again
      !    if size(this%d) == 3; start over for n+1 and continue until
      !    n == size(f,dim=1).
      return_value(n) = return_value(n) + this%d(l)%dot( n , f(:,l) )
    end do
  end do
  !
end function divergence_derivative_matrices
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function gradient_for_point_rank1(this,n,f) result(return_value)
  !
  !.. Formal Arguments ..
  class(derivative_matrices_t), intent(in) :: this
  integer,                      intent(in) :: n
  real(wp),                     intent(in) :: f(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(this%d))
  !
  !.. Local Scalars ..
  integer :: l
  !
continue
  !
  return_value = zero
  !
  do l = 1,size(this%d)
    return_value(l) = this%d(l)%dot( n , f(:) )
  end do
  !
end function gradient_for_point_rank1
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function gradient_for_point_rank2(this,n,f) result(return_value)
  !
  !.. Formal Arguments ..
  class(derivative_matrices_t), intent(in) :: this
  integer,                      intent(in) :: n
  real(wp),                     intent(in) :: f(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(this%d),1:size(f,dim=2))
  !
  !.. Local Scalars ..
  integer :: l,m
  !
continue
  !
  return_value = zero
  !
  do m = 1,size(f,dim=2)
    do l = 1,size(this%d)
      return_value(l,m) = this%d(l)%dot( n , f(:,m) )
    end do
  end do
  !
end function gradient_for_point_rank2
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function gradient_for_cell_rank1(this,f) result(return_value)
  !
  !.. Formal Arguments ..
  class(derivative_matrices_t), intent(in) :: this
  real(wp),                     intent(in) :: f(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(f,dim=1),1:size(this%d))
  !
  !.. Local Scalars ..
  integer :: l
  !
continue
  !
  return_value = zero
  !
  do l = 1,size(this%d) ! Derivative reference direction
    return_value(:,l) = this%d(l)%vm_mult( f(:) )
  end do
  !
end function gradient_for_cell_rank1
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function gradient_for_cell_rank2(this,f) result(return_value)
  !
  !.. Formal Arguments ..
  class(derivative_matrices_t), intent(in) :: this
  real(wp),                     intent(in) :: f(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(f,dim=1),1:size(f,dim=2),1:size(this%d))
  !
  !.. Local Scalars ..
  integer :: l,m
  !
continue
  !
  return_value = zero
  !
  do l = 1,size(this%d) ! Derivative reference direction
    do m = 1,size(f,dim=2) ! variable to differentiate in direction 'l'
      return_value(:,m,l) = this%d(l)%vm_mult( f(:,m) )
    end do
  end do
  !
end function gradient_for_cell_rank2
!
!###############################################################################
!
subroutine init_derivative_matrices()
  !
  !.. Use Statements ..
  use order_mod,        only : n_min_geom,n_max_geom
  use order_mod,        only : n_min_order,n_max_order
  use order_mod,        only : n_order
  use ovar,             only : loc_solution_pts
  use quadrature_mod,   only : geom_is_used
  use element_type_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,this_order,ndim
  integer :: n,n1,n2,i
  character(len=4)   :: gname
  character(len=200) :: array_name
  character(len=300) :: timer_message
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_derivative_matrices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  ! Allocate the vand array
  !
  allocate ( deriv(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"deriv",1,__LINE__,__FILE__,ierr,error_message)
  !
  call check_for_bad_number(__LINE__,__FILE__)
  geom_loop: do this_geom = gmin,gmax
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_Valid)) cycle geom_loop
    if (this_geom == Geom_Node) cycle geom_loop
    !
    ! Get the number of dimensions for this cell type
    !
    ndim = geom_dimen(this_geom)
    !
    order_loop: do this_order = omin,omax
      write (timer_message,100) trim(adjustl(Geom_Name(this_geom))),this_order
      100 format ("Geometry = ",a,"  ;  Order = ",i0)
      call debug_timer(start_timer,timer_message)
      !
      ! Allocate the number of dimensional derivative matrices
      !
      allocate ( deriv(this_geom,this_order)%d(1:ndim) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                       error_message,skip_alloc_pause)
      !
      ! Compute the derivative matrices depending on geometry type
      !
      select case (this_geom)
        case (Geom_Edge)
          !
          deriv(this_geom,this_order)%d(:) = Dmatrix_Edge( this_order , &
                                                           loc_solution_pts )
#ifdef SPECIAL_OUTPUT_SPARSITY
          write (array_name,10) "edge","x",this_order
          call output_sparsity( deriv(this_geom,this_order)%d(1) , array_name )
#endif
          !
          error_message = "dx Derivative"
          call output_matrix(deriv(this_geom,this_order)%d(1), &
                             error_message,this_geom,this_order)
          !
        case (Geom_Quad)
          !
          deriv(this_geom,this_order)%d(:) = &
                            Dmatrix_Quad( deriv(Geom_Edge,this_order)%d(1)%mat )
#ifdef SPECIAL_OUTPUT_SPARSITY
          write (array_name,10) "quad","x",this_order
          call output_sparsity( deriv(this_geom,this_order)%d(1) , array_name )
          write (array_name,10) "quad","y",this_order
          call output_sparsity( deriv(this_geom,this_order)%d(2) , array_name )
#endif
          !
          error_message = "dx Derivative"
          call output_matrix(deriv(this_geom,this_order)%d(1), &
                             error_message,this_geom,this_order)
          !
          error_message = "dy Derivative"
          call output_matrix(deriv(this_geom,this_order)%d(2), &
                             error_message,this_geom,this_order)
          !
        case (Geom_Hexa)
          !
          deriv(this_geom,this_order)%d(:) = &
                            Dmatrix_Hexa( deriv(Geom_Edge,this_order)%d(1)%mat )
#ifdef SPECIAL_OUTPUT_SPARSITY
          write (array_name,10) "hexa","x",this_order
          call output_sparsity( deriv(this_geom,this_order)%d(1) , array_name )
          write (array_name,10) "hexa","y",this_order
          call output_sparsity( deriv(this_geom,this_order)%d(2) , array_name )
          write (array_name,10) "hexa","z",this_order
          call output_sparsity( deriv(this_geom,this_order)%d(3) , array_name )
#endif
          !
          error_message = "dx Derivative"
          call output_matrix(deriv(this_geom,this_order)%d(1), &
                             error_message,this_geom,this_order)
          !
          error_message = "dy Derivative"
          call output_matrix(deriv(this_geom,this_order)%d(2), &
                             error_message,this_geom,this_order)
          !
          error_message = "dz Derivative"
          call output_matrix(deriv(this_geom,this_order)%d(3), &
                             error_message,this_geom,this_order)
          !
        case (Geom_Tria,Geom_Tetr,Geom_Pyra,Geom_Pris)
          !
          deriv(this_geom,this_order)%d(:) = &
                            Dmatrix_Unstructured(this_geom,this_order)
#ifdef SPECIAL_OUTPUT_SPARSITY
          select case (this_geom)
            case (Geom_Tria)
              gname = "tria"
            case (Geom_Tetr)
              gname = "tetr"
            case (Geom_Tetr)
              gname = "pyra"
            case (Geom_Tetr)
              gname = "pris"
            case default
              gname = "unst"
          end select
          write (array_name,10) gname,"x",this_order
          call output_sparsity( deriv(this_geom,this_order)%d(1) , array_name )
          write (array_name,10) gname,"y",this_order
          call output_sparsity( deriv(this_geom,this_order)%d(2) , array_name )
          if (this_geom /= Geom_Tria) then
            write (array_name,10) gname,"z",this_order
            call output_sparsity( deriv(this_geom,this_order)%d(3) , array_name )
          end if
#endif
          !
        case default
          !
          write (error_message,2)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
      end select
      !
      call debug_timer(stop_timer,timer_message)
      call check_for_bad_number(__LINE__,__FILE__)
      !
    end do order_loop
    !
  end do geom_loop
  !
  ! Check for sparsity
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(deriv(this_geom,this_order)%d)) then
        write (timer_message,101) trim(adjustl(Geom_Name(this_geom))),this_order
        101 format ("Sparsity check:  Geometry = ",a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
       !call deriv(this_geom,this_order)%d%check_for_sparsity
        call create_matrix_sparsity( deriv(this_geom,this_order)%d )
        call debug_timer(stop_timer,timer_message)
        call check_for_bad_number(__LINE__,__FILE__)
      end if
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("deriv(",a,",",i0,")%d")
  2 format (" A grid cell of an unknown geometry type was found!")
  !
  10 format (a4,".deriv.d",a1,".p",i0,".tec.dat")
  !
end subroutine init_derivative_matrices
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function Dmatrix_Edge(n,location) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : derivative_matrix_legendre_gauss
  use polynomial_mod, only : derivative_matrix_legendre_gauss_lobatto
  !
  !.. Formal Arguments ..
  integer, intent(in) :: n
  integer, intent(in) :: location
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: nsp,ierr
  !
  !.. Local Arrays ..
  real(mp) :: dmat(1:n+1,1:n+1)
  !
continue
  !
  nsp = n+1
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:1) , stat=ierr )
#ifdef DISABLE_PURE
  call check_for_bad_number(__LINE__,__FILE__)
#endif
  !
  ! Compute the derivative matrix using quad precision
  !
  select case (location)
    !
    case (Legendre_Gauss)
      !
      ! Derivative matrix for Legendre Gauss solution points
      !
      dmat(:,:) = derivative_matrix_legendre_gauss(n)
#ifdef DISABLE_PURE
      call check_for_bad_number(__LINE__,__FILE__)
#endif
      !
    case (Legendre_Gauss_Lobatto)
      !
      ! Derivative matrix for Legendre Gauss Lobatto solution points
      !
      dmat(:,:) = derivative_matrix_legendre_gauss_lobatto(n)
#ifdef DISABLE_PURE
      call check_for_bad_number(__LINE__,__FILE__)
#endif
      !
    case default
      !
      ! Default to Legendre Gauss solution points
      !
      dmat(:,:) = derivative_matrix_legendre_gauss(n)
#ifdef DISABLE_PURE
      call check_for_bad_number(__LINE__,__FILE__)
#endif
      !
  end select
  !
  !!!!!
  !!!!! Alternative Method if this is ever needed
  !!!!!
  !!!!return_value = invert_matrix(V)
  !!!!return_value = matmul(Vr,return_value)
  !
  ! We need to transpose this derivative matrix to make its use more efficient
  !
  dmat(:,:) = transpose( dmat )
#ifdef DISABLE_PURE
  call check_for_bad_number(__LINE__,__FILE__)
#endif
  !
  ! Use the chop function to copy dmat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value(1)%mat = chop( dmat ) ! F2003 AUTO-REALLOC
#ifdef DISABLE_PURE
  call check_for_bad_number(__LINE__,__FILE__)
#endif
  !
end function Dmatrix_Edge
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function Dmatrix_Quad(d1_mat) result(return_value)
  !
  ! Create the derivative matrices for quad elements.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: d1_mat(:,:)
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: i,j,np,nsp,ierr
  !
  !.. Local Arrays ..
  integer :: idx(1:size(d1_mat,dim=1))
  !
continue
  !
  ! ###########################################################################
  ! ###########################################################################
  ! ####  NOTE: No need to use chopping function in this procedure. The    ####
  ! ####        quad derivative matrix is really just the 1D derivative    ####
  ! ####        matrix applied in each direction of the standard quad      ####
  ! ####        element. The 1D derivative matrix is obtained by the       ####
  ! ####        function Dmatrix_Edge which applies the chopping function  ####
  ! ####        to its result so there is no need to repeat it again.      ####
  ! ###########################################################################
  ! ###########################################################################
  !
  np = size(d1_mat,dim=1)
  nsp = np*np
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:2) , stat=ierr )
  !
  ! Allocate each individual derivative matrix within return_value
  !
  allocate ( return_value(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the quadrilateral derivative matrix in the x/xi direction
  !
  do j = 1,np
    !
    ! Create the indices within the derivative matrix in the x/xi direction
    ! that correspond to the locations within the 1D derivative matrix
    !
    idx(:) = [ ( (j-1)*np+i , i=1,np ) ]
    !
    ! Copy the 1D derivative matrix to the corresponding locations given by idx
    ! NOTE: We need to copy the transpose due to the
    !       way the derivative matrices are stored
    !
   !return_value(1)%mat(idx,idx) = transpose( d1_mat )
    return_value(1)%mat(idx,idx) = d1_mat
    !
  end do
  !
  ! Create the quadrilateral derivative matrix in the y/eta direction
  !
  do i = 1,np
    !
    ! Create the indices within the derivative matrix in the y/eta direction
    ! that correspond to the locations within the 1D derivative matrix
    !
    idx(:) = [ ( (j-1)*np+i , j=1,np ) ]
    !
    ! Copy the 1D derivative matrix to the corresponding locations given by idx
    ! NOTE: We need to copy the transpose due to the
    !       way the derivative matrices are stored
    !
   !return_value(2)%mat(idx,idx) = transpose( d1_mat )
    return_value(2)%mat(idx,idx) = d1_mat
    !
  end do
  !
end function Dmatrix_Quad
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function Dmatrix_Hexa(d1_mat) result(return_value)
  !
  ! Create the derivative matrices for hexahedral elements.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: d1_mat(:,:)
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,np,nsp,ierr
  !
  !.. Local Arrays ..
  integer :: idx(1:size(d1_mat,dim=1))
  !
continue
  !
  ! ###########################################################################
  ! ###########################################################################
  ! ####  NOTE: No need to use chopping function in this procedure. The    ####
  ! ####        hexa derivative matrix is really just the 1D derivative    ####
  ! ####        matrix applied in each direction of the standard hexa      ####
  ! ####        element. The 1D derivative matrix is obtained by the       ####
  ! ####        function Dmatrix_Edge which applies the chopping function  ####
  ! ####        to its result so there is no need to repeat it again.      ####
  ! ###########################################################################
  ! ###########################################################################
  !
  np = size(d1_mat,dim=1)
  nsp = np*np*np
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:3) , stat=ierr )
  !
  ! Allocate each individual derivative matrix within return_value
  !
  allocate ( return_value(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(3)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the hexahedral derivative matrix in the x/xi direction
  !
  do k = 1,np
    do j = 1,np
      !
      ! Create the indices within the derivative matrix in the x/xi direction
      ! that correspond to the locations within the 1D derivative matrix
      !
      idx = [ ( (k-1)*np*np + (j-1)*np + i , i=1,np ) ]
      !
      ! Copy the 1D derivative matrix to the
      ! corresponding locations given by idx
      ! NOTE: We need to copy the transpose due to the
      !       way the derivative matrices are stored
      !
      return_value(1)%mat(idx,idx) = d1_mat
      !
    end do
  end do
  !
  ! Create the hexahedral derivative matrix in the y/eta direction
  !
  do k = 1,np
    do i = 1,np
      !
      ! Create the indices within the derivative matrix in the y/eta direction
      ! that correspond to the locations within the 1D derivative matrix
      !
      idx = [ ( (k-1)*np*np + (j-1)*np + i , j=1,np ) ]
      !
      ! Copy the 1D derivative matrix to the
      ! corresponding locations given by idx
      ! NOTE: We need to copy the transpose due to the
      !       way the derivative matrices are stored
      !
      return_value(2)%mat(idx,idx) = d1_mat
      !
    end do
  end do
  !
  ! Create the hexahedral derivative matrix in the z/zeta direction
  !
  do j = 1,np
    do i = 1,np
      !
      ! Create the indices within the derivative matrix in the z/zeta direction
      ! that correspond to the locations within the 1D derivative matrix
      !
      idx = [ ( (k-1)*np*np + (j-1)*np + i , k=1,np ) ]
      !
      ! Copy the 1D derivative matrix to the
      ! corresponding locations given by idx
      ! NOTE: We need to copy the transpose due to the
      !       way the derivative matrices are stored
      !
      return_value(3)%mat(idx,idx) = d1_mat
      !
    end do
  end do
  !
end function Dmatrix_Hexa
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function Dmatrix_Unstructured(this_geom,this_order) result(return_value)
  !
  !.. Use Statements ..
  use element_type_mod, only : std_elem
  use vandermonde_mod,  only : vand
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_geom
  integer, intent(in) :: this_order
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: n,ndim,ierr
  !
  !.. Local Allocatable Arrays ..
  type(matrix), allocatable :: grad_vand(:)
  real(mp),     allocatable :: dmat(:,:)
  real(mp),     allocatable :: dvand(:,:)
  real(mp),     allocatable :: vand_inv(:,:)
  !
continue
  !
  ndim = geom_dimen(this_geom)
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:ndim) , stat=ierr )
  !
  ! Compute the derivative of the Vandermonde matrix
  !
  select case (this_geom)
    case (Geom_Tria)
      ! F2003 AUTO-REALLOC
      grad_vand = GradVandermonde_Tria(std_elem(this_geom,this_order)%pts)
    case (Geom_Tetr)
      ! F2003 AUTO-REALLOC
      grad_vand = GradVandermonde_Tetr(std_elem(this_geom,this_order)%pts)
    case (Geom_Pyra)
      ! F2003 AUTO-REALLOC
      grad_vand = GradVandermonde_Pyra(std_elem(this_geom,this_order)%pts)
    case (Geom_Pris)
      ! F2003 AUTO-REALLOC
      grad_vand = GradVandermonde_Pris(std_elem(this_geom,this_order)%pts)
  end select
  !
  ! Create a quad precision copy of the inverse of the Vandermonde matrix
  !
  ! F2003 AUTO-REALLOC
  vand_inv = real( vand(this_geom,this_order)%nodal2modal%mat , &
                      kind=kind(vand_inv) )
  !
  ! Loop to create the derivative matrix for each coordinate direction
  ! [ r/x/xi , s/y/eta ]            for 2D (tria)
  ! [ r/x/xi , s/y/eta , t/z/zeta ] for 3D (tetr,pyra,pris)
  !
  do n = 1,ndim
    !
    ! Compute the derivative matrix in the r/x/xi direction
    !
    dvand = real( grad_vand(n)%mat , kind=kind(dvand) ) ! F2003 AUTO-REALLOC
    dmat = coef_dmat * matmul( dvand , vand_inv )       ! F2003 AUTO-REALLOC
    !
    ! Use the chop function to copy dmat into return_value while also
    ! setting any values that are less than working machine epsilon to zero.
    ! NOTE: Transpose this matrix to get it in the right ordering
    !
    return_value(n)%mat = transpose( chop( dmat ) ) ! F2003 AUTO-REALLOC
    !
  end do
  !
end function Dmatrix_Unstructured
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function GradVandermonde_Edge(r) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: r(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(r),1:size(r))
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(mp) :: coef,rj,val
  !
  !.. Local Arrays ..
  real(mp) :: xi(1:size(r))
  real(mp) :: grad_vand(1:size(r),1:size(r))
  !
  !.. Local Parameters ..
  real(mp), parameter :: one = 1.0_mp
  !
continue
  !
  return_value(:,:) = zero
  !
  xi(:) = real(r(:),kind=kind(xi))
  !
  do j = 2,size(xi)
    !
    rj = real(j-1,kind=kind(rj))
    coef = sqrt(rj*(rj+one))
    !
    do i = 1,size(xi)
      val = normalized_jacobi_poly(j-2,one,one,xi(i))
      return_value(i,j) = real(coef*val,kind=kind(return_value))
    end do
    !
  end do
  !
end function GradVandermonde_Edge
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function GradVandermonde_Tria(rs) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  use polynomial_mod, only : grad_normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rs(:,:)
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: i,im1,j,k,m,this_order,nsp,ierr
  real(mp) :: r,s,ri,twoiph,alpha,alphb
  real(mp) :: a,b,fa,gb,dfa,dgb,h1mb,dmdr,dmds
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: dvdr(:,:)
  real(mp), allocatable :: dvds(:,:)
  !
  !.. Local Parameters ..
  real(mp), parameter :: zero = 0.0_mp
  real(mp), parameter :: half = 0.5_mp
  real(mp), parameter :: one  = 1.0_mp
  real(mp), parameter :: two  = 2.0_mp
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = (nint(sqrt(real(9-8*(1-nsp))))-3)/2
  !
  ! Allocate the dimensions of the derivative of the Vandermonde matrices
  !
  allocate ( return_value(1:2) , stat=ierr )
  !
  ! Create the derivatives of the Vandermonde matrices using quad precision
  !
  dvdr = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  dvds = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  !
  m = 1
  do i = 0,this_order
    !
    ri   = real( i , kind=mp )
    im1  = max(i-1,0)
    !
    twoiph = sqrt(two) * two**i
    !
    alpha = zero
    alphb = real( 2*i+1 , kind=mp )
    !
    do j = 0,this_order-i
      !
      do k = 1,nsp
        !
        ! Get a,b from r,s
        !
        r = real( rs(1,k) , kind=mp )
        s = real( rs(2,k) , kind=mp )
        !
        ! NOTE: If s=1, r=-1, therefore (one+r)=zero and reciprocal(one-s)
        !       will return zero, thus a=-1. If we didnt use reciprocal for
        !       1/(one-s), we would have to create logic to prevent against
        !       division by zero when s=1, even though r=-1.
        !
        a = two * (one+r) * reciprocal(one-s) - one
        b = s
        !
        fa = normalized_jacobi_poly(i,alpha,zero,a)
        gb = normalized_jacobi_poly(j,alphb,zero,b)
        !
        dfa = grad_normalized_jacobi_poly(i,alpha,zero,a)
        dgb = grad_normalized_jacobi_poly(j,alphb,zero,b)
        !
        h1mb = half*(one-b)
        !
        dmdr = dfa * gb * h1mb**im1
        dmds = half*dmdr*(one+a) + fa*(dgb*h1mb**i - half*ri*gb*h1mb**im1)
        !
        dvdr(k,m) = twoiph * dmdr
        dvds(k,m) = twoiph * dmds
        !
      end do
      !
      m = m + 1
      !
    end do
  end do
  !
  ! Use the chop function to copy dvdr and dvds into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value(1)%mat = chop( dvdr ) ! F2003 AUTO-REALLOC
  return_value(2)%mat = chop( dvds ) ! F2003 AUTO-REALLOC
  !
end function GradVandermonde_Tria
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function GradVandermonde_Tetr(rst) result(return_value)
  !
  !.. Use Statements ..
  use order_mod,      only : geom_solpts
  use polynomial_mod, only : normalized_jacobi_poly
  use polynomial_mod, only : grad_normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rst(:,:)
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: i,im1,ipjm1,j,k,m,n,this_order,nsp,ierr
  real(mp) :: ri,rj,norm_coef,h1pa,h1mb,h1mc,tmps,tmpt
  real(mp) :: h1mb_to_im1,h1mb_to_i,h1mc_to_ipjm1,h1mc_to_ipj
  real(mp) :: r,alpha,a,fa,dfa,dmdr
  real(mp) :: s,alphb,b,gb,dgb,dmds
  real(mp) :: t,alphc,c,hc,dhc,dmdt
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: dvdr(:,:)
  real(mp), allocatable :: dvds(:,:)
  real(mp), allocatable :: dvdt(:,:)
  !
  !.. Local Parameters ..
  real(mp), parameter :: zero = 0.0_mp
  real(mp), parameter :: half = 0.5_mp
  real(mp), parameter :: one  = 1.0_mp
  real(mp), parameter :: two  = 2.0_mp
  !
continue
  !
  nsp = size(rst,dim=2)
  do n = 0,nsp
    if (geom_solpts(Geom_Tetr,n) >= nsp) then
      this_order = n
      exit
    end if
  end do
  !
  ! Allocate the dimensions of the derivative of the Vandermonde matrices
  !
  allocate ( return_value(1:3) , stat=ierr )
  !
  ! Create the derivatives of the Vandermonde matrices using quad precision
  !
  dvdr = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  dvds = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  dvdt = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  !
  m = 1
  do i = 0,this_order
    !
    ri   = real( i , kind=mp )
    im1  = max(i-1,0)
    !
    do j = 0,this_order-i
      !
      rj = real( j , kind=mp )
      ipjm1 = max(i+j-1,0)
      norm_coef = sqrt(two) * two**(2*i+j+1)
      !
      alpha = zero
      alphb = real( 2*i+1 , kind=mp )
      alphc = real( 2*i+2*j+2 , kind=mp )
      !
      do k = 0,this_order-i-j
        !
        do n = 1,nsp
          !
          ! Get a,b from r,s
          !
          r = real( rst(1,n) , kind=mp )
          s = real( rst(2,n) , kind=mp )
          t = real( rst(3,n) , kind=mp )
          !
          ! NOTE: If s+t=0, r=-1, therefore (one+r)=zero and reciprocal(-s-t)
          !       will return zero, thus a=-1. If we didnt use reciprocal for
          !       1/(-s-t), we would have to create logic to prevent against
          !       division by zero when s+t=0, even though r=-1.
          !
          ! NOTE: If t=1, s=-1, therefore (one+s)=zero and reciprocal(one-t)
          !       will return zero, thus b=-1. If we didnt use reciprocal for
          !       1/(one-t), we would have to create logic to prevent against
          !       division by zero when t=1, even though s=-1.
          !
          a = two * (one+r) * reciprocal( -s-t) - one
          b = two * (one+s) * reciprocal(one-t) - one
          c = t
          !
          fa = normalized_jacobi_poly(i,alpha,zero,a)
          gb = normalized_jacobi_poly(j,alphb,zero,b)
          hc = normalized_jacobi_poly(k,alphc,zero,c)
          !
          dfa = grad_normalized_jacobi_poly(i,alpha,zero,a)
          dgb = grad_normalized_jacobi_poly(j,alphb,zero,b)
          dhc = grad_normalized_jacobi_poly(k,alphc,zero,c)
          !
          h1pa = half*(one+a)
          h1mb = half*(one-b)
          h1mc = half*(one-c)
         !!
         !dmdr = dfa * gb * h1mb**im1
         !dmds = half*dmdr*(one+a) + fa*(dgb*h1mb**i - half*ri*gb*h1mb**im1)
         !dmdt = zero
          !
          h1mb_to_im1 = h1mb**im1
          h1mb_to_i   = h1mb**i
          !
          h1mc_to_ipjm1 = h1mc**ipjm1
          h1mc_to_ipj   = h1mc**(i+j)
          !
          dmdr = dfa*gb*hc * h1mb_to_im1 * h1mc_to_ipjm1
          !
          tmps = fa*hc*h1mc_to_ipjm1 * (dgb*h1mb_to_i - half*ri*gb*h1mb_to_im1)
          dmds = h1pa*dmdr + tmps
          !
          tmpt = fa*gb*h1mb_to_i* (dhc*h1mc_to_ipj - &
                                  half*(ri+rj)*hc*h1mc_to_ipjm1)
          dmdt = h1pa*dmdr + half*(one+b)*tmps + tmpt
          !
          dvdr(n,m) = norm_coef * dmdr
          dvds(n,m) = norm_coef * dmds
          dvdt(n,m) = norm_coef * dmdt
          !
        end do
        !
        m = m + 1
        !
      end do
      !
    end do
  end do
  !
  ! Use the chop function to copy dvdr and dvds into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value(1)%mat = chop( dvdr ) ! F2003 AUTO-REALLOC
  return_value(2)%mat = chop( dvds ) ! F2003 AUTO-REALLOC
  return_value(3)%mat = chop( dvdt ) ! F2003 AUTO-REALLOC
  !
end function GradVandermonde_Tetr
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function GradVandermonde_Pyra(rst) result(return_value)
  !
  !.. Use Statements ..
  use order_mod,      only : geom_solpts
  use polynomial_mod, only : normalized_jacobi_poly
  use polynomial_mod, only : grad_normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rst(:,:)
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: i,im1,j,k,m,n,this_order,nsp,ierr
  real(mp) :: ri,twoiph,beta,h1mb
  real(mp) :: r,alpha,a,fa,dfa,dmdr
  real(mp) :: s,alphb,b,gb,dgb,dmds
  real(mp) :: t,alphc,c,hc,dhc,dmdt
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: dvdr(:,:)
  real(mp), allocatable :: dvds(:,:)
  real(mp), allocatable :: dvdt(:,:)
  !
  !.. Local Parameters ..
  real(mp), parameter :: zero = 0.0_mp
  real(mp), parameter :: half = 0.5_mp
  real(mp), parameter :: one  = 1.0_mp
  real(mp), parameter :: two  = 2.0_mp
  !
continue
  !
  nsp = size(rst,dim=2)
  do n = 0,nsp
    if (geom_solpts(Geom_Pyra,n) >= nsp) then
      this_order = n
      exit
    end if
  end do
  !
  ! Allocate the dimensions of the derivative of the Vandermonde matrices
  !
  allocate ( return_value(1:3) , stat=ierr )
  !
  ! Create the derivatives of the Vandermonde matrices using quad precision
  !
  dvdr = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  dvds = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  dvdt = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  !
  m = 1
  do i = 0,this_order
    !
    ri   = real( i , kind=mp )
    im1  = max(i-1,0)
    !
    twoiph = sqrt(two) * two**i
    !
    alpha = zero
    alphb = real( 2*i+1 , kind=mp )
    beta  = zero
    !
    do j = 0,this_order-i
      !
      alphc = zero
      !
      do k = 0,this_order-i-j
        !
        do n = 1,nsp
          !
          ! Get a,b from r,s
          !
          r = real( rst(1,n) , kind=mp )
          s = real( rst(2,n) , kind=mp )
          t = real( rst(3,n) , kind=mp )
          !
          ! NOTE: If s=1, r=-1, therefore (one+r)=zero and reciprocal(one-s)
          !       will return zero, thus a=-1. If we didnt use reciprocal for
          !       1/(one-s), we would have to create logic to prevent against
          !       division by zero when s=1, even though r=-1.
          !
          a = two * (one+r) * reciprocal(one-s) - one
          b = s
          c = zero
          !
          fa = normalized_jacobi_poly(i,alpha,beta,a)
          gb = normalized_jacobi_poly(j,alphb,beta,b)
          hc = normalized_jacobi_poly(k,alphc,beta,c)
          !
          dfa = grad_normalized_jacobi_poly(i,alpha,beta,a)
          dgb = grad_normalized_jacobi_poly(j,alphb,beta,b)
          dhc = grad_normalized_jacobi_poly(k,alphc,beta,c)
          !
          h1mb = half*(one-b)
          !
          dmdr = dfa * gb * h1mb**im1
          dmds = half*dmdr*(one+a) + fa*(dgb*h1mb**i - half*ri*gb*h1mb**im1)
          dmdt = zero
          !
          dvdr(n,m) = twoiph * dmdr
          dvds(n,m) = twoiph * dmds
          dvdt(n,m) = twoiph * dmdt
          !
        end do
        !
        m = m + 1
        !
      end do
      !
    end do
  end do
  !
  ! Use the chop function to copy dvdr and dvds into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value(1)%mat = chop( dvdr ) ! F2003 AUTO-REALLOC
  return_value(2)%mat = chop( dvds ) ! F2003 AUTO-REALLOC
  return_value(3)%mat = chop( dvdt ) ! F2003 AUTO-REALLOC
  !
end function GradVandermonde_Pyra
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function GradVandermonde_Pris(rst) result(return_value)
  !
  !.. Use Statements ..
  use order_mod,      only : geom_solpts
  use polynomial_mod, only : normalized_jacobi_poly
  use polynomial_mod, only : grad_normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rst(:,:)
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: i,im1,j,k,m,n,this_order,nsp,ierr
  real(mp) :: ri,twoiph,beta,h1mb
  real(mp) :: r,alpha,a,fa,dfa,dmdr
  real(mp) :: s,alphb,b,gb,dgb,dmds
  real(mp) :: t,alphc,c,hc,dhc,dmdt
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: dvdr(:,:)
  real(mp), allocatable :: dvds(:,:)
  real(mp), allocatable :: dvdt(:,:)
  !
  !.. Local Parameters ..
  real(mp), parameter :: zero = 0.0_mp
  real(mp), parameter :: half = 0.5_mp
  real(mp), parameter :: one  = 1.0_mp
  real(mp), parameter :: two  = 2.0_mp
  !
continue
  !
  nsp = size(rst,dim=2)
  do n = 0,nsp
    if (geom_solpts(Geom_Pris,n) >= nsp) then
      this_order = n
      exit
    end if
  end do
  !
  ! Allocate the dimensions of the derivative of the Vandermonde matrices
  !
  allocate ( return_value(1:3) , stat=ierr )
  !
  ! Create the derivatives of the Vandermonde matrices using quad precision
  !
  dvdr = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  dvds = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  dvdt = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  !
  m = 1
  do i = 0,this_order
    !
    ri   = real( i , kind=mp )
    im1  = max(i-1,0)
    !
    twoiph = sqrt(two) * two**i
    !
    alpha = zero
    alphb = real( 2*i+1 , kind=mp )
    beta  = zero
    !
    do j = 0,this_order-i
      !
      alphc = zero
      !
      do k = 0,this_order-i-j
        !
        do n = 1,nsp
          !
          ! Get a,b from r,s
          !
          r = real( rst(1,n) , kind=mp )
          s = real( rst(2,n) , kind=mp )
          t = real( rst(3,n) , kind=mp )
          !
          ! NOTE: If s=1, r=-1, therefore (one+r)=zero and reciprocal(one-s)
          !       will return zero, thus a=-1. If we didnt use reciprocal for
          !       1/(one-s), we would have to create logic to prevent against
          !       division by zero when s=1, even though r=-1.
          !
          a = two * (one+r) * reciprocal(one-s) - one
          b = s
          c = zero
          !
          fa = normalized_jacobi_poly(i,alpha,beta,a)
          gb = normalized_jacobi_poly(j,alphb,beta,b)
          hc = normalized_jacobi_poly(k,alphc,beta,c)
          !
          dfa = grad_normalized_jacobi_poly(i,alpha,beta,a)
          dgb = grad_normalized_jacobi_poly(j,alphb,beta,b)
          dhc = grad_normalized_jacobi_poly(k,alphc,beta,c)
          !
          h1mb = half*(one-b)
          !
          dmdr = dfa * gb * h1mb**im1
          dmds = half*dmdr*(one+a) + fa*(dgb*h1mb**i - half*ri*gb*h1mb**im1)
          dmdt = zero
          !
          dvdr(n,m) = twoiph * dmdr
          dvds(n,m) = twoiph * dmds
          dvdt(n,m) = twoiph * dmdt
          !
        end do
        !
        m = m + 1
        !
      end do
      !
    end do
  end do
  !
  ! Use the chop function to copy dvdr and dvds into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value(1)%mat = chop( dvdr ) ! F2003 AUTO-REALLOC
  return_value(2)%mat = chop( dvds ) ! F2003 AUTO-REALLOC
  return_value(3)%mat = chop( dvdt ) ! F2003 AUTO-REALLOC
  !
end function GradVandermonde_Pris
!
!###############################################################################
!
subroutine derivatives_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
#ifndef DISABLE_DTIO
  use matrix_types_mod, only : matrix_memory_usage
#endif
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou,i,j,k
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: total
  !
#ifndef DISABLE_DTIO
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  ! Check the memory of the arrays within this module
  !
  if (allocated(deriv)) then
    !
    call dt%set( storage_size(deriv,kind=inttype) , &
                         size(deriv,kind=inttype) )
    !
    do j = lbound(deriv,dim=2),ubound(deriv,dim=2)
      do i = lbound(deriv,dim=1),ubound(deriv,dim=1)
        !
        if (allocated(deriv(i,j)%d)) then
          !
          call array%set( storage_size(deriv(i,j)%d,kind=inttype) , &
                                  size(deriv(i,j)%d,kind=inttype) )
          !
          write (array_name,4) "deriv",Geom_Name(i),j
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(deriv(i,j)%d,dim=1),ubound(deriv(i,j)%d,dim=1)
            !
            array = matrix_memory_usage( deriv(i,j)%d(k) )
            !
            write (array_name,4) "deriv",Geom_Name(i),j,k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "deriv"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "deriv"
    !
  end if
  !
  if (allocated(dqmat)) then
    !
    call dt%set( storage_size(dqmat,kind=inttype) , &
                         size(dqmat,kind=inttype) )
    !
    do j = lbound(dqmat,dim=2),ubound(dqmat,dim=2)
      do i = lbound(dqmat,dim=1),ubound(dqmat,dim=1)
        !
        if (allocated(dqmat(i,j)%d)) then
          !
          call array%set( storage_size(dqmat(i,j)%d,kind=inttype) , &
                                  size(dqmat(i,j)%d,kind=inttype) )
          !
          write (array_name,4) "dqmat",Geom_Name(i),j
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(dqmat(i,j)%d,dim=1),ubound(dqmat(i,j)%d,dim=1)
            !
            array = matrix_memory_usage( dqmat(i,j)%d(k) )
            !
            write (array_name,4) "dqmat",Geom_Name(i),j,k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "dqmat"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "dqmat"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module derivatives_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",a,",",i0,")%d",:,"(",i0,")%mat")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine derivatives_memory_usage
!
!###############################################################################
!
subroutine move_deriv_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,k,k1,k2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_deriv_to_device"
  !
continue
  !
  if (allocated(deriv)) then
    !$ACC ENTER DATA COPYIN(deriv)
    i1 = lbound(deriv,dim=1)
    i2 = ubound(deriv,dim=1)
    j1 = lbound(deriv,dim=2)
    j2 = ubound(deriv,dim=2)
    !$ACC ENTER DATA COPYIN(deriv(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        if (allocated(deriv(i,j)%d)) then
          !$ACC ENTER DATA COPYIN(deriv(i,j)%d)
          k1 = lbound(deriv(i,j)%d,dim=1)
          k2 = ubound(deriv(i,j)%d,dim=1)
          !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k1:k2))
          do k = k1,k2
            !
            if (allocated(deriv(i,j)%d(k)%mat)) then
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%mat)
              m1 = lbound(deriv(i,j)%d(k)%mat,dim=1)
              m2 = ubound(deriv(i,j)%d(k)%mat,dim=1)
              n1 = lbound(deriv(i,j)%d(k)%mat,dim=2)
              n2 = ubound(deriv(i,j)%d(k)%mat,dim=2)
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(deriv(i,j)%d(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%row_ptr)
              m1 = lbound(deriv(i,j)%d(k)%row_ptr,dim=1)
              m2 = ubound(deriv(i,j)%d(k)%row_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(deriv(i,j)%d(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%row_idx)
              m1 = lbound(deriv(i,j)%d(k)%row_idx,dim=1)
              m2 = ubound(deriv(i,j)%d(k)%row_idx,dim=1)
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(deriv(i,j)%d(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%row_val)
              m1 = lbound(deriv(i,j)%d(k)%row_val,dim=1)
              m2 = ubound(deriv(i,j)%d(k)%row_val,dim=1)
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%row_val(m1:m2))
            end if
            !
            if (allocated(deriv(i,j)%d(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%col_ptr)
              m1 = lbound(deriv(i,j)%d(k)%col_ptr,dim=1)
              m2 = ubound(deriv(i,j)%d(k)%col_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(deriv(i,j)%d(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%col_idx)
              m1 = lbound(deriv(i,j)%d(k)%col_idx,dim=1)
              m2 = ubound(deriv(i,j)%d(k)%col_idx,dim=1)
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(deriv(i,j)%d(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%col_val)
              m1 = lbound(deriv(i,j)%d(k)%col_val,dim=1)
              m2 = ubound(deriv(i,j)%d(k)%col_val,dim=1)
              !$ACC ENTER DATA COPYIN(deriv(i,j)%d(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
      end do
    end do
  end if
  !
end subroutine move_deriv_to_device
!
!###############################################################################
!
subroutine move_dqmat_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,k,k1,k2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_dqmat_to_device"
  !
continue
  !
  if (allocated(dqmat)) then
    !$ACC ENTER DATA COPYIN(dqmat)
    i1 = lbound(dqmat,dim=1)
    i2 = ubound(dqmat,dim=1)
    j1 = lbound(dqmat,dim=2)
    j2 = ubound(dqmat,dim=2)
    !$ACC ENTER DATA COPYIN(dqmat(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        if (allocated(dqmat(i,j)%d)) then
          !$ACC ENTER DATA COPYIN(dqmat(i,j)%d)
          k1 = lbound(dqmat(i,j)%d,dim=1)
          k2 = ubound(dqmat(i,j)%d,dim=1)
          !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k1:k2))
          do k = k1,k2
            !
            if (allocated(dqmat(i,j)%d(k)%mat)) then
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%mat)
              m1 = lbound(dqmat(i,j)%d(k)%mat,dim=1)
              m2 = ubound(dqmat(i,j)%d(k)%mat,dim=1)
              n1 = lbound(dqmat(i,j)%d(k)%mat,dim=2)
              n2 = ubound(dqmat(i,j)%d(k)%mat,dim=2)
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(dqmat(i,j)%d(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%row_ptr)
              m1 = lbound(dqmat(i,j)%d(k)%row_ptr,dim=1)
              m2 = ubound(dqmat(i,j)%d(k)%row_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(dqmat(i,j)%d(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%row_idx)
              m1 = lbound(dqmat(i,j)%d(k)%row_idx,dim=1)
              m2 = ubound(dqmat(i,j)%d(k)%row_idx,dim=1)
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(dqmat(i,j)%d(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%row_val)
              m1 = lbound(dqmat(i,j)%d(k)%row_val,dim=1)
              m2 = ubound(dqmat(i,j)%d(k)%row_val,dim=1)
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%row_val(m1:m2))
            end if
            !
            if (allocated(dqmat(i,j)%d(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%col_ptr)
              m1 = lbound(dqmat(i,j)%d(k)%col_ptr,dim=1)
              m2 = ubound(dqmat(i,j)%d(k)%col_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(dqmat(i,j)%d(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%col_idx)
              m1 = lbound(dqmat(i,j)%d(k)%col_idx,dim=1)
              m2 = ubound(dqmat(i,j)%d(k)%col_idx,dim=1)
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(dqmat(i,j)%d(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%col_val)
              m1 = lbound(dqmat(i,j)%d(k)%col_val,dim=1)
              m2 = ubound(dqmat(i,j)%d(k)%col_val,dim=1)
              !$ACC ENTER DATA COPYIN(dqmat(i,j)%d(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
      end do
    end do
  end if
  !
end subroutine move_dqmat_to_device
!
!###############################################################################
!
#include "Functions/check_for_bad_number.f90"
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module derivatives_mod
!
!###############################################################################
!
