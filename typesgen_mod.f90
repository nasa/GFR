module generic_types_mod
  !
  use module_kind_types, only : wp,zero,dot,intseq
  !
  implicit none
  !
  private
  !
  ! matrix : derived type that stores a matrix quantity
  !            Examples : interpolation matrix, derivative matrix, etc.
  !
  type, public :: matrix
    !
    real(wp), allocatable :: mat(:,:)
    !
    integer,  allocatable :: row_ptr(:)
    integer,  allocatable :: row_idx(:)
    real(wp), allocatable :: row_val(:)
    !
    integer,  allocatable :: col_ptr(:)
    integer,  allocatable :: col_idx(:)
    real(wp), allocatable :: col_val(:)
    !
    !
    procedure(dot_function), pointer :: dot
    procedure(mat_vec_mult), pointer :: mv_mult
    procedure(vec_mat_mult), pointer :: vm_mult
    !
  contains
    !
    procedure :: check_for_sparsity
    !
  end type matrix
  !
  !
  !
  abstract interface
    !
    function dot_function(this,idx,vector) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      integer, intent(in) :: idx
      real(wp), dimension(:), intent(in) :: vector
      real(wp) :: return_value
    end function dot_function
    !
    function mat_vec_mult(this,vector) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      real(wp), dimension(:), intent(in) :: vector
      real(wp), dimension(1:size(this%mat,dim=1)) :: return_value
    end function mat_vec_mult
    !
    function vec_mat_mult(this,vector) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      real(wp), dimension(:), intent(in) :: vector
      real(wp), dimension(1:size(this%mat,dim=2)) :: return_value
    end function vec_mat_mult
    !
  end interface
  !
contains
!
!###############################################################################
!
elemental subroutine check_for_sparsity(this)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(inout) :: this
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,nsz_dim1,nsz_dim2,total_nz,ierr
  real(wp) :: sparsity
  !
  !.. Local Parameters ..
  real(wp), parameter :: sparsity_limit = 0.5_wp
  !
continue
  !
  ! Default assignments of procedure pointer
  !
  this%dot => dot_intrinsic
  this%mv_mult => mat_vec_mult_intrinsic
  this%vm_mult => vec_mat_mult_intrinsic
  !
  if (allocated(this%mat)) then
    !
    nsz_dim1 = size(this%mat,dim=1)
    nsz_dim2 = size(this%mat,dim=2)
    !
    total_nz = count(this%mat /= zero)
    !
    sparsity = real(total_nz,kind=wp) / real(nsz_dim1*nsz_dim2,kind=wp)
    !
    if (sparsity <= sparsity_limit) then
      !
      ! Convert the mat array into the sparse
      ! arrays within the matrix derived-type
      !
      ! ######################################
      ! ######################################
      ! ######################################
      ! #####                            #####
      ! #####   Compressed Row Storage   #####
      ! #####                            #####
      ! ######################################
      ! ######################################
      ! ######################################
      !
      allocate ( this%row_ptr(1:nsz_dim2+1) , source=0 , stat=ierr )
      allocate ( this%row_idx(1:total_nz) , stat=ierr )
      allocate ( this%row_val(1:total_nz) , stat=ierr )
      !
      !
      !
      do n = 2,nsz_dim2+1
        this%row_ptr(n) = this%row_ptr(n-1) + count(this%mat(:,n-1) /= zero)
      end do
      !
      !
      !
      do n = 1,nsz_dim2
        !
        n1 = this%row_ptr(n)+1
        n2 = this%row_ptr(n+1)
        !
        !
        this%row_idx(n1:n2) = pack( intseq(1,nsz_dim1) , this%mat(:,n) /= zero )
        !
        this%row_val(n1:n2) = this%mat( this%row_idx(n1:n2) , n )
        !
      end do
      !
      ! Assign the procedure pointers to the sparse versions
      !
      this%dot => dot_sparse
      !
      ! #########################################
      ! #########################################
      ! #########################################
      ! #####                               #####
      ! #####   Compressed Column Storage   #####
      ! #####                               #####
      ! #########################################
      ! #########################################
      ! #########################################
      !
      allocate ( this%col_ptr(1:nsz_dim1+1) , source=0 , stat=ierr )
      allocate ( this%col_idx(1:total_nz) , stat=ierr )
      allocate ( this%col_val(1:total_nz) , stat=ierr )
      !
      !
      !
      do n = 2,nsz_dim1+1
        this%col_ptr(n) = this%col_ptr(n-1) + count(this%mat(n-1,:) /= zero)
      end do
      !
      !
      !
      do n = 1,nsz_dim1
        !
        n1 = this%col_ptr(n)+1
        n2 = this%col_ptr(n+1)
        !
        !
        this%col_idx(n1:n2) = pack( intseq(1,nsz_dim2) , this%mat(n,:) /= zero )
        !
        this%col_val(n1:n2) = this%mat( n , this%col_idx(n1:n2) )
        !
      end do
      !
      ! FOR NOW, KEEP THE MATRIX*VECTOR AND VECTOR*MATRIX POINTERS
      ! ASSOCIATED WITH THE INTRINSIC MATMUL FUNCTION. STORAGE FOR
      ! THE SPARSE DATA STILL NEEDS TO FIGURED OUT FOR THE ENTIRE
      ! MATRIX AND NOT JUST FOR THE COLUMN VECTORS THAT ARE NEEDED
      ! FOR THE DOT PRODUCT FUNCTIONS.
      !
     !this%mv_mult => mat_vec_mult_sparse
     !this%vm_mult => vec_mat_mult_sparse
      !
    end if
    !
  end if
  !
end subroutine check_for_sparsity
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function dot_sparse(this,idx,vector) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix),          intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: idx
  real(wp), dimension(:), intent(in) :: vector
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer :: n,n1,n2
  !
#if MATVEC_VERSION == 3
  !.. Local Arrays ..
  integer, dimension(1:num_nz(this,idx)) :: indices
#endif
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = this%row_ptr(idx+1)
  !
#if MATVEC_VERSION == 2
  ! VERSION 2
  !
  return_value = dot( this%row_val(n1:n2) , vector(this%row_idx(n1:n2)) )
  !
#elif MATVEC_VERSION == 3 
  ! VERSION 3
  !
  indices = this%row_idx(n1:n2)
  return_value = dot( this%row_val(n1:n2) , vector(indices) )
  !
#else
  ! VERSION 1
  !
  return_value = zero
  do n = n1,n2
    return_value = return_value + this%row_val(n)*vector(this%row_idx(n))
  end do
  !
#endif
  !
end function dot_sparse
!
!###############################################################################
!
pure function dot_intrinsic(this,idx,vector) result(return_value)
  !
  class(matrix),          intent(in) :: this
  integer,                intent(in) :: idx
  real(wp), dimension(:), intent(in) :: vector
  !
  real(wp) :: return_value
  !
  integer :: np
  !
continue
  !
#if MATVEC_VERSION == 1
  np = size(vector)
  return_value = dot( this%mat(1:np,idx) , vector(1:np) )
#else
  return_value = dot( this%mat(:,idx) , vector )
#endif
  !
end function dot_intrinsic
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function mat_vec_mult_sparse(this,vector) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix),          intent(in) :: this
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: vector
  !
  !.. Function Result ..
  real(wp), dimension(1:size(this%mat,dim=1)) :: return_value
  !
  !.. Local Scalars ..
  integer :: k,n,n1,n2
  real(wp) :: tmp
  !
continue
  !
  do k = 1,size(this%mat,dim=1)
    !
    n1 = this%col_ptr(k)+1
    n2 = this%col_ptr(k+1)
    !
#if MATVEC_VERSION == 1
    ! VERSION 1
    !
    tmp = zero
    do n = n1,n2
      tmp = tmp + this%col_val(n)*vector(this%col_idx(n))
    end do
    return_value(k) = tmp
    !
#else
    ! VERSION 2
    !
    return_value(k) = dot( this%col_val(n1:n2) , vector(this%col_idx(n1:n2)) )
    !
#endif
    !
  end do
  !
end function mat_vec_mult_sparse
!
!###############################################################################
!
pure function mat_vec_mult_intrinsic(this,vector) result(return_value)
  !
  class(matrix),          intent(in) :: this
  real(wp), dimension(:), intent(in) :: vector
  !
  real(wp), dimension(1:size(this%mat,dim=1)) :: return_value
  !
continue
  !
  return_value = matmul( this%mat , vector )
  !
end function mat_vec_mult_intrinsic
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function vec_mat_mult_sparse(this,vector) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix),          intent(in) :: this
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: vector
  !
  !.. Function Result ..
  real(wp), dimension(1:size(this%mat,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: k,n,n1,n2
  real(wp) :: tmp
  !
continue
  !
  do k = 1,size(this%mat,dim=2)
    !
    n1 = this%row_ptr(k)+1
    n2 = this%row_ptr(k+1)
    !
#if MATVEC_VERSION == 1
    ! VERSION 1
    !
    tmp = zero
    do n = n1,n2
      tmp = tmp + this%row_val(n)*vector(this%row_idx(n))
    end do
    return_value(k) = tmp
    !
#else
    ! VERSION 2
    !
    return_value(k) = dot( this%row_val(n1:n2) , vector(this%row_idx(n1:n2)) )
    !
#endif
    !
  end do
  !
end function vec_mat_mult_sparse
!
!###############################################################################
!
pure function vec_mat_mult_intrinsic(this,vector) result(return_value)
  !
  class(matrix),          intent(in) :: this
  real(wp), dimension(:), intent(in) :: vector
  !
  real(wp), dimension(1:size(this%mat,dim=2)) :: return_value
  !
continue
  !
  return_value = matmul( vector , this%mat )
  !
end function vec_mat_mult_intrinsic
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Might need to add mat_mat_mul functions. Probably need to create two versions:
! 1) mat x this%mat
! 2) this%mat x mat
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function num_nz(this,idx) result(return_value)
  !
  class(matrix), intent(in) :: this
  integer, intent(in) :: idx
  !
  integer :: return_value
  !
continue
  !
  return_value = this%row_ptr(idx+1) - this%row_ptr(idx)
  !
end function num_nz
!
!###############################################################################
!
end module generic_types_mod
!
!###############################################################################
!
