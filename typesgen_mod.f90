module generic_types_mod
  !
  use module_kind_types, only : wp,zero,half,dot,intseq
  use module_kind_types, only : alloc_error,io_error
  use module_kind_types, only : error_message
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
   !procedure(dot_vec_function), pointer, private :: dot_vec
   !procedure(dot_mat_function), pointer, private :: dot_mat
    procedure(dot_vec_function), pointer :: dot
    procedure(dot_mat_function), pointer :: dot_mat
    procedure(dot_mat2_function), pointer :: dot_mat2
    procedure(mat_vec_mult), pointer :: mv_mult
    procedure(vec_mat_mult), pointer :: vm_mult
   !procedure(mat_mat_mult),        pointer :: mm_mult
   !procedure(mat_mat_rev_mult),    pointer :: mm_rev_mult
   !procedure(mat_mat_mult_rt),     pointer :: mmrt_mult
   !procedure(mat_mat_rev_mult_rt), pointer :: mmrt_rev_mult
   !procedure(mat_mat_mult_ut),     pointer :: mmut_mult
   !procedure(mat_mat_rev_mult_ut), pointer :: mmut_rev_mult
    !
  contains
    !
    procedure :: check_for_sparsity
    procedure :: output_sparsity
    procedure :: memory_usage
    procedure :: mm_mult       => mat_mat_mult_intrinsic
    procedure :: mm_rev_mult   => mat_mat_rev_mult_intrinsic
    procedure :: mmrt_mult     => mm_mult_intsc_rtrn_transpose
    procedure :: mmrt_rev_mult => mm_rev_mult_intsc_rtrn_transpose
    procedure :: mmut_mult     => mm_mult_intsc_use_transpose
    procedure :: mmut_rev_mult => mm_rev_mult_intsc_use_transpose
   !generic :: dot => dot_vec,dot_mat
    !
  end type matrix
  !
  !
  !
  abstract interface
    !
    pure function dot_vec_function(this,idx,vector) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      integer, intent(in) :: idx
      real(wp), dimension(:), intent(in) :: vector
      real(wp) :: return_value
    end function dot_vec_function
    !
    pure function dot_mat_function(this,idx,array) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      integer, intent(in) :: idx
      real(wp), dimension(:,:), intent(in) :: array
      real(wp) :: return_value(1:size(array,dim=2))
    end function dot_mat_function
    !
    pure function dot_mat2_function(this,idx,array) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      integer, intent(in) :: idx
      real(wp), dimension(:,:,:), intent(in) :: array
      real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
    end function dot_mat2_function
    !
    pure function mat_vec_mult(this,vector) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      real(wp), dimension(:), intent(in) :: vector
      real(wp), dimension(1:size(this%mat,dim=1)) :: return_value
    end function mat_vec_mult
    !
    pure function vec_mat_mult(this,vector) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      real(wp), dimension(:), intent(in) :: vector
      real(wp), dimension(1:size(this%mat,dim=2)) :: return_value
    end function vec_mat_mult
    !
    !
    !
    !
    !
    !
    pure function mat_mat_mult(this,mat) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      real(wp), dimension(:,:), intent(in) :: mat
      real(wp), dimension(1:size(     mat,dim=1), &
                          1:size(this%mat,dim=2)) :: return_value
    end function mat_mat_mult
    !
    pure function mat_mat_rev_mult(this,mat) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      real(wp), dimension(:,:), intent(in) :: mat
      real(wp), dimension(1:size(this%mat,dim=1), &
                          1:size(     mat,dim=2)) :: return_value
    end function mat_mat_rev_mult
    !
    pure function mat_mat_mult_rt(this,mat) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      real(wp), dimension(:,:), intent(in) :: mat
      real(wp), dimension(1:size(this%mat,dim=2), &
                          1:size(     mat,dim=1)) :: return_value
    end function mat_mat_mult_rt
    !
    pure function mat_mat_rev_mult_rt(this,mat) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      real(wp), dimension(:,:), intent(in) :: mat
      real(wp), dimension(1:size(     mat,dim=2), &
                          1:size(this%mat,dim=1)) :: return_value
    end function mat_mat_rev_mult_rt
    !
    pure function mat_mat_mult_ut(this,mat) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      real(wp), dimension(:,:), intent(in) :: mat
      real(wp), dimension(1:size(this%mat,dim=2), &
                          1:size(     mat,dim=2)) :: return_value
    end function mat_mat_mult_ut
    !
    pure function mat_mat_rev_mult_ut(this,mat) result(return_value)
      import :: wp, matrix
      class(matrix), intent(in) :: this
      real(wp), dimension(:,:), intent(in) :: mat
      real(wp), dimension(1:size(     mat,dim=1), &
                          1:size(this%mat,dim=1)) :: return_value
    end function mat_mat_rev_mult_ut
    !
  end interface
  !
contains
!
!###############################################################################
!
elemental subroutine check_for_sparsity(this)
  !
  !.. Use Statements ..
  use ovar, only : use_unrolled_dot_products
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
  this%dot      => dot_intrinsic
  this%dot_mat  => dot_mat_intrinsic
  this%dot_mat2 => dot_mat2_intrinsic
  this%mv_mult  => mat_vec_mult_intrinsic
  this%vm_mult  => vec_mat_mult_intrinsic
  !
 !this%mm_mult       => mat_mat_mult_intrinsic
 !this%mm_rev_mult   => mat_mat_rev_mult_intrinsic
 !this%mmrt_mult     => mm_mult_intsc_rtrn_transpose
 !this%mmrt_rev_mult => mm_rev_mult_intsc_rtrn_transpose
 !this%mmut_mult     => mm_mult_intsc_use_transpose
 !this%mmut_rev_mult => mm_rev_mult_intsc_use_transpose
  !
  if (allocated(this%mat)) then
    !
    nsz_dim1 = size(this%mat,dim=1)
    nsz_dim2 = size(this%mat,dim=2)
    !
    total_nz = count(this%mat /= zero)
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
    ! Compute the sparsity of the matrix as the ratio of non-zero
    ! elements divided by the total number of matrix elements
    !
    sparsity = real(total_nz,kind=wp) / real(nsz_dim1*nsz_dim2,kind=wp)
    !
    ! If the matrix sparsity is less than a given limit, assign
    ! the procedure pointers to the sparse vector operations.
    !
    if (sparsity <= sparsity_limit) then
      !
      ! Assign the procedure pointers to the sparse versions
      !
      if (use_unrolled_dot_products) then
        !
        n1 = minval(this%row_ptr(2:nsz_dim2+1)-this%row_ptr(1:nsz_dim2))
        n2 = maxval(this%row_ptr(2:nsz_dim2+1)-this%row_ptr(1:nsz_dim2))
        !
        if (n1 == n2) then
          !
          select case (n1)
            case (1)
              this%dot => dot_sparse_unrolled_p0
              this%dot_mat => dot_mat_sparse_unrolled_p0
              this%dot_mat2 => dot_mat2_sparse_unrolled_p0
            case (2)
              this%dot => dot_sparse_unrolled_p1
              this%dot_mat => dot_mat_sparse_unrolled_p1
              this%dot_mat2 => dot_mat2_sparse_unrolled_p1
            case (3)
              this%dot => dot_sparse_unrolled_p2
              this%dot_mat => dot_mat_sparse_unrolled_p2
              this%dot_mat2 => dot_mat2_sparse_unrolled_p2
            case (4)
              this%dot => dot_sparse_unrolled_p3
              this%dot_mat => dot_mat_sparse_unrolled_p3
              this%dot_mat2 => dot_mat2_sparse_unrolled_p3
            case (5)
              this%dot => dot_sparse_unrolled_p4
              this%dot_mat => dot_mat_sparse_unrolled_p4
              this%dot_mat2 => dot_mat2_sparse_unrolled_p4
            case (6)
              this%dot => dot_sparse_unrolled_p5
              this%dot_mat => dot_mat_sparse_unrolled_p5
              this%dot_mat2 => dot_mat2_sparse_unrolled_p5
            case (7)
              this%dot => dot_sparse_unrolled_p6
              this%dot_mat => dot_mat_sparse_unrolled_p6
              this%dot_mat2 => dot_mat2_sparse_unrolled_p6
            case (8)
              this%dot => dot_sparse_unrolled_p7
              this%dot_mat => dot_mat_sparse_unrolled_p7
              this%dot_mat2 => dot_mat2_sparse_unrolled_p7
            case (9)
              this%dot => dot_sparse_unrolled_p8
              this%dot_mat => dot_mat_sparse_unrolled_p8
              this%dot_mat2 => dot_mat2_sparse_unrolled_p8
            case (10)
              this%dot => dot_sparse_unrolled_p9
              this%dot_mat => dot_mat_sparse_unrolled_p9
              this%dot_mat2 => dot_mat2_sparse_unrolled_p9
            case default
              this%dot => dot_sparse
              this%dot_mat => dot_mat_sparse
              this%dot_mat2 => dot_mat2_sparse
          end select
          !
        else
          !
          this%dot => dot_sparse
          this%dot_mat => dot_mat_sparse
          this%dot_mat2 => dot_mat2_sparse
          !
        end if
        !
      else
        !
        this%dot => dot_sparse
        this%dot_mat => dot_mat_sparse
        this%dot_mat2 => dot_mat2_sparse
        !
      end if
      !
      ! FOR NOW, KEEP THE MATRIX*VECTOR AND VECTOR*MATRIX POINTERS
      ! ASSOCIATED WITH THE INTRINSIC MATMUL FUNCTION. STORAGE FOR
      ! THE SPARSE DATA STILL NEEDS TO FIGURED OUT FOR THE ENTIRE
      ! MATRIX AND NOT JUST FOR THE COLUMN VECTORS THAT ARE NEEDED
      ! FOR THE DOT PRODUCT FUNCTIONS.
      !
      this%mv_mult => mat_vec_mult_sparse
      this%vm_mult => vec_mat_mult_sparse
     !this%mm_mult       => mat_mat_mult_intrinsic
     !this%mm_rev_mult   => mat_mat_rev_mult_intrinsic
     !this%mmrt_mult     => mm_mult_intsc_rtrn_transpose
     !this%mmrt_rev_mult => mm_rev_mult_intsc_rtrn_transpose
     !this%mmut_mult     => mm_mult_intsc_use_transpose
     !this%mmut_rev_mult => mm_rev_mult_intsc_use_transpose
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
#ifdef USE_INTEL_MKL
  use blas95, only : doti
#endif
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
  integer :: n1,n2
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = this%row_ptr(idx+1)
  !
#ifdef USE_INTEL_MKL
  return_value = doti( this%row_val(n1:n2) , this%row_idx(n1:n2) , vector )
#else
  return_value = dot( this%row_val(n1:n2) , vector(this%row_idx(n1:n2)) )
#endif
  !
end function dot_sparse
!
!###############################################################################
!
pure function dot_sparse_unrolled_p0(this,idx,vector) result(return_value)
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
  integer :: n
  !
continue
  !
  n = this%row_ptr(idx)+1
  !
  return_value = this%row_val(n) * vector( this%row_idx(n) )
  !
end function dot_sparse_unrolled_p0
!
!###############################################################################
!
pure function dot_sparse_unrolled_p1(this,idx,vector) result(return_value)
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
  integer :: n1,n2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2) :: a,b
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  !
  return_value = this%row_val(n1) * vector( this%row_idx(n1) ) + &
                 this%row_val(n2) * vector( this%row_idx(n2) )
  !
end function dot_sparse_unrolled_p1
!
!###############################################################################
!
pure function dot_sparse_unrolled_p2(this,idx,vector) result(return_value)
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
  integer :: n1,n2,n3
  !
  !.. Local Arrays ..
  real(wp), dimension(1:3) :: a,b
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  !
  return_value = this%row_val(n1) * vector( this%row_idx(n1) ) + &
                 this%row_val(n2) * vector( this%row_idx(n2) ) + &
                 this%row_val(n3) * vector( this%row_idx(n3) )
  !
end function dot_sparse_unrolled_p2
!
!###############################################################################
!
pure function dot_sparse_unrolled_p3(this,idx,vector) result(return_value)
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
  integer :: n1,n2,n3,n4
  !
  !.. Local Arrays ..
  real(wp), dimension(1:4) :: a,b
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  n4 = n3+1
  !
  return_value = this%row_val(n1) * vector( this%row_idx(n1) ) + &
                 this%row_val(n2) * vector( this%row_idx(n2) ) + &
                 this%row_val(n3) * vector( this%row_idx(n3) ) + &
                 this%row_val(n4) * vector( this%row_idx(n4) )
  !
end function dot_sparse_unrolled_p3
!
!###############################################################################
!
pure function dot_sparse_unrolled_p4(this,idx,vector) result(return_value)
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
  integer :: n1,n2,n3,n4,n5
  !
  !.. Local Arrays ..
  real(wp), dimension(1:5) :: a,b
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  n4 = n3+1
  n5 = n4+1
  !
  return_value = this%row_val(n1) * vector( this%row_idx(n1) ) + &
                 this%row_val(n2) * vector( this%row_idx(n2) ) + &
                 this%row_val(n3) * vector( this%row_idx(n3) ) + &
                 this%row_val(n4) * vector( this%row_idx(n4) ) + &
                 this%row_val(n5) * vector( this%row_idx(n5) )
  !
end function dot_sparse_unrolled_p4
!
!###############################################################################
!
pure function dot_sparse_unrolled_p5(this,idx,vector) result(return_value)
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
  integer :: n1,n2,n3,n4,n5,n6
  !
  !.. Local Arrays ..
  real(wp), dimension(1:6) :: a,b
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  n4 = n3+1
  n5 = n4+1
  n6 = n5+1
  !
  return_value = this%row_val(n1) * vector( this%row_idx(n1) ) + &
                 this%row_val(n2) * vector( this%row_idx(n2) ) + &
                 this%row_val(n3) * vector( this%row_idx(n3) ) + &
                 this%row_val(n4) * vector( this%row_idx(n4) ) + &
                 this%row_val(n5) * vector( this%row_idx(n5) ) + &
                 this%row_val(n6) * vector( this%row_idx(n6) )
  !
end function dot_sparse_unrolled_p5
!
!###############################################################################
!
pure function dot_sparse_unrolled_p6(this,idx,vector) result(return_value)
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
  integer :: n1,n2,n3,n4,n5,n6,n7
  !
  !.. Local Arrays ..
  real(wp), dimension(1:7) :: a,b
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  n4 = n3+1
  n5 = n4+1
  n6 = n5+1
  n7 = n6+1
  !
  return_value = this%row_val(n1) * vector( this%row_idx(n1) ) + &
                 this%row_val(n2) * vector( this%row_idx(n2) ) + &
                 this%row_val(n3) * vector( this%row_idx(n3) ) + &
                 this%row_val(n4) * vector( this%row_idx(n4) ) + &
                 this%row_val(n5) * vector( this%row_idx(n5) ) + &
                 this%row_val(n6) * vector( this%row_idx(n6) ) + &
                 this%row_val(n7) * vector( this%row_idx(n7) )
  !
end function dot_sparse_unrolled_p6
!
!###############################################################################
!
pure function dot_sparse_unrolled_p7(this,idx,vector) result(return_value)
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
  integer :: n1,n2,n3,n4,n5,n6,n7,n8
  !
  !.. Local Arrays ..
  real(wp), dimension(1:8) :: a,b
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  n4 = n3+1
  n5 = n4+1
  n6 = n5+1
  n7 = n6+1
  n8 = n7+1
  !
  return_value = this%row_val(n1) * vector( this%row_idx(n1) ) + &
                 this%row_val(n2) * vector( this%row_idx(n2) ) + &
                 this%row_val(n3) * vector( this%row_idx(n3) ) + &
                 this%row_val(n4) * vector( this%row_idx(n4) ) + &
                 this%row_val(n5) * vector( this%row_idx(n5) ) + &
                 this%row_val(n6) * vector( this%row_idx(n6) ) + &
                 this%row_val(n7) * vector( this%row_idx(n7) ) + &
                 this%row_val(n8) * vector( this%row_idx(n8) )
  !
end function dot_sparse_unrolled_p7
!
!###############################################################################
!
pure function dot_sparse_unrolled_p8(this,idx,vector) result(return_value)
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
  integer :: n1,n2,n3,n4,n5,n6,n7,n8,n9
  !
  !.. Local Arrays ..
  real(wp), dimension(1:9) :: a,b
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  n4 = n3+1
  n5 = n4+1
  n6 = n5+1
  n7 = n6+1
  n8 = n7+1
  n9 = n8+1
  !
  return_value = this%row_val(n1) * vector( this%row_idx(n1) ) + &
                 this%row_val(n2) * vector( this%row_idx(n2) ) + &
                 this%row_val(n3) * vector( this%row_idx(n3) ) + &
                 this%row_val(n4) * vector( this%row_idx(n4) ) + &
                 this%row_val(n5) * vector( this%row_idx(n5) ) + &
                 this%row_val(n6) * vector( this%row_idx(n6) ) + &
                 this%row_val(n7) * vector( this%row_idx(n7) ) + &
                 this%row_val(n8) * vector( this%row_idx(n8) ) + &
                 this%row_val(n9) * vector( this%row_idx(n9) )
  !
end function dot_sparse_unrolled_p8
!
!###############################################################################
!
pure function dot_sparse_unrolled_p9(this,idx,vector) result(return_value)
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
  integer :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
  !
  !.. Local Arrays ..
  real(wp), dimension(1:10) :: a,b
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  n7  = n6+1
  n8  = n7+1
  n9  = n8+1
  n10 = n9+1
  !
  return_value = this%row_val(n1 ) * vector( this%row_idx(n1 ) ) + &
                 this%row_val(n2 ) * vector( this%row_idx(n2 ) ) + &
                 this%row_val(n3 ) * vector( this%row_idx(n3 ) ) + &
                 this%row_val(n4 ) * vector( this%row_idx(n4 ) ) + &
                 this%row_val(n5 ) * vector( this%row_idx(n5 ) ) + &
                 this%row_val(n6 ) * vector( this%row_idx(n6 ) ) + &
                 this%row_val(n7 ) * vector( this%row_idx(n7 ) ) + &
                 this%row_val(n8 ) * vector( this%row_idx(n8 ) ) + &
                 this%row_val(n9 ) * vector( this%row_idx(n9 ) ) + &
                 this%row_val(n10) * vector( this%row_idx(n10) )
  !
end function dot_sparse_unrolled_p9
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
continue
  !
  return_value = dot( this%mat(:,idx) , vector )
  !
end function dot_intrinsic
!
!###############################################################################
!
pure function dot_mat_intrinsic(this,idx,array) result(return_value)
  !
  class(matrix), intent(in) :: this
  integer,       intent(in) :: idx
  real(wp),      intent(in) :: array(:,:)
  !
  real(wp) :: return_value(1:size(array,dim=2))
  !
continue
  !
  return_value(:) = matmul( this%mat(:,idx) , array )
 !do n = 1,size(array,dim=2)
 !  return_value(n) = dot( this%mat(:,idx) , array(:,n) )
 !end do
  !
end function dot_mat_intrinsic
!
!###############################################################################
!
pure function dot_mat_sparse(this,idx,array) result(return_value)
  !
#ifdef USE_INTEL_MKL
  use blas95, only : doti
#endif
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1,n2
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = this%row_ptr(idx+1)
  !
  do n = 1,size(array,dim=2)
#ifdef USE_INTEL_MKL
    return_value(n) = doti( this%row_val(n1:n2) , &
                            array( this%row_idx(n1:n2) , n ) )
#else
    return_value(n) = dot( this%row_val(n1:n2) , &
                           array( this%row_idx(n1:n2) , n ) )
#endif
  end do
  !
end function dot_mat_sparse
!
!###############################################################################
!
pure function dot_mat_sparse_unrolled_p0(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  !
  do n = 1,size(array,dim=2)
    return_value(n) = this%row_val(n1) * array( this%row_idx(n1) , n )
  end do
  !
end function dot_mat_sparse_unrolled_p0
!
!###############################################################################
!
pure function dot_mat_sparse_unrolled_p1(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1,n2
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  !
  do n = 1,size(array,dim=2)
    return_value(n) = this%row_val(n1) * array( this%row_idx(n1) , n ) + &
                      this%row_val(n2) * array( this%row_idx(n2) , n )
  end do
  !
end function dot_mat_sparse_unrolled_p1
!
!###############################################################################
!
pure function dot_mat_sparse_unrolled_p2(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,n3
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  !
  do n = 1,size(array,dim=2)
    return_value(n) = this%row_val(n1) * array( this%row_idx(n1) , n ) + &
                      this%row_val(n2) * array( this%row_idx(n2) , n ) + &
                      this%row_val(n3) * array( this%row_idx(n3) , n )
  end do
  !
end function dot_mat_sparse_unrolled_p2
!
!###############################################################################
!
pure function dot_mat_sparse_unrolled_p3(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,n3,n4
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  n4 = n3+1
  !
  do n = 1,size(array,dim=2)
    return_value(n) = this%row_val(n1) * array( this%row_idx(n1) , n ) + &
                      this%row_val(n2) * array( this%row_idx(n2) , n ) + &
                      this%row_val(n3) * array( this%row_idx(n3) , n ) + &
                      this%row_val(n4) * array( this%row_idx(n4) , n )
  end do
  !
end function dot_mat_sparse_unrolled_p3
!
!###############################################################################
!
pure function dot_mat_sparse_unrolled_p4(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,n3,n4,n5
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  !
  do n = 1,size(array,dim=2)
    return_value(n) = this%row_val(n1 ) * array( this%row_idx(n1 ) , n ) + &
                      this%row_val(n2 ) * array( this%row_idx(n2 ) , n ) + &
                      this%row_val(n3 ) * array( this%row_idx(n3 ) , n ) + &
                      this%row_val(n4 ) * array( this%row_idx(n4 ) , n ) + &
                      this%row_val(n5 ) * array( this%row_idx(n5 ) , n )
  end do
  !
end function dot_mat_sparse_unrolled_p4
!
!###############################################################################
!
pure function dot_mat_sparse_unrolled_p5(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,n3,n4,n5,n6
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  !
  do n = 1,size(array,dim=2)
    return_value(n) = this%row_val(n1 ) * array( this%row_idx(n1 ) , n ) + &
                      this%row_val(n2 ) * array( this%row_idx(n2 ) , n ) + &
                      this%row_val(n3 ) * array( this%row_idx(n3 ) , n ) + &
                      this%row_val(n4 ) * array( this%row_idx(n4 ) , n ) + &
                      this%row_val(n5 ) * array( this%row_idx(n5 ) , n ) + &
                      this%row_val(n6 ) * array( this%row_idx(n6 ) , n )
  end do
  !
end function dot_mat_sparse_unrolled_p5
!
!###############################################################################
!
pure function dot_mat_sparse_unrolled_p6(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,n3,n4,n5,n6,n7
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  n7  = n6+1
  !
  do n = 1,size(array,dim=2)
    return_value(n) = this%row_val(n1 ) * array( this%row_idx(n1 ) , n ) + &
                      this%row_val(n2 ) * array( this%row_idx(n2 ) , n ) + &
                      this%row_val(n3 ) * array( this%row_idx(n3 ) , n ) + &
                      this%row_val(n4 ) * array( this%row_idx(n4 ) , n ) + &
                      this%row_val(n5 ) * array( this%row_idx(n5 ) , n ) + &
                      this%row_val(n6 ) * array( this%row_idx(n6 ) , n ) + &
                      this%row_val(n7 ) * array( this%row_idx(n7 ) , n )
  end do
  !
end function dot_mat_sparse_unrolled_p6
!
!###############################################################################
!
pure function dot_mat_sparse_unrolled_p7(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,n3,n4,n5,n6,n7,n8
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  n7  = n6+1
  n8  = n7+1
  !
  do n = 1,size(array,dim=2)
    return_value(n) = this%row_val(n1 ) * array( this%row_idx(n1 ) , n ) + &
                      this%row_val(n2 ) * array( this%row_idx(n2 ) , n ) + &
                      this%row_val(n3 ) * array( this%row_idx(n3 ) , n ) + &
                      this%row_val(n4 ) * array( this%row_idx(n4 ) , n ) + &
                      this%row_val(n5 ) * array( this%row_idx(n5 ) , n ) + &
                      this%row_val(n6 ) * array( this%row_idx(n6 ) , n ) + &
                      this%row_val(n7 ) * array( this%row_idx(n7 ) , n ) + &
                      this%row_val(n8 ) * array( this%row_idx(n8 ) , n )
  end do
  !
end function dot_mat_sparse_unrolled_p7
!
!###############################################################################
!
pure function dot_mat_sparse_unrolled_p8(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,n3,n4,n5,n6,n7,n8,n9
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  n7  = n6+1
  n8  = n7+1
  n9  = n8+1
  !
  do n = 1,size(array,dim=2)
    return_value(n) = this%row_val(n1 ) * array( this%row_idx(n1 ) , n ) + &
                      this%row_val(n2 ) * array( this%row_idx(n2 ) , n ) + &
                      this%row_val(n3 ) * array( this%row_idx(n3 ) , n ) + &
                      this%row_val(n4 ) * array( this%row_idx(n4 ) , n ) + &
                      this%row_val(n5 ) * array( this%row_idx(n5 ) , n ) + &
                      this%row_val(n6 ) * array( this%row_idx(n6 ) , n ) + &
                      this%row_val(n7 ) * array( this%row_idx(n7 ) , n ) + &
                      this%row_val(n8 ) * array( this%row_idx(n8 ) , n ) + &
                      this%row_val(n9 ) * array( this%row_idx(n9 ) , n )
  end do
  !
end function dot_mat_sparse_unrolled_p8
!
!###############################################################################
!
pure function dot_mat_sparse_unrolled_p9(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2))
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  n7  = n6+1
  n8  = n7+1
  n9  = n8+1
  n10 = n9+1
  !
  do n = 1,size(array,dim=2)
    return_value(n) = this%row_val(n1 ) * array( this%row_idx(n1 ) , n ) + &
                      this%row_val(n2 ) * array( this%row_idx(n2 ) , n ) + &
                      this%row_val(n3 ) * array( this%row_idx(n3 ) , n ) + &
                      this%row_val(n4 ) * array( this%row_idx(n4 ) , n ) + &
                      this%row_val(n5 ) * array( this%row_idx(n5 ) , n ) + &
                      this%row_val(n6 ) * array( this%row_idx(n6 ) , n ) + &
                      this%row_val(n7 ) * array( this%row_idx(n7 ) , n ) + &
                      this%row_val(n8 ) * array( this%row_idx(n8 ) , n ) + &
                      this%row_val(n9 ) * array( this%row_idx(n9 ) , n ) + &
                      this%row_val(n10) * array( this%row_idx(n10) , n )
  end do
  !
end function dot_mat_sparse_unrolled_p9
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function dot_mat2_intrinsic(this,idx,array) result(return_value)
  !
  class(matrix), intent(in) :: this
  integer,       intent(in) :: idx
  real(wp),      intent(in) :: array(:,:,:)
  !
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  integer :: n
  !
continue
  !
  do n = 1,size(array,dim=3)
    return_value(:,n) = matmul( this%mat(:,idx) , array(:,:,n) )
  end do
  !
end function dot_mat2_intrinsic
!
!###############################################################################
!
pure function dot_mat2_sparse(this,idx,array) result(return_value)
  !
#ifdef USE_INTEL_MKL
  use blas95, only : doti
#endif
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1,n2
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = this%row_ptr(idx+1)
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
#ifdef USE_INTEL_MKL
      return_value(n,m) = doti( this%row_val(n1:n2) , &
                              array( this%row_idx(n1:n2) , n,m ) )
#else
      return_value(n,m) = dot( this%row_val(n1:n2) , &
                             array( this%row_idx(n1:n2) , n,m ) )
#endif
    end do
  end do
  !
end function dot_mat2_sparse
!
!###############################################################################
!
pure function dot_mat2_sparse_unrolled_p0(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
      return_value(n,m) = this%row_val(n1) * array( this%row_idx(n1) , n,m )
    end do
  end do
  !
end function dot_mat2_sparse_unrolled_p0
!
!###############################################################################
!
pure function dot_mat2_sparse_unrolled_p1(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1,n2
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
      return_value(n,m) = this%row_val(n1) * array( this%row_idx(n1) , n,m ) + &
                          this%row_val(n2) * array( this%row_idx(n2) , n,m )
    end do
  end do
  !
end function dot_mat2_sparse_unrolled_p1
!
!###############################################################################
!
pure function dot_mat2_sparse_unrolled_p2(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1,n2,n3
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
      return_value(n,m) = this%row_val(n1) * array( this%row_idx(n1) , n,m ) + &
                          this%row_val(n2) * array( this%row_idx(n2) , n,m ) + &
                          this%row_val(n3) * array( this%row_idx(n3) , n,m )
    end do
  end do
  !
end function dot_mat2_sparse_unrolled_p2
!
!###############################################################################
!
pure function dot_mat2_sparse_unrolled_p3(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1,n2,n3,n4
  !
continue
  !
  n1 = this%row_ptr(idx)+1
  n2 = n1+1
  n3 = n2+1
  n4 = n3+1
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
      return_value(n,m) = this%row_val(n1) * array( this%row_idx(n1) , n,m ) + &
                          this%row_val(n2) * array( this%row_idx(n2) , n,m ) + &
                          this%row_val(n3) * array( this%row_idx(n3) , n,m ) + &
                          this%row_val(n4) * array( this%row_idx(n4) , n,m )
    end do
  end do
  !
end function dot_mat2_sparse_unrolled_p3
!
!###############################################################################
!
pure function dot_mat2_sparse_unrolled_p4(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1,n2,n3,n4,n5
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
      return_value(n,m) = this%row_val(n1) * array( this%row_idx(n1) , n,m ) + &
                          this%row_val(n2) * array( this%row_idx(n2) , n,m ) + &
                          this%row_val(n3) * array( this%row_idx(n3) , n,m ) + &
                          this%row_val(n4) * array( this%row_idx(n4) , n,m ) + &
                          this%row_val(n5) * array( this%row_idx(n5) , n,m )
    end do
  end do
  !
end function dot_mat2_sparse_unrolled_p4
!
!###############################################################################
!
pure function dot_mat2_sparse_unrolled_p5(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1,n2,n3,n4,n5,n6
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
      return_value(n,m) = this%row_val(n1) * array( this%row_idx(n1) , n,m ) + &
                          this%row_val(n2) * array( this%row_idx(n2) , n,m ) + &
                          this%row_val(n3) * array( this%row_idx(n3) , n,m ) + &
                          this%row_val(n4) * array( this%row_idx(n4) , n,m ) + &
                          this%row_val(n5) * array( this%row_idx(n5) , n,m ) + &
                          this%row_val(n6) * array( this%row_idx(n6) , n,m )
    end do
  end do
  !
end function dot_mat2_sparse_unrolled_p5
!
!###############################################################################
!
pure function dot_mat2_sparse_unrolled_p6(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1,n2,n3,n4,n5,n6,n7
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  n7  = n6+1
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
      return_value(n,m) = this%row_val(n1) * array( this%row_idx(n1) , n,m ) + &
                          this%row_val(n2) * array( this%row_idx(n2) , n,m ) + &
                          this%row_val(n3) * array( this%row_idx(n3) , n,m ) + &
                          this%row_val(n4) * array( this%row_idx(n4) , n,m ) + &
                          this%row_val(n5) * array( this%row_idx(n5) , n,m ) + &
                          this%row_val(n6) * array( this%row_idx(n6) , n,m ) + &
                          this%row_val(n7) * array( this%row_idx(n7) , n,m )
    end do
  end do
  !
end function dot_mat2_sparse_unrolled_p6
!
!###############################################################################
!
pure function dot_mat2_sparse_unrolled_p7(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1,n2,n3,n4,n5,n6,n7,n8
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  n7  = n6+1
  n8  = n7+1
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
      return_value(n,m) = this%row_val(n1) * array( this%row_idx(n1) , n,m ) + &
                          this%row_val(n2) * array( this%row_idx(n2) , n,m ) + &
                          this%row_val(n3) * array( this%row_idx(n3) , n,m ) + &
                          this%row_val(n4) * array( this%row_idx(n4) , n,m ) + &
                          this%row_val(n5) * array( this%row_idx(n5) , n,m ) + &
                          this%row_val(n6) * array( this%row_idx(n6) , n,m ) + &
                          this%row_val(n7) * array( this%row_idx(n7) , n,m ) + &
                          this%row_val(n8) * array( this%row_idx(n8) , n,m )
    end do
  end do
  !
end function dot_mat2_sparse_unrolled_p7
!
!###############################################################################
!
pure function dot_mat2_sparse_unrolled_p8(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1,n2,n3,n4,n5,n6,n7,n8,n9
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  n7  = n6+1
  n8  = n7+1
  n9  = n8+1
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
      return_value(n,m) = this%row_val(n1) * array( this%row_idx(n1) , n,m ) + &
                          this%row_val(n2) * array( this%row_idx(n2) , n,m ) + &
                          this%row_val(n3) * array( this%row_idx(n3) , n,m ) + &
                          this%row_val(n4) * array( this%row_idx(n4) , n,m ) + &
                          this%row_val(n5) * array( this%row_idx(n5) , n,m ) + &
                          this%row_val(n6) * array( this%row_idx(n6) , n,m ) + &
                          this%row_val(n7) * array( this%row_idx(n7) , n,m ) + &
                          this%row_val(n8) * array( this%row_idx(n8) , n,m ) + &
                          this%row_val(n9) * array( this%row_idx(n9) , n,m )
    end do
  end do
  !
end function dot_mat2_sparse_unrolled_p8
!
!###############################################################################
!
pure function dot_mat2_sparse_unrolled_p9(this,idx,array) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: idx
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(array,dim=2),1:size(array,dim=3))
  !
  !.. Local Scalars ..
  integer :: n,m,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
  !
continue
  !
  n1  = this%row_ptr(idx)+1
  n2  = n1+1
  n3  = n2+1
  n4  = n3+1
  n5  = n4+1
  n6  = n5+1
  n7  = n6+1
  n8  = n7+1
  n9  = n8+1
  n10 = n9+1
  !
  do m = 1,size(array,dim=3)
    do n = 1,size(array,dim=2)
      return_value(n,m) = &
                        this%row_val(n1 ) * array( this%row_idx(n1 ) , n,m ) + &
                        this%row_val(n2 ) * array( this%row_idx(n2 ) , n,m ) + &
                        this%row_val(n3 ) * array( this%row_idx(n3 ) , n,m ) + &
                        this%row_val(n4 ) * array( this%row_idx(n4 ) , n,m ) + &
                        this%row_val(n5 ) * array( this%row_idx(n5 ) , n,m ) + &
                        this%row_val(n6 ) * array( this%row_idx(n6 ) , n,m ) + &
                        this%row_val(n7 ) * array( this%row_idx(n7 ) , n,m ) + &
                        this%row_val(n8 ) * array( this%row_idx(n8 ) , n,m ) + &
                        this%row_val(n9 ) * array( this%row_idx(n9 ) , n,m ) + &
                        this%row_val(n10) * array( this%row_idx(n10) , n,m )
    end do
  end do
  !
end function dot_mat2_sparse_unrolled_p9
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
  integer :: k,n1,n2
  !
continue
  !
  do k = 1,size(this%mat,dim=1)
    !
    n1 = this%col_ptr(k)+1
    n2 = this%col_ptr(k+1)
    !
    return_value(k) = dot( this%col_val(n1:n2) , vector(this%col_idx(n1:n2)) )
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
  integer :: k,n1,n2
  !
continue
  !
  do k = 1,size(this%mat,dim=2)
    !
    n1 = this%row_ptr(k)+1
    n2 = this%row_ptr(k+1)
    !
    return_value(k) = dot( this%row_val(n1:n2) , vector(this%row_idx(n1:n2)) )
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
! Might need to add additional mat_mat_mul functions.
! Probably need to create two versions:
! 1) mat x this%mat
! 2) this%mat x mat
!
!###############################################################################
!
pure function mat_mat_mult_intrinsic(this,mat) result(return_value)
  !
  class(matrix),            intent(in) :: this
  real(wp), dimension(:,:), intent(in) :: mat
  !
  real(wp), dimension(1:size(     mat,dim=1), &
                      1:size(this%mat,dim=2)) :: return_value
  !
continue
  !
  return_value = matmul( mat , this%mat )
  !
end function mat_mat_mult_intrinsic
!
!###############################################################################
!
pure function mat_mat_rev_mult_intrinsic(this,mat) result(return_value)
  !
  class(matrix),            intent(in) :: this
  real(wp), dimension(:,:), intent(in) :: mat
  !
  real(wp), dimension(1:size(this%mat,dim=1), &
                      1:size(     mat,dim=2)) :: return_value
  !
continue
  !
  return_value = matmul( this%mat , mat )
  !
end function mat_mat_rev_mult_intrinsic
!
!###############################################################################
!
pure function mm_mult_intsc_rtrn_transpose(this,mat) result(return_value)
  !
  class(matrix),            intent(in) :: this
  real(wp), dimension(:,:), intent(in) :: mat
  !
  real(wp), dimension(1:size(this%mat,dim=2), &
                      1:size(     mat,dim=1)) :: return_value
  !
continue
  !
  return_value = transpose( matmul( mat , this%mat ) )
  !
end function mm_mult_intsc_rtrn_transpose
!
!###############################################################################
!
pure function mm_rev_mult_intsc_rtrn_transpose(this,mat) result(return_value)
  !
  class(matrix),            intent(in) :: this
  real(wp), dimension(:,:), intent(in) :: mat
  !
  real(wp), dimension(1:size(     mat,dim=2), &
                      1:size(this%mat,dim=1)) :: return_value
  !
continue
  !
  return_value = transpose( matmul( this%mat , mat ) )
  !
end function mm_rev_mult_intsc_rtrn_transpose
!
!###############################################################################
!
pure function mm_mult_intsc_use_transpose(this,mat) result(return_value)
  !
  class(matrix),            intent(in) :: this
  real(wp), dimension(:,:), intent(in) :: mat
  !
  real(wp), dimension(1:size(this%mat,dim=2), &
                      1:size(     mat,dim=2)) :: return_value
  !
continue
  !
  return_value = transpose( matmul( transpose(mat) , this%mat ) )
  !
end function mm_mult_intsc_use_transpose
!
!###############################################################################
!
pure function mm_rev_mult_intsc_use_transpose(this,mat) result(return_value)
  !
  class(matrix),            intent(in) :: this
  real(wp), dimension(:,:), intent(in) :: mat
  !
  real(wp), dimension(1:size(     mat,dim=1), &
                      1:size(this%mat,dim=1)) :: return_value
  !
continue
  !
  return_value = transpose( matmul( this%mat , transpose(mat) ) )
  !
end function mm_rev_mult_intsc_use_transpose
!
!###############################################################################
!###############################################################################
!###############################################################################
!
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
subroutine output_sparsity(this,fname)
  !
  !.. Passed Argument from Invoking Object ..
  class(matrix), intent(inout) :: this
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: fname
  !
  !.. Local Scalars ..
  integer :: ifu,ierr
  integer :: i,j,n1,n2,np1,np2
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xy(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_sparsity"
  !
continue
  !
  n1 = size(this%mat,dim=1)
  n2 = size(this%mat,dim=2)
  !
  np1 = n1 + 1
  np2 = n2 + 1
  !
  allocate ( xy(1:np1,1:np2,1:2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xy",1,__LINE__,__FILE__,ierr,error_message)
  !
  do j = 1,np2
    do i = 1,np1
      xy(i,j,1) = real(i-1,kind=wp) + half
      xy(i,j,2) = real(j-1,kind=wp) + half
    end do
  end do
  !
  open (newunit=ifu,file=trim(adjustl(fname)),status="replace",action="write", &
        position="rewind",form="formatted",iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (ifu,1) np1,np2
  !
  write (ifu,2) ((xy(i,j,1), i=1,np1), j=1,np2)
  write (ifu,2) ((xy(i,j,2), i=1,np1), j=1,np2)
  !
 !write (ifu,3) ((merge( 1 , 0 , abs(this%mat(i,j)) > zero ), i=1,n1), &
 !                                                            j=n2,1,-1)
  write (ifu,3) ((merge( 1 , 0 , abs(this%mat(i,j)) > zero ), i=1,n1), &
                                                              j=1,n2)
  !
  close (ifu,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( xy , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xy",2,__LINE__,__FILE__,ierr,error_message)
  !
  1 format ('VARIABLES = "Vector" "Node" "Value"',/, &
            'Zone T="Matrix Sparsity", I=',i0,', J=',i0, &
            ', DATAPACKING=BLOCK, VARLOCATION=([3]=CELLCENTERED)')
  2 format (12f7.1)
  3 format (42i2)
  !
end subroutine output_sparsity
!
!###############################################################################
!
pure function memory_usage(this) result(return_value)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  class(matrix), intent(in) :: this
  !
  !.. Function Result ..
  type(memory) :: return_value
  !
  !.. Local Scalars ..
  type(memory) :: array
  !
#ifndef DISABLE_DTIO
continue
  !
  call return_value%reset
  !
  if (allocated(this%mat)) then
    call array%set( storage_size(this%mat,kind=inttype) , &
                    size(this%mat,kind=inttype) )
    call return_value%add(array)
  end if
  !
  if (allocated(this%row_ptr)) then
    call array%set( storage_size(this%row_ptr,kind=inttype) , &
                    size(this%row_ptr,kind=inttype) )
    call return_value%add(array)
  end if
  !
  if (allocated(this%row_idx)) then
    call array%set( storage_size(this%row_idx,kind=inttype) , &
                    size(this%row_idx,kind=inttype) )
    call return_value%add(array)
  end if
  !
  if (allocated(this%row_val)) then
    call array%set( storage_size(this%row_val,kind=inttype) , &
                    size(this%row_val,kind=inttype) )
    call return_value%add(array)
  end if
  !
  if (allocated(this%col_ptr)) then
    call array%set( storage_size(this%col_ptr,kind=inttype) , &
                    size(this%col_ptr,kind=inttype) )
    call return_value%add(array)
  end if
  !
  if (allocated(this%col_idx)) then
    call array%set( storage_size(this%col_idx,kind=inttype) , &
                    size(this%col_idx,kind=inttype) )
    call return_value%add(array)
  end if
  !
  if (allocated(this%col_val)) then
    call array%set( storage_size(this%col_val,kind=inttype) , &
                    size(this%col_val,kind=inttype) )
    call return_value%add(array)
  end if
  !
#endif
end function memory_usage
!
!###############################################################################
!
end module generic_types_mod
!
!###############################################################################
!
