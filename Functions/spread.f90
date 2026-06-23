!
!###############################################################################
!###############################################################################
!###############################################################################
! Spread functions with input argument a as 32-bit integer
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function spread_i4_rank0(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => int32
  !
  integer(lp), intent(in) :: a
  integer,     intent(in) :: dim,ncopies
  !
  integer(lp), allocatable :: return_value(:)
  !
  integer :: i
  !
continue
  !
  allocate ( return_value(1:ncopies) , source=a )
  !
end function spread_i4_rank0
!
pure function spread_i4_rank1(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => int32
  !
  integer(lp), intent(in) :: a(:)
  integer,     intent(in) :: dim,ncopies
  !
  integer(lp), allocatable :: return_value(:,:)
  !
  integer :: i,j
  integer :: n1
  !
continue
  !
  n1 = size(a)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1) )
    do i = 1,ncopies
      return_value(i,1:n1) = a(1:n1)
    end do
  else
    allocate ( return_value(1:n1,1:ncopies) )
    do j = 1,ncopies
      return_value(1:n1,j) = a(1:n1)
    end do
  end if
  !
end function spread_i4_rank1
!
pure function spread_i4_rank2(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => int32
  !
  integer(lp), intent(in) :: a(:,:)
  integer,     intent(in) :: dim,ncopies
  !
  integer(lp), allocatable :: return_value(:,:,:)
  !
  integer :: i,j,k
  integer :: n1,n2
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2) = a(1:n1,1:n2)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2) = a(1:n1,1:n2)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:ncopies) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k) = a(1:n1,1:n2)
    end do
  end if
  !
end function spread_i4_rank2
!
pure function spread_i4_rank3(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => int32
  !
  integer(lp), intent(in) :: a(:,:,:)
  integer,     intent(in) :: dim,ncopies
  !
  integer(lp), allocatable :: return_value(:,:,:,:)
  !
  integer :: i,j,k,l
  integer :: n1,n2,n3
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  n3 = size(a,dim=3)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2,1:n3) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2,1:n3) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 3) then
    allocate ( return_value(1:n1,1:n2,1:ncopies,1:n3) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:n3,1:ncopies) )
    do l = 1,ncopies
      return_value(1:n1,1:n2,1:n3,l) = a(1:n1,1:n2,1:n3)
    end do
  end if
  !
end function spread_i4_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! Spread functions with input argument a as 64-bit integer
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function spread_i8_rank0(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => int64
  !
  integer(lp), intent(in) :: a
  integer,     intent(in) :: dim,ncopies
  !
  integer(lp), allocatable :: return_value(:)
  !
  integer :: i
  !
continue
  !
  allocate ( return_value(1:ncopies) , source=a )
  !
end function spread_i8_rank0
!
pure function spread_i8_rank1(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => int64
  !
  integer(lp), intent(in) :: a(:)
  integer,     intent(in) :: dim,ncopies
  !
  integer(lp), allocatable :: return_value(:,:)
  !
  integer :: i,j
  integer :: n1
  !
continue
  !
  n1 = size(a)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1) )
    do i = 1,ncopies
      return_value(i,1:n1) = a(1:n1)
    end do
  else
    allocate ( return_value(1:n1,1:ncopies) )
    do j = 1,ncopies
      return_value(1:n1,j) = a(1:n1)
    end do
  end if
  !
end function spread_i8_rank1
!
pure function spread_i8_rank2(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => int64
  !
  integer(lp), intent(in) :: a(:,:)
  integer,     intent(in) :: dim,ncopies
  !
  integer(lp), allocatable :: return_value(:,:,:)
  !
  integer :: i,j,k
  integer :: n1,n2
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2) = a(1:n1,1:n2)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2) = a(1:n1,1:n2)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:ncopies) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k) = a(1:n1,1:n2)
    end do
  end if
  !
end function spread_i8_rank2
!
pure function spread_i8_rank3(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => int64
  !
  integer(lp), intent(in) :: a(:,:,:)
  integer,     intent(in) :: dim,ncopies
  !
  integer(lp), allocatable :: return_value(:,:,:,:)
  !
  integer :: i,j,k,l
  integer :: n1,n2,n3
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  n3 = size(a,dim=3)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2,1:n3) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2,1:n3) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 3) then
    allocate ( return_value(1:n1,1:n2,1:ncopies,1:n3) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:n3,1:ncopies) )
    do l = 1,ncopies
      return_value(1:n1,1:n2,1:n3,l) = a(1:n1,1:n2,1:n3)
    end do
  end if
  !
end function spread_i8_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! Spread functions with input argument a as 32-bit real
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function spread_r4_rank0(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real32
  !
  real(lp), intent(in) :: a
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:)
  !
  integer :: i
  !
continue
  !
  allocate ( return_value(1:ncopies) , source=a )
  !
end function spread_r4_rank0
!
pure function spread_r4_rank1(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real32
  !
  real(lp), intent(in) :: a(:)
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:,:)
  !
  integer :: i,j
  integer :: n1
  !
continue
  !
  n1 = size(a)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1) )
    do i = 1,ncopies
      return_value(i,1:n1) = a(1:n1)
    end do
  else
    allocate ( return_value(1:n1,1:ncopies) )
    do j = 1,ncopies
      return_value(1:n1,j) = a(1:n1)
    end do
  end if
  !
end function spread_r4_rank1
!
pure function spread_r4_rank2(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real32
  !
  real(lp), intent(in) :: a(:,:)
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:,:,:)
  !
  integer :: i,j,k
  integer :: n1,n2
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2) = a(1:n1,1:n2)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2) = a(1:n1,1:n2)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:ncopies) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k) = a(1:n1,1:n2)
    end do
  end if
  !
end function spread_r4_rank2
!
pure function spread_r4_rank3(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real32
  !
  real(lp), intent(in) :: a(:,:,:)
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:,:,:,:)
  !
  integer :: i,j,k,l
  integer :: n1,n2,n3
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  n3 = size(a,dim=3)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2,1:n3) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2,1:n3) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 3) then
    allocate ( return_value(1:n1,1:n2,1:ncopies,1:n3) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:n3,1:ncopies) )
    do l = 1,ncopies
      return_value(1:n1,1:n2,1:n3,l) = a(1:n1,1:n2,1:n3)
    end do
  end if
  !
end function spread_r4_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! Spread functions with input argument a as 64-bit real
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function spread_r8_rank0(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real64
  !
  real(lp), intent(in) :: a
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:)
  !
  integer :: i
  !
continue
  !
  allocate ( return_value(1:ncopies) , source=a )
  !
end function spread_r8_rank0
!
pure function spread_r8_rank1(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real64
  !
  real(lp), intent(in) :: a(:)
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:,:)
  !
  integer :: i,j
  integer :: n1
  !
continue
  !
  n1 = size(a)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1) )
    do i = 1,ncopies
      return_value(i,1:n1) = a(1:n1)
    end do
  else
    allocate ( return_value(1:n1,1:ncopies) )
    do j = 1,ncopies
      return_value(1:n1,j) = a(1:n1)
    end do
  end if
  !
end function spread_r8_rank1
!
pure function spread_r8_rank2(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real64
  !
  real(lp), intent(in) :: a(:,:)
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:,:,:)
  !
  integer :: i,j,k
  integer :: n1,n2
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2) = a(1:n1,1:n2)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2) = a(1:n1,1:n2)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:ncopies) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k) = a(1:n1,1:n2)
    end do
  end if
  !
end function spread_r8_rank2
!
pure function spread_r8_rank3(a,dim,ncopies) result(return_value)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real64
  !
  real(lp), intent(in) :: a(:,:,:)
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:,:,:,:)
  !
  integer :: i,j,k,l
  integer :: n1,n2,n3
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  n3 = size(a,dim=3)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2,1:n3) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2,1:n3) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 3) then
    allocate ( return_value(1:n1,1:n2,1:ncopies,1:n3) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:n3,1:ncopies) )
    do l = 1,ncopies
      return_value(1:n1,1:n2,1:n3,l) = a(1:n1,1:n2,1:n3)
    end do
  end if
  !
end function spread_r8_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! Spread functions with input argument a as 128-bit real
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function spread_r16_rank0(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = QP
  !
  real(lp), intent(in) :: a
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:)
  !
  integer :: i
  !
continue
  !
  allocate ( return_value(1:ncopies) , source=a )
  !
end function spread_r16_rank0
!
pure function spread_r16_rank1(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = QP
  !
  real(lp), intent(in) :: a(:)
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:,:)
  !
  integer :: i,j
  integer :: n1
  !
continue
  !
  n1 = size(a)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1) )
    do i = 1,ncopies
      return_value(i,1:n1) = a(1:n1)
    end do
  else
    allocate ( return_value(1:n1,1:ncopies) )
    do j = 1,ncopies
      return_value(1:n1,j) = a(1:n1)
    end do
  end if
  !
end function spread_r16_rank1
!
pure function spread_r16_rank2(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = QP
  !
  real(lp), intent(in) :: a(:,:)
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:,:,:)
  !
  integer :: i,j,k
  integer :: n1,n2
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2) = a(1:n1,1:n2)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2) = a(1:n1,1:n2)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:ncopies) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k) = a(1:n1,1:n2)
    end do
  end if
  !
end function spread_r16_rank2
!
pure function spread_r16_rank3(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = QP
  !
  real(lp), intent(in) :: a(:,:,:)
  integer,  intent(in) :: dim,ncopies
  !
  real(lp), allocatable :: return_value(:,:,:,:)
  !
  integer :: i,j,k,l
  integer :: n1,n2,n3
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  n3 = size(a,dim=3)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2,1:n3) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2,1:n3) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 3) then
    allocate ( return_value(1:n1,1:n2,1:ncopies,1:n3) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:n3,1:ncopies) )
    do l = 1,ncopies
      return_value(1:n1,1:n2,1:n3,l) = a(1:n1,1:n2,1:n3)
    end do
  end if
  !
end function spread_r16_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! Spread functions with input argument a as default logical kind
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function spread_l4_rank0(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = ldk
  !
  logical(lp), intent(in) :: a
  integer,     intent(in) :: dim,ncopies
  !
  logical(lp), allocatable :: return_value(:)
  !
  integer :: i
  !
continue
  !
  allocate ( return_value(1:ncopies) , source=a )
  !
end function spread_l4_rank0
!
pure function spread_l4_rank1(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = ldk
  !
  logical(lp), intent(in) :: a(:)
  integer,     intent(in) :: dim,ncopies
  !
  logical(lp), allocatable :: return_value(:,:)
  !
  integer :: i,j
  integer :: n1
  !
continue
  !
  n1 = size(a)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1) )
    do i = 1,ncopies
      return_value(i,1:n1) = a(1:n1)
    end do
  else
    allocate ( return_value(1:n1,1:ncopies) )
    do j = 1,ncopies
      return_value(1:n1,j) = a(1:n1)
    end do
  end if
  !
end function spread_l4_rank1
!
pure function spread_l4_rank2(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = ldk
  !
  logical(lp), intent(in) :: a(:,:)
  integer,     intent(in) :: dim,ncopies
  !
  logical(lp), allocatable :: return_value(:,:,:)
  !
  integer :: i,j,k
  integer :: n1,n2
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2) = a(1:n1,1:n2)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2) = a(1:n1,1:n2)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:ncopies) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k) = a(1:n1,1:n2)
    end do
  end if
  !
end function spread_l4_rank2
!
pure function spread_l4_rank3(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = ldk
  !
  logical(lp), intent(in) :: a(:,:,:)
  integer,     intent(in) :: dim,ncopies
  !
  logical(lp), allocatable :: return_value(:,:,:,:)
  !
  integer :: i,j,k,l
  integer :: n1,n2,n3
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  n3 = size(a,dim=3)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2,1:n3) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2,1:n3) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 3) then
    allocate ( return_value(1:n1,1:n2,1:ncopies,1:n3) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:n3,1:ncopies) )
    do l = 1,ncopies
      return_value(1:n1,1:n2,1:n3,l) = a(1:n1,1:n2,1:n3)
    end do
  end if
  !
end function spread_l4_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! Spread functions with input argument a as logical kind of least precision
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function spread_lk_rank0(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = lk
  !
  logical(lp), intent(in) :: a
  integer,     intent(in) :: dim,ncopies
  !
  logical(lp), allocatable :: return_value(:)
  !
  integer :: i
  !
continue
  !
  allocate ( return_value(1:ncopies) , source=a )
  !
end function spread_lk_rank0
!
pure function spread_lk_rank1(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = lk
  !
  logical(lp), intent(in) :: a(:)
  integer,     intent(in) :: dim,ncopies
  !
  logical(lp), allocatable :: return_value(:,:)
  !
  integer :: i,j
  integer :: n1
  !
continue
  !
  n1 = size(a)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1) )
    do i = 1,ncopies
      return_value(i,1:n1) = a(1:n1)
    end do
  else
    allocate ( return_value(1:n1,1:ncopies) )
    do j = 1,ncopies
      return_value(1:n1,j) = a(1:n1)
    end do
  end if
  !
end function spread_lk_rank1
!
pure function spread_lk_rank2(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = lk
  !
  logical(lp), intent(in) :: a(:,:)
  integer,     intent(in) :: dim,ncopies
  !
  logical(lp), allocatable :: return_value(:,:,:)
  !
  integer :: i,j,k
  integer :: n1,n2
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2) = a(1:n1,1:n2)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2) = a(1:n1,1:n2)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:ncopies) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k) = a(1:n1,1:n2)
    end do
  end if
  !
end function spread_lk_rank2
!
pure function spread_lk_rank3(a,dim,ncopies) result(return_value)
  !
  integer, parameter :: lp = lk
  !
  logical(lp), intent(in) :: a(:,:,:)
  integer,     intent(in) :: dim,ncopies
  !
  logical(lp), allocatable :: return_value(:,:,:,:)
  !
  integer :: i,j,k,l
  integer :: n1,n2,n3
  !
continue
  !
  n1 = size(a,dim=1)
  n2 = size(a,dim=2)
  n3 = size(a,dim=3)
  !
  if (dim == 1) then
    allocate ( return_value(1:ncopies,1:n1,1:n2,1:n3) )
    do i = 1,ncopies
      return_value(i,1:n1,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 2) then
    allocate ( return_value(1:n1,1:ncopies,1:n2,1:n3) )
    do j = 1,ncopies
      return_value(1:n1,j,1:n2,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else if (dim == 3) then
    allocate ( return_value(1:n1,1:n2,1:ncopies,1:n3) )
    do k = 1,ncopies
      return_value(1:n1,1:n2,k,1:n3) = a(1:n1,1:n2,1:n3)
    end do
  else
    allocate ( return_value(1:n1,1:n2,1:n3,1:ncopies) )
    do l = 1,ncopies
      return_value(1:n1,1:n2,1:n3,l) = a(1:n1,1:n2,1:n3)
    end do
  end if
  !
end function spread_lk_rank3
