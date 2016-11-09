subroutine reallocate1_i4(array,new_size)
  integer, intent(in) :: new_size
  integer(i4), allocatable, dimension(:), intent(inout) :: array
  integer(i4), allocatable, dimension(:) :: temp_array
  integer :: ierr,num_copy
  character(len=*), parameter :: pname = "reallocate1_i4"
continue
  ! Allocate the temporary array to store the old contents of array
  allocate ( temp_array(1:new_size) , source=0_i4 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_array",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  ! If array is allocated, copy the old contents of array into the temporary
  ! array. If decreasing the size of array, make sure to limit the amount of
  ! information copied over from array to the size of temp_array.
  if (allocated(array)) then
    num_copy = min( size(temp_array) , size(array) )
    temp_array(1:num_copy) = array(1:num_copy)
  end if
  ! move_alloc does:
  !     1. deallocates array if it is allocated
  !     2. allocates array with the same type, type parameters,
  !        array bounds, and value as temp_array
  !     3. deallocates temp_array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_i4
!
!###############################################################################
!
subroutine reallocate1_i8(array,new_size)
  integer, intent(in) :: new_size
  integer(i8), allocatable, dimension(:), intent(inout) :: array
  integer(i8), allocatable, dimension(:) :: temp_array
  integer :: ierr,num_copy
  character(len=*), parameter :: pname = "reallocate1_i8"
continue
  allocate ( temp_array(1:new_size) , source=0_i8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_array",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  if (allocated(array)) then
    num_copy = min( size(temp_array) , size(array) )
    temp_array(1:num_copy) = array(1:num_copy)
  end if
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_i8
!
!###############################################################################
!
subroutine reallocate1_r4(array,new_size)
  integer, intent(in) :: new_size
  real(r4), allocatable, dimension(:), intent(inout) :: array
  real(r4), allocatable, dimension(:) :: temp_array
  integer :: ierr,num_copy
  character(len=*), parameter :: pname = "reallocate1_r4"
continue
  allocate ( temp_array(1:new_size) , source=0.0_r4 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_array",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  if (allocated(array)) then
    num_copy = min( size(temp_array) , size(array) )
    temp_array(1:num_copy) = array(1:num_copy)
  end if
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_r4
!
!###############################################################################
!
subroutine reallocate1_r8(array,new_size)
  integer, intent(in) :: new_size
  real(r8), allocatable, dimension(:), intent(inout) :: array
  real(r8), allocatable, dimension(:) :: temp_array
  integer :: ierr,num_copy
  character(len=*), parameter :: pname = "reallocate1_r8"
continue
  allocate ( temp_array(1:new_size) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_array",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  if (allocated(array)) then
    num_copy = min( size(temp_array) , size(array) )
    temp_array(1:num_copy) = array(1:num_copy)
  end if
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_r8
!
!###############################################################################
!
subroutine reallocate1_c(array,new_size)
  integer, intent(in) :: new_size
  character(len=*), allocatable, dimension(:), intent(inout) :: array
  character(len=len(array)), allocatable, dimension(:) :: temp_array
  integer :: ierr,num_copy
  character(len=*), parameter :: pname = "reallocate1_c"
continue
 !write (iout,1) size(array),new_size
  allocate ( temp_array(1:new_size) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_array",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  if (allocated(array)) then
    num_copy = min( size(temp_array) , size(array) )
    temp_array(1:num_copy) = array(1:num_copy)
  end if
  call move_alloc( from=temp_array , to=array )
  1 format ("size(array) = ",i0,"; new_size = ",i0)
end subroutine reallocate1_c
!
!###############################################################################
!
subroutine reallocate2_i4(array,new_size1,new_size2)
  !
  !.. Formal Arguments ..
  integer, optional,                       intent(in) :: new_size1
  integer, optional,                       intent(in) :: new_size2
  integer(i4), allocatable, dimension(:,:), intent(inout) :: array
  !
  !.. Local Scalars ..
  integer :: temp_size1,temp_size2,num_copy1,num_copy2,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(i4), allocatable, dimension(:,:) :: temp_array
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reallocate2_i4"
  !
continue
  !
  ! Get the old dimensions of array
  !
  temp_size1 = 1
  temp_size2 = 1
  if (allocated(array)) then
    temp_size1 = size(array,1)
    temp_size2 = size(array,2)
  end if
  !
  if (present(new_size1)) temp_size1 = new_size1
  if (present(new_size2)) temp_size2 = new_size2
  !
  ! Allocate the temporary array to store the old contents of array
  !
  allocate ( temp_array(1:temp_size1,1:temp_size2) , source=0_i4 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_array",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  !
  ! Copy the old contents of array into the temporary array
  !
  if (allocated(array)) then
    num_copy1 = min( size(temp_array,dim=1) , size(array,dim=1) )
    num_copy2 = min( size(temp_array,dim=2) , size(array,dim=2) )
    temp_array(1:num_copy1,1:num_copy2) = array(1:num_copy1,1:num_copy2)
  end if
  !
  ! move_alloc does:
  !     1. deallocates array
  !     2. allocates array with the same type, type parameters,
  !        array bounds, and value as temp_array
  !     3. deallocates temp_array
  !
  call move_alloc( from=temp_array , to=array )
  !
end subroutine reallocate2_i4
!
!###############################################################################
!
subroutine reallocate2_i8(array,new_size1,new_size2)
  !
  !.. Formal Arguments ..
  integer, optional,                       intent(in) :: new_size1
  integer, optional,                       intent(in) :: new_size2
  integer(i8), allocatable, dimension(:,:), intent(inout) :: array
  !
  !.. Local Scalars ..
  integer :: temp_size1,temp_size2,num_copy1,num_copy2,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(i8), allocatable, dimension(:,:) :: temp_array
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reallocate2_i8"
  !
continue
  !
  ! Get the old dimensions of array
  !
  temp_size1 = 1
  temp_size2 = 1
  if (allocated(array)) then
    temp_size1 = size(array,1)
    temp_size2 = size(array,2)
  end if
  !
  if (present(new_size1)) temp_size1 = new_size1
  if (present(new_size2)) temp_size2 = new_size2
  !
  ! Allocate the temporary array to store the old contents of array
  !
  allocate ( temp_array(1:temp_size1,1:temp_size2) , source=0_i8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_array",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  !
  ! Copy the old contents of array into the temporary array
  !
  if (allocated(array)) then
    num_copy1 = min( size(temp_array,dim=1) , size(array,dim=1) )
    num_copy2 = min( size(temp_array,dim=2) , size(array,dim=2) )
    temp_array(1:num_copy1,1:num_copy2) = array(1:num_copy1,1:num_copy2)
  end if
  !
  ! move_alloc does:
  !     1. deallocates array
  !     2. allocates array with the same type, type parameters,
  !        array bounds, and value as temp_array
  !     3. deallocates temp_array
  !
  call move_alloc( from=temp_array , to=array )
  !
end subroutine reallocate2_i8
!
!###############################################################################
!
subroutine reallocate2_r4(array,new_size1,new_size2)
  !
  !.. Formal Arguments ..
  integer,  optional,                       intent(in) :: new_size1
  integer,  optional,                       intent(in) :: new_size2
  real(r4), allocatable, dimension(:,:), intent(inout) :: array
  !
  !.. Local Scalars ..
  integer :: temp_size1,temp_size2,num_copy1,num_copy2,ierr
  !
  !.. Local Allocatable Arrays ..
  real(r4), allocatable, dimension(:,:) :: temp_array
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reallocate2_r4"
  !
continue
  !
  ! Get the old dimensions of array
  !
  temp_size1 = 1
  temp_size2 = 1
  if (allocated(array)) then
    temp_size1 = size(array,1)
    temp_size2 = size(array,2)
  end if
  !
  if (present(new_size1)) temp_size1 = new_size1
  if (present(new_size2)) temp_size2 = new_size2
  !
  ! Allocate the temporary array to store the old contents of array
  !
  allocate ( temp_array(1:temp_size1,1:temp_size2) , source=0.0_r4 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_array",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  !
  ! Copy the old contents of array into the temporary array
  !
  if (allocated(array)) then
    num_copy1 = min( size(temp_array,dim=1) , size(array,dim=1) )
    num_copy2 = min( size(temp_array,dim=2) , size(array,dim=2) )
    temp_array(1:num_copy1,1:num_copy2) = array(1:num_copy1,1:num_copy2)
  end if
  !
  ! move_alloc does:
  !     1. deallocates array
  !     2. allocates array with the same type, type parameters,
  !        array bounds, and value as temp_array
  !     3. deallocates temp_array
  !
  call move_alloc( from=temp_array , to=array )
  !
end subroutine reallocate2_r4
!
!###############################################################################
!
subroutine reallocate2_r8(array,new_size1,new_size2)
  !
  !.. Formal Arguments ..
  integer,  optional,                       intent(in) :: new_size1
  integer,  optional,                       intent(in) :: new_size2
  real(r8), allocatable, dimension(:,:), intent(inout) :: array
  !
  !.. Local Scalars ..
  integer :: temp_size1,temp_size2,num_copy1,num_copy2,ierr
  !
  !.. Local Allocatable Arrays ..
  real(r8), allocatable, dimension(:,:) :: temp_array
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reallocate2_r8"
  !
continue
  !
  ! Get the old dimensions of array
  !
  temp_size1 = 1
  temp_size2 = 1
  if (allocated(array)) then
    temp_size1 = size(array,1)
    temp_size2 = size(array,2)
  end if
  !
  if (present(new_size1)) temp_size1 = new_size1
  if (present(new_size2)) temp_size2 = new_size2
  !
  ! Allocate the temporary array to store the old contents of array
  !
  allocate ( temp_array(1:temp_size1,1:temp_size2) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_array",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  !
  ! Copy the old contents of array into the temporary array
  !
  if (allocated(array)) then
    num_copy1 = min( size(temp_array,dim=1) , size(array,dim=1) )
    num_copy2 = min( size(temp_array,dim=2) , size(array,dim=2) )
    temp_array(1:num_copy1,1:num_copy2) = array(1:num_copy1,1:num_copy2)
  end if
  !
  ! move_alloc does:
  !     1. deallocates array
  !     2. allocates array with the same type, type parameters,
  !        array bounds, and value as temp_array
  !     3. deallocates temp_array
  !
  call move_alloc( from=temp_array , to=array )
  !
end subroutine reallocate2_r8
!
!###############################################################################
!
