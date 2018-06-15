subroutine heap_sort_integer(n,iarrin,indx)
  !
  ! Routine to perform a heap sort on an array of integers
  !
  ! N      : number of items in list to be sorted
  ! IARRIN : integer array to be sorted
  ! INDX   : map from unsorted to sorted variables
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n
  integer,  intent(in) :: iarrin(1:n)
  integer, intent(out) :: indx(1:n)
  !
  !.. Local Scalars ..
  integer :: j,l,ir,indxt,i,iq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "heap_sort_integer"
  !
continue
  !
  if (n <= 0) return
  !
  do j = 1,n
    indx(j) = j
  end do
  !
  if (n == 1) return
  !
  l = n/2 + 1
  ir= n
  !
  do
    !
    if(l > 1)then
      l=l-1
      indxt = indx(l)
      iq    = iarrin(indxt)
    else
      indxt = indx(ir)
      iq    = iarrin(indxt)
      indx(ir) = indx(1)
      ir = ir - 1
      if(ir == 1)then
        indx(1) = indxt
        return
      end if
    end if
    !
    i=l
    j=l+l
    !
    do
      !
      if (j > ir) exit
      !
      if (j < ir) then ! separate if to avoid array bounds issues
        if (iarrin(indx(j)) < iarrin(indx(j+1))) j = j+1
      end if
      !
      if (iq < iarrin(indx(j))) then
        indx(i) = indx(j)
        i = j
        j = j+j
      else
        j = ir+1
      end if
      !
    end do
    !
    indx(i) = indxt
    !
  end do
  !
  ! Check for sort failure
  !
  do i = 1,n-1
    if (iarrin(indx(i+1)) < iarrin(indx(i))) then
      write (error_message,1)
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
    end if
  end do
  !
  ! Format Statements
  !
  1 format ("MERGE SORT ALGORITHM FAILED TO WORK PROPERLY!")
  !
end subroutine heap_sort_integer
!
!###############################################################################
!
subroutine merge_sort_integer(nsort,isort,array)
  !
  ! Routine to sort an array of integers in ascending order
  ! (via the Merge-Sort algorithm) using the values in the ISORT array
  !
  ! NSORT  : no. of items to be sorted
  ! ISORT  : the array to be used for the sort
  ! ARRAY  : integer array to be sorted
  !
  !.. Formal Arguments ..
  integer,    intent(in) :: nsort
  integer, intent(inout) :: isort(1:nsort)
  integer, intent(inout) :: array(1:nsort)
  !
  !.. Local Scalars ..
  integer     :: ierr,n
  logical(lk) :: list_sorted
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: tsort(:)
  integer, allocatable :: tarray(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "merge_sort_integer"
  !
continue
  !
  ! Check for a list that is already sorted
  !
  list_sorted = true
  do n = 1,nsort-1
    if (isort(n+1) < isort(n)) then
      list_sorted = fals
      exit
    end if
  end do
  if (list_sorted) return
  !
  allocate ( tsort(1:(nsort+1)/2) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tsort",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( tarray(1:(nsort+1)/2) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tarray",1,__LINE__,__FILE__,ierr,error_message)
  !
  call merge_sort_int(nsort,isort,tsort,array,tarray)
  !
  if (allocated(tsort))  deallocate( tsort )
  if (allocated(tarray)) deallocate( tarray )
  !
  ! Check for sort failure
  !
  do n = 1,nsort-1
    if (isort(n+1) < isort(n)) then
      write (error_message,1)
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
    end if
  end do
  !
  ! Format Statements
  !
  1 format ("MERGE SORT ALGORITHM FAILED TO WORK PROPERLY!")
  !
end subroutine merge_sort_integer
!
!###############################################################################
!
recursive subroutine merge_sort_int(nsort,isort,tsort,array,tarray)
  !
  ! Merge-Sort routine to sort an array of integers in ascending order
  ! using the values in the ISORT array
  !
  ! NSORT  : no. of items to be sorted
  ! ISORT  : the array to be used for the sort
  ! TSORT  : temporary array used for storage
  ! ARRAY  : integer array to be sorted
  ! TARRAY : temporary integer array to be sorted
  !
  !.. Formal Arguments ..
  integer,    intent(in) :: nsort
  integer, intent(inout) :: isort(1:nsort)
  integer,   intent(out) :: tsort(1:(nsort+1)/2)
  integer, intent(inout) :: array(1:nsort)
  integer,   intent(out) :: tarray(1:(nsort+1)/2)
  !
  !.. Local Scalars ..
  integer :: na,nb,ntmp
  !
continue
  !
  if (nsort <= 1) return
  !
  if (nsort == 2) then
    if (isort(1) > isort(2)) then
      ntmp = isort(1)
      isort(1) = isort(2)
      isort(2) = ntmp
      ntmp = array(1)
      array(1) = array(2)
      array(2) = ntmp
    end if
  end if
  !
  na = (nsort+1)/2
  nb = nsort - na
  !
  call merge_sort_int(na,isort,tsort,array,tarray)
  call merge_sort_int(nb,isort(na+1),tsort,array(na+1),tarray)
  !
  if (isort(na) > isort(na+1)) then
    tsort(1:na) = isort(1:na)
    tarray(1:na) = array(1:na)
    call merge_int(tsort,na,isort(na+1),nb,isort,nsort, &
                   tarray,array(na+1),array)
  end if
  !
end subroutine merge_sort_int
!
!###############################################################################
!
subroutine merge_int(a,na,b,nb,c,nc,a_array,b_array,c_array)
  !
  !.. Formal Arguments ..
  integer,    intent(in) :: na
  integer,    intent(in) :: nb
  integer,    intent(in) :: nc
  integer,    intent(in) :: a(1:na)
  integer,    intent(in) :: b(1:nb)
  integer, intent(inout) :: c(1:nc)
  integer,    intent(in) :: a_array(1:na)
  integer,    intent(in) :: b_array(1:nb)
  integer, intent(inout) :: c_array(1:nc)
  !
  !.. Local Scalars ..
  integer :: i,j,k
  !
continue
  !
  i = 1
  j = 1
  k = 1
  !
  do while (i <= na .and. j <= nb)
    if (a(i) <= b(j)) then
      c(k) = a(i)
      c_array(k) = a_array(i)
      i = i + 1
    else
      c(k) = b(j)
      c_array(k) = b_array(j)
      j = j + 1
    end if
    k = k + 1
  end do
  !
  do while (i <= na)
    c(k) = a(i)
    c_array(k) = a_array(i)
    i = i + 1
    k = k + 1
  end do
  !
end subroutine merge_int
!
!###############################################################################
!
