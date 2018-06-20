module octtree_example1_mod

! Example1 demonstrating use of octtree search for 1d array.

  use module_kind_types, only : wp
  use geovar, only : nr
  use octtree_mod

! Source data is 1d array:  xyzdata(1:nr,1:npts).
  type ,extends(octedat_type) :: octedat_example1_type
    real(wp) ,dimension(:,:) ,pointer :: xyzdata
  contains
    procedure :: dealloc     => octedat_example1_dealloc
    procedure :: evalxyz_pt  => octedat_example1_evalxyz_pt
    procedure :: evalxyz_pts => octedat_example1_evalxyz_pts
  end type

contains

! =========================================================================
  subroutine octedat_example1_dealloc(this,status)

! Deallocate memory used within octedat.

  class(octedat_example1_type) ,intent(inout) :: this
  integer                      ,intent(out)   :: status

  status = OCT_STATUS_OK

! If pointer octedat%xyzdata is allocated rather than assigned,
! then deallocate using:
!   if ( associated(this%xyzdata) ) deallocate(this%xyzdata)
! Since the user is responsible for manually setting octedat%xyzdata,
! let him manually dealloc it as appropriate.

  return
  end subroutine octedat_example1_dealloc
! =========================================================================
  subroutine octedat_example1_evalxyz_pt(this,idx,xyz,status)

! Evaluate xyz at a single point given by the index values in idx.

  class(octedat_example1_type) ,intent(in)  :: this
  integer  ,dimension(:)       ,intent(in)  :: idx
  real(wp) ,dimension(:)       ,intent(out) :: xyz
  integer                      ,intent(out) :: status

  integer :: k

  status = OCT_STATUS_DATA_ERROR

! Decode the index info.
  k  = idx(1) ! face pt.

! Get xyz.
  if ( k .ge. 1 .and. k .le. size(this%xyzdata,dim=2) ) then
    xyz = this%xyzdata(:,k)
    status = OCT_STATUS_OK
  end if

  return
  end subroutine octedat_example1_evalxyz_pt
! =========================================================================
  subroutine octedat_example1_evalxyz_pts(this,npts,idx,xyz,status)
    
! Evaluate xyz at multiple points given by the index values in idx.

  class(octedat_example1_type) ,intent(in)  :: this
  integer                      ,intent(in)  :: npts
  integer  ,dimension(:,:)     ,intent(in)  :: idx
  real(wp) ,dimension(:,:)     ,intent(out) :: xyz
  integer                      ,intent(out) :: status

  integer :: k

  status = OCT_STATUS_OK

  if ( npts .ge. 1 .and. npts .le. size(idx,dim=2) ) then
    do n = 1,npts
!     Decode the index info.
      k  = idx(1,n) ! face pt.

!     Get xyz.
      if ( k .ge. 1 .and. k .le. size(this%xyzdata,dim=2) ) then
        xyz(:,n) = this%xyzdata(:,k)
      else
        status = OCT_STATUS_DATA_ERROR
      end if

    end do ! n-loop
  else
    status = OCT_STATUS_DATA_ERROR
  end if

  return
  end subroutine octedat_example1_evalxyz_pts
! =========================================================================
  subroutine octtree_example1(src_xyzdata)

! Given an array of src_xyzdata, fill octtree and show how to use it.

  use module_kind_types
  use octstat_mod
  use octedat_mod
  use octnode_mod
  use octtree_mod
  implicit none

  real(wp)  ,target  ,dimension(:,:)  ,intent(in)    :: src_xyzdata

  type(octedat_example1_type)         ,target        :: octedat
  type(octtree_type)                                 :: octtree
  integer                                            :: status
  integer            ,dimension(:)    ,allocatable   :: ixyz_pt
  real(wp)           ,dimension(1:nr)                ::  xyz_pt
  real(wp)           ,dimension(1:nr,1:2)            :: bbox
  integer :: k,l,npts
  logical :: dealloc_octedat

! Init the encapsulated data.
  octedat%varname = "src_xyzdata(1:nr,1:npts)"
  octedat%vardesc = "Example of 1d vector array."
  octedat%xyzdata => src_xyzdata

! Compute bounding box of data.
  bbox(1:nr,1) =  huge(octedat%xyzdata(1,1))
  bbox(1:nr,2) = -huge(octedat%xyzdata(1,1))
  do k = 1,size(octedat%xyzdata,dim=2)
    do l = 1,nr
      bbox(l,1) = min( bbox(l,1), octedat%xyzdata(l,k) )
      bbox(l,2) = max( bbox(l,2), octedat%xyzdata(l,k) )
    end do
  end do

! Init the empty octtree.
  octtree%octedat => octedat
  call octtree_init(octtree,bbox,status)

!
! Fill the octtree with data.
!
  allocate(ixyz_pt(1))  ! 1d array only tracks 1 index
!
! Loop over all points in octtree%octedat.
!
  if ( .true. ) then
! Loop over all points in local octedat%xyzdata since we know
! that is what octtree%octedat points to (per above init).
  do k = 1,size(octedat%xyzdata,dim=2)
    ixyz_pt = k
!-- Compute xyz_pt directly or via function call.
    xyz_pt = octedat%xyzdata(:,k)
!   call octedat%evalxyz_pt(ixyz_pt,xyz_pt,status)
!-- Push ixyz_pt and xyz_pt data for this point to the octree
!-- or directly to the head octnode.
    if ( .true. ) then
    call octtree_addxyz(octtree,ixyz_pt,xyz_pt,status)
    else
    call octnode_addxyz(octtree%octnode,octtree%octroot, &
           octtree%octedat,ixyz_pt,xyz_pt,status)
    end  if
    if ( status .ne. OCT_STATUS_OK ) exit
  end do ! k-loop

  else
! Loop over points in octtree%octedat%xyzdata.
! This is harder because we need to tell the compiler that
! octtree%octedat is of octedat_example1_type (not octedat_type).
! Further, the select type construct doesn't work with type
! components, so we need to associate it with a local name.
  associate ( octtree_octedat => octtree%octedat )
  select type (octtree_octedat)
  class is (octedat_example1_type)
    do k = 1,size(octtree_octedat%xyzdata,dim=2)
      ixyz_pt = k
!--   Compute xyz_pt directly or via function call.
!     xyz_pt = octtree%octedat%xyzdata(:,k)  ! Fails.  Unknown field xyzdata.
      xyz_pt = octtree_octedat%xyzdata(:,k)  ! Works.
!     call octtree%octedat%evalxyz_pt(ixyz_pt,xyz_pt,status)
!--   Push ixyz_pt and xyz_pt data for this point to the octree
!--   or directly to the head octnode.
      if ( .true. ) then
      call octtree_addxyz(octtree,ixyz_pt,xyz_pt,status)
      else
      call octnode_addxyz(octtree%octnode,octtree%octroot, &
             octtree_octedat,ixyz_pt,xyz_pt,status)
      end if
      if ( status .ne. OCT_STATUS_OK ) exit
    enddo ! k-loop
  class default
    write(*,*) 'Internal error.  Should never reach this point.'
  end select
  end associate

!
!   Now ready to search the octtree.
!

! This section only uses the local octedat%xyzdata.  To use
! octtree%octedat%xyzdata follow similar steps to those above.
  do
    write(*,*) 'Enter x,y,z and I will find the nearest point in octtree:'
    read (*,*) xyz_pt(1),xyz_pt(2),xyz_pt(3)
    call octtree_find_nearest_xyz(octtree,xyz_pt, ixyz_pt,status)
    if ( status .eq. OCT_STATUS_OK ) then
!--   Compute xyz_pt directly or via function call.
      xyz_pt = octedat%xyzdata(:,k)
!     call octedat%evalxyz_pt(idyx_pt,xyz_pt,status)
      write(*,*) 'index = ',(ixyz_pt(k), k=1,size(ixyz_pt))
      write(*,*) 'xyz   = ',( xyz_pt(l), l=1,nr)
    else if ( status .eq. OCT_STATUS_NOT_FOUND ) then
      write(*,*) 'Point not found in octtree.'
    else if ( status .eq. OCT_STATUS_MEMORY_ERROR ) then
      write(*,*) 'Memory error.'
    else if ( status .eq. OCT_STATUS_DATA_ERROR ) then
    write(*,*) 'Data error.'
    end if
    write(*,*) 'Enter another point? (0/1=no/yes)'
    read (*,*) l
    if ( l .eq. 0 ) exit
  end do ! infinite-loop

  end if

!
! Cleanup.
!
  dealloc_octedat = .false.         ! Do this manually instead.
  call octtree_dealloc(octtree,dealloc_octedat,status)
  nullify(octtree%octedat)          ! Manually clear octtree%octedat
! This part is not necessary, since octedat goes out of scope,
! but if octedat were of type pointer with allocated memory
! then would need to deallocate.
  octedat%varname = ""
  octedat%vardesc = ""
  nullify(octedat%xyzdata)
  

  return
  end subroutine octtree_example1
! =========================================================================
end module octtree_example1_mod
