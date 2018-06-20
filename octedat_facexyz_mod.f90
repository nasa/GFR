module octedat_facexyz_mod

! octedat type for passing facexyz data to octtree search.

  use module_kind_types, only : wp
#ifdef DEBUG_OCTEDAT_FACEXYZ
  use module_kind_types, only : iout  ! DEBUG
#endif
  use flowvar,           only : face_xyz
  use geovar,     only : cell_t
  use octedat_mod

! Source data based on facexyz and cell array structures.
  type ,extends(octedat_type) :: octedat_facexyz_type
!   We will use the global variables:  facexyz, cell
!   rather than pointer components to these global variables.
!   type(face_xyz) ,dimension(:) ,pointer :: facexyz
!   type(cell_t)   ,dimension(:) ,pointer :: cell
!
!   Cannot redefine MAXINDICES here.
!   Must set manually after instance created or implement a constructor.
!   integer :: MAXINDICES = 3
!   Encode the index info using 3 values.
!     nc  = idx(1)  ! cell no.
!     ncf = idx(2)  ! cell face no.
!     ifp = idx(3)  ! cell face facept no.
  contains
    procedure :: dealloc     => octedat_facexyz_dealloc
    procedure :: evalxyz_pt  => octedat_facexyz_evalxyz_pt
    procedure :: evalxyz_pts => octedat_facexyz_evalxyz_pts
  end type

contains

! =========================================================================
  subroutine octedat_facexyz_dealloc(this,status)

! Deallocate memory used within octedat.

  class(octedat_facexyz_type) ,intent(inout) :: this
  integer                     ,intent(out)   :: status

  status = OCT_STATUS_OK

! If pointer octedat% pointers are allocated rather than assigned,
! then deallocate using:
!   if ( associated(this%facexyz) ) deallocate(this%facexyz)
!   if ( associated(this%cell   ) ) deallocate(this%cell)
! Since the user is responsible for manually setting octedat%xyzdata,
! let him manually dealloc it as appropriate.

  return
  end subroutine octedat_facexyz_dealloc
! =========================================================================
  subroutine octedat_facexyz_evalxyz_pt(this,idx,xyz,status)

! Evaluate xyz at a single point given by the index values in idx.

  use geovar,  only : cell
  use flowvar, only : facexyz
  implicit none

  class(octedat_facexyz_type)  ,intent(in)  :: this
  integer  ,dimension(:)       ,intent(in)  :: idx
  real(wp) ,dimension(:)       ,intent(out) :: xyz
  integer                      ,intent(out) :: status

  integer :: ifp,nc,ncf,nf

  status = OCT_STATUS_DATA_ERROR
  xyz(:) = -9999

! Decode the index info.
  nc  = idx(1)  ! cell no.
  ncf = idx(2)  ! cell face no.
  ifp = idx(3)  ! cell face facept no.

! Get xyz.
  if ( nc .ge. 1 .and. nc .le. size(cell) ) then
    if ( ncf .ge. 1 .and. ncf .le. size(cell(nc)%face) ) then
      nf = cell(nc)%face(ncf)%idx  ! global face index
      xyz(:) = facexyz(nf)%v(:,ifp)
      status = OCT_STATUS_OK
    end if
  end if

  return
  end subroutine octedat_facexyz_evalxyz_pt
! =========================================================================
  subroutine octedat_facexyz_evalxyz_pts(this,npts,idx,xyz,status)

! Evaluate xyz at multiple points given by the index values in idx.

  use geovar,  only : cell
  use flowvar, only : facexyz
  implicit none

  class(octedat_facexyz_type)  ,intent(in)  :: this
  integer                      ,intent(in)  :: npts
  integer  ,dimension(:,:)     ,intent(in)  :: idx
  real(wp) ,dimension(:,:)     ,intent(out) :: xyz
  integer                      ,intent(out) :: status

  integer :: ifp,n,nc,ncf,nf

  status = OCT_STATUS_OK

!-DEBUG_BEG:
#ifdef DEBUG_OCTEDAT_FACEXYZ
  write(iout,*) 'DEBUG_BEG[octedat_facexyz_evalxyz_pts]'
  write(iout,*) '  npts      = ',npts
  write(iout,*) '  size(idx,dim=1) = ',size(idx,dim=1)
  write(iout,*) '  size(idx,dim=2) = ',size(idx,dim=2)
  write(iout,*) 'DEBUG_END[octedat_facexyz_evalxyz_pts]'
  flush(iout)
#endif
!-DEBUG_END:

  if ( npts .ge. 1 .and. npts .le. size(idx,dim=2) ) then
    do n = 1,npts
!     Decode the index info.
      nc  = idx(1,n)  ! cell no.
      ncf = idx(2,n)  ! cell face no.
      ifp = idx(3,n)  ! cell face facept no.

!     Get xyz.
      if ( nc .ge. 1 .and. nc .le. size(cell) ) then
        if ( ncf .ge. 1 .and. ncf .le. size(cell(nc)%face) ) then
          nf = cell(nc)%face(ncf)%idx  ! global face index
          xyz(:,n) = facexyz(nf)%v(:,ifp)
        else
          status = OCT_STATUS_DATA_ERROR
          xyz(:,n) = -9999
        end if
      else
        status = OCT_STATUS_DATA_ERROR
        xyz(:,n) = -9999
      end if
    end do ! n-loop
  else
    status = OCT_STATUS_DATA_ERROR
    xyz(:,:) = -9999
  end if

  return
  end subroutine octedat_facexyz_evalxyz_pts
! =========================================================================
! =========================================================================
end module octedat_facexyz_mod
