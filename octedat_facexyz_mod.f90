module octedat_facexyz_mod
  !
  ! octedat type for passing facexyz data to octtree search.
  !
  !.. Use Statements ..
  use module_kind_types, only : wp
#ifdef DEBUG_OCTEDAT_FACEXYZ
  use module_kind_types, only : iout  ! DEBUG
#endif
 !use flowvar_mod,       only : face_xyz
 !use geovar,            only : cell_t
  use octedat_mod
  use octstat_mod
  !
  implicit none
  !
  private
  !
  !
  ! Source data based on facexyz and cell array structures.
  type, public, extends(octedat_type) :: octedat_facexyz_type
  ! We will use the global variables:  facexyz, cell
  ! rather than pointer components to these global variables.
 !type(face_xyz), pointer :: facexyz(:)
 !type(cell_t),   pointer :: cell(:)
  !
  ! Cannot redefine MAXINDICES here.
  ! Must set manually after instance created or implement a constructor.
  ! integer :: MAXINDICES = 3
  ! Encode the index info using 3 values.
  !   nc  = idx(1)  ! cell no.
  !   ncf = idx(2)  ! cell face no.
  !   ifp = idx(3)  ! cell face facept no.
  contains
    procedure :: dealloc     => octedat_facexyz_dealloc
    procedure :: evalxyz_pt  => octedat_facexyz_evalxyz_pt
    procedure :: evalxyz_pts => octedat_facexyz_evalxyz_pts
  end type octedat_facexyz_type
  !
  ! MODULE PARAMETERS
  !
  real(wp), parameter :: ERROR_FINDING_XYZ = real(-9999,kind=wp)
  !
contains
!
!###############################################################################
!
subroutine octedat_facexyz_dealloc(this,istat)
  !
  ! Deallocate memory used within octedat.
  !
  !.. Passed Argument from Invoking Object ..
  class(octedat_facexyz_type), intent(inout) :: this
  !
  !.. Formal Arguments ..
  integer, intent(out) :: istat
  !
continue
  !
  istat = OCT_STATUS_OK
  !
  ! If pointer octedat% pointers are allocated rather than assigned,
  ! then deallocate using:
  !   if ( associated(this%facexyz) ) deallocate(this%facexyz)
  !   if ( associated(this%cell   ) ) deallocate(this%cell)
  ! Since the user is responsible for manually setting octedat%xyzdata,
  ! let him manually dealloc it as appropriate.
  !
end subroutine octedat_facexyz_dealloc
!
!###############################################################################
!
subroutine octedat_facexyz_evalxyz_pt(this,idx,xyz,istat)
  !
  ! Evaluate xyz at a single point given by the index values in idx.
  !
  !.. Use Statements ..
  use geovar,  only : cell,face
 !use flowvar_mod, only : facexyz_dt
  use flowvar_mod, only : facexyz
  !
  !.. Passed Argument from Invoking Object ..
  class(octedat_facexyz_type), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: idx(:)
  real(wp), intent(out) :: xyz(:)
  integer,  intent(out) :: istat
  !
  !.. Local Scalars ..
  integer :: ifp,nc,ncf,nf,f0
  !
continue
  !
  istat = OCT_STATUS_DATA_ERROR
  xyz(:) = ERROR_FINDING_XYZ
  !
  ! Decode the index info.
  !
  nc  = idx(1)  ! cell no.
  ncf = idx(2)  ! cell face no.
  ifp = idx(3)  ! cell face facept no.
  !
  ! Get xyz.
  !
  if (nc >= 1 .and. nc <= size(cell)) then
    if (ncf >= 1 .and. ncf <= size(cell(nc)%face)) then
      nf = cell(nc)%face(ncf)%idx  ! global face index
      f0 = face(nf)%beg_fp - 1     ! flxpt offset for this face
      xyz(:) = facexyz(:,f0+ifp)
     !xyz(:) = facexyz_dt(nf)%v(:,ifp)
      istat = OCT_STATUS_OK
    end if
  end if
  !
end subroutine octedat_facexyz_evalxyz_pt
!
!###############################################################################
!
subroutine octedat_facexyz_evalxyz_pts(this,npts,idx,xyz,istat)
  !
  ! Evaluate xyz at multiple points given by the index values in idx.
  !
  !.. Use Statements ..
  use geovar,  only : cell,face
 !use flowvar_mod, only : facexyz_dt
  use flowvar_mod, only : facexyz
  !
  !.. Passed Argument from Invoking Object ..
  class(octedat_facexyz_type), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: npts
  integer,   intent(in) :: idx(:,:)
  real(wp), intent(out) :: xyz(:,:)
  integer,  intent(out) :: istat
  !
  !.. Local Scalars ..
  integer :: ifp,n,nc,ncf,nf,f0
  !
continue
  !
  istat = OCT_STATUS_OK
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTEDAT_FACEXYZ
  write (iout,*) 'DEBUG_BEG[octedat_facexyz_evalxyz_pts]'
  write (iout,*) '  npts      = ',npts
  write (iout,*) '  size(idx,dim=1) = ',size(idx,dim=1)
  write (iout,*) '  size(idx,dim=2) = ',size(idx,dim=2)
  write (iout,*) 'DEBUG_END[octedat_facexyz_evalxyz_pts]'
  flush (iout)
#endif
  !-DEBUG_END:
  !
  !
  if (npts >= 1 .and. npts <= size(idx,dim=2)) then
    !
    do n = 1,npts
      !
      ! Decode the index info.
      !
      nc  = idx(1,n)  ! cell no.
      ncf = idx(2,n)  ! cell face no.
      ifp = idx(3,n)  ! cell face facept no.
      !
      ! Get xyz.
      !
      if (nc >= 1 .and. nc <= size(cell)) then
        if (ncf >= 1 .and. ncf <= size(cell(nc)%face)) then
          nf = cell(nc)%face(ncf)%idx  ! global face index
          f0 = face(nf)%beg_fp - 1     ! flxpt offset for this face
          xyz(:,n) = facexyz(:,f0+ifp)
         !xyz(:,n) = facexyz_dt(nf)%v(:,ifp)
        else
          istat = OCT_STATUS_DATA_ERROR
          xyz(:,n) = ERROR_FINDING_XYZ
        end if
      else
        istat = OCT_STATUS_DATA_ERROR
        xyz(:,n) = ERROR_FINDING_XYZ
      end if
      !
    end do ! n-loop
    !
  else
    !
    istat = OCT_STATUS_DATA_ERROR
    xyz(:,:) = ERROR_FINDING_XYZ
    !
  end if
  !
end subroutine octedat_facexyz_evalxyz_pts
!
!###############################################################################
!
end module octedat_facexyz_mod
