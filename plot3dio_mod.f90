module plot3dio_mod

! Module for reading generic plot3d files.
!   grid/solution/function
!   single/multi-block
!   2d/3d
!   formatted/unformatted
!   native/lil/big-endian
!   single/double-precision
!   iblanking or not

! -------------------------------------------------------------------------
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: plot3dio_read
  public :: plot3dio_ijk_matches
  !
! -------------------------------------------------------------------------

! Plot3d "container" for storing all details for a particular file.
  type, public :: type_plot3dinfo
    character(len=150) :: file
    character(len= 15) :: form
    character(len= 15) :: convert
    character(len= 15) :: type
    character(len= 15) :: blank
    character(len= 15) :: mgrid
    character(len= 15) :: dim
    character(len= 15) :: iprec
    character(len= 15) :: rprec
!   integer ,allocatable,dimension(:,:),intent(inout) :: ijk
!   integer ,allocatable,dimension(:)  ,intent(inout) :: nvar
!   integer ,allocatable,dimension(:)  ,intent(inout) :: iblank
!   real(wp),allocatable,dimension(:,:),intent(inout) :: data
  end type type_plot3dinfo

  character(len=*), parameter :: ind = "                    "
  !
  integer, save :: indx = 2
  !
  character(len=*), parameter :: &
     errstr = " ######################  ERROR!  ######################"

  ! Create a more friendly name to the front-end routine.
  interface plot3dio_read
    module procedure plot3dio_read_unknown_file
  end interface

contains
!
!###############################################################################
!
subroutine plot3dio_read_unknown_file(file,form,convert, &
                                      type,blank, &
                                      mgrid,dim, &
                                      iprec,rprec, &
                                      ijk,nvar,iblank,data, &
                                      iostat,iomsg,plot3dinfo)

! Generic plot3d file reader for when plot3d file parameters are unknown
! or only partially known.
! Allowed options shown below.
! read_plot3d(file       = fname, &
!             form       = "UNFORMATTED|FORMATTED", &
!             convert    = "NATIVE|LITTLE_ENDIAN|BIG_ENDIAN", &
!             type       = "XYZ|SOL|FUN", &
!             blank      = "NO|YES", &
!             mgrid      = "YES|NO", &
!             dim        = "3D|2D|1D", &
!             iprec      = "SINGLE|DOUBLE", &
!             rprec      = "SINGLE|DOUBLE", &
!             ijk        = allocatable_ijk(1:3,1:nblock), &
!             nvar       = allocatable_nvar(1:nblock), &
!             iblank     = allocatable_iblank(1:nblock) &
!             data       = allocatable_data(1:maxvar,1:nblock), &
!             iostat     = ierr, &
!             iomsg      = error_message )

  !.. Argument Variables ..
  character(len=*)                   ,intent(in)    :: file
  character(len=*)         ,optional ,intent(in)    :: form
  character(len=*)         ,optional ,intent(in)    :: convert
  character(len=*)         ,optional ,intent(in)    :: type
  character(len=*)         ,optional ,intent(in)    :: blank
  character(len=*)         ,optional ,intent(in)    :: mgrid
  character(len=*)         ,optional ,intent(in)    :: dim
  character(len=*)         ,optional ,intent(in)    :: iprec
  character(len=*)         ,optional ,intent(in)    :: rprec
  integer ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer ,allocatable,dimension(:)  ,intent(inout) ,optional :: iblank
  real(wp),allocatable,dimension(:,:),intent(inout) :: data
  integer                            ,intent(out)   :: iostat
  character(len=*)                   ,intent(out)   :: iomsg
  type(type_plot3dinfo)    ,optional ,intent(out)   :: plot3dinfo

  !.. Local Variables ..
  character(len=15) ,allocatable ,dimension(:) :: &
    p3d_forms, p3d_converts, p3d_types, p3d_blanks, p3d_mgrids, &
    p3d_dims,  p3d_iprecs,   p3d_rprecs
  character(len=15) :: &
    p3d_form, p3d_convert, p3d_type, p3d_blank, p3d_mgrid, &
    p3d_dim,  p3d_iprec,   p3d_rprec
  integer :: &
    ip3d_form, ip3d_convert, ip3d_type, ip3d_blank, ip3d_mgrid, &
    ip3d_dim,  ip3d_iprec,   ip3d_rprec
  integer :: iop3d, iostat_header, iostat_body
  integer :: loop_main_var, loop_header_var, loop_body_var
  integer ,allocatable,dimension(:) :: iblank_tmp
  logical :: file_exists

! -------------------------------------------------------------------------
continue

  ! Initialize the default options in the order to be tried.
  ! Some options must be tried in a certain order to detect things.
  !   Must try 3d before 2d.
  !   For unformatted, must try double before single.
  ! Other options are ordered based on liklihood.

  ! For formatted files, it is difficult to distinguish precision.
  ! So values will be read as double to preserve accuracy.

  ! This section uses F2003 Auto reallocation to resize char arrays.
  !                 123456789012345   123456789012345   123456789012345
  p3d_forms    = [ 'UNFORMATTED    ','FORMATTED      '                   ]
  p3d_converts = [ 'NATIVE         ','LITTLE_ENDIAN  ','BIG_ENDIAN     ' ]
  p3d_types    = [ 'XYZ            ','SOL            ','FUN            ' ]
  p3d_blanks   = [ 'YES            ','NO             '                   ]
  p3d_mgrids   = [ 'YES            ','NO             '                   ]
  p3d_dims     = [ '3D             ','2D             ','1D             ' ]
  p3d_iprecs   = [ 'DOUBLE         ','SINGLE         '                   ]
  p3d_rprecs   = [ 'DOUBLE         ','SINGLE         '                   ]

  ! Override defaults with any settings provided via argument variables.
  p3d_form    = ''; if (present(form   )) p3d_form    = trim(uppercase(form   ))
  p3d_convert = ''; if (present(convert)) p3d_convert = trim(uppercase(convert))
  p3d_type    = ''; if (present(type   )) p3d_type    = trim(uppercase(type   ))
  p3d_blank   = ''; if (present(blank  )) p3d_blank   = trim(uppercase(blank  ))
  p3d_mgrid   = ''; if (present(mgrid  )) p3d_mgrid   = trim(uppercase(mgrid  ))
  p3d_dim     = ''; if (present(dim    )) p3d_dim     = trim(uppercase(dim    ))
  p3d_iprec   = ''; if (present(iprec  )) p3d_iprec   = trim(uppercase(iprec  ))
  p3d_rprec   = ''; if (present(rprec  )) p3d_rprec   = trim(uppercase(rprec  ))

  if (len_trim(p3d_form)    > 0) p3d_forms    = [p3d_form   ]
  if (len_trim(p3d_convert) > 0) p3d_converts = [p3d_convert]
  if (len_trim(p3d_type)    > 0) p3d_types    = [p3d_type   ]
  if (len_trim(p3d_blank)   > 0) p3d_blanks   = [p3d_blank  ]
  if (len_trim(p3d_mgrid)   > 0) p3d_mgrids   = [p3d_mgrid  ]
  if (len_trim(p3d_dim)     > 0) p3d_dims     = [p3d_dim    ]
  if (len_trim(p3d_iprec)   > 0) p3d_iprecs   = [p3d_iprec  ]
  if (len_trim(p3d_rprec)   > 0) p3d_rprecs   = [p3d_rprec  ]

  ! Return Plot3d settings if requested.
  ! Set what is known at this point:  either user arg or empty.
  if ( present(plot3dinfo) ) then
    plot3dinfo%file    = file
    plot3dinfo%form    = p3d_form
    plot3dinfo%convert = p3d_convert
    plot3dinfo%type    = p3d_type
    plot3dinfo%mgrid   = p3d_mgrid
    plot3dinfo%dim     = p3d_dim
    plot3dinfo%iprec   = p3d_iprec
    plot3dinfo%blank   = p3d_blank
    plot3dinfo%rprec   = p3d_rprec
  end if

  ! Defensively set the default to some unknown error.
  iostat_header = 1
  iostat_body = 1
  iostat = 1
  iomsg  = 'Unknown error while reading plot3d file.'

  ! Going to wipe iblank anyway, so do it now.
  if ( present(iblank) ) then
    if ( allocated(iblank) ) deallocate (iblank)
  end if

  ! Main loop is just a container, so we can exit out.
  loop_main: &
  do loop_main_var = 1,1

    ! Make sure the file exists before trying to read it multiple times.
    inquire (file=trim(file),exist=file_exists)
    if ( .not. file_exists ) then
      iostat = 1
      iomsg = 'Plot3d file does not exist.'
      exit loop_main
    end if

    ! Header loop is just a container, so we can exit out.
    loop_header: &
    do loop_header_var = 1,1

    ! Loop over various plot3d option combos.
    do ip3d_convert = 1,size(p3d_converts)
    do ip3d_form    = 1,size(p3d_forms)

      p3d_convert = p3d_converts( ip3d_convert )
      p3d_form    = p3d_forms(    ip3d_form    )

!-DEBUG_BEG:
!-DEBUG_CMT:  Print info on what we are trying to read.
        if (mypnum == glb_root) then
          write (iout,*) ind(1:indx-2),'Reading plot3d header from: ', &
                                       trim(adjustl(file))
          write (iout,*) ind(1:indx),'Form    = ',trim(adjustl(p3d_form   ))
          write (iout,*) ind(1:indx),'Convert = ',trim(adjustl(p3d_convert))
          flush (iout)
        end if
!-DEBUG_END:

      ! Warning:  convert= is a non-standard fortran extension to the
      !           open command but is pretty common among compilers.
      open (newunit=iop3d,file=trim(file), &
            form=trim(p3d_form), &
            status="old",action="read",position="rewind", &
            convert=trim(p3d_convert), &
            iostat=iostat, iomsg=iomsg)
!-DEBUG_BEG:
!-DEBUG_CMT:  Print info on what we are trying to read.
        if (mypnum == glb_root) then
          if ( iostat == 0 ) then
            write (iout,*) ind(1:indx),'Open SUCCEEDED'
          else
            write (iout,*) ind(1:indx),'Open FAILED'
            write (iout,*) ind(1:indx+2),'iostat = ',iostat
            write (iout,*) ind(1:indx+2),'errmsg = ',trim(adjustl(iomsg))
          end if
          flush (iout)
        end if
!-DEBUG_END:
      ! If file opened ok, proceed.
      if ( iostat == 0 ) then
        indx = indx + 2

        ! Loop over remaining plot3d options for header.
        ! Loop order based on likelihood of guessing format (ie, popularity).
        ! Hardest to guess option should be inner-most loop.
        do ip3d_type  = 1,size(p3d_types)
        do ip3d_dim   = 1,size(p3d_dims)
        do ip3d_iprec = 1,size(p3d_iprecs)
        do ip3d_mgrid = 1,size(p3d_mgrids)

        p3d_type  = p3d_types(  ip3d_type  )
        p3d_dim   = p3d_dims(   ip3d_dim   )
        p3d_iprec = p3d_iprecs( ip3d_iprec )
        p3d_mgrid = p3d_mgrids( ip3d_mgrid )

!-DEBUG_BEG:
!-DEBUG_CMT:  Print info on what we are trying to read.
        if (mypnum == glb_root) then
!         write (iout,*) 'Reading plot3d header from: ',trim(file)
!         write (iout,*) '  Form    = ',trim(p3d_form   )
!         write (iout,*) '  Convert = ',trim(p3d_convert)
          write (iout,*) ind(1:indx),'Type    = ',trim(adjustl(p3d_type   ))
          write (iout,*) ind(1:indx),'Mgrid   = ',trim(adjustl(p3d_mgrid  ))
          write (iout,*) ind(1:indx),'Dim     = ',trim(adjustl(p3d_dim    ))
          write (iout,*) ind(1:indx),'Iprec   = ',trim(adjustl(p3d_iprec  ))
          flush (iout)
        end if
!-DEBUG_END:
        ! Read Plot3d header info.
        call plot3dio_read_header(unit=iop3d,file=file,form=p3d_form, &
                                  type=p3d_type, &
                                  mgrid=p3d_mgrid,dim=p3d_dim, &
                                  iprec=p3d_iprec, &
                                  ijk=ijk,nvar=nvar, &
                                  iostat=iostat,iomsg=iomsg)
!-DEBUG_BEG:
!-DEBUG_CMT:  Print info on what we are trying to read.
        if (mypnum == glb_root) then
          write (iout,*) ind(1:indx+2),'iostat  = ',iostat
          write (iout,*) ind(1:indx+2),'iomsg   = ',trim(adjustl(iomsg))
          flush (iout)
        end if
!-DEBUG_END:
        ! If header read ok, exit loop_header.
        if ( iostat == 0 ) exit loop_header
        ! Else, rewind and try again on next pass.
        rewind (iop3d)
        end do
        end do
        end do
        end do
        ! Must close file before reopening on next loop pass.
        close (iop3d)
      end if
    end do
    end do
    end do loop_header
    indx = indx - 2
    iostat_header = iostat
    ! If could not read header, then skip reading the body.
    if ( iostat /= 0 ) then
      close (iop3d)
      exit loop_main
    end if

    loop_body: &
    do loop_body_var = 1,1
    do ip3d_rprec = 1,size(p3d_rprecs)
    do ip3d_blank = 1,size(p3d_blanks)

      p3d_rprec = p3d_rprecs( ip3d_rprec )
      p3d_blank = p3d_blanks( ip3d_blank )

      call plot3dio_read_data(unit=iop3d,file=file,form=p3d_form, &
                              type=p3d_type,blank=p3d_blank, &
                              mgrid=p3d_mgrid,dim=p3d_dim, &
                              iprec=p3d_iprec,rprec=p3d_rprec, &
                              ijk=ijk,nvar=nvar,iblank=iblank_tmp, &
                              data=data, &
                              iostat=iostat,iomsg=iomsg)
!-DEBUG_BEG:
!-DEBUG_CMT:  Print info on what we tried to read.
      if (mypnum == glb_root) then
        write (iout,*) ind(1:indx-2),'Reading plot3d body from: ', &
                                     trim(adjustl(file))
        write (iout,*) ind(1:indx),'Blank   = ',trim(adjustl(p3d_blank  ))
        write (iout,*) ind(1:indx),'Rprec   = ',trim(adjustl(p3d_rprec  ))
        write (iout,*) ind(1:indx),'iostat  = ',iostat
        write (iout,'(1x,a,a,a)') ind(1:indx),'iomsg   = ',trim(adjustl(iomsg))
        flush (iout)
      end if
!-DEBUG_END:
      ! If body read ok, exit loop_body.
      if ( iostat == 0 ) then
        close (iop3d)
        exit loop_body
      end if
      ! Else, rewind, read header, and try again on next pass.
      rewind (iop3d)
      call plot3dio_read_header(unit=iop3d,file=file,form=p3d_form, &
                                type=p3d_type, &
                                mgrid=p3d_mgrid,dim=p3d_dim, &
                                iprec=p3d_iprec, &
                                ijk=ijk,nvar=nvar, &
                                iostat=iostat,iomsg=iomsg)
      ! Should never have trouble re-reading the header.
      ! But if so, exit out of here.
      if ( iostat /= 0 ) then
        close (iop3d)
        exit loop_body
      end if

    end do
    end do
    end do loop_body
    iostat_body = iostat
  end do loop_main

  if ( iostat == 0 ) then
    ! File was read successfully.
    iomsg = ''
    if ( present(iblank) ) call move_alloc(from=iblank_tmp , to=iblank)
    ! Make sure data will be returned allocated
    if (.not. allocated(data)) then
      iostat = 1
      write (iomsg,100)
    end if
  else
    ! Read failed.  Modify the error message to include the file name.
!-DEBUG_BEG:
    if (mypnum == glb_root) then
      write (iout,*) 'DEBUG_BEG[plot3dio_read_unknown_file]'
      write (iout,*) '  Unable to read Plot3d file'
      write (iout,*) '  ',trim(adjustl(file))
      write (iout,*) '  ',trim(adjustl(iomsg))
      write (iout,*) 'DEBUG_END[plot3dio_read_unknown_file]'
      flush (iout)
    end if
!-DEBUG_CMT: gcc "Fortran runtime error: End of file" on next regular line.
!-DEBUG_END:
    write (iomsg,'(a,a,/,a)') 'Unable to read Plot3d file ' // &
                              trim(adjustl(file)) // '.', trim(adjustl(iomsg))
  end if

  ! Return Plot3d settings if requested.
  ! Set what was finally found to work.
  if ( present(plot3dinfo) ) then
    plot3dinfo%file    = trim(adjustl(file))
    if ( iostat_header == 0 ) then
      plot3dinfo%form    = p3d_form
      plot3dinfo%convert = p3d_convert
      plot3dinfo%type    = p3d_type
      plot3dinfo%mgrid   = p3d_mgrid
      plot3dinfo%dim     = p3d_dim
      plot3dinfo%iprec   = p3d_iprec
    end if
    if ( iostat_body == 0 ) then
      plot3dinfo%blank   = p3d_blank
      plot3dinfo%rprec   = p3d_rprec
    end if
  end if

  return
  !
  ! Format Statements
  !
  100 format ("Plot3D read seemed to complete successfully, ", &
              "but data array returned deallocated!")
  !
end subroutine plot3dio_read_unknown_file
!
!###############################################################################
!
subroutine plot3dio_read_known_file(file,form,convert, &
                                  type,blank, &
                                  mgrid,dim, &
                                  iprec,rprec, &
                                  ijk,nvar,iblank,data, &
                                  iostat,iomsg)

! Generic plot3d file reader for when plot3d file parameters are known.
! Allowed options shown below.  Default setting is listed first.
! read_plot3d(file       = fname, &
!             form       = "UNFORMATTED|FORMATTED", &
!             convert    = "NATIVE|LITTLE_ENDIAN|BIG_ENDIAN", &
!             type       = "XYZ|SOL|FUN", &
!             blank      = "NO|YES", &
!             mgrid      = "YES|NO", &
!             dim        = "3D|2D|1D", &
!             iprec      = "SINGLE|DOUBLE", &
!             rprec      = "SINGLE|DOUBLE", &
!             ijk        = allocatable_ijk(1:3,1:nblock), &
!             nvar       = allocatable_nvar(1:nblock), &
!             iblank     = allocatable_iblank(1:nblock) &
!             data       = allocatable_data(1:maxvar,1:nblock), &
!             iostat     = ierr, &
!             iomsg      = error_message )

  !.. Argument Variables ..
  character(len=*)                   ,intent(in)    :: file
  character(len=*)         ,optional ,intent(in)    :: form
  character(len=*)         ,optional ,intent(in)    :: convert
  character(len=*)         ,optional ,intent(in)    :: type
  character(len=*)         ,optional ,intent(in)    :: blank
  character(len=*)         ,optional ,intent(in)    :: mgrid
  character(len=*)         ,optional ,intent(in)    :: dim
  character(len=*)         ,optional ,intent(in)    :: iprec
  character(len=*)         ,optional ,intent(in)    :: rprec
  integer ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer ,allocatable,dimension(:)  ,intent(inout) ,optional :: iblank
  real(wp),allocatable,dimension(:,:),intent(inout) :: data
  integer                            ,intent(out)   :: iostat
  character(len=*)                   ,intent(out)   :: iomsg

  !.. Local Variables ..
  character(len=:) ,allocatable :: p3d_form
  character(len=:) ,allocatable :: p3d_convert
  character(len=:) ,allocatable :: p3d_type
  character(len=:) ,allocatable :: p3d_blank
  character(len=:) ,allocatable :: p3d_mgrid
  character(len=:) ,allocatable :: p3d_dim
  character(len=:) ,allocatable :: p3d_iprec
  character(len=:) ,allocatable :: p3d_rprec
  character(len=1000)           :: p3d_msg
  integer                       :: p3d_ierr
  integer                       :: iop3d
  integer ,allocatable,dimension(:) :: iblank_tmp

! -------------------------------------------------------------------------
continue

  ! Default options, if none provided.
  p3d_form       = 'UNFORMATTED'
  p3d_convert    = 'NATIVE'
  p3d_type       = 'XYZ'
  p3d_blank      = 'NO'
  p3d_mgrid      = 'YES'
  p3d_dim        = '3D'
  p3d_iprec      = 'SINGLE'
  p3d_rprec      = 'SINGLE'

  ! Override defaults with provided options.
  if ( present(form      ) ) p3d_form       = trim(uppercase(form))
  if ( present(convert   ) ) p3d_convert    = trim(uppercase(convert))
  if ( present(type      ) ) p3d_type       = trim(uppercase(type))
  if ( present(blank     ) ) p3d_blank      = trim(uppercase(blank))
  if ( present(mgrid     ) ) p3d_mgrid      = trim(uppercase(mgrid))
  if ( present(dim       ) ) p3d_dim        = trim(uppercase(dim))
  if ( present(iprec     ) ) p3d_iprec      = trim(uppercase(iprec))
  if ( present(rprec     ) ) p3d_rprec      = trim(uppercase(rprec))

  ! Check for valid options.
  ! Skip form + convert.  Let open statement handle those errors.
  ! Skip iprec + rprec.  Let header and body read routines handle those.
  iostat = 1  ! Assume an error until proven wrong.
  iomsg  = 'Unknown error.'
  if      ( index('|XYZ|SOL|FUN|','|'//trim(p3d_type)//'|') == 0 ) then
    iomsg = 'Plot3d type option must be one of "XYZ|SOL|FUN".'
  else if ( index('|YES|NO|','|'//trim(p3d_blank)//'|') == 0 ) then
    iomsg  = 'Plot3d blank option must be one of "YES|NO".'
  else if ( index('|YES|NO|','|'//trim(p3d_mgrid)//'|') == 0 ) then
    iomsg  = 'Plot3d mgrid option must be one of "YES|NO".'
  else if ( index('|1D|2D|3D|','|'//trim(p3d_dim)//'|') == 0 ) then
    iomsg  = 'Plot3d ndim option must be one of "1D|2D|3D".'
!  else if ( index('|SINGLE|DOUBLE|','|'//trim(p3d_iprec)//'|') == 0 ) then
!    iomsg  = 'Plot3d iprec option must be one of "SINGLE|DOUBLE".'
!  else if ( index('|SINGLE|DOUBLE|','|'//trim(p3d_rprec)//'|') == 0 ) then
!    iomsg  = 'Plot3d rprec option must be one of "SINGLE|DOUBLE".'
  else
    !
    ! Options are ok.  Proceed with file read.
    !
    ! Warning:  convert= is a non-standard fortran extension to the
    !           open command but is pretty common among compilers.
    open (newunit=iop3d,file=trim(file), &
          form=trim(p3d_form), &
          status="old",action="read",position="rewind", &
          convert=trim(p3d_convert), &
          iostat=iostat, iomsg=iomsg)
    if ( iostat == 0 ) then
      ! Read Plot3d header info.
      if (mypnum == glb_root) then
        write (iout,*) '  Reading Plot3d header.'  ! DEBUG OUTPUT
      end if
      call plot3dio_read_header(unit=iop3d,file=file,form=p3d_form, &
                                type=p3d_type, &
                                mgrid=p3d_mgrid,dim=p3d_dim, &
                                iprec=p3d_iprec, &
                                ijk=ijk,nvar=nvar, &
                                iostat=iostat,iomsg=iomsg)
      if ( iostat == 0 ) then
        ! Deallocate iblank to help keep free memory contiguous.
        if ( present(iblank) ) then
          if ( allocated(iblank) ) deallocate (iblank)
        end if
        ! Read Plot3d body data.
        if (mypnum == glb_root) then
          write (iout,*) '  Reading Plot3d body.'  ! DEBUG OUTPUT
        end if
        call plot3dio_read_data(unit=iop3d,file=file,form=p3d_form, &
                                type=p3d_type,blank=p3d_blank, &
                                mgrid=p3d_mgrid,dim=p3d_dim, &
                                iprec=p3d_iprec,rprec=p3d_rprec, &
                                ijk=ijk,nvar=nvar,iblank=iblank_tmp, &
                                data=data, &
                                iostat=p3d_ierr,iomsg=p3d_msg)
        close (iop3d,iostat=iostat,iomsg=iomsg)
        ! Ignore the close status.  Use Plot3d status instead.
        iostat = p3d_ierr
        iomsg  = p3d_msg
        ! Copy out iblank, if requested.
        if ( p3d_ierr == 0 .and. present(iblank) ) &
          call move_alloc(from=iblank_tmp , to=iblank)
      end if
    end if
  end if

  return
end subroutine plot3dio_read_known_file
!
!###############################################################################
!
subroutine plot3dio_read_header(unit,file,form, &
                                type, &
                                mgrid,dim, &
                                iprec, &
                                ijk,nvar, &
                                iostat,iomsg)

! Read plot3d header info.  Front-end to more specific routines.

  !.. Argument Variables ..
  integer          ,intent(in)  :: unit
  character(len=*) ,intent(in)  :: file
  character(len=*) ,intent(in)  :: form
  character(len=*) ,intent(in)  :: type
  character(len=*) ,intent(in)  :: mgrid
  character(len=*) ,intent(in)  :: dim
  character(len=*) ,intent(in)  :: iprec
  integer ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer          ,intent(out) :: iostat
  character(len=*) ,intent(out) :: iomsg

  !.. Local Variables ..
  integer(i8) ,allocatable,dimension(:,:) :: ijk_dp
  integer(i8) ,allocatable,dimension(:)   :: nvar_dp
! -------------------------------------------------------------------------
continue

  iostat = 1
  iomsg  = 'Unknown error.'

  if ( iprec == 'SINGLE' ) then
    ! Read header info directly into i4 arrays.
    if ( form == 'FORMATTED' ) then
      call plot3dio_sp_read_header_form(unit,file,form, &
                                        type, &
                                        mgrid,dim, &
                                        iprec,ijk,nvar, &
                                        iostat,iomsg)
    else
      call plot3dio_sp_read_header_ufrm(unit,file,form, &
                                        type, &
                                        mgrid,dim, &
                                        iprec,ijk,nvar, &
                                        iostat,iomsg)
    end if
  else if ( iprec == 'DOUBLE' ) then
    ! Read header info into i8 temp arrays.
    if ( form == 'FORMATTED' ) then
      call plot3dio_dp_read_header_form(unit,file,form, &
                                        type, &
                                        mgrid,dim, &
                                        iprec,ijk_dp,nvar_dp, &
                                        iostat,iomsg)
    else
      call plot3dio_dp_read_header_ufrm(unit,file,form, &
                                        type, &
                                        mgrid,dim, &
                                        iprec,ijk_dp,nvar_dp, &
                                        iostat,iomsg)
    end if
    ! Copy header data from i8 tmp arrays into return arrays.
    if ( iostat == 0 ) then
      ijk  = int(ijk_dp)  ! Using F2003 Auto-Reallocation
      nvar = int(nvar_dp) ! Using F2003 Auto-Reallocation
    end if
    if ( allocated(ijk_dp) ) deallocate (ijk_dp)
    if ( allocated(nvar_dp) ) deallocate (nvar_dp)
  else
    iostat = 1
    iomsg  = 'Plot3d iprec option must be one of "SINGLE|DOUBLE".'
  end if

  return
end subroutine plot3dio_read_header
!
!###############################################################################
!
subroutine plot3dio_read_data(unit,file,form, &
                              type,blank, &
                              mgrid,dim, &
                              iprec,rprec, &
                              ijk,nvar,data,iblank, &
                              iostat,iomsg)

! Read plot3d body data info.  Front-end to more specific routines.

  !.. Argument Variables ..
  integer                            ,intent(in)    :: unit
  character(len=*)                   ,intent(in)    :: file
  character(len=*)                   ,intent(in)    :: form
  character(len=*)                   ,intent(in)    :: type
  character(len=*)                   ,intent(in)    :: blank
  character(len=*)                   ,intent(in)    :: mgrid
  character(len=*)                   ,intent(in)    :: dim
  character(len=*)                   ,intent(in)    :: iprec
  character(len=*)                   ,intent(in)    :: rprec
  integer ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer ,allocatable,dimension(:)  ,intent(inout) :: iblank
  real(wp),allocatable,dimension(:,:),intent(inout) :: data
  integer                            ,intent(out)   :: iostat
  character(len=*)                   ,intent(out)   :: iomsg

  !.. Local Variables ..
  real(sp),allocatable,dimension(:,:)               :: data_sp
  real(dp),allocatable,dimension(:,:)               :: data_dp
  integer :: ierr,mxvar,nblk,ndim,nnode

  !.. Local Parameters ..
  character(len=*), parameter :: pname = "plot3dio_read_data"

! -------------------------------------------------------------------------
continue

  iostat = 0  ! Assume no error until proven otherwise.
  iomsg  = ''

  if ( allocated(ijk) .and. allocated(nvar) ) then
    mxvar = maxval(nvar)
    ndim  = size(ijk,dim=1)
    nblk  = size(ijk,dim=2)
    nnode = sum( product( ijk(1:ndim,1:nblk) , dim=1 ) )

    ! Allocate storage for the return data, if not correct size.
    if ( allocated(data) ) then
      if ( size(data,dim=1) /= mxvar .or. &
           size(data,dim=2) /= nblk  ) then
        deallocate ( data, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"data",2,__LINE__,__FILE__,ierr,error_message)
      end if
    end if
    if ( .not. allocated(data) ) then
      allocate ( data(1:mxvar,1:nnode), source=0.0_wp, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"data",1,__LINE__,__FILE__,ierr,error_message)
    end if

    ! Storage for data_sp|data_dp will be allocated within the called subs.
    ! Copy data_sp|data_dp back into the return data.

    if ( rprec == 'SINGLE' ) then
      if ( form == 'FORMATTED' ) then
        call plot3dio_sp_read_data_form(unit,file,form,type,blank,mgrid,dim, &
                                        iprec,rprec,ijk,nvar,iblank,data_sp, &
                                        iostat,iomsg)
      else
        call plot3dio_sp_read_data_ufrm(unit,file,form,type,blank,mgrid,dim, &
                                        iprec,rprec,ijk,nvar,iblank,data_sp, &
                                        iostat,iomsg)
      end if
      if ( allocated(data_sp) ) then
        if ( iostat == 0 ) data(:,:) = real(data_sp(:,:),kind(wp))
        deallocate (data_sp)
      end if
    else if ( rprec == 'DOUBLE' ) then
      if ( form == 'FORMATTED' ) then
        call plot3dio_dp_read_data_form(unit,file,form,type,blank,mgrid,dim, &
                                        iprec,rprec,ijk,nvar,iblank,data_dp, &
                                        iostat,iomsg)
      else
        call plot3dio_dp_read_data_ufrm(unit,file,form,type,blank,mgrid,dim, &
                                        iprec,rprec,ijk,nvar,iblank,data_dp, &
                                        iostat,iomsg)
      end if
      if ( allocated(data_dp) ) then
        if ( iostat == 0 ) data(:,:) = real(data_dp(:,:),kind(wp))
        deallocate (data_dp)
      end if
    else
      iostat = 1
      iomsg  = 'Plot3d rprec option must be one of "SINGLE|DOUBLE".'
    end if
  else
    iostat = 1
    iomsg  = 'Plot3d header info must be set before reading body data.'
  end if

  return
end subroutine plot3dio_read_data
!
!###############################################################################
!
subroutine plot3dio_ijk_matches(ijk_xyz,ijk_sol,iostat,iomsg)

! Check whether two sets of ijk values match.

  !.. Argument Variables ..
  integer ,allocatable,dimension(:,:),intent(in)    :: ijk_xyz
  integer ,allocatable,dimension(:,:),intent(in)    :: ijk_sol
  integer                            ,intent(out)   :: iostat
  character(len=*)                   ,intent(out)   :: iomsg

  !.. Local Variables ..
! -------------------------------------------------------------------------
continue

  iostat = 1
  iomsg = 'Mismatch in size of Plot3d grid and solution files.'

  if ( allocated(ijk_xyz) .and. allocated(ijk_sol) ) then
    if ( size(ijk_xyz) == size(ijk_sol) ) then
      if ( all(ijk_xyz == ijk_sol) ) then
        iostat = 0
        iomsg = ''
      end if
    end if
  end if

  return
end subroutine plot3dio_ijk_matches
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! Subroutines for reading formatted single-precision Plot3D files
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine plot3dio_sp_read_header_form(unit,file,form, &
                                        type, &
                                        mgrid,dim, &
                                        iprec,ijk,nvar, &
                                        iostat,iomsg)

! Read plot3d header info for formatted single-precision files.

  !.. Argument Variables ..
  integer          ,intent(in)  :: unit
  character(len=*) ,intent(in)  :: file
  character(len=*) ,intent(in)  :: form
  character(len=*) ,intent(in)  :: type
  character(len=*) ,intent(in)  :: mgrid
  character(len=*) ,intent(in)  :: dim
  character(len=*) ,intent(in)  :: iprec
  integer(i4) ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer(i4) ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer          ,intent(out) :: iostat
  character(len=*) ,intent(out) :: iomsg

  !.. Local Variables ..
  integer(i4) nblk_i4
  integer nblk
  integer l,n,ndim
  integer ierr

  !.. Local Parameters ..
  character(len=*), parameter :: pname = "plot3dio_sp_read_header_form"
  integer :: huge_nblk = 100000
  integer :: huge_ijk  = 100000

! -------------------------------------------------------------------------
continue

  iostat = 0  ! Assume no error until proven otherwise.
  iomsg  = ''

  if ( form == 'FORMATTED' ) then
    read (dim,'(i1)') ndim
    if ( mgrid == 'NO' ) then
      nblk = 1
    else
      read (unit,*,iostat=iostat,iomsg=iomsg) nblk_i4
      nblk = int(nblk_i4)
    end if
    if ( iostat /= 0 ) then
      ! Already have iomsg from read command.
    else if ( nblk <= 0 ) then
      iostat = 1
      iomsg  = 'Number of blocks in header of Plot3d file must be positive.'
    else if ( nblk > huge_nblk ) then
      iostat = 1
      iomsg  = 'Number of blocks in header of Plot3d file is too large.'
    else
      if ( allocated(ijk) ) then
        deallocate ( ijk, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"ijk",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( ijk(1:3,1:nblk), source=1_i4, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"ijk",1,__LINE__,__FILE__,ierr,error_message)

      if ( allocated(nvar) ) then
        deallocate ( nvar, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"nvar",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( nvar(1:nblk), source=1_i4, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"nvar",1,__LINE__,__FILE__,ierr,error_message)

      if ( type == 'XYZ' ) then
        read (unit,*,iostat=iostat,iomsg=iomsg) &
          ((ijk(l,n),             l=1,ndim), n=1,nblk)
        nvar(:) = ndim
        if (ndim == 3) then
          if (nblk == 1) then
            if (ijk(1,nblk) == 1 .and. all(ijk(2:3,nblk)>1)) then
              iostat = 1
              iomsg  = 'This is most likely a 2D multi-block grid and ' // &
                       'not a 3D single-block grid with only 1 i-plane.'
            end if
          end if
        end if
      else if ( type == 'SOL' ) then
        read (unit,*,iostat=iostat) &
          ((ijk(l,n),             l=1,ndim), n=1,nblk)
        nvar(:) = ndim+2
      else if ( type == 'FUN' ) then
        read (unit,*,iostat=iostat) &
          ((ijk(l,n),nvar(n), l=1,ndim), n=1,nblk)
        if ( any(nvar(:) <= 0 ) ) then
          iostat = 1
          iomsg  = 'Number of variables in header of ' // &
                   'Plot3d function file must be positive.'
        else
          iostat = 0  ! header values read ok
          iomsg  = ''
        end if
      end if  ! file type

      if ( iostat /= 0 ) then
        ! Read failed.  Error already saved in iostat and iomsg.
      else if ( any(ijk(:,:) <= 0 ) ) then
        iostat = 1
        iomsg  = 'Array bounds in header of Plot3d file must be positive.'
      else if ( any(ijk(:,:) > huge_ijk ) ) then
        iostat = 1
        iomsg  = 'Array bounds in header of Plot3d file too large.'
      else
        iostat = 0  ! header values read ok
        iomsg  = ''
      end if
    end if  ! nblk check
  else
    iostat = 1
    iomsg  = 'Routine '//pname//' cannot read '//trim(form)//' files.'
  end if  ! plot3d format

9 continue
  return

end subroutine plot3dio_sp_read_header_form
!
!###############################################################################
!
subroutine plot3dio_sp_read_data_form(unit,file,form, &
                                      type,blank, &
                                      mgrid,dim, &
                                      iprec,rprec, &
                                      ijk,nvar,iblank,data, &
                                      iostat,iomsg)

! Read plot3d body data info for formatted double-precision files.

  !.. Argument Variables ..
  integer                            ,intent(in)    :: unit
  character(len=*)                   ,intent(in)    :: file
  character(len=*)                   ,intent(in)    :: form
  character(len=*)                   ,intent(in)    :: type
  character(len=*)                   ,intent(in)    :: blank
  character(len=*)                   ,intent(in)    :: mgrid
  character(len=*)                   ,intent(in)    :: dim
  character(len=*)                   ,intent(in)    :: iprec
  character(len=*)                   ,intent(in)    :: rprec
  integer ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer ,allocatable,dimension(:)  ,intent(inout) :: iblank
  real(sp),allocatable,dimension(:,:),intent(inout) :: data
  integer                            ,intent(out)   :: iostat
  character(len=*)                   ,intent(out)   :: iomsg

  !.. Local variables ..
  integer ierr
  integer i0,i,l,iblk,ndim,nblk,nnode,npts,mxvar
  integer(i8) ,allocatable,dimension(:) :: iblank_i8
  real(sp) :: fsmach,alpha,re,time

  !.. Local Parameters ..
  character(len=*), parameter :: pname = "plot3dio_sp_read_data_form"

! -------------------------------------------------------------------------
continue

  iostat = 0  ! Assume no error until proven otherwise.
  iomsg  = ''

  if ( form == 'FORMATTED' ) then

    ndim  = size(ijk,dim=1)
    nblk  = size(ijk,dim=2)
    nnode = sum( product( ijk(1:ndim,1:nblk) , dim=1 ) )
    mxvar = maxval(nvar(:))

    if ( allocated(data) ) then
      deallocate ( data, stat=ierr, errmsg=error_message )
      call alloc_error(pname,"data",2,__LINE__,__FILE__,ierr,error_message)
    end if
    allocate ( data(1:mxvar,1:nnode), source=0.0_sp, &
              stat=ierr, errmsg=error_message )
    call alloc_error(pname,"data",1,__LINE__,__FILE__,ierr,error_message)

    if ( type == 'XYZ' .and. blank == 'YES' ) then
      if ( allocated(iblank) ) then
        deallocate ( iblank, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"iblank",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( iblank(1:nnode), source=1, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"iblank",1,__LINE__,__FILE__,ierr,error_message)

      if ( iprec == 'SINGLE' ) then
        i0 = 0
        do iblk = 1,nblk
          npts = product(ijk(1:ndim,iblk))
          read (unit,*,iostat=iostat,iomsg=iomsg,end=9) &
            ((data(l,i0+i), i=1,npts),l=1,nvar(iblk)), &
             (iblank(i0+i), i=1,npts)
          i0 = i0 + npts
        end do
      else
        allocate ( iblank_i8(1:nnode), source=1_i8, &
                  stat=ierr, errmsg=error_message )
        call alloc_error(pname,"iblank_i8",1,__LINE__,__FILE__, &
                         ierr,error_message)
        i0 = 0
        do iblk = 1,nblk
          npts = product(ijk(1:ndim,iblk))
          read (unit,*,iostat=iostat,iomsg=iomsg,end=9) &
            ((data(l,i0+i), i=1,npts),l=1,nvar(iblk)), &
             (iblank_i8(i0+i), i=1,npts)
          i0 = i0 + npts
        end do
        iblank(:) = int(iblank_i8(:))
        deallocate (iblank_i8)
      end if

    else
      i0 = 0
      do iblk = 1,nblk
        npts = product(ijk(1:ndim,iblk))
        if ( type == 'SOL' ) then
          read (unit,iostat=iostat,iomsg=iomsg,end=9) fsmach,alpha,re,time
        end if
        read (unit,*,iostat=iostat,iomsg=iomsg,end=9) &
          ((data(l,i0+i), i=1,npts),l=1,nvar(iblk))
        i0 = i0 + npts
      end do
    end if
  else
    iostat = 1
    iomsg  = 'Routine '//pname//' cannot read '//trim(form)//' files.'

  end if

9 continue
  return
end subroutine plot3dio_sp_read_data_form
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! Subroutines for reading unformatted single-precision Plot3D files
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine plot3dio_sp_read_header_ufrm(unit,file,form, &
                                        type, &
                                        mgrid,dim, &
                                        iprec,ijk,nvar, &
                                        iostat,iomsg)

! Read plot3d header info for unformatted single-precision files.

  !.. Argument Variables ..
  integer          ,intent(in)  :: unit
  character(len=*) ,intent(in)  :: file
  character(len=*) ,intent(in)  :: form
  character(len=*) ,intent(in)  :: type
  character(len=*) ,intent(in)  :: mgrid
  character(len=*) ,intent(in)  :: dim
  character(len=*) ,intent(in)  :: iprec
  integer(i4) ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer(i4) ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer          ,intent(out) :: iostat
  character(len=*) ,intent(out) :: iomsg

  !.. Local Variables ..
  integer(i4) nblk_i4
  integer nblk
  integer l,n,ndim
  integer ierr

  !.. Local Parameters ..
  character(len=*), parameter :: pname = "plot3dio_sp_read_header_ufrm"
  integer :: huge_nblk = 100000
  integer :: huge_ijk  = 100000

! -------------------------------------------------------------------------
continue

  iostat = 0  ! Assume no error until proven otherwise.
  iomsg  = ''

  if ( form /= 'FORMATTED' ) then
    read (dim,'(i1)') ndim
    if ( mgrid == 'NO' ) then
      nblk = 1
    else
      read (unit,iostat=iostat,iomsg=iomsg) nblk_i4
      nblk = int(nblk_i4)
    end if
    if ( iostat /= 0 ) then
      ! Already have iomsg from read command.
    else if ( nblk <= 0 ) then
      iostat = 1
      iomsg  = 'Number of blocks in header of Plot3d file must be positive.'
    else if ( nblk > huge_nblk ) then
      iostat = 1
      iomsg  = 'Number of blocks in header of Plot3d file is too large.'
    else
      if ( allocated(ijk) ) then
        deallocate ( ijk, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"ijk",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( ijk(1:3,1:nblk), source=1_i4, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"ijk",1,__LINE__,__FILE__,ierr,error_message)

      if ( allocated(nvar) ) then
        deallocate ( nvar, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"nvar",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( nvar(1:nblk), source=1_i4, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"nvar",1,__LINE__,__FILE__,ierr,error_message)

      if ( type == 'XYZ' ) then
        read (unit,iostat=iostat,iomsg=iomsg) &
          ((ijk(l,n),             l=1,ndim), n=1,nblk)
        nvar(:) = ndim
      else if ( type == 'SOL' ) then
        read (unit,iostat=iostat) &
          ((ijk(l,n),             l=1,ndim), n=1,nblk)
        nvar(:) = ndim+2
      else if ( type == 'FUN' ) then
        read (unit,iostat=iostat) &
          ((ijk(l,n),nvar(n), l=1,ndim), n=1,nblk)
        if ( any(nvar(:) <= 0 ) ) then
          iostat = 1
          iomsg  = 'Number of variables in header of ' // &
                   'Plot3d function file must be positive.'
        else
          iostat = 0  ! header values read ok
          iomsg  = ''
        end if
      end if  ! file type

      if ( iostat /= 0 ) then
        ! Read failed.  Error already saved in iostat and iomsg.
      else if ( any(ijk(:,:) <= 0 ) ) then
        iostat = 1
        iomsg  = 'Array bounds in header of Plot3d file must be positive.'
      else if ( any(ijk(:,:) > huge_ijk ) ) then
        iostat = 1
        iomsg  = 'Array bounds in header of Plot3d file too large.'
      else
        iostat = 0  ! header values read ok
        iomsg  = ''
      end if
    end if  ! nblk check
  else
    iostat = 1
    iomsg  = 'Routine '//pname//' cannot read '//trim(form)//' files.'
  end if  ! plot3d format

9 continue
  return

end subroutine plot3dio_sp_read_header_ufrm
!
!###############################################################################
!
subroutine plot3dio_sp_read_data_ufrm(unit,file,form, &
                                      type,blank, &
                                      mgrid,dim, &
                                      iprec,rprec, &
                                      ijk,nvar,iblank,data, &
                                      iostat,iomsg)

! Read plot3d body data info for formatted double-precision files.

  !.. Argument Variables ..
  integer                            ,intent(in)    :: unit
  character(len=*)                   ,intent(in)    :: file
  character(len=*)                   ,intent(in)    :: form
  character(len=*)                   ,intent(in)    :: type
  character(len=*)                   ,intent(in)    :: blank
  character(len=*)                   ,intent(in)    :: mgrid
  character(len=*)                   ,intent(in)    :: dim
  character(len=*)                   ,intent(in)    :: iprec
  character(len=*)                   ,intent(in)    :: rprec
  integer ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer ,allocatable,dimension(:)  ,intent(inout) :: iblank
  real(sp),allocatable,dimension(:,:),intent(inout) :: data
  integer                            ,intent(out)   :: iostat
  character(len=*)                   ,intent(out)   :: iomsg

  !.. Local variables ..
  integer ierr
  integer i0,i,l,iblk,ndim,nblk,nnode,npts,mxvar
  integer(i8) ,allocatable,dimension(:) :: iblank_i8
  real(sp) :: fsmach,alpha,re,time

  !.. Local Parameters ..
  character(len=*), parameter :: pname = "plot3dio_sp_read_data_ufrm"

! -------------------------------------------------------------------------
continue

  iostat = 0  ! Assume no error until proven otherwise.
  iomsg  = ''

  if ( form /= 'FORMATTED' ) then

    ndim  = size(ijk,dim=1)
    nblk  = size(ijk,dim=2)
    nnode = sum( product( ijk(1:ndim,1:nblk) , dim=1 ) )
    mxvar = maxval(nvar(:))

    if ( allocated(data) ) then
      deallocate ( data, stat=ierr, errmsg=error_message )
      call alloc_error(pname,"data",2,__LINE__,__FILE__,ierr,error_message)
    end if
    allocate ( data(1:mxvar,1:nnode), source=0.0_sp, &
              stat=ierr, errmsg=error_message )
    call alloc_error(pname,"data",1,__LINE__,__FILE__,ierr,error_message)

    if ( type == 'XYZ' .and. blank == 'YES' ) then
      if ( allocated(iblank) ) then
        deallocate ( iblank, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"iblank",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( iblank(1:nnode), source=1, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"iblank",1,__LINE__,__FILE__,ierr,error_message)

      if ( iprec == 'SINGLE' ) then
        i0 = 0
        do iblk = 1,nblk
          npts = product(ijk(1:ndim,iblk))
          read (unit,iostat=iostat,iomsg=iomsg,end=9) &
            ((data(l,i0+i), i=1,npts),l=1,nvar(iblk)), &
             (iblank(i0+i), i=1,npts)
          i0 = i0 + npts
        end do
      else
        allocate ( iblank_i8(1:nnode), source=1_i8, &
                  stat=ierr, errmsg=error_message )
        call alloc_error(pname,"iblank_i8",1,__LINE__,__FILE__, &
                         ierr,error_message)
        i0 = 0
        do iblk = 1,nblk
          npts = product(ijk(1:ndim,iblk))
          read (unit,iostat=iostat,iomsg=iomsg,end=9) &
            ((data(l,i0+i), i=1,npts),l=1,nvar(iblk)), &
             (iblank_i8(i0+i), i=1,npts)
          i0 = i0 + npts
        end do
        iblank(:) = int(iblank_i8(:))
        deallocate (iblank_i8)
      end if

    else
      i0 = 0
      do iblk = 1,nblk
        npts = product(ijk(1:ndim,iblk))
        if ( type == 'SOL' ) then
          read (unit,iostat=iostat,iomsg=iomsg,end=9) fsmach,alpha,re,time
        end if
        read (unit,iostat=iostat,iomsg=iomsg,end=9) &
          ((data(l,i0+i), i=1,npts),l=1,nvar(iblk))
        i0 = i0 + npts
      end do
    end if
  else
    iostat = 1
    iomsg  = 'Routine '//pname//' cannot read '//trim(form)//' files.'

  end if

9 continue
  return
end subroutine plot3dio_sp_read_data_ufrm
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! Subroutines for reading formatted double-precision Plot3D files
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine plot3dio_dp_read_header_form(unit,file,form, &
                                        type, &
                                        mgrid,dim, &
                                        iprec,ijk,nvar, &
                                        iostat,iomsg)

! Read plot3d header info for formatted double-precision files.

  !.. Argument Variables ..
  integer          ,intent(in)  :: unit
  character(len=*) ,intent(in)  :: file
  character(len=*) ,intent(in)  :: form
  character(len=*) ,intent(in)  :: type
  character(len=*) ,intent(in)  :: mgrid
  character(len=*) ,intent(in)  :: dim
  character(len=*) ,intent(in)  :: iprec
  integer(i8) ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer(i8) ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer          ,intent(out) :: iostat
  character(len=*) ,intent(out) :: iomsg

  !.. Local Variables ..
  integer(i8) nblk_i8
  integer nblk
  integer l,n,ndim
  integer ierr

  !.. Local Parameters ..
  character(len=*), parameter :: pname = "plot3dio_dp_read_header_form"
  integer :: huge_nblk = 100000
  integer :: huge_ijk  = 100000

! -------------------------------------------------------------------------
continue

  iostat = 0  ! Assume no error until proven otherwise.
  iomsg  = ''

  if ( form == 'FORMATTED' ) then
    read (dim,'(i1)') ndim
    if ( mgrid == 'NO' ) then
      nblk = 1
    else
      read (unit,*,iostat=iostat,iomsg=iomsg) nblk_i8
      nblk = int(nblk_i8)
    end if
    if ( iostat /= 0 ) then
      ! Already have iomsg from read command.
    else if ( nblk <= 0 ) then
      iostat = 1
      iomsg  = 'Number of blocks in header of Plot3d file must be positive.'
    else if ( nblk > huge_nblk ) then
      iostat = 1
      iomsg  = 'Number of blocks in header of Plot3d file is too large.'
    else
      if ( allocated(ijk) ) then
        deallocate ( ijk, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"ijk",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( ijk(1:3,1:nblk), source=1_i8, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"ijk",1,__LINE__,__FILE__,ierr,error_message)

      if ( allocated(nvar) ) then
        deallocate ( nvar, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"nvar",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( nvar(1:nblk), source=1_i8, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"nvar",1,__LINE__,__FILE__,ierr,error_message)

      if ( type == 'XYZ' ) then
        read (unit,*,iostat=iostat,iomsg=iomsg) &
          ((ijk(l,n),             l=1,ndim), n=1,nblk)
        nvar(:) = ndim
        if (ndim == 3) then
          if (nblk == 1) then
            if (ijk(1,nblk) == 1 .and. all(ijk(2:3,nblk)>1)) then
              iostat = 1
              iomsg  = 'This is most likely a 2D multi-block grid and ' // &
                       'not a 3D single-block grid with only 1 i-plane.'
            end if
          end if
        end if
      else if ( type == 'SOL' ) then
        read (unit,*,iostat=iostat) &
          ((ijk(l,n),             l=1,ndim), n=1,nblk)
        nvar(:) = ndim+2
      else if ( type == 'FUN' ) then
        read (unit,*,iostat=iostat) &
          ((ijk(l,n),nvar(n), l=1,ndim), n=1,nblk)
        if ( any(nvar(:) <= 0 ) ) then
          iostat = 1
          iomsg  = 'Number of variables in header of ' // &
                   'Plot3d function file must be positive.'
        else
          iostat = 0  ! header values read ok
          iomsg  = ''
        end if
      end if  ! file type

      if ( iostat /= 0 ) then
        ! Read failed.  Error already saved in iostat and iomsg.
      else if ( any(ijk(:,:) <= 0 ) ) then
        iostat = 1
        iomsg  = 'Array bounds in header of Plot3d file must be positive.'
      else if ( any(ijk(:,:) > huge_ijk ) ) then
        iostat = 1
        iomsg  = 'Array bounds in header of Plot3d file too large.'
      else
        iostat = 0  ! header values read ok
        iomsg  = ''
      end if
    end if  ! nblk check
  else
    iostat = 1
    iomsg  = 'Routine '//pname//' cannot read '//trim(form)//' files.'
  end if  ! plot3d format

9 continue
  return

end subroutine plot3dio_dp_read_header_form
!
!###############################################################################
!
subroutine plot3dio_dp_read_data_form(unit,file,form, &
                                      type,blank, &
                                      mgrid,dim, &
                                      iprec,rprec, &
                                      ijk,nvar,iblank,data, &
                                      iostat,iomsg)

! Read plot3d body data info for formatted double-precision files.

  !.. Argument Variables ..
  integer                            ,intent(in)    :: unit
  character(len=*)                   ,intent(in)    :: file
  character(len=*)                   ,intent(in)    :: form
  character(len=*)                   ,intent(in)    :: type
  character(len=*)                   ,intent(in)    :: blank
  character(len=*)                   ,intent(in)    :: mgrid
  character(len=*)                   ,intent(in)    :: dim
  character(len=*)                   ,intent(in)    :: iprec
  character(len=*)                   ,intent(in)    :: rprec
  integer ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer ,allocatable,dimension(:)  ,intent(inout) :: iblank
  real(dp),allocatable,dimension(:,:),intent(inout) :: data
  integer                            ,intent(out)   :: iostat
  character(len=*)                   ,intent(out)   :: iomsg

  !.. Local variables ..
  integer ierr
  integer i0,i,l,iblk,ndim,nblk,nnode,npts,mxvar
  integer(i8) ,allocatable,dimension(:) :: iblank_i8
  real(dp) :: fsmach,alpha,re,time

  !.. Local Parameters ..
  character(len=*), parameter :: pname = "plot3dio_dp_read_data_form"

! -------------------------------------------------------------------------
continue

  iostat = 0  ! Assume no error until proven otherwise.
  iomsg  = ''

  if ( form == 'FORMATTED' ) then

    ndim  = size(ijk,dim=1)
    nblk  = size(ijk,dim=2)
    nnode = sum( product( ijk(1:ndim,1:nblk) , dim=1 ) )
    mxvar = maxval(nvar(:))

    if ( allocated(data) ) then
      deallocate ( data, stat=ierr, errmsg=error_message )
      call alloc_error(pname,"data",2,__LINE__,__FILE__,ierr,error_message)
    end if
    allocate ( data(1:mxvar,1:nnode), source=0.0_dp, &
              stat=ierr, errmsg=error_message )
    call alloc_error(pname,"data",1,__LINE__,__FILE__,ierr,error_message)

    if ( type == 'XYZ' .and. blank == 'YES' ) then
      if ( allocated(iblank) ) then
        deallocate ( iblank, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"iblank",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( iblank(1:nnode), source=1, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"iblank",1,__LINE__,__FILE__,ierr,error_message)

      if ( iprec == 'SINGLE' ) then
        i0 = 0
        do iblk = 1,nblk
          npts = product(ijk(1:ndim,iblk))
          read (unit,*,iostat=iostat,iomsg=iomsg,end=9) &
            ((data(l,i0+i), i=1,npts),l=1,nvar(iblk)), &
             (iblank(i0+i), i=1,npts)
          i0 = i0 + npts
        end do
      else
        allocate ( iblank_i8(1:nnode), source=1_i8, &
                  stat=ierr, errmsg=error_message )
        call alloc_error(pname,"iblank_i8",1,__LINE__,__FILE__, &
                         ierr,error_message)
        i0 = 0
        do iblk = 1,nblk
          npts = product(ijk(1:ndim,iblk))
          read (unit,*,iostat=iostat,iomsg=iomsg,end=9) &
            ((data(l,i0+i), i=1,npts),l=1,nvar(iblk)), &
             (iblank_i8(i0+i), i=1,npts)
          i0 = i0 + npts
        end do
        iblank(:) = int(iblank_i8(:))
        deallocate (iblank_i8)
      end if

    else
      i0 = 0
      do iblk = 1,nblk
        npts = product(ijk(1:ndim,iblk))
        if ( type == 'SOL' ) then
          read (unit,iostat=iostat,iomsg=iomsg,end=9) fsmach,alpha,re,time
        end if
        read (unit,*,iostat=iostat,iomsg=iomsg,end=9) &
          ((data(l,i0+i), i=1,npts),l=1,nvar(iblk))
        i0 = i0 + npts
      end do
    end if
  else
    iostat = 1
    iomsg  = 'Routine '//pname//' cannot read '//trim(form)//' files.'

  end if

9 continue
  return
end subroutine plot3dio_dp_read_data_form
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! Subroutines for reading unformatted double-precision Plot3D files
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine plot3dio_dp_read_header_ufrm(unit,file,form, &
                                        type, &
                                        mgrid,dim, &
                                        iprec,ijk,nvar, &
                                        iostat,iomsg)

! Read plot3d header info for unformatted double-precision files.

  !.. Argument Variables ..
  integer          ,intent(in)  :: unit
  character(len=*) ,intent(in)  :: file
  character(len=*) ,intent(in)  :: form
  character(len=*) ,intent(in)  :: type
  character(len=*) ,intent(in)  :: mgrid
  character(len=*) ,intent(in)  :: dim
  character(len=*) ,intent(in)  :: iprec
  integer(i8) ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer(i8) ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer          ,intent(out) :: iostat
  character(len=*) ,intent(out) :: iomsg

  !.. Local Variables ..
  integer(i8) nblk_i8
  integer nblk
  integer l,n,ndim
  integer ierr

  !.. Local Parameters ..
  character(len=*), parameter :: pname = "plot3dio_dp_read_header_ufrm"
  integer :: huge_nblk = 100000
  integer :: huge_ijk  = 100000

! -------------------------------------------------------------------------
continue

  iostat = 0  ! Assume no error until proven otherwise.
  iomsg  = ''

  if ( form /= 'FORMATTED' ) then
    read (dim,'(i1)') ndim
    if ( mgrid == 'NO' ) then
      nblk = 1
    else
      read (unit,iostat=iostat,iomsg=iomsg) nblk_i8
      nblk = int(nblk_i8)
    end if
    if ( iostat /= 0 ) then
      ! Already have iomsg from read command.
    else if ( nblk <= 0 ) then
      iostat = 1
      iomsg  = 'Number of blocks in header of Plot3d file must be positive.'
    else if ( nblk > huge_nblk ) then
      iostat = 1
      iomsg  = 'Number of blocks in header of Plot3d file is too large.'
    else
      if ( allocated(ijk) ) then
        deallocate ( ijk, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"ijk",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( ijk(1:3,1:nblk), source=1_i8, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"ijk",1,__LINE__,__FILE__,ierr,error_message)

      if ( allocated(nvar) ) then
        deallocate ( nvar, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"nvar",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( nvar(1:nblk), source=1_i8, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"nvar",1,__LINE__,__FILE__,ierr,error_message)

      if ( type == 'XYZ' ) then
        read (unit,iostat=iostat,iomsg=iomsg) &
          ((ijk(l,n),             l=1,ndim), n=1,nblk)
        nvar(:) = ndim
      else if ( type == 'SOL' ) then
        read (unit,iostat=iostat) &
          ((ijk(l,n),             l=1,ndim), n=1,nblk)
        nvar(:) = ndim+2
      else if ( type == 'FUN' ) then
        read (unit,iostat=iostat) &
          ((ijk(l,n),nvar(n), l=1,ndim), n=1,nblk)
        if ( any(nvar(:) <= 0 ) ) then
          iostat = 1
          iomsg  = 'Number of variables in header of ' // &
                   'Plot3d function file must be positive.'
        else
          iostat = 0  ! header values read ok
          iomsg  = ''
        end if
      end if  ! file type

      if ( iostat /= 0 ) then
        ! Read failed.  Error already saved in iostat and iomsg.
      else if ( any(ijk(:,:) <= 0 ) ) then
        iostat = 1
        iomsg  = 'Array bounds in header of Plot3d file must be positive.'
      else if ( any(ijk(:,:) > huge_ijk ) ) then
        iostat = 1
        iomsg  = 'Array bounds in header of Plot3d file too large.'
      else
        iostat = 0  ! header values read ok
        iomsg  = ''
      end if
    end if  ! nblk check
  else
    iostat = 1
    iomsg  = 'Routine '//pname//' cannot read '//trim(form)//' files.'
  end if  ! plot3d format

9 continue
  return

end subroutine plot3dio_dp_read_header_ufrm
!
!###############################################################################
!
subroutine plot3dio_dp_read_data_ufrm(unit,file,form, &
                                      type,blank, &
                                      mgrid,dim, &
                                      iprec,rprec, &
                                      ijk,nvar,iblank,data, &
                                      iostat,iomsg)

! Read plot3d body data info for formatted double-precision files.

  !.. Argument Variables ..
  integer                            ,intent(in)    :: unit
  character(len=*)                   ,intent(in)    :: file
  character(len=*)                   ,intent(in)    :: form
  character(len=*)                   ,intent(in)    :: type
  character(len=*)                   ,intent(in)    :: blank
  character(len=*)                   ,intent(in)    :: mgrid
  character(len=*)                   ,intent(in)    :: dim
  character(len=*)                   ,intent(in)    :: iprec
  character(len=*)                   ,intent(in)    :: rprec
  integer ,allocatable,dimension(:,:),intent(inout) :: ijk
  integer ,allocatable,dimension(:)  ,intent(inout) :: nvar
  integer ,allocatable,dimension(:)  ,intent(inout) :: iblank
  real(dp),allocatable,dimension(:,:),intent(inout) :: data
  integer                            ,intent(out)   :: iostat
  character(len=*)                   ,intent(out)   :: iomsg

  !.. Local variables ..
  integer ierr
  integer i0,i,l,iblk,ndim,nblk,nnode,npts,mxvar
  integer(i8) ,allocatable,dimension(:) :: iblank_i8
  real(dp) :: fsmach,alpha,re,time

  !.. Local Parameters ..
  character(len=*), parameter :: pname = "plot3dio_dp_read_data_ufrm"

! -------------------------------------------------------------------------
continue

  iostat = 0  ! Assume no error until proven otherwise.
  iomsg  = ''

  if ( form /= 'FORMATTED' ) then

    ndim  = size(ijk,dim=1)
    nblk  = size(ijk,dim=2)
    nnode = sum( product( ijk(1:ndim,1:nblk) , dim=1 ) )
    mxvar = maxval(nvar(:))

    if ( allocated(data) ) then
      deallocate ( data, stat=ierr, errmsg=error_message )
      call alloc_error(pname,"data",2,__LINE__,__FILE__,ierr,error_message)
    end if
    allocate ( data(1:mxvar,1:nnode), source=0.0_dp, &
              stat=ierr, errmsg=error_message )
    call alloc_error(pname,"data",1,__LINE__,__FILE__,ierr,error_message)

    if ( type == 'XYZ' .and. blank == 'YES' ) then
      if ( allocated(iblank) ) then
        deallocate ( iblank, stat=ierr, errmsg=error_message )
        call alloc_error(pname,"iblank",2,__LINE__,__FILE__,ierr,error_message)
      end if
      allocate ( iblank(1:nnode), source=1, &
                stat=ierr, errmsg=error_message )
      call alloc_error(pname,"iblank",1,__LINE__,__FILE__,ierr,error_message)

      if ( iprec == 'SINGLE' ) then
        i0 = 0
        do iblk = 1,nblk
          npts = product(ijk(1:ndim,iblk))
          read (unit,iostat=iostat,iomsg=iomsg,end=9) &
            ((data(l,i0+i), i=1,npts),l=1,nvar(iblk)), &
             (iblank(i0+i), i=1,npts)
          i0 = i0 + npts
        end do
      else
        allocate ( iblank_i8(1:nnode), source=1_i8, &
                  stat=ierr, errmsg=error_message )
        call alloc_error(pname,"iblank_i8",1,__LINE__,__FILE__, &
                         ierr,error_message)
        i0 = 0
        do iblk = 1,nblk
          npts = product(ijk(1:ndim,iblk))
          read (unit,iostat=iostat,iomsg=iomsg,end=9) &
            ((data(l,i0+i), i=1,npts),l=1,nvar(iblk)), &
             (iblank_i8(i0+i), i=1,npts)
          i0 = i0 + npts
        end do
        iblank(:) = int(iblank_i8(:))
        deallocate (iblank_i8)
      end if

    else
      i0 = 0
      do iblk = 1,nblk
        npts = product(ijk(1:ndim,iblk))
        if ( type == 'SOL' ) then
          read (unit,iostat=iostat,iomsg=iomsg,end=9) fsmach,alpha,re,time
        end if
        read (unit,iostat=iostat,iomsg=iomsg,end=9) &
          ((data(l,i0+i), i=1,npts),l=1,nvar(iblk))
        i0 = i0 + npts
      end do
    end if
  else
    iostat = 1
    iomsg  = 'Routine '//pname//' cannot read '//trim(form)//' files.'

  end if

9 continue
  return
end subroutine plot3dio_dp_read_data_ufrm
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
end module plot3dio_mod
