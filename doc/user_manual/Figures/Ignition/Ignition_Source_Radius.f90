!==========================================================================
!> Generate Tecplot files with data for the decay of
!> the chemistry ignition source with distance.
!==========================================================================
program ignition_source_radius
  implicit none
  integer, parameter :: nmax = 201
  integer :: n, io
  real    :: r           ! (r/R)
  real    :: rmin = 0.0  ! (r/R)min
  real    :: rmax = 2.0  ! (r/R)min
  real    :: source

continue

  open(newunit=io,file='Ignition_Source_Radius.tec',form='formatted',&
       status='unknown')
    write(io,1) 'Variables="r/R","Source"'
    !
    write(io,2) 'Ignition_Source_Radius < 0'
    do n = 1, nmax
      r = rmin + (rmax-rmin) * real(n-1) / real(nmax-1)
      source = exp(-36.*r**2)
      write(io,3) r,source
    end do
    !
    write(io,2) 'Ignition_Source_Radius > 0'
    do n = 1, nmax
      r = rmin + (rmax-rmin) * real(n-1) / real(nmax-1)
      source = 1.0
      if ( r > 1.0 ) source = exp(-36.*(r-1.0)**3)
      write(io,3) r,source
    end do
    !
  1 format(a)
  2 format('Zone T="',a,'" F=POINT')
  3 format(2x,f8.4,e14.5)
    !
  close(unit=io)

end program ignition_source_radius
