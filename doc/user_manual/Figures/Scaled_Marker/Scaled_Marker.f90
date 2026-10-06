!> Compute data for scaled marker plot.
program main

  implicit none

  ! Number of Delta Psi intervals
  real,    parameter :: pi        =  acos(-1.)
  real,    parameter :: zeta_min  = -pi/2. * 3.
  real,    parameter :: zeta_max  =  pi/2. * 3.
  integer, parameter :: nmax      = 1 + 100* 3

  integer :: io,n
  real    :: zeta

continue

  open(newunit=io,file='Scaled_Marker.tec',form='formatted',status='unknown')
    write(io,1) 'Variables="zeta","zeta/(pi/2)","S"'
    write(io,1) 'Zone T="Scaled Marker" F=POINT'
    do n = 1, nmax
      zeta = zeta_min + (zeta_max-zeta_min) * real(n-1) / real(nmax-1)
      write(io,2) zeta, zeta/(pi/2.), scaled_marker(zeta)
    enddo
  close(io)
  1 format(a)
  2 format(2x,3f10.5)

contains
!==========================================================================
function scaled_marker(zeta) result(return_value)

  real, intent(in) :: zeta
  real             :: return_value

continue

  if ( zeta .le. -pi/2. ) then
    return_value = 0
  elseif ( abs(zeta) .le. pi/2. ) then
    return_value = ( 1. + sin(zeta) ) / 2.
  else
    return_value = 1
  endif

end function scaled_marker
!==========================================================================
end program main
