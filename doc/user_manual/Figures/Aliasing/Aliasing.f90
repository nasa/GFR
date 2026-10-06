!==========================================================================
!> Generate data to illustrate aliasing error.
!==========================================================================
program main

  implicit none

  integer, parameter :: nmax = 201
  real   , parameter :: xmin = -1.0
  real   , parameter :: xmax =  1.0
  integer :: io,n
  real    :: x
  real    :: xpts(5) = [ -0.906180, &
                         -0.538469, &
                          0.      , &
                          0.538469, &
                          0.906180  ] ! Gauss points

continue

  open(newunit=io,file='Aliasing.tec',form='formatted',status='unknown')
    write(io,1) 'Variables="x","u"'

    ! Exact solution (u)
    write(io,2) 'Exact Solution, u'
    do n = 1,nmax
      x = xmin + (xmax-xmin) * real(n-1) / real(nmax-1)
      write(io,3) x,u(x)
    enddo

    ! Solution at Gauss points (u)
    write(io,2) 'Polynomial Solution, u'
    do n = 1,size(xpts)
      write(io,3) xpts(n),u(xpts(n))
    enddo

    ! Exact flux (u^2)
    write(io,2) 'Exact Flux, u^2'
    do n = 1,nmax
      x = xmin + (xmax-xmin) * real(n-1) / real(nmax-1)
      write(io,3) x,u(x)**2
    enddo

    ! Flux at Gauss points (u^2)
    write(io,2) 'Polynomial Flux, u^2'
    do n = 1,size(xpts)
      write(io,3) xpts(n),u(xpts(n))**2
    enddo

  close(io)

  1 format(a)
  2 format('Zone T="',a,'" F=POINT')
  3 format(2x,f10.5,f10.5)

contains

!==========================================================================
! Define the solution as a 4th order polynomial.
!   u(x) = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4
!==========================================================================
function u(x) result(return_value)
  real             :: return_value
  real, intent(in) :: x

  ! Use Tecplot to get coefficients of 4th order polynomial fit of points:
  !   (x1,y1) = (-1.0, 1.0)
  !   (x2,y2) = (-0.5,-0.5)
  !   (x3,y3) = ( 0.0, 0.3)
  !   (x4,y4) = ( 0.5,-0.5)
  !   (x5,y5) = ( 1.0, 1.0)
  real, parameter :: c0 =  0.3
  real, parameter :: c1 =  0.
  real, parameter :: c2 = -4.5
  real, parameter :: c3 =  0.
  real, parameter :: c4 =  5.2

continue
  return_value = c0 + c1*x + c2*x**2 + c3*x**3 + c4*x**4
end function u
!==========================================================================
end program main
