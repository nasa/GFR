!> Program to compute the modal and nodal decomposition of a polynomial.
program poly_form
  !
  implicit none
  !
  character(len=80)  :: filename = 'Poly_Form.tec'
  integer, parameter :: P = 3
  integer, parameter :: NP = P+1
  integer, parameter :: imax = 201
  real,    parameter :: xmin = -1.0
  real,    parameter :: xmax =  1.0
  !
  character(len=20), dimension(NP) :: modal_name, nodal_name
  integer :: io,i,m,n
  real, dimension(NP) :: u_modal, u_nodal, x_nodal, lagrange
  real, dimension(NP) :: modal_term, nodal_term
  real                :: x,px
  !
continue
  !
  ! Define the solution polynomial in modal form.
  u_modal(:) = [1.0, 2.0, 3.0, 4.0]
  !
  ! Define the nodal points.
  ! Evenly spaced for now. Should use Gauss points instead.
! do n = 1,NP
!   x_nodal(n) = xmin + (xmax-xmin) * real(n-1) / real(NP-1)
! end do
  ! Get NP=P+1 Gauss points for polynomial degree P.
  x_nodal(:) = get_gauss_points(P)
  !
  ! Compute the nodal solution values.
  do n = 1,NP
    x = x_nodal(n)
    u_nodal(n) = 0.
    do m = 1,NP
      u_nodal(n) = u_nodal(n) + u_modal(m) * psi(m,x)
    end do ! m
  end do ! n
  !
  ! Names of the modal and nodal output terms.
  do n = 1,NP
    write(modal_name(n),101) n,n
    write(nodal_name(n),102) n,n
  end do
101 format('uhat_',i0,'*psi_',i0)
102 format(   'u_',i0,  '*l_',i0)
  !
  ! Output the result.
  open(newunit=io,file=trim(filename),form='formatted',status='unknown')
    write(io,201) 'r', &
                  'zero', &
                  'p(r)', &
                  ( trim(modal_name(n)), n=1,NP ), &
                  ( trim(nodal_name(n)), n=1,NP )
    write(io,202) 'P3 Solution'
!   do i = 1,imax
!     x = xmin + (xmax-xmin) * real(i-1) / real(imax-1)
    do i = 1,NP
      x = x_nodal(i)
      !
      ! Compute the solution polynomial based on the modal form.
      px = 0.
      do n = 1,NP
        px = px + u_modal(n) * psi(n,x)
      end do ! n
      !
      ! Compute the Lagrange polynomial at  this x:  l_n(x)
      do n = 1, NP
        lagrange(n) = 1.
        do m = 1,NP
          if ( m .ne. n ) then
            lagrange(n) = lagrange(n) * ( x          - x_nodal(m) ) &
                                      / ( x_nodal(n) - x_nodal(m) )
          end if
        end do ! m
      end do ! n
      !
      ! Compute the modal and nodal decomposition at this x.
      do n = 1, NP
        modal_term(n) = u_modal(n) * psi(n,x)
        nodal_term(n) = u_nodal(n) * lagrange(n)
      end do ! n
      !
      ! Output the array values.
      write(io,203) x, &
                    0., &
                    px, &
                    ( modal_term(n), n=1,NP ), &
                    ( nodal_term(n), n=1,NP )
    end do ! i
    !
    ! Output a zone to show y=0
    !
    if ( .false. ) then
    write(io,202) 'y=0'
    do i = 1,2 
      x = xmin + (xmax-xmin) * real(i-1) / real(2-1)
      write(io,203) x, &
                    0., &
                    ( 0.*n, n=1,NP ), &
                    ( 0.*n, n=1,NP )
      write(io,203) x, &
                    0., &
                    ( 0.*n, n=1,NP ), &
                    ( 0.*n, n=1,NP )
    end do ! i
  close(unit=io)
201 format('Variables=',99('"',a,'"',:,','))
202 format('Zone T="',a,'" F=POINT')
203 format(99f15.10)
    end if
!
contains
!
function psi(n,x) result(return_value)
  !
  implicit none
  !
  integer, intent(in) :: n
  real,    intent(in) :: x
  real                :: return_value
  !
  return_value = x**(n-1)
  !
end function psi
!
function get_gauss_points(n) result(return_value)
  !
  implicit none
  !
  integer, intent(in) :: n                         !< Polynomial order
  real, allocatable, dimension(:) :: return_value  !< n+1 gauss points
  !
continue
  !
  select case (n)
    case (1) ; return_value = [ -0.577350, &
                                 0.577350  ]
    case (2) ; return_value = [ -0.774597, &
                                 0.      , &
                                 0.774597  ]
    case (3) ; return_value = [ -0.861136, &
                                -0.339981, &
                                 0.339981, &
                                 0.861136  ]
    case (4) ; return_value = [ -0.906180, &
                                -0.538469, &
                                 0.      , &
                                 0.538469, &
                                 0.906180  ]
    case default
      write(*,*) 'ERROR: Gauss points not programmed for p=',n
      stop
  end select
end function get_gauss_points
!
end program poly_form
