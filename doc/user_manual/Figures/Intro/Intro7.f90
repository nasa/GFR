!> Program to generate plots for introduction.
program main
  !
  implicit none
  !
  integer, parameter :: ncell = 5
  integer, parameter :: order = 3
  integer, parameter :: nsp   = order + 1  ! num sol  points
  integer, parameter :: nnp   = order + 3  ! num node points (sol+2face)
  !
  real,    parameter :: zero  =  0.0
  real,    parameter :: one   =  1.0
  real,    parameter :: two   =  2.0
  real,    parameter :: three =  3.0
  real,    parameter :: five  =  5.0
  real,    parameter :: six   =  6.0
  real,    parameter :: seven =  7.0
  real,    parameter :: ten   = 10.0
  !
  real, dimension(:), allocatable :: xx,yy,zz
  real, dimension(nsp,ncell) :: r_sp,x_sp,u_sp,fdisc_sp
  real, dimension(nnp,ncell) :: r_np,x_np,u_np,fdisc_np,fcomm_np
  real, dimension(nnp,ncell) :: unit_gl_np,unit_gr_np
  real, dimension(nnp,ncell) :: gl_np,gr_np,g_np,fcont_np
  real, dimension(  2,ncell) :: fcomm_fp
  real, dimension(    ncell) :: xmin,xmax,x0
  real                       :: dom_width,dom_xmin,dom_xmax
  real                       :: fdl,fdr,fluxjumpl,fluxjumpr
  real                       :: xf,xpt,upt,fpt
  integer                    :: n,nmax,nc,lc,rc,io,o
  !
continue
  !
  ! Domain of all cells (each is unit cell [-1:1])
  !
  dom_width =      ncell * two
  dom_xmin  = -dom_width / two
  dom_xmax  =  dom_width / two
  !
  ! Solution points
  !
  r_sp(:,:) = -9999.
  x_sp(:,:) = -9999.
  r_np(:,:) = -9999.
  x_np(:,:) = -9999.
  do nc = 1,ncell
!   xmin(nc) = -one - two*real(nc-1)          ! -1,-3,-5,-7 for nc=1:4
!   xmax(nc) = xmin(nc) + two                 !  1,-1,-3,-5 for nc=1:4
!   x0(nc) = -two * real(nc-1)
    xmin(nc) = dom_xmin + two*real(nc-1)
    xmax(nc) = xmin(nc) + two
    x0(nc) = ( xmin(nc) + xmax(nc) ) / two
    ! Solution points
    r_sp(:,nc) =          GaussPoints(order)  ! comput   coord, [-1:1]
    x_sp(:,nc) = x0(nc) + r_sp(:,nc)          ! physical coord
    ! Node points (solution points + boundary points)
    r_np(:,nc) = [-one, GaussPoints(order), one]
    x_np(:,nc) = x0(nc) + r_np(:,nc)
  end do
  !
  ! Solution values
  !
  u_sp(:,1) = [1.40, 1.50, 1.60, 1.60]
  u_sp(:,2) = [1.40, 1.45, 1.40, 1.30]
  u_sp(:,3) = [1.60, 1.65, 1.55, 1.65]
  u_sp(:,4) = [1.50, 1.45, 1.50, 1.60]
  u_sp(:,5) = [1.40, 1.35, 1.45, 1.40]
  !
  !
  ! Discontinuous flux values (for Burgers equation)
  !
  fdisc_sp(:,:) = u_sp(:,:)**2
  !
  ! Common flux values on interfaces
  !   Not a true common flux, simply weighted toward upstream value.
  !
  fcomm_fp(:,:) = zero
  fcomm_np(:,:) = zero
  do nc = 1, ncell
    ! left face
    lc = max(1,nc-1)      ! left  cell
    rc = nc               ! right cell
    xf = xmin(nc)  ! face location
    fpt = 0.7*Lagrange(x_sp(:,lc),fdisc_sp(:,lc),xf) &
        + 0.3*Lagrange(x_sp(:,rc),fdisc_sp(:,rc),xf)
    fcomm_fp(1,nc) = fpt
    fcomm_np(1,nc) = fpt
    ! right face
    lc = nc               ! left  cell
    rc = min(nc+1,ncell)  ! right cell
    xf = xmax(nc)  ! face location
    fpt =  0.7*Lagrange(x_sp(:,lc),fdisc_sp(:,lc),xf) &
        + 0.3*Lagrange(x_sp(:,rc),fdisc_sp(:,rc),xf)
    fcomm_fp(  2,nc) = fpt
    fcomm_np(nnp,nc) = fpt
  end do
  !
  ! Interpolate from solution points to node points.
  ! Evaluate at the node points. This will over-define the polynomials
  ! by 2 points, but they will remain the same order.
  !
  do nc = 1, ncell
    do n = 1, nnp
      xpt = r_np(n,nc)
          u_np(n,nc) = Lagrange(r_sp(:,nc),    u_sp(:,nc),xpt)
      fdisc_np(n,nc) = Lagrange(r_sp(:,nc),fdisc_sp(:,nc),xpt)
    end do
  end do
  !
  ! Correction functions (1 degree higher than solution)
  ! Evaluate at the node points. This will over-define the polynomials
  ! by 2 points, but they will remain the same order.
  !
  gl_np(:,:) = zero
  gr_np(:,:) = zero
  do nc = 1, ncell
    ! Discontinuous flux function at left/right boundaries of cell.
    fdl = Lagrange(x_sp(:,nc),fdisc_sp(:,nc),xmin(nc))
    fdr = Lagrange(x_sp(:,nc),fdisc_sp(:,nc),xmax(nc))
    ! Flux jumps at left/right boundaries of cell needed to
    ! obtain the common interface flux there.
    fluxjumpl = fcomm_fp(1,nc) - fdl
    fluxjumpr = fcomm_fp(2,nc) - fdr
    ! Left/right correction functions, and their sum.
    unit_gl_np(:,nc) = RightRadau(order+1,r_np(:,nc))
    unit_gr_np(:,nc) =  LeftRadau(order+1,r_np(:,nc))
         gl_np(:,nc) = RightRadau(order+1,r_np(:,nc)) * fluxjumpl
         gr_np(:,nc) =  LeftRadau(order+1,r_np(:,nc)) * fluxjumpr
          g_np(:,nc) = gl_np(:,nc) + gr_np(:,nc)
  end do
  !
  ! Continuous flux
  !
  do nc = 1, ncell
    fcont_np(:,nc) = fdisc_np(:,nc) + g_np(:,nc)
  end do
  !
  !
  ! Output Tecplot file for plotting.
  !
  !
  ! Output Data for Flux Reconstuction Steps.
  !
  open(newunit=io,file='Intro.tec',form='formatted',status='unknown')
  write(io,901) 'x','u','fdisc','fcomm', &
                'unit_gl','unit_gr','Gl','Gr','Gl+Gr','fcont'
  ! Output volume data.
  do nc = 1, ncell
    write(io,902) 'Cell ',nc
    do n = 1, nnp
      write(io,903) x_np(n,nc),u_np(n,nc),fdisc_np(n,nc),fcomm_np(n,nc), &
                    unit_gl_np(n,nc),unit_gr_np(n,nc), &
                    gl_np(n,nc),gr_np(n,nc),g_np(n,nc),fcont_np(n,nc)
    end do
  end do
  ! Output face data.
  if ( 0 == 1 ) then
  do nc = 1, ncell
    write(io,902) 'Face ',nc
    write(io,903) xmin(nc),zero,zero,fcomm_fp(1,nc), &
                  zero,zero,zero,zero,zero,zero
  end do
  nc = ncell
  write(io,902) 'Face ',ncell+1
  write(io,903) xmax(nc),zero,zero,fcomm_fp(2,nc), &
                zero,zero,zero,zero,zero,zero
  end if
  close(unit=io)
  !
  !
  ! Output Legendre polynomial data.
  !
  open(newunit=io,file='Legendre.tec',form='formatted',status='unknown')
    write(io,901) 'x','y'
  ! Points for plotting polynomial.
  nmax = 201
  allocate(xx(nmax))
  do n = 1, nmax
    xx(n) = -one + two * real(n-1) / real(nmax-1)
  end do
  ! Legendre Polynomial
  do nc = 0, 4 ! loop order
    write(io,902) 'Poly Leg',nc
    yy = Legendre(nc,xx)
    do n = 1, nmax
      write(io,903) xx(n),yy(n)
    end do
  end do
  ! Gauss Points
  do nc = 1, 4 ! loop order
    write(io,902) 'Zeros Leg',nc
    xx = GaussPoints(nc-1)
    yy = Legendre(nc,xx)
    do n = 1, size(xx)
      write(io,903) xx(n),yy(n)
    end do
  end do
  close(unit=io)
  deallocate(xx)
  !
  !
  ! Output Radau polynomial data.
  !
  open(newunit=io,file='Radau.tec',form='formatted',status='unknown')
  write(io,901) 'x','y'
  ! Points for plotting polynomial.
  nmax = 201
  allocate(xx(nmax))
  do n = 1, nmax
    xx(n) = -one + two * real(n-1) / real(nmax-1)
  end do
  ! Left Radau Polynomial
  do nc = 1, 4  ! loop order
    write(io,902) 'Poly R_L',nc
    yy = LeftRadau(nc,xx)
    do n = 1, nmax
      write(io,903) xx(n),yy(n)
    end do
  end do
  ! Left Radau Points
  do nc = 1, 4  ! loop order
    write(io,902) 'Points R_L',nc
    zz = [       LeftRadauZeros(nc), one]
    yy = [zero * LeftRadauZeros(nc), one]
    do n = 1, size(zz)
      write(io,903) zz(n),yy(n)
    end do
  end do
  ! Right Radau Polynomial
  do nc = 1, 4  ! loop order
    write(io,902) 'Poly R_R',nc
    yy = RightRadau(nc,xx)
    do n = 1, nmax
      write(io,903) xx(n),yy(n)
    end do
  end do
  ! Right Radau Points
  do nc = 1, 4  ! loop order
    write(io,902) 'Points R_R',nc
    zz = [-one,        RightRadauZeros(nc)]
    yy = [ one, zero * RightRadauZeros(nc)]
    do n = 1, size(zz)
      write(io,903) zz(n),yy(n)
    end do
  end do
  close(unit=io)
  deallocate(xx)
  !
901 format('Variables=',99('"',a,'"',:,','))
902 format('Zone T="',a,i0,'" F=POINT')
903 format(99es15.5)
  !
  call OutputEqualPoints(   '1D_Points_Equal.tec'  )
  call OutputGaussPoints(   '1D_Points_Gauss.tec'  )
  call OutputLobattoPoints( '1D_Points_Lobatto.tec')
  !
contains
!==========================================================================
!==========================================================================
function Lagrange(xpt,ypt,x) result(return_value)
  !
  real, dimension(:), intent(in) :: xpt,ypt
  real,               intent(in) :: x
  real                           :: return_value
  !
  integer :: i,n
  real    :: val
  real,   dimension(size(xpt)) :: unit_lagrange
  !
continue
  !
  do n = 1, size(xpt)
    ! Compute the unit Lagrange function l_i(x) at this x.
    val = one
    do i = 1, size(xpt)
      if ( i /= n ) then
        val = val * ( x - xpt(i) ) / ( xpt(n) - xpt(i) )
      end if
    end do
    unit_lagrange(n) = val
    !
  end do
  !
  ! Compute the Lagrange polynomial at this x.
  ! Sum( u_i*l_i(x) )
  return_value = sum( ypt(:) * unit_lagrange(:) )
  !
end function Lagrange
!==========================================================================
!> Compute the Legendre polynomial values at a given set of r locations.
!> The r locations should be between [-1:1].
!>
!> Recursive formula:
!>   P[-1] = 0
!>   P[ 0] = 1
!>   P[ k] = (2k-1)/k * r * P[k-1] - (k-1)/k * P[k-2]
!> which yields:
!>   P[ 1] = r
!>   P[ 2] = (3*x^2-1)/2
!>   P[ 3] = (5*r^3-3*r) / 2
!>   P[ 4] = (35*r^4-30*r^2+3) / 8
!==========================================================================
recursive function Legendre(order,r) result(return_value)
  !
  integer,               intent(in) :: order
  real,    dimension(:), intent(in) :: r
  real,    dimension(size(r))       :: return_value
  !
continue
  !
  if     ( order == -1 ) then
    return_value(:) = 0
  elseif ( order ==  0 ) then
    return_value(:) = 1
  elseif ( order >   0 ) then
    return_value(:) = &
      real(2*order-1) / real(order) * r(:) * Legendre(order-1,r(:)) &
    - real(  order-1) / real(order) *        Legendre(order-2,r(:))
  else ! invalid order
    write(*,*) 'Fun Legendre: Invalid order, ',order
    return_value(:) = 0
  end if
  !
end function Legendre
!==========================================================================
!> Compute the derivative of the Legendre polynomial at a given set of r
!> locations. The r locations should be between [-1:1].
!>
!> Recursive formula:
!>   P'[-1] = 0
!>   P'[ 0] = 0
!>   P'[ k] = (2k-1)/k * P[k-1] + (2k-1)/k * r * P'[k-1] - (k-1)/k * P'[k-2]
!> which yields:
!>   P'[ 1] = 1
!>   P'[ 2] = 3*r
!>   P'[ 3] = 15/2*r^2 - 3/2
!>   P'[ 4] = 35/2*r^3 - 15/2*r
!==========================================================================
recursive function LegendreDeriv(order,r) result(return_value)
  !
  integer,               intent(in) :: order
  real,    dimension(:), intent(in) :: r
  real,    dimension(size(r))       :: return_value
  !
continue
  !
  if     ( order == -1 ) then
    return_value(:) = 0
  elseif ( order ==  0 ) then
    return_value(:) = 0
  elseif ( order >   0 ) then
    return_value(:) = &
      real(2*order-1) / real(order) *        Legendre(     order-1,r(:)) &
    + real(2*order-1) / real(order) * r(:) * LegendreDeriv(order-1,r(:)) &
    - real(  order-1) / real(order) *        LegendreDeriv(order-2,r(:))
  else ! invalid order
    write(*,*) 'Fun LegendreDeriv: Invalid order, ',order
    return_value(:) = 0
  end if
  !
end function LegendreDeriv
!==========================================================================
!> Compute the left Radau polynomial values at a given set of r locations.
!> The r locations should be between [-1:1].
!==========================================================================
function LeftRadau(order,r) result(return_value)
  !
  integer,               intent(in) :: order
  real,    dimension(:), intent(in) :: r
  real,    dimension(size(r))       :: return_value
  !
  return_value(:) = ( Legendre(order,r(:)) + Legendre(order-1,r(:)) ) / two
  !
end function LeftRadau
!==========================================================================
!> Compute the right Radau polynomial values at a given set of r locations.
!> The r locations should be between [-1:1].
!==========================================================================
function RightRadau(order,r) result(return_value)
  !
  integer,               intent(in) :: order
  real,    dimension(:), intent(in) :: r
  real,    dimension(size(r))       :: return_value
  !
  return_value(:) = ( Legendre(order,r(:)) - Legendre(order-1,r(:)) ) / two &
                  * (-one)**order
  !
end function RightRadau
!==========================================================================
!> Compute zeros of the left Radau polynomial function.
!==========================================================================
function LeftRadauZeros(order) result(return_value)
  !
  integer,               intent(in) :: order
  real,    dimension(max(1,order))  :: return_value
  !
  real :: r1,r2,r3,r4
  !
  select case (order)
    case (1)
      r1 = -one
      return_value = [r1]
    case (2)
      r1 = -one
      r2 =  one / three
      return_value = [r1,r2]
    case (3)
      r1 = -one
      r2 = (one - sqrt(six)) / five
      r3 = (one + sqrt(six)) / five
      return_value = [r1,r2,r3]
    case (4)
      r1 = -one
      r2 = -0.575318923521694
      r3 =  0.181066271118531
      r4 =  0.822824080974592
      return_value = [r1,r2,r3,r4]
    case default
      write(*,*) 'Fun LeftRadauZeros: Invalid order, ',order
      return_value(:) = zero
  end select
end function LeftRadauZeros
!==========================================================================
!> Compute zeros of the right Radau polynomial function.
!==========================================================================
function RightRadauZeros(order) result(return_value)
  !
  integer,               intent(in) :: order
  real,    dimension(max(1,order))  :: return_value
  !
  real :: r1,r2,r3,r4
  !
  select case (order)
    case (1)
      r1 =  one
      return_value = [r1]
    case (2)
      r1 = -one / three
      r2 =  one
      return_value = [r1,r2]
    case (3)
      r1 = -(one + sqrt(six)) / five
      r2 = -(one - sqrt(six)) / five
      r3 =  one
      return_value = [r1,r2,r3]
    case (4)
      r1 = -0.822824080974592
      r2 = -0.181066271118531
      r3 =  0.575318923521694
      r4 =  one
      return_value = [r1,r2,r3,r4]
    case default
      write(*,*) 'Fun RightRadauZeros: Invalid order, ',order
      return_value(:) = zero
  end select
end function RightRadauZeros
!==========================================================================
!==========================================================================
function EqualPoints(order) result(return_value)
  !
  integer, intent(in) :: order
  real,    dimension(max(1,order+1)) :: return_value
  !
  real, parameter :: rmin = -1., rmax = 1.
  integer ::n,nmax
  !
  nmax = order+1
  do n = 1, nmax
    return_value(n) = rmin + (rmax-rmin) * real(n-1) / real(nmax-1)
  end do
  !
end function EqualPoints
!==========================================================================
!> Compute the Gauss points in the interval r=[-1:1].
!>
!> These are the zeros of the Legendre polynomial of specified order.
!> The number of Gauss points is the same as the order.
!==========================================================================
function GaussPoints(order) result(return_value)
  !
  integer, intent(in)              :: order
  real,    dimension(max(1,order+1)) :: return_value
  !
  real :: r1,r2
  !
continue
  !
  ! Exact values taken from:
  ! https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
  select case (order)
    case (0) ; return_value = [0]
    case (1) ; r1 = one/sqrt(three)
               return_value = [-r1,r1]
    case (2) ; r1 = sqrt(three/five)
               return_value = [-r1,zero,r1]
    case (3) ; r1 = sqrt( (three - two*sqrt(six/five)) / seven )
               r2 = sqrt( (three + two*sqrt(six/five)) / seven )
               return_value = [-r2,-r1,r1,r2]
    case (4) ; r1 = sqrt( five - two*sqrt(ten/seven) ) / three
               r2 = sqrt( five + two*sqrt(ten/seven) ) / three
               return_value = [-r2,-r1,zero,r1,r2]
    case default
      write(*,*) 'Fun GaussPoints: Invalid order, ',order
      return_value(:) = zero
  end select
  !
end function GaussPoints
!==========================================================================
!> Compute the Lobatto points in the interval r=[-1:1].
!> These are the points where left Radau polynomials RL_k & RL_k-1 cross.
!>
!> Since Radau polynomials are based on Legendre polynomials,
!> these are the points where Legendre polynomials P_k & P_k-2 cross.
!==========================================================================
function LobattoPoints(order) result(return_value)
  !
  integer, intent(in)              :: order
  real,    dimension(max(1,order+1)) :: return_value
  !
  real :: r1,r2
  !
continue
  !
  ! Exact values computed by hand from Legendre polynomials.
  select case (order)
    case (0) ; return_value = [0]
    case (1) ; r1 = one
               return_value = [-r1,r1]
    case (2) ; r1 = one
               return_value = [-r1,zero,r1]
    case (3) ; r1 = two * sqrt( two / 35. )
               r2 = one
               return_value = [-r2,-r1,r1,r2]
    case (4) ; r1 = sqrt( 3. / 7. )
               r2 = one
               return_value = [-r2,-r1,zero,r1,r2]
    case default
      write(*,*) 'Fun LobattoPoints: Invalid order, ',order
      return_value(:) = zero
  end select
  !
end function LobattoPoints
!==========================================================================
!> Output 1D equidistant points to Tecplot file.
!==========================================================================
subroutine OutputEqualPoints(filename)
  ! 
  character(len=*), intent(in) :: filename
  ! 
  integer :: order,n,nmax
  real, allocatable, dimension(:) :: r,y
  ! 
  open(newunit=io,file=trim(filename),form='formatted',status='unknown')
    write(io,901) 'r','y'
    do order = 1, 4
      write(io,902) 'Equal Points, Order=',order
      nmax = order+1
      r = EqualPoints(order)
      allocate(y(nmax),source=real(order))
      do n = 1, nmax
        write(io,903) r(n),y(n)
      end do
      deallocate(r)
      deallocate(y)
    end do
  close(unit=io)
  ! 
901 format('Variables=',99('"',a,'"',:,','))
902 format('Zone T="',a,i0,'" F=POINT')
903 format(99es15.5)
end subroutine OutputEqualPoints
!==========================================================================
!> Output 1D Gauss points to Tecplot file.
!==========================================================================
subroutine OutputGaussPoints(filename)
  ! 
  character(len=*), intent(in) :: filename
  !
  integer :: order,n,nmax
  real, allocatable, dimension(:) :: r,y
  real, parameter :: rmin = -1., rmax = 1.
  ! 
  open(newunit=io,file=trim(filename),form='formatted',status='unknown')
    write(io,901) 'r','y'
    do order = 1, 4
      write(io,902) 'Gauss Points, Order=',order
      nmax = order+1
      r = GaussPoints(order)
      allocate(y(nmax),source=real(order))
      do n = 1, nmax
        write(io,903) r(n),y(n)
      end do
      deallocate(r)
      deallocate(y)
    end do
  close(unit=io)
  ! 
901 format('Variables=',99('"',a,'"',:,','))
902 format('Zone T="',a,i0,'" F=POINT')
903 format(99es15.5)
end subroutine OutputGaussPoints
!==========================================================================
!> Output 1D Lobatto points to Tecplot file.
!==========================================================================
subroutine OutputLobattoPoints(filename)
  !
  character(len=*), intent(in) :: filename
  !
  integer :: order,n,nmax
  real, allocatable, dimension(:) :: r,y
  real, parameter :: rmin = -1., rmax = 1.
  ! 
  open(newunit=io,file=trim(filename),form='formatted',status='unknown')
    write(io,901) 'r','y'
    do order = 1, 4
      write(io,902) 'Equal Points, Order=',order
      nmax = order+1
      r = LobattoPoints(order)
      allocate(y(nmax),source=real(order))
      do n = 1, nmax
        write(io,903) r(n),y(n)
      end do
      deallocate(r)
      deallocate(y)
    end do
  close(unit=io)
  ! 
901 format('Variables=',99('"',a,'"',:,','))
902 format('Zone T="',a,i0,'" F=POINT')
903 format(99es15.5)
end subroutine OutputLobattoPoints
!==========================================================================
!
end program main
