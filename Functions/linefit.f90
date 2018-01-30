program linefit

 implicit none

 real(8), allocatable :: x(:)
 real(8), allocatable :: y(:)
 real(8), allocatable :: s(:)

 integer :: i,n
 real(8) :: a,b,ea,eb,x2

 open(1,file='xy.dat',status='old')
 n=0
 do
   read(1,*,end=1)a
   n=n+1
 end do
 1 rewind(1)

 allocate(x(n))
 allocate(y(n))
 allocate(s(n))

 do i=1,n
   read(1,*)x(i),y(i),s(i)
 end do
 close(1)
 call lfit(n,x,y,s,a,b,ea,eb,x2)

 print*,'Intercept and expected error : ',a,ea
 print*,'Slope and expected error     : ',b,eb
 print*,'Number of data points        : ',n
 print*,'X2 per degree of freedom     : ',x2

contains
!--------------------------------------------------------------------!
!Fits a straight line through n points with x-values in x(mdat),     !
!y-values in y(mdat) and statistical error in s(mdat). The returned  !
!line-parameters are a and b (y=a+b*x), and their statistical errors !
!(one standard deviation) are ea and eb. The chi-squared value (per  !
!degrees of freedom = n-2) is returned as x2.                        !
!--------------------------------------------------------------------!
subroutine lfit(n,x,y,s,a,b,ea,eb,x2)
!-------------------------------------!

 integer :: i,n
 real(8) ::  a,b,ea,eb,da,db,w,wx,wy,wxx,wxy,x2,wgt
 real(8) ::  x(n),y(n),s(n)

 w=0.d0
 wx=0.d0
 wy=0.d0
 wxx=0.d0
 wxy=0.d0
 do i=1,n
    wgt=1.d0/s(i)**2
    w=w+wgt
    wx=wx+wgt*x(i)
    wy=wy+wgt*y(i)
    wxx=wxx+wgt*x(i)*x(i)
    wxy=wxy+wgt*x(i)*y(i)
 end do
 b=(wxy*w-wx*wy)/(wxx*w-wx**2)
 a=(wy-b*wx)/w
 x2=0.d0
 do i=1,n
    x2=x2+(a+b*x(i)-y(i))**2/(s(i)**2)
 end do
 x2=x2/dble(n-2)

 ea=0.d0
 eb=0.d0
 do i=1,n
    db=(x(i)*w-wx)/((wxx*w-wx**2)*s(i)**2)
    da=(1.d0/s(i)**2-db*wx)/w
    ea=ea+(da*s(i))**2
    eb=eb+(db*s(i))**2
 end do
 ea=sqrt(ea)
 eb=sqrt(eb)

end subroutine lfit
!-------------------!

end program linefit
