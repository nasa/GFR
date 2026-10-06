program main

  implicit none

  integer, parameter :: npts = 101
  real    ::  etacs(5) = [.0,.2,.4,.6,.8]
  real    :: orders(5) = [2,4,8,16,32]
  real    :: alphas(5) = [2,4,8,16,32]
  real    :: etac,order,alpha
  integer :: io,i
  character(len=80) :: zonename

continue
  open(newunit=io,file='Filter.tec',form='formatted',status='unknown')
    write(io,1) 'Variables="eta","filter"'
    ! Sweep etac
    do i = 1, size(etacs)
      etac = etacs(i)
      write(zonename,2) '<greek>h</greek><sub>c</sub>',etac
      call write_filter(zonename,etac=etac,alpha=32.0,order=8.0)
    end do
    ! Sweep alpha
    do i = 1, size(alphas)
      alpha = alphas(i)
      write(zonename,3) '<greek>a</greek>',int(alpha)
      call write_filter(zonename,etac=0.6,alpha=alpha,order=8.0)
    end do
    ! Sweep order
    do i = 1, size(orders)
      order = orders(i)
      write(zonename,3) 's',int(order)
      call write_filter(zonename,etac=0.6,alpha=32.0,order=order)
    end do
  close(io)
  1 format(a)
  2 format(a,' = ',f4.1)
  3 format(a,' = ',i0)
contains
!==========================================================================
!> Exponential filter.
function filter(eta,etac,alpha,order) result(return_value)
  real             :: return_value
  real, intent(in) :: eta    ! modal fraction [0:1]
  real, intent(in) :: etac   ! modal fraction cutoff [0:1]
  real, intent(in) :: alpha  ! filter strength [0:32]
  real, intent(in) :: order  ! filter order, typically [0,2,4,8,16,32]
continue
  return_value = 1.
  if ( eta > etac ) then
    return_value = exp( -alpha * ((eta-etac)/(1.-etac))**order )
  end if
end function filter
!==========================================================================
subroutine write_filter(zonename,etac,alpha,order)
  character(len=*) :: zonename  ! Tecplot zone name
  real, intent(in) :: etac      ! modal fraction cutoff [0:1]
  real, intent(in) :: alpha     ! filter strength [0:32]
  real, intent(in) :: order     ! filter order, typically [0,2,4,8,16,32]
  real    :: etamin = 0.
  real    :: etamax = 1.
  real    :: eta
  integer :: n
continue
  write(io,1) trim(zonename)
  write(io,2) 'etac' ,etac
  write(io,3) 'alpha',int(alpha)
  write(io,3) 'order',int(order)
  do n = 1, npts
    eta = etamin + (etamax-etamin) * real(n-1) / real(npts-1)
    write(io,5) eta,filter(eta,etac,alpha,order)
  end do
  1 format('Zone T="',a,'" F=POINT')
  2 format('AuxData ',a,'="',f3.1,'"')
  3 format('AuxData ',a,'="',i0,'"')
  5 format(2x,f10.8,e20.10)
end subroutine write_filter
!==========================================================================
end program main
