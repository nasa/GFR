elemental function solpts2order(cell_geom,npts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: cell_geom
  integer, intent(in) :: npts
  !
  !.. Function Result ..
  integer :: return_value
  !
continue
  !
  select case (cell_geom)
    case (Geom_Node)
      return_value = solpts2order_node(npts)
    case (Geom_Edge)
      return_value = solpts2order_edge(npts)
    case (Geom_Tria)
      return_value = solpts2order_tria(npts)
    case (Geom_Quad)
      return_value = solpts2order_quad(npts)
    case (Geom_Tetr)
      return_value = solpts2order_tetr(npts)
    case (Geom_Pyra)
      return_value = solpts2order_pyra(npts)
    case (Geom_Pris)
      return_value = solpts2order_pris(npts)
    case (Geom_Hexa)
      return_value = solpts2order_hexa(npts)
    case default
      return_value = 0
  end select
  !
end function solpts2order
!
!###############################################################################
!
!elemental function solpts2flxpts(cell_geom,npts) result(return_value)
!  !
!  !.. Formal Arguments ..
!  integer, intent(in) :: cell_geom
!  integer, intent(in) :: npts
!  !
!  !.. Function Result ..
!  integer :: return_value
!  !
!continue
!  !
!  return_value = cell_flxpts( cell_geom , solpts2order(npts) )
!  !
!end function solpts2flxpts
!
!###############################################################################
!
elemental function solpts2order_node(npts) result(return_value)
  integer, intent(in) :: npts
  integer :: return_value
continue
  return_value = 1
end function solpts2order_node
!
!###############################################################################
!
elemental function solpts2order_edge(npts) result(return_value)
  integer, intent(in) :: npts
  integer :: return_value
continue
  return_value = npts - 1
end function solpts2order_edge
!
!###############################################################################
!
elemental function solpts2order_tria(npts) result(return_value)
  integer, intent(in) :: npts
  integer :: return_value
continue
  return_value = (nint(sqrt(real(9-8*(1-npts),kind=wp)))-3)/2
end function solpts2order_tria
!
!###############################################################################
!
elemental function solpts2order_quad(npts) result(return_value)
  integer, intent(in) :: npts
  integer :: return_value
continue
  return_value = nint( sqrt(real(npts,kind=wp)) ) - 1
end function solpts2order_quad
!
!###############################################################################
!
elemental function solpts2order_tetr(npts) result(return_value)
  integer, intent(in) :: npts
  integer :: return_value
  integer :: i
continue
  return_value = -1
  do i = 0,256
    if ( ((i+1)*(i+2)*(i+3))/6 == npts) then
      return_value = i
      exit
    end if
  end do
end function solpts2order_tetr
!
!###############################################################################
!
elemental function solpts2order_pyra(npts) result(return_value)
  integer, intent(in) :: npts
  integer :: return_value
  integer :: i
continue
  return_value = -1
  do i = 0,256
    if ( ((i+1)*(i+2)*(2*i+3))/6 == npts) then
      return_value = i
      exit
    end if
  end do
end function solpts2order_pyra
!
!###############################################################################
!
elemental function solpts2order_pris(npts) result(return_value)
  integer, intent(in) :: npts
  integer :: return_value
  integer :: i
continue
  return_value = -1
  do i = 0,256
    if ( ((i+1)*(i+1)*(i+2))/2 == npts) then
      return_value = i
      exit
    end if
  end do
end function solpts2order_pris
!
!###############################################################################
!
elemental function solpts2order_hexa(npts) result(return_value)
  integer, intent(in) :: npts
  integer :: return_value
continue
  return_value = nint( real(npts,kind=wp)**one3 ) - 1
end function solpts2order_hexa
