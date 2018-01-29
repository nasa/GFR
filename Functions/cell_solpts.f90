elemental function cell_solpts(cell_geom,cell_order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: cell_geom
  integer, intent(in) :: cell_order
  !
  !.. Function Result ..
  integer :: return_value
  !
continue
  !
  select case (cell_geom)
    case (Geom_Node)
      return_value = node_solpts(cell_order)
    case (Geom_Edge)
      return_value = edge_solpts(cell_order)
    case (Geom_Tria)
      return_value = tria_solpts(cell_order)
    case (Geom_Quad)
      return_value = quad_solpts(cell_order)
    case (Geom_Tetr)
      return_value = tetr_solpts(cell_order)
    case (Geom_Pyra)
      return_value = pyra_solpts(cell_order)
    case (Geom_Pris)
      return_value = pris_solpts(cell_order)
    case (Geom_Hexa)
      return_value = hexa_solpts(cell_order)
    case default
      return_value = 0
  end select
  !
end function cell_solpts
!
!###############################################################################
!
elemental function cell_flxpts(cell_geom,cell_order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: cell_geom
  integer, intent(in) :: cell_order
  !
  !.. Function Result ..
  integer :: return_value
  !
continue
  !
  select case (cell_geom)
    case (Geom_Node)
      return_value = node_solpts(cell_order)
    case (Geom_Edge)
      return_value = 2 * node_solpts(cell_order)
    case (Geom_Tria)
      return_value = 3 * edge_solpts(cell_order)
    case (Geom_Quad)
      return_value = 4 * edge_solpts(cell_order)
    case (Geom_Tetr)
      return_value = 4 * tria_solpts(cell_order)
    case (Geom_Pyra)
      return_value = 4 * tria_solpts(cell_order) + 1 * quad_solpts(cell_order)
    case (Geom_Pris)
      return_value = 2 * tria_solpts(cell_order) + 3 * quad_solpts(cell_order)
    case (Geom_Hexa)
      return_value = 6 * quad_solpts(cell_order)
    case default
      return_value = 0
  end select
  !
end function cell_flxpts
!
!###############################################################################
!
elemental function cell_edgpts(cell_geom,cell_order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: cell_geom
  integer, intent(in) :: cell_order
  !
  !.. Function Result ..
  integer :: return_value
  !
continue
  !
  select case (cell_geom)
    case (Geom_Tetr)
      return_value =  6 * edge_solpts(cell_order)
    case (Geom_Pyra)
      return_value =  8 * edge_solpts(cell_order)
    case (Geom_Pris)
      return_value =  9 * edge_solpts(cell_order)
    case (Geom_Hexa)
      return_value = 12 * edge_solpts(cell_order)
    case default
      return_value = 0
  end select
  !
end function cell_edgpts
!
!###############################################################################
!
elemental function node_solpts(order) result(return_value)
  integer, intent(in) :: order
  integer :: return_value
continue
  return_value = 1
end function node_solpts
!
!###############################################################################
!
elemental function edge_solpts(order) result(return_value)
  integer, intent(in) :: order
  integer :: return_value
continue
  return_value = order + 1
end function edge_solpts
!
!###############################################################################
!
elemental function tria_solpts(order) result(return_value)
  integer, intent(in) :: order
  integer :: return_value
continue
  return_value = ( (order + 1) * (order + 2) ) / 2
end function tria_solpts
!
!###############################################################################
!
elemental function quad_solpts(order) result(return_value)
  integer, intent(in) :: order
  integer :: return_value
continue
  return_value = (order + 1)**2
end function quad_solpts
!
!###############################################################################
!
elemental function tetr_solpts(order) result(return_value)
  integer, intent(in) :: order
  integer :: return_value
continue
  return_value = ( (order + 1) * (order + 2) * (order + 3) ) / 6
end function tetr_solpts
!
!###############################################################################
!
elemental function pyra_solpts(order) result(return_value)
  integer, intent(in) :: order
  integer :: return_value
continue
  return_value = ( (order + 1) * (order + 2) * (2*order + 3) ) / 6
end function pyra_solpts
!
!###############################################################################
!
elemental function pris_solpts(order) result(return_value)
  integer, intent(in) :: order
  integer :: return_value
continue
  return_value = ( (order + 1) * (order + 1) * (order + 2) ) / 2
end function pris_solpts
!
!###############################################################################
!
elemental function hexa_solpts(order) result(return_value)
  integer, intent(in) :: order
  integer :: return_value
continue
  return_value = (order + 1)**3
end function hexa_solpts
