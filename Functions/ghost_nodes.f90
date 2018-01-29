pure function ghost_nodes(nodes,geom,side) result(return_value)
  !
  integer, intent(in) :: geom
  integer, intent(in) :: side
  integer, dimension(:), intent(in) :: nodes
  !
  integer, dimension(1:geom_nodes(geom)) :: return_value
  !
continue
  !
  ! Set the untouched nodes to a weird negative value
  !
  return_value = -99
  !
  if (any(geom == Geom_2D)) then
    !
    ! The same switch applies for both triangles and quadrilaterals
    !
    return_value(1) = nodes(2)
    return_value(2) = nodes(1)
    !
  else if (geom == Geom_Pyra) then
    !
    ! Special node ordering for pyramids depending on which side of the host
    ! cell is on the boundary face
    !
    if (side == 5) then
      return_value(1:size(nodes)) = nodes
    else
      return_value(1) = nodes(2)
      return_value(2) = nodes(1)
      return_value(5) = nodes(3)
    end if
    !
  else if (geom == Geom_Pris) then
    !
    ! Special node ordering for prisms depending on which side of the host
    ! cell is on the boundary face
    !
    if (any(side == [1,2])) then
      return_value(1:size(nodes)) = nodes
    else
      return_value(1) = nodes(1)
      return_value(3) = nodes(2)
      return_value(6) = nodes(3)
      return_value(4) = nodes(4)
    end if
    !
  else
    !
    ! For all other geometries (Geom_Edge,Geom_Tetr,Geom_Hexa), just copy
    ! the nodes over to return_value keeping them in the same order.
    !
    return_value(1:size(nodes)) = nodes
    !
  end if
  !
end function ghost_nodes
