pure function bc_string_to_integer(bc_text) result(return_value)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: bc_text
  !
  !.. Function Return Value ..
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: n,nb,ne
  character(len=len_trim(bc_text)) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: cbnd = "BOUNDARY"
  !
continue
  !
  ! Copy the input BC string into text so we can manipulate it
  !
  text = uppercase( trim(adjustl(bc_text)) )
  !
  nb = 1
  ne = len_trim(text)
  !
  ! Remove all spaces from the character string
  !
  do
    n = index(text(nb:ne)," ")
    if (n == 0) exit
    text = text(nb:n-1) // text(n+1:ne)
    ne = len_trim(text)
  end do
  !
  ! Remove the substring "BOUNDARY" if it exists
  !
  n = index(text(nb:ne),cbnd)
  if (n > 0) then
    text = text(nb:n-1) // text(n+len_trim(cbnd):ne)
    ne = len_trim(text)
  end if
  !
  ! Use the character string given in the PhysicalName
  ! section to set the value of the boundary condition\
  !
  select case (trim(adjustl(text)))
    !
    case ("MESHINTERIOR")
      return_value = not_a_bc
      !
    case ("INFLOW","INLET")
      return_value = bc_generic_inflow
    case ("SUBINFLOW")
      return_value = bc_sub_inflow
    case ("SUPINFLOW")
      return_value = bc_sup_inflow
      !
    case ("OUTFLOW","EXIT")
      return_value = bc_generic_outflow
    case ("SUBOUTFLOW")
      return_value = bc_sub_outflow
    case ("SUPOUTFLOW","EXTRAP")
      return_value = bc_sup_outflow
      !
    case ("FREEFLOW")
      return_value = bc_generic_freeflow
    case ("CHARACTERISTIC","CHAR","CHARACTER")
      return_value = bc_characteristic
    case ("FREESTREAM","FARFIELD","REFERENCE")
      return_value = bc_freestream
    case ("FIXED")
      return_value = bc_fixed
      !
    case ("WALL")
      return_value = bc_default_wall
    case ("SWALL","SLIP","SLIPWALL")
      return_value = bc_slip_wall
    case ("EWALL","EULER","EULERWALL")
      return_value = bc_euler_wall
    case ("AWALL","ADIABATIC")
      return_value = bc_adiabatic_wall
    case ("IWALL","ISOTHERMAL")
      return_value = bc_isothermal_wall
      !
    case ("SYMM","SYMMETRY")
      return_value = bc_symmetry
    case ("PERIODIC")
      return_value = bc_periodic
      !
    case default
      return_value = bc_unknown
      !
  end select
  !
end function bc_string_to_integer
!
!###############################################################################
!
