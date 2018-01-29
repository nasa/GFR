elemental function bc_integer_to_string(ibc) result(return_value)
  !
  integer, intent(in) :: ibc
  !
  character(len=28) :: return_value
  !
continue
  !
  select case (ibc)
    case (bc_generic_inflow)
      return_value = "generic_inflow"
    case (bc_sub_inflow)
      return_value = "sub_inflow"
    case (bc_sup_inflow)
      return_value = "sup_inflow"
    case (bc_mdot_inflow)
      return_value = "mdot_inflow"
    case (bc_generic_outflow)
      return_value = "generic_outflow"
    case (bc_sub_outflow)
      return_value = "sub_outflow"
    case (bc_sup_outflow)
      return_value = "sup_outflow"
    case (bc_mdot_outflow)
      return_value = "mdot_outflow"
    case (bc_generic_freeflow)
      return_value = "generic_freeflow"
    case (bc_characteristic)
      return_value = "characteristic"
    case (bc_freestream)
      return_value = "freestream"
    case (bc_fixed)
      return_value = "fixed"
    case (bc_slip_wall)
      return_value = "slip_wall"
    case (bc_euler_wall)
      return_value = "euler_wall"
    case (bc_adiabatic_wall)
      return_value = "adiabatic_wall"
    case (bc_isothermal_wall)
      return_value = "isothermal_wall"
    case (not_a_bc)
      return_value = "not_a_bc"
    case (bc_symmetry)
      return_value = "symmetry"
    case (bc_periodic)
      return_value = "periodic"
    case (bc_mms_dirichlet)
      return_value = "MMS_dirichlet"
    case (bc_cpu_bnd)
      return_value = "cpu_bnd"
    case (bc_unknown)
      return_value = "unknown"
    case default
      return_value = "BC_INTEGER_TO_STRING ERROR!!"
  end select
  !
  return_value = uppercase(return_value)
  !
end function bc_integer_to_string
