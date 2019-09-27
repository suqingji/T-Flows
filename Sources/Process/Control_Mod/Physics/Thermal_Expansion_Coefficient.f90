!==============================================================================!
  subroutine Control_Mod_Thermal_Expansion_Coefficient(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('THERMAL_EXPANSION_COEFFICIENT', 1.0,  &
                                   val, verbose)

  end subroutine
