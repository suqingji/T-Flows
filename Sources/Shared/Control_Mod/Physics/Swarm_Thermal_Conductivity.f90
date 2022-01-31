!==============================================================================!
  subroutine Control_Mod_Swarm_Thermal_Conductivity(therm_cond, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: therm_cond
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def 
!==============================================================================!

  def = 1.38

  call Control_Mod_Read_Real_Item('SWARM_THERMAL_CONDUCTIVITY',  &
                                   def, therm_cond, verbose)

  end subroutine
