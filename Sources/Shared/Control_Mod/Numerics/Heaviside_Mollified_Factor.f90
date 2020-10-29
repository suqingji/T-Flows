!==============================================================================!
  subroutine Control_Mod_Heaviside_Mollified_Factor(val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real              :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Real_Item('HEAVISIDE_MOLLIFIED_FACTOR',  &
                                   1.5, val, verbose)

  end subroutine
