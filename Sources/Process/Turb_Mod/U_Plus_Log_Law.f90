!==============================================================================!
  real function U_Plus_Log_Law(turb, wall_dist, y_plus, z_o)
!------------------------------------------------------------------------------!
!   Calculates U+ by using log law. 
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: TINY
  use Turb_Mod,  only: Turb_Type, kappa, e_log
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type) :: turb
  real            :: y_plus, wall_dist, z_o
!==============================================================================!

  if(z_o > tiny) then 

    U_Plus_Log_Law = log( (wall_dist + z_o) / z_o)  &
                      / (kappa + TINY) + TINY
  else

    U_Plus_Log_Law = log( max(y_plus, 1.05) * e_log ) / kappa
  end if

  end function
