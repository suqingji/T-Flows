!==============================================================================!
  real function Tau_Wall_Log_Law(turb, dens, u_tau, u_tan, wall_dist, & 
                                 y_plus, z_o)
!------------------------------------------------------------------------------!
!   Calculates wall shear stress by using log law.                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turb_Mod,  only: Turb_Type, kappa, e_log
  use Const_Mod, only: TINY
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type) :: turb
  real            :: dens, u_tau, u_tan, wall_dist, y_plus, z_o
!==============================================================================!

  if(z_o > tiny) then

    Tau_Wall_Log_Law = dens * kappa * u_tau * u_tan  &
                       / log(((wall_dist + z_o) / z_o))
  else

    Tau_Wall_Log_Law = dens * kappa * u_tau * u_tan   &   
                       / log(e_log * max(y_plus, 1.05))
  end if

  end function
