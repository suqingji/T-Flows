!==============================================================================!
  subroutine Load_Physical_Properties(grid)
!------------------------------------------------------------------------------!
!   Reads physical properties from control file.                               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Field_Mod
  use Control_Mod
  use Rans_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  call Control_Mod_Dynamic_Viscosity    (viscosity)
  call Control_Mod_Heat_Capacity        (capacity)
  call Control_Mod_Mass_Density         (density)
  call Control_Mod_Thermal_Conductivity (conductivity)
<<<<<<< HEAD
  call Control_Mod_Thermal_Expansion_Coefficient (beta_tec)
=======
  call Control_Mod_Species_Diffusivity  (diffusivity)
>>>>>>> f1d9682d818b16c7eaf3998e4932fc283d2d2a4a
  call Control_Mod_Roughness_Coefficient(z_o)

  end subroutine
