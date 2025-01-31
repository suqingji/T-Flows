!==============================================================================!
  subroutine Numerics_Mod_Min_Max(phi)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c1, c2, s
!==============================================================================!

  ! Take alias to Grid
  Grid => phi % pnt_grid

  phi % min(:) = phi % n(:)
  phi % max(:) = phi % n(:)

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    phi % min(c1) = min(phi % min(c1), phi % n(c2))
    phi % min(c2) = min(phi % min(c2), phi % n(c1))

    phi % max(c1) = max(phi % max(c1), phi % n(c2))
    phi % max(c2) = max(phi % max(c2), phi % n(c1))
  end do

  call Grid % Exchange_Cells_Real(phi % min)
  call Grid % Exchange_Cells_Real(phi % max)

  end subroutine
