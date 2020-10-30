!==============================================================================!
  subroutine Multiphase_Mod_Vof_Heaviside_Function(mult)
!------------------------------------------------------------------------------!
!   Computes the Heaviside function, necessary to thicken curvature when using !
!   the distance function                                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vof
  type(Var_Type),   pointer :: dist_func
  integer                   :: c, s, c1, c2, nb, nc
  real                      :: eps_grid
!==============================================================================!

  flow      => mult % pnt_flow
  grid      => mult % pnt_grid
  vof       => mult % vof
  dist_func => mult % dist_func

  nb = grid % n_bnd_cells
  nc = grid % n_cells


  do c = 1, grid % n_cells

    eps_grid = mult % heaviside_mollified_factor * grid % vol(c) ** ONE_THIRD

    if (dist_func % n(c) > eps_grid) then
      mult % heaviside_func(c) = 1.0
    else if (dist_func % n(c) < -eps_grid) then
      mult % heaviside_func(c) = 0.0
    else
      mult % heaviside_func(c) = 0.5 * ( 1.0 + dist_func % n(c) / eps_grid     &
                                       + 1.0 / PI * sin(PI * dist_func % n(c)  &
                                       / eps_grid))
    end if

  end do

  ! Values at boundaries
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (c2 < 0) then
      mult % heaviside_func(c2) = mult % heaviside_func(c1)
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, mult % heaviside_func(-nb:nc))

  end subroutine
