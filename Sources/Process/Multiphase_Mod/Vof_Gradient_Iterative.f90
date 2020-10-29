!==============================================================================!
  subroutine Multiphase_Mod_Vof_Gradient_Iterative(mult, var, g_x, g_y, g_z)
!------------------------------------------------------------------------------!
!   Computes de gradient of scalr var. using Gauss theorem in a iterative way  !
!   to take into account skewness and to produce values at faces that can be   !
!   reconstructed using gradients at cells. This technique can be found at     !
!   https://spiral.imperial.ac.uk/handle/10044/1/28101                         !
!                                                                              !
!   Arguments:                                                                 !
!   - g_x, g_y, g_z                 : Gradient components of var               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: aux_x  => r_cell_27,  &
                      aux_y  => r_cell_28,  &
                      aux_z  => r_cell_29,  &
                      tol_x  => r_cell_24,  &
                      tol_y  => r_cell_25,  &
                      tol_z  => r_cell_26
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: g_x(-mult % pnt_grid % n_bnd_cells :   &
                                        mult % pnt_grid % n_cells),       &
                                   g_y(-mult % pnt_grid % n_bnd_cells :   &
                                        mult % pnt_grid % n_cells),       &
                                   g_z(-mult % pnt_grid % n_bnd_cells :   &
                                        mult % pnt_grid % n_cells),       &
                                   var(-mult % pnt_grid % n_bnd_cells :   &
                                        mult % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),     pointer :: flow
  type(Grid_Type),      pointer :: grid
  integer                       :: s, c, c1, c2, n, i, n_iter, nb, nc, bc
  real                          :: varf, toler_grad, epsloc, fs
!==============================================================================!

  flow => mult % pnt_flow
  grid => mult % pnt_grid

  nb = grid % n_bnd_cells
  nc = grid % n_cells
  bc = grid % comm % n_buff_cells

  epsloc = epsilon(epsloc)

  ! Calculate gradients using Least Squares as a guess
  call Field_Mod_Grad_Component(flow, var(-nb:nc), 1, g_x(-nb:nc))
  call Field_Mod_Grad_Component(flow, var(-nb:nc), 2, g_y(-nb:nc))
  call Field_Mod_Grad_Component(flow, var(-nb:nc), 3, g_z(-nb:nc))

  toler_grad = HUGE
  n_iter = 100
  i = 0

  do while (i < n_iter .and. toler_grad > 1.0e-10)

    i = i + 1

    aux_x(-nb:nc) = 0.0
    aux_y(-nb:nc) = 0.0
    aux_z(-nb:nc) = 0.0

    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      ! Add new velocity component correction
      varf = var(c1)

      varf = varf + dot_product( (/g_x(c1), g_y(c1), g_z(c1)/),            &
                    (/grid % xr(s), grid % yr(s), grid % zr(s)/) )

      aux_x(c1) = aux_x(c1) + varf * grid % sx(s) / grid % vol(c1)
      aux_y(c1) = aux_y(c1) + varf * grid % sy(s) / grid % vol(c1)
      aux_z(c1) = aux_z(c1) + varf * grid % sz(s) / grid % vol(c1)
    end do

    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      ! Add new velocity component correction
      varf = fs * var(c1) + (1.0 - fs) * var(c2)

      varf = varf +         fs * dot_product( (/g_x(c1), g_y(c1), g_z(c1)/),            &
                                 (/grid % xr(s), grid % yr(s), grid % zr(s)/) )

      varf = varf + (1.0 - fs) * dot_product( (/g_x(c2), g_y(c2), g_z(c2)/),            &
                                 (/grid % xr(s), grid % yr(s), grid % zr(s)/) )

      aux_x(c1) = aux_x(c1) + varf * grid % sx(s) / grid % vol(c1)
      aux_y(c1) = aux_y(c1) + varf * grid % sy(s) / grid % vol(c1)
      aux_z(c1) = aux_z(c1) + varf * grid % sz(s) / grid % vol(c1)

      aux_x(c2) = aux_x(c2) - varf * grid % sx(s) / grid % vol(c2)
      aux_y(c2) = aux_y(c2) - varf * grid % sy(s) / grid % vol(c2)
      aux_z(c2) = aux_z(c2) - varf * grid % sz(s) / grid % vol(c2)
    end do

    call Grid_Mod_Exchange_Cells_Real(grid, aux_x(-nb:nc))
    call Grid_Mod_Exchange_Cells_Real(grid, aux_y(-nb:nc))
    call Grid_Mod_Exchange_Cells_Real(grid, aux_z(-nb:nc))

    do c = 1, grid % n_cells
      tol_x(c) = abs(g_x(c) - aux_x(c))
      tol_y(c) = abs(g_y(c) - aux_y(c))
      tol_z(c) = abs(g_z(c) - aux_z(c))
    end do

    g_x(-nb:nc) = aux_x(-nb:nc)
    g_y(-nb:nc) = aux_y(-nb:nc)
    g_z(-nb:nc) = aux_z(-nb:nc)

    toler_grad = max(maxval(tol_x(1:nc-bc)),   &
                     maxval(tol_y(1:nc-bc)),   &
                     maxval(tol_z(1:nc-bc)))

    call Comm_Mod_Global_Max_Real(toler_grad)

  end do

  end subroutine
