!==============================================================================!
  subroutine Multiphase_Mod_Vof_Solver_Dist_Function_Cell_Loop(flow,        &
                                                               a,           &
                                                               dist_curr,   &
                                                               r_phi,       &
                                                               grad_i,      &
                                                               grad_j,      &
                                                               grad_k)
!------------------------------------------------------------------------------!
!    Loop on cells for calculation of the Distance function. This function     !
!    only reduces the size of subroutine Compute_Distance_Function             !
!                                                                              !
!    Arguments:                                                                !
!    - dist_curr                       : distance function at current iteration!
!    - r_phi                           : norm of gradient of dist_curr         !
!    - grad_i, grad_j, grad_k          : gradient components of dist_curr      !
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: gra_x  => r_node_05,  &
                      gra_y  => r_node_06,  &
                      gra_z  => r_node_07,  &
                      dist_n => r_node_08
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: flow
  type(Matrix_Type), target :: a
  real                      :: grad_i   (-flow % pnt_grid % n_bnd_cells    &
                                        : flow % pnt_grid % n_cells),      &
                               grad_j   (-flow % pnt_grid % n_bnd_cells    &
                                        : flow % pnt_grid % n_cells),      &
                               grad_k   (-flow % pnt_grid % n_bnd_cells    &
                                        : flow % pnt_grid % n_cells),      &
                               dist_curr(-flow % pnt_grid % n_bnd_cells    &
                                        : flow % pnt_grid % n_cells),      &
                               r_phi    (-flow % pnt_grid % n_bnd_cells    &
                                        : flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  integer :: c, s, c1, c2, nb, nc, nn
  integer :: i_fac, i_nod, n
  real    :: fs
  real    :: dot_p1, dot_p2, grad_fn(3)
  real    :: g_fx, g_fy, g_fz
  real    :: a_i(8), b_i(8), c_i(8), r_nod(3)
  real    :: corr_x, corr_y, corr_z
!==============================================================================!

  grid => flow % pnt_grid

  nb = grid % n_bnd_cells
  nc = grid % n_cells
  nn = grid % n_nodes

  ! Interpolate to nodes

  call Field_Mod_Interpolate_Cells_To_Nodes(flow,                     &
                                      dist_curr(-nb:nc), dist_n(1:nn))

  ! Find gradient at nodal positions
  call Field_Mod_Grad_Component_Cells_To_Nodes(flow,                     &
                                               dist_curr, dist_n,  &
                                            1, gra_x(1:nn))
  call Field_Mod_Grad_Component_Cells_To_Nodes(flow,                     &
                                               dist_curr, dist_n,  &
                                            2, gra_y(1:nn))
  call Field_Mod_Grad_Component_Cells_To_Nodes(flow,                     &
                                               dist_curr, dist_n,  &
                                            3, gra_z(1:nn))

  do c = 1, grid % n_cells

    a_i = -HUGE; b_i = - HUGE; c_i = -HUGE

    ! Loop on nodes
    do i_nod = 1, grid % cells_n_nodes(c)
      n  = grid % cells_n(i_nod, c)
      r_nod = (/grid % xn(n) - grid % xc(c),   &
                grid % yn(n) - grid % yc(c),   &
                grid % zn(n) - grid % zc(c)/)

      if( (dist_curr(c) > 0.0 .and. r_nod(1) > 0.0) .or.     &
          (dist_curr(c) < 0.0 .and. r_nod(1) < 0.0) ) then
        a_i(i_nod) = min(0.0, gra_x(n)) ** 2
      else
        if( (dist_curr(c) > 0.0 .and. r_nod(1) < 0.0) .or.     &
            (dist_curr(c) < 0.0 .and. r_nod(1) > 0.0) ) then
          a_i(i_nod) = max(0.0, gra_x(n)) ** 2
        end if
      end if

      if( (dist_curr(c) > 0.0 .and. r_nod(2) > 0.0) .or.     &
          (dist_curr(c) < 0.0 .and. r_nod(2) < 0.0) ) then
        b_i(i_nod) = min(0.0, gra_y(n)) ** 2
      else
        if( (dist_curr(c) > 0.0 .and. r_nod(2) < 0.0) .or.     &
            (dist_curr(c) < 0.0 .and. r_nod(2) > 0.0) ) then
          b_i(i_nod) = max(0.0, gra_y(n)) ** 2
        end if
      end if

      if( (dist_curr(c) > 0.0 .and. r_nod(3) > 0.0) .or.     &
          (dist_curr(c) < 0.0 .and. r_nod(3) < 0.0) ) then
        c_i(i_nod) = min(0.0, gra_z(n)) ** 2
      else
        if( (dist_curr(c) > 0.0 .and. r_nod(3) < 0.0) .or.     &
            (dist_curr(c) < 0.0 .and. r_nod(3) > 0.0) ) then
          c_i(i_nod) = max(0.0, gra_z(n)) ** 2
        end if
      end if
    end do ! End loop nodes

    r_phi(c) = sqrt( maxval(a_i(1:grid % cells_n_nodes(c)))      &
                   + maxval(b_i(1:grid % cells_n_nodes(c)))      &
                   + maxval(c_i(1:grid % cells_n_nodes(c))) )
  end do

  end subroutine
