!==============================================================================!
  subroutine Field_Mod_Interpolate_Cells_To_Nodes(flow, var_cell, var_node)
!------------------------------------------------------------------------------!
!   Interpolates a variable from cell centers to nodes                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)  :: flow
  real, intent(in)  :: var_cell( -flow % pnt_grid % n_bnd_cells  &
                                 :flow % pnt_grid % n_cells)
  real, intent(out) :: var_node(1:flow % pnt_grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: n, s, c, c1, c2, i_cell, i_fac
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! It might need a call to: call Grid_Mod_Exchange_Cells_Real(grid, var_cell)

  ! Interpolate from all cells to nodes
  do n = 1, grid % n_nodes

    var_node(n) = 0.0

    ! Loop on cells
    do i_cell = 1, grid % nodes_n_cells(n)

      c = grid % nodes_c(i_cell, n)
      var_node(n) = var_node(n)  &
                  + grid % nodes_weight_c(i_cell, n) * var_cell(c)

    end do

  end do

  call Grid_Mod_Exchange_Nodes_Real(grid, var_node)

  end subroutine
