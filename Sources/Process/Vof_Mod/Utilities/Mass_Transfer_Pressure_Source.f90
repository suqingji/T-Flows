!==============================================================================!
  subroutine Mass_Transfer_Pressure_Source(Vof, b)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  real                    :: b(Vof % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  integer                   :: c, e
!==============================================================================!

  ! Take aliases
  Grid => Vof % pnt_grid
  Flow => Vof % pnt_flow

  if(.not. Flow % mass_transfer) return

  ! Add volume over all cells, avoiding buffer cells
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    e = Vof % Front % cell_at_elem(c)  ! front element
    if(e .ne. 0) then
      b(c) = b(c)            &
           + Vof % m_dot(c) * Vof % Front % elem(e) % area
    end if
  end do

  end subroutine
