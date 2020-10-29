!==============================================================================!
  subroutine Grid_Mod_Find_Cells_Cells(grid)
!------------------------------------------------------------------------------!
!   Find cells surrounding each cell; needed for VOF for now.                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s   ! counters
!==============================================================================!

  grid % cells_n_cells = 0

  !---------------------------------------!
  !   Browse through all faces to form:   !
  !        cells_n_faces, cells_f         !
  !---------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    grid % cells_n_cells(c1) = grid % cells_n_cells(c1) + 1
    grid % cells_c(grid % cells_n_cells(c1), c1) = c2
    if (c2 > 0) then
      grid % cells_n_cells(c2) = grid % cells_n_cells(c2) + 1
      grid % cells_c(grid % cells_n_cells(c2), c2) = c1
    end if
  end do

  end subroutine
