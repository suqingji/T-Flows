!==============================================================================!
  subroutine Form_Maps(Grid)
!------------------------------------------------------------------------------!
!   Forms maps for parallel backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, c1g, c2g, s, cnt
!==============================================================================!
!   There is an issue with this procedure, but it's more related to MPI/IO     !
!   functions than T-Flows.  In cases a subdomain has no physical boundary     !
!   cells, variable "nb_sub" turns out to be zero.  This, per se, should not   !
!   be an issue if MPI/IO functions could handle calls to:                     !
!   "Mpi_Type_Create_Indexed_Block(...)" and "Mpi_File_Write(...)" with zero   !
!   length.  But they don't.  Therefore, I avoid allocation with zero size     !
!   (max(nb_sub,1)) here and creation of new types with zero size in           !
!   "Comm_Mod_Create_New_Types".  It is a bit of a dirty trick :-(             !
!------------------------------------------------------------------------------!

  ! Initialize number of cells in subdomain
  Grid % Comm % nc_sub = Grid % n_cells - Grid % Comm % n_buff_cells

  ! Initialize number of boundary cells in subdomain
  Grid % Comm % nb_sub = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % Comm % cell_proc(c) .eq. this_proc) then
      Grid % Comm % nb_sub = Grid % Comm % nb_sub + 1
    end if
  end do

  ! Initialize number of faces in subdomain
  Grid % Comm % nf_sub = 0
  Grid % Comm % nf_tot = 0
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    c1g = Grid % Comm % cell_glo(c1)
    c2g = Grid % Comm % cell_glo(c2)
    ! Count boundary faces as well as inside whose c1g < c2g for total count
    if(c2g < 0 .or. c1g < c2g) then
      Grid % Comm % nf_tot = Grid % Comm % nf_tot + 1
    end if
    Grid % Comm % nf_sub = Grid % Comm % nf_sub + 1
  end do

  ! First and last cell to send
  Grid % Comm % nb_f = Grid % n_bnd_cells
  Grid % Comm % nb_l = Grid % n_bnd_cells - Grid % Comm % nb_sub + 1

  ! Initialize total number of cells
  Grid % Comm % nc_tot = Grid % Comm % nc_sub
  Grid % Comm % nb_tot = Grid % Comm % nb_sub

  !--------------------------------!
  !                                !
  !   For run with one processor   !
  !                                !
  !--------------------------------!
  if(n_proc < 2) then

    !-------------------------------------!
    !   Global cell numbers for T-Flows   !
    !-------------------------------------!
    do c = -Grid % Comm % nb_tot, Grid % Comm % nc_tot
      Grid % Comm % cell_glo(c) = c
    end do

    !-----------------------------!
    !   Create mapping matrices   !
    !-----------------------------!
    allocate(Grid % Comm % cell_map    (Grid % Comm % nc_sub))
    allocate(Grid % Comm % bnd_cell_map(Grid % Comm % nb_sub))
    allocate(Grid % Comm % face_map    (Grid % Comm % nf_sub))

    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, Grid % Comm % nc_tot
      Grid % Comm % cell_map(c) = int(c - 1, SP)
    end do

    ! -1 is to start from zero, as needed by MPI functions
    do c = 1, Grid % Comm % nb_tot
      Grid % Comm % bnd_cell_map(c) = int(c - 1, SP)
    end do

    ! -1 is to start from zero, as needed by MPI functions
    do s = 1, Grid % Comm % nf_tot
      Grid % Comm % face_map(s) = int(s - 1, SP)
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)
      c1g = Grid % Comm % cell_glo(c1)
      c2g = Grid % Comm % cell_glo(c2)
    end do

  !-----------------------!
  !                       !
  !   For parallel runs   !
  !                       !
  !-----------------------!
  else

    call Comm_Mod_Global_Sum_Int(Grid % Comm % nc_tot)
    call Comm_Mod_Global_Sum_Int(Grid % Comm % nb_tot)
    call Comm_Mod_Global_Sum_Int(Grid % Comm % nf_tot)

    !-----------------------------!
    !   Create mapping matrices   !
    !-----------------------------!
    allocate(Grid % Comm % cell_map        (Grid % Comm % nc_sub))
    allocate(Grid % Comm % bnd_cell_map(max(Grid % Comm % nb_sub,1)))
    allocate(Grid % Comm % face_map        (Grid % Comm % nf_sub))
    Grid % Comm % cell_map(:)     = 0
    Grid % Comm % bnd_cell_map(:) = 0
    Grid % Comm % face_map(:)     = 0

    !---------------------!
    !   Inside cell map   !
    !---------------------!
    do c = 1, Grid % Comm % nc_sub
      ! Take cell mapping to be the same as global cell numbers but start from 0
      Grid % Comm % cell_map(c) = int(Grid % Comm % cell_glo(c) - 1, SP)
    end do

    !-----------------------!
    !   Boundary cell map   !
    !-----------------------!
    cnt = 0
    do c = -Grid % Comm % nb_f, -Grid % Comm % nb_l
      cnt = cnt + 1
      Grid % Comm % bnd_cell_map(cnt) = int(  Grid % Comm % cell_glo(c)  &
                                            + Grid % Comm % nb_tot, SP)
    end do

    !--------------!
    !   Face map   !
    !--------------!
    cnt = 0
    do s = 1, Grid % n_faces
      cnt = cnt + 1
      Grid % Comm % face_map(cnt) = int(Grid % Comm % face_glo(s) - 1, SP)
    end do
  end if

  end subroutine
