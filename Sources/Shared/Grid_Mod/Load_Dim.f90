!==============================================================================!
  subroutine Load_Dim(Grid, this_proc, domain)
!------------------------------------------------------------------------------!
!   Reads file with grid dimensions (.dim, used to be .geo)                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
  integer, intent(in) :: this_proc
  integer, optional   :: domain
!-----------------------------------[Locals]-----------------------------------!
  integer       :: c, n, s, fu
  character(SL) :: name_in
!==============================================================================!

  !----------------------------!
  !                            !
  !     Read the file with     !
  !   geometrical dimensions   !
  !                            !
  !----------------------------!
  call File % Set_Name(name_in, processor=this_proc, extension='.dim',  &
                       domain=domain)
  call File % Open_For_Reading_Binary(name_in, fu, this_proc)

  read(fu) (Grid % xn(n), n = 1, Grid % n_nodes)
  read(fu) (Grid % yn(n), n = 1, Grid % n_nodes)
  read(fu) (Grid % zn(n), n = 1, Grid % n_nodes)

  read(fu) (Grid % xc(c), c = -Grid % n_bnd_cells, Grid % n_cells)
  read(fu) (Grid % yc(c), c = -Grid % n_bnd_cells, Grid % n_cells)
  read(fu) (Grid % zc(c), c = -Grid % n_bnd_cells, Grid % n_cells)

  read(fu) (Grid % wall_dist(c), c = -Grid % n_bnd_cells, Grid % n_cells)
  read(fu) (Grid % vol(c), c = 1, Grid % n_cells)

  read(fu) (Grid % sx(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % sy(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % sz(s), s = 1, Grid % n_faces + Grid % n_shadows)

  read(fu) (Grid % dx(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % dy(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % dz(s), s = 1, Grid % n_faces + Grid % n_shadows)

  read(fu) (Grid % f(s), s = 1, Grid % n_faces + Grid % n_shadows)

  read(fu) (Grid % xf(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % yf(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % zf(s), s = 1, Grid % n_faces + Grid % n_shadows)

  read(fu) (Grid % rx(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % ry(s), s = 1, Grid % n_faces + Grid % n_shadows)
  read(fu) (Grid % rz(s), s = 1, Grid % n_faces + Grid % n_shadows)

  read(fu) Grid % per_x
  read(fu) Grid % per_y
  read(fu) Grid % per_z

  close(fu)

  end subroutine
