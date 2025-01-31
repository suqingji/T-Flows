!==============================================================================!
  subroutine Allocate_Front(Front, Flow)
!------------------------------------------------------------------------------!
!   Surface genesis                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  type(Field_Type),  target :: Flow
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: nb, nc, nf
!==============================================================================!

  ! Take aliases to object vertex Flow around
  Front % pnt_flow => Flow
  Front % pnt_grid => Flow % pnt_grid

  ! Front divides the surface mesh among processors
  Front % mesh_divided = .true.

  ! Allocate memory
  allocate(Front % Elem(MAX_SURFACE_ELEMENTS))
  allocate(Front % Vert(MAX_SURFACE_VERTICES))
  allocate(Front % side(MAX_SURFACE_ELEMENTS * 3))

  if(Flow % mass_transfer) then
    Grid => Flow % pnt_grid
    nb   =  Grid % n_bnd_cells
    nc   =  Grid % n_cells
    nf   =  Grid % n_faces

    allocate(Front % cell_at_elem(-nb:nc)); Front % cell_at_elem(-nb:nc) = 0
    allocate(Front % face_at_elem(  2,nf)); Front % face_at_elem(  :, :) = 0
  end if

  ! Initialize front's local variables
  call Front % Initialize_Front()

  end subroutine
