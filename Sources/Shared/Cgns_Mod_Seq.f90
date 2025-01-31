!=============================================================================!
  module Cgns_Mod
!-----------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
  use Sort_Mod
!-----------------------------------------------------------------------------!
  implicit none
!-----------------------------------------------------------------------------!
  include "cgns_io_f.h"
  include "cgnslib_f.h"
!=============================================================================!

  !----------!
  !   File   ! -> contains base
  !----------!
    !
    !----------!
    !   Base   ! -> contains blocks
    !----------!
      !
      !------------!
      !   Blocks   ! -> contains coordinates, sections, bnd_conds, solution
      !------------!
        !
        !-----------------!
        !   Coordinates   !
        !-----------------!
        !
        !----------------------!
        !   Element sections   !
        !----------------------!
        !
        !----------------------!
        !   Solution section   !
        !----------------------!
          !
          !------------------!
          !   Field section  !
          !------------------!

  ! File
  integer       :: file_id
  character(SL) :: file_name
  integer       :: file_mode
  logical       :: verbose = .false.

  ! Solution section
  type Cgns_Solution_Type
    character(SL) :: name
    integer       :: sol_type
  end type

  ! Element section
  type Cgns_Section_Type
    character(SL) :: name
    integer       :: cell_type
    integer       :: first_cell
    integer       :: last_cell
  end type

  ! Blocks
  type Cgns_Block_Type
    character(SL)                         :: name
    integer                               :: type
    integer                               :: mesh_info(3)
    integer                               :: n_sects
    type(Cgns_Section_Type), allocatable  :: section(:)
    integer                               :: n_coords
    character(SL)                         :: coord_name(3)
    integer                               :: n_solutions
    type(Cgns_Solution_Type), allocatable :: solution(:)
  end type

  ! Base
  integer :: n_bases
  type Cgns_Base_Type
    character(SL)                      :: name
    integer                            :: cell_dim
    integer                            :: phys_dim
    integer                            :: n_blocks
    type(Cgns_Block_Type), allocatable :: block(:)
  end type
  type(Cgns_Base_Type), allocatable :: cgns_base(:)

  integer :: cnt_cells  ! number of cells (except boundary and interface)

  ! Cells (3d)
  integer :: cnt_hex
  integer :: cnt_pyr
  integer :: cnt_wed
  integer :: cnt_tet

  ! If actual grid was written, further saves have just a link to that grid
  logical       :: mesh_written = .false.
  character(SL) :: file_with_mesh

  contains

  include 'Cgns_Mod/Initialize_Counters.f90'

  include 'Cgns_Mod/Write_Link_To_Mesh_In_File.f90'
  include 'Cgns_Mod/Write_Link_To_Field.f90'
  include 'Cgns_Mod/Write_Dimensions_Info.f90'

  ! Seq and Par
  include 'Cgns_Mod/Sequential/Open_File.f90'
  include 'Cgns_Mod/Sequential/Close_File.f90'
  include 'Cgns_Mod/Sequential/Write_Base_Info.f90'
  include 'Cgns_Mod/Sequential/Write_Block_Info.f90'
  include 'Cgns_Mod/Sequential/Write_Coordinate_Array.f90'
  include 'Cgns_Mod/Sequential/Write_Section_Connections.f90'
  include 'Cgns_Mod/Sequential/Write_Solution_Info.f90'
  include 'Cgns_Mod/Sequential/Write_Field.f90'

  end module
