!==============================================================================!
  module Vector_Mod
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !-----------------!
  !   Vector type   !
  !-----------------!
  type Vector_Type

    type(Grid_Type), pointer :: pnt_grid

    real, allocatable :: val(:)    ! value
  end type

  contains

  include 'Vector_Mod/Allocate.f90'

  end module 
