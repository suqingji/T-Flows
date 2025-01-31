!==============================================================================!
  module Results_Mod
!------------------------------------------------------------------------------!
!   Module containig functions for saving numerical results for visualization. !
!   It comes in two flavors: "Vtu" and "Cgns", depending on the file format    !
!   one wants to save.  It has (and uses) a sister module "Save_Grid_Mod",     !
!   which is in the directory "Shared".                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Backup_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Results type   !
  !------------------!
  type Results_Type

    integer :: interval  ! result save interval
    logical :: initial   ! save intial condition or not
    logical :: boundary  ! save results at boundaries or not

    contains
      procedure :: Main_Results

      procedure, private :: Save_Results
      procedure, private :: Save_Scalar_Int
      procedure, private :: Save_Scalar_Real
      procedure, private :: Save_Front
      procedure, private :: Save_Surf
      procedure, private :: Save_Swarm
      procedure, private :: Save_Vector_Real
      procedure, private :: Time_To_Save

  end type

  type(Results_Type) :: Results

  contains

  include 'Results_Mod/Main_Results.f90'
  include 'Results_Mod/Save_Results.f90'      ! binary
  include 'Results_Mod/Save_Front.f90'
  include 'Results_Mod/Save_Scalar_Int.f90'   ! binary
  include 'Results_Mod/Save_Scalar_Real.f90'  ! binary
  include 'Results_Mod/Save_Swarm.f90'
  include 'Results_Mod/Save_Vector_Real.f90'  ! binary
  include 'Results_Mod/Save_Surf.f90'
  include 'Results_Mod/Time_To_Save.f90'

  end module 
