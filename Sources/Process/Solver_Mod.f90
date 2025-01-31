!==============================================================================!
  module Solver_Mod
!------------------------------------------------------------------------------!
!   Module used for native linear solvers.                                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Matrix_Mod
  use Vector_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Solvers type   !
  !------------------!
  type Solver_Type

    type(Grid_Type), pointer :: pnt_grid

    ! Matrix for all variables except momentum
    type(Matrix_Type) :: A

    ! Matrix for discretized momentum equations
    type(Matrix_Type) :: M

    ! Preconditioning matrix for all variables (used inside solvers only)
    type(Matrix_Type) :: D

    ! Right-hand side for all variables
    ! (used in solvers and during discretization)
    type(Vector_Type) :: b

    contains
      procedure          :: Alias_Solver
      procedure          :: Bicg                 ! bicg solver
      procedure          :: Cg                   ! cg solver
      procedure          :: Cgs                  ! cgs solver
      procedure          :: Create_Solver
      procedure, private :: Normalized_Root_Mean_Square
      procedure, private :: Prec_Form
      procedure, private :: Prec_Solve
      procedure, private :: Residual_Vector
      procedure, private :: Root_Mean_Square

  end type

  contains

  include 'Solver_Mod/Alias_Solver.f90'
  include 'Solver_Mod/Bicg.f90'
  include 'Solver_Mod/Cg.f90'
  include 'Solver_Mod/Cgs.f90'
  include 'Solver_Mod/Create_Solver.f90'
  include 'Solver_Mod/Normalized_Root_Mean_Square.f90'
  include 'Solver_Mod/Prec_Form.f90'
  include 'Solver_Mod/Prec_Solve.f90'
  include 'Solver_Mod/Residual_Vector.f90'
  include 'Solver_Mod/Root_Mean_Square.f90'

  end module 
