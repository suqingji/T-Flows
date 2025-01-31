!==============================================================================!
  subroutine Solve_System(Vof, Sol, b)
!------------------------------------------------------------------------------!
!   Solves linear system for VOF                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),   target :: Vof
  type(Solver_Type), target :: Sol
  real, contiguous,  target :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),    pointer :: fun
  type(Field_Type),  pointer :: flow
  type(Matrix_Type), pointer :: A
  character(SL)              :: solver
!==============================================================================!

  ! Take aliases
  fun  => Vof % fun
  flow => Vof % pnt_flow
  A    => Sol % A

  ! Get solver
  call Control_Mod_Solver_For_Vof(solver)

  call Cpu_Timer % Start('Linear_Solver_For_Vof')
  call Sol % Bicg(A,              &
                  fun % n,        &
                  b,              &
                  fun % precond,  &
                  fun % mniter,   &      ! max number of iterations
                  fun % eniter,   &      ! executed number of iterations
                  fun % tol,      &
                  fun % res)
  call Cpu_Timer % Stop('Linear_Solver_For_Vof')

  if(.not. flow % heat_transfer) then
    call Info_Mod_Iter_Fill_At(1, 6, fun % name, fun % eniter, fun % res)
  else
    call Info_Mod_Iter_Fill_At(2, 1, fun % name, fun % eniter, fun % res)
  end if

  end subroutine
