!==============================================================================!
  subroutine Grad_Least_Pressure_Correction(Flow, pp)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure correction.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
  type(Var_Type)    :: pp
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: s, c1, c2
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  ! Refresh buffers for variable
  call Grid % Exchange_Cells_Real(pp % n)

  ! No correction at boundaries
  ! Tried to extrapolate to boundaries in previous revision, but didn't
  ! work for jet.  It blew in the first time step.  Maybe extrapolations
  ! should be done at later stages of simulation, I am not sure
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid % Bnd_Cond_Type(c2) .ne. PRESSURE) then
        pp % n(c2) = 0.0
      end if
    end if
  end do

  ! Compute individual gradients without refreshing buffers
  call Flow % Grad_Component_No_Refresh(pp % n, 1, pp % x)  ! dp/dx
  call Flow % Grad_Component_No_Refresh(pp % n, 2, pp % y)  ! dp/dy
  call Flow % Grad_Component_No_Refresh(pp % n, 3, pp % z)  ! dp/dz

  ! Refresh buffers for gradient components
  call Grid % Exchange_Cells_Real(pp % x)
  call Grid % Exchange_Cells_Real(pp % y)
  call Grid % Exchange_Cells_Real(pp % z)

  end subroutine
