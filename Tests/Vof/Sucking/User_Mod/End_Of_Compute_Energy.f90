!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Energy(Flow, turb, Vof, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Energy function.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt  ! current time step
  integer, intent(in)         :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: s, fu
  logical, save            :: first_entry = .true.
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  if(first_entry) then
    call File % Open_For_Writing_Ascii('sucking_solution.dat', fu)
  else
    call File % Append_For_Writing_Ascii('sucking_solution.dat', fu)
  end if

  do s = 1, Grid % n_faces

    if(any(Vof % Front % face_at_elem(1:2,s) .ne. 0)) then

      ! Write down Stefan's solution
      if(ini .eq. 1                            .and.  &
         Math % Approx_Real(grid % ys(s), 0.0) .and.  &
         Math % Approx_Real(grid % zs(s), 0.0)) then
        write(fu,  '(99(es12.4))') curr_dt * Flow % dt, Grid % xs(s)
      end if

    end if

  end do

  first_entry = .false.

  close(fu)

  end subroutine
