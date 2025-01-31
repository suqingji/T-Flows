!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, turb, Vof, swarm,  &
                                       n, n_stat_t, n_stat_p, time)
!------------------------------------------------------------------------------!
!   This function is computing benchmark for rising bubble.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: swarm
  integer, intent(in)      :: n         ! current time step
  integer, intent(in)      :: n_stat_t  ! 1st t.s. statistics turbulence
  integer, intent(in)      :: n_stat_p  ! 1st t.s. statistics particles
  real,    intent(in)      :: time      ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  integer                  :: c, last_cell, fu
  real                     :: b_volume, surface, rise_velocity
  real                     :: circularity, c_position
  real, save :: w_in, c_position_old
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  if(n .eq. 1) c_position_old = 0.0

  !-------------------!
  !   Bubble volume   !
  !-------------------!
  b_volume = 0.0
  surface = 0.0
  c_position = 0.0
  rise_velocity = 0.0
  call Flow % Grad_Variable(fun)

  do c = 1, Grid % n_cells - Grid % comm % n_buff_cells
    b_volume = b_volume + Grid % vol(c) * fun % n(c)
    if (norm2((/fun % x(c),fun % y(c),fun % z(c)/)) > 1.0) then

      surface = surface + sqrt(fun % x(c) ** 2                    &
                             + fun % y(c) ** 2                    &
                             + fun % z(c) ** 2) * Grid % vol(c)
    end if
    c_position = c_position + Grid % zc(c) * fun % n(c) * Grid % vol(c)
    rise_velocity = rise_velocity + Flow % w % n(c) * fun % n(c) * Grid % vol(c)
  end do

  call Comm_Mod_Global_Sum_Real(b_volume)
  call Comm_Mod_Global_Sum_Real(surface)
  call Comm_Mod_Global_Sum_Real(c_position)
  call Comm_Mod_Global_Sum_Real(rise_velocity)
  c_position   = c_position    / b_volume
  rise_velocity = rise_velocity / b_volume
  w_in = (c_position - c_position_old) / Flow % dt

  IF(N > 1) THEN
    IF(THIS_PROC < 2) THEN
      PRINT *, 'C_POSITION        = ', C_POSITION
      PRINT *, 'C_POSITION_OLD    = ', C_POSITION_OLD
      PRINT *, 'RISE_VELOCITY (1) = ', RISE_VELOCITY
      PRINT *, 'RISE_VELOCITY (2) = ', W_IN
    END IF
  END IF
  c_position_old = c_position


  IF(N > 1) THEN
    if(w_in > 0.0 .and. c_position > 0.0) then
      do c = -Grid % n_bnd_cells, -1
        if (Grid % Bnd_Cond_Name( c) .eq. 'IN') then
          Flow % w % n(c) = -w_in
          Flow % w % b(c) = Flow % w % n(c)
        end if
      end do
      IF(THIS_PROC < 2) PRINT *, 'INCREASING THE PUMP; NEW W_IN = ', W_IN
    end if
  END IF

  ! Write to file
  if (this_proc < 2) then
    call File % Append_For_Writing_Ascii('bench-data.dat', fu)
    ! With circularity 2D:
!    write(fu,'(5(2x,e16.10e2))') time, b_volume,                    &
!                                2.0*PI/surface*sqrt(b_volume/PI),  &
!                                c_position,                        &
!                                rise_velocity
    ! With sphericity 3D:
    write(fu,'(6(2x,e16.10e2))') time, b_volume,                         &
                                PI**(1.0/3.0)*(6.0*b_volume)**(2.0/3.0)  &
                                /surface,                                &
                                c_position,                              &
                                rise_velocity,                           &
                                w_in
    close(fu)
  end if

  end subroutine
