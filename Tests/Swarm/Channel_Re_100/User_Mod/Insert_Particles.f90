!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, turb, Vof, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: swarm
  integer, intent(in)      :: n     ! time step
  real,    intent(in)      :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  integer :: i, j, k, n_parts_in_buffers
  real    :: x, y, z, dy, dz
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NJ = 10
  integer, parameter :: NZ = 10
!==============================================================================!

  !-------------------!
  !   1st time step   !
  !-------------------!
  if(n .eq. 1) then

    dy = 1.0 / NJ
    dz = 1.0 / NZ

    ! Place 100 particles where you want them
    do j = 1, NJ
      do k = 1, NZ
        i = (k-1)*NJ + j  ! particle number

        ! Placing particles (only at the 1st time step)
        x = 0.05
        y = dy * 0.5 + (j-1) * dy
        z = dz * 0.5 + (k-1) * dz

        swarm % Particle(i) % x_n = x
        swarm % Particle(i) % y_n = y
        swarm % Particle(i) % z_n = z

        ! Searching for the closest cell and node to place the moved Particle
        call swarm % Particle(i) % Find_Nearest_Cell(n_parts_in_buffers)
        call swarm % Particle(i) % Find_Nearest_Node()
      end do
    end do

  end if

  end subroutine
