!==============================================================================!
  subroutine Swarm_Mod_Check_Periodicity(swarm, k, n_parts_in_buffers)
!------------------------------------------------------------------------------!
!   Check if particle left the the domain in a periodic direction              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k
  integer                  :: n_parts_in_buffers
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: Part
  integer                      :: s, c1, c2
  real                         :: xc1, xc2, yc1, yc2, zc1, zc2, d_sq_1, d_sq_2
!==============================================================================!

  ! Take aliases
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid

  !-------------------------------------!
  !                                     !
  !   Mark cells which have particles   !
  !                                     !
  !-------------------------------------!
  swarm % cell_has_particles(:) = .false.  ! assume no cells has particles

  Part => swarm % particle(k)
  if(Part % proc .eq. this_proc) then
    swarm % cell_has_particles(Part % cell) = .true.
  end if

  !-----------------------------------!
  !                                   !
  !   Browse through periodic faces   !
  !                                   !
  !-----------------------------------!
  do s = 1, grid % n_faces

    if(grid % faces_s(s) .ne. 0) then

      c1 = grid % faces_c(1, s)
      c2 = grid % faces_c(2, s)

      !------------------------------!
      !   If c1 contains particles   !
      !------------------------------!
      if( swarm % cell_has_particles(c1) ) then

        Part => swarm % particle(k)

        if(Part % proc .eq. this_proc) then

          xc1 = grid % xc(c1)
          yc1 = grid % yc(c1)
          zc1 = grid % zc(c1)

          ! Take position of c2 at its shadow position
          xc2 = grid % xc(c1) + grid % dx(s)
          yc2 = grid % yc(c1) + grid % dy(s)
          zc2 = grid % zc(c1) + grid % dz(s)

          d_sq_1 = (xc1 - Part % x_n)**2  &
                 + (yc1 - Part % y_n)**2  &
                 + (zc1 - Part % z_n)**2

          d_sq_2 = (xc2 - Part % x_n)**2  &
                 + (yc2 - Part % y_n)**2  &
                 + (zc2 - Part % z_n)**2

          ! Particle is indeed closer to c2's shadow, mark it as "lost"
          if( d_sq_2 < d_sq_1 ) then

            ! Move c1 to its shadow position
            xc1 = grid % xc(c2) - grid % dx(s)
            yc1 = grid % yc(c2) - grid % dy(s)
            zc1 = grid % zc(c2) - grid % dz(s)

            Part % x_n = xc1 + Part % x_n - grid % xc(c1)
            Part % y_n = yc1 + Part % y_n - grid % yc(c1)
            Part % z_n = zc1 + Part % z_n - grid % zc(c1)
            Part % x_o = xc1 + Part % x_o - grid % xc(c1)
            Part % y_o = yc1 + Part % y_o - grid % yc(c1)
            Part % z_o = zc1 + Part % z_o - grid % zc(c1)
                              !<---- dist(Part, c1) ---->!
            Part % cell = c2
            ! call Swarm_Mod_Find_Nearest_Node(swarm, k)
            call Part % Find_Nearest_Node()

            ! If c2 is in the buffer, tell that particle wants to go there
            if(grid % comm % cell_proc(c2) .ne.  &
               grid % comm % cell_proc(c1)) then
              Part % buff = grid % comm % cell_proc(c2)
              n_parts_in_buffers = n_parts_in_buffers + 1
            end if
          end if
        end if

      end if    ! c1 has particles

      !------------------------------!
      !   If c2 contains particles   !
      !------------------------------!
      if( swarm % cell_has_particles(c2) ) then

        Part => swarm % particle(k)

        if(Part % proc .eq. this_proc) then

          ! Take position of c1 at its shadow position
          xc1 = grid % xc(c2) - grid % dx(s)
          yc1 = grid % yc(c2) - grid % dy(s)
          zc1 = grid % zc(c2) - grid % dz(s)

          xc2 = grid % xc(c2)
          yc2 = grid % yc(c2)
          zc2 = grid % zc(c2)

          d_sq_1 = (xc1 - Part % x_n)**2  &
                 + (yc1 - Part % y_n)**2  &
                 + (zc1 - Part % z_n)**2

          d_sq_2 = (xc2 - Part % x_n)**2  &
                 + (yc2 - Part % y_n)**2  &
                 + (zc2 - Part % z_n)**2

          ! Particle is indeed closer to c1's shadow, mark it as "lost"
          if( d_sq_1 < d_sq_2 ) then

            ! Move c2 to its shadow position
            xc2 = grid % xc(c1) + grid % dx(s)
            yc2 = grid % yc(c1) + grid % dy(s)
            zc2 = grid % zc(c1) + grid % dz(s)

            Part % x_n = xc2 + Part % x_n - grid % xc(c2)
            Part % y_n = yc2 + Part % y_n - grid % yc(c2)
            Part % z_n = zc2 + Part % z_n - grid % zc(c2)
            Part % x_o = xc2 + Part % x_o - grid % xc(c2)
            Part % y_o = yc2 + Part % y_o - grid % yc(c2)
            Part % z_o = zc2 + Part % z_o - grid % zc(c2)
                              !<---- dist(Part, c2) ---->!
            Part % cell = c1
            ! call Swarm_Mod_Find_Nearest_Node(swarm, k)
            call Part % Find_Nearest_Node()

            ! If c1 is in the buffer, tell that particle wants to go there
            if(grid % comm % cell_proc(c1) .ne.  &
               grid % comm % cell_proc(c2)) then
              Part % buff = grid % comm % cell_proc(c1)
              n_parts_in_buffers = n_parts_in_buffers + 1
            end if
          end if
        end if

      end if    ! c2 has particles

    end if  ! is a shadow face

  end do

  end subroutine
