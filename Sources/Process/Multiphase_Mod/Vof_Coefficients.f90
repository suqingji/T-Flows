!==============================================================================!
  subroutine Multiphase_Mod_Vof_Coefficients(mult, a, b, dt, beta_f)
!------------------------------------------------------------------------------!
!   Computes matrix coefficients for volume fraction equation                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Matrix_Type),     target :: a
  real,                  target :: b(:)
  real                          :: dt
  real                          :: beta_f(:)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Face_Type),  pointer :: v_flux
  type(Var_Type),   pointer :: vof
  integer                   :: c, c1, c2, s
  real                      :: upwd1, upwd2, upwd3, a0
!==============================================================================!

  ! Take aliases
  flow   => mult % pnt_flow
  grid   => flow % pnt_grid
  v_flux => flow % v_flux
  vof    => mult % vof

  ! Initialize matrix and right hand side
  b       = 0.0
  a % val = 0.0

  !-------------------------!
  !   Matrix Coefficients   !
  !-------------------------!
  if(vof % adv_scheme .eq. UPWIND) then

    ! At boundaries
    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
        b(c1) = b(c1) - v_flux % n(s) * vof % n(c2)
      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + v_flux % n(s)
      end if

    end do

    ! Interior faces
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      upwd1 = 0.5 * max( v_flux % n(s), 0.0)
      upwd2 = 0.5 * max(-v_flux % n(s), 0.0)

      a % val(a % dia(c1)) = a % val(a % dia(c1)) + upwd1
      a % val(a % dia(c2)) = a % val(a % dia(c2)) + upwd2

      a % val(a % pos(1,s)) =  - upwd2
      a % val(a % pos(2,s)) =  - upwd1

      b(c1) = b(c1) - ( upwd1 * vof % o(c1) - upwd2 * vof % o(c2) )
      b(c2) = b(c2) - ( upwd2 * vof % o(c2) - upwd1 * vof % o(c1) )

    end do

  else if(vof % adv_scheme .eq. CICSAM .or. &
           vof % adv_scheme .eq. STACS) then

    ! At boundaries
    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
        a % val(a % dia(c1)) = a % val(a % dia(c1)) + v_flux % n(s)
      else
        b(c1) = b(c1) - v_flux % n(s) * vof % n(c2)
      end if
    end do

    ! Interior faces
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      upwd1 = (0.5 - beta_f(s)) * max(-v_flux % n(s), 0.0)       &
             - 0.5 * beta_f(s)  * v_flux % n(s)
      upwd2 = (0.5 - beta_f(s)) * max(v_flux % n(s), 0.0)       &
             + 0.5 * beta_f(s)  * v_flux % n(s)
      upwd3 = 0.5 * v_flux % n(s)

      a % val(a % dia(c1)) = a % val(a % dia(c1)) + upwd1 + upwd3
      a % val(a % dia(c2)) = a % val(a % dia(c2)) + upwd2 - upwd3

      a % val(a % pos(1,s)) = - upwd1
      a % val(a % pos(2,s)) = - upwd2

      b(c1) = b(c1) - (upwd1 + upwd3) * vof % o(c1) + upwd1 * vof % o(c2)
      b(c2) = b(c2) - (upwd2 - upwd3) * vof % o(c2) + upwd2 * vof % o(c1)
    end do

  end if

  !------------------------------------------------------------!
  !                Calculate Main coefficient and              !
  !   Calculate Source term (part 2 : Temporal contribution)   !
  !------------------------------------------------------------!

  ! Two time levels; Linear interpolation
  if(vof % td_scheme .eq. LINEAR) then
    do c = 1, grid % n_cells
      a0 = grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + a0
      b(c) = b(c) + a0 * vof % o(c)
    end do
  end if

  ! Three time levels; parabolic interpolation
  if(vof % td_scheme .eq. PARABOLIC) then
    do c = 1, grid % n_cells
      a0 = grid % vol(c) / dt
      a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
      b(c) = b(c) + 2.0 * a0 * vof % o(c) - 0.5 * a0 * vof % oo(c)
    end do
  end if

  ! Phase change
  if(mult % phase_change) then

    !do s = 1, grid % n_bnd_faces
    !  c1 = grid % faces_c(1,s)
    !  c2 = grid % faces_c(2,s)

    !  if(mult % ic(c1) == 1) then
    !    !b(c1) = b(c1) + v_flux % n(s) * vof % n(c1)
    !    a % val(a % dia(c1)) = a % val(a % dia(c1)) - v_flux % n(s)
    !  end if
    !end do

    !! Interior faces
    !do s = grid % n_bnd_faces + 1, grid % n_faces
    !  c1 = grid % faces_c(1,s)
    !  c2 = grid % faces_c(2,s)

    !  if(mult % ic(c1) == 1) then
    !    !b(c1) = b(c1) + v_flux % n(s) * vof % n(c1)
    !    a % val(a % dia(c1)) = a % val(a % dia(c1)) - v_flux % n(s)
    !  end if

    !  if(mult % ic(c2) == 1) then
    !    !b(c2) = b(c2) - v_flux % n(s) * vof % n(c2)
    !    a % val(a % dia(c2)) = a % val(a % dia(c2)) + v_flux % n(s)
    !  end if
    !end do

    do c = 1, grid % n_cells
      b(c) = b(c) + mult % flux_rate(c) * grid % vol(c)    &
                  / mult % phase_dens(1)
      !b(c) = b(c) + mult % flux_rate(c) * grid % vol(c) * vof % n(c)    &
      !          * (1.0 / mult % phase_dens(2) - 1.0 / mult % phase_dens(1))
      !b(c) = b(c) + mult % flux_rate(c) * grid % vol(c)                       &
      !            * ( 1.0 / mult % phase_dens(1)                              &
      !            - vof % n(c)                                                &
      !            * (1.0 / mult % phase_dens(1) - 1.0 / mult % phase_dens(2)) )
      !b(c) = b(c) + mult % flux_rate(c) * grid % vol(c)    &
      !            / mult % phase_dens(1)
      !a % val(a % dia(c)) = a % val(a % dia(c)) + mult % flux_rate(c) * grid % vol(c)   &
      !                    * (1.0 / mult % phase_dens(1) - 1.0 / mult % phase_dens(2))
    end do
  end if

  end subroutine
