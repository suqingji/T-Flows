!==============================================================================!
  subroutine Multiphase_Mod_Vof_Correct_Beta(mult, grid, beta_f, c_d)
!------------------------------------------------------------------------------!
!   Step 2 of CICSAM: Correct beta for computation of volume fraction          !
!                                                                              !
!   Arguments                                                                  !
!   - beta_f                     : Coefficient beta at faces, governs the      !
!                                  advection of vof using CICSAM               !
!   - c_d                        : Courant number at cells                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Grid_Type)               :: grid
  type(Var_Type)                :: phi
  real                          :: beta_f(grid % n_faces)
  real                          :: c_d   (-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),     pointer :: flow
  type(Var_Type),       pointer :: vof
  type(Face_Type),      pointer :: v_flux
  integer                       :: s, c1, c2, donor, accept
  real                          :: fs, e_plus, e_minus, cf, delta_alfa, bcorr
  real                          :: epsloc
!==============================================================================!

  epsloc = epsilon(epsloc)

  ! Take aliases
  flow   => mult % pnt_flow
  vof    => mult % vof
  v_flux => flow % v_flux

  ! Interior faces
  do s = grid % n_bnd_faces + 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)
    if (abs(v_flux % n(s)) > epsloc) then

      if (v_flux % n(s) > 0.0) then
        donor = c1
        accept = c2
      else
        donor = c2
        accept = c1
      end if

      !--------------------!
      !   Correct beta_f   !
      !--------------------!
      bcorr = 0.0
      delta_alfa = 0.5 * (vof % n(accept) + vof % o(accept)      &
                       - (vof % n(donor) + vof % o(donor)))

      cf = c_d(donor)

      if (vof % n(donor) < 0.0) then
        e_minus = max(-vof % n(donor), 0.0)
        ! Donor value < 0.0 Ex: sd = -0.1 -> e_minus = +0.1
        if (e_minus > epsloc .and. cf > epsloc) then
          if (delta_alfa > e_minus) then
            bcorr = e_minus * (2.0 + cf - 2.0 * cf * beta_f(s))     &
                            / (2.0 * cf * (delta_alfa - e_minus))
            bcorr = min(bcorr, beta_f(s))
          end if
        end if
      end if

      if (vof % n(donor) > 1.0) then
        e_plus = max(vof % n(donor) - 1.0, 0.0)
        ! Donor value > 1.0 Ex: sd = 1.1 -> e_plus = +0.1
        if (e_plus > epsloc .and. cf > epsloc) then
          if (delta_alfa < - e_plus) then
            bcorr = e_plus * (2.0 + cf - 2.0 * cf * beta_f(s))     &
                           / (2.0 * cf * (-delta_alfa - e_plus))
            bcorr = min(bcorr, beta_f(s))
          end if
        end if

      end if

      beta_f(s) = beta_f(s) - bcorr
      beta_f(s) = max(beta_f(s),0.0)

    end if

  end do

  end subroutine
