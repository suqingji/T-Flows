!==============================================================================!
  subroutine Turb_Mod_Src_Kin_K_Eps_Zeta_F(turb, Sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation.                       !
!------------------------------------------------------------------------------!
!   In kinetic energy equation there are two source terms:                     !
!                                                                              !
!     /
!    |                                                                         !
!    | (density (p_kin - eps)) dV                                              !
!    |                                                                         !
!   /                                                                          !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: Sol
!---------------------------------[Calling]------------------------------------!
  real :: Roughness_Coefficient
  real :: Tau_Wall_Log_Law
  real :: Y_Plus
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Var_Type),    pointer :: kin, eps, zeta, f, ut, vt, wt, t2
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: u_tan, u_tau
  real                       :: lf, ebf, p_kin_int, p_kin_wf, l_rans_d, l_rans_v
  real                       :: kin_vis
  real                       :: z_o, alpha_d, alpha_v, l_sgs_d, l_sgs_v
  real                       :: ut_log_law, vt_log_law, wt_log_law
  real                       :: nx, ny, nz, qx, qy, qz, g_buoy_wall
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f)
  call Turb_Mod_Alias_Heat_Fluxes (turb, ut, vt, wt)
  call Sol % Alias_Solver         (A, b)
  call Turb_Mod_Alias_T2          (turb, t2)

  ! Production source:
  do c = 1, Grid % n_cells
    turb % p_kin(c) = max(turb % vis_t(c) * Flow % shear(c)**2, TINY)
    b(c) = b(c) + turb % p_kin(c) * Grid % vol(c)
  end do

  if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
    do c = 1, Grid % n_cells
      turb % g_buoy(c) = -Flow % beta                    &
                        * ( Flow % grav_x * ut % n(c)    &
                          + Flow % grav_y * vt % n(c)    &
                          + Flow % grav_z * wt % n(c) )  &
                        * Flow % density(c)

! In general, this clipping should be avoided.  
!      if(turb % g_buoy(c) + turb % p_kin(c) < 0.0) then
!        turb % g_buoy(c) = 0.0
!      end if

      b(c) = b(c) + max(0.0, turb % g_buoy(c) * Grid % vol(c))
      A % val(A % dia(c)) = A % val(A % dia(c))         &
                          + max(0.0,-turb % g_buoy(c)   &
                          * Grid % vol(c)               &
                          / (kin % n(c) + TINY))
    end do
  end if

  if(turb % model .eq. HYBRID_LES_RANS) then
    do c = 1, Grid % n_cells

      lf = Grid % vol(c)**ONE_THIRD

      ! Distance switch
      l_sgs_d  = 0.8 * lf
      l_rans_d = 0.41 * Grid % wall_dist(c)
      alpha_d  = max(1.0,l_rans_d/l_sgs_d)

      ! Velocity switch
      l_sgs_v  = lf * Flow % shear(c)
      l_rans_v = sqrt(kin % n(c) * zeta % n(c))
      alpha_v  = l_rans_v / (l_sgs_v + TINY)

      if( (turb % hybrid_les_rans_switch .eq. SWITCH_DISTANCE)  &
          .and. (alpha_d < 1.05)                                &
          .or.                                                  &
          (turb % hybrid_les_rans_switch .eq. SWITCH_VELOCITY)  &
          .and. (alpha_v < 0.5 .or. alpha_d < 1.05) ) then
        A % val(A % dia(c)) = A % val(A % dia(c))             &
                            + Flow % density(c) * eps % n(c)  &
                            / (kin % n(c) + TINY) * Grid % vol(c)
      else
        A % val(A % dia(c)) = A % val(A % dia(c))                        &
          + Flow % density(c)                                            &
          * min(alpha_d**1.4 * eps % n(c), kin % n(c)**1.5 / (lf*0.01))  &
          / (kin % n(c) + TINY) * Grid % vol(c)
      end if
    end do
  else  ! turbuence model will be K_EPS_ZETA_F
    do c = 1, Grid % n_cells
      A % val(A % dia(c)) = A % val(A % dia(c))             &
                          + Flow % density(c) * eps % n(c)  &
                          / (kin % n(c) + TINY) * Grid % vol(c)

    end do
  end if

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then
      ! Kinematic viscosities
      kin_vis = Flow % viscosity(c1) / Flow % density(c1)

      if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or. &
         Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

        ! Set up roughness coefficient 
        z_o = Roughness_Coefficient(turb, turb % z_o_f(c1))
        if(turb % rough_walls) then
          z_o = max(Grid % wall_dist(c1)   &
              / (e_log * max(turb % y_plus(c1), 1.0)), z_o)
        end if

        ! Compute tangential velocity component
        u_tan = Flow % U_Tan(s)

        u_tau = c_mu25 * sqrt(kin % n(c1))
     
        turb % y_plus(c1) = Y_Plus(turb,                 &
                                   u_tau,                &
                                   Grid % wall_dist(c1), &
                                   kin_vis,              &
                                   z_o)

        turb % tau_wall(c1) = Tau_Wall_Log_Law(turb,                &
                                              Flow % density(c1),   &
                                              u_tau,                &
                                              u_tan,                &
                                              Grid % wall_dist(c1), &
                                              turb % y_plus(c1),    &
                                              z_o)

        ebf = Turb_Mod_Ebf_Momentum(turb, c1)

        p_kin_wf  = turb % tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))  &
                  / ((Grid % wall_dist(c1) + z_o) * kappa)

        p_kin_int = turb % vis_t(c1) * Flow % shear(c1)**2

        turb % p_kin(c1) = exp(-1.0 * ebf) * p_kin_int   &
                         + exp(-1.0 / ebf) * p_kin_wf

        b(c1) = b(c1) + (turb % p_kin(c1)  &
              - turb % vis_t(c1) * Flow % shear(c1)**2) * Grid % vol(c1)

        ! Implementation of wall function for buoyancy-driven flows
        if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then

          nx = Grid % sx(s) / Grid % s(s)
          ny = Grid % sy(s) / Grid % s(s)
          nz = Grid % sz(s) / Grid % s(s)
          qx = t % q(c2) * nx
          qy = t % q(c2) * ny
          qz = t % q(c2) * nz

          ut_log_law = - turb % con_w(c1)  &
                     / (Flow % density(c1) * Flow % capacity(c1))   &
                     * (t % n(c2) - t % n(c1))/Grid % wall_dist(c1) * nx
          vt_log_law = - turb % con_w(c1)  &
                     / (Flow % density(c1) * Flow % capacity(c1))   &
                     * (t % n(c2) - t % n(c1))/Grid % wall_dist(c1) * ny
          wt_log_law = - turb % con_w(c1)  &
                     / (Flow % density(c1) * Flow % capacity(c1))   &
                     * (t % n(c2) - t % n(c1))/Grid % wall_dist(c1) * nz

          ut % n(c1) = ut % n(c1) * exp(-1.0 * ebf)  &
                     + ut_log_law * exp(-1.0 / ebf)
          vt % n(c1) = vt % n(c1) * exp(-1.0 * ebf)  &
                     + vt_log_law * exp(-1.0 / ebf)
          wt % n(c1) = wt % n(c1) * exp(-1.0 * ebf)  &
                     + wt_log_law * exp(-1.0 / ebf)

          if(Grid % Bnd_Cond_Type(c2) .eq. WALL)             &
            t % q(c2) = turb % con_w(c1) * (t % n(c1) - t % n(c2))  &
                      / Grid % wall_dist(c1)

          g_buoy_wall = Flow % density(c1)                                    &
                      * Flow % beta * abs(  Flow % grav_x                     &
                                          + Flow % grav_y                     &
                                          + Flow % grav_z)                    &
                      * sqrt(abs(  t % q(c2)                                  &
                                 / (Flow % density(c1)*Flow % capacity(c1)))  &
                             * c_mu_theta5                                    &
                             * sqrt(abs(t2 % n(c1) * kin % n(c1))))

          ! Clean up b(c) from old values of g_buoy
          b(c1)      = b(c1) - turb % g_buoy(c1) * Grid % vol(c1)

          turb % g_buoy(c1) = turb % g_buoy(c1) * exp(-1.0 * ebf) &
                            + g_buoy_wall * exp(-1.0 / ebf)

          ! Add new values of g_buoy based on wall function approach
          b(c1)      = b(c1) + turb % g_buoy(c1) * Grid % vol(c1)

        end if ! Flow % buoyancy .eq. THERMALLY_DRIVEN

      end if  ! Grid % Bnd_Cond_Type(c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  end subroutine
