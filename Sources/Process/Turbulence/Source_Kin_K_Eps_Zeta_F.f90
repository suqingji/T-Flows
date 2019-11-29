!==============================================================================!
  subroutine Source_Kin_K_Eps_Zeta_F(flow, sol)
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
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Field_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod,   only: Grid_Type
  use Solver_Mod, only: Solver_Type
  use Matrix_Mod, only: Matrix_Type
  use Field_Mod,  only: Field_Type, beta_tec
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
  real :: Y_Plus_Rough_Walls
  real :: Roughness_Coefficient
  real :: Thermal_Expansion_Coefficient
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: t, u, v, w
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: u_tan, u_nor_sq, u_nor, u_tot_sq
  real                       :: lf, ebf, p_kin_int, p_kin_wf
  real                       :: alpha1, l_rans, l_sgs, kin_vis
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
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t
  a    => sol % a
  b    => sol % b % val

  ! Production source:
  do c = 1, grid % n_cells
    p_kin(c) = vis_t(c) * shear(c)**2
    b(c)     = b(c) + p_kin(c) * grid % vol(c)
  end do

  if(turbulence_model .eq. HYBRID_LES_RANS) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD
      l_sgs  = 0.8*lf
      l_rans = 0.41*grid % wall_dist(c)
      alpha1 = max(1.0,l_rans/l_sgs)

      if(alpha1 < 1.05) then
        a % val(a % dia(c)) = a % val(a % dia(c))   &
                            + density * eps % n(c)  &
                            / (kin % n(c) + TINY) * grid % vol(c)
      else
        a % val(a % dia(c)) = a % val(a % dia(c))   &
          + density                                 &
          * min(alpha1**1.45 * eps % n(c), kin % n(c)**1.5 / (lf*0.01))  &
          / (kin % n(c) + TINY) * grid % vol(c)
      end if
    end do
  else  ! turbuence model will be K_EPS_ZETA_F
    do c = 1, grid % n_cells
      a % val(a % dia(c)) = a % val(a % dia(c))   &
                          + density * eps % n(c)  &
                          / (kin % n(c) + TINY) * grid % vol(c)

      if(buoyancy) then
        g_buoy(c) = -beta_tec * (grav_x * ut % n(c) +  &
                                 grav_y * vt % n(c) +  &
                                 grav_z * wt % n(c)) * density

        g_buoy(c) = max(g_buoy(c) ,0.0)

        b(c) = b(c) + g_buoy(c) * grid % vol(c)
      end if
    end do
  end if

  ! Kinematic viscosities
  kin_vis = viscosity / density

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or. &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        ! Compute tangential velocity component
        u_tot_sq = u % n(c1) * u % n(c1) &
                 + v % n(c1) * v % n(c1) &
                 + w % n(c1) * w % n(c1)
        u_nor  = ( u % n(c1) * grid % sx(s)     &
                 + v % n(c1) * grid % sy(s)     &
                 + w % n(c1) * grid % sz(s) )   &
                 / sqrt(  grid % sx(s)*grid % sx(s)  &
                        + grid % sy(s)*grid % sy(s)  &
                        + grid % sz(s)*grid % sz(s))
        u_nor_sq = u_nor**2

        if( u_tot_sq  > u_nor_sq) then
          u_tan = sqrt(u_tot_sq - u_nor_sq)
        else
          u_tan = TINY
        end if

        u_tau(c1) = c_mu25 * sqrt(kin % n(c1))
        y_plus(c1) = Y_Plus_Low_Re(u_tau(c1), grid % wall_dist(c1), kin_vis)

        tau_wall(c1) = density*kappa*u_tau(c1)*u_tan  &
                     / log(e_log*max(y_plus(c1),1.05))

        ebf = max(0.01 * y_plus(c1)**4 / (1.0 + 5.0*y_plus(c1)), TINY)

        p_kin_wf  = tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))  &
                    / (grid % wall_dist(c1) * kappa)

        p_kin_int = vis_t(c1) * shear(c1)**2

        p_kin(c1) = p_kin_wf

        if(rough_walls) then
          z_o = Roughness_Coefficient(grid, z_o_f(c1), c1)    
          y_plus(c1) = Y_Plus_Rough_Walls(u_tau(c1), &
                       grid % wall_dist(c1), kin_vis) 

          tau_wall(c1) = density*kappa*u_tau(c1)*u_tan  &
                       / log(((grid % wall_dist(c1)+z_o) / z_o))

          p_kin(c1) = tau_wall(c1) * c_mu25 * sqrt(kin % n(c1)) &
                      / (kappa*(grid % wall_dist(c1)+z_o))
        end if ! rough_walls

        b(c1) = b(c1) + (p_kin(c1) - vis_t(c1) * shear(c1)**2) * grid % vol(c1)

        ! Implementation of wall function for buoyancy-driven flows
        if(buoyancy) then
          
          nx = grid % sx(s) / grid % s(s)
          ny = grid % sy(s) / grid % s(s)
          nz = grid % sz(s) / grid % s(s)
          qx = t % q(c2) * nx
          qy = t % q(c2) * ny
          qz = t % q(c2) * nz

          ut_log_law = - con_wall(c1) &
                     * (t % n(c2) - t % n(c1))/grid % wall_dist(c1) * nx
          vt_log_law = - con_wall(c1) &
                     * (t % n(c2) - t % n(c1))/grid % wall_dist(c1) * ny
          wt_log_law = - con_wall(c1) &
                     * (t % n(c2) - t % n(c1))/grid % wall_dist(c1) * nz

          ut % n(c1) = ut %n(c1)  * exp(-1.0 * ebf) &
                     + ut_log_law * exp(-1.0 / ebf)
          vt % n(c1) = vt %n(c1)  * exp(-1.0 * ebf) &
                     + vt_log_law * exp(-1.0 / ebf)
          wt % n(c1) = wt %n(c1)  * exp(-1.0 * ebf) &
                     + wt_log_law * exp(-1.0 / ebf)

          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) &
          t % q(c2) = con_wall(c1)*(t % n(c1) &
                      - t % n(c2)) / grid % wall_dist(c1)

          g_buoy_wall = beta_tec*abs(grav_z)*sqrt(abs(t % q(c2))*       &
                        c_mu_theta5*sqrt(abs(t2 % n(c1) * kin % n(c1))))
         
          ! Clean up b(c) from old values of g_buoy         
          b(c1)      = b(c1) - g_buoy(c1) * grid % vol(c1)

          g_buoy(c1) = g_buoy(c1) * exp(-1.0 * ebf) &
                     + g_buoy_wall * exp(-1.0 / ebf)

          ! Add new values of g_buoy based on wall function approach          
          b(c1)      = b(c1) + g_buoy(c1) * grid % vol(c1)
        end if ! buoyancy

      end if   ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if     ! c2 < 0
  end do

  end subroutine
