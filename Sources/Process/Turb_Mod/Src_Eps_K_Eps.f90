!==============================================================================!
  subroutine Turb_Mod_Src_Eps_K_Eps(turb, Sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in the eps transport equation,                   !
!   wall shear stress (wall function approuch)                                 !
!------------------------------------------------------------------------------!
!   int( density (c_1e eps/kin Gk - c_2e eps^2/kin) )dV                        !
!                                                                              !
!   assigns epsilon from the wall function:                                    !
!                                                                              !
!   Eps_w = Cmu^(3/4)* Kin^(3/2)/(Kappa/y)                                     !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: Sol
!---------------------------------[Calling]------------------------------------!
  real :: Roughness_Coefficient
  real :: Tau_Wall_Log_Law
  real :: Y_Plus
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, eps
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2, j
  real                       :: u_tan, u_tau
  real                       :: re_t, f_mu, u_tau_new, fa, kin_vis
  real                       :: eps_wf, eps_int, y_star, dia_coef_tmp
  real                       :: p_kin_int, p_kin_wf, ebf, z_o
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  call Turb_Mod_Alias_K_Eps    (turb, kin, eps)
  call Sol % Alias_Solver      (A, b)

  do c = 1, Grid % n_cells
    kin_vis =  Flow % viscosity(c) / Flow % density(c)

    ! Positive contribution:
    b(c) = b(c) + &
            c_1e * turb % p_kin(c) * eps % n(c) / kin % n(c) * Grid % vol(c)

    ! Negative contribution:
    re_t = kin % n(c)*kin % n(c)/(kin_vis*eps % n(c))
    y_star = sqrt(sqrt(kin_vis * eps % n(c))) *     &
             Grid % wall_dist(c)/kin_vis
    f_mu = (1.0 - exp(-y_star/3.1))**2              &
         * (1.0 - 0.3*exp(-(re_t/6.5)*(re_t/6.5)))

    f_mu = min(f_mu,1.0)

    A % val(A % dia(c)) = A % val(A % dia(c))                             &
                        +    Flow % density(c) * f_mu* c_2e * eps % n(c)  &
                           / kin % n(c) * Grid % vol(c)

    ! Buoyancy contribution
    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      b(c) = b(c) + max(0.0, c_1e * turb % g_buoy(c) &
                    * eps % n(c) / kin % n(c) * Grid % vol(c))
      A % val(A % dia(c)) = A % val(A % dia(c))                &
                          + max(0.0,(-c_1e * turb % g_buoy(c)  &
                          * eps % n(c)                         &
                          / kin % n(c) * Grid % vol(c))        &
                          / (eps % n(c) + TINY))
    end if

  end do

  ! Imposing a boundary condition on wall for eps
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 < 0) then
      if( Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
          Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

        ! Set up roughness coefficient
        z_o = Roughness_Coefficient(turb, turb % z_o_f(c1))
        if(turb % rough_walls) then
          z_o = max(Grid % wall_dist(c1)  &
              / (e_log * max(turb % y_plus(c1), 1.0)), z_o)
        end if 

        ! Compute tangential velocity component
        u_tan = Flow % U_Tan(s)

        u_tau = c_mu25 * sqrt(kin % n(c1))

        turb % y_plus(c1) = Y_Plus(turb,                  &
                                   u_tau,                 &
                                   Grid % wall_dist(c1),  &
                                   kin_vis,               &
                                   z_o)

        turb % tau_wall(c1) = Tau_Wall_Log_Law(turb,                 &
                                              Flow % density(c1),    &
                                              u_tau,                 &
                                              u_tan,                 &
                                              Grid % wall_dist(c1),  &
                                              turb % y_plus(c1),     &
                                              z_o)

        u_tau_new = sqrt(turb % tau_wall(c1) / Flow % density(c1))

        turb % y_plus(c1) = Y_Plus(turb,                  &
                                   u_tau_new,             &
                                   Grid % wall_dist(c1),  &
                                   kin_vis,               &
                                   z_o)

        eps_int = 2.0* kin_vis * kin % n(c1)  &
                / Grid % wall_dist(c1)**2

        eps_wf  = c_mu75 * kin % n(c1)**1.5   &
                / ((Grid % wall_dist(c1) + z_o) * kappa)

        ebf = Turb_Mod_Ebf_Momentum(turb, c1)

        p_kin_wf  = turb % tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))  &
                / ((Grid % wall_dist(c1) + z_o) * kappa)

        p_kin_int = turb % vis_t(c1) * Flow % shear(c1)**2

        turb % p_kin(c1) = p_kin_int * exp(-1.0 * ebf) + p_kin_wf  &
                           * exp(-1.0 / ebf)

        fa = min( p_kin_wf * exp(-1.0 / ebf) / (turb % p_kin(c1) + TINY), 1.0)

        eps % n(c1) = (1.0 - fa)**0.5 * eps_int + fa**0.5 * eps_wf

        if(turb % y_plus(c1) > 3) then
    
          dia_coef_tmp = A % val(A % dia(c1))

          ! Adjusting coefficient to fix eps value in near wall calls
          do j = A % row(c1), A % row(c1 + 1) - 1 
            A % val(j) = 0.0 
          end do

          b(c1) = eps % n(c1) * dia_coef_tmp
          A % val(A % dia(c1)) = dia_coef_tmp

        else

          eps % n(c2) = 2.0* kin_vis * kin % n(c1)  &
                      / Grid % wall_dist(c1)**2
        end if  ! y_plus(c1) < 3
      end if     ! wall or wall_flux
    end if       ! c2 < 0
  end do

  end subroutine
