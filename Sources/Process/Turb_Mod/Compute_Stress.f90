!==============================================================================!
  subroutine Turb_Mod_Compute_Stress(turb, Sol, curr_dt, ini, phi)
!------------------------------------------------------------------------------!
!   Discretizes and solves transport equation for Re stresses for RSM.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: phi_x       => r_cell_01,  &
                      phi_y       => r_cell_02,  &
                      phi_z       => r_cell_03,  &
                      u1uj_phij   => r_cell_04,  &
                      u2uj_phij   => r_cell_05,  &
                      u3uj_phij   => r_cell_06,  &
                      u1uj_phij_x => r_cell_07,  &
                      u2uj_phij_y => r_cell_08,  &
                      u3uj_phij_z => r_cell_09
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence should be outlined                   !
!                                                                              !
!   Main_Pro                                  (allocates Work_Mod)             !
!     |                                                                        !
!     +----> Turb_Mod_Main                    (does not use Work_Mod)          !
!              |                                                               !
!              +---> Turb_Mod_Compute_Stress  (safe to use r_cell_01..09)      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: Sol
  integer, intent(in)       :: curr_dt
  integer, intent(in)       :: ini
  type(Var_Type)            :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, eps, zeta, f22, ut, vt, wt
  type(Var_Type),    pointer :: uu, vv, ww, uv, uw, vw
  real, contiguous,  pointer :: flux(:)
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2, nc, nb
  real                       :: f_ex, f_im
  real                       :: a0, a12, a21
  real                       :: vis_eff
  real                       :: phix_f, phiy_f, phiz_f
  real                       :: vis_t_f
  real                       :: dt
  real                       :: visc_f
!==============================================================================!
!                                                                              !
!   The form of equations which are being solved:                              !
!                                                                              !
!     /               /                /                     /                 !
!    |     dphi      |                | mu_eff              |                  !
!    | rho ---- dV + | rho u phi dS = | ------ DIV phi dS + | G dV             !
!    |      dt       |                |  sigma              |                  !
!   /               /                /                     /                   !
!                                                                              !
!------------------------------------------------------------------------------!

  call Cpu_Timer % Start('Compute_Turbulence (without solvers)')

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  nc   =  Grid % n_cells
  nb   =  Grid % n_bnd_cells
  dt   =  Flow % dt
  flux => Flow % v_flux % n
  call Flow % Alias_Momentum(u, v, w)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_Heat_Fluxes (turb, ut, vt, wt)
  call Sol % Alias_Solver         (A, b)

  ! Initialize advection and cross diffusion sources, matrix and right hand side
  phi % a(:) = 0.0
  phi % c(:) = 0.0
  A % val(:) = 0.0
  b      (:) = 0.0

  ! Old values (o) and older than old (oo)
  if(ini .eq. 1) then
    do c = 1, Grid % n_cells
      phi % oo(c) = phi % o(c)
      phi % o (c) = phi % n(c)
    end do
  end if

  ! Gradients
  call Flow % Grad(phi % n, phi_x(-nb:nc),  &
                            phi_y(-nb:nc),  &
                            phi_z(-nb:nc))

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(phi, Flow % density, flux, b)

  !------------------!
  !                  !
  !     Difusion     !
  !                  !
  !------------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! vis_tur is used to make diaginal element more dominant.
    ! This contribution is later substracted.
    vis_t_f =      Grid % fw(s)  * turb % vis_t(c1)  &
            + (1.0-Grid % fw(s)) * turb % vis_t(c2)

    visc_f =        Grid % fw(s)  * Flow % viscosity(c1)  &
           + (1.0 - Grid % fw(s)) * Flow % viscosity(c2)

    vis_eff = visc_f + vis_t_f

    if(turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
      if(turb % model_variant .ne. STABILIZED) then
        vis_eff = 1.5 * visc_f + vis_t_f
      end if
    end if

    phix_f = Grid % fw(s) * phi_x(c1) + (1.0-Grid % fw(s)) * phi_x(c2)
    phiy_f = Grid % fw(s) * phi_y(c1) + (1.0-Grid % fw(s)) * phi_y(c2)
    phiz_f = Grid % fw(s) * phi_z(c1) + (1.0-Grid % fw(s)) * phi_z(c2)


    ! Total (exact) diffusive flux plus turb. diffusion
    f_ex = vis_eff * (  phix_f * Grid % sx(s)  &
                      + phiy_f * Grid % sy(s)  &
                      + phiz_f * Grid % sz(s) ) 

    a0 = vis_eff * A % fc(s)

    ! Implicit diffusive flux
    ! (this is a very crude approximation: f_coef is
    !  not corrected at interface between materials)
    f_im=( phix_f*Grid % dx(s)       &
          +phiy_f*Grid % dy(s)       &
          +phiz_f*Grid % dz(s))*a0

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + f_ex - f_im
    if(c2  > 0) then
      phi % c(c2) = phi % c(c2) - f_ex + f_im
    end if

    ! Compute coefficients for the sysytem matrix
    a12 = a0
    a21 = a0

    a12 = a12  - min(flux(s), 0.0) * Flow % density(c1)
    a21 = a21  + max(flux(s), 0.0) * Flow % density(c2)

    ! Fill the system matrix
    if(c2  > 0) then
      A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - a12
      A % val(A % dia(c1))  = A % val(A % dia(c1))  + a12
      A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - a21
      A % val(A % dia(c2))  = A % val(A % dia(c2))  + a21
    else if(c2  < 0) then

      ! Outflow is not included because it was causing problems     
      ! Convect is commented because for turbulent scalars convect 
      ! outflow is treated as classic outflow.
      if((Grid % Bnd_Cond_Type(c2) .eq. INFLOW).or.     &
         (Grid % Bnd_Cond_Type(c2) .eq. WALL).or.       &
!!!      (Grid % Bnd_Cond_Type(c2) .eq. CONVECT).or.    &
         (Grid % Bnd_Cond_Type(c2) .eq. WALLFL) ) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
        b(c1) = b(c1) + a12 * phi % n(c2)
      end if
    end if

  end do  ! through faces

  !------------------------------!
  !   Turbulent diffusion term   !
  !------------------------------!
  if(phi % name .eq. 'EPS') then
    c_mu_d = 0.18
  else
    c_mu_d = 0.22
  end if

  if(turb % model_variant .ne. STABILIZED) then
    if(turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
      do c = 1, Grid % n_cells
        u1uj_phij(c) = Flow % density(c) * c_mu_d / phi % sigma        &
                     * kin % n(c)                                      &
                     / max(eps % n(c), TINY)                           &
                     * (  uu % n(c) * phi_x(c)                         &
                        + uv % n(c) * phi_y(c)                         &
                        + uw % n(c) * phi_z(c))                        &
                     - Flow % viscosity(c) * phi_x(c)

        u2uj_phij(c) = Flow % density(c) * c_mu_d / phi % sigma        &
                     * kin % n(c)                                      &
                     / max(eps % n(c), TINY)                           &
                     * (  uv % n(c) * phi_x(c)                         &
                        + vv % n(c) * phi_y(c)                         &
                        + vw % n(c) * phi_z(c))                        &
                     - Flow % viscosity(c) * phi_y(c)

        u3uj_phij(c) = Flow % density(c) * c_mu_d / phi % sigma        &
                     * kin % n(c)                                      &
                     / max(eps % n(c), TINY)                           &
                     * (  uw % n(c) * phi_x(c)                         &
                        + vw % n(c) * phi_y(c)                         &
                        + ww % n(c) * phi_z(c))                        &
                     - Flow % viscosity(c) * phi_z(c)
      end do
    else if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
      do c = 1, Grid % n_cells
        u1uj_phij(c) = Flow % density(c) * c_mu_d / phi % sigma            &
                     * turb % t_scale(c)                                   &
                     * (  uu % n(c) * phi_x(c)                             &
                        + uv % n(c) * phi_y(c)                             &
                        + uw % n(c) * phi_z(c))

        u2uj_phij(c) = Flow % density(c) * c_mu_d / phi % sigma            &
                     * turb % t_scale(c)                                   &
                     * (  uv % n(c) * phi_x(c)                             &
                        + vv % n(c) * phi_y(c)                             &
                        + vw % n(c) * phi_z(c))

        u3uj_phij(c) = Flow % density(c) * c_mu_d / phi % sigma            &
                     * turb % t_scale(c)                                   &
                     * (  uw % n(c) * phi_x(c)                             &
                        + vw % n(c) * phi_y(c)                             &
                        + ww % n(c) * phi_z(c))
      end do
    end if

    call Flow % Grad_Component(u1uj_phij(-nb:nc), 1, u1uj_phij_x(-nb:nc))
    call Flow % Grad_Component(u2uj_phij(-nb:nc), 2, u2uj_phij_y(-nb:nc))
    call Flow % Grad_Component(u3uj_phij(-nb:nc), 3, u3uj_phij_z(-nb:nc))

    do c = 1, Grid % n_cells
      b(c) = b(c) + (  u1uj_phij_x(c)  &
                     + u2uj_phij_y(c)  &
                     + u3uj_phij_z(c) ) * Grid % vol(c)
    end do
  end if

  !------------------------------------------------------------------!
  !   Here we clean up transport equation from the false diffusion   !
  !------------------------------------------------------------------!
  if(turb % model_variant .ne. STABILIZED) then
    do s = 1, Grid % n_faces

      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      vis_eff = (Grid % fw(s)      * turb % vis_t(c1)  &
              + (1.0-Grid % fw(s)) * turb % vis_t(c2))

      phix_f = Grid % fw(s) * phi_x(c1) + (1.0-Grid % fw(s)) * phi_x(c2)
      phiy_f = Grid % fw(s) * phi_y(c1) + (1.0-Grid % fw(s)) * phi_y(c2)
      phiz_f = Grid % fw(s) * phi_z(c1) + (1.0-Grid % fw(s)) * phi_z(c2)
      f_ex = vis_eff * (  phix_f * Grid % sx(s)  &
                        + phiy_f * Grid % sy(s)  &
                        + phiz_f * Grid % sz(s))
      a0 = vis_eff * A % fc(s)
      f_im = (   phix_f * Grid % dx(s)        &
               + phiy_f * Grid % dy(s)        &
               + phiz_f * Grid % dz(s)) * a0

      b(c1) = b(c1)                                             &
             - vis_eff * (phi % n(c2) - phi%n(c1)) * A % fc(s)  &
             - f_ex + f_im
      if(c2  > 0) then
        b(c2) = b(c2)                                            &
              + vis_eff * (phi % n(c2) - phi%n(c1)) * A % fc(s)  &
              + f_ex - f_im
      end if
    end do
  end if

  !------------------------------------------------!
  !   Source term contains difference between      !
  !   explicity and implicitly treated advection   !
  !------------------------------------------------!
  do c = 1, Grid % n_cells
    b(c) = b(c) + phi % c(c)
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!
  call Numerics_Mod_Inertial_Term(phi, Flow % density, a, b, dt)

  !-------------------------------------!
  !                                     !
  !   Source terms and wall function    !
  !                                     !
  !-------------------------------------!
  if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
    call Flow % Grad_Variable(f22)

    call Turb_Mod_Src_Rsm_Manceau_Hanjalic(turb, Sol, phi % name)
  else if(turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Turb_Mod_Src_Rsm_Hanjalic_Jakirlic(turb, Sol, phi % name)
  end if

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(phi, a, b)

  ! Call linear solver to solve the equations
  call Cpu_Timer % Start('Linear_Solver_For_Turbulence')
  call Sol % Bicg(A,              &
                  phi % n,        &
                  b,              &
                  phi % precond,  &
                  phi % mniter,   &
                  phi % eniter,   &
                  phi % tol,      &
                  phi % res)
  call Cpu_Timer % Stop('Linear_Solver_For_Turbulence')

  ! Print info on the screen
  if( phi % name .eq. 'UU' )   &
    call Info_Mod_Iter_Fill_At(3, 1, phi % name, phi % eniter, phi % res)
  if( phi % name .eq. 'VV' )   &
    call Info_Mod_Iter_Fill_At(3, 2, phi % name, phi % eniter, phi % res)
  if( phi % name .eq. 'WW' )   &
    call Info_Mod_Iter_Fill_At(3, 3, phi % name, phi % eniter, phi % res)
  if( phi % name .eq. 'UV' )   &
    call Info_Mod_Iter_Fill_At(3, 4, phi % name, phi % eniter, phi % res)
  if( phi % name .eq. 'UW' )   &
    call Info_Mod_Iter_Fill_At(3, 5, phi % name, phi % eniter, phi % res)
  if( phi % name .eq. 'VW' )   &
    call Info_Mod_Iter_Fill_At(3, 6, phi % name, phi % eniter, phi % res)
  if( phi % name .eq. 'EPS' )  &
    call Info_Mod_Iter_Fill_At(4, 1, phi % name, phi % eniter, phi % res)

  if(phi % name .eq. 'EPS') then
    do c= 1, Grid % n_cells
      phi % n(c) = phi % n(c)
     if( phi % n(c) < 0.) then
       phi % n(c) = phi % o(c)
     end if
    end do
  end if

  if(phi % name .eq. 'UU' .or.  &
     phi % name .eq. 'VV' .or.  &
     phi % name .eq. 'WW') then
    do c = 1, Grid % n_cells
      phi % n(c) = phi % n(c)
      if(phi % n(c) < 0.) then
        phi % n(c) = phi % o(c)
      end if
    end do
  end if

  call Flow % Grad_Variable(phi)

  call Cpu_Timer % Stop('Compute_Turbulence (without solvers)')

  end subroutine
