!==============================================================================!
  subroutine Compute_Scalar(Flow, turb, Vof, Sol, curr_dt, ini, sc)
!------------------------------------------------------------------------------!
!   Purpose: Solve transport equation for use scalar.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
  use Work_Mod, only: u1uj_phij   => r_cell_01,  &
                      u2uj_phij   => r_cell_02,  &
                      u3uj_phij   => r_cell_03,  &
                      u1uj_phij_x => r_cell_04,  &
                      u2uj_phij_y => r_cell_05,  &
                      u3uj_phij_z => r_cell_06
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence should be outlined                   !
!                                                                              !
!   Main_Pro                (allocates Work_Mod)                               !
!     |                                                                        |
!     +----> Compute_Scalar (safe to use r_cell_01..06)
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Arguments]--------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt
  integer, intent(in)         :: ini
  integer, intent(in)         :: sc
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: uu, vv, ww, uv, uw, vw
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  type(Face_Type),   pointer :: v_flux
  type(Var_Type),    pointer :: phi
  integer                    :: c, s, c1, c2, row, col
  real                       :: a12, a21
  real                       :: ns, dt
  real                       :: dif_eff1, f_ex1, f_im1
  real                       :: dif_eff2, f_ex2, f_im2
  real                       :: phix_f1, phiy_f1, phiz_f1
  real                       :: phix_f2, phiy_f2, phiz_f2
!------------------------------------------------------------------------------!
!
!  The form of equations which are solved:
!
!     /                /                /
!    |     d phi      |                |
!    | rho ----- dV   | rho u phi dS = | gamma DIV phi dS
!    |      dt        |                |
!   /                /                /
!
!==============================================================================!

  call Cpu_Timer % Start('Compute_Scalars (without solvers)')

  ! Take aliases
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  phi    => Flow % scalar(sc)
  dt     =  Flow % dt
  call Turb_Mod_Alias_Stresses(turb, uu, vv, ww, uv, uw, vw)
  call Sol % Alias_Solver     (A, b)

  ! User function
  call User_Mod_Beginning_Of_Compute_Scalar(Flow, turb, Vof, Sol,  &
                                            curr_dt, ini, sc)

  ! Initialize matrix and right hand side
  A % val(:) = 0.0
  b      (:) = 0.0

  !-------------------------------------!
  !   Initialize variables and fluxes   !
  !-------------------------------------!

  ! Old values (o and oo)
  if(ini.lt.2) then
    do c = 1, Grid % n_cells
      phi % oo(c) = phi % o(c)
      phi % o (c) = phi % n(c)
    end do
  end if

  ! Gradients
  call Flow % Grad_Variable(phi)

  !---------------!
  !               !
  !   Advection   !
  !               !
  !---------------!
  call Numerics_Mod_Advection_Term(phi, Flow % density, v_flux % n, b)

  !--------------!
  !              !
  !   Difusion   !
  !              !
  !--------------!

  !----------------------------!
  !   Spatial discretization   !
  !----------------------------!
  call Control_Mod_Turbulent_Schmidt_Number(sc_t)  ! get default sc_t (0.9)

  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

!   For species, we don't have a function for turbulent Schmidt number
!   if(turb % model .ne. LES_SMAGORINSKY    .or.  &
!      turb % model .ne. LES_DYNAMIC        .or.  &
!      turb % model .ne. HYBRID_LES_PRANDTL .or.  &
!      turb % model .ne. LES_WALE           .or.  &
!      turb % model .ne. DNS) then
!     sc_t1 = Turb_Mod_Prandtl_Number(turb, c1)
!     sc_t2 = Turb_Mod_Prandtl_Number(turb, c2)
!     sc_t  = Grid % fw(s) * sc_t1 + (1.0-Grid % fw(s)) * sc_t2
!   end if

    ! Gradients on the cell face 
    if(c2 > 0) then
      phix_f1 = Grid % fw(s)*phi % x(c1) + (1.0-Grid % fw(s))*phi % x(c2)
      phiy_f1 = Grid % fw(s)*phi % y(c1) + (1.0-Grid % fw(s))*phi % y(c2)
      phiz_f1 = Grid % fw(s)*phi % z(c1) + (1.0-Grid % fw(s))*phi % z(c2)
      phix_f2 = phix_f1
      phiy_f2 = phiy_f1
      phiz_f2 = phiz_f1
      dif_eff1 = Grid % f(s) *(Flow % diffusivity)  &
           + (1.-Grid % f(s))*(Flow % diffusivity)
      if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
         turb % model .ne. DNS) then
        dif_eff1 = Grid % f(s) *(Flow % diffusivity + turb % vis_t(c1)/sc_t)  &
             + (1.-Grid % f(s))*(Flow % diffusivity + turb % vis_t(c2)/sc_t)
        if(turb % model .eq. HYBRID_LES_RANS) then
          dif_eff1 = Grid % f(s)   &
                  * (Flow % diffusivity + turb % vis_t_eff(c1) / sc_t)  &
               + (1.-Grid % f(s))  &
                  * (Flow % diffusivity + turb % vis_t_eff(c2) / sc_t)
        end if
      end if
      dif_eff2 = dif_eff1
    else
      phix_f1 = phi % x(c1)
      phiy_f1 = phi % y(c1)
      phiz_f1 = phi % z(c1)
      phix_f2 = phix_f1
      phiy_f2 = phiy_f1
      phiz_f2 = phiz_f1
      dif_eff1 = Flow % diffusivity
      if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
         turb % model .ne. DNS) then
        dif_eff1 = Flow % diffusivity + turb % vis_t(c1) / sc_t
        if(turb % model .eq. HYBRID_LES_RANS) then
          dif_eff1 = Flow % diffusivity + turb % vis_t_eff(c1) / sc_t
        end if
      end if
      dif_eff2 = dif_eff1
    end if

    ! Wall diffusivity for species
    if(turb % model .eq. K_EPS .or.        &
       turb % model .eq. K_EPS_ZETA_F .or. &
       turb % model .eq. HYBRID_LES_RANS) then
      if(c2 < 0) then
        if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALL .or.  &
           Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALLFL) then
          dif_eff1 = turb % diff_w(c1)
          dif_eff2 = dif_eff1
        end if
      end if
    end if

    ! Total (exact) diffusive flux
    f_ex1 = dif_eff1 * (  phix_f1 * Grid % sx(s)  &
                        + phiy_f1 * Grid % sy(s)  &
                        + phiz_f1 * Grid % sz(s))
    f_ex2 = dif_eff2 * (  phix_f2 * Grid % sx(s)  &
                        + phiy_f2 * Grid % sy(s)  &
                        + phiz_f2 * Grid % sz(s))

    ! Implicit diffusive flux
    f_im1 = dif_eff1 * A % fc(s)           &
          * (  phix_f1 * Grid % dx(s)      &
             + phiy_f1 * Grid % dy(s)      &
             + phiz_f1 * Grid % dz(s) )
    f_im2 = dif_eff2 * A % fc(s)           &
          * (  phix_f2 * Grid % dx(s)      &
             + phiy_f2 * Grid % dy(s)      &
             + phiz_f2 * Grid % dz(s) )

    ! Cross diffusion part
    phi % c(c1) = phi % c(c1) + f_ex1 - f_im1
    if(c2 .gt. 0) then
      phi % c(c2) = phi % c(c2) - f_ex2 + f_im2
    end if

    ! Calculate the coefficients for the sysytem matrix

    a12 = dif_eff1 * A % fc(s)
    a21 = dif_eff2 * A % fc(s)

    a12 = a12  - min(v_flux % n(s), 0.0) * Flow % density(c1)
    a21 = a21  + max(v_flux % n(s), 0.0) * Flow % density(c2)

    ! Fill the system matrix
    if(c2 > 0) then
      A % val(A % dia(c1))  = A % val(A % dia(c1)) + a12
      A % val(A % dia(c2))  = A % val(A % dia(c2)) + a21
      A % val(A % pos(1,s)) = A % val(A % pos(1,s)) - a12
      A % val(A % pos(2,s)) = A % val(A % pos(2,s)) - a21
    else if(c2 < 0) then

      ! Outflow is included because of the flux
      ! corrections which also affects velocities
      if( (Var_Mod_Bnd_Cond_Type(phi,c2) .eq. INFLOW) .or.  &
          (Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALL)   .or.  &
          (Var_Mod_Bnd_Cond_Type(phi,c2) .eq. CONVECT) ) then
        A % val(A % dia(c1)) = A % val(A % dia(c1)) + a12
        b(c1)  = b(c1)  + a12 * phi % n(c2)

      ! In case of wallflux 
      else if(Var_Mod_Bnd_Cond_Type(phi,c2) .eq. WALLFL) then
        b(c1) = b(c1) + Grid % s(s) * phi % q(c2)
      end if

    end if

  end do  ! through sides

  ! Implicit treatment for cross difusive terms
  do c = 1, Grid % n_cells
    if(phi % c(c) >= 0) then
      b(c)  = b(c) + phi % c(c)
    else
      A % val(A % dia(c)) = A % val(A % dia(c))  &
                          - phi % c(c) / (phi % n(c) + MICRO)
    end if
  end do

  !--------------------!
  !                    !
  !   Inertial terms   !
  !                    !
  !--------------------!
  call Numerics_Mod_Inertial_Term(phi, Flow % density, A, b, dt)

  !-------------------------------------!
  !                                     !
  !   Source terms and wall function    !
  !                                     !
  !-------------------------------------!
  if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    if(turb % model_variant .ne. STABILIZED) then
      do c = 1, Grid % n_cells
        u1uj_phij(c) = -0.22 * turb % t_scale(c) *  &
                   (  uu % n(c) * phi % x(c)          &
                    + uv % n(c) * phi % y(c)          &
                    + uw % n(c) * phi % z(c))
        u2uj_phij(c) = -0.22 * turb % t_scale(c) *  &
                   (  uv % n(c) * phi % x(c)          &
                    + vv % n(c) * phi % y(c)          &
                    + vw % n(c) * phi % z(c))
        u3uj_phij(c) = -0.22 * turb % t_scale(c) *  &
                   (  uw % n(c) * phi % x(c)          &
                    + vw % n(c) * phi % y(c)          &
                    + ww % n(c) * phi % z(c))
      end do
      call Flow % Grad_Component(u1uj_phij, 1, u1uj_phij_x)
      call Flow % Grad_Component(u2uj_phij, 2, u2uj_phij_y)
      call Flow % Grad_Component(u3uj_phij, 3, u3uj_phij_z)
      do c = 1, Grid % n_cells
        b(c) = b(c) - (  u1uj_phij_x(c)  &
                       + u2uj_phij_y(c)  &
                       + u3uj_phij_z(c) ) * Grid % vol(c)
      end do

      !------------------------------------------------------------------!
      !   Here we clean up transport equation from the false diffusion   !
      !------------------------------------------------------------------!
      do s = 1, Grid % n_faces

        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

!   For species, we don't have a function for turbulent Schmidt number
!       sc_t1 = Turb_Mod_Prandtl_Number(turb, c1)
!       sc_t2 = Turb_Mod_Prandtl_Number(turb, c2)
!       sc_t  = Grid % fw(s) * sc_t1 + (1.0-Grid % fw(s)) * sc_t2

        if(c2 > 0) then
          phix_f1 = Grid % fw(s)*phi % x(c1) + (1.0-Grid % fw(s))*phi % x(c2)
          phiy_f1 = Grid % fw(s)*phi % y(c1) + (1.0-Grid % fw(s))*phi % y(c2)
          phiz_f1 = Grid % fw(s)*phi % z(c1) + (1.0-Grid % fw(s))*phi % z(c2)
          phix_f2 = phix_f1 
          phiy_f2 = phiy_f1 
          phiz_f2 = phiz_f1 
          dif_eff1 =      Grid % f(s)  * (turb % vis_t(c1)/sc_t )  &
                  + (1. - Grid % f(s)) * (turb % vis_t(c2)/sc_t )
          dif_eff2 = dif_eff1 
        else
          phix_f1 = phi % x(c1)
          phiy_f1 = phi % y(c1)
          phiz_f1 = phi % z(c1)
          phix_f2 = phix_f1
          phiy_f2 = phiy_f1
          phiz_f2 = phiz_f1
          dif_eff1 = turb % vis_t(c1) / sc_t
          dif_eff2 = dif_eff1
        end if

        ! Total (exact) diffusive flux
        f_ex1 = dif_eff1 * (  phix_f1 * Grid % sx(s)  &
                            + phiy_f1 * Grid % sy(s)  &
                            + phiz_f1 * Grid % sz(s))
        f_ex2 = dif_eff2 * (  phix_f2 * Grid % sx(s)  &
                            + phiy_f2 * Grid % sy(s)  &
                            + phiz_f2 * Grid % sz(s))

        ! Implicit diffusive flux
        f_im1 = dif_eff1 * A % fc(s) *         &
                (  phix_f1 * Grid % dx(s)      &
                 + phiy_f1 * Grid % dy(s)      &
                 + phiz_f1 * Grid % dz(s) )
        f_im2 = dif_eff2 * A % fc(s) *         &
                (  phix_f2 * Grid % dx(s)      &
                 + phiy_f2 * Grid % dy(s)      &
                 + phiz_f2 * Grid % dz(s) )

        b(c1) = b(c1) - dif_eff1 * (phi % n(c2) - phi % n(c1)) * A % fc(s)  &
              - f_ex1 + f_im1
        if(c2  > 0) then
          b(c2) = b(c2) + dif_eff1 * (phi % n(c2) - phi % n(c1)) * A % fc(s)  &
                + f_ex2 - f_im2
        end if
      end do
    end if
  end if

  call User_Mod_Source(Flow, phi, A, b)

  !---------------------------------!
  !                                 !
  !   Solve the equations for phi   !
  !                                 !
  !---------------------------------!

  ! Under-relax the equations
  call Numerics_Mod_Under_Relax(phi, A, b)

  ! Call linear solver to solve them
  call Cpu_Timer % Start('Linear_Solver_For_Scalars')
  call Sol % Bicg(A,              &
                  phi % n,        &
                  b,              &
                  phi % precond,  &
                  phi % mniter,   &
                  phi % eniter,   &
                  phi % tol,      &
                  phi % res)
  call Cpu_Timer % Stop('Linear_Solver_For_Scalars')

  read(phi % name(3:4), *) ns  ! reterive the number of scalar
  row = ceiling(ns/6)          ! will be 1 (scal. 1-6), 2 (scal. 6-12), etc.
  col = nint(ns) - (row-1)*6   ! will be in range 1 - 6

  call Info_Mod_Iter_Fill_User_At(row, col, phi % name, phi % eniter, phi % res)

  call Flow % Grad_Variable(phi)

  ! User function
  call User_Mod_End_Of_Compute_Scalar(Flow, turb, Vof, Sol, curr_dt, ini, sc)

  call Cpu_Timer % Stop('Compute_Scalars (without solvers)')

  end subroutine
