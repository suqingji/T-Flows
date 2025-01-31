!==============================================================================!
  subroutine Compute_Momentum(Flow, turb, Vof, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Discretizes and solves momentum conservation equations                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt
  integer, intent(in)         :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Bulk_Type),   pointer :: bulk
  type(Matrix_Type), pointer :: M
  type(Var_Type),    pointer :: ui, uj, uk, t, p
  type(Face_Type),   pointer :: v_flux
  real, contiguous,  pointer :: b(:)
  real, contiguous,  pointer :: ui_i(:), ui_j(:), ui_k(:), uj_i(:), uk_i(:)
  real, contiguous,  pointer :: si(:), sj(:), sk(:), di(:), dj(:), dk(:)
  real, contiguous,  pointer :: fi(:), p_i(:), cell_fi(:), st_i(:)
  integer                    :: s, c, c1, c2, i
  real                       :: f_ex, f_im, f_stress
  real                       :: vel_max, dt
  real                       :: m0, m12, m21
  real                       :: vis_eff
  real                       :: ui_i_f, ui_j_f, ui_k_f, uj_i_f, uk_i_f
  real                       :: grav_i, p_drop_i
  real                       :: ui_si, ui_di
!------------------------------------------------------------------------------!
!                                                                              !
!  Stress tensor on the face s:                                                !
!                                                                              !
!    t = mu * [    2*du/dx     du/dy+dv/dx   du/dz+dw/dx  ]                    !
!             [  du/dy+dv/dx     2*dv/dy     dv/dz+dw/dy  ]                    !
!             [  du/dz+dw/dx   dv/dz+dw/dy     2*dw/dz    ]                    !
!                                                                              !
!  The forces, acting on the cell face are:                                    !
!                                                                              !
!    fx = t11*sx + t12*sy + t13*sz                                             !
!    fy = t21*sx + t22*sy + t23*sz                                             !
!    fz = t31*sx + t32*sy + t33*sz                                             !
!                                                                              !
!  which could also be written in the compact form:                            !
!                                                                              !
!    {f} = [t]{s}                                                              !
!                                                                              !
!  or in expended form:                                                        !
!                                                                              !
!    fx = txx*sx + txy*sy + txz*sz                                             !
!    fy = tyx*sx + tyy*sy + tyz*sz                                             !
!    fz = tzx*sx + tzy*sy + tzz*sz                                             !
!                                                                              !
!------------------------------------------------------------------------------!
!                                                                              !
!  The form of equations which I am solving:                                   !
!                                                                              !
!     /             /              /               /             /             !
!    |     du      |              |               |             |              !
!    | rho -- dV + | rho u u dS = | mu DIV u dS - | GRAD p dV + | f dV         !
!    |     dt      |              |               |             |              !
!   /             /              /               /             /               !
!                                                                              !
!  Dimension of the system under consideration                                 !
!                                                                              !
!     [M]{u} = {b}   [kgm/s^2]   [N]                                           !
!                                                                              !
!  Dimensions of certain variables:                                            !
!                                                                              !
!     M              [kg/s]                                                    !
!     u, v, w        [m/s]                                                     !
!     bu, bv, bw     [kgm/s^2]      [N]                                        !
!     p, pp          [kg/(m s^2)]   [N/m^2]                                    !
!     v_flux         [m^3/s]                                                   !
!     au*, av*, aw*  [kgm/s^2]      [N]                                        !
!     du*, dv*, dw*  [kgm/s^2]      [N]                                        !
!     cu*, cv*, cw*  [kgm/s^2]      [N]                                        !
!==============================================================================!

  call Cpu_Timer % Start('Compute_Momentum (without solvers)')

  ! Take aliases
  Grid   => Flow % pnt_grid
  bulk   => Flow % bulk
  v_flux => Flow % v_flux
  t      => Flow % t
  p      => Flow % p
  dt     =  Flow % dt
  M      => Sol % M
  b      => Sol % b % val

  ! User function
  call User_Mod_Beginning_Of_Compute_Momentum(Flow, turb, Vof, Sol,  &
                                              curr_dt, ini)

  !-------------------------------------------------------!
  !   Store the old volume fluxes for Choi's correction   !
  !-------------------------------------------------------!
  if (Flow % piso_status .eqv. .false.) then  ! check about this
    if(ini .eq. 1) then
      do s = 1, Grid % n_faces
        v_flux % oo(s) = v_flux % o(s)
        v_flux % o (s) = v_flux % n(s)
      end do
    end if
  end if

  !--------------------------------------------!
  !                                            !
  !                                            !
  !   Browse through all velocity components   !
  !                                            !
  !                                            !
  !--------------------------------------------!
  do i = 1, 3

    if(i .eq. 1) then
      ui   => Flow % u;   uj   => Flow % v;   uk   => Flow % w
      ui_i => ui % x;     ui_j => ui % y;     ui_k => ui % z
      si   => Grid % sx;  sj   => Grid % sy;  sk   => Grid % sz
      di   => Grid % dx;  dj   => Grid % dy;  dk   => Grid % dz
      p_i  => p % x;      uj_i => uj % x;     uk_i => uk % x
      fi       => Flow % fx
      cell_fi  => Flow % cell_fx
      grav_i   =  Flow % grav_x
      p_drop_i =  bulk % p_drop_x
      st_i     => Vof % surf_fx
    end if
    if(i .eq. 2) then
      ui   => Flow % v;   uj   => Flow % w;   uk   => Flow % u
      ui_i => ui % y;     ui_j => ui % z;     ui_k => ui % x
      si   => Grid % sy;  sj   => Grid % sz;  sk   => Grid % sx
      di   => Grid % dy;  dj   => Grid % dz;  dk   => Grid % dx
      p_i  => p % y;      uj_i => uj % y;     uk_i => uk % y
      fi       => Flow % fy
      cell_fi  => Flow % cell_fy
      grav_i   =  Flow % grav_y
      p_drop_i =  bulk % p_drop_y
      st_i     => Vof % surf_fy
    end if
    if(i .eq. 3) then
      ui   => Flow % w;   uj   => Flow % u;   uk   => Flow % v
      ui_i => ui % z;     ui_j => ui % x;     ui_k => ui % y
      si   => Grid % sz;  sj   => Grid % sx;  sk   => Grid % sy
      di   => Grid % dz;  dj   => Grid % dx;  dk   => Grid % dy
      p_i  => p % z;      uj_i => uj % z;     uk_i => uk % z
      fi       => Flow % fz
      cell_fi  => Flow % cell_fz
      grav_i   =  Flow % grav_z
      p_drop_i =  bulk % p_drop_z
      st_i     => Vof % surf_fz
    end if

    ! Initialize advection, cross diffusion, forces, matrix and right hand side
    ui % a (:) = 0.0
    ui % c (:) = 0.0
    fi     (:) = 0.0  ! all "internal" forces acting on this component
    f_stress   = 0.0  ! this is presumably not needed
    M % val(:) = 0.0
    b      (:) = 0.0

    ! Calculate velocity magnitude for normalization
    vel_max = MICRO
    do c = -Grid % n_bnd_cells, Grid % n_cells
      vel_max = max(vel_max, sqrt(ui % n(c)**2 + uj % n(c)**2 + uk % n(c)**2))
    end do
    call Comm_Mod_Global_Max_Real(vel_max)

    ! Old values (o) and older than old (oo)
    if (Flow % piso_status .eqv. .false.) then
      if(ini .eq. 1) then
        do c = 1, Grid % n_cells
          ui % oo(c) = ui % o(c)
          ui % o (c) = ui % n(c)
        end do
      end if
    end if

    !--------------------------------------------------------!
    !   Compute buoyancy force for this velocity component   !
    !--------------------------------------------------------!
    call Flow % Buoyancy_Forces(i)

    !---------------!
    !               !
    !   Advection   !
    !               !
    !---------------!
    call Numerics_Mod_Advection_Term(ui, Flow % density, v_flux % n, fi)

    !---------------!
    !               !
    !   Diffusion   !
    !               !
    !---------------!

    !----------------------------!
    !   Spatial discretization   !
    !----------------------------!
    do s = 1, Grid % n_faces

      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      call Turb_Mod_Face_Vis   (turb, vis_eff,  s)
      call Turb_Mod_Face_Stress(turb, ui, f_stress, s)

      ui_i_f = Grid % fw(s)*ui_i(c1) + (1.0-Grid % fw(s))*ui_i(c2)
      ui_j_f = Grid % fw(s)*ui_j(c1) + (1.0-Grid % fw(s))*ui_j(c2)
      ui_k_f = Grid % fw(s)*ui_k(c1) + (1.0-Grid % fw(s))*ui_k(c2)
      uj_i_f = Grid % fw(s)*uj_i(c1) + (1.0-Grid % fw(s))*uj_i(c2)
      uk_i_f = Grid % fw(s)*uk_i(c1) + (1.0-Grid % fw(s))*uk_i(c2)

      ui_si = (  (ui_i_f + ui_i_f) * si(s)    &
               + (ui_j_f + uj_i_f) * sj(s)    &
               + (ui_k_f + uk_i_f) * sk(s) )
      ui_di = (  ui_i_f * di(s)  &
               + ui_j_f * dj(s)  &
               + ui_k_f * dk(s))

      ! Total (exact) viscous stress
      f_ex = vis_eff * ui_si

      ! Implicit viscous stress
      m0 = vis_eff * M % fc(s)
      f_im = ui_di * m0

      ! Cross diffusion part
      ui % c(c1) = ui % c(c1) + f_ex - f_im + f_stress * Flow % density(c1)
      if(c2  > 0) then
        ui % c(c2) = ui % c(c2) - f_ex + f_im - f_stress * Flow % density(c2)
      end if

      ! Compute the coefficients for the sysytem matrix
      m12 = m0 - min(v_flux % n(s), 0.0) * Flow % density(c1)
      m21 = m0 + max(v_flux % n(s), 0.0) * Flow % density(c2)

      ! Fill the system matrix
      if(c2 > 0) then
        M % val(M % pos(1,s)) = M % val(M % pos(1,s)) - m12
        M % val(M % dia(c1))  = M % val(M % dia(c1))  + m12
        M % val(M % pos(2,s)) = M % val(M % pos(2,s)) - m21
        M % val(M % dia(c2))  = M % val(M % dia(c2))  + m21
      else if(c2  < 0) then
        ! Outflow is not included because it was causing problems
        if((Grid % Bnd_Cond_Type(c2) .eq. INFLOW)  .or.  &
           (Grid % Bnd_Cond_Type(c2) .eq. WALL)    .or.  &
           (Grid % Bnd_Cond_Type(c2) .eq. CONVECT) .or.  &
           (Grid % Bnd_Cond_Type(c2) .eq. WALLFL)) then
           ! (Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW) ) then
          M % val(M % dia(c1)) = M % val(M % dia(c1)) + m12
          fi(c1) = fi(c1) + m12 * ui % n(c2)
        end if
      end if

      ! Here we clean up momentum from the false diffusion
      call Turb_Mod_Substract_Face_Stress(turb, ui_si, ui_di,            &
                                                ui % n(c1), ui % n(c2),  &
                                                M % fc(s), fi, s)

    end do  ! through faces

    ! Explicit treatment for cross diffusion terms
    ! (Shouldn't theese, in an ideal world,
    !  also be treated in Rhie and Chow?)
    do c = 1, Grid % n_cells
      fi(c) = fi(c) + ui % c(c)
    end do

    !--------------------!
    !                    !
    !   Inertial terms   !
    !                    !
    !--------------------!
    call Numerics_Mod_Inertial_Term(ui, Flow % density, M, fi, dt)

    !---------------------------------!
    !                                 !
    !   Various force contributions   !
    !                                 !
    !---------------------------------!

    !--------------------------!
    !   Global pressure drop   !
    !--------------------------!
    do c = 1, Grid % n_cells
      fi(c) = fi(c) + p_drop_i * Grid % vol(c)
    end do

    !--------------------!
    !   Buoyancy force   !
    !--------------------!
    do c = 1, Grid % n_cells
      fi(c) = fi(c) + cell_fi(c) * Grid % vol(c)
    end do

    !----------------------------------------!
    !   All other terms defined by the user  !
    !----------------------------------------!
    call User_Mod_Force(Flow, ui, M, fi)

    !-----------------------------------------------------------!
    !   Copy forces from current component to right hand side   !
    !   (Note: pressure gradients are not with other forces.    !
    !    Same is true for surface tension, see nex comments)    !
    !-----------------------------------------------------------!
    do c = 1, Grid % n_cells
      b(c) = fi(c) - p_i(c) * Grid % vol(c)
    end do

    !----------------------------------------------------------------!
    !   In case of vof simulations, add the surface tension forces   !
    !   (Note: they are treated like pressure everywhere in the      !
    !    discretized form, meaning separater from other forces)      !
    !----------------------------------------------------------------!
    if(Flow % with_interface) then
      call Vof % Surface_Tension_Force(i)
      do c = 1, Grid % n_cells
        b(c) = b(c) + st_i(c) * Grid % vol(c)
      end do
    end if

    !------------------------------------------------!
    !   Save the coefficients from the discretized   !
    !   momentum equation before under-relaxation    !
    !                                                !
    !   If you save them like this, before the un-   !
    !   der-relaxation, results are independent      !
    !   from under-relaxation factors and Majumdar   !
    !   correction in Rhie_And_Chow is not needed    !
    !------------------------------------------------!
    do c = 1, Grid % n_cells
      M % sav(c) = M % val(M % dia(c))
    end do

    !----------------------------------------------!
    !   Explicit solution for the PISO algorithm   !
    !----------------------------------------------!
    call Compute_Momentum_Explicit(Flow, ui, Sol)

    !-----------------------------------!
    !                                   !
    !   Solve the equations for u,v,w   !
    !                                   !
    !-----------------------------------!

    !--------------------------------------------------------!
    !   If not inside the PRIME part of the PISO algorithm   !
    !--------------------------------------------------------!
    if(Flow % piso_status .eqv. .false.) then

      ! Under-relax the equations
      call Numerics_Mod_Under_Relax(ui, M, b)

      ! Call linear solver
      call Cpu_Timer % Start('Linear_Solver_For_Momentum')

      call Sol % Bicg(M,             &
                      ui % n,        &
                      b,             &
                      ui % precond,  &
                      ui % mniter,   &
                      ui % eniter,   &
                      ui % tol,      &
                      ui % res,      &
                      norm = vel_max)
      call Cpu_Timer % Stop('Linear_Solver_For_Momentum')

      ! Fill the info screen up
      if (Flow % p_m_coupling == SIMPLE) then
        call Info_Mod_Iter_Fill_At(1, i, ui % name, ui % eniter, ui % res)
      end if

    end if

  end do  ! browsing through components

  ! Refresh buffers for M % sav before discretizing for pressure
  call Grid % Exchange_Cells_Real(M % sav)

  ! User function
  call User_Mod_End_Of_Compute_Momentum(Flow, turb, Vof, Sol, curr_dt, ini)

  call Cpu_Timer % Stop('Compute_Momentum (without solvers)')

  end subroutine
