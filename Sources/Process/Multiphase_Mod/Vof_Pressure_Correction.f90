!==============================================================================!
  subroutine Multiphase_Mod_Vof_Pressure_Correction(mult, sol, ini, mass_err)
!------------------------------------------------------------------------------!
!   Correct fluxes on pressure equation due to surface tension, temporal terms,!
!   under relaxation and skewness                                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: curr_colour => r_cell_01,  &
                      go_x        => r_cell_02,  &
                      go_y        => r_cell_03,  &
                      go_z        => r_cell_04
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer                       :: ini
  real                          :: mass_err
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Face_Type),   pointer :: v_flux
  type(Var_Type),    pointer :: vof
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: new_var
  type(Matrix_Type), pointer :: a
  real, contiguous,  pointer :: b(:), old_var(:), si(:)
  real,              pointer :: u_relax, dt_corr
  real                       :: gf_x, gf_y, gf_z
  integer                    :: i_dir
  integer                    :: c, c1, c2, s, nb, nc
  real                       :: a12, fs
  real                       :: u_fo, v_fo, w_fo, tf
  real                       :: stens_source, dotprod
  real                       :: factor2, correction, dens_h, curv_f
!==============================================================================!

  ! Take aliases
  grid    => mult % pnt_grid
  flow    => mult % pnt_flow
  vof     => mult % vof
  u_relax => flow % u_rel_corr
  dt_corr => flow % dt_corr
  v_flux  => flow % v_flux
  a       => sol % a
  b       => sol % b % val

  nb = grid % n_bnd_cells
  nc = grid % n_cells

  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! Correct for Surface tension
  if(mult % surface_tension > TINY) then

    curr_colour(-nb:nc) = vof % n(-nb:nc)

    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      ! Interpolate VOF gradients
      dens_h = 2.0 / ( 1.0 / flow % density(c1) + 1.0 / flow % density(c2) )

      ! Unit for dotprod: [1/m]
      dotprod = 0.5 * dens_h                                                 &
                    * ( mult % curv(c1) * vof % x(c1) / flow % density(c1)   &
                      + mult % curv(c2) * vof % x(c2) / flow % density(c2) ) &
                    * grid % dx(s)                                           &
              + 0.5 * dens_h                                                 &
                    * ( mult % curv(c1) * vof % y(c1) / flow % density(c1)   &
                      + mult % curv(c2) * vof % y(c2) / flow % density(c2) ) &
                    * grid % dy(s)                                           &
              + 0.5 * dens_h                                                 &
                    * ( mult % curv(c1) * vof % z(c1) / flow % density(c1)   &
                      + mult % curv(c2) * vof % z(c2) / flow % density(c2) ) &
                    * grid % dz(s)

      ! Unit for a12: [m^4s/kg]
      a12 = u_relax * 0.5 * ( grid % vol(c1) / a % sav(c1)     &
                            + grid % vol(c2) / a % sav(c2) ) * a % fc(s)

      ! Curvature at the face; unit: [1/m]
      curv_f = 0.5 * ( mult % curv(c1) + mult % curv(c2) )

      ! Unit for stens_source: [kg/s^2 * m^4s/kg * 1/m = m^3/s]
      stens_source = mult % surface_tension * a12               &
                   * ( curv_f * (curr_colour(c2) -  curr_colour(c1)) - dotprod )

      v_flux % n(s) = v_flux % n(s) + stens_source

      b(c1) = b(c1) - stens_source
      b(c2) = b(c2) + stens_source

    end do

  end if

  ! Introduce temporal correction and subrelaxation
  ! (See equation 3.61 in Denner's thesis)
  if (flow % temp_corr) then
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)
      u_fo = fs * u % o(c1) + (1.0 - fs) * u % o(c2)
      v_fo = fs * v % o(c1) + (1.0 - fs) * v % o(c2)
      w_fo = fs * w % o(c1) + (1.0 - fs) * w % o(c2)

      dens_h = 2.0 / ( 1.0 / flow % density(c1) + 1.0 / flow % density(c2) )

      ! Unit for factor2: [1]
      factor2 = u_relax * 0.5 * dens_h * ( grid % vol(c1)               &
                                         / (a % sav(c1) * dt_corr)      &
                                         + grid % vol(c2)               &
                                         / (a % sav(c2) * dt_corr) )

      ! Unit for correction [m^3/s]
      correction = (1.0 - u_relax) * ( v_flux % star(s) - v_flux % avg(s) )   &
                 + factor2 * ( v_flux % o(s) - ( u_fo * grid % sx(s)          &
                                               + v_fo * grid % sy(s)          &
                                               + w_fo * grid % sz(s) ) )

      v_flux % n(s) = v_flux % n(s) + correction

      b(c1) = b(c1) - correction
      b(c2) = b(c2) + correction

    end do
  end if

  if (mult % skew_corr ) then
    ! loop on each direction

    do i_dir = 1, 3
      select case(i_dir)
      case(1)
        old_var => u % o
        new_var => u
        si      => grid % sx
      case(2)
        old_var => v % o
        new_var => v
        si      => grid % sy
      case(3)
        old_var => w % o
        new_var => w
        si      => grid % sz
      end select
      call Grid_Mod_Exchange_Cells_Real(grid, old_var)

      call Field_Mod_Grad(flow, old_var(-nb:nc), go_x(-nb:nc),  &
                                                 go_y(-nb:nc),  &
                                                 go_z(-nb:nc))
      do s = grid % n_bnd_faces + 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        fs = grid % f(s)

        ! Add new velocity component correction
        gf_x = fs * new_var % x(c1) + (1.0 - fs) * new_var % x(c2)
        gf_y = fs * new_var % y(c1) + (1.0 - fs) * new_var % y(c2)
        gf_z = fs * new_var % z(c1) + (1.0 - fs) * new_var % z(c2)

        correction = dot_product( (/gf_x, gf_y, gf_z/),                     &
                     (/grid % xr(s), grid % yr(s), grid % zr(s)/)) * si(s)

        ! Add old velocity component correction
        gf_x = fs * go_x(c1) + (1.0 - fs) * go_x(c2)
        gf_y = fs * go_y(c1) + (1.0 - fs) * go_y(c2)
        gf_z = fs * go_z(c1) + (1.0 - fs) * go_z(c2)

        dens_h = 2.0 / ( 1.0 / flow % density(c1) + 1.0 / flow % density(c2) )

        ! Unit for factor2: [1]
        factor2 = u_relax * 0.5 * dens_h * ( grid % vol(c1)               &
                                           / (a % sav(c1) * dt_corr)      &
                                           + grid % vol(c2)               &
                                           / (a % sav(c2) * dt_corr) )

        correction = correction                                             &
                   - dot_product( (/gf_x, gf_y, gf_z/),                     &
                     (/grid % xr(s), grid % yr(s), grid % zr(s)/)) * si(s)  &
                   * factor2

        v_flux % n(s) = v_flux % n(s) + correction

        b(c1) = b(c1) - correction
        b(c2) = b(c2) + correction

      end do
    end do
  end if

  if (mult % phase_change) then
    do c = 1, grid % n_cells
      b(c) = b(c) + mult % flux_rate(c) * grid % vol(c)                    &
                                        * ( 1.0 / mult % phase_dens(1)     &
                                          - 1.0 / mult % phase_dens(2) )
    end do
  end if

  end subroutine
