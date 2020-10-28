!============================================================================!
  subroutine Multiphase_Mod_Vof_Momentum_Contribution(mult, sol, ui, i)
!----------------------------------------------------------------------------!
!   Computes Surface tension, Gravity and phase change sources for Momentum  !
!   Equation if a two-phase flow calculation is performed. Additionally and  !
!   for the moment, PISO calculations are run here                           !
!----------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]--------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  type(Var_Type),        target :: ui
  integer                       :: i
!-----------------------------------[Locals]---------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: vof, avg_var
  type(Face_Type),   pointer :: v_flux
  type(Matrix_Type), pointer :: a
  real, contiguous,  pointer :: b(:)
  real, contiguous,  pointer :: surf_fx(:), surf_fy(:), surf_fz(:)
  real, contiguous,  pointer :: si(:)
  real,              pointer :: u_relax, dt_corr
  integer                    :: s, c, c1, c2, nt, ni, i_dir
  real                       :: gf_x, gf_y, gf_z, correction
  real                       :: dotprod, epsloc, fs
  real                       :: corr_x, corr_y, corr_z
  real                       :: u_f, v_f, w_f
!============================================================================!

  ! Take aliases
  flow    => mult % pnt_flow
  grid    => mult % pnt_grid
  vof     => mult % vof
  surf_fx => mult % surf_fx
  surf_fy => mult % surf_fy
  surf_fz => mult % surf_fz
  v_flux  => flow % v_flux
  a       => sol % a
  b       => sol % b % val
  u_relax => flow % u_rel_corr

  epsloc = epsilon(epsloc)

  ! Surface tension contribution
  if (mult % surface_tension > TINY) then

    select case(i)
      case(1)
        do c = 1, grid % n_cells
          surf_fx(c) = mult % surface_tension  &
                     * mult % curv(c)          &
                     * vof % x(c)              &
                     * grid % vol(c)
          b(c) = b(c) + surf_fx(c)
         end do
      case(2)
        do c = 1, grid % n_cells
          surf_fy(c) = mult % surface_tension  &
                     * mult % curv(c)          &
                     * vof % y(c)              &
                     * grid % vol(c)
          b(c) = b(c) + surf_fy(c)
        end do
      case(3)
        do c = 1, grid % n_cells
          surf_fz(c) = mult % surface_tension  &
                     * mult % curv(c)          &
                     * vof % z(c)              &
                     * grid % vol(c)
          b(c) = b(c) + surf_fz(c)
        end do

    end select

  end if

  ! Momentum variables for pressure correction
  ! This is here because they need to be collected before
  ! u, v, w are calculated

  if (flow % temp_corr) then
    ! Guessed face velocity
    if (i == 1) then
      do s = grid % n_bnd_faces + 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        fs = grid % f(s)
        u_f = fs * flow % u % n(c1) + (1.0 - fs) * flow % u % n(c2)
        v_f = fs * flow % v % n(c1) + (1.0 - fs) * flow % v % n(c2)
        w_f = fs * flow % w % n(c1) + (1.0 - fs) * flow % w % n(c2)
        v_flux % avg(s) = ( u_f * grid % sx(s)     &
                          + v_f * grid % sy(s)     &
                          + w_f * grid % sz(s) )
      end do

      do s = grid % n_bnd_faces + 1, grid % n_faces
        v_flux % star(s) = v_flux % n(s)
      end do

      if (mult % skew_corr ) then
        do i_dir = 1, 3
          select case(i_dir)
          case(1)
            avg_var => flow % u
            si      => grid % sx
          case(2)
            avg_var => flow % v
            si      => grid % sy
          case(3)
            avg_var => flow % w
            si      => grid % sz
          end select

          call Field_Mod_Grad_Variable(flow, avg_var)

          do s = grid % n_bnd_faces + 1, grid % n_faces
            c1 = grid % faces_c(1,s)
            c2 = grid % faces_c(2,s)
            fs = grid % f(s)

            ! Add new velocity component correction
            gf_x = fs * avg_var % x(c1) + (1.0 - fs) * avg_var % x(c2)
            gf_y = fs * avg_var % y(c1) + (1.0 - fs) * avg_var % y(c2)
            gf_z = fs * avg_var % z(c1) + (1.0 - fs) * avg_var % z(c2)

            correction = dot_product( (/gf_x, gf_y, gf_z/),                    &
                         (/grid % xr(s), grid % yr(s), grid % zr(s)/)) * si(s)

            v_flux % avg(s) = v_flux % avg(s) + correction

          end do
        end do
      end if
    end if
  end if

  end subroutine
