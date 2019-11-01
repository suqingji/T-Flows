!==============================================================================!
  subroutine Source_T2(flow, sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in t2 transport equation for k-eps_t2 model      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grad_Mod
  use Grid_Mod,   only: Grid_Type
  use Solver_Mod, only: Solver_Type
  use Matrix_Mod, only: Matrix_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
  real :: Roughness_Coefficient
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: EBF, kin_vis, p_t2_wall
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  a        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
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

  call Grad_Mod_Array(grid, t % n, t % x, t % y, t % z, .true.)

  !-----------------------------------------!
  !   Compute the sources in all the cells  !
  !-----------------------------------------!

  ! Production source:
  do c = 1, grid % n_cells

    p_t2(c) = - 2.0 * (  ut % n(c) * t % x(c)   &
                       + vt % n(c) * t % y(c)   &
                       + wt % n(c) * t % z(c))

    b(c) = b(c) + p_t2(c) * grid % vol(c)

   ! Negative contribution
   a % val(a % dia(c)) = a % val(a % dia(c)) +  &
         2.0 * density * eps % n(c) / (kin % n(c) + TINY) * grid % vol(c)

  end do

  ! Kinematic viscosities
  kin_vis = viscosity / density

  ! Implementation of wall function approach for buoyancy-driven flows

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or. &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        y_plus(c1) = Y_Plus_Low_Re(u_tau(c1),           &
                     grid % wall_dist(c1), kin_vis)

        EBF  = 0.01*y_plus(c1)**4.0/(1.0+5.0*y_plus(c1))

        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) & 
        t % q(c2) = abs(con_wall(c1)*(t % n(c1) &
                    - t % n(c2))/grid % wall_dist(c1))

        p_t2_wall  = t % q(c2)*c_mu_theta5*sqrt(abs(t2 % n(c1))) &
                     /(kappa_theta*c_mu25*grid % wall_dist(c1))

        b(c1) = b(c1) - p_t2(c1) * grid % vol(c1)

        if(y_plus(c1) > 11.0) then
          b(c1) = b(c1) + p_t2_wall * grid % vol(c1)
        else  
          b(c1) = b(c1) + (p_t2(c1) * exp(-1.0 * EBF) + &
                  p_t2_wall * exp(-1.0/EBF)) * grid % vol(c1)
        end if

        t2 % n(c2) = 0.0
         
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  end subroutine

