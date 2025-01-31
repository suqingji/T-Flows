!==============================================================================!
  subroutine Turb_Mod_Calculate_Mean(turb, n0, n1)
!------------------------------------------------------------------------------!
!   Calculates time averaged velocity and velocity fluctuations.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),  target :: turb
  integer                  :: n0, n1
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w, p, t, phi
  type(Var_Type),   pointer :: kin, eps, f22, zeta, vis, t2
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),   pointer :: ut, vt, wt
  integer                   :: c, n, sc
  real, contiguous, pointer :: u_mean(:), v_mean(:), w_mean(:),  &
                               p_mean(:), t_mean(:), q_mean(:)
  real, contiguous, pointer :: kin_mean (:), eps_mean(:),  &
                               zeta_mean(:), f22_mean(:)
  real, contiguous, pointer :: uu_res(:), vv_res(:), ww_res(:),  &
                               uv_res(:), vw_res(:), uw_res(:)
  real, contiguous, pointer :: ut_res(:), vt_res(:), wt_res(:), t2_res(:)
  real, contiguous, pointer :: uu_mean(:), vv_mean(:), ww_mean(:)
  real, contiguous, pointer :: uv_mean(:), vw_mean(:), uw_mean(:)
  real, contiguous, pointer :: ut_mean(:), vt_mean(:), wt_mean(:), t2_mean(:)
  real, contiguous, pointer :: phi_mean(:,:)
!==============================================================================!

  if(.not. turb % statistics) return

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  p    => Flow % p
  vis  => turb % vis
  t2   => turb % t2
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_Heat_Fluxes (turb, ut, vt, wt)

  ! Time averaged momentum and energy equations
  u_mean => turb % u_mean;  v_mean => turb % v_mean;  w_mean => turb % w_mean
  p_mean => turb % p_mean;  t_mean => turb % t_mean;  q_mean => turb % q_mean

  ! Time averaged modeled quantities
  kin_mean  => turb % kin_mean;   eps_mean  => turb % eps_mean
  zeta_mean => turb % zeta_mean;  f22_mean  => turb % f22_mean

  ! Time-averaged modelled Reynolds stresses and heat fluxes
  uu_mean => turb % uu_mean;  vv_mean => turb % vv_mean
  ww_mean => turb % ww_mean;  uv_mean => turb % uv_mean
  vw_mean => turb % vw_mean;  uw_mean => turb % uw_mean
  ut_mean => turb % ut_mean;  vt_mean => turb % vt_mean
  wt_mean => turb % wt_mean;  t2_mean => turb % t2_mean

  ! Resolved Reynolds stresses and heat fluxes
  uu_res => turb % uu_res;  vv_res => turb % vv_res;  ww_res => turb % ww_res
  uv_res => turb % uv_res;  vw_res => turb % vw_res;  uw_res => turb % uw_res
  ut_res => turb % ut_res;  vt_res => turb % vt_res;  wt_res => turb % wt_res
  t2_res => turb % t2_res

  n = n1 - n0

  if(n > -1) then

    do c = -Grid % n_bnd_cells, Grid % n_cells

      ! Mean velocities (and temperature)
      u_mean(c) = (u_mean(c) * real(n) + u % n(c)) / real(n+1)
      v_mean(c) = (v_mean(c) * real(n) + v % n(c)) / real(n+1)
      w_mean(c) = (w_mean(c) * real(n) + w % n(c)) / real(n+1)
      p_mean(c) = (p_mean(c) * real(n) + p % n(c)) / real(n+1)

      if(Flow % heat_transfer) then
        t_mean(c) = (t_mean(c) * real(n) + t % n(c)) / real(n+1)
        q_mean(c) = (q_mean(c) * real(n) + t % q(c)) / real(n+1)
      end if

      ! Resolved Reynolds stresses
      uu_res(c) = (uu_res(c) * real(n) + u % n(c) * u % n(c)) / real(n+1)
      vv_res(c) = (vv_res(c) * real(n) + v % n(c) * v % n(c)) / real(n+1)
      ww_res(c) = (ww_res(c) * real(n) + w % n(c) * w % n(c)) / real(n+1)

      uv_res(c) = (uv_res(c) * real(n) + u % n(c) * v % n(c)) / real(n+1)
      uw_res(c) = (uw_res(c) * real(n) + u % n(c) * w % n(c)) / real(n+1)
      vw_res(c) = (vw_res(c) * real(n) + v % n(c) * w % n(c)) / real(n+1)

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        t2_res(c) = (t2_res(c) * real(n) + t % n(c) * t % n(c)) / real(n+1)
        ut_res(c) = (ut_res(c) * real(n) + u % n(c) * t % n(c)) / real(n+1)
        vt_res(c) = (vt_res(c) * real(n) + v % n(c) * t % n(c)) / real(n+1)
        wt_res(c) = (wt_res(c) * real(n) + w % n(c) * t % n(c)) / real(n+1)
      end if

      if(turb % model .eq. K_EPS                 .or.  &
         turb % model .eq. K_EPS_ZETA_F          .or.  &
         turb % model .eq. HYBRID_LES_RANS       .or.  &
         turb % model .eq. RSM_HANJALIC_JAKIRLIC .or.  &
         turb % model .eq. RSM_MANCEAU_HANJALIC ) then

        ! Resolved turbulent heat fluxes
        if(Flow % heat_transfer) then
          t2_mean(c) = (t2_mean(c) * real(n) + t2 % n(c)) / real(n+1)
          ut_mean(c) = (ut_mean(c) * real(n) + ut % n(c)) / real(n+1)
          vt_mean(c) = (vt_mean(c) * real(n) + vt % n(c)) / real(n+1)
          wt_mean(c) = (wt_mean(c) * real(n) + wt % n(c)) / real(n+1)
        end if
      end if

      !-----------------!
      !   K-eps model   !
      !-----------------!
      if(turb % model .eq. K_EPS) then

        ! Time-averaged modeled quantities
        kin_mean(c) = (kin_mean(c) * real(n) + kin % n(c)) / real(n+1)
        eps_mean(c) = (eps_mean(c) * real(n) + eps % n(c)) / real(n+1)
      end if

      !------------------!
      !   K-eps-zeta-f   !
      !------------------!
      if(turb % model .eq. K_EPS_ZETA_F .or.  &
         turb % model .eq. HYBRID_LES_RANS) then

        ! Time-averaged modeled quantities
        kin_mean (c) = (kin_mean (c) * real(n) + kin  % n(c)) / real(n+1)
        eps_mean (c) = (eps_mean (c) * real(n) + eps  % n(c)) / real(n+1)
        zeta_mean(c) = (zeta_mean(c) * real(n) + zeta % n(c)) / real(n+1)
        f22_mean (c) = (f22_mean (c) * real(n) + f22  % n(c)) / real(n+1)
      end if

      !----------------------------!
      !   Reynolds stress models   !
      !----------------------------!
      if(turb % model .eq. RSM_HANJALIC_JAKIRLIC .or.  &
         turb % model .eq. RSM_MANCEAU_HANJALIC) then

        ! Time-averaged modeled quantities (modelled Reynolds stresses)
        uu_mean (c) = (uu_mean (c) * real(n) + uu  % n(c)) / real(n+1)
        vv_mean (c) = (vv_mean (c) * real(n) + vv  % n(c)) / real(n+1)
        ww_mean (c) = (ww_mean (c) * real(n) + ww  % n(c)) / real(n+1)
        uv_mean (c) = (uv_mean (c) * real(n) + uv  % n(c)) / real(n+1)
        uw_mean (c) = (uw_mean (c) * real(n) + uw  % n(c)) / real(n+1)
        vw_mean (c) = (vw_mean (c) * real(n) + vw  % n(c)) / real(n+1)
        kin_mean(c) = (kin_mean(c) * real(n) + kin % n(c)) / real(n+1)
        eps_mean(c) = (eps_mean(c) * real(n) + eps % n(c)) / real(n+1)
        if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
          f22_mean(c) = (f22_mean(c) * real(n) + f22 % n(c)) / real(n+1)
        end if
      end if

      !-------------!
      !   Scalars   !
      !-------------!
      do sc = 1, Flow % n_scalars
        phi      => Flow % scalar(sc)
        phi_mean => turb % scalar_mean
        phi_mean(sc, c) = (phi_mean(sc, c) * real(n) + phi % n(c)) / real(n+1)
      end do
    end do

  end if

  end subroutine
