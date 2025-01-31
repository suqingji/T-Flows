!==============================================================================!
  subroutine Turb_Mod_Const_K_Eps_Zeta_F(turb)
!------------------------------------------------------------------------------!
!   Initializes constants for k-eps-zeta-f and hybrid k-eps-zeta-f             !
!   turbulence models.                                                         ! 
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
!==============================================================================!

  c_1e    =  1.4
  c_2e    =  1.9
  c_mu    =  0.09
  c_mu_d  =  0.22
  c_mu25  = sqrt(sqrt(c_mu))
  c_mu75  = c_mu25**3
  kappa   =  0.41
  e_log   =  8.342
  c_l     =  0.36
  c_t     =  6.0
  c_nu    = 85.0
  alpha   =  0.012
  c_f1    =  1.4
  c_f2    =  0.3

  ! Roughness coefficient
  turb % z_o = 0.0 

  ! Constants for buoyancy wall function
  c_mu_theta   =  0.1225
  c_mu_theta5  =  0.35
  kappa_theta  =  0.38         ! von Karman constant for temperature

  ! Constants for GGDH turbulent flux
  c_theta =  0.22

  ! Constants for AFM turbulent flux
  afm_psi      = 0.1
  afm_eta      = 0.1

  turb % kin  % sigma = 1.0
  turb % eps  % sigma = 1.3
  turb % zeta % sigma = 1.2
  turb % t2   % sigma = 1.0

  end subroutine
