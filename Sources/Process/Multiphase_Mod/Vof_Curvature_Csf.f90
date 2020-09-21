!==============================================================================!
  subroutine Multiphase_Mod_Vof_Curvature_Csf(grid, mult,                  &
                                              grad_kx, grad_ky, grad_kz,   &
                                              curr_colour)
!------------------------------------------------------------------------------!
!   Computes the Curvature based on Brackbill's CSF using Gauss theorem        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: div_x => r_cell_10,  &
                      div_y => r_cell_11,  &
                      div_z => r_cell_12,  &
                      vof_n => r_node_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Grid_Type)               :: grid
  real                          :: grad_kx    (-grid % n_bnd_cells    &
                                              : grid % n_cells),      &
                                   grad_ky    (-grid % n_bnd_cells    &
                                              : grid % n_cells),      &
                                   grad_kz    (-grid % n_bnd_cells    &
                                              : grid % n_cells),      &
                                   curr_colour(-grid % n_bnd_cells    &
                                              : grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),     pointer :: flow
  type(Var_Type),       pointer :: vof
  integer                       :: s, c, c1, c2, n, i_fac,i_nod, tot_cells,sub
  integer                       :: c_inte, fu, nb, nc
  real, contiguous,     pointer :: fs_x(:), fs_y(:), fs_z(:)
  real                          :: vol_face, grad_face(3), d_n(3)
  real                          :: dotprod, sxyz_mod, sxyz_control, fs, epsloc
  real                          :: dotprod2, stabilize
  real                          :: n_0(3), n_f(3), n_w(3), reflex(3)
  real                          :: theta, theta0, a, b, s_vector(3)
  real                          :: vof_fx, vof_fy, vof_fz, vof_c1, vof_c2, voff
  real                          :: res1, res2, resul, term_c, sumtot
  real                          :: sumx, sumy, sumz, norm_grad, coeff
  real                          :: v1(3), v2(3), v3(3), v4(3)
  real                          :: gf_x, gf_y, gf_z, curv_loc
  real                          :: costheta0, costheta, a1, a2
!==============================================================================!

  vof  => mult % vof
  flow => mult % pnt_flow

  nb = grid % n_bnd_cells
  nc = grid % n_cells

  epsloc = epsilon(epsloc)

  ! Tangent vector to walls/symmetries
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.   &
        Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then

      norm_grad = norm2((/grad_kx(c1),grad_ky(c1),grad_kz(c1)/))

      if (norm_grad > epsloc) then
        v1 = (/grid % sx(s), grid % sy(s), grid % sz(s)/)
        v2 = (/grad_kx(c1), grad_ky(c1), grad_kz(c1)/)
        v3 = Math_Mod_Cross_Product(v1, v2)
        v4 = Math_Mod_Cross_Product(v3, v1)
        ! projection on v4
        norm_grad = norm2(v4)
        if (norm_grad > epsloc) then
          grad_kx(c2) = v4(1) / norm_grad
          grad_ky(c2) = v4(2) / norm_grad
          grad_kz(c2) = v4(3) / norm_grad
        end if
      end if
    else
      grad_kx(c2) = grad_kx(c1)
      grad_ky(c2) = grad_ky(c1)
      grad_kz(c2) = grad_kz(c1)
    end if
  end do

  mult % curv = 0.0

  ! For contact angle
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    curv_loc = 0.0
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
      dotprod = dot_product((/grid % dx(s), grid % dy(s), grid % dz(s)/),   &
                            (/grid % sx(s), grid % sy(s), grid % sz(s)/))

      dotprod2 = grad_kx(c1) * grid % dx(s)   &
               + grad_ky(c1) * grid % dy(s)   &
               + grad_kz(c1) * grid % dz(s)

      gf_x = grad_kx(c1) + grid % sx(s) / dotprod    &
           * (curr_colour(c2) - curr_colour(c1) - dotprod2)
      gf_y = grad_ky(c1) + grid % sy(s) / dotprod    &
           * (curr_colour(c2) - curr_colour(c1) - dotprod2)
      gf_z = grad_kz(c1) + grid % sz(s) / dotprod    &
           * (curr_colour(c2) - curr_colour(c1) - dotprod2)

      norm_grad = norm2((/gf_x, gf_y, gf_z/))

      if(norm_grad > epsloc) then

        gf_x = gf_x / (norm_grad + epsloc)
        gf_y = gf_y / (norm_grad + epsloc)
        gf_z = gf_z / (norm_grad + epsloc)

        costheta0 = dot_product((/gf_x, gf_y, gf_z/),   &
                                 (/grid % sx(s), grid % sy(s), grid % sz(s)/)) &
                   / grid % s(s)
        theta0 = acos(costheta0)

        theta = vof % q(c2) * PI /180.0
        costheta = cos(theta)

        a1 = cos(theta0 - theta)
        a2 = 1.0 - costheta0 * costheta0

        a = (costheta - costheta0 * a1) / (a2 + epsloc)
        b = (a1 - costheta0 * costheta) / (a2 + epsloc)

        grad_kx(c1) = b * gf_x + a * grid % sx(s) / grid % s(s)
        grad_ky(c1) = b * gf_y + a * grid % sy(s) / grid % s(s)
        grad_kz(c1) = b * gf_z + a * grid % sz(s) / grid % s(s)
        grad_kx(c2) = grad_kx(c1)
        grad_ky(c2) = grad_ky(c1)
        grad_kz(c2) = grad_kz(c1)
      end if

    end if

  end do

  call Grid_Mod_Exchange_Cells_Real(grid, grad_kx(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, grad_ky(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, grad_kz(-nb:nc))

  ! Normalize vector at cells
  do c = 1, grid % n_cells
    norm_grad = sqrt(grad_kx(c) ** 2 + grad_ky(c) ** 2 + grad_kz(c) ** 2)
    grad_kx(c) = grad_kx(c) / (norm_grad + epsloc)
    grad_ky(c) = grad_ky(c) / (norm_grad + epsloc)
    grad_kz(c) = grad_kz(c) / (norm_grad + epsloc)
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, grad_kx(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, grad_ky(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, grad_kz(-nb:nc))

  !--------------------!
  !   Find Curvature   !
  !--------------------!

  mult % curv = 0.0

  ! Find divergence of normals
  call Field_Mod_Grad_Component(flow, grad_kx(-nb:nc),  &
                                1,    div_x  (-nb:nc),  &
                                impose_symmetry = .false.)
  call Field_Mod_Grad_Component(flow, grad_ky(-nb:nc),  &
                                2,    div_y  (-nb:nc),  &
                                impose_symmetry = .false.)
  call Field_Mod_Grad_Component(flow, grad_kz(-nb:nc),  &
                                3,    div_z  (-nb:nc),  &
                                impose_symmetry = .false.)

  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_x(-nb:nc)
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_y(-nb:nc)
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_z(-nb:nc)

  call Grid_Mod_Exchange_Cells_Real(grid, mult % curv)

  ! At boundaries
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    mult % curv(c2) = mult % curv(c1)
  end do

  call Multiphase_Mod_Vof_Smooth_Curvature(grid, mult,                  &
                          grad_kx(-nb:nc), grad_ky(-nb:nc), grad_kz(-nb:nc))

  ! At boundaries
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
      if(abs(mult % curv(c1)) > epsloc) then
        mult % curv(c1) = mult % curv(c2)
      endif
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, mult % curv)

  end subroutine
