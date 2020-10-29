!==============================================================================!
  subroutine Multiphase_Mod_Vof_Curvature_Nodal(grid, mult, smooth_k,   &
                                                var_node_k)
!------------------------------------------------------------------------------!
!   Computes the Curvature based on Brackbill's CSF using Least Squares method !
!   and nodal values                                                           !
!                                                                              !
!   Arguments:                                                                 !
!   - smooth_k, var_node_k          : Both hold the distance function or vof at!
!                                     cell centers and nodes respectively.     !
!                                     In any case, they have been smoothed out !
!                                     previously to enhance curvature          !
!                                     calculation                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: div_x  => r_cell_11,  &
                      div_y  => r_cell_12,  &
                      div_z  => r_cell_13,  &
                      div_xx => r_cell_14,  &
                      div_yy => r_cell_15,  &
                      div_zz => r_cell_16,  &
                      grad_x => r_node_05,  &
                      grad_y => r_node_06,  &
                      grad_z => r_node_07,  &
                      corr_x => r_node_01,  &
                      corr_y => r_node_02,  &
                      corr_z => r_node_03,  &
                      summ_n => r_node_12,  &
                      mark   => i_node_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Grid_Type)               :: grid
  real                          :: smooth_k(-grid % n_bnd_cells :   &
                                             grid % n_cells),       &
                                   var_node_k(1:grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),     pointer :: flow
  type(Var_Type),       pointer :: vof
  integer                       :: s, c, c1, c2, n, i_fac, i_nod, tot_cells,sub
  integer                       :: c_inte, fu, nb, nc, nn
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
  real                          :: c_c
  real                          :: gf_x, gf_y, gf_z, a1, a2, costheta, costheta0
!==============================================================================!

  vof  => mult % vof
  flow => mult % pnt_flow

  nb = grid % n_bnd_cells
  nc = grid % n_cells
  nn = grid % n_nodes

  epsloc = epsilon(epsloc)

  ! Calculate gradients at nodes
  call Field_Mod_Grad_Component_Cells_To_Nodes(flow,                  &
                                               smooth_k, var_node_k,  &
                                            1, grad_x(1:nn))
  call Field_Mod_Grad_Component_Cells_To_Nodes(flow,                  &
                                               smooth_k, var_node_k,  &
                                            2, grad_y(1:nn))
  call Field_Mod_Grad_Component_Cells_To_Nodes(flow,                  &
                                               smooth_k, var_node_k,  &
                                            3, grad_z(1:nn))

  ! Tangent vector to symmetries
  mark(1:nn) = 0

  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(     (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY)    &
       .or. (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) ) then

      do i_nod = 1, grid % faces_n_nodes(s)
        n = grid % faces_n(i_nod,s)
        if(mark(n) == 0) then
          v1 = (/grid % sx(s), grid % sy(s), grid % sz(s)/)
          v2 = (/grad_x(n), grad_y(n), grad_z(n)/)
          v3 = Math_Mod_Cross_Product(v1, v2)
          v4 = Math_Mod_Cross_Product(v3, v1)
          ! projection on v4
          norm_grad = norm2(v4)
          if (norm_grad > epsloc) then
            grad_x(n) = v4(1) / norm_grad
            grad_y(n) = v4(2) / norm_grad
            grad_z(n) = v4(3) / norm_grad
            mark(n) = 1
          end if
        end if
      end do

    end if
  end do

  ! Normalize vectors at nodes
  do n = 1, grid % n_nodes
    norm_grad = sqrt(grad_x(n) ** 2 + grad_y(n) ** 2 + grad_z(n) ** 2)
    if (norm_grad >= epsloc) then
      grad_x(n) = grad_x(n) / (norm_grad + epsloc)
      grad_y(n) = grad_y(n) / (norm_grad + epsloc)
      grad_z(n) = grad_z(n) / (norm_grad + epsloc)
    else
      grad_x(n) = 0.0
      grad_y(n) = 0.0
      grad_z(n) = 0.0
    end if
  end do


  call Grid_Mod_Exchange_Nodes_Real(grid, grad_x(1:nn))
  call Grid_Mod_Exchange_Nodes_Real(grid, grad_y(1:nn))
  call Grid_Mod_Exchange_Nodes_Real(grid, grad_z(1:nn))

  ! Interpolate node values to cells
  call Field_Mod_Interpolate_Nodes_To_Cells(flow, grad_x(1:nn), div_x(-nb:nc))
  call Field_Mod_Interpolate_Nodes_To_Cells(flow, grad_y(1:nn), div_y(-nb:nc))
  call Field_Mod_Interpolate_Nodes_To_Cells(flow, grad_z(1:nn), div_z(-nb:nc))

  ! Correct for contact angle at walls
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
      dotprod = dot_product((/grid % dx(s), grid % dy(s), grid % dz(s)/),   &
                            (/grid % sx(s), grid % sy(s), grid % sz(s)/))

      dotprod2 = div_x(c1) * grid % dx(s)   &
               + div_y(c1) * grid % dy(s)   &
               + div_z(c1) * grid % dz(s)

      gf_x = div_x(c1) + grid % sx(s) / dotprod    &
           * (smooth_k(c2) - smooth_k(c1) - dotprod2)
      gf_y = div_y(c1) + grid % sy(s) / dotprod    &
           * (smooth_k(c2) - smooth_k(c1) - dotprod2)
      gf_z = div_z(c1) + grid % sz(s) / dotprod    &
           * (smooth_k(c2) - smooth_k(c1) - dotprod2)

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

        div_x(c1) = b * gf_x + a * grid % sx(s) / grid % s(s)
        div_y(c1) = b * gf_y + a * grid % sy(s) / grid % s(s)
        div_z(c1) = b * gf_z + a * grid % sz(s) / grid % s(s)
        div_x(c2) = div_x(c1)
        div_y(c2) = div_y(c1)
        div_z(c2) = div_z(c1)
      end if

    end if

  end do

  call Grid_Mod_Exchange_Cells_Real(grid, div_x(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, div_y(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, div_z(-nb:nc))
  ! Normalize vector at cells
  do c = 1, grid % n_cells
    norm_grad = sqrt(div_x(c) ** 2 + div_y(c) ** 2 + div_z(c) ** 2)
    if (norm_grad >= epsloc) then
      div_x(c) = div_x(c) / (norm_grad + epsloc)
      div_y(c) = div_y(c) / (norm_grad + epsloc)
      div_z(c) = div_z(c) / (norm_grad + epsloc)
    else
      div_x(c) = 0.0
      div_y(c) = 0.0
      div_z(c) = 0.0
    end if
  end do

  ! Correct nodal values at walls
  mark(1:nn)   = 0
  corr_x(1:nn) = 0.0
  corr_y(1:nn) = 0.0
  corr_z(1:nn) = 0.0
  summ_n(1:nn) = 0.0
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
      do i_nod = 1, grid % faces_n_nodes(s)
        n = grid % faces_n(i_nod,s)
        mark(n) = 1
        summ_n(n) = summ_n(n) + 1.0
        corr_x(n) = corr_x(n) + div_x(c2)
        corr_y(n) = corr_y(n) + div_y(c2)
        corr_z(n) = corr_z(n) + div_z(c2)
      end do
    end if
  end do

  do n = 1, nn
    if (mark(n) == 1) then
      grad_x(n) = corr_x(n) / (summ_n(n) + epsloc)
      grad_y(n) = corr_y(n) / (summ_n(n) + epsloc)
      grad_z(n) = corr_z(n) / (summ_n(n) + epsloc)
      ! Normalize
      norm_grad = norm2((/grad_x(n), grad_y(n), grad_z(n)/))
      grad_x(n) = grad_x(n) / (norm_grad + epsloc)
      grad_y(n) = grad_y(n) / (norm_grad + epsloc)
      grad_z(n) = grad_z(n) / (norm_grad + epsloc)
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, div_x(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, div_y(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, div_z(-nb:nc))

  call Grid_Mod_Exchange_Nodes_Real(grid, grad_x)
  call Grid_Mod_Exchange_Nodes_Real(grid, grad_y)
  call Grid_Mod_Exchange_Nodes_Real(grid, grad_z)

  !--------------------!
  !   Find Curvature   !
  !--------------------!

  mult % curv = 0.0

  !! Derivatives of normals using nodes

  !call Field_Mod_Grad_Component_Nodes_To_Cells(flow,                        &
  !                                            div_x(-nb:nc), grad_x(1:nn),  &
  !                                         1, div_xx(-nb:nc))
  !mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_xx(-nb:nc)

  !call Field_Mod_Grad_Component_Nodes_To_Cells(flow,                        &
  !                                            div_y(-nb:nc), grad_y(1:nn),  &
  !                                         2, div_yy(-nb:nc))
  !mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_yy(-nb:nc)

  !call Field_Mod_Grad_Component_Nodes_To_Cells(flow,                        &
  !                                            div_z(-nb:nc), grad_z(1:nn),  &
  !                                         3, div_zz(-nb:nc))
  !mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_zz(-nb:nc)

  call Multiphase_Mod_Vof_Gradient_Iterative(mult, div_x(-nb:nc),    &
                                                   div_xx(-nb:nc),   &
                                                   div_yy(-nb:nc),   &
                                                   div_zz(-nb:nc))
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_xx(-nb:nc)

  call Multiphase_Mod_Vof_Gradient_Iterative(mult, div_y(-nb:nc),    &
                                                   div_xx(-nb:nc),   &
                                                   div_yy(-nb:nc),   &
                                                   div_zz(-nb:nc))
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_yy(-nb:nc)

  call Multiphase_Mod_Vof_Gradient_Iterative(mult, div_z(-nb:nc),    &
                                                   div_xx(-nb:nc),   &
                                                   div_yy(-nb:nc),   &
                                                   div_zz(-nb:nc))
  mult % curv(-nb:nc) = mult % curv(-nb:nc) - div_zz(-nb:nc)

  call Grid_Mod_Exchange_Cells_Real(grid, mult % curv(-nb:nc))

  ! At boundaries
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    mult % curv(c2) = mult % curv(c1)
  end do

  call Multiphase_Mod_Vof_Smooth_Curvature(grid, mult,              &
                      div_x(-nb:nc), div_y(-nb:nc), div_z(-nb:nc))

  call Grid_Mod_Exchange_Cells_Real(grid, mult % curv(-nb:nc))

  end subroutine
