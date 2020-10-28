!==============================================================================!
  subroutine Multiphase_Mod_Vof_Smooth_Curvature(grid, mult,                  &
                                                 norm_nx, norm_ny, norm_nz)
!------------------------------------------------------------------------------!
!   Smoothes curvature in two steps: first a smoothing curvature around the    !
!   Interface and second in the direction of the normal. This technique can    !
!   be found at https://spiral.imperial.ac.uk/handle/10044/1/28101
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: k_star       => r_cell_14,    &
                      sum_k_weight => r_cell_18,    &
                      sum_weight   => r_cell_19,    &
                      wint         => i_cell_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Grid_Type)               :: grid
  real                          :: norm_nx(-grid % n_bnd_cells:grid % n_cells)
  real                          :: norm_ny(-grid % n_bnd_cells:grid % n_cells)
  real                          :: norm_nz(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  real          , pointer :: vof(:)
  integer                 :: i, s, c, c1, c2, c_iter, iter, i_fac, nb, nc
  integer                 :: face_init, face_end, face_step
  real                    :: fs, w_v1, w_v2, w_m1, w_m2
  real                    :: weight_s, weight_n, sum1, sum2
  real                    :: norma, epsloc, curvf, dotprod
!==============================================================================!

  ! Take aliases
  if (mult % d_func) then
    vof => mult % dist_func % oo
  else
    vof => mult % vof % n
  end if

  nb = grid % n_bnd_cells
  nc = grid % n_cells

  epsloc = epsilon(epsloc)

  sum_k_weight(-nb:nc) = 0.0
  sum_weight  (-nb:nc) = 0.0
  weight_s             = 8.0
  weight_n             = 8.0

  !-------------------------!
  !   Smoothing curvature   !
  !-------------------------!

  ! What is the curvature at boundaries??? For now zero gradient
  ! Preliminary results using wetting (with symmetry in some boundaries) show
  ! it is better not to take into aacount the boundaries

  ! Interior faces
  do s = grid % n_bnd_faces + 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    w_v1 = (1.0 - 2.0 * abs(0.5 - vof(c1))) ** weight_s
    w_v2 = (1.0 - 2.0 * abs(0.5 - vof(c2))) ** weight_s
    sum_k_weight(c1) = sum_k_weight(c1) + mult % curv(c2) * w_v2
    sum_weight(c1) = sum_weight(c1) + w_v2

    sum_k_weight(c2) = sum_k_weight(c2) + mult % curv(c1) * w_v1
    sum_weight(c2) = sum_weight(c2) + w_v1
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, sum_k_weight(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, sum_weight  (-nb:nc))

  do c = 1, grid % n_cells
    w_v1 = (1.0 - 2.0 * abs(0.5 - vof(c))) ** weight_s
    k_star(c) = (w_v1 * mult % curv(c) + sum_k_weight(c))    &
              / (w_v1 + sum_weight(c) + epsloc)
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, k_star(-nb:nc))

  !-------------------------------------------------------------------------!
  !   Smoothing curvature in the direction of the normal to the interface   !
  !-------------------------------------------------------------------------!

  sum_k_weight(-nb:nc) = 0.0
  sum_weight  (-nb:nc) = 0.0

  ! Interior faces
  do s = grid % n_bnd_faces + 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    w_v1 = (1.0 - 2.0 * abs(0.5 - vof(c1))) ** weight_n
    w_v2 = (1.0 - 2.0 * abs(0.5 - vof(c2))) ** weight_n

    w_m1 = abs(dot_product((/norm_nx(c1), norm_ny(c1), norm_nz(c1)/),     &
                           (/grid % dx(s), grid % dy(s), grid % dz(s)/)   &
                           / grid % d(s))) ** weight_n

    w_m2 = abs(dot_product((/norm_nx(c2), norm_ny(c2), norm_nz(c2)/),     &
                           (/grid % dx(s), grid % dy(s), grid % dz(s)/)   &
                           / (-grid % d(s)))) ** weight_n

    sum_k_weight(c1) = sum_k_weight(c1) + k_star(c2) * w_v2 * w_m2
    sum_weight(c1) = sum_weight(c1) + w_v2 * w_m2

    sum_k_weight(c2) = sum_k_weight(c2) + k_star(c1) * w_v1 * w_m1
    sum_weight(c2) = sum_weight(c2) + w_v1 * w_m1
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, sum_k_weight(-nb:nc))
  call Grid_Mod_Exchange_Cells_Real(grid, sum_weight  (-nb:nc))

  do c = 1, grid % n_cells
    w_v1 = (1.0 - 2.0 * abs(0.5 - vof(c))) ** weight_n
    mult % curv(c) = (w_v1 * k_star(c) + sum_k_weight(c))    &
                   / (w_v1 + sum_weight(c) + epsloc)
  end do

  !--------------------------------------------------------!
  !   Expand curvature to two additional layers of cells   !
  !--------------------------------------------------------!
  call Grid_Mod_Exchange_Cells_Real(grid, mult % curv)

  wint = 0
  do c = 1, grid % n_cells
    if(vof(c) > 0.0 .and. vof(c) < 1.0) then
      wint(c) = 1
    end if
  end do

  call Grid_Mod_Exchange_Cells_Int(grid, wint)

  iter = 2
  do i = 1, iter
    do c = 1, grid % n_cells
      if (wint(c) .ne. i) then
        sum1 = 0.0
        sum2 = 0.0
        do i_fac = 1, grid % cells_n_faces(c)
          s  =  grid % cells_f(i_fac, c)
          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)
          if (all((/c1,c2/) > 0)) then
            if (c .eq. c1) then
              if (wint(c2) == i) then
                w_m2 = abs(dot_product((/norm_nx(c2), norm_ny(c2), norm_nz(c2)/),     &
                                       (/grid % dx(s), grid % dy(s), grid % dz(s)/)   &
                                       / (-grid % d(s)))) ** weight_n

                sum1 = sum1 + mult % curv(c2) * w_m2
                sum2 = sum2 + w_m2
              end if
            else
              if (wint(c1) == i) then
                w_m1 = abs(dot_product((/norm_nx(c1), norm_ny(c1), norm_nz(c1)/),     &
                                       (/grid % dx(s), grid % dy(s), grid % dz(s)/)   &
                                       / grid % d(s))) ** weight_n

                sum1 = sum1 + mult % curv(c1) * w_m1
                sum2 = sum2 + w_m1
              end if
            end if
          end if
        end do
        if ((sum2 > epsloc) .and. wint(c) < i) then
          mult % curv(c) = sum1 / sum2
          wint(c) = i + 1
        else if(wint(c) == 0) then
          mult % curv(c) = 0.0
        end if
      end if
    end do

    call Grid_Mod_Exchange_Cells_Int(grid, wint)
    call Grid_Mod_Exchange_Cells_Real(grid, mult % curv)
  end do

  end subroutine
