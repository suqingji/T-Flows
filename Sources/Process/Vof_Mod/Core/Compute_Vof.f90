!==============================================================================!
  subroutine Compute_Vof(Vof, Sol, dt, curr_dt)
!------------------------------------------------------------------------------!
!   Solves Volume Fraction equation using UPWIND ADVECTION and CICSAM          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),   target :: Vof
  type(Solver_Type), target :: Sol
  real                      :: dt
  integer, intent(in)       :: curr_dt  ! current time step
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: fun
  type(Face_Type),   pointer :: v_flux
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  real, contiguous,  pointer :: beta_f(:)
  real, contiguous,  pointer :: beta_c(:)
  real, contiguous,  pointer :: c_d(:)
  real                       :: courant_max
  integer                    :: i_sub, n_sub, wrong_vf, n_wrong_vf0, n_wrong_vf1
  integer                    :: s, c, c1, c2, fu, corr
!==============================================================================!

  call Cpu_Timer % Start('Compute_Vof (without solvers)')

  call User_Mod_Beginning_Of_Compute_Vof(Vof, Sol, curr_dt)

  ! Take aliases
  Flow   => Vof % pnt_flow
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  fun    => Vof % fun
  beta_f => Vof % beta_f
  beta_c => Vof % beta_c
  c_d    => Vof % c_d
  A      => Sol % A
  b      => Sol % b % val

  if(fun % adv_scheme .eq. CICSAM .or. &
     fun % adv_scheme .eq. STACS) then

    if(fun % adv_scheme .eq. CICSAM) then
      ! Compute courant Number close to the interface:
      call Vof % Max_Courant_Number(dt, 1, courant_max)

      n_sub = min(max(ceiling(courant_max / Vof % courant_max_param), 1),  &
                  Vof % n_sub_param)

      ! Warning if Courant Number is exceeded
      if(n_sub > 1) then
        if(this_proc < 2) then
          call File % Append_For_Writing_Ascii('alert-dt-vof.dat', fu)
          write(fu,*) 'Courant Number was exceded at iteration: ', curr_dt
          write(fu,*) 'Co_max = ', courant_max
          write(fu,*) 'Try reducing time step'
          close(fu)
        end if
      end if
    else
      n_sub = 1
    endif
  else
    ! Old volume fraction:
    fun % o(:) = fun % n(:)
  end if

  if(fun % adv_scheme .eq. UPWIND) then

    !------------------------------------------------!
    !   Discretize the system for the vof function   !
    !------------------------------------------------!
    call Vof % Discretize(A, b, dt)

    ! Solve System
    call Vof % Solve_System(Sol, b)

    call Grid % Exchange_Cells_Real(fun % n)

    !-----------------------------!
    !   Correct Volume Fraction   !
    !-----------------------------!
    do c = 1, Grid % n_cells
      fun % n(c) = max(min(fun % n(c),1.0),0.0)
    end do

  else if(fun % adv_scheme .eq. CICSAM .or. &
           fun % adv_scheme .eq. STACS) then

    do i_sub = 1, n_sub

      ! Courant number full domain:
      call Vof % Max_Courant_Number(dt / real(n_sub), 0, courant_max)

      !---------------------------!
      !   Predict beta at faces   !
      !---------------------------!

      ! Impose zero gradient at boundaries
      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        if(c2 < 0) then
          if(Grid % Bnd_Cond_Type(c2) .ne. INFLOW) then
            fun % n(c2) = fun % n(c1)
          end if
        end if
      end do

      ! Old volume fraction:
      fun % o(:) = fun % n(:)

      ! Compute gradient:
      call Flow % Grad_Variable(fun)

      call Vof % Predict_Beta()

      do corr = 1, Vof % corr_num_max

        !------------------------------------------------!
        !   Discretize the system for the vof function   !
        !------------------------------------------------!
        call Vof % Discretize(A, b, dt / real(n_sub))

        ! Solve System
        call Vof % Solve_System(Sol, b)

        do s = 1, Grid % n_faces
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)
          if(c2 < 0) then
            if(Grid % Bnd_Cond_Type(c2) .ne. INFLOW) then
              fun % n(c2) = fun % n(c1)
            end if
          end if
        end do

        n_wrong_vf0 = 0
        n_wrong_vf1 = 0
        wrong_vf = 0

        ! Determine if 0 <= fun <= 1.0
        do c = 1, Grid % n_cells
          if(fun % n(c) < -fun % tol) then
            n_wrong_vf0 = n_wrong_vf0 + 1
          end if

          if(fun % n(c) - 1.0 > fun % tol) then
            n_wrong_vf1 = n_wrong_vf1 + 1
          end if
        end do

        call Grid % Exchange_Cells_Real(fun % n)

        !---------------------------!
        !   Correct beta at faces   !
        !---------------------------!

        if(n_wrong_vf0 > 0 .or. n_wrong_vf1 > 0) then
          wrong_vf = 1
        end if

        call Comm_Mod_Global_Sum_Int(wrong_vf)

        if(wrong_vf == 0) then
          goto 1
        else
          call Vof % Correct_Beta()
        end if

      end do
1     continue

      !------------------------!
      !   Correct boundaries   !
      !------------------------!
      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        if(c2 < 0) then
          if(Grid % Bnd_Cond_Type( c2) .ne. INFLOW) then
            if(fun % n(c2) < FEMTO) then
               fun % n(c2) = 0.0
            end if
            if(fun % n(c2) - 1.0 >= FEMTO) then
               fun % n(c2) = 1.0
            end if
          end if
        end if

      end do

      !--------------------------------------!
      !   Correct Interior Volume Fraction   !
      !--------------------------------------!
      do c = 1, Grid % n_cells
        if(fun % n(c) < FEMTO) then
          fun % n(c) = 0.0
        end if
        if(fun % n(c) - 1.0 >= FEMTO) then
          fun % n(c) = 1.0
        end if
      end do

      call Grid % Exchange_Cells_Real(fun % n)
    end do


  end if

  !------------------------------!
  !   Volume fraction at faces   !
  !------------------------------!

  if(fun % adv_scheme .eq. UPWIND) then

    ! At boundaries
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid % Bnd_Cond_Type(c2) .ne. INFLOW) then
          fun % n(c2) = fun % n(c1)
        end if
      end if
    end do

  else if(fun % adv_scheme .eq. CICSAM .or. fun % adv_scheme .eq. STACS) then

    ! At boundaries
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW) then
          fun % n(c2) = fun % n(c1)
        else if(Grid % Bnd_Cond_Type(c2) .eq. PRESSURE) then
          fun % n(c2) = fun % n(c1)
        else if(Grid % Bnd_Cond_Type(c2) .eq. CONVECT) then
          fun % n(c2) = fun % n(c1)
        else if(Grid % Bnd_Cond_Type(c2) .eq. INFLOW) then
        else
          fun % n(c2) = fun % n(c1)
        end if
      end if
    end do

  end if

  call Flow % Grad_Variable(fun)

  call User_Mod_End_Of_Compute_Vof(Vof, Sol, curr_dt)

  call Cpu_Timer % Stop('Compute_Vof (without solvers)')

  end subroutine
