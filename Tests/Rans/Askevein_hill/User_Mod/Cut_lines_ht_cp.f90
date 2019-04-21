!==============================================================================!
   subroutine User_Mod_Cut_lines_ht_cp(flow, save_name) 
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for scalar.        !
!   It is called from "Compute_Scalar" function, just before calling the       !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent should the user want to stabilize the         !
!   system for always positive variables, for example.                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod, only: Field_Type, heat_transfer, heat_flux,  &
                       density, viscosity, capacity, conductivity
  use Grid_Mod,   only: Grid_Type
  use Var_Mod,    only: Var_Type
  use Matrix_Mod, only: Matrix_Type
  use Bulk_Mod,   only: Bulk_Type  
  use Rans_Mod
  use Comm_Mod                       ! parallel stuff
  use Name_Mod,  only: problem_name
  use Const_Mod                      ! constants
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  character(len=*)         :: save_name
!   type(Var_Type),   target :: phi
!   type(Matrix_Type)        :: a_matrix
!   real, dimension(:)       :: b_vector
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  real,            pointer :: flux(:)
 !------------------------------------------------------------------------------! 
    INTERFACE
    LOGICAL FUNCTION Approx_Real(A,B,tol)
      REAL           :: A,B
      REAL, OPTIONAL :: tol
    END FUNCTION Approx_Real
  END INTERFACE

!-----------------------------[Parameters]-----------------------------!

  REAL :: Rad_2, Ufric, Rad_1 
!------------------------------[Calling]-------------------------------!
  INTEGER             :: Nprob, pl, c, i, count, k, c1, c2, s, N_hor, n
  INTEGER             :: Nfound
  CHARACTER           :: name_ht*12, name_cp*12, name_inflow*12
  character(len=80)   :: res_name_cp, res_name_ht
  character(len=80)   :: store_name
  REAL,ALLOCATABLE    :: y_p(:), x_p(:), z_p(:), Ump(:), Vmp(:), Wmp(:), & 
                                 uup(:), vvp(:), wwp(:), &
                                 uvp(:), uwp(:), vwp(:), &
                                 Tmp(:), TTp(:),         &
                                 uTp(:), vTp(:), wTp(:), &
                                 Ksgsp(:),               & 
                                 var_1(:), var_2(:), var_3(:), Rad_mp(:), &
                                 var_4(:), var_5(:), R_p(:)  
  INTEGER,ALLOCATABLE :: Np(:), Ncount(:)
  REAL                :: R, Urad_n, Utan_n, R1, R2, Urad, Utan, dummy 
  REAL                :: xc_max, yc_max, zc_max, range
  logical             :: there
!==============================================================================!
  ! Take aliases
  grid => flow % pnt_grid
  flux => flow % flux
  bulk => flow % bulk
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t

!--------------------------------------------------------------------! 

  ! Store the name
  store_name = problem_name
  problem_name = save_name

  call Name_File(0, res_name_ht, "-ht.dat")
  call Name_File(0, res_name_cp, "-cp.dat")
 
  N_hor = 2
  Nfound = 0
  zc_max = 0.0
  range = HUGE

  allocate(x_p(N_hor))
  allocate(y_p(N_hor))
  allocate(R_p(N_hor))

  do c = -1, -grid % n_bnd_cells, -1
    if( Grid_Mod_Bnd_Cond_Type(grid, c) .eq. WALL .or.  &
        Grid_Mod_Bnd_Cond_Type(grid, c) .eq. WALLFL) then 
      if(grid % zc(c) > zc_max) then
        zc_max = grid % zc(c)
      end if
    end if
  end do
  
  call Comm_Mod_Global_Max_Real(zc_max)
  
  if(this_proc< 2) write(*,*) 'Hill top height is ',  zc_max

  call Comm_Mod_Wait

  xc_max = 0.0
  yc_max = 0.0

  do c = 1, grid % n_cells
    range = min(range, grid % delta(c))
  end do

  call Comm_Mod_Global_Min_Real(range)

  do c = -1, -grid % n_bnd_cells, -1
    if( Grid_Mod_Bnd_Cond_Type(grid, c) .eq. WALL .or. &
        Grid_Mod_Bnd_Cond_Type(grid, c) .eq. WALLFL) then 
      if(Approx_Real( grid % zc(c), zc_max)) then
        xc_max = grid % xc(c)
        yc_max = grid % yc(c)
      end if
    end if
  end do
  
  call Comm_Mod_Global_Sum_Real(xc_max)
  call Comm_Mod_Global_Sum_Real(yc_max)
  
  if(this_proc< 2) write(*,*)'Coordinate of hill top (HT)', xc_max,  yc_max , zc_max
  
  call Comm_Mod_Wait
 
  ! HT 
  x_p(1) = xc_max
  y_p(1) = yc_max
  
 !CP                                   !
 !distance between HT and CP is 400 m  !
 
!  x_p(2) = xc_max - 410.0*cos(78.0*3.14159/180.0)
!  y_p(2) = yc_max - 410.0*sin(73.0*3.14159/180.0)

  x_p(2) = xc_max - .5*cos(78.0*3.14159/180.0)
  y_p(2) = yc_max - .5*sin(73.0*3.14159/180.0)
!-Radius for averaging --adjust radius according to number of counts 
  R_p(1) = range
  R_p(2) = range

  !------------------!
  !   Read 1d file   !
  !------------------!
  inquire(file='Z_coordinate.dat', exist=there)
  if(.not. there) then
    if(this_proc < 2) then
      print *, '==============================================================='
      print *, 'In order to extract profiles and write them in ascii files'
      print *, 'the code has to read cell-faces coordinates '
      print *, 'in wall-normal direction in the ascii file ''case_name.1d.'''
      print *, 'The file format should be as follows:'
      print *, '10  ! number of cells + 1'
      print *, '1 0.0'
      print *, '2 0.1'
      print *, '3 0.2'
      print *, '... '
      print *, '==============================================================='
    end if

    ! Restore the name and return
    problem_name = store_name
    return
  end if

  if(this_proc< 2) write(6, *) '# Now reading Z_coordinate.dat '  

  open(9, FILE='Z_coordinate.dat')
!---- write the number of probes 
  read(9,*) Nprob
  allocate(z_p(Nprob))

!---- write the probe coordinates out
  do pl=1,Nprob
    read(9,*) dummy, z_p(pl)
  end do
  close(9)
    
  allocate(Np(Nprob));     Np     = 0 
  allocate(Ump(Nprob));    Ump    = 0.0
  allocate(uwp(Nprob));    uwp    = 0.0
  allocate(Rad_mp(Nprob)); Rad_mp = 0.0
  allocate(var_1(Nprob));  var_1  = 0.0
  allocate(var_2(Nprob));  var_2  = 0.0
  allocate(var_3(Nprob));  var_3  = 0.0
  allocate(var_4(Nprob));  var_4  = 0.0

  allocate(Ncount(Nprob)); Ncount=0
  count = 0

!+++++++++++++++++++++++++++++!
!     average the results     !
!+++++++++++++++++++++++++++++!

  do k = 1, 2 
    do i = 1, Nprob-1
      101 continue
      Ncount(i) = 0 
      do c = 1, grid % n_cells - grid % comm % n_buff_cells
        Rad_1 = sqrt((grid % xc(c)-x_p(k))**2 + (grid % yc(c)-y_p(k))**2)
        Rad_2 = grid % wall_dist(c) !grid % zc(c) 
        if(Rad_1 < R_p(k)) then
          if(Rad_2 > z_p(i) .and. Rad_2 < z_p(i+1)) then
            Rad_mp(i) = Rad_mp(i) + Rad_2 
            Ncount(i) = Ncount(i) + 1
          end if
        end if
      end do  

    call Comm_Mod_Wait
    call Comm_Mod_Global_Sum_Int(Ncount(i))

    if(Ncount(i) /= 0) then
      Nfound    = Ncount(i) !max(Nfound, Ncount(i)) 
    end if

    call Comm_Mod_Wait

    if(Nfound > 6) then
      R_p(k) = R_p(k) * 0.9 
      write(*,*) 'vise od 6', Nfound, R_p(k), k
      goto 101  
    else  
      continue
    end if
!    write(*,*) 'It found ', Nfound, ' points!'

    end do
  end do

  do i = 1, Nprob
    Ncount(i) = 0
    Rad_mp(i) = 0.0 
  end do

!  write(*,*) x_p(2), x_p(1)
!  write(*,*) y_p(2), y_p(1)
   
  do k = 1, 2 
    do i = 1, Nprob-1
      do c = 1, grid % n_cells
        Rad_1 = sqrt((grid % xc(c)-x_p(k))**2 + (grid % yc(c)-y_p(k))**2)
        Rad_2 = grid % wall_dist(c) 
        if(Rad_1 < R_p(k)) then
          if(Rad_2 > z_p(i) .and. Rad_2 < z_p(i+1)) then
            Ncount(i) = Ncount(i) + 1
            Ump(i)   = Ump(i) + u % n(c)
            uwp(i)   = uwp(i) + vis_t(c)*(u % z(c) + w % x(c))
            var_1(i) = var_1(i) + kin % n(c)
            var_2(i) = var_2(i) + eps % n(c) 
            var_3(i) = var_3(i) + 5.64 + 1.41*log(grid % wall_dist(c))
            var_4(i) = var_4(i) + grid % wall_dist(c)
            Rad_mp(i) = Rad_mp(i) + Rad_2 
          end if
        end if
      end do 
    end do 

!---- average over all processors
    do pl=1, Nprob
      call Comm_Mod_Global_Sum_Int(Ncount(pl))
      call Comm_Mod_Global_Sum_Real(Ump(pl))
      call Comm_Mod_Global_Sum_Real(Rad_mp(pl))
      call Comm_Mod_Global_Sum_Real(uwp(pl))
      call Comm_Mod_Global_Sum_Real(var_1(pl))
      call Comm_Mod_Global_Sum_Real(var_2(pl))
      call Comm_Mod_Global_Sum_Real(var_3(pl))
      call Comm_Mod_Global_Sum_Real(var_4(pl))

      count =  count + Ncount(pl) 
    end do
   
    if(k == 1) then
      open(k,FILE=res_name_ht)
    else if(k == 2) then
      open(k,FILE=res_name_cp)
    end if        

    do i = 1, Nprob
      if(Ncount(i) /= 0) then
        Ump(i)    =  Ump(i)/Ncount(i)
        uwp(i)    =  uwp(i)/Ncount(i)
        var_1(i)  =  var_1(i)/Ncount(i)
        var_2(i)  =  var_2(i)/Ncount(i)
        var_3(i)  =  var_3(i)/Ncount(i)
        var_4(i)  =  var_4(i)/Ncount(i)
        Rad_mp(i) =  Rad_mp(i)/Ncount(i)
        
!--------------Write variables ------------------------------------------------------------!

      write(k,'(3E15.7, I7)') var_4(i), (Ump(i) - var_3(i))/var_3(i), &
                                   (var_1(i)/(8.9**2)), Ncount(i)

!------------------------------------------------------------------------------------------!

        Ump(i)    = 0.0 
        uwp(i)    = 0.0 
        var_1(i)  = 0.0 
        var_2(i)  = 0.0 
        var_3(i)  = 0.0
        var_4(i)  = 0.0
        Rad_mp(i) = 0.0 
        
        Ncount(i) = 0
      end if
    end do 

    close(k)
  end do
 
  deallocate(Np)
  deallocate(z_p)
  deallocate(Ump)
  deallocate(uwp)
  deallocate(var_1)
  deallocate(var_2)
  deallocate(var_3)
  deallocate(var_4)
  deallocate(Rad_mp)
  deallocate(Ncount)
  
  call Comm_Mod_Wait
   
  if(this_proc< 2) write(*,*) 'Finished with user function for Askervein hill ! '

  ! Restore the name
  problem_name = store_name

  END SUBROUTINE User_Mod_Cut_lines_ht_cp
