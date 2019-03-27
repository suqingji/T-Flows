!==============================================================================!
   subroutine User_Mod_Tecplot_plane(flow, save_name) 
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for scalar.        !
!   It is called from "Compute_Scalar" function, just before calling the       !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent should the user want to stabilize the         !
!   system for always positive variables, for example.                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod,  only: Field_Type, heat_flux
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
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  real,            pointer :: flux(:)

  INTEGER             :: n
  REAL                :: x1_dis, x0_dis, y1_dis, y0_dis, z1_dis, z0_dis
  REAL                :: TauWup, TauWdown
  REAL    :: Utan, UnorSq, Unor, UtotSq, dely, Stot, R, UtanX
!-----------------------------------[Locals]-----------------------------------!
  INTEGER             :: c, s, c1, c2
  CHARACTER*80        :: namout, res_name, name_out_9
  CHARACTER*39        :: path
  character(len=80)   :: store_name
!------------------------------------------------------------------------------!
  ! Take aliases
  grid => flow % pnt_grid
  flux => flow % flux
  bulk => flow % bulk
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t
!------------------------------------------------------------------------------!

  ! Store the name
  store_name = problem_name
  problem_name = save_name

  call Name_File(this_proc, res_name, "-tec360.dat")
!  call Name_File(this_proc, name_out_9, '-tec360.dat')
!  open(9, FILE='slice.dat')

!  if(this_proc < 2) write(*,*) '# Now reading the file: slice.dat '  
!  read(9,*) z1_dis, z0_dis
!  if(this_proc < 2) write(*,*) '# X:[ ', x0_dis, " ;", x1_dis, "]", &
!  ' Y:[ ', y0_dis, " ;", y1_dis, "]", ' Z:[ ', z0_dis, " ;", z1_dis, "]"

!    inflowfile = 'tecplt-10m-xxxxxx-xx.dat'

!    write(inflowfile(12:17),'(I6.6)') n
!    write(inflowfile(19:20),'(I2.2)') this_proc
  open(500+this_proc, FILE=res_name)

  if( this_proc < 2 ) then
    write(*,*) 'Capturing field..  ', res_name
  end if


  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or. &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
  
        write (500+this_proc,'(8E17.7E3)') grid % xc(c1), grid % yc(c1),&
            grid % zc(c1), U % n(c1), V % n(c1), W % n(c1), Kin % n(c1), &
            sqrt( (U % n(c1))**2+ (V % n(c1))**2+ (W % n(c1))**2)
                 
      endif
    endif
  end do

  call Comm_Mod_Wait
   
  close(9)
  close(500+this_proc)
    
  call Comm_Mod_Wait
   
  if ( this_proc < 2 ) then
    write(*,*) 'It is done'
  end if

  ! Restore the name
  problem_name = store_name

  END SUBROUTINE 
