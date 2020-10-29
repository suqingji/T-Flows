!==============================================================================!
  subroutine Control_Mod_Least_Squares_Curvature(least_squares_curvature,   &
                                                 verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: least_squares_curvature
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item( 'LEAST_SQUARES_CURVATURE', 'yes',    &
                                   val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    least_squares_curvature = .true.

  else if( val .eq. 'NO' ) then
    least_squares_curvature = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for curvature calculation method: ', &
                trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
