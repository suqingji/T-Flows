!==============================================================================!
  subroutine Control_Mod_Save_Results_At_Boundaries(save_results_bnd, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: save_results_bnd
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('SAVE_RESULTS_AT_BOUNDARIES', 'yes',  &
                                   val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    save_results_bnd = .true.

  else if( val .eq. 'NO' ) then
    save_results_bnd = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for save_results_bnd: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
