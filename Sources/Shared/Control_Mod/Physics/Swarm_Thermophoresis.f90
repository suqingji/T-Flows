!==============================================================================!
  subroutine Control_Mod_Swarm_Thermophoresis(thermophoresis, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: thermophoresis
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val 
!==============================================================================!

  call Control_Mod_Read_Char_Item('THERMOPHORESIS', 'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    thermophoresis = .true.

  else if( val .eq. 'NO' ) then
    thermophoresis = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for thermophoresis: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
