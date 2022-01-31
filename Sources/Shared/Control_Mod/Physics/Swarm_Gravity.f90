!==============================================================================!
  subroutine Control_Mod_Swarm_Gravity(swarm_gravity, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical           :: swarm_gravity
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val 
!==============================================================================!

  call Control_Mod_Read_Char_Item('SWARM_GRAVITY', 'no', val, verbose)
  call To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    swarm_gravity = .true.

  else if( val .eq. 'NO' ) then
    swarm_gravity = .false.

  else
    if(this_proc < 2) then
      print *, '# ERROR!  Unknown state for swarm gravity: ', trim(val)
      print *, '# Exiting!'
    end if
    call Comm_Mod_End

  end if

  end subroutine
