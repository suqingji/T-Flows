!==============================================================================!
  subroutine Logo_Gen()
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!==============================================================================!

  print *,'#===================================' // &
          '===================================='
  print *,'#                                                                   '
  print *,'#    ______________________ ____    ________  __      __  _________ '
  print *,'#    \__    ___/\_   _____/|    |   \_____  \/  \    /  \/   _____/ '
  print *,'#      |    |    |    __)  |    |    /   |   \   \/\/   /\_____  \  '
  print *,'#      |    |    |     \   |    |___/    |    \        / /        \ '
  print *,'#      |____|    \___  /   |_______ \_______  /\__/\  / /_______  / '
  print *,'#                    \/            \/       \/      \/          \/  '
  print *,'#                   _____                      __                   '
  print *,'#                  / ___/__ ___  ___ _______ _/ /____               '
  print *,'#                 / (_ / -_) _ \/ -_) __/ _ `/ __/ -_)              '
  print *,'#                 \___/\__/_//_/\__/_/  \_,_/\__/\__/               '
  print *,'#                                                                   '
  if(RP .eq. DP) then
  print *,'#                        Double precision mode'
  else
  print *,'#                        Single precision mode'
  end if
  print *,'#-----------------------------------' // &
          '------------------------------------'

  end subroutine

