include '../User_Mod/Tecplot_plane.f90'
include '../User_Mod/Cut_lines_ht_cp.f90'

!==============================================================================!
  subroutine User_Mod_Save_Results(flow, save_name)
!------------------------------------------------------------------------------!
!   Calls User_Impinging_Jet_Nu and User_Impinging_Jet_Profile functions.      !
!------------------------------------------------------------------------------!
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  character(len=*) :: save_name
!==============================================================================!

  call User_Mod_Tecplot_plane  (flow, save_name)
  call User_Mod_Cut_lines_ht_cp(flow, save_name)

  end subroutine
