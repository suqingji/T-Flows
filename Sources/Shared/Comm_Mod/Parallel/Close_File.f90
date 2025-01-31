!==============================================================================!
  subroutine Close_File(Comm, fh)
!------------------------------------------------------------------------------!
!   Close file for parallel runs.                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Comm
  integer          :: fh  ! file handle
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  ! Close the file
  call Mpi_File_Close(fh, error)

  end subroutine
