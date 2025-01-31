!==============================================================================!
  module Elem_Mod
!------------------------------------------------------------------------------!
!   Storage for Elem_Type used by Front_Mod and Surf_mod                       !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !   Elem type   !
  !---------------!
  type Elem_Type

    integer :: nne         ! number of neighbouring element
    integer :: nv          ! number of vertices for the element
    integer :: ns          ! number of sides for the element
    integer :: v(24)
    integer :: e(24)
    integer :: s(24)
    integer :: cell        ! cell at which the surface resides
    integer :: face        ! face which the surface intersects
    real    :: nx, ny, nz  ! surface normal vector
    real    :: xc, yc, zc  ! center of a sphere
    real    :: xe, ye, ze  ! center of element
    real    :: sx, sy, sz  ! surface area components
    real    :: area
    real    :: curv

    contains
      procedure :: Initialize_Elem

  end type

  contains

  include 'Elem_Mod/Initialize_Elem.f90'

  end module
