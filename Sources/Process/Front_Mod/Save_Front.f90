!==============================================================================!
  subroutine Save_Front(Front, time_step)
!------------------------------------------------------------------------------!
!   Writes surface vertices in VTU file format (for VisIt and Paraview)        !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Front_Type), target :: Front
  integer                   :: time_step
!----------------------------------[Locals]------------------------------------!
  type(Vert_Type), pointer :: Vert
  integer                  :: v, e     ! vertex and element counters
  integer                  :: offset, fu
  character(SL)            :: name_out
!-----------------------------[Local parameters]-------------------------------!
  integer, parameter :: VTK_TRIANGLE = 5  ! cell shapes in VTK format
  character(len= 0)  :: IN_0 = ''         ! indentation levels
  character(len= 2)  :: IN_1 = '  '
  character(len= 4)  :: IN_2 = '    '
  character(len= 6)  :: IN_3 = '      '
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
!==============================================================================!

  if(Front % n_verts < 1) return

  !----------------------------!
  !                            !
  !   Create .front.vtu file   !
  !                            !
  !----------------------------!

! if(this_proc < 2) then

    call File % Set_Name(name_out,               &
                         processor = this_proc,  &
                         time_step = time_step,  &
                         appendix  = '-front',   &
                         extension = '.vtu')
    call File % Open_For_Writing_Ascii(name_out, fu)

    !------------!
    !            !
    !   Header   !
    !            !
    !------------!
    write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" '//&
                            'byte_order="LittleEndian">'
    write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'

    write(fu,'(a,a,i0.0,a,i0.0,a)')   &
                IN_2, '<Piece NumberOfPoints="', Front % n_verts,  &
                           '" NumberOfCells ="', Front % n_elems, '">'

    !------------------------!
    !                        !
    !   Vertex coordinates   !
    !                        !
    !------------------------!
    write(fu,'(a,a)') IN_3, '<Points>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents' //  &
                            '="3" format="ascii">'
    do v = 1, Front % n_verts
      Vert => Front % Vert(v)
      write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
                  IN_5, Vert % x_n, Vert % y_n, Vert % z_n
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_3, '</Points>'

    !----------------!
    !                !
    !   Point data   !
    !                !
    !----------------!
    write(fu,'(a,a)') IN_3, '<PointData Scalars="scalars" vectors="velocity">'

    !--------------------!
    !   Particle i.d.s   !
    !--------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="Index" ' // &
                            'format="ascii">'
    do v = 1, Front % n_verts
      write(fu,'(a,i9)') IN_5, v
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !--------------------------!
    !   Number of neighbours   !
    !--------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="Neighbours" ' // &
                            'format="ascii">'
    do v = 1, Front % n_verts
      write(fu,'(a,i9)') IN_5, Front % Vert(v) % nne
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-----------------------------!
    !   Curvatures at the nodes   !
    !-----------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Float64" Name="NodeCurv" ' // &
                           ' format="ascii">'
    do v = 1, Front % n_verts
      Vert => Front % Vert(v)
      write(fu,'(a,1pe16.6e4)') IN_5, Vert % curv
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    write(fu,'(a,a)') IN_3, '</PointData>'

    !-----------!
    !           !
    !   Cells   !
    !           !
    !-----------!
    write(fu,'(a,a)') IN_3, '<Cells>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                            ' format="ascii">'
    ! Cell topology
    do e = 1, Front % n_elems
      write(fu,'(a,99i9)') IN_5, Front % Elem(e) % v(1:Front % Elem(e) % nv)-1
    end do

    ! Cell offsets
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets"' //  &
                            ' format="ascii">'
    offset = 0
    do e = 1, Front % n_elems
      offset = offset + Front % Elem(e) % nv
      write(fu,'(a,i9)') IN_5, offset
    end do

    ! Cell types
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="types"' //  &
                            ' format="ascii">'
    do e = 1, Front % n_elems
      write(fu,'(a,i9)') IN_5, VTK_POLYGON
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_3, '</Cells>'

    !---------------!
    !               !
    !   Cell data   !
    !               !
    !---------------!

    ! Beginning of cell data
    write(fu,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

    !-------------------------------------!
    !   Number of neighbouring elements   !
    !-------------------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="Neighbours"' //  &
                            ' format="ascii">'
    do e = 1, Front % n_elems
      write(fu,'(a,i9)') IN_5, Front % Elem(e) % nne
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !---------------------!
    !   Surface normals   !
    !---------------------!
    write(fu,'(4a)') IN_4,                                                &
                   '<DataArray type="Float64" Name="ElementNormals" ' //  &
                   ' NumberOfComponents="3" format="ascii">'
    do e = 1, Front % n_elems
      write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
                IN_5, Front % Elem(e) % nx,           &
                      Front % Elem(e) % ny,           &
                      Front % Elem(e) % nz
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------!
    !   Element areas   !
    !-------------------!
    write(fu,'(4a)') IN_4,                                                &
                   '<DataArray type="Float64" Name="ElementArea" ' //  &
                   ' format="ascii">'
    do e = 1, Front % n_elems
      write(fu,'(a,1pe16.6e4)') IN_5, Front % Elem(e) % area
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------------!
    !   Element coordinates   !
    !-------------------------!
    write(fu,'(4a)') IN_4,                                                    &
                   '<DataArray type="Float64" Name="ElementCoordinates" ' //  &
                   ' NumberOfComponents="3" format="ascii">'
    do e = 1, Front % n_elems
      write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
                IN_5, Front % Elem(e) % xe,           &
                      Front % Elem(e) % ye,           &
                      Front % Elem(e) % ze
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !------------------------!
    !   Surface curvatures   !
    !------------------------!
    write(fu,'(4a)') IN_4,                                                &
                   '<DataArray type="Float64" Name="ElementCurv" ' //  &
                   ' format="ascii">'
    do e = 1, Front % n_elems
      write(fu,'(a,1pe16.6e4)') IN_5, Front % Elem(e) % curv
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    ! End of cell data
    write(fu,'(a,a)') IN_3, '</CellData>'

    !------------!
    !            !
    !   Footer   !
    !            !
    !------------!
    write(fu,'(a,a)') IN_2, '</Piece>'
    write(fu,'(a,a)') IN_1, '</UnstructuredGrid>'
    write(fu,'(a,a)') IN_0, '</VTKFile>'
    close(fu)
! end if

  end subroutine
