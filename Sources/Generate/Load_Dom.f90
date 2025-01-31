!==============================================================================!
  subroutine Load_Dom(dom, smr, ref, Grid)
!------------------------------------------------------------------------------!
!   Reads: .dom file                                                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Domain_Mod
  use Smooths_Mod
  use Refines_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Domain_Type)  :: dom
  type(Smooths_Type) :: smr
  type(Refines_Type) :: ref
  type(Grid_Type)    :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer       :: b, i, l, s, i_fac, n, n1, n2, n3, n4, dumi, fu
  integer       :: n_faces_check, n_nodes_check
  integer       :: ni, nj, nk, npnt, nsurf
  character(SL) :: dum
  character(SL) :: domain_name
  character(SL) :: answer
  real          :: xt(8), yt(8), zt(8)
  integer       :: fn(6,4)
!------------------------------------------------------------------------------!
  include 'Block_Numbering.f90'
!==============================================================================!

  ! Copy face-node numbering for blocks
  fn = hex_block

  !---------------------------------------------!
  !   The Grid will not have polyhedral cells   !
  !---------------------------------------------!
  Grid % polyhedral = .false.

  print *, '#========================================'
  print *, '# Input problem name: (without extension)'
  print *, '#----------------------------------------'
  call File % Read_Line(5)
  read(line % tokens(1), *) problem_name(1)

  Grid % name = problem_name(1)
  call To_Upper_Case(Grid % name)

  call File % Set_Name(domain_name, extension='.dom')
  call File % Open_For_Reading_Ascii(domain_name, fu)

  !-----------------------------------------------------------------!
  !   Max. number of nodes (cells), boundary faces and cell faces   !
  !-----------------------------------------------------------------!
  call File % Read_Line(fu)
  read(line % tokens(1), *) Grid % max_n_nodes
  read(line % tokens(2), *) Grid % max_n_bnd_cells
  read(line % tokens(3), *) Grid % max_n_faces

  !---------------------!
  !   Allocate memory   !
  !---------------------!
  print *, '# Allocating memory for: '
  print *, '#', Grid % max_n_nodes,     ' nodes and cells'
  print *, '#', Grid % max_n_bnd_cells, ' boundary cells'
  print *, '#', Grid % max_n_faces,     ' cell faces'

  allocate (Grid % bnd_cond % color(-Grid % max_n_bnd_cells-1:-1))
  Grid % bnd_cond % color = 0

  ! Variables in Grid_Mod
  call Grid % Allocate_Nodes(Grid % max_n_nodes)

  call Grid % Allocate_Cells(Grid % max_n_nodes,      &
                             Grid % max_n_bnd_cells)

  call Grid % Allocate_Faces(Grid % max_n_faces, 0)

  call Refines_Mod_Allocate_Cells(ref,                     &
                                  Grid % max_n_bnd_cells,  &
                                  Grid % max_n_nodes)

  ! Variables still declared in Gen_Mod.h90:
  allocate (face_c_to_c(Grid % max_n_faces,2))
  face_c_to_c = 0
  allocate (twin_n(Grid % max_n_nodes,0:8))
  twin_n (:,:) = 0

  print *, '# Allocation successfull !'

  !-------------!
  !   Corners   !
  !-------------!
  call File % Read_Line(fu)
  read(line % tokens(1), *) dom % n_points  ! number of points

  call Domain_Mod_Allocate_Points(dom, dom % n_points)

  do i = 1, dom % n_points
    call File % Read_Line(fu)
    read(line % tokens(2),*) dom % points(i) % x
    read(line % tokens(3),*) dom % points(i) % y
    read(line % tokens(4),*) dom % points(i) % z
  end do

  !------------!
  !   Blocks   !
  !------------!
  call File % Read_Line(fu)
  read(line % tokens(1), *) dom % n_blocks  ! number of blocks

  call Domain_Mod_Allocate_Blocks(dom, dom % n_blocks)

  ! Initialize weights
  do b=1, dom % n_blocks
    dom % blocks(b) % weights      = 1.0
    dom % blocks(b) % face_weights = 1.0
  end do

  do b = 1, dom % n_blocks
    dom % blocks(b) % corners(0)=1       ! suppose it is properly oriented

    call File % Read_Line(fu)
    read(line % tokens(2),*) dom % blocks(b) % resolutions(1)
    read(line % tokens(3),*) dom % blocks(b) % resolutions(2)
    read(line % tokens(4),*) dom % blocks(b) % resolutions(3)

    call File % Read_Line(fu)
    read(line % whole, *)               &  ! block weights
         dom % blocks(b) % weights(1),  &
         dom % blocks(b) % weights(2),  &
         dom % blocks(b) % weights(3)

    call File % Read_Line(fu)
    read(line % whole, *)                                             &
         dom % blocks(b) % corners(1), dom % blocks(b) % corners(2),  &
         dom % blocks(b) % corners(3), dom % blocks(b) % corners(4),  &
         dom % blocks(b) % corners(5), dom % blocks(b) % corners(6),  &
         dom % blocks(b) % corners(7), dom % blocks(b) % corners(8)

    !---------------------------!
    !   Check if the block is   !
    !     properly oriented     !
    !---------------------------!
    do n=1,8
      xt(n) = dom % points(dom % blocks(b) % corners(n)) % x
      yt(n) = dom % points(dom % blocks(b) % corners(n)) % y
      zt(n) = dom % points(dom % blocks(b) % corners(n)) % z
    end do

    if(Math % Tet_Volume( xt(2),yt(2),zt(2), xt(5),yt(5),zt(5),  &
                          xt(3),yt(3),zt(3), xt(1),yt(1),zt(1) )  < 0) then
      dom % blocks(b) % corners(0)=-1            !  It's nor properly oriented
      call Swap_Int(dom % blocks(b) % corners(2),  &
                    dom % blocks(b) % corners(3))
      call Swap_Int(dom % blocks(b) % corners(6),  &
                    dom % blocks(b) % corners(7))
      call Swap_Real(dom % blocks(b) % weights(1),  &
                     dom % blocks(b) % weights(2))
      dom % blocks(b) % weights(1) = 1.0 / dom % blocks(b) % weights(1)
      dom % blocks(b) % weights(2) = 1.0 / dom % blocks(b) % weights(2)
      call Swap_Int(dom % blocks(b) % resolutions(1),  &
                    dom % blocks(b) % resolutions(2))
      print *, 'Warning: Block ',b,' was not properly oriented'
    end if
  end do                 ! through dom % blocks

  !-----------------------------!
  !   Set the corners of each   !
  !      face of the block      !
  !-----------------------------!
  do b = 1, dom % n_blocks
    do i_fac = 1, 6
      do n = 1, 4
        dom % blocks(b) % faces(i_fac, n) =  &
        dom % blocks(b) % corners(fn(i_fac,n))
      end do
    end do
  end do

  !----------------------------------------------!
  !   Lines                                      !
  !----------------------------------------------!
  !   Lines can be prescribed point by point     !
  !   or with just a weighting factor.           !
  !----------------------------------------------!
  call File % Read_Line(fu)
  read(line % tokens(1), *) dom % n_lines  ! number of defined dom % lines

  call Domain_Mod_Allocate_Lines(dom, dom % n_lines)

  do l=1, dom % n_lines
    call File % Read_Line(fu)

    read(line % tokens(1),*) npnt
    read(line % tokens(2),*) dom % lines(l) % points(1)
    read(line % tokens(3),*) dom % lines(l) % points(2)

    call Domain_Mod_Find_Line(dom,                         &
                              dom % lines(l) % points(1),  &
                              dom % lines(l) % points(2),  &
                              dom % lines(l) % resolution)

    ! Does this need a more elegant solution?
    allocate(dom % lines(l) % x( dom % lines(l) % resolution ))
    allocate(dom % lines(l) % y( dom % lines(l) % resolution ))
    allocate(dom % lines(l) % z( dom % lines(l) % resolution ))

    ! Point by point
    if(npnt > 0) then
      do n=1,dom % lines(l) % resolution
        call File % Read_Line(fu)
        read(line % tokens(2),*) dom % lines(l) % x(n)
        read(line % tokens(3),*) dom % lines(l) % y(n)
        read(line % tokens(4),*) dom % lines(l) % z(n)
      end do

    ! Weight factor
    else
      call File % Read_Line(fu)
      read(line % tokens(1), *) dom % lines(l) % weight
    end if

  end do

  !----------------------------------------!
  !   Copy block weights to face weights   !
  !----------------------------------------!
  do b = 1, dom % n_blocks
    do i_fac = 1,6                          !  face of the block
      dom % blocks(b) % face_weights(i_fac, 1) = dom % blocks(b) % weights(1)
      dom % blocks(b) % face_weights(i_fac, 2) = dom % blocks(b) % weights(2)
      dom % blocks(b) % face_weights(i_fac, 3) = dom % blocks(b) % weights(3)
    end do
  end do

  !--------------!
  !   Surfaces   !
  !--------------!
  call File % Read_Line(fu)
  read(line % tokens(1), *) nsurf     ! number of defined surfaces

  do s = 1, nsurf
    call File % Read_Line(fu)
    read(line % whole,*) dum, n1, n2, n3, n4
    call Domain_Mod_Find_Surface(dom, n1, n2, n3, n4, b, i_fac)
    print *, '# block: ', b, ' surf: ', i_fac
    n = (b-1)*6 + i_fac         ! surface number

    call File % Read_Line(fu)
    read(line % whole, *) dom % blocks(b) % face_weights(i_fac,1),  &
                          dom % blocks(b) % face_weights(i_fac,2),  &
                          dom % blocks(b) % face_weights(i_fac,2)
  end do

  !---------------------------------------!
  !   Is there enough allocated memory?   !
  !---------------------------------------!

  ! Nodes & faces
  n_nodes_check = 0
  n_faces_check = 0
  do b = 1, dom % n_blocks
    ni = dom % blocks(b) % resolutions(1)
    nj = dom % blocks(b) % resolutions(2)
    nk = dom % blocks(b) % resolutions(3)
    n_nodes_check=n_nodes_check + ni*nj*nk
    n_faces_check=n_faces_check + ni*nj*nk + 2*( (ni*nj)+(nj*nk)+(ni*nk) )
  end do

  if( (n_faces_check > Grid % max_n_faces) .or.  &
      (n_nodes_check > Grid % max_n_nodes) ) then
    print *, '# Error message from T-Flows:'
  end if

  if( n_faces_check  > Grid % max_n_faces ) then
    print *, '# The estimated number of faces is :', n_faces_check
    print *, '# There is space available only for:', Grid % max_n_faces
    print *, '# Increase the parameter Grid % max_n_faces in the input file'
    print *, '# and re-run the code !'
  end if

  if( n_nodes_check  > Grid % max_n_nodes ) then
    print *, '# The estimated number of nodes is :', n_nodes_check
    print *, '# There is space available only for:', Grid % max_n_nodes
    print *, '# Increase the parameter Grid % max_n_nodes in the input file'
    print *, '# and re-run the code !'
  end if

  if( (n_faces_check > Grid % max_n_faces) .or.  &
      (n_nodes_check > Grid % max_n_nodes) ) then
    stop
  end if

  !---------------------------------------!
  !   Boundary conditions and materials   !
  !---------------------------------------!
  call File % Read_Line(fu)
  read(line % tokens(1), *) dom % n_regions  ! number of regions (can be bnd.
                                             ! conditions or materials)

  call Domain_Mod_Allocate_Regions(dom, dom % n_regions)

  do n = 1, dom % n_regions
    dom % regions(n) % face=''

    call File % Read_Line(fu)
    if(line % n_tokens .eq. 7) then
      read(line % whole,*)  dum,            &
                   dom % regions(n) % is,   &
                   dom % regions(n) % js,   &
                   dom % regions(n) % ks,   &
                   dom % regions(n) % ie,   &
                   dom % regions(n) % je,   &
                   dom % regions(n) % ke
    else if(line % n_tokens .eq. 2) then
      read(line % tokens(1),*)       dum
      read(line % tokens(2),'(A4)')  &
           dom % regions(n) % face
      call To_Upper_Case(dom % regions(n) % face)
    end if

    call File % Read_Line(fu)
    read(line % tokens(1), *) dom % regions(n) % block
    read(line % tokens(2), *) dom % regions(n) % name
    call To_Upper_Case(dom % regions(n) % name)

    ! if( dom % blocks(b_cond(n,7)) % points(0) .eq. -1 ) then
    !   call Swap_Int( dom % regions(n) % is,dom % regions(n) % js )
    !   call Swap_Int( dom % regions(n) % ie,dom % regions(n) % je )
    ! end if

  end do

  !-------------------------!
  !   Periodic boundaries   !
  !-------------------------!
  call File % Read_Line(fu)
  read(line % tokens(1), *)  n_periodic_cond  ! number of periodic boundaries
  print '(a38,i7)', '# Number of periodic boundaries:     ', n_periodic_cond

  allocate (periodic_cond(n_periodic_cond,8))

  do n=1,n_periodic_cond
    call File % Read_Line(fu)
    read(line % whole, *) dum, periodic_cond(n,1), periodic_cond(n,2),  &
                               periodic_cond(n,3), periodic_cond(n,4)
    call File % Read_Line(fu)
    read(line % whole, *)      periodic_cond(n,5), periodic_cond(n,6),  &
                               periodic_cond(n,7), periodic_cond(n,8)
  end do

  !---------------------!
  !   Copy boundaries   !
  !---------------------!
  call File % Read_Line(fu)
  read(line % tokens(1), *)  dumi  ! used to be number of copy boundaries

  !-----------------------------------!
  !   Refinement levels and regions   !
  !-----------------------------------!
  call File % Read_Line(fu)
  read(line % tokens(1), *) ref % n_levels     ! number of refinement levels
  print '(a38,i7)', '# Number of refinement levels:       ', ref % n_levels

  call Refines_Mod_Allocate_Levels(ref, ref % n_levels)

  do l = 1, ref % n_levels
    print *, '# Level: ', l
    call File % Read_Line(fu)
    read(line % tokens(2), *) ref % n_regions(l)

    ! Browse through regions in level "l"
    do n = 1, ref % n_regions(l)
      call File % Read_Line(fu)
      read(line % tokens(3),*) answer
      call To_Upper_Case(answer)
      ref % region(l,n) % shape = -1
      if(answer .eq. 'RECTANGLE') ref % region(l,n) % shape = RECTANGLE
      if(answer .eq. 'ELIPSOID')  ref % region(l,n) % shape = ELIPSOID
      if(answer .eq. 'PLANE')     ref % region(l,n) % shape = PLANE
      if(ref % region(l,n) % shape .eq. -1) then
        print *, 'ERROR!  Refinement shape not specified well by: ', answer
        stop
      end if

      call File % Read_Line(fu)
      read(line % whole, *)                      &
                ref % region(l,n) % pnt(1) % x,  &
                ref % region(l,n) % pnt(1) % y,  &
                ref % region(l,n) % pnt(1) % z,  &
                ref % region(l,n) % pnt(2) % x,  &
                ref % region(l,n) % pnt(2) % y,  &
                ref % region(l,n) % pnt(2) % z
    end do
  end do

  !-----------------------!
  !   Smoothing regions   !
  !-----------------------!
  call File % Read_Line(fu)
  read(line % tokens(1), *) smr % n_smooths  ! number of smoothing regions

  print '(a38,i7)', '# Number of (non)smoothing regions:  ', smr % n_smooths

  call Smooths_Mod_Allocate_Smooths(smr, smr % n_smooths)

  do n = 1, smr % n_smooths
    smr % in_x(n) = .false.
    smr % in_y(n) = .false.
    smr % in_z(n) = .false.
    call File % Read_Line(fu)
    read(line % tokens(1), *) dumi    ! this line is probably not needed
    if(line % n_tokens .eq. 4) then   ! smoothing in three directions
      smr % in_x(n) = .true.
      smr % in_y(n) = .true.
      smr % in_z(n) = .true.
    else if(line % n_tokens .eq. 3) then
      call To_Upper_Case(line % tokens(2))
      call To_Upper_Case(line % tokens(3))
      if( line % tokens(2) .eq. 'X' )  smr % in_x(n) = .true.
      if( line % tokens(3) .eq. 'X' )  smr % in_x(n) = .true.
      if( line % tokens(2) .eq. 'Y' )  smr % in_y(n) = .true.
      if( line % tokens(3) .eq. 'Y' )  smr % in_y(n) = .true.
      if( line % tokens(2) .eq. 'Z' )  smr % in_z(n) = .true.
      if( line % tokens(3) .eq. 'Z' )  smr % in_z(n) = .true.
    else if(line % n_tokens .eq. 2) then
      call To_Upper_Case(line % tokens(2))
      if( line % tokens(2) .eq. 'X' )  smr % in_x(n) = .true.
      if( line % tokens(2) .eq. 'Y' )  smr % in_y(n) = .true.
      if( line % tokens(2) .eq. 'Z' )  smr % in_z(n) = .true.
    end if

    ! Read the coordinates of the (non)smoothed region
    call File % Read_Line(fu)
    read(line % whole, *) smr % iters(n),  &
                          smr % relax(n)

    call File % Read_Line(fu)
    read(line % whole, *) smr % x_min(n),   &
                          smr % y_min(n),   &
                          smr % z_min(n),   &
                          smr % x_max(n),   &
                          smr % y_max(n),   &
                          smr % z_max(n)
  end do

  close(fu)

  end subroutine
