!==============================================================================!
  subroutine Print_Surf_Statistics(Surf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type),  target :: Surf
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: elem(:)
  type(Side_Type), pointer :: side(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: item, i, j, k, c, d, e, s, v, si, sj, sk
  integer                  :: nne_s, nne_e
  real, allocatable        :: nne(:)
  real                     :: a(3), b(3), tri_v(3)
  real                     :: max_rat, min_rat, max_l, min_l
  character(len=160)       :: line
!-----------------------------[Local parameters]-------------------------------!
  integer, parameter :: T=33  ! indent
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ns   => Surf % n_sides
  ne   => Surf % n_elems
  Vert => Surf % Vert
  side => Surf % side
  elem => Surf % elem

  !--------------------------!
  !   Compute side lengths   !
  !--------------------------!
  do s = 1, ns
    c = side(s) % c
    d = side(s) % d
    side(s) % length = Math_Mod_Distance(                 &
            Vert(c) % x_n, Vert(c) % y_n, Vert(c) % z_n,  &
            Vert(d) % x_n, Vert(d) % y_n, Vert(d) % z_n)
  end do

  !-----------------------------!
  !   Compute elements' areas   !
  !-----------------------------!
  do e = 1, ne
    i = elem(e) % v(1)
    j = elem(e) % v(2)
    k = elem(e) % v(3)
    a(1) = Vert(j) % x_n - Vert(i) % x_n
    a(2) = Vert(j) % y_n - Vert(i) % y_n
    a(3) = Vert(j) % z_n - Vert(i) % z_n
    b(1) = Vert(k) % x_n - Vert(i) % x_n
    b(2) = Vert(k) % y_n - Vert(i) % y_n
    b(3) = Vert(k) % z_n - Vert(i) % z_n
    tri_v = Math_Mod_Cross_Product(a, b)
    elem(e) % area = sqrt(dot_product(tri_v, tri_v)) * 0.5
  end do

  !------------------------------!
  !   Extreme side size ratios   !
  !------------------------------!
  max_rat = -HUGE
  min_rat = +HUGE
  do e = 1, ne
    si = elem(e) % s(1)
    sj = elem(e) % s(2)
    sk = elem(e) % s(3)
    max_l = max(side(si) % length, side(sj) % length, side(sk) % length)
    min_l = min(side(si) % length, side(sj) % length, side(sk) % length)
    max_rat = max(max_rat, max_l/min_l)
    min_rat = min(min_rat, max_l/min_l)
  end do

  !--------------------------------!
  !   Count number of neighbours   !
  !--------------------------------!
  call Surf % Find_Vertex_Elements()
  nne_s = minval(Vert(1:nv) % nne)
  nne_e = maxval(Vert(1:nv) % nne)
  allocate(nne(nne_s:nne_e)); nne = 0.0
  do v = 1, nv
    nne(Vert(v) % nne) = nne(Vert(v) % nne) + 1.0
  end do

  if(this_proc < 2) then

    line( 1:160) = ' '
    line( 1+T:63+T) =   &
               '#=============================================================#'
    print *, trim(line)
    line( 1+T:63+T) =   &
               '#                   Surface mesh statistics                   #'
    print *, trim(line)
    line( 1+T:63+T) =   &
               '#-------------------------------------------------------------#'
    print *, trim(line)

    ! Number of elements (1), vertices (2) and sides (3)
    do item = 1, 3
      line( 1:160) = ' '
      line( 1+T: 1+T) = '#'
      line(63+T:63+T) = '#'
      line( 3+T: 3+T) = '-'
      if(item.eq.1) line( 5+T: 5+T+20) = 'Number of elements: '
      if(item.eq.1) write(line(32+T:37+T), '(i6)') Surf % n_elems
      if(item.eq.2) line( 5+T: 5+T+20) = 'Number of vertices: '
      if(item.eq.2) write(line(32+T:37+T), '(i6)') Surf % n_verts
      if(item.eq.3) line( 5+T: 5+T+20) = 'Number of sides:    '
      if(item.eq.3) write(line(32+T:37+T), '(i6)') Surf % n_sides
      print *, trim(line)
    end do

    line( 1+T:63+T) =  &
               '#-------------------------------------------------------------#'
    print *, trim(line)

    ! Maximum (1) and minimum (2) element area
    ! Maximum (3) and minimum (4) side length
    ! Maximum (5) and minimum (6) side length ratio
    do item = 1, 6
      line( 1:160) = ' '
      line( 1+T: 1+T) = '#'
      line(63+T:63+T) = '#'
      line( 3+T: 3+T) = '-'
      if(item.eq.1) line( 5+T: 5+T+33) = 'Maximum element area:          '
      if(item.eq.1) write(line(36+T:47+T), '(1pe12.5)')  &
                          maxval(elem(1:ne) % area)
      if(item.eq.2) line( 5+T: 5+T+33) = 'Minimum element area:          '
      if(item.eq.2) write(line(36+T:47+T), '(1pe12.5)')  &
                          minval(elem(1:ne) % area)
      if(item.eq.3) line( 5+T: 5+T+33) = 'Maximum side length:           '
      if(item.eq.3) write(line(36+T:47+T), '(1pe12.5)')  &
                          maxval(side(1:ns) % length)
      if(item.eq.4) line( 5+T: 5+T+33) = 'Minimum side length:           '
      if(item.eq.4) write(line(36+T:47+T), '(1pe12.5)')  &
                          minval(side(1:ns) % length)
      if(item.eq.5) line( 5+T: 5+T+33) = 'Maximum side ratio in element: '
      if(item.eq.5) write(line(36+T:47+T), '(1pe12.5)') max_rat
      if(item.eq.6) line( 5+T: 5+T+33) = 'Minimum side ratio in element: '
      if(item.eq.6) write(line(36+T:47+T), '(1pe12.5)') min_rat
      print *, trim(line)
    end do

    line( 1+T:63+T) =  &
               '#-------------------------------------------------------------#'
    print *, trim(line)

    ! Number of neighbours
    do item = nne_s, nne_e
      line( 1:160) = ' '
      line( 1+T: 1+T) = '#'
      line(63+T:63+T) = '#'
      line( 3+T: 3+T) = '-'
      line( 5+T: 5+T+43) = 'Percentage of vertices with XX neighbours: '
      write(line(33+T:34+T), '(i2)') item
      write(line(48+T:53+T), '(f6.2)') nne(item) / nv * 100.0
      line(55+T:55+T) = '%'
      print *, trim(line)
    end do

    line( 1+T:63+T) =  &
               '#-------------------------------------------------------------#'
    print *, trim(line)
    print *, ''

  end if

  end subroutine
