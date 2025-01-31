!==============================================================================!
  subroutine Swap_Side(Surf, s)
!------------------------------------------------------------------------------!
!  This function calculates radii of inscribed and circumscribed circle        !
!  for a given element (int e)                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  integer                  :: s
!-----------------------------------[Locals]-----------------------------------!
  integer                  :: a, b, c, d, ea, eb
  integer                  :: eac, ead, ebc, ebd
  integer                  :: sad, sac, sbc, sbd
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  type(Side_Type), pointer :: side(:)
!==============================================================================!

  ! Take aliases
  Vert => Surf % Vert
  Elem => Surf % Elem
  side => Surf % side

  ea = side(s) % ea
  eb = side(s) % eb
  a  = side(s) % a
  b  = side(s) % b
  c  = side(s) % c
  d  = side(s) % d

  eac = 0; ead = 0; ebc = 0; ebd = 0
  sad = 0; sac = 0; sbc = 0; sbd = 0

  if(Elem(ea) % e(1) .eq. eb) then
    ead = Elem(ea) % e(2);  eac = Elem(ea) % e(3)
    sad = Elem(ea) % s(2);  sac = Elem(ea) % s(3)
  end if
  if(Elem(ea) % e(2) .eq. eb) then
    ead = Elem(ea) % e(3);  eac = Elem(ea) % e(1)
    sad = Elem(ea) % s(3);  sac = Elem(ea) % s(1)
  end if
  if(Elem(ea) % e(3) .eq. eb) then
    ead = Elem(ea) % e(1);  eac = Elem(ea) % e(2)
    sad = Elem(ea) % s(1);  sac = Elem(ea) % s(2)
  end if

  if(Elem(eb) % e(1) .eq. ea) then
    ebc = Elem(eb) % e(2);  ebd = Elem(eb) % e(3)
    sbc = Elem(eb) % s(2);  sbd = Elem(eb) % s(3)
  end if
  if(Elem(eb) % e(2) .eq. ea) then
    ebc = Elem(eb) % e(3);  ebd = Elem(eb) % e(1)
    sbc = Elem(eb) % s(3);  sbd = Elem(eb) % s(1)
  end if
  if(Elem(eb) % e(3) .eq. ea) then
    ebc = Elem(eb) % e(1);  ebd = Elem(eb) % e(2)
    sbc = Elem(eb) % s(1);  sbd = Elem(eb) % s(2)
  end if

  !--------------------------------------------!
  !   Change the orientation of the elements   !
  !--------------------------------------------!
  Elem(ea) % v(1) = a;    Elem(ea) % v(2) = b;    Elem(ea) % v(3) = d
  Elem(ea) % e(1) = ebd;  Elem(ea) % e(2) = ead;  Elem(ea) % e(3) = eb
  Elem(ea) % s(1) = sbd;  Elem(ea) % s(2) = sad;  Elem(ea) % s(3) = s

  Elem(eb) % v(1) = a;    Elem(eb) % v(2) = c;    Elem(eb) % v(3) = b
  Elem(eb) % e(1) = ebc;  Elem(eb) % e(2) = ea;   Elem(eb) % e(3) = eac
  Elem(eb) % s(1) = sbc;  Elem(eb) % s(2) = s;    Elem(eb) % s(3) = sac

  !-------------------------------------------------!
  !   Change the orientation of the side properly   !
  !-------------------------------------------------!
  if(a < b) then
    side(s) % c  = a;   side(s) % d  = b
    side(s) % a  = d;   side(s) % b  = c
    side(s) % ea = ea;  side(s) % eb = eb
  else
    side(s) % c  = b;   side(s) % d  = a
    side(s) % a  = c;   side(s) % b  = d
    side(s) % ea = eb;  side(s) % eb = ea
  end if

  !-----------------------------------------------------!
  !   Update information on the neighbouring elements   !
  !-----------------------------------------------------!
  if(eac .ne. 0) then
    if(Elem(eac) % e(1) .eq. ea) Elem(eac) % e(1) = eb
    if(Elem(eac) % e(2) .eq. ea) Elem(eac) % e(2) = eb
    if(Elem(eac) % e(3) .eq. ea) Elem(eac) % e(3) = eb
  end if

  if(ebd .ne. 0) then
    if(Elem(ebd) % e(1) .eq. eb) Elem(ebd) % e(1) = ea
    if(Elem(ebd) % e(2) .eq. eb) Elem(ebd) % e(2) = ea
    if(Elem(ebd) % e(3) .eq. eb) Elem(ebd) % e(3) = ea
  end if

  !--------------------------------------------------!
  !   Update information on the neighbouring sides   !
  !--------------------------------------------------!
  if(side(sad) % ea .eq. ea) side(sad) % a = b
  if(side(sad) % eb .eq. ea) side(sad) % b = b

  if(side(sbc) % ea .eq. eb) side(sbc) % a = a
  if(side(sbc) % eb .eq. eb) side(sbc) % b = a

  if(side(sbd) % ea .eq. eb) then
    side(sbd) % ea = ea
    side(sbd) % a  = a
  end if
  if(side(sbd) % eb .eq. eb) then
    side(sbd) % eb = ea
    side(sbd) % b  = a
  end if

  if(side(sac) % ea .eq. ea) then
    side(sac) % ea = eb
    side(sac) % a  = b
  end if
  if(side(sac) % eb .eq. ea) then
    side(sac) % eb = eb
    side(sac) % b  = b
  end if

  end subroutine
