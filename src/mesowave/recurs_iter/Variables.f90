!****************************************************
!  this part of definition module defines all global
!  variables in this data assimilation.
!
!  history: jan. 2004 by yuanfu xie.
!****************************************************

!===================
!  file names:
!===================

character*60 :: datafile
integer      :: namelens

!===================
!  parameters:
!===================

integer      :: maxitr

!===================
!  arrays:
!===================

integer      :: l(4),n(4),np(1:3,1:mv),nrf(mv)
real         :: a(mx,my,mt,mv),dm(2,3),d(3),al(1:3,1:mv)
real         :: s(mx,my,mt,mv),qc_cons(mv)
