!***********************************************************
!  this header file defines large arrays for observations to
!  avoid the limitation of parameter passing.
!
!  history: jan. 2004 by yuanfu xie.
!***********************************************************

! maximum number of obs:
integer, parameter :: mobs = 200000

! observations:
integer :: nobs
integer :: vid(mobs),idx(3,mobs)
real    :: o(4,mobs),coe(3,mobs),w(mobs)
common /obsblock/nobs,vid,idx,o,coe,w
