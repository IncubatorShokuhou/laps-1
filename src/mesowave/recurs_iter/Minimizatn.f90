module minimizatn

!*********************************************************
!  this module defines a minimization procedure analyzing
!  surface observations in time sequence.
!
!  history: jan. 2004 by yuanfu xie.
!*********************************************************

  use definition

contains
  
  include 'functn.f90'
  include 'functn_ad.f90'
  ! include 'functndiv.f90'
  ! include 'functndiv_ad.f90'
  include 'iterates.f90'
  include 'minimize.f90'
  include 'rf1d.f90'
  include 'rf1d_ad.f90'
  include 'rf3d.f90'
  include 'rf3d_ad.f90'

end module minimizatn
