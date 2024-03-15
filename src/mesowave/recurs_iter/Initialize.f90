module initialize

!****************************************************
!  this module initializes the gps data assimilation
!  package.
!
!  history: jan. 2004 by yuanfu xie.
!****************************************************

  use definition
  use util_tools

contains

  include 'namelist.f90'
  include 'readobsn.f90'
  include 'grid2obs.f90'

end module initialize
