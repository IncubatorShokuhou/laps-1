subroutine readobsn

!*********************************************************
!  this routine reads in all observations avaiable on this
!  surface analysis.
!
!  history: jan. 2004 by yuanfu xie.
!*********************************************************

  implicit none

  integer :: m

  open(unit=10,file=datafile(1:namelens),status='old')

  ! domain x time interval:
  read(10,*) dm(1:2,1:3)
  d(1:3) = (dm(2,1:3)-dm(1,1:3))/float(n(1:3)-1)

  nobs = 1		! count observations
  m = 0			! count variables to be analyzed

1 format(i2,5e14.6)
2 continue
  read(10,1,end=9) vid(nobs),o(1:4,nobs),w(nobs)
  if (m .lt. vid(nobs)) m = vid(nobs)
  nobs = nobs+1
  goto 2

9 continue

  ! close open file:
  close(10)

  nobs = nobs-1
  
  ! check dimension for variables:
  if (m .gt. mv) then
     print*,'readobsn: too much variables'
     stop
  endif
  if (m .le. 0) then
     print*,'readobsn: no variables to analyze'
     stop
  endif
  ! check if number of obs is fit into the array:
  if (nobs .gt. mobs) then
     print*,'readobsn: too many obs'
     stop
  endif
  if (nobs .le. 0) then
     print*,'readobsn: no obs to analyze'
     stop
  endif

end subroutine readobsn
