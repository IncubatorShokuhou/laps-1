subroutine writeout

!*************************************************
!  this routine writes out the solution to a file.
!
!  history: feb. 2004 by yuanfu xie.
!*************************************************

  use definition

  implicit none

  open(unit=10,file='analysis.dat')

  ! write(10,*) n

  ! write(10,*) s(1:n(1),1:n(2),1:n(3),1:n(4))

  write(10,*) n(1),n(2),1,n(4)

  write(10,*) s(1:n(1),1:n(2),n(3),1:n(4))

  close(10)

  !open(unit=10,file='../dat/analysis.bin',form='unformatted')

  !write(10) n

  !write(10) s(1:n(1),1:n(2),1:n(3),1:n(4))

  !close(10)

end subroutine writeout
