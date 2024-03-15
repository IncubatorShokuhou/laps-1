subroutine namelist

!*********************************************************
!  this routine reads in namelist for the surface analysis
!
!  history: jan. 2004 by yuanfu xie.
!*********************************************************

  use definition

  integer :: i

  open(unit=10,file='run/namelist.txt',status='old')
  read(10,*) ! skip headline

  l(1:4) = (/mx,my,mt,mv/)

  ! grid dimensions:
  read(10,*) n(1),n(2),n(3),n(4)

  ! check:
  if ((n(1) .gt. mx) .or. (n(2) .gt. my) .or. &
      (n(3) .gt. mt) .or. (n(4) .gt. mv)) then
     print*,'namelist: analysis array is too small!'
     stop
  endif

  ! data filename:
  read(10,*) namelens,datafile

  ! recursive filters:
  read(10,*) ! skip a specification line
  do i=1,n(4)
     read(10,*) al(1:3,i),np(1:3,i)
  enddo

  ! number of minimization iterations:
  read(10,*) maxitr

  ! number of recursive filter iterations:
  read(10,*) nrf(1:n(4))

  close(10)

end subroutine namelist
