program main

!*********************************************************
!  this program analyzes a set of surface in time using a
!  recursive filter.
!
!  history: jan. 2004 by yuanfu xie.
!*********************************************************

  use definition
  use initialize
  use minimizatn
  use configlaps

  implicit none

  ! local variable:
  integer :: id,ierr
  real    :: ds(3)

  call lapsinfo
  call lso_data

  ! call namelist

  ! call readobsn

  call grid2obs

  !call minimize
  ds(1) = grid_spacingx
  ds(2) = grid_spacingy
  ds(3) = d(3)
  do id=1,n(4)
     if (id .ne. 4) then  ! do not analyze station pressure
        call iterates(id,bkgd,ldf,nx,ny,ds,ncycles,nvlaps,nfic)
        print*,'variable ',id,' has been analyzed'
     endif
  enddo

  ! release memory of ldf:
  deallocate(bkgd,ldf,stat=ierr)

  ! call writeout
  ! call writeout
  call writeanalysis(s(1:n(1),1:n(2),1:n(3),1:n(4)),n)

end program main
