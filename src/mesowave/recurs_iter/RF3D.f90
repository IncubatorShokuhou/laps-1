subroutine rf3d(u,l,n,a,np)

!****************************************************
!  this routine applies the one-d recursive filter in
!  x, y, and z direction np times with alpha value of
!  a(1),a(2),a(3) to a three-d array u.
!
!  history: apr. 2003 by yuanfu xie.
!****************************************************

  implicit none

  integer, intent(in) :: l(3),n(3),np(3)
  real,    intent(in) :: a(3)
  real                :: u(l(1),l(2),l(3))

  ! local variables:
  integer :: j,k

  ! x direction:
  do k=1,n(3)
     do j=1,n(2)
        call rf1d(u(1:n(1),j,k),n(1),a(1),np(1))
     enddo
  enddo

  ! y direction:
  do k=1,n(3)
     do j=1,n(1)
        call rf1d(u(j,1:n(2),k),n(2),a(2),np(2))
     enddo
  enddo

  ! z direction:
  do k=1,n(2)
     do j=1,n(1)
        call rf1d(u(j,k,1:n(3)),n(3),a(3),np(3))
     enddo
  enddo

end subroutine rf3d
