subroutine intplt3d(x,l,n,d,dm,idx,coe,ier)

!*************************************************
!  this routine computes an interpolation indices
!  and coefficients using bi-linear interpolation.
!
!  history: jan. 2004 by yuanfu xie.
!*************************************************

  implicit none

  integer, intent(in) :: l(3),n(3)
  real,    intent(in) :: x(3),d(3),dm(2,3)

  integer, intent(out) :: idx(3),ier
  real,    intent(out) :: coe(3)

  ! local variables:
  integer :: i

  ier = 0
  ! check in box?
  do i=1,3
     if (x(i) .lt. dm(1,i)) ier = -i
     if (x(i) .gt. dm(2,i)) ier = i
  enddo

  ! indices:
  idx = (x-dm(1,1:3))/d

  ! coefficients:
  coe = (x-idx*d-dm(1,1:3))/d

  idx = idx+1

end subroutine intplt3d
