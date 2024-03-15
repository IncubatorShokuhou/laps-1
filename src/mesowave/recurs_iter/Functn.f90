! 1 "functn.f"
! 1 "<built-in>"
! 1 "<command line>"
! 1 "functn.f"
subroutine functn(f,v,l,n,id,np,al)

!***************************************************************
! this routine evaluates the cost function of a surface data
! analysis problem.
!
! history: jan. 2004 by yuanfu xie.
!***************************************************************

  implicit none

  double precision, intent(out) :: f

  integer, intent(in) :: l(4),n(4) ! grid dimensions
  integer, intent(in) :: id ! variable id
  integer, intent(in) :: np(3,l(4)) ! rf passes

  real, intent(in) :: al(3,l(4)) ! rf alpha values
  real, intent(in) :: v(l(1),l(2),l(3)) ! control grid

  ! local variables:
  integer :: iobs,i,j,k
  real :: x(l(1),l(2),l(3)),vo,a(2,3)

! 1 "../../obscommon.f90" 1
!***********************************************************
! this header file defines large arrays for observations to
! avoid the limitation of parameter passing.
!
! history: jan. 2004 by yuanfu xie.
!***********************************************************

! maximum number of obs:
integer, parameter :: mobs = 200000

! observations:
integer :: nobs
integer :: vid(mobs),idx(3,mobs)
real :: o(4,mobs),coe(3,mobs),w(mobs)
common /obsblock/nobs,vid,idx,o,coe,w
! 26 "functn.f" 2

  !------------------------------------------------------------
  ! recursive filter:
  !------------------------------------------------------------
  x = v
  call rf3d(x(1,1,1),l,n,al(1,id),np(1,id))

  !------------------------------------------------------------
  ! evaluate the cost function:
  !------------------------------------------------------------

  f = 0.0d0

  do iobs=1,nobs

     if (id .eq. vid(iobs)) then

        ! h(x) = y
        vo = 0.0
        a(1,1:3) = 1.0-coe(1:3,iobs)
        a(2,1:3) = coe(1:3,iobs)
        do k=1,2
           if ((idx(3,iobs)+k-1 .ge. 1) .and. &
               (idx(3,iobs)+k-1 .ge. n(3))) then
           do j=1,2
              do i=1,2
                 vo = vo + x(idx(1,iobs)+i-1,idx(2,iobs)+j-1, &
                             idx(3,iobs)+k-1)* &
                             a(i,1)*a(j,2)*a(k,3)
              enddo
           enddo
           endif
        enddo

        ! y - y0
        f = f + w(iobs)*(vo-o(1,iobs))**2

     endif

  enddo

  f = 0.5*f

end subroutine functn
