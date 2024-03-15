subroutine rf1d(u,n,a,np)

!****************************************************
!  this routine applies the recursive filter to u for
!  np times with alpha value of a.
!
!  history: apr. 2003 by yuanfu xie.
!****************************************************

  implicit none

  integer, intent(in) :: n,np
  real,    intent(in) :: a
  real                :: u(n)

  ! local variables:
  integer :: i,ip
  real    :: one_a,r(n)

  one_a = 1.0-a

  ! recurisve filter number of np times:
  do ip=1,np
        
     ! left to right:
     if (ip .eq. 1) then
        r(1) = one_a*u(1) 
     else if (ip == 2 ) then
        r(1) = u(1) / ( 1.0 + a )
     else
        r(1) = one_a * ( u(1) - a**3 * u(2) ) / &
             (1.0-a**2)**2
     end if
     do i=2,n
        r(i) = a*r(i-1)+one_a*u(i)
     enddo
        
     ! right to left:
     if (ip .eq. 1) then
        u(n) = r(n)/(1.0+a) 
     else
        u(n) = one_a * ( r(n) - a**3 * r(n-1) ) / &
             (1.0-a**2)**2
     end if
     do i=n-1,1,-1
        u(i) = a*u(i+1)+one_a*r(i)
     enddo
  enddo

end subroutine rf1d
