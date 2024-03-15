subroutine grid2obs

!*************************************************
!  this routine maps grid functions to observation
!  sites.
!
!  history: jan. 2004 by yuanfu xie.
!*************************************************

  implicit none

  integer :: iobs,ier

  do iobs=1,nobs

     call intplt3d(o(2,iobs),l,n,d,dm, &
                   idx(1,iobs),coe(1,iobs),ier)

     if (ier .ne. 0) then
        w(iobs) = 0.0
        idx(1:3,iobs) = 1
        coe(1:3,iobs) = 0.0
     endif

  enddo

end subroutine grid2obs
