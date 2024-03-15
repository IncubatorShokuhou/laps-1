subroutine gridbarnes(f,l,n,b)

!==========================================================
!  this routine applied a barnes analysis to a gridded two
!  dimensional function as a filter.
!
!  history: jan. 2005 by yuanfu xie.
!==========================================================

   implicit none

   integer, intent(in) :: l(2),n(2)
   real,    intent(in) :: f(l(1),l(2))
   real,    intent(out) :: b(l(1),l(2))

   ! local variables:
   integer :: i,j,ir,ni,nj,ip
   real    :: gama,s,w,ws,t(l(1),l(2)),tmp

   gama = 0.2
   s = 2.0e3/gama  ! kapa_0
   ir = 60	   ! neighbors
   
   ! clean analysis field:
   b = 0.0
   t = 0.0

   ! barnes analysis:
   do ip=1,1

      s = s*gama

      ! for every gridpoint:
      do j=1,n(2)
         do i=1,n(1)
	    
	    ! for each neighbor within ir:
	    ws = 0.0
	    tmp = 0.0
	    do nj=-ir,ir,2
	       do ni=-ir,ir,2
		  w = exp(-(float(ni)**2+float(nj)**2)/s)

		  if ((i+ni .ge. 1) .and. (i+ni .le. n(1)) .and. &
		      (j+nj .ge. 1) .and. (j+nj .le. n(2))) then
		     tmp = tmp+(f(i+ni,j+nj)-t(i+ni,j+nj))*w
		     ws = ws + w
		  endif

	       enddo
	    enddo

	    b(i,j) = b(i,j)+tmp/ws
         enddo
      enddo

      t = b

   enddo

end subroutine gridbarnes
