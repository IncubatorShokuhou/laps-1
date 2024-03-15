module drawcountour
!*************************************************
! draw countour picture for sufer software
! history: february 2008, separated from the input_bg_obs module by zhongjie he.
!*************************************************

   use prmtrs_stmas

   public drcontour, drcontour_2d

contains

   subroutine drcontour(ox, ex, oy, ey, nd, ng, uv, fn)
!*************************************************
! draw radial wind (affiliate)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: nd, ng(numdims), iu, i, j, k
      real  :: mn, mx, ox, ex, oy, ey
      real  :: uv(ng(1), ng(2), ng(3))
      character(len=6) :: fn
! --------------------
      k = 2 !(ng(3)+1)/2
      iu = 2
      open (iu, file=fn, status='unknown')
      mx = -100000000.0
      mn = 100000000.0
      do i = 1, ng(1)
      do j = 1, ng(2)
         if (uv(i, j, k) .lt. 9000.0) then
            if (uv(i, j, k) .gt. mx) mx = uv(i, j, k)
            if (uv(i, j, k) .lt. mn) mn = uv(i, j, k)
         end if
      end do
      end do
      write (iu, '(a4)') 'dsaa'
      write (iu, '(2i4)') ng(1), ng(2)
      write (iu, *) ox, ex
      write (iu, *) oy, ey
      write (iu, *) mn, mx
      do i = 1, ng(1)
      do j = 1, ng(2)
         if (uv(i, j, k) .gt. 9000.0) uv(i, j, k) = 2.e38
      end do
      end do
      do j = 1, ng(2)
         write (iu, *) (uv(i, j, k), i=1, ng(1))
      end do
      close (iu)
      return
   end subroutine drcontour

   subroutine drcontour_2d(ox, ex, oy, ey, im, jm, uv, fn)
!*************************************************
! draw radial wind (affiliate)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: iu, i, j, k, im, jm
      real  :: mn, mx, ox, ex, oy, ey
      real  :: uv(im, jm)
      character(len=6) :: fn
! --------------------
      iu = 2
      open (iu, file=fn, status='unknown')
      mx = -100000000.0
      mn = 100000000.0
      do i = 1, im
      do j = 1, jm
         if (uv(i, j) .lt. 9000.0) then
            if (uv(i, j) .gt. mx) mx = uv(i, j)
            if (uv(i, j) .lt. mn) mn = uv(i, j)
         end if
      end do
      end do
      write (iu, '(a4)') 'dsaa'
      write (iu, '(2i4)') im, jm
      write (iu, *) ox, ex
      write (iu, *) oy, ey
      write (iu, *) mn, mx
      do i = 1, im
      do j = 1, jm
         if (uv(i, j) .gt. 9000.0) uv(i, j) = 2.e38
      end do
      end do
      do j = 1, jm
         write (iu, *) (uv(i, j), i=1, im)
      end do
      close (iu)
      return
   end subroutine drcontour_2d

end module drawcountour

