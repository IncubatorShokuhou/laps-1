module generaltools

contains

   subroutine interpltn(nd, ng, co, ac, oc)
!*************************************************
! calculate the interpolating coefficent (general)
! history: august 2007, coded by wei li.
!*************************************************
! 2     , 3     , 5     , 9
! 2**0+1, 2**1+1, 2**2+1, 2**3+1
      implicit none
! --------------------
      integer  :: i, j
      integer, intent(in)  :: nd, ng
      real, intent(in)  :: ac(nd, ng), oc(nd)
      real, intent(out) :: co(ng)
      real                  :: dd(nd, ng), dc(nd)
! --------------------
! declare:
!        'nd' is the number of dimensions of the field
!        'ng' is the number of gridpoints neighbor to the position, which is to be interpolated, ng is equal to 2**nd
!        'co' is the interpolating coefficients.
!        'oc' is the position where to be interpolated
!        'ac' is the index of the gridpoints, which is neighbor to the positioin
! --------------------
      do i = 1, nd
         if ((oc(i) - ac(i, 1))*(oc(i) - ac(i, 2**(i - 1) + 1)) .gt. 0.0) then
            print *, 'interpltn: oc is not in the domain ', i, oc(i), ac(i, 1), ac(i, 2**(i - 1) + 1)
            stop
         end if
      end do
      do i = 1, nd
         do j = 1, ng
            dd(i, j) = ac(i, j) - oc(i)
         end do
         dc(i) = abs(ac(i, 2**(i - 1) + 1) - ac(i, 1))
!    if(dc(i).le.0.0)stop 'grid arrange wrong!!!'            ! omitted by zhongjie he
      end do
      do i = 1, nd
         if (dc(i) .ne. 0) then
            do j = 1, ng
               dd(i, j) = dd(i, j)/dc(i)
               if (dd(i, j) .gt. 0.0) then
                  dd(i, j) = dd(i, j) - 1.0
               else
                  dd(i, j) = dd(i, j) + 1.0
               end if
               dd(i, j) = abs(dd(i, j))
            end do
         else                                     !     added by zhongjie he
            do j = 1, ng
               dd(i, j) = 1.
            end do
         end if
      end do
      do j = 1, ng
         co(j) = 1.0
         do i = 1, nd
            co(j) = co(j)*dd(i, j)
         end do
      end do
      return
   end subroutine interpltn

   subroutine interpltn_xie(nd, ng, co, ac, oc, ip, np, pv)
!*************************************************
! calculate the interpolating coefficent (general)
! history: august 2007, coded by wei li.
!
!          march 31, 2009 by yuanfu xie for log(p)
!          interpolation.
!*************************************************
! 2     , 3     , 5     , 9
! 2**0+1, 2**1+1, 2**2+1, 2**3+1
      implicit none
! --------------------
      integer  :: i, j, k, k1
      integer, intent(in)  :: nd, ng, ip, np
      real, intent(in)  :: ac(nd, ng), oc(nd), pv(np)
      real, intent(out) :: co(ng)
      real                  :: dd(nd, ng), dc(nd), p
! --------------------
! declare:
!        'nd' is the number of dimensions of the field
!        'ng' is the number of gridpoints neighbor to the position, which is to be interpolated, ng is equal to 2**nd
!        'co' is the interpolating coefficients.
!        'oc' is the position where to be interpolated
!        'ac' is the index of the gridpoints, which is neighbor to the positioin
! --------------------
      do i = 1, nd
         if ((oc(i) - ac(i, 1))*(oc(i) - ac(i, 2**(i - 1) + 1)) .gt. 0.0) then
            print *, 'interpltn: oc is not in the domain ', i, oc(i), ac(i, 1), ac(i, 2**(i - 1) + 1)
            stop
         end if
      end do
      do i = 1, nd
         do j = 1, ng
            dd(i, j) = ac(i, j) - oc(i)
         end do
         dc(i) = abs(ac(i, 2**(i - 1) + 1) - ac(i, 1))
      end do
      do i = 1, nd
         if (dc(i) .ne. 0) then
            do j = 1, ng
               dd(i, j) = dd(i, j)/dc(i)
               if (dd(i, j) .gt. 0.0) then
                  ! log(p):
                  if (i .eq. ip) then
                     k1 = ac(i, j)
                     if (k1 .le. 1) then
                        dd(i, j) = 0.0
                     else
                        k = k1 - 1
                        p = pv(k)*dd(i, j) + (1.0 - dd(i, j))*pv(k + 1)
                        dd(i, j) = (log(p) - log(pv(k)))/(log(pv(k1)) - log(pv(k)))
                     end if
                  else
                     dd(i, j) = dd(i, j) - 1.0
                  end if
               else
                  if (i .eq. ip) then
                     k = ac(i, j)
                     if (k .ge. np) then
                        dd(i, j) = 1.0
                     else
                        k1 = k + 1
                        p = (1.0 + dd(i, j))*pv(k) - dd(i, j)*pv(k1)
                        dd(i, j) = (log(pv(k1)) - log(p))/(log(pv(k1)) - log(pv(k)))
                     end if
                  else
                     dd(i, j) = dd(i, j) + 1.0
                  end if
               end if
               dd(i, j) = abs(dd(i, j))
            end do
         else                                     !     added by zhongjie he
            do j = 1, ng
               dd(i, j) = 1.
            end do
         end if
      end do

      do j = 1, ng
         co(j) = 1.0
         do i = 1, nd
            co(j) = co(j)*dd(i, j)
         end do
      end do
      return
   end subroutine interpltn_xie

   subroutine gderiveit(z1, z2, z3, a, b, c)
!*************************************************
! general derivative of interior (general)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      real  :: x, y
      real, intent(in)  :: z1, z2, z3
      real, intent(out) :: a, b, c
! --------------------
      x = z2 - z1
      y = z3 - z2
      a = -y/(x*x + x*y)
      b = (-x + y)/(x*y)
      c = x/(x*y + y*y)
      return
   end subroutine gderiveit

   subroutine gderivelb(z1, z2, z3, a, b, c)
!*************************************************
! general derivative of left boundary (general)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      real  :: x, y
      real, intent(in)  :: z1, z2, z3
      real, intent(out) :: a, b, c
! --------------------
      x = z2 - z1
      y = z3 - z2
      a = (-2.0*x - y)/(x*x + x*y)
      b = (x + y)/(x*y)
      c = -x/(x*y + y*y)
      return
   end subroutine gderivelb

   subroutine gderiverb(z1, z2, z3, a, b, c)
!*************************************************
! general derivative of right boundary (general)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      real  :: x, y
      real, intent(in)  :: z1, z2, z3
      real, intent(out) :: a, b, c
! --------------------
      x = z2 - z1
      y = z3 - z2
      a = y/(x*x + x*y)
      b = (-x - y)/(x*y)
      c = (x + 2.0*y)/(x*y + y*y)
      return
   end subroutine gderiverb

   subroutine g2orderit(z1, z2, z3, a, b, c)
!*************************************************
! general 2-order derivative of interior (general)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      real  :: x, y
      real, intent(in)  :: z1, z2, z3
      real, intent(out) :: a, b, c
! --------------------
      x = z2 - z1
      y = z3 - z2
      a = 2.0/(x*x + x*y)
      b = -2.0/(x*y)
      c = 2.0/(x*y + y*y)
      return
   end subroutine g2orderit

   subroutine pstn2numb(nn, np, nm, nc)
!*************************************************
! transform position to number (general)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer               :: pt, m, n
      integer, intent(in)  :: nn
      integer, intent(in)  :: np(nn), nm(nn)
      integer, intent(out) :: nc
! --------------------
! declare:
!        'nn' is the number of dimensions of the field
!        'np' is the index of the gridpoint at each dimension
!        'nm' is the total grid number of each dimensioin
!        'nc' is the index of the gridpoint, when the field is translated to a 1 dimensional filed
! --------------------
      nc = np(1)
      if (nn .eq. 1) return
      do n = 2, nn
         pt = np(n) - 1
         do m = 1, n - 1
            pt = pt*nm(m)
         end do
         nc = nc + pt
      end do
      return
   end subroutine pstn2numb

   subroutine numb2pstn(nn, np, nm, nc)
!*************************************************
! transform number to position (general) (not used)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer               :: n, n0
      integer, intent(in)  :: nn, nc
      integer, intent(in)  :: nm(nn)
      integer, intent(out) :: np(nn)
! --------------------
! declare:
!        'nn' is the number of dimensions of the field
!        'np' is the index of the gridpoint at each dimension, when the field is translated to a nn dimensional filed
!        'nm' is the total grid number of each dimensioin
!        'nc' is the index of the gridpoint in the 1 dimensional field
! --------------------
      n0 = nc
      np(1) = mod(n0, nm(1))
      if (np(1) .eq. 0) np(1) = nm(1)
      if (nn .eq. 1) return
      do n = 2, nn
         n0 = (n0 - np(n - 1))/nm(n - 1) + 1
         np(n) = mod(n0, nm(n))
         if (nm(n) .eq. 1) np(n) = 1
         if (np(n) .eq. 0) np(n) = nm(n)
      end do
      return
   end subroutine numb2pstn

   subroutine vrtclpstn8(km, lm, pp, p, ps, is)
!*************************************************
! get vertical position (general)
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: k
      integer, intent(in)  :: km
      real, intent(in)  :: pp(km)
      real, intent(in)  :: p, lm
      real, intent(out) :: ps
      integer, intent(out) :: is
! --------------------
! declare:
!        'km' is the number of levels of the field
!        'pp' is the presure of each level
!        'p'  is the presure of the position, which we want to get its index in the vertical coordinate
!        'ps' is the index of the point in the vertical coordinate
!        'is' is the state of the result of this subroutine, 0 means wrong, 1 means ok
!        'lm' is a limitation to decide whether the observation is available (just used for the case: km=1)
! --------------------
      is = 0
      if (km .ge. 2) then
         do k = 1, km - 1
            if ((p - pp(k))*(p - pp(k + 1)) .gt. 0.0) cycle
            ps = k + abs(p - pp(k))/abs(pp(k + 1) - pp(k))
            is = 1
         end do
      elseif (km .eq. 1) then
         if (abs(p - pp(1)) .le. lm) then
            ps = 1
            is = 1
         end if
      end if
      return
   end subroutine vrtclpstn8

   subroutine vrtclpstn(km, lm, pp, p, ps, is)
!*************************************************
! get vertical position (general)
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: k
      integer, intent(in)  :: km
      real, intent(in)  :: pp(km), p, lm
      real, intent(out) :: ps
      integer, intent(out) :: is
! --------------------
! declare:
!        'km' is the number of levels of the field
!        'pp' is the presure of each level
!        'p'  is the presure of the position, which we want to get its index in the vertical coordinate
!        'ps' is the index of the point in the vertical coordinate
!        'is' is the state of the result of this subroutine, 0 means wrong, 1 means ok
!        'lm' is a limitation of preesure to decide whether the observation is available (just used for the case: km=1)
! --------------------
      is = 0
      if (km .ge. 2) then
         do k = 1, km - 1
            if ((p - pp(k))*(p - pp(k + 1)) .gt. 0.0) cycle
            ps = k + abs(p - pp(k))/abs(pp(k + 1) - pp(k))
            is = 1
         end do
      elseif (km .eq. 1) then
         if (abs(p - pp(1)) .le. lm) then
            ps = 1
            is = 1
         end if
      end if

      return
   end subroutine vrtclpstn

   subroutine zpconvert(km, zz, pp, z, p, is)
!*************************************************
! find pressure from height (general)
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer              :: k
      integer, intent(in)  :: km
      real, intent(in)  :: zz(km), z
      real, intent(in)  :: pp(km)
      real, intent(out) :: p
      integer, intent(out) :: is
! --------------------
! declare:
!        'km' is the number of levels of the field
!        'zz' is the height of each level
!        'pp' is the presure of each level
!        'z'  is the height of point, which we want to get its presure
!        'p'  is the presure calculated from the zz, pp, and z
!        'is' is the state of the result of this subroutine, 0 means wrong, 1 means ok
! --------------------
      is = 0
      do k = 1, km - 1
         if ((z - zz(k))*(z - zz(k + 1)) .gt. 0.0) cycle
         p = pp(k)*abs(zz(k + 1) - z)/abs(zz(k + 1) - zz(k)) + pp(k + 1)*abs(z - zz(k))/abs(zz(k + 1) - zz(k))
         is = 1
      end do
      return
   end subroutine zpconvert

   subroutine fcst2bkgd(nd, ng, ns, nf, xf, yf, zf, tf, nb, xb, yb, zb, tb, fc, bk)
!*************************************************
! interpolate model forecast to background used in analyzing (general)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer, parameter  :: md = 4
      integer, intent(in) :: nd, ng, ns
      integer, intent(in) :: nf(md), nb(md)
      real, intent(in) :: xf(nf(1)), yf(nf(2)), zf(nf(3)), tf(nf(4))
      real, intent(in) :: xb(nb(1)), yb(nb(2)), zb(nb(3)), tb(nb(4))
      real, intent(in) :: fc(nf(1), nf(2), nf(3), nf(4), ns)
      integer              :: i, j, k, l, i0, j0, k0, l0, i1, j1, k1, l1, s, m, n, nn
      real                 :: pp(md)
      real                 :: ac(nd, ng), oc(nd), co(ng)
      real, intent(out):: bk(nb(1), nb(2), nb(3), nb(4), ns)
      integer              :: lm(nd)
! --------------------
! declare:
!        'nd' is the number of dimensions of the field
!        'ng' is the number of gridpoints neighbor to the position, which is to be interpolated, ng is equal to 2**nd
!        'ns' is the number of variables that to be interpolated
!        'nf' the grid number of the original field
!        'nb' the grid number of the target field
!        'xf', 'yf', 'zf' and 'tf' are the scales of the original filed in the x, y, z and t coordinates respectively
!        'xb', 'yb', 'zb' and 'tb' are the scales of the target field n the x, y, z and t coordinates
!        'fc' the original field
!        'bk' the target field
!        'co' is the interpolating coefficients.
!        'oc' is the position where to be interpolated
!        'ac' is the index of the gridpoints, which is neighbor to the positioin
! --------------------
      do l = 1, nb(4)
      do k = 1, nb(3)
      do j = 1, nb(2)
      do i = 1, nb(1)
         call checkpstn(nf(1), xf, xb(i), i0)
         call checkpstn(nf(2), yf, yb(j), j0)
         call checkpstn(nf(3), zf, zb(k), k0)
         call checkpstn(nf(4), tf, tb(l), l0)
         if (i0 .eq. 0 .or. j0 .eq. 0 .or. k0 .eq. 0 .or. l0 .eq. 0) stop 'backgound can not be generated'
         m = 0
!===============================================
!    do l1=l0,min0(l0+1,nf(4))
!    do k1=k0,min0(k0+1,nf(3))
!    do j1=j0,min0(j0+1,nf(2))
!    do i1=i0,min0(i0+1,nf(1))
!      m=m+1
!      pp(1)=xf(i1)
!      pp(2)=yf(j1)
!      pp(3)=zf(k1)
!      pp(4)=tf(l1)
!      nn=0
!      do n=1,md
!        if(nf(n).ge.2)then
!          nn=nn+1
!          ac(nn,m)=pp(n)
!        endif
!      enddo
!    enddo
!    enddo
!    enddo
!    enddo
!================ modified by zhongjie he =====
         lm(1) = i0 + 1
         lm(2) = j0 + 1
         lm(3) = k0 + 1
         lm(4) = l0 + 1
         do n = 1, nd
            if (nd .lt. n) lm(n) = lm(n) - 1
         end do
         do l1 = l0, lm(4)
         do k1 = k0, lm(3)
         do j1 = j0, lm(2)
         do i1 = i0, lm(1)
            m = m + 1
            pp(1) = xf(min0(i1, nf(1)))
            pp(2) = yf(min0(j1, nf(3)))
            pp(3) = zf(min0(k1, nf(2)))
            pp(4) = tf(min0(l1, nf(1)))
            nn = 0
            do n = 1, md
               if (nf(n) .ge. 2) then
                  nn = nn + 1
                  ac(nn, m) = pp(n)
               end if
            end do
         end do
         end do
         end do
         end do
!===============================================
         pp(1) = xb(i)
         pp(2) = yb(j)
         pp(3) = zb(k)
         pp(4) = tb(l)
         nn = 0
         do n = 1, md
            if (nf(n) .ge. 2) then
               nn = nn + 1
               oc(nn) = pp(n)
            end if
         end do
         call interpltn(nd, ng, co, ac, oc)
         do s = 1, ns
            m = 0
            bk(i, j, k, l, s) = 0.0
            do l1 = l0, min0(l0 + 1, nf(4))
            do k1 = k0, min0(k0 + 1, nf(3))
            do j1 = j0, min0(j0 + 1, nf(2))
            do i1 = i0, min0(i0 + 1, nf(1))
               m = m + 1
               bk(i, j, k, l, s) = bk(i, j, k, l, s) + co(m)*fc(i1, j1, k1, l1, s)
            end do
            end do
            end do
            end do
         end do
      end do
      end do
      end do
      end do
      return
   end subroutine fcst2bkgd

   subroutine checkpstn(km, zz, z, kz)
!*************************************************
! find the position (general)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: k
      integer, intent(in)  :: km
      real, intent(in)  :: zz(km), z
      integer, intent(out) :: kz
! --------------------
      kz = 0
      if (km .eq. 1 .and. zz(1) .eq. z) then
         kz = 1
         return
      end if
      do k = 1, km - 1
         if ((z - zz(k))*(z - zz(k + 1)) .gt. 0.0) cycle
         kz = k
         exit
      end do
      return
   end subroutine checkpstn

   subroutine gridtrans(nv, n1, x1, y1, z1, t1, n2, x2, y2, z2, t2, g1, g2)
!*************************************************
!  this routine maps grid functions to other grid.
!
!  history: apr. 2008, coded by yuanfu xie.
!*************************************************

      implicit none

      integer, intent(in) :: nv        ! number of var to map
      integer, intent(in) :: n1(4)        ! numbers of input gridpts
      integer, intent(in) :: n2(4)        ! numbers of output gridpts

      real, intent(in) :: x1(n1(1)), y1(n1(2)), z1(n1(3)), t1(n1(4))
      ! bkg grid positions
      real, intent(in) :: x2(n2(1)), y2(n2(2)), z2(n2(3)), t2(n2(4))
      ! fine grid positions
      real, intent(in) :: g1(n1(1), n1(2), n1(3), n1(4), nv)
      ! background variables
      real, intent(out) :: g2(n2(1), n2(2), n2(3), n2(4), nv)

      ! local variables:
      integer :: i, j, k, l
      integer :: idx(2, n2(1)), idy(2, n2(2)), idz(2, n2(3)), idt(2, n2(4))
      real    :: cox(2, n2(1)), coy(2, n2(2)), coz(2, n2(3)), cot(2, n2(4)), s

      ! coefficients and indices:
      ! left:
      idx(1:2, 1) = 1
      idy(1:2, 1) = 1
      idz(1:2, 1) = 1
      idt(1:2, 1) = 1
      cox(1, 1) = 1.0
      coy(1, 1) = 1.0
      coz(1, 1) = 1.0
      cot(1, 1) = 1.0
      cox(2, 1) = 0.0
      coy(2, 1) = 0.0
      coz(2, 1) = 0.0
      cot(2, 1) = 0.0
      ! right:
      idx(1:2, n2(1)) = n1(1)
      idy(1:2, n2(2)) = n1(2)
      idz(1:2, n2(3)) = n1(3)
      idt(1:2, n2(4)) = n1(4)
      cox(1, n2(1)) = 1.0
      coy(1, n2(2)) = 1.0
      coz(1, n2(3)) = 1.0
      cot(1, n2(4)) = 1.0
      cox(2, n2(1)) = 0.0
      coy(2, n2(2)) = 0.0
      coz(2, n2(3)) = 0.0
      cot(2, n2(4)) = 0.0

      ! x: uniform grid
      s = x1(2) - x1(1)
      do i = 2, n2(1) - 1
         idx(1, i) = int((x2(i) - x1(1))/s) + 1
         idx(2, i) = idx(1, i) + 1
         cox(2, i) = (x2(i) - x1(idx(1, i)))/s
         cox(1, i) = 1.0 - cox(2, i)
      end do

      ! y: uniform grid
      s = y1(2) - y1(1)
      do j = 2, n2(2) - 1
         idy(1, j) = int((y2(j) - y1(1))/s) + 1
         idy(2, j) = idy(1, j) + 1
         coy(2, j) = (y2(j) - y1(idy(1, j)))/s
         coy(1, j) = 1.0 - coy(2, j)
      end do

      ! z: non-uniform grid
      do k = 2, n2(3) - 1
         idz(1, k) = idz(1, k - 1)
         idz(2, k) = idz(1, k) + 1
         do while ((idz(1, k) .lt. n1(3)) .and. ( &
                   ((z1(1) .lt. z1(2)) .and. (z2(k) .ge. z1(idz(2, k)))) .or. &
                   ((z1(1) .gt. z1(2)) .and. (z2(k) .le. z1(idz(2, k))))))
            idz(1:2, k) = idz(1:2, k) + 1
         end do
      end do

      ! t: uniform grid:
      if (n1(4) .gt. 1) s = t1(2) - t1(1)
      do l = 2, n2(4) - 1
         idt(1, l) = int((t2(l) - t1(1))/s) + 1
         idt(2, l) = idt(1, l) + 1
         cot(2, l) = (t2(l) - t1(idt(1, l)))/s
         cot(1, l) = 1.0 - cot(2, l)
      end do

      ! interpolations:
      do l = 1, n2(4)
      do k = 1, n2(3)
      do j = 1, n2(2)
      do i = 1, n2(1)
         g2(i, j, k, l, 1:nv) = &
            g1(idx(1, i), idy(1, j), idz(1, k), idt(1, l), 1:nv)*cox(1, i)*coy(1, j)*coz(1, k)*cot(1, l) + &
            g1(idx(2, i), idy(1, j), idz(1, k), idt(1, l), 1:nv)*cox(2, i)*coy(1, j)*coz(1, k)*cot(1, l) + &
            g1(idx(1, i), idy(2, j), idz(1, k), idt(1, l), 1:nv)*cox(1, i)*coy(2, j)*coz(1, k)*cot(1, l) + &
            g1(idx(1, i), idy(1, j), idz(2, k), idt(1, l), 1:nv)*cox(1, i)*coy(1, j)*coz(2, k)*cot(1, l) + &
            g1(idx(1, i), idy(1, j), idz(1, k), idt(2, l), 1:nv)*cox(1, i)*coy(1, j)*coz(1, k)*cot(2, l) + &
            g1(idx(2, i), idy(2, j), idz(1, k), idt(1, l), 1:nv)*cox(2, i)*coy(2, j)*coz(1, k)*cot(1, l) + &
            g1(idx(2, i), idy(1, j), idz(2, k), idt(1, l), 1:nv)*cox(2, i)*coy(1, j)*coz(2, k)*cot(1, l) + &
            g1(idx(2, i), idy(1, j), idz(1, k), idt(2, l), 1:nv)*cox(2, i)*coy(1, j)*coz(1, k)*cot(2, l) + &
            g1(idx(1, i), idy(2, j), idz(2, k), idt(1, l), 1:nv)*cox(1, i)*coy(2, j)*coz(2, k)*cot(1, l) + &
            g1(idx(1, i), idy(2, j), idz(1, k), idt(2, l), 1:nv)*cox(1, i)*coy(2, j)*coz(1, k)*cot(2, l) + &
            g1(idx(1, i), idy(1, j), idz(2, k), idt(2, l), 1:nv)*cox(1, i)*coy(1, j)*coz(2, k)*cot(2, l) + &
            g1(idx(2, i), idy(2, j), idz(2, k), idt(1, l), 1:nv)*cox(2, i)*coy(2, j)*coz(2, k)*cot(1, l) + &
            g1(idx(2, i), idy(2, j), idz(1, k), idt(2, l), 1:nv)*cox(2, i)*coy(2, j)*coz(1, k)*cot(2, l) + &
            g1(idx(2, i), idy(1, j), idz(2, k), idt(2, l), 1:nv)*cox(2, i)*coy(1, j)*coz(2, k)*cot(2, l) + &
            g1(idx(1, i), idy(2, j), idz(2, k), idt(2, l), 1:nv)*cox(1, i)*coy(2, j)*coz(2, k)*cot(2, l) + &
            g1(idx(2, i), idy(2, j), idz(2, k), idt(2, l), 1:nv)*cox(2, i)*coy(2, j)*coz(2, k)*cot(2, l)
      end do
      end do
      end do
      end do

   end subroutine gridtrans

   subroutine bkgtofine(ns, nb, xb, yb, zb, tb, nf, xf, yf, zf, tf, bk, fg)
!*************************************************
!  this routine maps model background variables to
!  the finest grid of stmas multigrid.
!
!  history: apr. 2008, coded by yuanfu xie.
!*************************************************

      implicit none

      integer, intent(in) :: ns        ! number of var to map
      integer, intent(in) :: nb(4)        ! numbers of bk gridpoint
      integer, intent(in) :: nf(4)        ! numbers of fine gridpts

      real, intent(in) :: xb(nb(1)), yb(nb(2)), zb(nb(3)), tb(nb(4))
      ! bkg grid positions
      real, intent(in) :: xf(nf(1)), yf(nf(2)), zf(nf(3)), tf(nf(4))
      ! fine grid positions
      real, intent(in) :: bk(nb(1), nb(2), nb(3), nb(4), ns)
      ! background variables
      real, intent(out) :: fg(nf(1), nf(2), nf(3), nf(4), ns)

      ! local variables:
      integer :: i, j, k, l, id(2, 4)
      real    :: co(2, 4)

      ! for all grid points over the finest grid:
      id(1, 4) = 1
      id(2, 4) = min0(2, nb(4))
      do l = 1, nf(4)

         ! fine grid points passes bk grid points:
         if (nb(4) .eq. 1) then        ! only one gridpoint case
            id(2, 4) = id(1, 4) + 1
         else
            do while ((id(1, 4) .lt. nb(4)) .and. (tf(l) .ge. tb(id(2, 4))))
               id(1, 4) = id(1, 4) + 1
               id(2, 4) = id(1, 4) + 1
            end do
         end if

         ! coefficients:
         if (id(2, 4) .gt. nb(4)) then
            id(1:2, 4) = nb(4)
            co(2, 4) = 0.0
         else
            co(2, 4) = (tf(l) - tb(id(1, 4)))/(tb(id(2, 4)) - tb(id(1, 4)))
         end if
         co(1, 4) = 1.0 - co(2, 4)

         ! initialize grid level 3:
         id(1, 3) = 1
         id(2, 3) = min0(2, nb(3))
         do k = 1, nf(3)

            ! fine grid points passes bk grid points:
            if (nb(3) .eq. 1) then                ! only one gridpoint case
               id(2, 3) = id(1, 3) + 1
            else
               do while ((id(1, 3) .lt. nb(3)) .and. ( &
                         ((zb(1) .gt. zb(2)) .and. (zf(k) .le. zb(id(2, 3)))) .or. &
                         ((zb(1) .lt. zb(2)) .and. (zf(k) .ge. zb(id(2, 3))))))
                  id(1, 3) = id(1, 3) + 1
                  id(2, 3) = id(1, 3) + 1
               end do
            end if

            ! coefficients:
            if (id(2, 3) .gt. nb(3)) then
               id(1:2, 3) = nb(3)
               co(2, 3) = 0.0
            else
               co(2, 3) = (zf(k) - zb(id(1, 3)))/(zb(id(2, 3)) - zb(id(1, 3)))
            end if
            co(1, 3) = 1.0 - co(2, 3)

            ! initialize grid level 2:
            id(1, 2) = 1
            id(2, 2) = min0(2, nb(2))
            do j = 1, nf(2)

               ! fine grid points passes bk grid points:
               if (nb(2) .eq. 1) then        ! only one gridpoint case
                  id(2, 2) = id(1, 2) + 1
               else
                  do while ((id(1, 2) .lt. nb(2)) .and. (yf(j) .ge. yb(id(2, 2))))
                     id(1, 2) = id(1, 2) + 1
                     id(2, 2) = id(1, 2) + 1
                  end do
               end if

               ! coefficients:
               if (id(2, 2) .gt. nb(2)) then
                  id(1:2, 2) = nb(2)
                  co(2, 2) = 0.0
               else
                  co(2, 2) = (yf(j) - yb(id(1, 2)))/(yb(id(2, 2)) - yb(id(1, 2)))
               end if
               co(1, 2) = 1.0 - co(2, 2)

               ! initialize grid level 2:
               id(1, 1) = 1
               id(2, 1) = min0(2, nb(1))
               do i = 1, nf(1)

                  ! fine grid points passes bk grid points:
                  if (nb(1) .eq. 1) then        ! only one gridpoint case
                     id(2, 1) = id(1, 1) + 1
                  else
                     do while ((id(1, 1) .lt. nb(1)) .and. (xf(i) .ge. xb(id(2, 1))))
                        id(1, 1) = id(1, 1) + 1
                        id(2, 1) = id(1, 1) + 1
                     end do
                  end if

                  ! coefficients:
                  if (id(2, 1) .gt. nb(1)) then
                     id(1:2, 1) = nb(1)
                     co(2, 1) = 0.0
                  else
                     co(2, 1) = (xf(i) - xb(id(1, 1)))/(xb(id(2, 1)) - xb(id(1, 1)))
                  end if
                  co(1, 1) = 1.0 - co(2, 1)

                  fg(i, j, k, l, 1:ns) = &
                     bk(id(1, 1), id(1, 2), id(1, 3), id(1, 4), 1:ns)*co(1, 1)*co(1, 2)*co(1, 3)*co(1, 4) + &
                     bk(id(2, 1), id(1, 2), id(1, 3), id(1, 4), 1:ns)*co(2, 1)*co(1, 2)*co(1, 3)*co(1, 4) + &
                     bk(id(1, 1), id(2, 2), id(1, 3), id(1, 4), 1:ns)*co(1, 1)*co(2, 2)*co(1, 3)*co(1, 4) + &
                     bk(id(1, 1), id(1, 2), id(2, 3), id(1, 4), 1:ns)*co(1, 1)*co(1, 2)*co(2, 3)*co(1, 4) + &
                     bk(id(1, 1), id(1, 2), id(1, 3), id(2, 4), 1:ns)*co(1, 1)*co(1, 2)*co(1, 3)*co(2, 4) + &
                     bk(id(2, 1), id(2, 2), id(1, 3), id(1, 4), 1:ns)*co(2, 1)*co(2, 2)*co(1, 3)*co(1, 4) + &
                     bk(id(2, 1), id(1, 2), id(2, 3), id(1, 4), 1:ns)*co(2, 1)*co(1, 2)*co(2, 3)*co(1, 4) + &
                     bk(id(2, 1), id(1, 2), id(1, 3), id(2, 4), 1:ns)*co(2, 1)*co(1, 2)*co(1, 3)*co(2, 4) + &
                     bk(id(1, 1), id(2, 2), id(2, 3), id(1, 4), 1:ns)*co(1, 1)*co(2, 2)*co(2, 3)*co(1, 4) + &
                     bk(id(1, 1), id(2, 2), id(1, 3), id(2, 4), 1:ns)*co(1, 1)*co(2, 2)*co(1, 3)*co(2, 4) + &
                     bk(id(1, 1), id(1, 2), id(2, 3), id(2, 4), 1:ns)*co(1, 1)*co(1, 2)*co(2, 3)*co(2, 4) + &
                     bk(id(2, 1), id(2, 2), id(2, 3), id(1, 4), 1:ns)*co(2, 1)*co(2, 2)*co(2, 3)*co(1, 4) + &
                     bk(id(2, 1), id(2, 2), id(1, 3), id(2, 4), 1:ns)*co(2, 1)*co(2, 2)*co(1, 3)*co(2, 4) + &
                     bk(id(2, 1), id(1, 2), id(2, 3), id(2, 4), 1:ns)*co(2, 1)*co(1, 2)*co(2, 3)*co(2, 4) + &
                     bk(id(1, 1), id(2, 2), id(2, 3), id(2, 4), 1:ns)*co(1, 1)*co(2, 2)*co(2, 3)*co(2, 4) + &
                     bk(id(2, 1), id(2, 2), id(2, 3), id(2, 4), 1:ns)*co(2, 1)*co(2, 2)*co(2, 3)*co(2, 4)
               end do
            end do
         end do
      end do
   end subroutine bkgtofine

   subroutine direction(u, v, dg)
!*************************************************
! calculate direction (general) (not used)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      real  :: u, v, dg, ra
! --------------------
      ra = atan(1.0d0)/45.0d0
      if (u .eq. 0.0d0) then
         if (v .gt. 0.0d0) then
            dg = 90.0d0
         else
            dg = 270.0d0
         end if
      else
         dg = atan2(v, u)/ra
         if (dg .lt. 0.0d0) dg = dg + 360.0d0
      end if
      return
   end subroutine direction

   subroutine getobdate(ad, dh, bd)
!*************************************************
! calculate observation date, modified from 'raob2dwl.f' (general) (not used)
! history: september 2007, coded by wei li.
!          march 2008, modified by zhongjie he.
!*************************************************
      implicit none
! --------------------
      integer  :: mn(12)
      integer  :: dt, iy, im, id
      real  :: ad, bd, hr, dh
      logical      :: fg
      data mn/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!  real ,external :: relday
! --------------------
      dt = nint(ad)
      iy = mod(dt/1000000, 100)
      im = mod(dt/10000, 100)
      id = mod(dt/100, 100)
      hr = ad - dt/100*100 + dh
      if (mod(iy, 4) .eq. 0) mn(2) = 29
      fg = .true.
      do while (fg)
         if (hr .lt. 0) then
            hr = hr + 24
            id = id - 1
            if (id .eq. 0) then
               im = im - 1
               if (im .eq. 0) then
                  im = 12
                  iy = iy - 1
                  if (iy .lt. 0) iy = 99
               end if
               id = mn(im)
            end if
         elseif (hr .ge. 24) then
            hr = hr - 24
            id = id + 1
            if (id .gt. mn(im)) then
               id = 1
               im = im + 1
               if (im .gt. 12) then
                  im = 1
                  iy = mod(iy + 1, 100)
               end if
            end if
         else
            fg = .false.
         end if
      end do
!  bd = iy*1000000 + im*10000 + id*100 + hr
!  if(abs(bd-1.0d0*idnint(bd)).le..01)bd=1.0d0*idnint(bd)
!  iy=iy+2000
!  bd=relday(0,0,0,id,im,iy)-relday(0,0,0,1,1,2000)
!  bd=bd+hr/24.0
      return
   end subroutine getobdate

   subroutine obstogrid(y0, x0, y00, x00, im, jm, x, y, is)
!*************************************************
! calculate the position of observation at the model grid.
! history: feburary 2008, coded by zhongjie he.
!*************************************************
      implicit none
      integer  :: i, j
      integer, intent(in)  :: im, jm
      real, intent(in)  :: x0, y0
      real, intent(in)  :: x00(im, jm), y00(im, jm)
      real, intent(out) :: x, y
      integer, intent(out) :: is
!--------------------
      is = 0
      x = 0
      y = 0

      if (x0 .le. x00(1, 1)) then
         if (abs(x00(1, 1) - x00(2, 1)) .le. 1e-10) return
         x = 1 - (x0 - x00(1, 1))/(x00(1, 1) - x00(2, 1))
      end if
      if (x0 .ge. x00(im, jm)) then
         if (abs(x00(im, jm) - x00(im - 1, jm)) .le. 1e-10) return
         x = im + (x0 - x00(im, jm))/(x00(im, jm) - x00(im - 1, jm))
      end if
      do i = 2, im
         if (x0 .ge. x00(i - 1, 1) .and. x0 .lt. x00(i, 1)) then
            if (abs(x00(i - 1, 1) - x00(i, 1)) .le. 1e-10) return
            x = i - 1 + (x0 - x00(i - 1, 1))/(x00(i, 1) - x00(i - 1, 1))
         end if
      end do

      if (y0 .le. y00(1, 1)) then
         if (abs(y00(1, 1) - y00(1, 2)) .le. 1e-10) return
         y = 1 - (y0 - y00(1, 1))/(y00(1, 1) - y00(1, 2))
      end if
      if (y0 .ge. y00(im, jm)) then
         if (abs(y00(im, jm) - y00(im, jm - 1)) .le. 1e-10) return
         y = jm + (y0 - y00(im, jm))/(y00(im, jm) - y00(im, jm - 1))
      end if
      do j = 2, jm
         if (y0 .ge. y00(1, j - 1) .and. y0 .lt. y00(1, j)) then
            if (abs(y00(1, j - 1) - y00(1, j)) .le. 1e-10) return
            y = j - 1 + (y0 - y00(1, j - 1))/(y00(1, j) - y00(1, j - 1))
         end if
      end do

      if (x .ge. 1 .and. x .le. im .and. y .ge. 1 .and. y .le. jm) is = 1

      return
   end subroutine obstogrid

end module generaltools
