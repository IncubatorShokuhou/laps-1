module stmas4d_core

   use prmtrs_stmas
   use generaltools
   use prep_stmas4d
   use post_stmas4d
   use costfun_grad
   use wcompt_gradt, only: wcompgernl

contains

   subroutine mganalyss
!*************************************************
! multi-grid analysis
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
!------------------------
      integer       :: it, fw
!                  fw=1: from coarse grid to fine grid
!                  fw=0: from fine grid to coarse grid
      integer  :: i, j, k, t, s

      integer  ::num, no, n, o   !shuyuan 20100901
      integer  :: np(maxdims)  !shuyuan
!------------------------

      print *, 'preprocss'
      call preprocss

      grdlevl = 0
      it = 0
      fw = 1
      num = 0
      do while (.true.)

         grdlevl = grdlevl + 1
         print *, 'number', grdlevl, 'level grid is in processing'
         call rdinitbgd
         print *, 'physcpstn'
         call physcpstn
         print *, 'getcoefft'
         call getcoefft

         if (numvars .eq. 0) exit

         print *, 'minimizer'
         call minimizer_xie
         if (grdlevl .eq. fnstgrd .and. it .eq. itrepet) exit

         if (ifrepet .eq. 0) then
            print *, 'coas2fine'
            call coas2fine_xie                ! use xie's as the old has error
         else
            if (grdlevl .eq. fnstgrd) fw = 0
            if (grdlevl .eq. 1) then
               fw = 1
               it = it + 1
            end if
            if (fw .eq. 1) then
               print *, 'coas2fine'
               call coas2fine_xie                ! use xie's as the old has error
            else
               grdlevl = grdlevl - 2
               print *, 'fine2coas'
               call fine2coas_xie                ! use xie's as the old has error
            end if
         end if

      end do
!   call getw
      print *, 'pstprocss'
      call pstprocss
      return
   end subroutine mganalyss

   subroutine tmpmemalc
!*************************************************
! memory allocate for tmp array
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, er
! --------------------
      allocate (tmpanals(ntmpgrd(1), ntmpgrd(2), ntmpgrd(3), ntmpgrd(4), numstat), stat=er)
      if (er .ne. 0) stop 'tmpanals allocate wrong'
      do t = 1, ntmpgrd(4)
      do k = 1, ntmpgrd(3)
      do j = 1, ntmpgrd(2)
      do i = 1, ntmpgrd(1)
         do s = 1, numstat
            tmpanals(i, j, k, t, s) = grdanals(i, j, k, t, s)
         end do
      end do
      end do
      end do
      end do
      return
   end subroutine tmpmemalc

   subroutine tmpmemrls
!*************************************************
! memory release for tmp array
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
      deallocate (tmpanals)
      return
   end subroutine tmpmemrls

   subroutine getcoefft
!*************************************************
! get grid point coefficent for every observation
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, o, m, n, uv, rr
      integer  :: lm(maxdims)
      integer  :: np(maxdims), nn(maxdims)
      real     :: ac(numdims, ngptobs), oc(numdims), co(ngptobs)
! --------------------
! get interpolation coefficent
      uv = 0
      rr = 0
      if (nallobs .eq. 0) return
      do o = 1, nallobs
         do n = 1, maxdims
            np(n) = 1
         end do
         do n = 1, numdims
            np(n) = int((obspostn(n, o) - oripstn(n))/grdspac(n)) + 1
            if (np(n) .eq. numgrid(n) .and. numgrid(n) .ne. 1) np(n) = numgrid(n) - 1
            if (numgrid(n) .eq. 1) np(n) = 1
         end do
         do n = 1, maxdims
            obsidxpc(n, o) = np(n)
         end do
         m = 0
!====================================================
!    do t=np(4),min0(np(4)+1,numgrid(4))
!    do k=np(3),min0(np(3)+1,numgrid(3))
!    do j=np(2),min0(np(2)+1,numgrid(2))
!    do i=np(1),min0(np(1)+1,numgrid(1))
!      nn(1)=i
!      nn(2)=j
!      nn(3)=k
!      nn(4)=t
!      m=m+1
!      do n=1,numdims
!        ac(n,m)=nn(n)*1.0
!      enddo
!    enddo
!    enddo
!    enddo
!    enddo
!    do n=1,numdims
!      oc(n)=(obspostn(n,o)-oripstn(n))/grdspac(n)+1
!    enddo
!================== modified by zhongjie he ========
         do n = 1, maxdims
            lm(n) = np(n) + 1
            if (numdims .lt. n) lm(n) = np(n)
         end do
         do t = np(4), lm(4)
         do k = np(3), lm(3)
         do j = np(2), lm(2)
         do i = np(1), lm(1)
            nn(1) = min0(i, numgrid(1))
            nn(2) = min0(j, numgrid(2))
            nn(3) = min0(k, numgrid(3))
            nn(4) = min0(t, numgrid(4))
            m = m + 1
            do n = 1, numdims
               ac(n, m) = nn(n)*1.0
            end do
         end do
         end do
         end do
         end do
         do n = 1, numdims
            if (numgrid(n) .ge. 2) then
               oc(n) = (obspostn(n, o) - oripstn(n))/grdspac(n) + 1
            else
               oc(n) = 1
            end if
         end do
!====================================================
         ! call interpltn_xie(numdims,ngptobs,co,ac,oc,3,numgrid(3),ppm)
         call interpltn(numdims, ngptobs, co, ac, oc)
         do m = 1, ngptobs
            obscoeff(m, o) = co(m)
         end do
      end do
      uv = nobstat(u_cmpnnt) + nobstat(v_cmpnnt)
      rr = nobstat(numstat + 1)
      obsradar = 1.0
      if (rr .ne. 0) obsradar = 100.0*uv/float(rr)
      rr = nobstat(numstat + 2)
      obs_sfmr = 1.0
      if (rr .ne. 0) obs_sfmr = 1.0/rr
!jhui
      rr = nobstat(numstat + 3)
      if (rr .ne. 0) obsref = 1.0/rr

      return
   end subroutine getcoefft

   subroutine minimizer_xie
!*************************************************
! minimize the cost function
! history: august 2007, coded by wei li.
!
!          modified by yuanfu for using a single
!          precision lbfgs routine
!
!          december 2013 modified by yuanfu for
!          a) changing automatic arrays to allocatables
!          b) switching to lbfgsb 3.0 with yuanfu's
!             modification for super large minimization
!          c) removing initialization of iw and wa for
!             efficiency
!*************************************************
      implicit none
! --------------------
      integer, parameter :: mm = 5
      character(len=60) :: ta, cs
      logical           :: ls(4)
      integer    :: i0, ic, ip, it, isbmn, n, o, s, t, k, j, i, no, er
      integer    :: is(44)
      integer, allocatable :: nb(:), iw(:)
      integer    :: nn(maxdims + 1), ng(maxdims + 1), nc
      real       :: mn(numstat + 3), mx(numstat + 3)

      real, allocatable :: lb(:), ub(:)
      real :: fa, pg, ds(29)
      real, allocatable :: wa(:)

      real  ::temp, temp1, temp2, dif1, dif2, temp0
      integer::  locx, locy, locz, loct, locv
! --------------------

      ! allocate working array:
      allocate (wa(2*mm*numvars + 5*numvars + 12*mm*mm + 12*mm), &
                lb(numvars), ub(numvars), nb(numvars), iw(3*numvars), stat=er)
      if (er .ne. 0) then
         print *, 'minimizer: cannot allocate enough memory for the working array'
         stop
      else
         print *, 'successfully allocate memory for lbfgs!'
      end if

      ip = 1
      fa = 1.0d+2
      pg = 1.0d-10
      isbmn = 1        ! subspace minimization
      do n = 1, numvars
         nb(n) = 0
         lb(n) = 0.0
         ub(n) = 0.0
      end do

      ! bound: if ifbound eq 1, use (min(obs), max(obs)) to bound
      !        analysis
      if (ifbound .eq. 1 .and. nallobs .ge. 1) then

         ! find the min(obs) and max(obs):

         ! bound at least allow zero increment: by yuanfu correcting zhongjie's setting
         mn = 0.0
         mx = 0.0

         ! conventional obs:
         do n = 1, maxdims
            ng(n) = numgrid(n)
         end do
         ng(maxdims + 1) = numstat + 1
         o = 0
         do s = 1, numstat
            do no = 1, nobstat(s)
               o = o + 1
               mn(s) = min(mn(s), obsvalue(o))
               mx(s) = max(mx(s), obsvalue(o))
            end do
         end do

         ! radar and sfmr:
         s = numstat + 1
         do no = 1, nobstat(s)
            o = o + 1
            mn(u_cmpnnt) = min(mn(u_cmpnnt), obsvalue(o))
            mx(u_cmpnnt) = max(mx(u_cmpnnt), obsvalue(o))
            mn(v_cmpnnt) = min(mn(v_cmpnnt), obsvalue(o))
            mx(v_cmpnnt) = max(mx(v_cmpnnt), obsvalue(o))
         end do
         s = numstat + 2
         do no = 1, nobstat(s)
            o = o + 1
            mn(u_cmpnnt) = min(mn(u_cmpnnt), -1.*obsvalue(o))
            mx(u_cmpnnt) = max(mx(u_cmpnnt), obsvalue(o))
            mn(v_cmpnnt) = min(mn(v_cmpnnt), -1.*obsvalue(o))
            mx(v_cmpnnt) = max(mx(v_cmpnnt), obsvalue(o))
         end do
!jhui
!----------------------------------
!changed by shuyuan 20100722
!      mn(6) =0.0
!      mx(6) =0.0

         s = numstat + 3
         mn(s) = 0.0
         mx(s) = 0.0
         do no = 1, nobstat(s)
            o = o + 1
            mn(s) = min(mn(s), obsvalue(o))
            mx(s) = max(mx(s), obsvalue(o))
         end do

         do s = 6, 7
            mn(s) = 0.0
            mx(s) = 0.0
            do no = 1, nobstat(s)
               o = o + 1
               mn(s) = 0.
               mx(s) = 1000.!just for test
            end do
         end do
!--------------------------------------------------
         do s = 1, 5!!numstat+1
            do t = 1, numgrid(4)
            do k = 1, numgrid(3)
            do j = 1, numgrid(2)
            do i = 1, numgrid(1)
               nn(1) = i
               nn(2) = j
               nn(3) = k
               nn(4) = t
               nn(maxdims + 1) = s
               call pstn2numb(maxdims + 1, nn, ng, nc)
               if (mn(s) .le. mx(s)) then
                  nb(nc) = 2
                  lb(nc) = mn(s)
                  ub(nc) = mx(s)
               else
                  nb(nc) = 0
                  lb(nc) = 0.0d0
                  ub(nc) = 0.0d0
               end if
            end do
            end do
            end do
            end do
         end do

         if (ifbkgnd .eq. 1) then
            do n = 1, numvars
               if (nb(n) .eq. 2) then
                  lb(n) = min(lb(n), 0.0)
                  ub(n) = max(ub(n), 0.0)
               end if
            end do
         end if

      end if

! bound
! for 3d radar
      if (w_cmpnnt .ne. 0) then
         do n = 1, maxdims
            ng(n) = numgrid(n)
         end do
         ng(maxdims + 1) = numstat
         s = w_cmpnnt
         k = 1
         do t = 1, numgrid(4)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            nn(1) = i
            nn(2) = j
            nn(3) = k
            nn(4) = t
            nn(maxdims + 1) = s
            call pstn2numb(maxdims + 1, nn, ng, nc)
            nb(nc) = 2
            lb(nc) = 0.0
            ub(nc) = 0.0
         end do
         end do
         end do
      end if

      ta = 'start'
      i0 = 0
      ic = 0
      write (*, 9001)

      if (w_cmpnnt .ne. 0) then         ! by zhongjie he
         call costfunct1
         call costgradt1
      else                           ! by zhongjie he
         call wcompgernl

         call costfunct2
         call costgradt2
      end if
      print *, 'minvalue of bk: ', minval(grdbkgnd(1:numgrid(1), 1:numgrid(2), 1:numgrid(3), 1:numgrid(4), 5))

      ! count number of iterations:
      it = 0
      iterloop: do while (.true.)

         ! add low bound for specific humidity:
         ng(1:4) = numgrid(1:4)
         ng(5) = numstat  ! the number of states in the 5th dimension
         do t = 1, numgrid(4)
            do k = 1, numgrid(3)
               do j = 1, numgrid(2)
                  do i = 1, numgrid(1)
                     nn(1) = i
                     nn(2) = j
                     nn(3) = k
                     nn(4) = t
                     nn(5) = 5  ! position of specific humidity

                     call pstn2numb(5, nn, ng, nc)

                     nb(nc) = 0
                     if ((maxgrid(1) - 1)/2 + 1 .le. numgrid(1)) nb(nc) = 1

                     ! low bound from reflectivity did not get scaled and so here it does:
                     lb(nc) = -grdbkgnd(i, j, k, t, 5) + grdbkgnd(i, j, k, t, numstat + 1)
                     ub(nc) = -grdbkgnd(i, j, k, t, 5) + grdbkgnd(i, j, k, t, numstat + 2)
                  end do
               end do
            end do
         end do
         ! add lower bounds for rain and snow:
         if (numstat .gt. 5) then ! rain and snow are 6th and 7th:
         do t = 1, numgrid(4)
            do k = 1, numgrid(3)
               do j = 1, numgrid(2)
                  do i = 1, numgrid(1)
                     nn(1) = i
                     nn(2) = j
                     nn(3) = k
                     nn(4) = t
                     nn(5) = rour_cmpnnt  ! position of rain

                     call pstn2numb(5, nn, ng, nc)

                     nb(nc) = 1
                     lb(nc) = -grdbkgnd(i, j, k, t, rour_cmpnnt)

                     nn(5) = rous_cmpnnt  ! position of snow

                     call pstn2numb(5, nn, ng, nc)

                     nb(nc) = 1
                     lb(nc) = -grdbkgnd(i, j, k, t, rous_cmpnnt)
                  end do
               end do
            end do
         end do
         end if

         call setulb(numvars, mm, grdanals, lb, ub, nb, costfun, gradint, fa, pg, wa, iw, &
                     ta, ip, cs, ls, is, ds)

         if (ta(1:2) .eq. 'fg') then
            if (w_cmpnnt .ne. 0) then         ! by zhongjie he
               call costfunct1
               call costgradt1
            else                           ! by zhongjie he
               call wcompgernl
               call costfunct2
               call costgradt2
            end if
            ic = ic + 1
            cycle iterloop
         elseif (ta(1:5) .eq. 'new_x') then
            if (i0 .eq. 0) then
               write (*, 1003) is(30), is(34), ds(13), costfun
            else
               write (*, 3001) is(30), is(34), ds(14), ds(13), costfun
            end if
            i0 = i0 + 1
            it = it + 1

            ! exceeds maximum iterations:
            if (grdlevl .le. midgrid .and. it .gt. cosstep) exit iterloop
            if (grdlevl .gt. midgrid .and. it .gt. finstep) exit iterloop

            cycle iterloop
         else
            if (ip .le. -1 .and. ta(1:4) .ne. 'stop') write (6, *) ta
            exit iterloop
         end if
      end do iterloop
1003  format(2(1x, i4), 5x, '-', 3x, 1p, 2(1x, d10.3))
9001  format(/, 3x, 'it', 3x, 'nf', 2x, 'stepl', 5x, 'projg', 8x, 'f')
3001  format(2(1x, i4), 1p, 2x, d8.1, 1p, 2(1x, d10.3))

      if (w_cmpnnt .ne. 0) then         ! by zhongjie he
         call costfunct1
!    call costgradt1
      else                           ! by zhongjie he
         call wcompgernl
         call costfunct2

!    call costgradt2
      end if

      ! deallocate working array:
      deallocate (wa, lb, ub, nb, iw, stat=er)
      if (er .ne. 0) then
         print *, 'minimizer: cannot deallocate enough memory for the working array'
         stop
      end if

      write (8, *) 'bottom function =', costfun
      return
   end subroutine minimizer_xie

   subroutine coas2fine_xie
!*************************************************
! interpolation from coarse grid to fine grid
! history: august 2007, coded by wei li.
!                modified by yuanfu from wei li's
!                coas2fine,which is incorrectly
!                interpolating coase to fine grid.
!
!        note: this routine assume fine grid halves
!                the coase resolution.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, n, i0, j0, k0, t0, sm
      integer  :: ic(maxdims)
! --------------------
      do n = 1, maxdims
         ntmpgrd(n) = numgrid(n)
      end do
      call tmpmemalc
      call grdmemrls
      do n = 1, maxdims
         ic(n) = 1
         if (ntmpgrd(n) .lt. maxgrid(n)) then
!      if(n.ne.3 .or. grdlevl.gt.3) then         ! added by zhongjie he
            ic(n) = 2
            numgrid(n) = 2*numgrid(n) - 1
!      endif
         end if
      end do

      ! pressure coordinate levels:
      i = (maxgrid(pressure) - 1)/(numgrid(pressure) - 1)
      do k = 1, numgrid(pressure)
         ppm(k) = pp0((k - 1)*i + 1)
      end do

      call grdmemalc
! get the grid spacing information for new grd array
      do n = 1, numdims
         if (ic(n) .eq. 2) grdspac(n) = grdspac(n)/2.0
      end do

! note fine grid doubles resolution from coase grid:

! 1. project coase grid onto fine grid:
      do t = 1, numgrid(4), ic(4)
         t0 = (t - 1)/ic(4) + 1
         do k = 1, numgrid(3), ic(3)
            k0 = (k - 1)/ic(3) + 1
            do j = 1, numgrid(2), ic(2)
               j0 = (j - 1)/ic(2) + 1
               do i = 1, numgrid(1), ic(1)
                  i0 = (i - 1)/ic(1) + 1

                  grdanals(i, j, k, t, 1:numstat) = tmpanals(i0, j0, k0, t0, 1:numstat)
               end do
            end do
         end do
      end do

      ! release temporary memory:
      call tmpmemrls

! 2. x direction:
      if (ic(1) .eq. 2) then
         do t = 1, numgrid(4), ic(4)
         do k = 1, numgrid(3), ic(3)
         do j = 1, numgrid(2), ic(2)
         do i = 2, numgrid(1), ic(1)

            grdanals(i, j, k, t, 1:numstat) = 0.5*(grdanals(i - 1, j, k, t, 1:numstat) &
                                                   + grdanals(i + 1, j, k, t, 1:numstat))
         end do
         end do
         end do
         end do
      end if

! 3. y direction:
      if (ic(2) .eq. 2) then
         do t = 1, numgrid(4), ic(4)
         do k = 1, numgrid(3), ic(3)
         do j = 2, numgrid(2), ic(2)
         do i = 1, numgrid(1)

            grdanals(i, j, k, t, 1:numstat) = 0.5*(grdanals(i, j - 1, k, t, 1:numstat) &
                                                   + grdanals(i, j + 1, k, t, 1:numstat))
         end do
         end do
         end do
         end do
      end if

! 4. z direction:
      if (ic(3) .eq. 2) then
         do t = 1, numgrid(4), ic(4)
         do k = 2, numgrid(3), ic(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)

            grdanals(i, j, k, t, 1:numstat) = 0.5*(grdanals(i, j, k - 1, t, 1:numstat) &
                                                   + grdanals(i, j, k + 1, t, 1:numstat))
         end do
         end do
         end do
         end do
      end if

! 5. t direction:
      if (ic(4) .eq. 2) then
         do t = 2, numgrid(4), ic(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)

            grdanals(i, j, k, t, 1:numstat) = 0.5*(grdanals(i, j, k, t - 1, 1:numstat) &
                                                   + grdanals(i, j, k, t + 1, 1:numstat))
         end do
         end do
         end do
         end do
      end if

      ! make sure interpolated sh analysis positive:
      do t = 1, numgrid(4)
      do k = 1, numgrid(3)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         grdanals(i, j, k, t, humidity) = &
            max(-grdbkgnd(i, j, k, t, humidity), grdanals(i, j, k, t, humidity))
      end do
      end do
      end do
      end do

      return
   end subroutine coas2fine_xie

   subroutine coas2fine
!*************************************************
! interpolation from coarse grid to fine grid
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, n, i0, j0, k0, t0, sm
      integer  :: ic(maxdims)
! --------------------
      do n = 1, maxdims
         ntmpgrd(n) = numgrid(n)
      end do
      call tmpmemalc
      call grdmemrls
      do n = 1, maxdims
         ic(n) = 1
         if (ntmpgrd(n) .lt. maxgrid(n)) then
!      if(n.ne.3 .or. grdlevl.gt.3) then         ! added by zhongjie he
            ic(n) = 2
            numgrid(n) = 2*numgrid(n) - 1
!      endif
         end if
      end do
      call grdmemalc
! get the grid spacing information for new grd array
      do n = 1, numdims
         if (ic(n) .eq. 2) grdspac(n) = grdspac(n)/2.0
      end do
! get the information for new grid array by interpolation
!======================== modified by zhongjie he
      do t = 1, numgrid(4), ic(4)
      do k = 1, numgrid(3), ic(3)
      do j = 1, numgrid(2), ic(2)
      do i = 1, numgrid(1), ic(1)
         i0 = 0.5*(i + 1)
         j0 = 0.5*(j + 1)
         k0 = 0.5*(k + 1)
         t0 = 0.5*(t + 1)
         if (ic(1) .eq. 1) i0 = i
         if (ic(2) .eq. 1) j0 = j
         if (ic(3) .eq. 1) k0 = k
         if (ic(4) .eq. 1) t0 = t

         do s = 1, numstat
            sm = 0
            if (ic(1) .eq. 2) then
               if (i0 .eq. 1) then
                  grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0 + 1, j0, k0, t0, s)
                  sm = sm + 1
               elseif (i0 .eq. ntmpgrd(1)) then
                  grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0 - 1, j0, k0, t0, s)
                  sm = sm + 1
               else
               grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0 - 1, j0, k0, t0, s) + tmpanals(i0 + 1, j0, k0, t0, s)
                  sm = sm + 2
               end if
            end if
            if (ic(2) .eq. 2) then
               if (j0 .eq. 1) then
                  grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0, j0 + 1, k0, t0, s)
                  sm = sm + 1
               elseif (j0 .eq. ntmpgrd(2)) then
                  grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0, j0 - 1, k0, t0, s)
                  sm = sm + 1
               else
               grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0, j0 - 1, k0, t0, s) + tmpanals(i0, j0 + 1, k0, t0, s)
                  sm = sm + 2
               end if
            end if
            if (ic(3) .eq. 2) then
               if (k0 .eq. 1) then
                  grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0, j0, k0 + 1, t0, s)
                  sm = sm + 1
               elseif (k0 .eq. ntmpgrd(3)) then
                  grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0, j0, k0 - 1, t0, s)
                  sm = sm + 1
               else
               grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0, j0, k0 - 1, t0, s) + tmpanals(i0, j0, k0 + 1, t0, s)
                  sm = sm + 2
               end if
            end if
            if (ic(4) .eq. 2) then
               if (t0 .eq. 1) then
                  grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0, j0, k0, t0 + 1, s)
                  sm = sm + 1
               elseif (t0 .eq. ntmpgrd(4)) then
                  grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0, j0, k0, t0 - 1, s)
                  sm = sm + 1
               else
               grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0, j0, k0, t0 - 1, s) + tmpanals(i0, j0, k0, t0 + 1, s)
                  sm = sm + 2
               end if
            end if

            if (sm .ge. 1) then
               grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s)/float(sm)
            else
               grdanals(i, j, k, t, s) = tmpanals(i0, j0, k0, t0, s)
            end if
         end do
      end do
      end do
      end do
      end do
!========================================================== end of modification by zhongjie he
! x direction
      if (ic(1) .eq. 2) then
         do t = 1, numgrid(4), ic(4)
         do k = 1, numgrid(3), ic(3)
         do j = 1, numgrid(2), ic(2)
         do i = 2, numgrid(1) - 1, ic(1)
            do s = 1, numstat
               grdanals(i, j, k, t, s) = 0.5*(grdanals(i - 1, j, k, t, s) + grdanals(i + 1, j, k, t, s))
            end do
         end do
         end do
         end do
         end do
      end if
! y direction
      if (ic(2) .eq. 2) then
         do t = 1, numgrid(4), ic(4)
         do k = 1, numgrid(3), ic(3)
         do j = 2, numgrid(2) - 1, ic(2)
         do i = 1, numgrid(1)
            do s = 1, numstat
               grdanals(i, j, k, t, s) = 0.5*(grdanals(i, j - 1, k, t, s) + grdanals(i, j + 1, k, t, s))
            end do
         end do
         end do
         end do
         end do
      end if
! z direction
      if (ic(3) .eq. 2) then
         do t = 1, numgrid(4), ic(4)
         do k = 2, numgrid(3) - 1, ic(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            do s = 1, numstat
               grdanals(i, j, k, t, s) = 0.5*(grdanals(i, j, k - 1, t, s) + grdanals(i, j, k + 1, t, s))
            end do
         end do
         end do
         end do
         end do
      end if
! t direction
      if (ic(4) .eq. 2) then
         do t = 2, numgrid(4) - 1, ic(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            do s = 1, numstat
               grdanals(i, j, k, t, s) = 0.5*(grdanals(i, j, k, t - 1, s) + grdanals(i, j, k, t + 1, s))
            end do
         end do
         end do
         end do
         end do
      end if
      call tmpmemrls
      return
   end subroutine coas2fine

   subroutine fine2coas_xie
!*************************************************
! interpolation from coarse grid to fine grid
! history: august 2007, coded by wei li.
!                modified by yuanfu from wei li's
!                coas2fine,which is incorrectly
!                interpolating coase to fine grid.
!                thus yuanfu rewrites the fine2coas.
!
!        note: this routine assume fine grid halves
!                the coase resolution.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, n, i0, j0, k0, t0, sm
      integer  :: ic(maxdims)
! --------------------
      do n = 1, maxdims
         ntmpgrd(n) = numgrid(n)
      end do
      call tmpmemalc                ! temporary memory and copy grdanals to here
      call grdmemrls
      do n = 1, maxdims
         ic(n) = 1
         if (ntmpgrd(n) .gt. inigrid(n)) then
            ic(n) = 2
            numgrid(n) = (numgrid(n) - 1)/2 + 1                        ! assume numgrid odd number
         end if
      end do
      call grdmemalc
! get the grid spacing information for new grd array
      do n = 1, numdims
         if (ic(n) .eq. 2) grdspac(n) = grdspac(n)*2.0
      end do

! note fine grid doubles resolution from coase grid:

! 1. project coase grid onto fine grid:
      do t = 1, numgrid(4)
         t0 = (t - 1)*ic(4) + 1
         do k = 1, numgrid(3)
            k0 = (k - 1)*ic(3) + 1
            do j = 1, numgrid(2)
               j0 = (j - 1)*ic(2) + 1
               do i = 1, numgrid(1)
                  i0 = (i - 1)*ic(1) + 1

                  grdanals(i, j, k, t, 1:numstat) = tmpanals(i0, j0, k0, t0, 1:numstat)
               end do
            end do
         end do
      end do

      ! release temporary memory:
      call tmpmemrls

      return
   end subroutine fine2coas_xie

   subroutine fine2coas
!*************************************************
! projection from fine grid to coarse grid (not used)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, n, m, i0, j0, k0, t0, sm
      integer  :: ic(maxdims)
! --------------------
      do n = 1, maxdims
         ntmpgrd(n) = numgrid(n)
      end do
      call tmpmemalc
      do n = 1, maxdims
         ic(n) = 1
         if (ntmpgrd(n) .gt. inigrid(n)) ic(n) = 2
      end do

!======================== modified by zhongjie he
      do t = 1, ntmpgrd(4)
      do k = 1, ntmpgrd(3)
      do j = 1, ntmpgrd(2)
      do i = 1, ntmpgrd(1)
         do s = 1, numstat
            sm = 0
            tmpanals(i, j, k, t, s) = 0
            if (ic(1) .eq. 2) then
               if (i .eq. 1) then
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i + 1, j, k, t, s)
                  sm = sm + 1
               elseif (i .eq. ntmpgrd(1)) then
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i - 1, j, k, t, s)
                  sm = sm + 1
               else
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i - 1, j, k, t, s) + grdanals(i + 1, j, k, t, s)
                  sm = sm + 2
               end if
            end if
            if (ic(2) .eq. 2) then
               if (j .eq. 1) then
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i, j + 1, k, t, s)
                  sm = sm + 1
               elseif (j .eq. ntmpgrd(2)) then
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i, j - 1, k, t, s)
                  sm = sm + 1
               else
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i, j - 1, k, t, s) + grdanals(i, j + 1, k, t, s)
                  sm = sm + 2
               end if
            end if
            if (ic(3) .eq. 2) then
               if (k .eq. 1) then
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i, j, k + 1, t, s)
                  sm = sm + 1
               elseif (k .eq. ntmpgrd(3)) then
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i, j, k - 1, t, s)
                  sm = sm + 1
               else
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i, j, k - 1, t, s) + grdanals(i, j, k + 1, t, s)
                  sm = sm + 2
               end if
            end if
            if (ic(4) .eq. 2) then
               if (t .eq. 1) then
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i, j, k, t + 1, s)
                  sm = sm + 1
               elseif (t .eq. ntmpgrd(4)) then
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i, j, k, t - 1, s)
                  sm = sm + 1
               else
                  tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s) + grdanals(i, j, k, t - 1, s) + grdanals(i, j, k, t + 1, s)
                  sm = sm + 2
               end if
            end if

            if (sm .ge. 1) then
               tmpanals(i, j, k, t, s) = tmpanals(i, j, k, t, s)/float(sm)
            else
               tmpanals(i, j, k, t, s) = grdanals(i, j, k, t, s)
            end if
         end do
      end do
      end do
      end do
      end do
!========================================================== end of modification by zhongjie he

      call grdmemrls
      do n = 1, maxdims
         if (ic(n) .eq. 2) numgrid(n) = 0.5*(numgrid(n) + 1)
      end do
      call grdmemalc
      do n = 1, numdims
         if (ic(n) .eq. 2) grdspac(n) = grdspac(n)*2.0
      end do
! ----
      if (.false.) then
! ----
         do t = 1, ntmpgrd(4), ic(4)
         do k = 1, ntmpgrd(3), ic(3)
         do j = 1, ntmpgrd(2), ic(2)
         do i = 1, ntmpgrd(1), ic(1)
            do s = 1, numstat
               grdanals(i, j, k, t, s) = 0.0
               m = 0
               do t0 = max0(t - 1, 1), min0(t + 1, ntmpgrd(4))
               do k0 = max0(k - 1, 1), min0(k + 1, ntmpgrd(3))
               do j0 = max0(j - 1, 1), min0(j + 1, ntmpgrd(2))
               do i0 = max0(i - 1, 1), min0(i + 1, ntmpgrd(1))
                  if (i0 .eq. i .and. j0 .eq. j .and. k0 .eq. k .and. t0 .eq. t) cycle
                  m = m + 1
                  grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s) + tmpanals(i0, j0, k0, t0, s)
               end do
               end do
               end do
               end do
               if (m .gt. 0) grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s)/m
            end do
         end do
         end do
         end do
         end do
! ----
      else
! ----
         do t = 1, ntmpgrd(4), ic(4)
         do k = 1, ntmpgrd(3), ic(3)
         do j = 1, ntmpgrd(2), ic(2)
         do i = 1, ntmpgrd(1), ic(1)
            i0 = 0.5*(i + 1)
            j0 = 0.5*(j + 1)
            k0 = 0.5*(k + 1)
            t0 = 0.5*(t + 1)
            if (ic(1) .eq. 1) i0 = i
            if (ic(2) .eq. 1) j0 = j
            if (ic(3) .eq. 1) k0 = k
            if (ic(4) .eq. 1) t0 = t
            do s = 1, numstat
               grdanals(i0, j0, k0, t0, s) = tmpanals(i, j, k, t, s)
            end do
         end do
         end do
         end do
         end do
! ----
      end if
! ----
      call tmpmemrls
      return
   end subroutine fine2coas

   subroutine check_f_g
!*************************************************
! check whether the gradient match cost function (affiliate)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s
      real     :: control(numgrid(1), numgrid(2), numgrid(3), numgrid(4), numstat)
      real     :: ed, f2, f1, gg
! --------------------
      ed = 0.00001
      do s = 1, numstat
      do t = 1, numgrid(4)
      do k = 1, numgrid(3)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         control(i, j, k, t, s) = grdanals(i, j, k, t, s)
      end do
      end do
      end do
      end do
      end do
      call wcompgernl
      call costgradt2
      do s = 1, numstat
      do t = 1, numgrid(4)
      do k = 1, numgrid(3)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         print *, '       '
         print *, '       '
         print *, '       '
         print *, '       '
         print *, '       '
         grdanals(i, j, k, t, s) = control(i, j, k, t, s) + ed
         call costfunct2
         f2 = costfun
         grdanals(i, j, k, t, s) = control(i, j, k, t, s) - ed
         call costfunct2
         f1 = costfun
         grdanals(i, j, k, t, s) = control(i, j, k, t, s)
         gg = (f2 - f1)/(2.0*ed)
!    print*,f2,f1,gg
         if (abs(gradint(i, j, k, t, s) - gg) .gt. 1.0e-19) &
            print *, i, j, k, t, s, gradint(i, j, k, t, s), gg, gradint(i, j, k, t, s) - gg
         if (abs(gradint(i, j, k, t, s) - gg) .gt. 1000.0) stop
      end do
      end do
      end do
      end do
      end do
      return
   end subroutine check_f_g

   subroutine getw
      implicit none
! --------------------
      integer  :: i, j, k, t
! --------------------
      call wcompgernl      ! modified by zhongjie he
      do t = 1, numgrid(4)
      do k = 1, numgrid(3)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         grdanals(i, j, k, t, 3) = www(i, j, k, t) !*scl(u_cmpnnt)
      end do
      end do
      end do
      end do
      return
   end subroutine getw

end module stmas4d_core
