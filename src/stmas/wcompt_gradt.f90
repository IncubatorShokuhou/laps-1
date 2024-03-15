module wcompt_gradt
!*************************************************
! calculate the vertical velocity and it's gradients to control variables
! history: january 2008, separated from input_bg_obs module by zhongjie he.
!*************************************************

   use prmtrs_stmas

   public wcompgernl, wwgradient

!***************************************************
!!comment:
!   this module is used by the module of costfun_grad to calculate the vertical velocity and gradients to control variables.
!   subroutines:
!      wcompgernl: calculate vertical velocities at each grid points.
!      wwgradient: calculate the gradient of vertical velocities to control variables.
!***************************************************
contains

   subroutine wcompgernl
!*************************************************
! calculate the vertical velocity www for general coordinate
! history : january 2008, modified from wei li's program by zhongjie he
!
!*************************************************

      implicit none
! --------------------
      integer  :: i, j, k, t, i1, i2, j1, j2, uu, vv, zz, k2, k1
      real     :: ht(numgrid(1), numgrid(2), numgrid(3), numgrid(4)), sz
! ---------------------
! declare :
!          'ht' is the temporary array to save height of each grid point
!          'uu', 'vv' and 'zz' are indexes of the contral variables for u, v and height respectively
! --------------------

      uu = u_cmpnnt
      vv = v_cmpnnt
      zz = pressure
      do k = 1, numgrid(3)
      do t = 1, numgrid(4)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         www(i, j, k, t) = 0.0
      end do
      end do
      end do
      end do

      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then         ! for sigma and height coordinate
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            ht(i, j, k, t) = zzz(i, j, k, t)
         end do
         end do
         end do
         end do
         sz = scp(psl)
      else                                            ! for pressure coordinate
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            ht(i, j, k, t) = ppp(k)
         end do
         end do
         end do
         end do
         sz = scp(psl)
      end if

      do k = 2, numgrid(3)
      do t = 1, numgrid(4)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         i1 = max0(i - 1, 1)
         i2 = min0(i + 1, numgrid(1))
         j1 = max0(j - 1, 1)
         j2 = min0(j + 1, numgrid(2))

         www(i, j, k, t) = www(i, j, k - 1, t) - (ht(i, j, k, t) - ht(i, j, k - 1, t))*sz/2.0 &
                           *((grdanals(i2, j, k - 1, t, uu) - grdanals(i1, j, k - 1, t, uu) &
                              + grdbkgnd(i2, j, k - 1, t, uu) - grdbkgnd(i1, j, k - 1, t, uu)) &
                             /(xxx(i2, j) - xxx(i1, j))/scp(xsl)*1.0 &
                             + (grdanals(i, j2, k - 1, t, vv) - grdanals(i, j1, k - 1, t, vv) &
                                + grdbkgnd(i, j2, k - 1, t, vv) - grdbkgnd(i, j1, k - 1, t, vv)) &
                             /(yyy(i, j2) - yyy(i, j1))/scp(ysl)*1.0 &
                             + (grdanals(i2, j, k, t, uu) - grdanals(i1, j, k, t, uu) &
                                + grdbkgnd(i2, j, k, t, uu) - grdbkgnd(i1, j, k, t, uu)) &
                             /(xxx(i2, j) - xxx(i1, j))/scp(xsl)*1.0 &
                             + (grdanals(i, j2, k, t, vv) - grdanals(i, j1, k, t, vv) &
                                + grdbkgnd(i, j2, k, t, vv) - grdbkgnd(i, j1, k, t, vv)) &
                             /(yyy(i, j2) - yyy(i, j1))/scp(ysl)*1.0 &
                             - (grdanals(i, j2, k - 1, t, uu) - grdanals(i, j1, k - 1, t, uu) &
                                + grdbkgnd(i, j2, k - 1, t, uu) - grdbkgnd(i, j1, k - 1, t, uu)) &
                             /(yyy(i, j2) - yyy(i, j1))/scp(ysl)*0.0 &
                             + (grdanals(i2, j, k - 1, t, vv) - grdanals(i1, j, k - 1, t, vv) &
                                + grdbkgnd(i2, j, k - 1, t, vv) - grdbkgnd(i1, j, k - 1, t, vv)) &
                             /(xxx(i2, j) - xxx(i1, j))/scp(xsl)*0.0 &
                             - (grdanals(i, j2, k, t, uu) - grdanals(i, j1, k, t, uu) &
                                + grdbkgnd(i, j2, k, t, uu) - grdbkgnd(i, j1, k, t, uu)) &
                             /(yyy(i, j2) - yyy(i, j1))/scp(ysl)*0.0 &
                             + (grdanals(i2, j, k, t, vv) - grdanals(i1, j, k, t, vv) &
                                + grdbkgnd(i2, j, k, t, vv) - grdbkgnd(i1, j, k, t, vv)) &
                             /(xxx(i2, j) - xxx(i1, j))/scp(xsl)*0.0)

         if (ifpcdnt .eq. 0) then                                          !     for sigma
            www(i, j, k, t) = www(i, j, k, t) - (ht(i, j, k, t) - ht(i, j, k - 1, t))*sz/2.0 &
                              *((grdanals(i, j, k, t, uu) - grdanals(i, j, k - 1, t, uu) &
                                 + grdbkgnd(i, j, k, t, uu) - grdbkgnd(i, j, k - 1, t, uu)) &
                                /(ht(i, j, k, t) - ht(i, j, k - 1, t)) &
                                *(ht(i2, j, k, t) - ht(i1, j, k, t) + ht(i2, j, k - 1, t) - ht(i1, j, k - 1, t)) &
                                /(xxx(i2, j) - xxx(i1, j))/scp(xsl)*1.0 &
                                + (grdanals(i, j, k, t, vv) - grdanals(i, j, k - 1, t, vv) &
                                   + grdbkgnd(i, j, k, t, vv) - grdbkgnd(i, j, k - 1, t, vv)) &
                                /(ht(i, j, k, t) - ht(i, j, k - 1, t)) &
                                *(ht(i, j2, k, t) - ht(i, j1, k, t) + ht(i, j2, k - 1, t) - ht(i, j1, k - 1, t)) &
                                /(yyy(i, j2) - yyy(i, j1))/scp(ysl)*1.0 &
                                - (grdanals(i, j, k, t, uu) - grdanals(i, j, k - 1, t, uu) &
                                   + grdbkgnd(i, j, k, t, uu) - grdbkgnd(i, j, k - 1, t, uu)) &
                                /(ht(i, j, k, t) - ht(i, j, k - 1, t)) &
                                *(ht(i, j2, k, t) - ht(i, j1, k, t) + ht(i, j2, k - 1, t) - ht(i, j1, k - 1, t)) &
                                /(yyy(i, j2) - yyy(i, j1))/scp(ysl)*0.0 &
                                + (grdanals(i, j, k, t, vv) - grdanals(i, j, k - 1, t, vv) &
                                   + grdbkgnd(i, j, k, t, vv) - grdbkgnd(i, j, k - 1, t, vv)) &
                                /(ht(i, j, k, t) - ht(i, j, k - 1, t)) &
                                *(ht(i2, j, k, t) - ht(i1, j, k, t) + ht(i2, j, k - 1, t) - ht(i1, j, k - 1, t)) &
                                /(xxx(i2, j) - xxx(i1, j))/scp(xsl)*0.0)
         end if
      end do
      end do
      end do
      end do

      if (ifpcdnt .eq. 1) then                                          !     for pressure
         do k = 2, numgrid(3)
         do t = 1, numgrid(4)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            k1 = k - 1
            k2 = min(k + 1, numgrid(3))
            www(i, j, k, t) = www(i, j, k, t) &
                              *(grdanals(i, j, k2, t, zz) - grdanals(i, j, k1, t, zz) &
                                + grdbkgnd(i, j, k2, t, zz) - grdbkgnd(i, j, k1, t, zz)) &
                              /(ht(i, j, k2, t) - ht(i, j, k1, t))*scl(zz)*z_trans/scp(psl)
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine wcompgernl

   subroutine wwgradient(cc, st)
!*************************************************
! calculate the gradient of dw/du, dw/dv and dw/dz
! history : january 25 2008, coded by zhongjie he
!*************************************************
! --------------------
      implicit none
! --------------------
      integer           :: i, j, k, t, i1, i2, j1, j2, k1, k2, uu, vv, zz, m, n, s, o, no
      real              :: z1, z2, z3, z4, sz
      integer           :: np(maxdims)
      real              :: ht(numgrid(1), numgrid(2), numgrid(3), numgrid(4))
      real, intent(in)  :: cc(ngptobs, nallobs)
      integer, intent(in):: st
     real              :: sg(numgrid(1), numgrid(2)), cg(numgrid(1), numgrid(2)), zp(numgrid(1), numgrid(2), numgrid(3), numgrid(4))
      real              :: xscle, yscle, xcoefc, ycoefc, xcoefs, ycoefs
! ---------------------
! declare :
!          'cc' is the coefficent used to calculate gw
!          'ht' is the temporary array to save height of each grid point
!          'np' is the index of the grid that involves the observation
!          'sg' and 'cg' are the cos and sin values of deg respectively
!          'uu', 'vv' and 'zz' are indexes of the contral variables for u, v and height respectively
!          'st' is the index of the observation state
! ---------------------

      uu = u_cmpnnt
      vv = v_cmpnnt
      zz = pressure
      if (nallobs .eq. 0) return

      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then           ! for sigma and height coordinate
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            ht(i, j, k, t) = zzz(i, j, k, t)
            zp(i, j, k, t) = 1.0
         end do
         end do
         end do
         end do
         sz = scp(psl)
      else                                              ! for pressure coordinate
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            ht(i, j, k, t) = ppp(k)
         end do
         end do
         end do
         end do
         sz = scp(psl)

!   the following is coded for the case of omiga
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            k1 = max(k - 1, 1)
            k2 = min(k + 1, numgrid(3))
            zp(i, j, k, t) = (grdanals(i, j, k2, t, zz) - grdanals(i, j, k1, t, zz) &
                              + grdbkgnd(i, j, k2, t, zz) - grdbkgnd(i, j, k1, t, zz)) &
                             /(ht(i, j, k2, t) - ht(i, j, k1, t))*scl(zz)*z_trans/scp(psl)
         end do
         end do
         end do
         end do

      end if

      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         sg(i, j) = 0.0
         cg(i, j) = 1.0
      end do
      end do

      xscle = sz/scp(xsl)/2.0
      yscle = sz/scp(ysl)/2.0

! this do loop is to calculate the gradient of vertical velocity to the control variable of u and v, the horizontal velocity.
      o = 0
      do s = 1, st - 1
         o = o + nobstat(s)
      end do
      s = st
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            xcoefc = cc(m, o)/(xxx(i2, j) - xxx(i1, j))*xscle*cg(i, j)*zp(i, j, k, t)
            ycoefc = cc(m, o)/(yyy(i, j2) - yyy(i, j1))*yscle*cg(i, j)*zp(i, j, k, t)
            xcoefs = cc(m, o)/(xxx(i2, j) - xxx(i1, j))*xscle*sg(i, j)*zp(i, j, k, t)
            ycoefs = cc(m, o)/(yyy(i, j2) - yyy(i, j1))*yscle*sg(i, j)*zp(i, j, k, t)
            do k1 = 2, k
               z2 = ht(i, j, k1, t)
               z1 = ht(i, j, k1 - 1, t)

               gradint(i2, j, k1 - 1, t, uu) = gradint(i2, j, k1 - 1, t, uu) - (z2 - z1)*xcoefc
               gradint(i1, j, k1 - 1, t, uu) = gradint(i1, j, k1 - 1, t, uu) + (z2 - z1)*xcoefc
               gradint(i, j2, k1 - 1, t, vv) = gradint(i, j2, k1 - 1, t, vv) - (z2 - z1)*ycoefc
               gradint(i, j1, k1 - 1, t, vv) = gradint(i, j1, k1 - 1, t, vv) + (z2 - z1)*ycoefc
               gradint(i2, j, k1, t, uu) = gradint(i2, j, k1, t, uu) - (z2 - z1)*xcoefc
               gradint(i1, j, k1, t, uu) = gradint(i1, j, k1, t, uu) + (z2 - z1)*xcoefc
               gradint(i, j2, k1, t, vv) = gradint(i, j2, k1, t, vv) - (z2 - z1)*ycoefc
               gradint(i, j1, k1, t, vv) = gradint(i, j1, k1, t, vv) + (z2 - z1)*ycoefc
               gradint(i, j2, k1 - 1, t, uu) = gradint(i, j2, k1 - 1, t, uu) + (z2 - z1)*ycoefs
               gradint(i, j1, k1 - 1, t, uu) = gradint(i, j1, k1 - 1, t, uu) - (z2 - z1)*ycoefs
               gradint(i2, j, k1 - 1, t, vv) = gradint(i2, j, k1 - 1, t, vv) - (z2 - z1)*xcoefs
               gradint(i1, j, k1 - 1, t, vv) = gradint(i1, j, k1 - 1, t, vv) + (z2 - z1)*xcoefs
               gradint(i, j2, k1, t, uu) = gradint(i, j2, k1, t, uu) + (z2 - z1)*ycoefs
               gradint(i, j1, k1, t, uu) = gradint(i, j1, k1, t, uu) - (z2 - z1)*ycoefs
               gradint(i2, j, k1, t, vv) = gradint(i2, j, k1, t, vv) - (z2 - z1)*xcoefs
               gradint(i1, j, k1, t, vv) = gradint(i1, j, k1, t, vv) + (z2 - z1)*xcoefs
            end do

            if (ifpcdnt .eq. 0) then
            do k1 = 2, k
               z1 = ht(i1, j, k1, t)
               z2 = ht(i2, j, k1, t)
               z3 = ht(i1, j, k1 - 1, t)
               z4 = ht(i2, j, k1 - 1, t)
               gradint(i, j, k1, t, uu) = gradint(i, j, k1, t, uu) - (z2 - z1 + z4 - z3)*xcoefc
               gradint(i, j, k1 - 1, t, uu) = gradint(i, j, k1 - 1, t, uu) + (z2 - z1 + z4 - z3)*xcoefc
               gradint(i, j, k1, t, vv) = gradint(i, j, k1, t, vv) - (z2 - z1 + z4 - z3)*xcoefs
               gradint(i, j, k1 - 1, t, vv) = gradint(i, j, k1 - 1, t, vv) + (z2 - z1 + z4 - z3)*xcoefs
            end do
            do k1 = 2, k
               z1 = ht(i, j1, k1, t)
               z2 = ht(i, j2, k1, t)
               z3 = ht(i, j1, k1 - 1, t)
               z4 = ht(i, j2, k1 - 1, t)
               gradint(i, j, k1, t, vv) = gradint(i, j, k1, t, vv) - (z2 - z1 + z4 - z3)*ycoefc
               gradint(i, j, k1 - 1, t, vv) = gradint(i, j, k1 - 1, t, vv) + (z2 - z1 + z4 - z3)*ycoefc
               gradint(i, j, k1, t, uu) = gradint(i, j, k1, t, uu) + (z2 - z1 + z4 - z3)*ycoefs
               gradint(i, j, k1 - 1, t, uu) = gradint(i, j, k1 - 1, t, uu) - (z2 - z1 + z4 - z3)*ycoefs
            end do
            end if

         end do
         end do
         end do
         end do
      end do

! the following statement is to calculate the gradient of vertical velocity to the variable of height, this is only necessary for the case of presure vertical coordinate.

      if (ifpcdnt .eq. 1) then
         o = 0
         do s = 1, st - 1
            o = o + nobstat(s)
         end do
         s = st
         do no = 1, nobstat(s)
            o = o + 1
            do n = 1, maxdims
               np(n) = obsidxpc(n, o)
            end do
            m = 0
            do t = np(4), min0(np(4) + 1, numgrid(4))
            do k = np(3), min0(np(3) + 1, numgrid(3))
            do j = np(2), min0(np(2) + 1, numgrid(2))
            do i = np(1), min0(np(1) + 1, numgrid(1))
               m = m + 1
               i1 = max0(i - 1, 1)
               i2 = min0(i + 1, numgrid(1))
               j1 = max0(j - 1, 1)
               j2 = min0(j + 1, numgrid(2))

               k1 = max(k - 1, 1)
               k2 = min(k + 1, numgrid(3))
   gradint(i, j, k1, t, zz) = gradint(i, j, k1, t, zz) + www(i, j, k, t)/zp(i, j, k, t)*(-1.0/(ppp(k2) - ppp(k1))/scp(psl))*cc(m, o)
    gradint(i, j, k2, t, zz) = gradint(i, j, k2, t, zz) + www(i, j, k, t)/zp(i, j, k, t)*(1.0/(ppp(k2) - ppp(k1))/scp(psl))*cc(m, o)

            end do
            end do
            end do
            end do
         end do
      end if

   end subroutine wwgradient

end module wcompt_gradt
