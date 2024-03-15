module gsbcost_grad
!*************************************************
! calculate cost function and the relevent gradients of gs-balance penalty term, and the penalty term of continual condition.
! history: modified and departed from stmas4d_core module by zhongjie he, january 2008.
!*************************************************

   use prmtrs_stmas
   use generaltools, only: gderivelb, gderiveit, gderiverb

  public  gsblncost_non_uniform_z, gsblngrad_non_uniform_z, gsblncost_p_hs, gsblngrad_p_hs, cntnscost, cntnsgrad, gsblncost_p, gsblngrad_p

!***************************************************
!!comment:
!   this module is used by the module of stmas4d_core to calculate the cost function and the relevent gradients of penalty terms.
!   subroutines:
!      gsblncost_non_uniform_z : calculate the costfunction of gs-balance for the sigma and height coordinate system.
!      gsblngrad_non_uniform_z : calculate gradients of costfunctions correspond to gsblncost_non_uniform_z.
!      gsblncost_p_hs          : calculate the costfunction for the case of pressure coordinate, the hydrostatic approximation is used.
!      gsblngrad_p_hs          : calculate gradients of costfunctions correspond to gsblncost_p_hs.
!      cntnscost               : calculate the penalty term of continual conditions.
!      cntnsgrad               : calculate the gradient of costfunctions correspond to cntnscost.
!      gsblncost_p             : just like gsblncost_p_hs, not use the hydrostatic approximation.
!      gsblngrad_p             : calculate gradients of costfunctions correspond to gsblncost_p.
contains

   subroutine gsblncost_non_uniform_z
!*************************************************
! geostrophic balance penalty term of cost function
! history: august 2007, coded by wei li.
!          february 2008, modified by zhongjie he
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, i1, i2, j1, j2, k1, k2, k3, uu, vv, pp
      real     :: z1, z2, z3, az, bz, cz
      real     :: gs
! --------------------

      uu = u_cmpnnt
      vv = v_cmpnnt
      pp = pressure
      if (pnlt0pu .lt. 1.0e-10 .and. pnlt0pv .lt. 1.0e-10) return

! calculate the term of dp/dx=den*cor*v
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            k1 = max0(k - 1, 1)
            k1 = min0(k1, numgrid(3) - 2)
            k2 = min0(k1 + 1, numgrid(3) - 1)
            k3 = min0(k2 + 1, numgrid(3))
            z1 = zzz(i, j, k1, t)
            z2 = zzz(i, j, k2, t)
            z3 = zzz(i, j, k3, t)
            if (k .eq. 1) call gderivelb(z1, z2, z3, az, bz, cz)
            if (k .ge. 2 .and. k .le. numgrid(3) - 1) call gderiveit(z1, z2, z3, az, bz, cz)
            if (k .eq. numgrid(3)) call gderiverb(z1, z2, z3, az, bz, cz)

            gs = ((grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                  /(xxx(i2, j) - xxx(i1, j)) &
                  - (grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                     + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                  *(zzz(i2, j, k, t) - zzz(i1, j, k, t))/(xxx(i2, j) - xxx(i1, j)))*scl(pp)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                - ((grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                    /(yyy(i, j2) - yyy(i, j1)) &
                    - (grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                  *(zzz(i, j2, k, t) - zzz(i, j1, k, t))/(yyy(i, j2) - yyy(i, j1)))*scl(pp)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                 - cor(i, j)*den(i, j, k, t)*(grdanals(i, j, k, t, vv) + grdbkgnd(i, j, k, t, vv))

            costfun = costfun + pnlt_pv*gs**2
         end do
         end do
         end do
         end do
      end if

! calculate the term of dp/dy=-den*cor*u
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            k1 = max0(k - 1, 1)
            k1 = min0(k1, numgrid(3) - 2)
            k2 = min0(k1 + 1, numgrid(3) - 1)
            k3 = min0(k2 + 1, numgrid(3))
            z1 = zzz(i, j, k1, t)
            z2 = zzz(i, j, k2, t)
            z3 = zzz(i, j, k3, t)
            if (k .eq. 1) call gderivelb(z1, z2, z3, az, bz, cz)
            if (k .ge. 2 .and. k .le. numgrid(3) - 1) call gderiveit(z1, z2, z3, az, bz, cz)
            if (k .eq. numgrid(3)) call gderiverb(z1, z2, z3, az, bz, cz)

            gs = ((grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                  /(yyy(i, j2) - yyy(i, j1)) &
                  - (grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                     + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                  *(zzz(i, j2, k, t) - zzz(i, j1, k, t))/(yyy(i, j2) - yyy(i, j1)))*scl(pp)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                + ((grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                    /(xxx(i2, j) - xxx(i1, j)) &
                    - (grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                  *(zzz(i2, j, k, t) - zzz(i1, j, k, t))/(xxx(i2, j) - xxx(i1, j)))*scl(pp)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                 + cor(i, j)*den(i, j, k, t)*(grdanals(i, j, k, t, uu) + grdbkgnd(i, j, k, t, uu))

            costfun = costfun + pnlt_pu*gs**2
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine gsblncost_non_uniform_z

   subroutine gsblngrad_non_uniform_z
!*************************************************
! geostrophic balance penalty term of gradient
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, i1, i2, j1, j2, k1, k2, k3, uu, vv, pp
      real     :: z1, z2, z3, az, bz, cz
      real     :: gs
! --------------------

      uu = u_cmpnnt
      vv = v_cmpnnt
      pp = pressure
      if (pnlt0pu .lt. 1.0e-10 .and. pnlt0pv .lt. 1.0e-10) return
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            k1 = max0(k - 1, 1)
            k1 = min0(k1, numgrid(3) - 2)
            k2 = min0(k1 + 1, numgrid(3) - 1)
            k3 = min0(k2 + 1, numgrid(3))
            z1 = zzz(i, j, k1, t)
            z2 = zzz(i, j, k2, t)
            z3 = zzz(i, j, k3, t)
            if (k .eq. 1) call gderivelb(z1, z2, z3, az, bz, cz)
            if (k .ge. 2 .and. k .le. numgrid(3) - 1) call gderiveit(z1, z2, z3, az, bz, cz)
            if (k .eq. numgrid(3)) call gderiverb(z1, z2, z3, az, bz, cz)

            gs = ((grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                  /(xxx(i2, j) - xxx(i1, j)) &
                  - (grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                     + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                  *(zzz(i2, j, k, t) - zzz(i1, j, k, t))/(xxx(i2, j) - xxx(i1, j)))*scl(pp)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                - ((grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                    /(yyy(i, j2) - yyy(i, j1)) &
                    - (grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                  *(zzz(i, j2, k, t) - zzz(i, j1, k, t))/(yyy(i, j2) - yyy(i, j1)))*scl(pp)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                 - cor(i, j)*den(i, j, k, t)*(grdanals(i, j, k, t, vv) + grdbkgnd(i, j, k, t, vv))
!
            gradint(i, j, k, t, vv) = gradint(i, j, k, t, vv) + 2.0*pnlt_pv*gs*cor(i, j)*den(i, j, k, t)*(-1.0)
!
            gradint(i2, j, k, t, pp) = gradint(i2, j, k, t, pp) + 2.0*pnlt_pv*gs*scl(pp)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                                       /(xxx(i2, j) - xxx(i1, j))
!
            gradint(i1, j, k, t, pp) = gradint(i1, j, k, t, pp) + 2.0*pnlt_pv*gs*scl(pp)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                                       /(xxx(i2, j) - xxx(i1, j))*(-1.0)
!
            gradint(i, j, k1, t, pp) = gradint(i, j, k1, t, pp) + 2.0*pnlt_pv*gs*scl(pp)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                                       *az*(-1.0)*(zzz(i2, j, k, t) - zzz(i1, j, k, t))/(xxx(i2, j) - xxx(i1, j))
!
            gradint(i, j, k2, t, pp) = gradint(i, j, k2, t, pp) + 2.0*pnlt_pv*gs*scl(pp)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                                       *bz*(-1.0)*(zzz(i2, j, k, t) - zzz(i1, j, k, t))/(xxx(i2, j) - xxx(i1, j))
!
            gradint(i, j, k3, t, pp) = gradint(i, j, k3, t, pp) + 2.0*pnlt_pv*gs*scl(pp)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                                       *cz*(-1.0)*(zzz(i2, j, k, t) - zzz(i1, j, k, t))/(xxx(i2, j) - xxx(i1, j))
!
            gradint(i, j2, k, t, pp) = gradint(i, j2, k, t, pp) - 2.0*pnlt_pu*gs*scl(pp)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                                       /(yyy(i, j2) - yyy(i, j1))
!
            gradint(i, j1, k, t, pp) = gradint(i, j1, k, t, pp) - 2.0*pnlt_pu*gs*scl(pp)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                                       /(yyy(i, j2) - yyy(i, j1))*(-1.0)
!
            gradint(i, j, k1, t, pp) = gradint(i, j, k1, t, pp) - 2.0*pnlt_pu*gs*scl(pp)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                                       *az*(-1.0)*(zzz(i, j2, k, t) - zzz(i, j1, k, t))/(yyy(i, j2) - yyy(i, j1))
!
            gradint(i, j, k2, t, pp) = gradint(i, j, k2, t, pp) - 2.0*pnlt_pu*gs*scl(pp)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                                       *bz*(-1.0)*(zzz(i, j2, k, t) - zzz(i, j1, k, t))/(yyy(i, j2) - yyy(i, j1))
!
            gradint(i, j, k3, t, pp) = gradint(i, j, k3, t, pp) - 2.0*pnlt_pu*gs*scl(pp)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                                       *cz*(-1.0)*(zzz(i, j2, k, t) - zzz(i, j1, k, t))/(yyy(i, j2) - yyy(i, j1))

         end do
         end do
         end do
         end do
      end if
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            k1 = max0(k - 1, 1)
            k1 = min0(k1, numgrid(3) - 2)
            k2 = min0(k1 + 1, numgrid(3) - 1)
            k3 = min0(k2 + 1, numgrid(3))
            z1 = zzz(i, j, k1, t)
            z2 = zzz(i, j, k2, t)
            z3 = zzz(i, j, k3, t)
            if (k .eq. 1) call gderivelb(z1, z2, z3, az, bz, cz)
            if (k .ge. 2 .and. k .le. numgrid(3) - 1) call gderiveit(z1, z2, z3, az, bz, cz)
            if (k .eq. numgrid(3)) call gderiverb(z1, z2, z3, az, bz, cz)
            gs = ((grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                  /(yyy(i, j2) - yyy(i, j1)) &
                  - (grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                     + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                  *(zzz(i, j2, k, t) - zzz(i, j1, k, t))/(yyy(i, j2) - yyy(i, j1)))*scl(pp)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                + ((grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                    /(xxx(i2, j) - xxx(i1, j)) &
                    - (grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                  *(zzz(i2, j, k, t) - zzz(i1, j, k, t))/(xxx(i2, j) - xxx(i1, j)))*scl(pp)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                 + cor(i, j)*den(i, j, k, t)*(grdanals(i, j, k, t, uu) + grdbkgnd(i, j, k, t, uu))
!
            gradint(i, j, k, t, uu) = gradint(i, j, k, t, uu) + 2.0*pnlt_pu*gs*cor(i, j)*den(i, j, k, t)
!
            gradint(i, j2, k, t, pp) = gradint(i, j2, k, t, pp) + 2.0*pnlt_pu*gs*scl(pp)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                                       /(yyy(i, j2) - yyy(i, j1))
!
            gradint(i, j1, k, t, pp) = gradint(i, j1, k, t, pp) + 2.0*pnlt_pu*gs*scl(pp)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                                       /(yyy(i, j2) - yyy(i, j1))*(-1.0)
!
            gradint(i, j, k1, t, pp) = gradint(i, j, k1, t, pp) + 2.0*pnlt_pu*gs*scl(pp)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                                       *az*(-1.0)*(zzz(i, j2, k, t) - zzz(i, j1, k, t))/(yyy(i, j2) - yyy(i, j1))
!
            gradint(i, j, k2, t, pp) = gradint(i, j, k2, t, pp) + 2.0*pnlt_pu*gs*scl(pp)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                                       *bz*(-1.0)*(zzz(i, j2, k, t) - zzz(i, j1, k, t))/(yyy(i, j2) - yyy(i, j1))
!
            gradint(i, j, k3, t, pp) = gradint(i, j, k3, t, pp) + 2.0*pnlt_pu*gs*scl(pp)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                                       *cz*(-1.0)*(zzz(i, j2, k, t) - zzz(i, j1, k, t))/(yyy(i, j2) - yyy(i, j1))
!
            gradint(i2, j, k, t, pp) = gradint(i2, j, k, t, pp) + 2.0*pnlt_pv*gs*scl(pp)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                                       /(xxx(i2, j) - xxx(i1, j))
!
            gradint(i1, j, k, t, pp) = gradint(i1, j, k, t, pp) + 2.0*pnlt_pv*gs*scl(pp)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                                       /(xxx(i2, j) - xxx(i1, j))*(-1.0)
!
            gradint(i, j, k1, t, pp) = gradint(i, j, k1, t, pp) + 2.0*pnlt_pv*gs*scl(pp)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                                       *az*(-1.0)*(zzz(i2, j, k, t) - zzz(i1, j, k, t))/(xxx(i2, j) - xxx(i1, j))
!
            gradint(i, j, k2, t, pp) = gradint(i, j, k2, t, pp) + 2.0*pnlt_pv*gs*scl(pp)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                                       *bz*(-1.0)*(zzz(i2, j, k, t) - zzz(i1, j, k, t))/(xxx(i2, j) - xxx(i1, j))
!
            gradint(i, j, k3, t, pp) = gradint(i, j, k3, t, pp) + 2.0*pnlt_pv*gs*scl(pp)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                                       *cz*(-1.0)*(zzz(i2, j, k, t) - zzz(i1, j, k, t))/(xxx(i2, j) - xxx(i1, j))
!
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine gsblngrad_non_uniform_z

   subroutine gsblncost_p_hs
!*************************************************
! geostrophic balance penalty term of cost function
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      real, parameter :: ga = 9.8d0
      integer  :: i, j, k, t, i1, i2, j1, j2, uu, vv, pp
      real     :: gs
! --------------------

      uu = u_cmpnnt
      vv = v_cmpnnt
      pp = pressure
      if (pnlt0pu .lt. 1.0e-10 .and. pnlt0pv .lt. 1.0e-10) return

! calculate the term of dp/dx=den*cor*v
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            gs = (grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                 /(xxx(i2, j) - xxx(i1, j))*(scl(pp)*z_trans)/scp(xsl)/scp(csl)/scl(vv)*ga*1.0 &
                 - (grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                 /(yyy(i, j2) - yyy(i, j1))*(scl(pp)*z_trans)/scp(ysl)/scp(csl)/scl(vv)*ga*0.0 &
                 - cor(i, j)*(grdanals(i, j, k, t, vv) + grdbkgnd(i, j, k, t, vv))
            costfun = costfun + pnlt_pv*gs**2
         end do
         end do
         end do
         end do
      end if

! calculate the term of dp/dy=-den*cor*u
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            gs = (grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                 /(yyy(i, j2) - yyy(i, j1))*(scl(pp)*z_trans)/scp(ysl)/scp(csl)/scl(uu)*ga*1.0 &
                 + (grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                 /(xxx(i2, j) - xxx(i1, j))*(scl(pp)*z_trans)/scp(xsl)/scp(csl)/scl(uu)*ga*0.0 &
                 + cor(i, j)*(grdanals(i, j, k, t, uu) + grdbkgnd(i, j, k, t, uu))
            costfun = costfun + pnlt_pu*gs**2
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine gsblncost_p_hs

   subroutine gsblngrad_p_hs
!*************************************************
! geostrophic balance penalty term of gradient
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      real, parameter :: ga = 9.8d0
      integer  :: i, j, k, t, i1, i2, j1, j2, uu, vv, pp
      real     :: gs
! --------------------

      uu = u_cmpnnt
      vv = v_cmpnnt
      pp = pressure
      if (pnlt0pu .lt. 1.0e-10 .and. pnlt0pv .lt. 1.0e-10) return
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            gs = (grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                 /(xxx(i2, j) - xxx(i1, j))*(scl(pp)*z_trans)/scp(xsl)/scp(csl)/scl(vv)*ga*1.0 &
                 - (grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                 /(yyy(i, j2) - yyy(i, j1))*(scl(pp)*z_trans)/scp(ysl)/scp(csl)/scl(vv)*ga*0.0 &
                 - cor(i, j)*(grdanals(i, j, k, t, vv) + grdbkgnd(i, j, k, t, vv))
!
            gradint(i, j, k, t, vv) = gradint(i, j, k, t, vv) + 2.0*pnlt_pv*gs*cor(i, j)*(-1.0)
!
            gradint(i2, j, k, t, pp) = gradint(i2, j, k, t, pp) + 2.0*pnlt_pv*gs*(scl(pp)*z_trans)/scp(xsl)/scp(csl)/scl(vv)*1.0 &
                                       /(xxx(i2, j) - xxx(i1, j))*ga
!
       gradint(i1, j, k, t, pp) = gradint(i1, j, k, t, pp) + 2.0*pnlt_pv*gs*(scl(pp)*z_trans)/scp(xsl)/scp(csl)/scl(vv)*(-1.0)*1.0 &
                                       /(xxx(i2, j) - xxx(i1, j))*ga
!
            gradint(i, j2, k, t, pp) = gradint(i, j2, k, t, pp) - 2.0*pnlt_pu*gs*(scl(pp)*z_trans)/scp(ysl)/scp(csl)/scl(vv)*0.0 &
                                       /(yyy(i, j2) - yyy(i, j1))*ga
!
       gradint(i, j1, k, t, pp) = gradint(i, j1, k, t, pp) - 2.0*pnlt_pu*gs*(scl(pp)*z_trans)/scp(ysl)/scp(csl)/scl(vv)*(-1.0)*0.0 &
                                       /(yyy(i, j2) - yyy(i, j1))*ga
!
         end do
         end do
         end do
         end do
      end if
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            gs = (grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                 /(yyy(i, j2) - yyy(i, j1))*(scl(pp)*z_trans)/scp(ysl)/scp(csl)/scl(uu)*ga*1.0 &
                 + (grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                 /(xxx(i2, j) - xxx(i1, j))*(scl(pp)*z_trans)/scp(xsl)/scp(csl)/scl(uu)*ga*0.0 &
                 + cor(i, j)*(grdanals(i, j, k, t, uu) + grdbkgnd(i, j, k, t, uu))
!
            gradint(i, j, k, t, uu) = gradint(i, j, k, t, uu) + 2.0*pnlt_pu*gs*cor(i, j)
!
            gradint(i, j2, k, t, pp) = gradint(i, j2, k, t, pp) + 2.0*pnlt_pu*gs*(scl(pp)*z_trans)/scp(ysl)/scp(csl)/scl(uu)*1.0 &
                                       /(yyy(i, j2) - yyy(i, j1))*ga
!
       gradint(i, j1, k, t, pp) = gradint(i, j1, k, t, pp) + 2.0*pnlt_pu*gs*(scl(pp)*z_trans)/scp(ysl)/scp(csl)/scl(uu)*(-1.0)*1.0 &
                                       /(yyy(i, j2) - yyy(i, j1))*ga
!
            gradint(i2, j, k, t, pp) = gradint(i2, j, k, t, pp) + 2.0*pnlt_pv*gs*(scl(pp)*z_trans)/scp(xsl)/scp(csl)/scl(uu)*0.0 &
                                       /(xxx(i2, j) - xxx(i1, j))*ga
!
       gradint(i1, j, k, t, pp) = gradint(i1, j, k, t, pp) + 2.0*pnlt_pv*gs*(scl(pp)*z_trans)/scp(xsl)/scp(csl)/scl(uu)*(-1.0)*0.0 &
                                       /(xxx(i2, j) - xxx(i1, j))*ga
!
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine gsblngrad_p_hs

   subroutine gsblncost_p
!*************************************************
! geostrophic balance penalty term of cost function
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, i1, i2, j1, j2, k1, k2, k3, uu, vv, pp
      real     :: z1, z2, z3, az, bz, cz
      real     :: gs
! --------------------

      uu = u_cmpnnt
      vv = v_cmpnnt
      pp = pressure
      if (pnlt0pu .lt. 1.0e-10 .and. pnlt0pv .lt. 1.0e-10) return
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            k1 = max0(k - 1, 1)
            k1 = min0(k1, numgrid(3) - 2)
            k2 = min0(k1 + 1, numgrid(3) - 1)
            k3 = min0(k2 + 1, numgrid(3))
            z1 = ppp(k1)
            z2 = ppp(k2)
            z3 = ppp(k3)
            if (k .eq. 1) call gderivelb(z1, z2, z3, az, bz, cz)
            if (k .ge. 2 .and. k .le. numgrid(3) - 1) call gderiveit(z1, z2, z3, az, bz, cz)
            if (k .eq. numgrid(3)) call gderiverb(z1, z2, z3, az, bz, cz)
            gs = (grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                 /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                   + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                 /(xxx(i2, j) - xxx(i1, j))*scp(psl)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                 - (grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                 /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                   + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                 /(yyy(i, j2) - yyy(i, j1))*scp(psl)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                 + cor(i, j)*den(i, j, k, t)*(grdanals(i, j, k, t, vv) + grdbkgnd(i, j, k, t, vv))
            costfun = costfun + pnlt_pv*gs**2
         end do
         end do
         end do
         end do
      end if
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            k1 = max0(k - 1, 1)
            k1 = min0(k1, numgrid(3) - 2)
            k2 = min0(k1 + 1, numgrid(3) - 1)
            k3 = min0(k2 + 1, numgrid(3))
            z1 = ppp(k1)
            z2 = ppp(k2)
            z3 = ppp(k3)
            if (k .eq. 1) call gderivelb(z1, z2, z3, az, bz, cz)
            if (k .ge. 2 .and. k .le. numgrid(3) - 1) call gderiveit(z1, z2, z3, az, bz, cz)
            if (k .eq. numgrid(3)) call gderiverb(z1, z2, z3, az, bz, cz)
            gs = (grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                 /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                   + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                 /(yyy(i, j2) - yyy(i, j1))*scp(psl)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                 + (grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                 /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                   + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                 /(xxx(i2, j) - xxx(i1, j))*scp(psl)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                 - cor(i, j)*den(i, j, k, t)*(grdanals(i, j, k, t, uu) + grdbkgnd(i, j, k, t, uu))
            costfun = costfun + pnlt_pu*gs**2
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine gsblncost_p

   subroutine gsblngrad_p
!*************************************************
! geostrophic balance penalty term of gradient
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, i1, i2, j1, j2, k1, k2, k3, uu, vv, pp
      real     :: z1, z2, z3, az, bz, cz
      real     :: gs
! --------------------

      uu = u_cmpnnt
      vv = v_cmpnnt
      pp = pressure
      if (pnlt0pu .lt. 1.0e-10 .and. pnlt0pv .lt. 1.0e-10) return
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            k1 = max0(k - 1, 1)
            k1 = min0(k1, numgrid(3) - 2)
            k2 = min0(k1 + 1, numgrid(3) - 1)
            k3 = min0(k2 + 1, numgrid(3))
            z1 = ppp(k1)
            z2 = ppp(k2)
            z3 = ppp(k3)
            if (k .eq. 1) call gderivelb(z1, z2, z3, az, bz, cz)
            if (k .ge. 2 .and. k .le. numgrid(3) - 1) call gderiveit(z1, z2, z3, az, bz, cz)
            if (k .eq. numgrid(3)) call gderiverb(z1, z2, z3, az, bz, cz)
            gs = (grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                 /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                   + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                 /(xxx(i2, j) - xxx(i1, j))*scp(psl)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                 - (grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                 /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                   + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                 /(yyy(i, j2) - yyy(i, j1))*scp(psl)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                 + cor(i, j)*den(i, j, k, t)*(grdanals(i, j, k, t, vv) + grdbkgnd(i, j, k, t, vv))
!
            gradint(i, j, k, t, vv) = gradint(i, j, k, t, vv) + 2.0*pnlt_pv*gs*cor(i, j)*den(i, j, k, t)
!
            gradint(i2, j, k, t, pp) = gradint(i2, j, k, t, pp) + 2.0*pnlt_pv*gs*scp(psl)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                                       /(xxx(i2, j) - xxx(i1, j)) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)
!
       gradint(i1, j, k, t, pp) = gradint(i1, j, k, t, pp) + 2.0*pnlt_pv*gs*scp(psl)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*(-1.0)*1.0 &
                                       /(xxx(i2, j) - xxx(i1, j)) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)
!
            gradint(i, j, k1, t, pp) = gradint(i, j, k1, t, pp) + 2.0*pnlt_pv*gs*scp(psl)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                                       *(grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) &
                                      + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp))/(xxx(i2, j) - xxx(i1, j))*az*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
            gradint(i, j, k2, t, pp) = gradint(i, j, k2, t, pp) + 2.0*pnlt_pv*gs*scp(psl)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                                       *(grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) &
                                      + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp))/(xxx(i2, j) - xxx(i1, j))*bz*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
            gradint(i, j, k3, t, pp) = gradint(i, j, k3, t, pp) + 2.0*pnlt_pv*gs*scp(psl)/scp(xsl)/scp(csl)/scl(vv)/scp(dsl)*1.0 &
                                       *(grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) &
                                      + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp))/(xxx(i2, j) - xxx(i1, j))*cz*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
            gradint(i, j2, k, t, pp) = gradint(i, j2, k, t, pp) - 2.0*pnlt_pu*gs*scp(psl)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                                       /(yyy(i, j2) - yyy(i, j1)) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)
!
       gradint(i, j1, k, t, pp) = gradint(i, j1, k, t, pp) - 2.0*pnlt_pu*gs*scp(psl)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*(-1.0)*0.0 &
                                       /(yyy(i, j2) - yyy(i, j1)) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)
!
            gradint(i, j, k1, t, pp) = gradint(i, j, k1, t, pp) - 2.0*pnlt_pu*gs*scp(psl)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                                       *(grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) &
                                      + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp))/(yyy(i, j2) - yyy(i, j1))*az*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
            gradint(i, j, k2, t, pp) = gradint(i, j, k2, t, pp) - 2.0*pnlt_pu*gs*scp(psl)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                                       *(grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) &
                                      + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp))/(yyy(i, j2) - yyy(i, j1))*bz*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
            gradint(i, j, k3, t, pp) = gradint(i, j, k3, t, pp) - 2.0*pnlt_pu*gs*scp(psl)/scp(ysl)/scp(csl)/scl(vv)/scp(dsl)*0.0 &
                                       *(grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) &
                                      + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp))/(yyy(i, j2) - yyy(i, j1))*cz*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
         end do
         end do
         end do
         end do
      end if
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            k1 = max0(k - 1, 1)
            k1 = min0(k1, numgrid(3) - 2)
            k2 = min0(k1 + 1, numgrid(3) - 1)
            k3 = min0(k2 + 1, numgrid(3))
            z1 = ppp(k1)
            z2 = ppp(k2)
            z3 = ppp(k3)
            if (k .eq. 1) call gderivelb(z1, z2, z3, az, bz, cz)
            if (k .ge. 2 .and. k .le. numgrid(3) - 1) call gderiveit(z1, z2, z3, az, bz, cz)
            if (k .eq. numgrid(3)) call gderiverb(z1, z2, z3, az, bz, cz)
            gs = (grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp)) &
                 /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                   + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                 /(yyy(i, j2) - yyy(i, j1))*scp(psl)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                 + (grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp)) &
                 /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                   + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz) &
                 /(xxx(i2, j) - xxx(i1, j))*scp(psl)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                 - cor(i, j)*den(i, j, k, t)*(grdanals(i, j, k, t, uu) + grdbkgnd(i, j, k, t, uu))
!
            gradint(i, j, k, t, uu) = gradint(i, j, k, t, uu) + 2.0*pnlt_pu*gs*cor(i, j)*den(i, j, k, t)*(-1.0)
!
            gradint(i, j2, k, t, pp) = gradint(i, j2, k, t, pp) + 2.0*pnlt_pu*gs*scp(psl)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                                       /(yyy(i, j2) - yyy(i, j1)) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)
!
       gradint(i, j1, k, t, pp) = gradint(i, j1, k, t, pp) + 2.0*pnlt_pu*gs*scp(psl)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*(-1.0)*1.0 &
                                       /(yyy(i, j2) - yyy(i, j1)) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)
!
            gradint(i, j, k1, t, pp) = gradint(i, j, k1, t, pp) + 2.0*pnlt_pu*gs*scp(psl)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                                       *(grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) &
                                      + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp))/(yyy(i, j2) - yyy(i, j1))*az*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
            gradint(i, j, k2, t, pp) = gradint(i, j, k2, t, pp) + 2.0*pnlt_pu*gs*scp(psl)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                                       *(grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) &
                                      + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp))/(yyy(i, j2) - yyy(i, j1))*bz*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
            gradint(i, j, k3, t, pp) = gradint(i, j, k3, t, pp) + 2.0*pnlt_pu*gs*scp(psl)/scp(ysl)/scp(csl)/scl(uu)/scp(dsl)*1.0 &
                                       *(grdanals(i, j2, k, t, pp) - grdanals(i, j1, k, t, pp) &
                                      + grdbkgnd(i, j2, k, t, pp) - grdbkgnd(i, j1, k, t, pp))/(yyy(i, j2) - yyy(i, j1))*cz*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
            gradint(i2, j, k, t, pp) = gradint(i2, j, k, t, pp) + 2.0*pnlt_pv*gs*scp(psl)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                                       /(xxx(i2, j) - xxx(i1, j)) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)
!
       gradint(i1, j, k, t, pp) = gradint(i1, j, k, t, pp) + 2.0*pnlt_pv*gs*scp(psl)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*(-1.0)*0.0 &
                                       /(xxx(i2, j) - xxx(i1, j)) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                       + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)
!
            gradint(i, j, k1, t, pp) = gradint(i, j, k1, t, pp) + 2.0*pnlt_pv*gs*scp(psl)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                                       *(grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) &
                                      + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp))/(xxx(i2, j) - xxx(i1, j))*az*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
            gradint(i, j, k2, t, pp) = gradint(i, j, k2, t, pp) + 2.0*pnlt_pv*gs*scp(psl)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                                       *(grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) &
                                      + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp))/(xxx(i2, j) - xxx(i1, j))*bz*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
            gradint(i, j, k3, t, pp) = gradint(i, j, k3, t, pp) + 2.0*pnlt_pv*gs*scp(psl)/scp(xsl)/scp(csl)/scl(uu)/scp(dsl)*0.0 &
                                       *(grdanals(i2, j, k, t, pp) - grdanals(i1, j, k, t, pp) &
                                      + grdbkgnd(i2, j, k, t, pp) - grdbkgnd(i1, j, k, t, pp))/(xxx(i2, j) - xxx(i1, j))*cz*(-1.0) &
                                      /(grdanals(i, j, k1, t, pp)*az + grdanals(i, j, k2, t, pp)*bz + grdanals(i, j, k3, t, pp)*cz &
                                    + grdbkgnd(i, j, k1, t, pp)*az + grdbkgnd(i, j, k2, t, pp)*bz + grdbkgnd(i, j, k3, t, pp)*cz)**2
!
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine gsblngrad_p

   subroutine cntnscost
!*************************************************
! continous equation
! history: october 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, i1, i2, j1, j2, k1, k2, k3, uu, vv, ww
      real     :: z1, z2, z3, az, bz, cz
      real     :: ce
! --------------------

      uu = u_cmpnnt
      vv = v_cmpnnt
      ww = w_cmpnnt
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            k1 = max0(k - 1, 1)
            k1 = min0(k1, numgrid(3) - 2)
            k2 = min0(k1 + 1, numgrid(3) - 1)
            k3 = min0(k2 + 1, numgrid(3))
            z1 = ppp(k1)
            z2 = ppp(k2)
            z3 = ppp(k3)
            if (k .eq. 1) call gderivelb(z1, z2, z3, az, bz, cz)
            if (k .ge. 2 .and. k .le. numgrid(3) - 1) call gderiveit(z1, z2, z3, az, bz, cz)
            if (k .eq. numgrid(3)) call gderiverb(z1, z2, z3, az, bz, cz)
            ce = (grdanals(i2, j, k, t, uu) - grdanals(i1, j, k, t, uu) + grdbkgnd(i2, j, k, t, uu) - grdbkgnd(i1, j, k, t, uu)) &
                 /(xxx(i2, j) - xxx(i1, j))*scl(uu)/scp(xsl)*scp(psl)/scl(ww)*1.0 &
                 - (grdanals(i, j2, k, t, uu) - grdanals(i, j1, k, t, uu) + grdbkgnd(i, j2, k, t, uu) - grdbkgnd(i, j1, k, t, uu)) &
                 /(yyy(i, j2) - yyy(i, j1))*scl(uu)/scp(ysl)*scp(psl)/scl(ww)*0.0 &
                 + (grdanals(i, j2, k, t, vv) - grdanals(i, j1, k, t, vv) + grdbkgnd(i, j2, k, t, vv) - grdbkgnd(i, j1, k, t, vv)) &
                 /(yyy(i, j2) - yyy(i, j1))*scl(vv)/scp(ysl)*scp(psl)/scl(ww)*1.0 &
                 + (grdanals(i2, j, k, t, vv) - grdanals(i1, j, k, t, vv) + grdbkgnd(i2, j, k, t, vv) - grdbkgnd(i1, j, k, t, vv)) &
                 /(xxx(i2, j) - xxx(i1, j))*scl(vv)/scp(xsl)*scp(psl)/scl(ww)*0.0 &
                 + (grdanals(i, j, k1, t, ww)*az + grdanals(i, j, k2, t, ww)*bz + grdanals(i, j, k3, t, ww)*cz &
                    + grdbkgnd(i, j, k1, t, ww)*az + grdbkgnd(i, j, k2, t, ww)*bz + grdbkgnd(i, j, k3, t, ww)*cz)
            costfun = costfun + pnlt0pu*ce**2
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine cntnscost

   subroutine cntnsgrad
!*************************************************
! continous equation
! history: october 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, i1, i2, j1, j2, k1, k2, k3, uu, vv, ww
      real     :: z1, z2, z3, az, bz, cz
      real     :: ce
! --------------------

      uu = u_cmpnnt
      vv = v_cmpnnt
      ww = w_cmpnnt
      if (numgrid(1) .ge. 2 .and. numgrid(2) .ge. 2 .and. numgrid(3) .ge. 3 .and. numstat .ge. 3) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            i1 = max0(i - 1, 1)
            i2 = min0(i + 1, numgrid(1))
            j1 = max0(j - 1, 1)
            j2 = min0(j + 1, numgrid(2))
            k1 = max0(k - 1, 1)
            k1 = min0(k1, numgrid(3) - 2)
            k2 = min0(k1 + 1, numgrid(3) - 1)
            k3 = min0(k2 + 1, numgrid(3))
            z1 = ppp(k1)
            z2 = ppp(k2)
            z3 = ppp(k3)
            if (k .eq. 1) call gderivelb(z1, z2, z3, az, bz, cz)
            if (k .ge. 2 .and. k .le. numgrid(3) - 1) call gderiveit(z1, z2, z3, az, bz, cz)
            if (k .eq. numgrid(3)) call gderiverb(z1, z2, z3, az, bz, cz)
            ce = (grdanals(i2, j, k, t, uu) - grdanals(i1, j, k, t, uu) + grdbkgnd(i2, j, k, t, uu) - grdbkgnd(i1, j, k, t, uu)) &
                 /(xxx(i2, j) - xxx(i1, j))*scl(uu)/scp(xsl)*scp(psl)/scl(ww)*1.0 &
                 - (grdanals(i, j2, k, t, uu) - grdanals(i, j1, k, t, uu) + grdbkgnd(i, j2, k, t, uu) - grdbkgnd(i, j1, k, t, uu)) &
                 /(yyy(i, j2) - yyy(i, j1))*scl(uu)/scp(ysl)*scp(psl)/scl(ww)*0.0 &
                 + (grdanals(i, j2, k, t, vv) - grdanals(i, j1, k, t, vv) + grdbkgnd(i, j2, k, t, vv) - grdbkgnd(i, j1, k, t, vv)) &
                 /(yyy(i, j2) - yyy(i, j1))*scl(vv)/scp(ysl)*scp(psl)/scl(ww)*1.0 &
                 + (grdanals(i2, j, k, t, vv) - grdanals(i1, j, k, t, vv) + grdbkgnd(i2, j, k, t, vv) - grdbkgnd(i1, j, k, t, vv)) &
                 /(xxx(i2, j) - xxx(i1, j))*scl(vv)/scp(xsl)*scp(psl)/scl(ww)*0.0 &
                 + (grdanals(i, j, k1, t, ww)*az + grdanals(i, j, k2, t, ww)*bz + grdanals(i, j, k3, t, ww)*cz &
                    + grdbkgnd(i, j, k1, t, ww)*az + grdbkgnd(i, j, k2, t, ww)*bz + grdbkgnd(i, j, k3, t, ww)*cz)
            gradint(i2, j, k, t, uu) = gradint(i2, j, k, t, uu) + 2.0*pnlt0pu*ce &
                                       /(xxx(i2, j) - xxx(i1, j))*scl(uu)/scp(xsl)*scp(psl)/scl(ww)*1.0
            gradint(i1, j, k, t, uu) = gradint(i1, j, k, t, uu) + 2.0*pnlt0pu*ce*(-1.0) &
                                       /(xxx(i2, j) - xxx(i1, j))*scl(uu)/scp(xsl)*scp(psl)/scl(ww)*1.0
            gradint(i, j2, k, t, uu) = gradint(i, j2, k, t, uu) + 2.0*pnlt0pu*ce*(-1.0) &
                                       /(yyy(i, j2) - yyy(i, j1))*scl(uu)/scp(ysl)*scp(psl)/scl(ww)*0.0
            gradint(i, j1, k, t, uu) = gradint(i, j1, k, t, uu) + 2.0*pnlt0pu*ce &
                                       /(yyy(i, j2) - yyy(i, j1))*scl(uu)/scp(ysl)*scp(psl)/scl(ww)*0.0
            gradint(i, j2, k, t, vv) = gradint(i, j2, k, t, vv) + 2.0*pnlt0pu*ce &
                                       /(yyy(i, j2) - yyy(i, j1))*scl(vv)/scp(ysl)*scp(psl)/scl(ww)*1.0
            gradint(i, j1, k, t, vv) = gradint(i, j1, k, t, vv) + 2.0*pnlt0pu*ce*(-1.0) &
                                       /(yyy(i, j2) - yyy(i, j1))*scl(vv)/scp(ysl)*scp(psl)/scl(ww)*1.0
            gradint(i2, j, k, t, vv) = gradint(i2, j, k, t, vv) + 2.0*pnlt0pu*ce &
                                       /(xxx(i2, j) - xxx(i1, j))*scl(vv)/scp(xsl)*scp(psl)/scl(ww)*0.0
            gradint(i1, j, k, t, vv) = gradint(i1, j, k, t, vv) + 2.0*pnlt0pu*ce*(-1.0) &
                                       /(xxx(i2, j) - xxx(i1, j))*scl(vv)/scp(xsl)*scp(psl)/scl(ww)*0.0
            gradint(i, j, k1, t, ww) = gradint(i, j, k1, t, ww) + 2.0*pnlt0pu*ce*az
            gradint(i, j, k2, t, ww) = gradint(i, j, k2, t, ww) + 2.0*pnlt0pu*ce*bz
            gradint(i, j, k3, t, ww) = gradint(i, j, k3, t, ww) + 2.0*pnlt0pu*ce*cz
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine cntnsgrad

end module gsbcost_grad
