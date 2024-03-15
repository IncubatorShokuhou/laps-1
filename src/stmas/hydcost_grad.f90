module hydcost_grad
!*************************************************
! calculate cost function and the relevent gradients
! of hydrostatic condition penalty term.
! history: corded by zhongjie he, may 2008.
!
! modified by yuanfu xie for adding scl to hydrograd
!          and adding hydrocast_xie and hydrograd
!          for dp/dz+pg/(rt)
!*************************************************

   use prmtrs_stmas

   public hydrocost, hydrograd, hydrocost_xie, hydrograd_xie
   public hydrocost_shuyuan, hydrograd_shuyuan

!  real , parameter ::  r=287
!  real , parameter ::  grav=9.8

!***************************************************
!!comment:
!   this module is used by the module of stmas4d_core to calculate the cost function and gradients of hydrostatic condition penalty terms.
!   subroutines:
!      hydrocost: costfunction of hydrostatic condition for pressure coordinate.
!      hydrograd: calculate gradients of hydrocost.
!      hydrograd_xie:
contains

   subroutine hydrocost
!*************************************************
! hydrostatic condition penalty term of cost function, for the case of pressure coordinate.
! history: may 2008, coded by zhongjie he.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, tt, zz
      real     :: hy
      real     :: r, grav
! --------------------

      r = 287
      grav = 9.8
      tt = temprtur
      zz = pressure
      if (pnlt0hy .lt. 1.0e-10) return
      if (numgrid(3) .ge. 2) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3) - 1
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            hy=grav*(grdanals(i,j,k+1,t,zz)-grdanals(i,j,k,t,zz)+grdbkgnd(i,j,k+1,t,zz)-grdbkgnd(i,j,k,t,zz))*scl(zz)  &
       + 0.5*r*(grdanals(i, j, k + 1, t, tt) + grdanals(i, j, k, t, tt) + grdbkgnd(i, j, k + 1, t, tt) + grdbkgnd(i, j, k, t, tt)) &
                *(log(ppp(k + 1)*scp(psl) + orivtcl) - log(ppp(k)*scp(psl) + orivtcl))*scl(tt)
            costfun = costfun + pnlt_hy*hy**2
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine hydrocost

   subroutine hydrocost_xie
!*************************************************
! hydrostatic condition penalty term of cost function, for the case of pressure coordinate.
! history: may 2008, coded by zhongjie he.
!
!          modified by yuanfu xie for penalizing
!          dp/dz+g*p/(rt).
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, tt, zz
      real     :: hy, dp, dz, ps, tm
      real     :: r, grav
! --------------------

      r = 287
      grav = 9.8
      tt = temprtur
      zz = pressure
      if (pnlt0hy .lt. 1.0e-10) return
      if (numgrid(3) .ge. 2) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3) - 1
            dp = (ppp(k + 1) - ppp(k))*scp(psl)
            ps = 0.5*(ppp(k + 1) + ppp(k))*scp(psl) + orivtcl
            do j = 1, numgrid(2)
            do i = 1, numgrid(1)

               dz = (grdanals(i, j, k + 1, t, zz) - grdanals(i, j, k, t, zz) + &
                     grdbkgnd(i, j, k + 1, t, zz) - grdbkgnd(i, j, k, t, zz))*scl(zz)
               tm = 0.5*(grdanals(i, j, k + 1, t, tt) + grdanals(i, j, k, t, tt) + &
                         grdbkgnd(i, j, k + 1, t, tt) + grdbkgnd(i, j, k, t, tt))*scl(tt)
               hy = dp/dz + grav*ps/(r*tm)

               costfun = costfun + pnlt_hy*hy*hy        ! divide weight by yuanfu
            end do
            end do
         end do
         end do
      end if
      return
   end subroutine hydrocost_xie

   subroutine hydrograd
!*************************************************
! gradient of hydrostatic condition, corresponding to subroutine hydrocost.
! history: may 2008, coded by zhongjie he.
!
!          modified by yuanfu xie for correcting
!          gradient computation. the code was missing
!          a multiplication of scls
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, tt, zz
      real     :: hy
      real     :: r, grav
! --------------------
      r = 287
      grav = 9.8

      tt = temprtur
      zz = pressure
      if (pnlt0hy .lt. 1.0e-10) return
      if (numgrid(3) .ge. 2) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3) - 1
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            hy=grav*(grdanals(i,j,k+1,t,zz)-grdanals(i,j,k,t,zz)+grdbkgnd(i,j,k+1,t,zz)-grdbkgnd(i,j,k,t,zz))*scl(zz)  &
       + 0.5*r*(grdanals(i, j, k + 1, t, tt) + grdanals(i, j, k, t, tt) + grdbkgnd(i, j, k + 1, t, tt) + grdbkgnd(i, j, k, t, tt)) &
                *(log(ppp(k + 1)*scp(psl) + orivtcl) - log(ppp(k)*scp(psl) + orivtcl))*scl(tt)
!
            gradint(i, j, k, t, tt) = gradint(i, j, k, t, tt) + 2.0*pnlt_hy*hy*scl(tt) &
                                      *0.5*r*(log(ppp(k + 1)*scp(psl) + orivtcl) - log(ppp(k)*scp(psl) + orivtcl))
!
            gradint(i, j, k + 1, t, tt) = gradint(i, j, k + 1, t, tt) + 2.0*pnlt_hy*hy*scl(tt) &
                                          *0.5*r*(log(ppp(k + 1)*scp(psl) + orivtcl) - log(ppp(k)*scp(psl) + orivtcl))
!
            gradint(i, j, k, t, zz) = gradint(i, j, k, t, zz) + 2.0*pnlt_hy*hy*(-grav)*scl(zz)
!
            gradint(i, j, k + 1, t, zz) = gradint(i, j, k + 1, t, zz) + 2.0*pnlt_hy*hy*grav*scl(zz)
!
         end do
         end do
         end do
         end do
      end if
      return
   end subroutine hydrograd

   subroutine hydrograd_xie
!*************************************************
! gradient of hydrostatic condition, corresponding to subroutine hydrocost.
! history: may 2008, coded by zhongjie he.
!
!          modified by yuanfu xie for penalizing
!          dp/dz+g*p/(rt).
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, tt, zz
      real     :: hy, dp, dz, ps, tm, d1
      real     :: r, grav
! --------------------
      r = 287
      grav = 9.8
      tt = temprtur
      zz = pressure
      if (pnlt0hy .lt. 1.0e-10) return
      if (numgrid(3) .ge. 2) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3) - 1
            dp = (ppp(k + 1) - ppp(k))*scp(psl)
            ps = 0.5*(ppp(k + 1) + ppp(k))*scp(psl) + orivtcl
            do j = 1, numgrid(2)
            do i = 1, numgrid(1)

               dz = (grdanals(i, j, k + 1, t, zz) - grdanals(i, j, k, t, zz) + &
                     grdbkgnd(i, j, k + 1, t, zz) - grdbkgnd(i, j, k, t, zz))*scl(zz)
               tm = 0.5*(grdanals(i, j, k + 1, t, tt) + grdanals(i, j, k, t, tt) + &
                         grdbkgnd(i, j, k + 1, t, tt) + grdbkgnd(i, j, k, t, tt))*scl(tt)
               hy = dp/dz + grav*ps/(r*tm)

               d1 = -grav*ps/(r*tm*tm)
!
               gradint(i, j, k, t, tt) = gradint(i, j, k, t, tt) + 2.0*pnlt_hy*hy*0.5*d1*scl(tt)
!
               gradint(i, j, k + 1, t, tt) = gradint(i, j, k + 1, t, tt) + 2.0*pnlt_hy*hy*0.5*d1*scl(tt)

               d1 = -dp/(dz*dz)
!
               gradint(i, j, k, t, zz) = gradint(i, j, k, t, zz) - 2.0*pnlt_hy*hy*d1*scl(zz)
!
               gradint(i, j, k + 1, t, zz) = gradint(i, j, k + 1, t, zz) + 2.0*pnlt_hy*hy*d1*scl(zz)
!
            end do
            end do
         end do
         end do
      end if

      return

   end subroutine hydrograd_xie
!=================================================
   subroutine hydrocost_shuyuan
!*************************************************
! hydrostatic condition penalty term of cost function, for the case of pressure coordinate.
! history: may 2008, coded by zhongjie he.
!          modified by yuanfu xie for penalizing
!          dp/dz=g*p/(rt).
!          modified by shuyuan liu for using density 01-09-2010
!          dp/dz=g*den,den=denv+denr+dens
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, tt, zz, rr, rs
      real     :: hy, dp, dz, ps, tm
      real     :: r, grav, tempr, temps, tempv
      real     :: den, denv, dens, denr! air density

! --------------------

      r = 287
      grav = 9.8
      tt = temprtur
      zz = pressure
      rr = rour_cmpnnt
      rs = rous_cmpnnt

      if (pnlt0hy .lt. 1.0e-10) return

      if (numgrid(3) .ge. 2) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3) - 1
            dp = (ppp(k + 1) - ppp(k))*scp(psl)
            ps = 0.5*(ppp(k + 1) + ppp(k))*scp(psl) + orivtcl
            do j = 1, numgrid(2)
            do i = 1, numgrid(1)

               tempr = 0.5*(grdanals(i, j, k, t, rr) + grdbkgnd(i, j, k, t, rr) + &
                            grdanals(i, j, k + 1, t, rr) + grdbkgnd(i, j, k + 1, t, rr))*scl(numstat + 1)

               dz = (grdanals(i, j, k + 1, t, zz) - grdanals(i, j, k, t, zz) + &
                     grdbkgnd(i, j, k + 1, t, zz) - grdbkgnd(i, j, k, t, zz))*scl(zz)
               tm = 0.5*(grdanals(i, j, k + 1, t, tt) + grdanals(i, j, k, t, tt) + &
                         grdbkgnd(i, j, k + 1, t, tt) + grdbkgnd(i, j, k, t, tt))*scl(tt)
               hy = dp/dz + grav*(tempr + ps/(r*tm))

               costfun = costfun + 0.5*pnlt_hy*hy*hy        ! divide weight by yuanfu
            end do
            end do
         end do
         end do
      end if
      return
   end subroutine hydrocost_shuyuan
!========================================
!========================================
   subroutine hydrograd_shuyuan
!*************************************************
! gradient of hydrostatic condition, corresponding to subroutine hydrocost.
! history: may 2008, coded by zhongjie he.
!
!          modified by yuanfu xie for penalizing
!          dp/dz+g*p/(rt).
!          modified by shuyuan liu for adding the rain water 2010/09
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, tt, zz, rr, rs
      real     :: hy, dp, dz, ps, tm, d1
      real     :: r, grav, tempr, temps, tempv
      real     :: den, denv, dens, denr! air density
! --------------------
      r = 287
      grav = 9.8
      tt = temprtur
      zz = pressure
      rr = rour_cmpnnt
      rs = rous_cmpnnt

      if (pnlt0hy .lt. 1.0e-10) return

      if (numgrid(3) .ge. 2) then
         do t = 1, numgrid(4)
         do k = 1, numgrid(3) - 1
            dp = (ppp(k + 1) - ppp(k))*scp(psl)
            ps = 0.5*(ppp(k + 1) + ppp(k))*scp(psl) + orivtcl
            do j = 1, numgrid(2)
            do i = 1, numgrid(1)
               tempr = 0.5*(grdanals(i, j, k, t, rr) + grdbkgnd(i, j, k, t, rr) + &
                            grdanals(i, j, k + 1, t, rr) + grdbkgnd(i, j, k + 1, t, rr))*scl(numstat + 1)

               dz = (grdanals(i, j, k + 1, t, zz) - grdanals(i, j, k, t, zz) + &
                     grdbkgnd(i, j, k + 1, t, zz) - grdbkgnd(i, j, k, t, zz))*scl(zz)
               tm = 0.5*(grdanals(i, j, k + 1, t, tt) + grdanals(i, j, k, t, tt) + &
                         grdbkgnd(i, j, k + 1, t, tt) + grdbkgnd(i, j, k, t, tt))*scl(tt)
               hy = dp/dz + grav*(tempr + ps/(r*tm))

               d1 = -grav*ps/(r*tm*tm)  !!tm
!
               gradint(i, j, k, t, tt) = gradint(i, j, k, t, tt) + pnlt_hy*hy*0.5*d1*scl(tt)
!
               gradint(i, j, k + 1, t, tt) = gradint(i, j, k + 1, t, tt) + pnlt_hy*hy*0.5*d1*scl(tt)

               d1 = -dp/(dz*dz)  !!dz
!
               gradint(i, j, k, t, zz) = gradint(i, j, k, t, zz) - pnlt_hy*hy*d1*scl(zz)
!
               gradint(i, j, k + 1, t, zz) = gradint(i, j, k + 1, t, zz) + pnlt_hy*hy*d1*scl(zz)

               d1 = grav   !!tempr
               gradint(i, j, k, t, rr) = gradint(i, j, k, t, rr) + pnlt_hy*hy*0.5*d1*scl(numstat + 1)
!
               gradint(i, j, k + 1, t, rr) = gradint(i, j, k + 1, t, rr) + pnlt_hy*hy*0.5*d1*scl(numstat + 1)
!
            end do
            end do
         end do
         end do
      end if

      return

   end subroutine hydrograd_shuyuan
!============================================

end module hydcost_grad

