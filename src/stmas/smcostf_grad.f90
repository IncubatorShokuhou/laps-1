module smcostf_grad
! calculate the penalty cost function of smooth and the relertive gradient
! history: january 2008, separated frome the module 'stmas4d_core' by zhongjie he

   use prmtrs_stmas
   use generaltools, only: g2orderit

   public smothcost, smothgrad

!**************************************************
!comment:
!   this module is used by the module of costfun_grad to calculate the penalty cost function of smooth and the relertive gradient
!   subroutines:
!      smothcost: calculate smooth penalty term of cost function.
!      smothgrad: calculate gradients of smooth term.
!**************************************************
contains

   subroutine smothcost
!*************************************************
! smooth penalty term of cost function
! history: august 2007, coded by wei li.
!
!          may 2015 modified by yuanfu xie
!          add a penalty parameter check
!          and laplacian check: skip if small enough
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s
      real     :: z1, z2, z3, az, bz, cz, sm
      double precision :: smx, smy, smz, smt
! --------------------
! penalty term

      print *, 'smoothing x: ', penal_x(1:numstat)
      print *, 'smoothing y: ', penal_y(1:numstat)
      print *, 'smoothing z: ', penal_z(1:numstat)
      print *, 'smoothing t: ', penal_t(1:numstat)
      smx = 0.0d0
      if (numgrid(1) .ge. 3) then
         do s = 1, numstat
            if (penal_x(s) .gt. 0.0) then
            do t = 1, numgrid(4)
            do k = 1, numgrid(3)
            do j = 1, numgrid(2)
            do i = 2, numgrid(1) - 1
               ! costfun=costfun+penal_x(s)* &
               sm = &
                  ((grdanals(i - 1, j, k, t, s) - 2*grdanals(i, j, k, t, s) + grdanals(i + 1, j, k, t, s))/ &
                   (grdspac(1)*grdspac(1)))**2
               ! yuanfu added a check to skip to small laplacian: may 2015
               if (sm .gt. penal_x(s)*0.1) smx = smx + penal_x(s)*sm
            end do
            end do
            end do
            end do
            end if
         end do
      end if

      smy = 0.0d0
      if (numgrid(2) .ge. 3) then
         do s = 1, numstat
            if (penal_y(s) .gt. 0.0) then
            do t = 1, numgrid(4)
            do k = 1, numgrid(3)
            do j = 2, numgrid(2) - 1
            do i = 1, numgrid(1)
               ! costfun=costfun+penal_y(s)* &
               sm = &
                  ((grdanals(i, j - 1, k, t, s) - 2*grdanals(i, j, k, t, s) + grdanals(i, j + 1, k, t, s))/ &
                   (grdspac(2)*grdspac(2)))**2
               ! yuanfu added a check to skip to small laplacian: may 2015
               if (sm .gt. penal_y(s)*0.1) smy = smy + penal_y(s)*sm
            end do
            end do
            end do
            end do
            end if
         end do
      end if

      smz = 0.0d0
      if (numgrid(3) .ge. 3) then
         do s = 1, numstat
            if (penal_z(s) .gt. 0.0) then
            do t = 1, numgrid(4)
            do k = 2, numgrid(3) - 1
            do j = 1, numgrid(2)
            do i = 1, numgrid(1)
               if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then       ! for sigma and height coordinate
                  z1 = zzz(i, j, k - 1, t)
                  z2 = zzz(i, j, k, t)
                  z3 = zzz(i, j, k + 1, t)
               elseif (ifpcdnt .eq. 1) then                     ! for pressure coordinate
                  z1 = ppp(k - 1)
                  z2 = ppp(k)
                  z3 = ppp(k + 1)
               end if

               call g2orderit(z1, z2, z3, az, bz, cz)

               az = az*(z2 - z1)*(z3 - z2)
               bz = bz*(z2 - z1)*(z3 - z2)
               cz = cz*(z2 - z1)*(z3 - z2)
               ! costfun=costfun+penal_z(s)* &
               sm = &
                  (az*grdanals(i, j, k - 1, t, s) + bz*grdanals(i, j, k, t, s) + cz*grdanals(i, j, k + 1, t, s))**2
               ! yuanfu added a check to skip to small laplacian: may 2015
               if (sm .gt. penal_z(s)*0.1) smz = smz + penal_z(s)*sm
            end do
            end do
            end do
            end do
            end if
         end do
      end if

      smt = 0.0d0
      if (numgrid(4) .ge. 3) then
         do s = 1, numstat
            if (penal_t(s) .gt. 0.0) then
            do t = 2, numgrid(4) - 1
            do k = 1, numgrid(3)
            do j = 1, numgrid(2)
            do i = 1, numgrid(1)
               ! costfun=costfun+penal_t(s)* &
               sm = &
                  ((grdanals(i, j, k, t - 1, s) - 2*grdanals(i, j, k, t, s) + grdanals(i, j, k, t + 1, s))/ &
                   (grdspac(4)*grdspac(4)))**2
               ! yuanfu added a check to skip to small laplacian: may 2015
               if (sm .gt. penal_t(s)*0.1) smt = smt + penal_t(s)*sm
            end do
            end do
            end do
            end do
            end if
         end do
      end if

      costfun = costfun + smx + smy + smz + smt

      return
   end subroutine smothcost

   subroutine smothgrad
!*************************************************
! smooth penalty term of gradient
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s
      real     :: z1, z2, z3, az, bz, cz
      real     :: sm
! --------------------
! penalty term
      if (numgrid(1) .ge. 3) then
         do s = 1, numstat
            if (penal_x(s) .gt. 0.0) then
            do t = 1, numgrid(4)
            do k = 1, numgrid(3)
            do j = 1, numgrid(2)
            do i = 2, numgrid(1) - 1
               sm = (grdanals(i - 1, j, k, t, s) - 2.0*grdanals(i, j, k, t, s) + grdanals(i + 1, j, k, t, s))/ &
                    (grdspac(1)*grdspac(1))
               ! yuanfu added a check to skip to small laplacian: may 2015
               if (sm .gt. penal_x(s)*0.1) then
                  gradint(i - 1, j, k, t, s) = gradint(i - 1, j, k, t, s) + 2.0*penal_x(s)/(grdspac(1)*grdspac(1))*sm
                  gradint(i, j, k, t, s) = gradint(i, j, k, t, s) + 2.0*penal_x(s)/(grdspac(1)*grdspac(1))*(-2.0)*sm
                  gradint(i + 1, j, k, t, s) = gradint(i + 1, j, k, t, s) + 2.0*penal_x(s)/(grdspac(1)*grdspac(1))*sm
               end if
            end do
            end do
            end do
            end do
            end if
         end do
      end if
      if (numgrid(2) .ge. 3) then
         do s = 1, numstat
            if (penal_y(s) .gt. 0.0) then
            do t = 1, numgrid(4)
            do k = 1, numgrid(3)
            do i = 1, numgrid(1)
            do j = 2, numgrid(2) - 1
               sm = (grdanals(i, j - 1, k, t, s) - 2.0*grdanals(i, j, k, t, s) + grdanals(i, j + 1, k, t, s))/ &
                    (grdspac(2)*grdspac(2))
               ! yuanfu added a check to skip to small laplacian: may 2015
               if (sm .gt. penal_y(s)*0.1) then
                  gradint(i, j - 1, k, t, s) = gradint(i, j - 1, k, t, s) + 2.0*penal_y(s)/(grdspac(2)*grdspac(2))*sm
                  gradint(i, j, k, t, s) = gradint(i, j, k, t, s) + 2.0*penal_y(s)/(grdspac(2)*grdspac(2))*(-2.0)*sm
                  gradint(i, j + 1, k, t, s) = gradint(i, j + 1, k, t, s) + 2.0*penal_y(s)/(grdspac(2)*grdspac(2))*sm
               end if
            end do
            end do
            end do
            end do
            end if
         end do
      end if
      if (numgrid(3) .ge. 3) then
         do s = 1, numstat
            if (penal_z(s) .gt. 0.0) then
            do t = 1, numgrid(4)
            do j = 1, numgrid(2)
            do i = 1, numgrid(1)
            do k = 2, numgrid(3) - 1
               if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then       ! for sigma and height coordinate
                  z1 = zzz(i, j, k - 1, t)
                  z2 = zzz(i, j, k, t)
                  z3 = zzz(i, j, k + 1, t)
               elseif (ifpcdnt .eq. 1) then                     ! for pressur coordinate
                  z1 = ppp(k - 1)
                  z2 = ppp(k)
                  z3 = ppp(k + 1)
               end if
               call g2orderit(z1, z2, z3, az, bz, cz)
               az = az*(z2 - z1)*(z3 - z2)
               bz = bz*(z2 - z1)*(z3 - z2)
               cz = cz*(z2 - z1)*(z3 - z2)
               sm = az*grdanals(i, j, k - 1, t, s) + bz*grdanals(i, j, k, t, s) + cz*grdanals(i, j, k + 1, t, s)
               ! yuanfu added a check to skip to small laplacian: may 2015
               if (sm .gt. penal_z(s)*0.1) then
                  gradint(i, j, k - 1, t, s) = gradint(i, j, k - 1, t, s) + 2.0*penal_z(s)*az*sm
                  gradint(i, j, k, t, s) = gradint(i, j, k, t, s) + 2.0*penal_z(s)*bz*sm
                  gradint(i, j, k + 1, t, s) = gradint(i, j, k + 1, t, s) + 2.0*penal_z(s)*cz*sm
               end if
            end do
            end do
            end do
            end do
            end if
         end do
      end if
      if (numgrid(4) .ge. 3) then
         do s = 1, numstat
            if (penal_t(s) .gt. 0.0) then
            do k = 1, numgrid(3)
            do j = 1, numgrid(2)
            do i = 1, numgrid(1)
            do t = 2, numgrid(4) - 1
               sm = (grdanals(i, j, k, t - 1, s) - 2.0*grdanals(i, j, k, t, s) + grdanals(i, j, k, t + 1, s))/ &
                    (grdspac(4)*grdspac(4))
               ! yuanfu added a check to skip to small laplacian: may 2015
               if (sm .gt. penal_t(s)*0.1) then
                  gradint(i, j, k, t - 1, s) = gradint(i, j, k, t - 1, s) + 2.0*penal_t(s)/(grdspac(4)*grdspac(4))*sm
                  gradint(i, j, k, t, s) = gradint(i, j, k, t, s) + 2.0*penal_t(s)/(grdspac(4)*grdspac(4))*(-2.0)*sm
                  gradint(i, j, k, t + 1, s) = gradint(i, j, k, t + 1, s) + 2.0*penal_t(s)/(grdspac(4)*grdspac(4))*sm
               end if
            end do
            end do
            end do
            end do
            end if
         end do
      end if
      return
   end subroutine smothgrad

end module smcostf_grad
