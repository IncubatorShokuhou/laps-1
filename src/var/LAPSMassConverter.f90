!dis
!dis    open source license/disclaimer, forecast systems laboratory
!dis    noaa/oar/fsl, 325 broadway boulder, co 80305
!dis
!dis    this software is distributed under the open source definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    in particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - all modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - if significant modifications or enhancements are made to this
!dis    software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    this software and its documentation are in the public domain
!dis    and are furnished "as is."  the authors, the united states
!dis    government, its instrumentalities, officers, employees, and
!dis    agents make no warranty, express or implied, as to the usefulness
!dis    of the software and documentation for any purpose.  they assume
!dis    no responsibility (1) for the use of the software and
!dis    documentation; or (2) to provide technical support to users.
!dis

subroutine laps2mass(vlaps, imax, jmax, kmax, plaps, psurf, eta, &
                     iorder, logp, vmass)

!==========================================================
!  this routine converts the laps grid function into gsi
!  grid over mass coordinate.
!
!  history: jul. 2006 by yuanfu xie.
!==========================================================

   ! this routine converts a laps variable (varlaps) on laps
   ! coordinate into one on a mass coordinate.

   implicit none

   integer, intent(in) :: imax, jmax, kmax          ! 3d array dimensions
   integer, intent(in) :: iorder                ! 2rd or 4th order
   ! interpolation (2 or 4)
   integer, intent(in) :: logp                        ! use of log(p) in
   ! interpolation
   real, intent(in) :: vlaps(imax, jmax, kmax)
   real, intent(in) :: eta(kmax)                ! eta=(p-pt)/(ps-pt)
   real, intent(in) :: plaps(kmax)                ! laps pressure level
   real, intent(in) :: psurf(imax, jmax)        ! surface pressure
   real, intent(out) :: vmass(imax, jmax, kmax)

   ! local variables:
   integer :: i, j, k, l, mp, istart, iend
   real :: a(4), x(4), p

   ! conversion:
   do j = 1, jmax
      do i = 1, imax

         ! use pressure decreasing property from bottom to top:
         mp = 1        ! memorize the previous level searched

         ! convert:
         do k = 1, kmax

            ! mass pressure:
            p = eta(k)*(psurf(i, j) - plaps(kmax)) + plaps(kmax)

            ! search the interval of laps pressure levels
            ! containing the mass level:
            do l = mp, kmax
               if (p .ge. plaps(l)) then
                  mp = l                        ! memorize

                  ! find indices for interpolation:
                  istart = l - iorder/2
                  if (l .le. iorder/2) istart = 1
                  if (l + iorder/4 .gt. kmax) istart = kmax - iorder + 1
                  iend = istart + iorder - 1

                  ! check if log p is needed:
                  if (logp .eq. 1) then
                     x(1:iorder) = log(plaps(istart:iend))
                     p = log(p)
                  else
                     x(1:iorder) = plaps(istart:iend)
                  end if

                  ! second order:
                  if (iorder .eq. 2) &
                     call intplt2(p, x, vlaps(i, j, istart:iend), vmass(i, j, k))

                  ! fourth order:
                  if (iorder .eq. 4) &
                     call intplt4(p, x, vlaps(i, j, istart:iend), vmass(i, j, k))

                  ! end of search for k
                  goto 10

               end if
            end do

            ! search next k:
10          continue

         end do
      end do
   end do

end subroutine laps2mass

subroutine mass2laps(vlaps, imax, jmax, kmax, plaps, psurf, eta, &
                     iorder, logp, vmass)

!==========================================================
!  this routine converts the laps grid function into gsi
!  grid over mass coordinate.
!
!  history: jul. 2006 by yuanfu xie.
!==========================================================

   ! this routine converts a laps variable (vlaps) on laps
   ! coordinate into one on a mass coordinate.

   implicit none

   integer, intent(in) :: imax, jmax, kmax          ! 3d array dimensions
   integer, intent(in) :: iorder                ! 2rd or 4th order
   ! interpolation (2 or 4)
   integer, intent(in) :: logp                        ! use of log(p) in
   ! interpolation
   real, intent(out) :: vlaps(imax, jmax, kmax)
   real, intent(in) :: eta(kmax)                ! eta=(p-pt)/(ps-pt)
   real, intent(in) :: plaps(kmax)                ! laps pressure level
   real, intent(in) :: psurf(imax, jmax)        ! surface pressure
   real, intent(in) :: vmass(imax, jmax, kmax)

   ! local variables:
   integer :: i, j, k, l, mp, m, istart, iend
   real :: a(4), x(4), p

   ! conversion:
   do j = 1, jmax
      do i = 1, imax

         ! use pressure decreasing property from bottom to top:
         mp = 1        ! memorize the previous level searched

         ! convert:
         do k = 1, kmax

            ! search the interval of mass levels
            ! containing the laps pressure level:
            do l = mp, kmax
               p = eta(l)*(psurf(i, j) - plaps(kmax)) + plaps(kmax)

               if (plaps(k) .ge. p) then

                  mp = l

                  ! find indices for interpolation:
                  istart = l - iorder/2
                  if (l .le. iorder/2) istart = 1
                  if (l + iorder/4 .gt. kmax) istart = kmax - iorder + 1
                  iend = istart + iorder - 1

                  ! mass values:
                  do m = 1, iorder
                     x(m) = eta(istart + m - 1)* &
                            (psurf(i, j) - plaps(kmax)) + plaps(kmax)
                  end do
                  p = plaps(k)

                  ! check if log p is needed:
                  if (logp .eq. 1) then
                     x(1:iorder) = log(x(1:iorder))
                     p = log(plaps(k))
                  end if

                  ! second order:
                  if (iorder .eq. 2) &
                     call intplt2(p, x, vmass(i, j, istart:iend), vlaps(i, j, k))

                  ! fourth order:
                  if (iorder .eq. 4) &
                     call intplt4(p, x, vmass(i, j, istart:iend), vlaps(i, j, k))

                  ! end of search for k
                  goto 10

               end if

            end do

10          continue

         end do

      end do
   end do

end subroutine mass2laps

subroutine dryairmass(dam, pdam, imax, jmax, kmax, pres_1d, &
                      heights_3d, p_laps_bkg, t_laps_bkg)

!==========================================================
!  this routine computes the dry air mass and perturbation
!  using temperature field.
!
!  history: mar. 2006 by yuanfu xie.
!==========================================================

   implicit none

   integer, intent(in) :: imax, jmax, kmax          ! 3d array dimensions
   real, intent(out) :: dam(imax, jmax)   ! dry air mass in column
   real, intent(out) :: pdam(imax, jmax)  ! perturbation of dam
   real, intent(in) :: pres_1d(kmax)    ! pressure values@levels
   real, intent(in) :: heights_3d(imax, jmax, kmax)   ! heights
   real, intent(in) :: p_laps_bkg(imax, jmax, kmax)   ! p bkg
   real, intent(in) :: t_laps_bkg(imax, jmax, kmax)   ! t bkg

   ! local variables:
   integer :: i, j, k

   ! use pressure equation: dp/p = - g/(rt) dz [pp 20 holton]
   ! integration between p_top and p_surface:
   ! int_z=h(surface)^h(top) d ln(p) = - g/r int_h(s)^h(t) dz/t.
   ! thus: p(h(s)) = p(h(t))*exp(g/r*int_h(s)^h(t) dz/t, where
   ! p(h(t)) = pressure at the highest level: pres_1d(kmax).
   ! integration is replaced by riemann sum.
   do j = 1, jmax
      do i = 1, imax
         ! riemann sum:
         dam(i, j) = 0.0                ! initial summation
         pdam(i, j) = p_laps_bkg(i, j, 1)        ! default value at bottom
         do k = kmax, 2, -1
            ! summation til height=0
            if (heights_3d(i, j, k - 1) .gt. 0.0) then
               ! summation approximates integral:
               dam(i, j) = dam(i, j) + 2.0/ &
                           (t_laps_bkg(i, j, k) + t_laps_bkg(i, j, k - 1))* &
                           (heights_3d(i, j, k) - heights_3d(i, j, k - 1))
            else
               ! interpolate the full pressure between +/- heights:
               ! note: the values of these heights are the weights
               ! so the indices are swapped:
               pdam(i, j) = (p_laps_bkg(i, j, k - 1)*heights_3d(i, j, k) - &
                             p_laps_bkg(i, j, k)*heights_3d(i, j, k - 1))/ &
                            (heights_3d(i, j, k) - heights_3d(i, j, k - 1))

               ! riemann sum over the partial grid:
               ! positive height * 0.5*(t_k + interpolated t at z =0):
               dam(i, j) = dam(i, j) + 2.0*heights_3d(i, j, k)/ &
                           (t_laps_bkg(i, j, k) + &
                            (t_laps_bkg(i, j, k - 1)*heights_3d(i, j, k) - &
                             t_laps_bkg(i, j, k)*heights_3d(i, j, k - 1))/ &
                            (heights_3d(i, j, k) - heights_3d(i, j, k - 1)))
               goto 10
            end if
         end do
         ! dry air mass: integral of t -> dry surface pressure:
10       dam(i, j) = exp(9.806/287.0*dam(i, j))*pres_1d(kmax)
      end do
   end do

   ! perturbation pressure:
   pdam(1:imax, 1:jmax) = pdam(1:imax, 1:jmax) - dam(1:imax, 1:jmax)

end subroutine dryairmass

