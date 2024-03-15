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
!dis
!dis
!dis
!dis

subroutine staggerz(vin, nx, ny, nz, ptop, psfc, znw, iorder, logp, vout)

!==========================================================
!  this routine computes a grid function on a vertically
!  stagger grid using a given uniform grid function with
!  options of 2nd/4th order accuracy and logp transfer.
!
!  input:
!        ptop:         top pressure;
!        psfc:        sfc pressure;
!        znw:        grid spacing ratio relative to psfc (nz);
!        vin:        input grid function on a uniform grid;
!                (nx*ny*nz)
!        nx:        number of uniform grid point in x;
!        ny:        number of uniform grid point in y;
!        nz:        number of uniform vertical points;
!        logp:        0 -- standard interpolation;
!                1 -- log(p) interpolation;
!       iorder: 2 or 4 for second or fourth order;
!
!  output:
!        vout:        output grid funcion on a stagger grid.
!                (nx*ny*(nz-1))
!
!  history: jul. 2006 by yuanfu xie.
!==========================================================

   implicit none

   integer, intent(in) :: nx, ny, nz, logp, iorder
   real, intent(in) ::  ptop, psfc(nx, ny), znw(nz)
   real, intent(in) ::    vin(nx, ny, nz)
   real, intent(out) ::   vout(nx, ny, nz - 1)

   ! local variables:
   integer :: i, j, k, l, m
   real :: p, x(4)

   ! lower boundary:
   do j = 1, ny
      do i = 1, nx

         do l = 1, iorder
            x(l) = znw(l)*(psfc(i, j) - ptop) + ptop
         end do

         ! first stagger grid point:
         p = 0.5*(x(1) + x(2))

         ! log(p):
         if (logp .eq. 1) then
            x = log(x)
            p = log(p)
         end if

         ! interpolation:
         if (iorder .eq. 2) &
            call intplt2(p, x, vin(i, j, 1:iorder), vout(i, j, 1))
         if (iorder .eq. 4) &
            call intplt4(p, x, vin(i, j, 1:iorder), vout(i, j, 1))

      end do
   end do

   ! interior:
   do k = 2, nz - 2

      m = k - iorder/2

      do j = 1, ny
         do i = 1, nx

            do l = 1, iorder
               x(l) = znw(m + l)*(psfc(i, j) - ptop) + ptop
            end do

            ! first stagger grid point:
            p = 0.5*(x(iorder/2) + x(iorder/2 + 1))

            ! log(p):
            if (logp .eq. 1) then
               x = log(x)
               p = log(p)
            end if

            ! interpolation:
            if (iorder .eq. 2) &
               call intplt2(p, x, vin(i, j, k:k + 1), vout(i, j, k))
            if (iorder .eq. 4) &
               call intplt4(p, x, vin(i, j, k - 1:k + 2), vout(i, j, k))

         end do
      end do

   end do

   ! upper boundary:
   m = nz - iorder
   do j = 1, ny
      do i = 1, nx

         do l = 1, iorder
            x(l) = znw(m + l)*(psfc(i, j) - ptop) + ptop
         end do

         ! first stagger grid point:
         p = 0.5*(x(iorder - 1) + x(iorder))

         ! log(p):
         if (logp .eq. 1) then
            x = log(x)
            p = log(p)
         end if

         ! interpolation:
         if (iorder .eq. 2) &
            call intplt2(p, x, vin(i, j, nz - 1:nz), vout(i, j, nz - 1))
         if (iorder .eq. 4) &
            call intplt4(p, x, vin(i, j, nz - 3:nz), vout(i, j, nz - 1))

      end do
   end do

end subroutine staggerz

subroutine unstaggerz(vin, nx, ny, nz, ptop, psfc, znu, iorder, logp, vout)

!==========================================================
!  this routine computes a grid function on a uniform grid
!  using a given stagger grid function with options of a
!  second/fourth order accuracy and log p transfer.
!
!  input:
!        ptop:         top pressure;
!        psfc:        sfc pressure;
!        znw:        eta value for a uniform grid;
!        vin:        input grid function on a stagger grid;
!                (nx*ny*(nz-1))
!        nx:        number of uniform grid point in x;
!        ny:        number of uniform grid point in y;
!        nz:        number of uniform vertical points;
!        logp:        0 -- standard interpolation;
!                1 -- log(p) interpolation;
!        iorder: 2 or 4 order accuracy;
!
!  output:
!        vout:        output grid funcion on a uniform grid.
!                (nx*ny*nz)
!
!  history: jul. 2006 by yuanfu xie.
!==========================================================

   implicit none

   integer, intent(in) :: nx, ny, nz, logp, iorder
   real, intent(in) ::  ptop, psfc(nx, ny), znu(nz - 1)
   real, intent(in) ::    vin(nx, ny, nz - 1)
   real, intent(out) ::   vout(nx, ny, nz)

   ! local variables:
   integer :: i, j, k, l, m
   real :: p, x(4)

   ! lower boundary:
   do j = 1, ny
      do i = 1, nx

         do l = 1, iorder
            x(l) = znu(l)*(psfc(i, j) - ptop) + ptop
         end do

         ! first stagger grid point:
         p = 1.5*x(1) - 0.5*x(2)

         ! log(p):
         if (logp .eq. 1) then
            x = log(x)
            p = log(p)
         end if

         ! interpolation:
         if (iorder .eq. 2) &
            call intplt2(p, x, vin(i, j, 1:2), vout(i, j, 1))
         if (iorder .eq. 4) &
            call intplt4(p, x, vin(i, j, 1:4), vout(i, j, 1))

         ! second stagger grid point:
         p = 0.5*(x(1) + x(2))

         ! log(p):
         if (logp .eq. 1) then
            p = log(p)
         end if

         ! interpolation:
         if (iorder .eq. 2) &
            call intplt2(p, x, vin(i, j, 1:2), vout(i, j, 2))
         if (iorder .eq. 4) &
            call intplt4(p, x, vin(i, j, 1:4), vout(i, j, 2))

      end do
   end do

   ! interior:
   do k = 3, nz - 2

      m = k - iorder/2 - 1

      do j = 1, ny
         do i = 1, nx

            do l = 1, iorder
               x(l) = znu(m + l)*(psfc(i, j) - ptop) + ptop
            end do

            ! first stagger grid point:
            p = 0.5*(x(iorder/2) + x(iorder/2 + 1))

            ! log(p):
            if (logp .eq. 1) then
               x = log(x)
               p = log(p)
            end if

            ! interpolation:
            if (iorder .eq. 2) &
               call intplt2(p, x, vin(i, j, k - 1:k), vout(i, j, k))
            if (iorder .eq. 4) &
               call intplt4(p, x, vin(i, j, k - 2:k + 1), vout(i, j, k))

         end do
      end do

   end do

   ! upper boundary:
   m = nz - iorder - 1
   do j = 1, ny
      do i = 1, nx

         do l = 1, iorder
            x(l) = znu(m + l)*(psfc(i, j) - ptop) + ptop
         end do

         ! first stagger grid point from the top:
         p = 1.5*x(iorder) - 0.5*x(iorder - 1)

         ! log(p):
         if (logp .eq. 1) then
            x = log(x)
            p = log(p)
         end if

         ! interpolation:
         if (iorder .eq. 2) &
            call intplt2(p, x, vin(i, j, nz - 2:nz - 1), vout(i, j, nz))
         if (iorder .eq. 4) &
            call intplt4(p, x, vin(i, j, nz - 4:nz - 1), vout(i, j, nz))

         ! second stagger grid point:
         p = 0.5*(x(iorder - 1) + x(iorder))

         ! log(p):
         if (logp .eq. 1) then
            p = log(p)
         end if

         ! interpolation:
         if (iorder .eq. 2) &
            call intplt2(p, x, vin(i, j, nz - 2:nz - 1), vout(i, j, nz - 1))
         if (iorder .eq. 4) &
            call intplt4(p, x, vin(i, j, nz - 4:nz - 1), vout(i, j, nz - 1))

      end do
   end do

end subroutine unstaggerz
