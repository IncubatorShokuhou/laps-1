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

subroutine intplt2(p, x, vin, vout)

!==========================================================
!  this routine interpolates a value using a second order
!  scheme.
!
!  input:
!        p: position where to interpolate;
!        x: positions of 2 grid points;
!        vin: input grid values;
!  output:
!        vout: interpolated value;
!
!  reference: intpl2.doc on mac notebook.
!
!  history:
!        jul. 2006 by yuanfu xie.
!==========================================================

   implicit none

   real, intent(in) :: p, x(2), vin(2)
   real, intent(out) :: vout

   integer :: i, j
   real :: e(2), r, b, c

   ! check the 4 grid points:
   do i = 1, 2
      do j = 1, 2
         if ((i .ne. j) .and. &
             (abs(x(i) - x(j)) .lt. 1.0e-5*x(i))) then
            print *, 'intplt2: there are overlap gridpoints, stop'
            stop
         end if
      end do
   end do

   do i = 1, 2
      if (abs(p - x(i)) .lt. 1.0e-5*maxval(x)) then
         e = 0.0
         e(i) = 1.0
         goto 1
      end if
   end do

   ! e2:
   e(2) = (p - x(1))/(x(2) - x(1))

   ! e1:
   e(1) = 1.0 - e(2)

1  continue

   ! interpolate:
   vout = 0.0
   do i = 1, 2
      vout = vout + e(i)*vin(i)
   end do

end subroutine intplt2

subroutine intplt4(p, x, vin, vout)

!==========================================================
!  this routine interpolates a value using a fourth order
!  scheme.
!
!  input:
!        p: position where to interpolate;
!        x: positions of 4 grid points;
!        vin: input grid values;
!  output:
!        vout: interpolated value;
!
!  reference: intpl2.doc on mac notebook.
!
!  history:
!        jul. 2006 by yuanfu xie.
!==========================================================

   implicit none

   real, intent(in) :: p, x(4), vin(4)
   real, intent(out) :: vout

   integer :: i, j
   real :: e(4), r, b, c

   ! check the 4 grid points:
   do i = 1, 4
      do j = 1, 4
         if ((i .ne. j) .and. &
             (abs(x(i) - x(j)) .lt. 1.0e-5*x(i))) then
            print *, 'intplt4: there are overlap gridpoints, stop'
            stop
         end if
      end do
   end do

   do i = 1, 4
      if (abs(p - x(i)) .lt. 1.0e-5*maxval(x)) then
         e = 0.0
         e(i) = 1.0
         goto 1
      end if
   end do

   ! e4:
   r = 1.0
   c = 1.0
   do i = 1, 3
      r = r*(p - x(i))
      c = c*(x(4) - x(i))
   end do
   e(4) = r/c

   ! e3:
   r = 1.0
   b = 1.0
   c = 1.0
   do i = 1, 2
      r = r*(p - x(i))
      c = c*(x(4) - x(i))
      b = b*(x(3) - x(i))
   end do
   e(3) = (r - c*e(4))/b

   ! e2:
   e(2) = (p - x(1) - (x(4) - x(1))*e(4) - (x(3) - x(1))*e(3))/(x(2) - x(1))

   ! e1:
   e(1) = 1.0 - e(2) - e(3) - e(4)

1  continue

   ! interpolate:
   vout = 0.0
   do i = 1, 4
      vout = vout + e(i)*vin(i)
   end do

end subroutine intplt4
