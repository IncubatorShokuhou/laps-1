cdis
cdis open source license/disclaimer, forecast systems laboratory
cdis noaa/oar/fsl, 325 broadway boulder, co 80305
cdis
cdis this software is distributed under the open source definition,
cdis which may be found at http://www.opensource.org/osd.html.
cdis
cdis in particular, redistribution and use in source and binary forms,
cdis with or without modification, are permitted provided that the
cdis following conditions are met:
cdis
cdis - redistributions of source code must retain this notice, this
cdis list of conditions and the following disclaimer.
cdis
cdis - redistributions in binary form must provide access to this
cdis notice, this list of conditions and the following disclaimer, and
cdis the underlying source code.
cdis
cdis - all modifications to this software must be clearly documented,
cdis and are solely the responsibility of the agent making the
cdis modifications.
cdis
cdis - if significant modifications or enhancements are made to this
cdis software, the fsl software policy manager
cdis(softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis this software and its documentation are in the public domain
cdis and are furnished "as is."the authors, the united states
cdis government, its instrumentalities, officers, employees, and
cdis agents make no warranty, express or implied, as to the usefulness
cdis of the software and documentation for any purpose.they assume
cdis no responsibility(1) for the use of the software and
cdis documentation; or(2) to provide technical support to users.
cdis
cdis
cdis
cdis
cdis

subroutine array_plot(a, b, imax, jmax, name, name_array, kmax, cld_hts
1 , scale)

!    ~1990        s. albers - ascii plots of cloud fields
!     1997 aug 01 k. dritz  - changed nx_l to imax and ny_l to jmax
!     1997 aug 01 k. dritz  - removed include of lapsparms.for
!     1997 oct    s. albers - make plots work with variable domain sizes.

integer imax, jmax
real a(imax, jmax), b(imax, jmax)
integer ia(imax, jmax)
character*1 c1a_array(imax, jmax), c1b_array(imax, jmax)
character*1 name_array(imax, jmax)
real cld_hts(kmax)
character name*10
character*1 c1_cov(-1:13)

integer max_plot
parameter(max_plot=65)
character c_mxp_a*(max_plot)
character c_mxp_b*(max_plot)

data c1_cov
1 /' ', '.', '1', '2', '3', '4', '5', '6', '7', '8', '9', '#', '*', ')', '.'/

c find max and min
amax = -1.e30
amin = 1.e30
!     sum=0.
cnt = 0.

do j = 1, jmax
do i = 1, imax
   ia(i, j) = 0
   if (a(i, j) .gt. 0.05 .or. b(i, j) .gt. 0.05) then
      cnt = cnt + 1.
   end if

   if (a(i, j) .ge. 0.00) then
      if (a(i, j) .gt. amax) amax = a(i, j)
      if (a(i, j) .lt. amin) amin = a(i, j)
!             sum=sum+a(i,j)
   end if

1  continue
end do
end do

if (cnt .eq. 0) then
   write (6, 1004) name
1004 format(1x, 'all values in array < 0.05 ', a10)
   return
else
!          write(6,1005) amax,amin,name
1005 format(1x, 2e12.3, ' max min ', a10)
end if

ihigh = imax

iskip = max(imax/48, 1)

nplot = (ihigh - 1)/iskip + 1

if (nplot .gt. max_plot) then ! prevent lines from getting too long
   nplot = max_plot
   ihigh = 1 + (nplot - 1)*iskip
end if

nspace = 3

jskip = max(jmax/32, 2)

write (6, *) 'iskip/jskip/nplot/ihigh', iskip, jskip, nplot, ihigh

do j = 1, jmax
do i = 1, ihigh, iskip
   c1a_array(i, j) = c1_cov(int(min(max(a(i, j)*10.*scale, -1.), 13.)))
   c1b_array(i, j) = c1_cov(int(min(max(b(i, j)*10.*scale, -1.), 13.)))

   if (name(1:4) .eq. 'horz') then ! horizontal section
      iil = max(i - iskip/2, 1)
      iih = min(i + (iskip - 1)/2, imax)
      jjl = max(j - jskip/2, 1)
      jjh = min(j + (jskip - 1)/2, jmax)

      do ii = iil, iih
      do jj = jjl, jjh
         if (name_array(ii, jj) .ne. ' ') then
            c1a_array(i, j) = name_array(ii, jj)
            c1b_array(i, j) = name_array(ii, jj)
         end if

      end do ! jj
      end do ! ii

   end if
end do
end do

if (name(1:4) .eq. 'vert') then     ! vertical section
   do j = jmax, 1, -2
      iarg = nint(cld_hts(j))

      do i1 = 1, max_plot
         i2 = 1 + (i1 - 1)*iskip
         if (i2 .le. ihigh) then
            c_mxp_a(i1:i1) = c1a_array(i2, j)
            c_mxp_b(i1:i1) = c1b_array(i2, j)
         end if
      end do ! i1

      write (6, 1001) iarg, c_mxp_a(1:nplot), c_mxp_b(1:nplot)
1001  format(i6, 2x, a, 2x, a)
   end do ! j

elseif (name(1:4) .eq. 'horz') then ! horizontal array plot
   do j = jmax, 1, -jskip

      do i1 = 1, max_plot          ! small grid
         i2 = 1 + (i1 - 1)*iskip ! large grid
         if (i2 .le. ihigh) then
            c_mxp_a(i1:i1) = c1a_array(i2, j)
            c_mxp_b(i1:i1) = c1b_array(i2, j)
         end if
      end do ! i1

      write (6, 1002) c_mxp_a(1:nplot), c_mxp_b(1:nplot)
1002  format(1x, a, 3x, a)
   end do ! j

end if

return
end

