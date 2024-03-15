cdis forecast systems laboratory
cdis noaa/oar/erl/fsl
cdis 325 broadway
cdis boulder, co 80303
cdis
cdis forecast research division
cdis local analysis and prediction branch
cdis laps
cdis
cdis this software and its documentation are in the public domain and
cdis are furnished "as is."the united states government, its
cdis instrumentalities, officers, employees, and agents make no
cdis warranty, express or implied, as to the usefulness of the software
cdis and documentation for any purpose.they assume no responsibility
cdis(1) for the use of the software and documentation; or(2) to provide
cdis technical support to users.
cdis
cdis permission to use, copy, modify, and distribute this software is
cdis hereby granted, provided that the entire disclaimer notice appears
cdis in all copies.all modifications to this software must be clearly
cdis documented, and are solely the responsibility of the agent making
cdis the modifications.if significant modifications or enhancements
cdis are made to this software, the fsl software policy manager
cdis(softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
c
c
subroutine regress_precip(num_sfc, radar, gauge, a_t, b_t, rbar, gbar
1 , istatus)
c
c*******************************************************************************
c
c routine to calculate regression coefficients from gauge/radar data
c(formerly t and elev data) .

c values returned are the coefficients for the regression equations
c for the temperatures and dew points.also calculates the mean
c elevation for the surface stations.
c
c changes:
c p.a.stamus 12 - 01 - 88 original(from j.mcginley)
c
c inputs/outputs:
c
c variable var type i/o description
c - -------------------------------------
c num_sfc i i number of surface stations.
c radar ra i station elevation.
c gauge ra i temperature.
c a_t r o 'a'"       ""    " (slope)
c b_t r o 'b'regression value for temp. (intrcp)
c rbar r o mean elevation of the stations.
c
c user notes:
c
c 1.units are not changed in this routine.
c
c*******************************************************************************
c
real radar(num_sfc), gauge(num_sfc)

c
badflag = -99.9
c
c.....set up storage variables.
c
cnt = 0.
sumht = 0.
sumh = 0.
sumt = 0.
sumh2 = 0.
sumt2 = 0.
c
c.....gather sums and then calculate the 'a'and 'b'for the regression
c.....equation y = az + b, for both the temperature and dew point.the
c.....'b'is the intercept with sea level, and the 'a'is the lapse rate.
c.....'y'represents the gauge value and 'z'is radar/analyzed
c.....also calculate the mean elevation of the stations.
c
istatus = 0

gmax = 0.
rmax = 0.
gmin = 999.
rmin = 999.

write (6, *) '     gauge analyzed'
do 10 i = 1, num_sfc
!         if(radar(i).le.badflag .or. gauge(i).le.badflag) go to 10
   write (6, 5) i, gauge(i), radar(i)
5  format(i3, 2f8.3)
   sumht = (radar(i)*gauge(i)) + sumht
   sumh = radar(i) + sumh
   sumh2 = (radar(i)*radar(i)) + sumh2
   sumt = gauge(i) + sumt
   cnt = cnt + 1.

   rmin = min(radar(i), rmin)
   rmax = max(radar(i), rmax)
   gmin = min(gauge(i), gmin)
   gmax = max(gauge(i), gmax)

   istatus = 1
10 continue

   if (cnt .eq. 0.) then
      write (6, *) ' count = 0'
      istatus = 0
      return
   end if

!       slope
   denominator = (cnt*sumh2 - sumh*sumh)
   if (denominator .ne. 0.) then
      a_t = (cnt*sumht - sumh*sumt)/denominator
   else
      a_t = 0.
      istatus = 0
   end if

!       intercept
   b_t = (sumt - a_t*sumh)/cnt
   c
   rbar = sumh/cnt
   gbar = sumt/cnt

   write (6, *) 'num_sfc,rbar,gbar,a_t,b_t'
   write (6, *) num_sfc, rbar, gbar, a_t, b_t

   write (6, 11) gmin, gmax
11 format('  gauge range = ', 2f9.3)
   write (6, 12) rmin, rmax
12 format('  radar range = ', 2f9.3)

   if (a_t .lt. 0.1 .or. a_t .gt. 10.) then
      write (6, *) ' warning, slope is ill conditioned'
      istatus = 0
   end if
   c
   c.....end of routine
   c
   return
end
