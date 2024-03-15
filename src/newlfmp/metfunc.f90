function esat(t)

! computes saturation vapor pressure (pa) from temperature (k).
!  note: computed with respect to liquid.

   use constants

   implicit none

   real :: esat, t

   esat = 611.21*exp((17.502*(t - t0))/(t - 32.18))

   return
end

!===============================================================================

function potential_temp(t, p)

! computes potential temperature from temp (k) and pressure (pa).

   use constants

   implicit none

   real :: t, p, potential_temp

   potential_temp = t*(p0/p)**kappa

   return
end

!===============================================================================

function relhum(t, mixrat, p)

! computes relative humidity (fraction form)

   implicit none

   real :: t       ! temperature (k)
   real :: p       ! pressure (pa)
   real :: mixrat  ! vapor mixing ratio
   real :: mixsat  ! saturation vapor mix. ratio
   real :: relhum  ! rh (fraction)

   relhum = mixrat/mixsat(t, p)

   return
end

!===============================================================================

function mixsat(t, p)

! computes saturation vapor mixing ratio as function of temp (k) and press( pa).

   use constants

   implicit none

   real :: esat, mixrat, p, satvpr, t, mixsat

   satvpr = esat(t)
   mixsat = (e*satvpr)/(p - satvpr)

   return
end

!===============================================================================

function eq_potential_temp(t, p, w, rh)

! calculates equivalent potential temperature given temperature,
! pressure, mixing ratio, and relative humidity via
! bolton's equation.  (mwr 1980, p 1052, eq. 43)

   implicit none

   real :: t   ! temp (k)
   real :: p   ! pressure (pa)
   real :: w   ! mixing ratio (kg/kg)
   real :: rh  ! relative humidity (fraction)
   real :: thtm, tlcl, eq_potential_temp

   thtm = t*(100000./p)**(2./7.*(1.-(0.28*w)))
   eq_potential_temp = thtm*exp((3.376/tlcl(t, rh) - 0.00254) &
                                *(w*1000.0*(1.0 + 0.81*w)))

   return
end

!===============================================================================

function tlcl(t, rh)

! computes the temperature of the lifting condensation level
! given surface t and rh using bolton's equation.  (mwr 1980, p 1048, #22)

   implicit none

   real :: denom, rh, t, term1, term2, tlcl

   term1 = 1.0/(t - 55.0)
   term2 = alog(rh/1.0)/2840.
   denom = term1 - term2
   tlcl = (1.0/denom) + 55.0

   return
end

!===============================================================================

function dewpt(t, rh)

! compute dew point from temperature (k) and rh (fraction).

   use constants

   implicit none

   real :: dewpt, rh, t

   dewpt = t/((-rvolv*alog(max(rh, 0.01))*t) + 1.0)
   dewpt = min(t, dewpt)

   return
end

!===============================================================================

function wobf(t)

! name: wobuf function
!
! this function calculates the difference of the wet bulb potential
! temperatures for saturated and dry air given the temperature.
!
! it was created by herman wobus of the navy weather research
! facility from data in the smithsonian meteorological tables.
!
!
! let wbpts = wet bulb potential temperature for saturated air
!             at temperature t in celsius
!
! let wbptd = wet bult potential temperature for dry air at
!             the same temperature.
!
! the wobus function wobf (in degrees celsius) is defined by:
!
!             wobf(t) = wbpts - wbptd
!
! although wbpts and wbptd are functions of both pressure and
! temperature, their difference is a function of temperature only.
!
! the wobus function is useful for evaluating several thermodynamic
! quantities.
!
! if t is at 1000 mb, then t is potential temperature pt and
! wbpts = pt.  thus,
!
!             wobf(pt) = pt - wbptd
!
! if t is at the condensation level, then t is the condensation
! temperature tc and wbpts is the wet bulb potential temperature
! wbpt.  thus,
!
!             wobf(tc) = wbpt - wbptd
!
! manipulating the above equations we get,
!
!             wbpt = pt - wobf(pt) + wobf(tc)   and
!
!             wbpts = pt - wobf(pt) + wobf(t)
!
! if t is equivalent potential temperature ept (implying that
! the air at 1000 mb is completely dry), then
!
!             wbpts = ept and wbptd = wbpt, thus,
!
!             wobf(ept) = ept - wbpt

   implicit none

   real :: pol, t, wobf, x

   x = t - 20.

   if (x <= 0.) then
      pol = 1.0 + x*(-8.8416605e-03 &
                     + x*(1.4714143e-04 + x*(-9.6719890e-07 &
                                             + x*(-3.2607217e-08 + x*(-3.8598073e-10)))))
      wobf = 15.130/(pol**4)
   else
      pol = 1.+x*(3.6182989e-03 &
                  + x*(-1.3603273e-05 + x*(4.9618922e-07 &
                                           + x*(-6.1059365e-09 + x*(3.9401551e-11 &
                                                                    + x*(-1.2588129e-13 + x*(1.6688280e-16)))))))
      wobf = (29.930/(pol**4)) + (0.96*x) - 14.8
   end if

   return
end

!===============================================================================

subroutine the2t(thetae, p, tparcel)

! calculates temperature k) at any pressure (mb) level along a saturation
! adiabat by iteratively solving eq 2.76 from wallace and hobbs (1977).

! adpated from usaf routine, b. shaw, noaa/fsl

   use constants

   implicit none

   integer           :: iter
   real, external    :: mixsat
   real, intent(in)  :: p
   real              :: tcheck
   real              :: theta
   real, intent(in)  :: thetae
   real              :: tovtheta
   real              :: arg
   real, intent(out) :: tparcel
   logical           :: converged

   converged = .false.
   tovtheta = (p/1000.)**kappa
   tparcel = thetae/exp(lv*0.012/(cp*295.0))*tovtheta

   do iter = 1, 50
      arg = lv*mixsat(tparcel, p*100.)/(cp*tparcel)
      if (arg .le. 80.) then
         theta = thetae/exp(arg)
      else
         print *, 'warning, large argument for exponent ', arg
         goto 999
      end if
      tcheck = theta*tovtheta
      if (abs(tparcel - tcheck) < 0.05) then
         converged = .true.
         exit
      end if
      tparcel = tparcel + (tcheck - tparcel)*0.3
   end do

999 continue
   if (.not. converged) then
      print *, 'warning: thetae to temp calc did not converge.'
      print *, '  thetae and p:', thetae, p
      print *, '  tparcel:', tparcel
      print *, '  theta,tovtheta:', theta, tovtheta
   end if

   return
end

!===============================================================================

function wdir(u, v, xlon, orient, conefact)

! computes wind direction from u/v components.  converts the direction
!   to true direction if conefactor \= 0.

   use constants

   implicit none

   real :: conefact, diff, orient, u, v, wdir, xlon

! handle case where u is very small to prevent divide by 0.

   if (abs(u) < 0.001) then
      if (v <= 0.) then
         wdir = 0.
      else
         wdir = 180.
      end if

! otherwise, standard trig problem.

   else
      wdir = 270.-(atan2(v, u)*rad2deg)
      if (wdir > 360.) wdir = wdir - 360.
   end if

! change to earth relative.

   diff = (orient - xlon)*conefact
   if (diff > 180.) diff = diff - 360.
   if (diff < -180.) diff = diff + 360.

   wdir = wdir - diff
   if (wdir > 360.) wdir = wdir - 360.
   if (wdir < 0.) wdir = wdir + 360.

   return
end

!===============================================================================

function wspd(u, v)

! computes wind velocity from u/v components.

   implicit none

   real :: u, v, wspd

   wspd = sqrt((u*u) + (v*v))

   return
end

!===============================================================================

function fahren(t)

! converts from celsius to fahrenheit

   implicit none

   real :: t, fahren

   fahren = (1.8*t) + 32.

   return
end

!===============================================================================

function celsius(tf)

   implicit none

   real, parameter :: factor = 5./9.
   real :: tf, celsius

   celsius = factor*(tf - 32.)

   return
end

!===============================================================================

function heatindex(temp_k, rh_pct)

! computes heat index from a temperature (k) and rh (%).

   use constants

   implicit none

   real :: celsius, fahren, heatindex, rh_pct, rh_pct_sqr, tf, tf_sqr, temp_k

   rh_pct_sqr = rh_pct*rh_pct
   tf = fahren(temp_k - t0)
   tf_sqr = tf*tf

   heatindex = -42.379 + (2.04901523*tf) &
               + (10.1433312*rh_pct) &
               - (0.22475541*tf*rh_pct) &
               - (6.83783e-03*tf_sqr) &
               - (5.481717e-02*rh_pct_sqr) &
               + (1.22874e-03*tf_sqr*rh_pct) &
               + (8.52e-04*rh_pct_sqr*tf) &
               - (1.99e-06*tf_sqr*rh_pct_sqr)

   heatindex = celsius(heatindex) + t0

   return
end

!===============================================================================

function tw(t, td, p)

! this function returns the wet-bulb temperature tw (k)
!  given the temperature t (k), dew point td (k)
!  and pressure p (pa).  see p.13 in stipanuk (1973), referenced
!  above, for a description of the technique.

! baker,schlatter        17-may-1982        original version

! determine the mixing ratio line thru td and p.

   aw = w(td, p)

! determine the dry adiabat thru t and p.

   ao = o(t, p)
   pi = p

! iterate to locate pressure pi at the intersection of the two
!  curves.  pi has been set to p for the initial guess.

   do i = 1, 10
      x = 0.02*(tmr(aw, pi) - tda(ao, pi))
      if (abs(x) < 0.01) goto 5
      pi = pi*(2.**(x))
   end do
5  continue

! find the temperature on the dry adiabat ao at pressure pi.

   ti = tda(ao, pi)

! the intersection has been located...now, find a saturation
!  adiabat thru this point.  function os returns the equivalent
!  potential temperature (k) of a parcel saturated at temperature
!  ti and pressure pi.

   aos = os(ti, pi)

! function tsa returns the wet-bulb temperature (k) of a parcel at
!  pressure p whose equivalent potential temperature is aos.

   tw = tsa(aos, p)

   return
end

!===============================================================================

function w(t, p)

! this function returns the mixing ratio (kg/kg) given the temperature t (k)
!  and pressure (pa).  the formula is quoted in most meteorological texts.

! baker,schlatter        17-may-1982        original version

   implicit none

   real :: t, p, x, esat, w

   x = esat(t)
   w = 0.622*x/(p - x)

   return
end

!===============================================================================

function o(t, p)

! g.s. stipanuk     1973                original version.
! reference stipanuk paper entitled:
!  "algorithms for generating a skew-t, log p
!     diagram and computing selected meteorological quantities."
!     atmospheric sciences laboratory
!     u.s. army electronics command
!     white sands missile range, new mexico 88002
!     33 pages
! baker, schlatter  17-may-1982

! this function returns potential temperature (k) given
!  temperature t (k) and pressure p (pa) by solving the poisson
!  equation.

   implicit none

   real :: t, p, o

   o = t*((100000./p)**.286)

   return
end

!===============================================================================

function tmr(w, p)

!        g.s. stipanuk     1973                original version.
!        reference stipanuk paper entitled:
!            "algorithms for generating a skew-t, log p
!             diagram and computing selected meteorological
!             quantities."
!             atmospheric sciences laboratory
!             u.s. army electronics command
!             white sands missile range, new mexico 88002
!             33 pages
!        baker, schlatter  17-may-1982

!   this function returns the temperature (k) on a mixing
!   ratio line w (kg/kg) at pressure p (pa). the formula is given in
!   table 1 on page 7 of stipanuk (1973).

!   initialize constants

   data c1/.0498646455/, c2/2.4082965/, c3/7.07475/
   data c4/38.9114/, c5/.0915/, c6/1.2035/

   x = alog10(w*p*0.01/(0.622 + w))
   tmr = 10.**(c1*x + c2) - c3 + c4*((10.**(c5*x) - c6)**2.)

   return
end

!===============================================================================

function tda(th, p)

!        g.s. stipanuk     1973                original version.
!        reference stipanuk paper entitled:
!            "algorithms for generating a skew-t, log p
!             diagram and computing selected meteorological
!             quantities."
!             atmospheric sciences laboratory
!             u.s. army electronics command
!             white sands missile range, new mexico 88002
!             33 pages
!        baker, schlatter  17-may-1982

! this function returns the temperature tda (k) on a dry adiabat
!  at pressure p (pa).  the dry adiabat is given by
!  potential temperature th (k).  the computation is based on
!  poisson's equation.

   implicit none

   real :: th, p, tda

   tda = th*((p*0.00001)**.286)

   return
end

!===============================================================================

function os(t, p)

!        g.s. stipanuk     1973                original version.
!        reference stipanuk paper entitled:
!            "algorithms for generating a skew-t, log p
!             diagram and computing selected meteorological
!             quantities."
!             atmospheric sciences laboratory
!             u.s. army electronics command
!             white sands missile range, new mexico 88002
!             33 pages
!        baker, schlatter  17-may-1982

! this function returns the equivalent potential temperature os
!  (k) for a parcel of air saturated at temperature t (k)
!  and pressure p (pa).

   implicit none

   real :: t, p, b, w, os

   data b/2651.8986/

! b is an empirical constant approximately equal to the latent heat
!  of vaporization for water divided by the specific heat at constant
!  pressure for dry air.

   os = t*((100000./p)**.286)*(exp(b*w(t, p)/t))

   return
end

!===============================================================================

function tsa(os, p)

!        g.s. stipanuk     1973                original version.
!        reference stipanuk paper entitled:
!            "algorithms for generating a skew-t, log p
!             diagram and computing selected meteorological
!             quantities."
!             atmospheric sciences laboratory
!             u.s. army electronics command
!             white sands missile range, new mexico 88002
!             33 pages
!        baker, schlatter  17-may-1982

! this function returns the temperature tsa (k) on a saturation
!  adiabat at pressure p (pa). os is the equivalent potential
!  temperature of the parcel (k). sign(a,b) replaces the
!  algebraic sign of a with that of b.
!  b is an empirical constant approximately equal to 0.001 of the latent
!  heat of vaporization for water divided by the specific heat at constant
!  pressure for dry air.

   data b/2651.8986/

! tq is the first guess for tsa.

   tsa = 253.15

! d is an initial value used in the iteration below.

   d = 120.

! iterate to obtain sufficient accuracy....see table 1, p.8
!  of stipanuk (1973) for equation used in iteration.

   do i = 1, 12
      d = d/2.
      x = os*exp(-b*w(tsa, p)/tsa) - tsa*((100000./p)**.286)
      if (abs(x) < 1.e-7) goto 2
      tsa = tsa + sign(d, x)
   end do
2  continue

   return
end
