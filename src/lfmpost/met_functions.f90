
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function esat(t)
     !computes saturation vapor pressure (pa) from temperature
     ! note: computed with respect to liquid!!!
     use constants
     real                       :: esat
     real, intent(in)           :: t

     esat = 611.21*exp((17.502*(t - t0))/(t - 32.18))
  end function esat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mixsat(t, p)
     ! computes saturation vapor mixing ratio as function of temp and press
     use constants
     implicit none

     real, external             :: esat
     real                      :: mixsat
     real, intent(in)          :: p
     real                      :: satvpr
     real, intent(in)          :: t

     satvpr = esat(t)
     mixsat = (e*satvpr)/(p - satvpr)

  end function mixsat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function potential_temp(t, p)

     ! computes potential temperature from temp (k) and pressure (pa)

     use constants
     implicit none
     real, intent(in)           :: t
     real, intent(in)           :: p
     real                       :: potential_temp

     potential_temp = t*(p0/p)**kappa

  end function potential_temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function eq_potential_temp(t, p, w, rh)

     ! calculates equivalent potential temperature given temperature,
     ! pressure, mixing ratio, and relative humidity via
     ! bolton's equation. (mwr 1980, p 1052, eq. 43)

     implicit none

     real, intent(in)             :: t  ! temp in k
     real, intent(in)             :: p  ! pressure in pa
     real, intent(in)             :: w  ! mixing ratio kg/kg
     real, intent(in)             :: rh ! relative humidity (fraction)
     ! real, external               :: potential_temp
     real                         :: thtm
     real, external               :: tlcl
     real                         :: eq_potential_temp

     thtm = t*(100000./p)**(2./7.*(1.-(0.28*w)))
     eq_potential_temp = thtm* &
                         exp((3.376/tlcl(t, rh) - 0.00254)* &
                             (w*1000.0*(1.0 + 0.81*w)))

  end function eq_potential_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function tlcl(t, rh)

     ! computes the temperature of the lifting condensation level
     ! given surface t and rh using bolton's equation (mwr 1980, p 1048, #22)

     implicit none
     real                               :: denom
     real, intent(in)                   :: rh
     real, intent(in)                   :: t
     real                               :: term1
     real                               :: term2
     real                               :: tlcl

     term1 = 1.0/(t - 55.0)
     term2 = alog(rh/1.0)/2840.
     denom = term1 - term2
     tlcl = (1.0/denom) + 55.0

  end function tlcl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function dewpt(t, rh)

     ! compute dew point from temperature (k) and rh (fraction)
     use constants
     implicit none

     real                :: dewpt
     real, intent(in)    :: rh
     real, intent(in)    :: t

     dewpt = t/((-rvolv*alog(max(rh, 0.01))*t) + 1.0)

  end function dewpt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function relhum(t, mixrat, p)
     ! computes relative humidity (fraction form)
     implicit none
     real, intent(in)                 :: p ! pressure in pa
     real, intent(in)                 :: mixrat ! vapor mixing ratio
     real, external                   :: mixsat ! saturation vapor mix. ratio
     real                             :: relhum  !(fraction)
     real, intent(in)                 :: t

     relhum = mixrat/mixsat(t, p)
  end function relhum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function fahren(t)

     ! converts from celsius to fahrenheit

     implicit none

     real, intent(in)          :: t
     real                      :: fahren

     fahren = (1.8*t) + 32.
  end function fahren
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function celsius(tf)

     implicit none
     real, intent(in)             :: tf
     real                         :: celsius
     ! converts fahrenheit to celsius

     celsius = (5./9.)*(tf - 32.0)

  end function celsius
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function wdir(u, v, xlon, orient, conefact)
     ! computes wind direction from u/v components.  converts the direction
     ! to true direction if conefactor != 0.

     use constants
     implicit none

     real, intent(in)               :: conefact
     real                           :: diff
     real, intent(in)               :: orient
     real, intent(in)               :: u
     real, intent(in)               :: v
     real                           :: wdir
     real, intent(in)               :: xlon

     ! handle case where u is very small to prevent divide by 0

     if (abs(u) .lt. 0.001) then
        if (v .le. 0.0) then
           wdir = 0.0
        else
           wdir = 180.
        end if

        ! otherwise, standard trig problem!

     else
        wdir = 270.-(atan2(v, u)*rad2deg)

        if (wdir .gt. 360.) then
           wdir = wdir - 360.
        end if
     end if

     ! change to earth relative

     diff = (orient - xlon)*conefact
     if (diff .gt. 180.0) diff = diff - 360.
     if (diff .lt. -180.) diff = diff + 360.

     wdir = wdir - diff
     if (wdir .gt. 360.0) wdir = wdir - 360.
     if (wdir .lt. 0.) wdir = wdir + 360.

  end function wdir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function wspd(u, v)

     ! computes wind velocity from u/v components
     implicit none
     real, intent(in)        :: u
     real, intent(in)        :: v
     real                    :: wspd

     wspd = sqrt((u*u) + (v*v))
  end function wspd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function wobf(t)

!-----------------------------------------------------------------------
!--
!--   name: wobuf function
!--
!--   this function calculates the difference of the wet bulb potential
!--   temperatures for saturated and dry air given the temperature.
!--
!--   it was created by herman wobus of the navy weather research
!--   facility from data in the smithsonian meteorological tables.
!--
!--
!--   let wbpts = wet bulb potential temperature for saturated air
!--               at temperature t in celsius
!--
!--   let wbptd = wet bult potential temperature for dry air at
!--               the same temperature.
!--
!--   the wobus function wobf (in degrees celsius) is defined by:
!--
!--               wobf(t) = wbpts - wbptd
!--
!--   although wbpts and wbptd are functions of both pressure and
!--   temperature, their difference is a function of temperature only.
!--
!--   the wobus function is useful for evaluating several thermodynamic
!--   quantities.
!--
!--   if t is at 1000 mb, then t is potential temperature pt and
!--   wbpts = pt.  thus,
!--
!--               wobf(pt) = pt - wbptd
!--
!--   if t is at the condensation level, then t is the condensation
!--   temperature tc and wbpts is the wet bulb potential temperature
!--   wbpt.  thus,
!--
!--               wobf(tc) = wbpt - wbptd
!--
!--   manipulating the above equations we get,
!--
!--               wbpt = pt - wobf(pt) + wobf(tc)   and
!--
!--               wbpts = pt - wobf(pt) + wobf(t)
!--
!--   if t is equivalent potential temperature ept (implying that
!--   the air at 1000 mb is completely dry), then
!--
!--               wbpts = ept and wbptd = wbpt, thus,
!--
!--               wobf(ept) = ept - wbpt
!--
!-----------------------------------------------------------------------

     implicit none

     real                        :: pol
     real, intent(in)         :: t
     real                        :: wobf
     real                        :: x

     x = t - 20.0

     if (x .le. 0.0) then

        pol = 1.0 + x*(-8.8416605e-03 &
                       + x*(1.4714143e-04 + x*(-9.6719890e-07 &
                                               + x*(-3.2607217e-08 + x*(-3.8598073e-10)))))

        wobf = 15.130/(pol**4)

     else

        pol = 1.0 + x*(3.6182989e-03 &
                       + x*(-1.3603273e-05 + x*(4.9618922e-07 &
                                                + x*(-6.1059365e-09 + x*(3.9401551e-11 &
                                                                         + x*(-1.2588129e-13 + x*(1.6688280e-16)))))))

        wobf = (29.930/(pol**4)) + (0.96*x) - 14.8

     end if

  end function wobf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function heatindex(temp_k, rh_pct)

     ! computes heat index from a temperature (k) and rh (%).

     use constants
     implicit none
     real, external                  :: celsius
     real, external                  :: fahren
     real                            :: heatindex
     real, intent(in)                :: rh_pct
     real                            :: rh_pct_sqr
     real                            :: tf
     real                            :: tf_sqr
     real, intent(in)                :: temp_k

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

  end function heatindex

  function tcvp(p, mr, z, rho, nz)
     use constants
     implicit none

     real                          :: tcvp
     integer, intent(in)           :: nz
     real, intent(in)              :: p(nz)
     real, intent(in)              :: mr(nz)
     real, intent(in)              :: z(nz)
     real, intent(in)              :: rho(nz)

     integer   :: k, kbot
     real      :: pvapor(nz)
     real      :: mrmean, dz

     ! set top vapor pressure to 0
     pvapor(nz) = 0.

     ! integrate moisture downward

     do kbot = nz - 1, 1, -1

        pvapor(kbot) = 0.

        do k = nz - 1, kbot, -1

           ! compute dz and mean qv for this layer
           dz = z(k + 1) - z(k)
           mrmean = (mr(k) + mr(k + 1))*0.5

           pvapor(kbot) = pvapor(kbot) + grav*mrmean*rho(k)*dz/(1.+mrmean)
        end do
     end do
     tcvp = pvapor(1)
     return
  end function tcvp

