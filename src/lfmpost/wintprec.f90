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
  subroutine wintprec(tsig, zsig, zprs, psfc, tsfc, &
                      terrain, prcpinc, imax, jmax, &
                      ksig, kprs, k700, k850, &
                      k1000, prcpconv, preciptype)

!-----------------------------------------------------------------------
!--   (winter) precipitation algorithm
!--
!--   purpose:  this product identifies areas of precipitation and the
!--   type expected, based on mm5 data.  the process is essentially two-
!--   fold.  first, areas of precipitation are identified when the mm5
!--   precipitation array (prcpinc) exceeds 0.01 inch.
!--
!--   second, thickness thresholds are used at the gridpoints for which
!--   it has been determined that precipitation is occurring and the
!--   surface pressure is greater than 850mb (i.e., non-mountainous
!--   regions).  the thickness thresholds utilized are based on
!--   meteorological research from the pertinent sources
!--   (e.g., mwr, w&f), and are as follows:
!--
!--                               (thick1)     (thick2)
!--                             1000mb-850mb  850mb-700mb  <--thickness
!--
!--   precipitation type:   rain   gt 1310      gt 1540
!--
!--                freezing rain   gt 1310      gt 1540 [sig1 t < 0co]
!--
!--                    ice/mixed   le 1310      gt 1540
!--
!--                         snow   le 1310      le 1540.
!--
!--   over mountainous terrain, precipitation type is limited to either
!--   rain or snow.  this is consistent with climatic data presented in
!--   "a regional climatology of freezing precipitation for the
!--   contiguous united states" (10th conf. on applied climatology, 20-
!--   24 oct 97).  where a precipitation occurrence has been determined,
!--   the temperatures in the lowest 1500 m are checked.  if all are
!--   below freezing, snow is forecasted; otherwise rain is predicted.
!--
!--   modification:  added ability to predict regions where thunderstorm
!--   activity may occur.  prior to exiting the main loop, a check is
!--   made:  where rain is predicted and the convective component of the
!--   precip exceeds 0.01", forecast for thunderstorms.
!--
!--   updates
!--   =======
!--   jan 2001 1998  initial version, adapted from usaf weather agency
!--       brent shaw, noaa/forecast systems lab
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------

     use constants
     implicit none

     integer                        :: i
     integer, intent(in)        :: imax
     integer                        :: j
     integer, intent(in)        :: jmax
     integer                        :: k
     integer, intent(in)        :: k700
     integer, intent(in)        :: k850
     integer, intent(in)        :: k1000
     integer                        :: kprs
     integer                        :: ksig
     integer                        :: k1500
     real, intent(in)        :: prcpconv(imax, jmax)
     real, intent(in)        :: prcpinc(imax, jmax)
     integer, intent(out)       :: preciptype(imax, jmax)
     real, intent(in)        :: psfc(imax, jmax)
     real, intent(in)        :: terrain(imax, jmax)
     real, intent(in)        :: tsfc(imax, jmax)
     real                           :: tsfcf
     real, intent(in)        :: tsig(imax, jmax, ksig)
     real                           :: thickhigh
     real                           :: thicklow
     real                           :: tsig1, tsig2, tsig3
     real, intent(in)        :: zprs(imax, jmax, kprs)
     real, intent(in)        :: zsig(imax, jmax, ksig)
     real, external          :: fahren

     do j = 1, jmax
        do i = 1, imax

!-----------------------------------------------------------------------
!--       the threshold for calculating precip type is 0.0001 meter
!--       per time period.
!-----------------------------------------------------------------------

           if (prcpinc(i, j) .le. 0.0001) then

              preciptype(i, j) = 0

           else

!-----------------------------------------------------------------------
!--         check the surface pressure to determine whether
!--         high or low-elevation logic is used. 850 mb is
!--         the cutoff.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!--         low elevation grid point.
!-----------------------------------------------------------------------

              if (psfc(i, j) .gt. 85000.0) then

!-----------------------------------------------------------------------
!--           calculate thicknesses that will be used to determine
!--           precip type.
!-----------------------------------------------------------------------

                 thicklow = zprs(i, j, k850) - zprs(i, j, k1000)

                 thickhigh = zprs(i, j, k700) - zprs(i, j, k850)

!-----------------------------------------------------------------------
!--           rain, or if surface temperature is below freezing,
!--           freezing rain.
!-----------------------------------------------------------------------

                 if ((thicklow .gt. 1310.0) .and. (thickhigh .gt. 1540.0)) then

                    if (tsfc(i, j) .ge. t0) then
                       preciptype(i, j) = 1
                    else
                       preciptype(i, j) = 3
                    end if

!-----------------------------------------------------------------------
!--           ice/mixed.
!-----------------------------------------------------------------------

                 elseif ((thicklow .le. 1310.0) .and. &
                         (thickhigh .gt. 1540.0)) then

                    preciptype(i, j) = 4

!-----------------------------------------------------------------------
!--           rain or snow.
!-----------------------------------------------------------------------

                 elseif ((thicklow .le. 1310.0) .and. &
                         (thickhigh .le. 1540.0)) then

                    tsfcf = fahren(tsfc(i, j) - t0)

                    if (tsfcf .ge. 37.0) then
                       preciptype(i, j) = 1
                    else
                       preciptype(i, j) = 5
                    end if

!-----------------------------------------------------------------------
!--           rain.
!-----------------------------------------------------------------------

                 else

                    preciptype(i, j) = 1

                 end if

!-----------------------------------------------------------------------
!--         high terrain grid point.
!-----------------------------------------------------------------------

              else

!-----------------------------------------------------------------------
!--           find top of 1500 m agl layer.
!-----------------------------------------------------------------------

                 do k = 1, ksig

                    if ((zsig(i, j, k) - terrain(i, j)) .ge. 1500.0) then

                       k1500 = k
                       exit

                    end if

                 end do
!-----------------------------------------------------------------------
!--           if the model top is ever lowered, the above code
!--           could fail on the top of a high mountain.
!-----------------------------------------------------------------------

                 k1500 = max(k1500, 1)

!-----------------------------------------------------------------------
!--           find temperature at the bottom, top and middle of
!--           1500 m agl layer.  for middle layer, recycle variable
!--           k1500.
!-----------------------------------------------------------------------

                 tsig1 = tsig(i, j, 1)

                 tsig3 = tsig(i, j, k1500)

                 k1500 = max(1, nint(float(k1500)/2.0))

                 tsig2 = tsig(i, j, k1500)

!----------------------------------------------------------------------
!--           snow.
!-----------------------------------------------------------------------

                 if ((tsig1 .lt. t0) .and. &
                     (tsig2 .lt. t0) .and. &
                     (tsig3 .lt. t0)) then

                    preciptype(i, j) = 5

!-----------------------------------------------------------------------
!--           rain & check for thunderstorms.
!-----------------------------------------------------------------------

                 else

                    preciptype(i, j) = 1

                 end if

              end if

           end if

           if ((preciptype(i, j) .eq. 1) .and. &
               (prcpconv(i, j) .ge. 0.001)) then

              preciptype(i, j) = 2         ! thunderstorm

           end if

        end do
     end do
     return
  end subroutine wintprec

