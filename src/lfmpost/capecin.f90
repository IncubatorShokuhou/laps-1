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

  subroutine capecin(psig, tsig, thetaesig, thetasig, &
                     rhsig, zsigfull, &
                     tprs, li, &
                     posbuoyen, negbuoyen, k500, imax, jmax, ksig, kprs)

!-----------------------------------------------------------------------
!--
!--   purpose
!--   =======
!--   calculate positive buoyant energy (or convective available
!--   energy) and negative buoyant energy (or convective inhibition).
!--
!--   also calculate lifted index.
!--
!--   reference
!--   =========
!--   doswell and rasmussen (1994), wea and fcsting, p 625.
!--
!--   updates
!--   =======
!--   4 jan 01  - adapted from usaf weather agency routine
!--               b. shaw, noaa/fsl
!-----------------------------------------------------------------------

     use constants
     implicit none

     integer                     :: i
     integer, intent(in)       :: imax
     integer                     :: j
     integer, intent(in)       :: jmax
     integer                     :: k
     integer, intent(in)       :: k500
     integer                     :: kmax
     integer, intent(in)       :: ksig
     integer, intent(in)       :: kprs
     real, intent(out)      :: li(imax, jmax)
     real, intent(out)      :: negbuoyen(imax, jmax)
     real, intent(out)      :: posbuoyen(imax, jmax)
     real                        :: prslcl
     real, intent(in)       :: psig(imax, jmax, ksig)
     real, intent(in)       :: rhsig(imax, jmax, ksig)
     real, intent(in)       :: thetasig(imax, jmax, ksig)
     real, intent(in)       :: thetaesig(imax, jmax, ksig)
     real                        :: thw
     real                        :: thwmax
     real, external         :: tlcl
     real                        :: tmplcl
     real                        :: tparcel
     real, intent(in)       :: tprs(imax, jmax, kprs)
     real, intent(in)       :: tsig(imax, jmax, ksig)
     real, external         :: wobf
     real, intent(in)       :: zsigfull(imax, jmax, ksig + 1)
     real                        :: deltaz, dtheta, thetaparcel
     real, external              :: potential_temp
     real, parameter             :: pref = 1000.
     real                        :: psave
     logical                     :: compute_cin
     real, allocatable           :: buoy(:)
     posbuoyen = 0.0
     negbuoyen = 0.0
     allocate (buoy(ksig))
     do j = 1, jmax
        do i = 1, imax

           thwmax = -9999.0

           find_most_unstable: do k = 1, ksig

!           ------------------------------------------------------------
!           pick the most unstable parcel in the lowest 50 mb as
!           indicated by the sigma level with the highest wet bulb
!           potential temperature.  store index in kmax.
!           ------------------------------------------------------------

              if (((psig(i, j, 1) - psig(i, j, k)) .lt. 50.0)) then

                 thw = thetaesig(i, j, k) - wobf(thetaesig(i, j, k) - t0)

                 if (thw .gt. thwmax) then
                    kmax = k
                    thwmax = thw
                 end if

              else

                 exit find_most_unstable

              end if

           end do find_most_unstable

!         --------------------------------------------------------------
!         calculate lifted index by lifting the most unstable
!         parcel.
!         --------------------------------------------------------------

           call the2t(thetaesig(i, j, kmax), 500.0, tparcel)
           li(i, j) = tprs(i, j, k500) - tparcel

!         --------------------------------------------------------------
!         calculate the temperature and pressure of the lifting
!         condensation level.
!         --------------------------------------------------------------

           tmplcl = tlcl(tsig(i, j, kmax), rhsig(i, j, kmax))
           prslcl = psig(i, j, kmax)*(tmplcl/tsig(i, j, kmax))**cpor

!         --------------------------------------------------------------
!         calculate the buoyancy.
!         --------------------------------------------------------------
           posbuoyen(i, j) = 0.
           negbuoyen(i, j) = 0.
           do k = kmax, ksig

!           ------------------------------------------------------------
!           above the lcl, calculate virtual temperature of the
!           parcel as it moves along a moist adiabat.  below the
!           lcl, lift parcel along a dry adiabat.
!           ------------------------------------------------------------

              if (psig(i, j, k) .le. prslcl) then

                 call the2t(thetaesig(i, j, kmax), psig(i, j, k), tparcel)
              else

                 tparcel = thetasig(i, j, kmax)/(pref/psig(i, j, k))**kappa

              end if

              ! compute the potential temperature of the parcel
              thetaparcel = potential_temp(tparcel, psig(i, j, k)*100.)
              dtheta = thetaparcel - thetasig(i, j, k)
              deltaz = zsigfull(i, j, k + 1) - zsigfull(i, j, k)
              buoy(k) = deltaz*dtheta/thetasig(i, j, k)
           end do

           ! now loop through the column again, partitioning the buoyency
           ! into positive (cape) and negative (cin) component.  we terminate
           ! the contribution to cin when/if a layer of cape greater than
           ! 150 mb deep is found.

           compute_cin = .true.
           psave = -100.
           do k = kmax, ksig
              if (buoy(k) .gt. 0.) then
                 if (psave .lt. 0) then
                    psave = psig(i, j, k)
                 else
                    if ((psave - psig(i, j, k)) .gt. 150.) then
                       compute_cin = .false.
                    end if
                 end if
                 posbuoyen(i, j) = posbuoyen(i, j) + buoy(k)
              else if (buoy(k) .lt. 0.) then
                 psave = -100.
                 if (compute_cin) then
                    negbuoyen(i, j) = negbuoyen(i, j) + buoy(k)
                 end if
              end if
           end do
        end do
     end do

     posbuoyen = grav*posbuoyen
     negbuoyen = grav*negbuoyen

     ! cap the negative buoyancy to a maximum value of 700 j/kg

     where (negbuoyen .lt. -700) negbuoyen = -700.
     deallocate (buoy)
  end subroutine capecin
