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
  subroutine helicity(usig, vsig, zsig, terrain, imax, jmax, ksig, srelhel)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--
!--   name: helicity
!--
!--   purpose
!--   =======
!--   calculate storm relative helicity.
!--
!--   remarks
!--   =======
!--   helicity is equal to the vertical integral of:
!--
!--   (v - c) x dv/dz
!--
!--   where v is the model predicted wind and c is the
!--   estimated storm motion vector.
!--
!--   updates
!--   =======
!--   ??? 97   initial version......................................dnxm
!--   dec 98   modified call to wdir function.  an error in wdir was
!--            corrected............................................dnxm
!--   jan 99   changed variable names usigcrs and vsigcrs to usig and
!--            and vsig as winds are now at the cross point.........dnxm
!--   jan 01   adapted for use at fsl by b. shaw, noaa/fsl
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

     use constants
     implicit none

     real                         :: dz
     integer                      :: i
     integer, intent(in)      :: imax
     integer                      :: j
     integer, intent(in)      :: jmax
     integer                      :: k
     integer                      :: k3km
     integer                      :: k10km
     integer, intent(in)      :: ksig
     real, intent(out)     :: srelhel(imax, jmax)
     real                         :: stormdir
     real                         :: stormspd
     real                         :: stormu
     real                         :: stormv
     real                         :: sumdz
     real                         :: sumhel
     real                         :: sumu
     real                         :: sumv
     real, intent(in)      :: terrain(imax, jmax)
     real, intent(in)      :: usig(imax, jmax, ksig)
     real, intent(in)      :: vsig(imax, jmax, ksig)
     real, external        :: wdir
     real, external        :: wspd
     real, intent(in)      :: zsig(imax, jmax, ksig)

!-----------------------------------------------------------------------
!--   determine indices of 3 and 10 km layers (above ground level)
!-----------------------------------------------------------------------

     do j = 1, jmax
        do i = 1, imax

           do k = 1, ksig

              if (((zsig(i, j, k) - terrain(i, j)) .ge. 3000.0)) then
                 k3km = k
                 exit
              end if

           end do

           do k = ksig, k3km, -1

              if (((zsig(i, j, k) - terrain(i, j)) .le. 10000.0)) then
                 k10km = k
                 exit
              end if

           end do

!-----------------------------------------------------------------------
!--       estimate storm motion vector. speed is 75 percent of the mean
!--       wind between 3 and 10 km.  direction is 30 degrees to the
!--       right of the mean wind.
!-----------------------------------------------------------------------

           sumu = 0.0
           sumv = 0.0
           sumdz = 0.0

           do k = k3km, k10km

              dz = zsig(i, j, k) - zsig(i, j, k - 1)
              sumdz = sumdz + dz
              sumu = sumu + (0.5*dz*(usig(i, j, k) + usig(i, j, k - 1)))
              sumv = sumv + (0.5*dz*(vsig(i, j, k) + vsig(i, j, k - 1)))

           end do

           stormu = sumu/sumdz
           stormv = sumv/sumdz

           stormspd = wspd(stormu, stormv)*0.75

!-----------------------------------------------------------------------
!--       when calling wdir, send in a cone factor of zero
!--       so that direction is grid relative.
!-----------------------------------------------------------------------

           stormdir = wdir(stormu, stormv, 0.0, 0.0, 0.0) + 30.0
           if (stormdir .gt. 360.0) stormdir = stormdir - 360.0

           stormu = -stormspd*sin(stormdir*deg2rad)
           stormv = -stormspd*cos(stormdir*deg2rad)

           sumhel = 0.0

!-----------------------------------------------------------------------
!--       calculate helicity.  integrate between the ground and 3 km,
!--       a depth that is frequently used in the literature.
!-----------------------------------------------------------------------

           do k = 1, k3km

              sumhel = sumhel + &
                       ((usig(i, j, k) - stormu)* &
                        (vsig(i, j, k + 1) - vsig(i, j, k))) - &
                       ((vsig(i, j, k) - stormv)* &
                        (usig(i, j, k + 1) - usig(i, j, k)))

           end do

           srelhel(i, j) = -sumhel

        end do
     end do
     return
  end subroutine helicity
