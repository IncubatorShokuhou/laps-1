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

  subroutine clouds(nx, ny, nz, grid_spacing, cldliqmr, cldicemr, snowmr, &
                    height, topo, cldbase, cldtop, ceiling, cldamt)

     ! this routine is used to compute cloud ceiling (agl), cloud base
     ! height (msl), cloud top height (msl), and coverage (fraction) given
     ! mixing ratios of the various microphysical species.

     ! adapted from the usaf weather agency mm5 post processer.
     !  brent shaw, noaa forecast systems lab, dec 2000

     implicit none

     ! arguments
     integer, intent(in)             :: nx      ! array x-dimension
     integer, intent(in)             :: ny      ! array y-dimension
     integer, intent(in)             :: nz      ! array z-dimension
     real, intent(in)                :: grid_spacing  ! in meters
     real, intent(in)                :: cldliqmr(nx, ny, nz) ! cloud liq mix rat
     real, intent(in)                :: cldicemr(nx, ny, nz) ! cloud ice mix rat
     real, intent(in)                :: snowmr(nx, ny, nz) ! snow mix rat
     real, intent(in)                :: height(nx, ny, nz) ! heights in meters
     real, intent(in)                :: topo(nx, ny)    ! topography in m
     real, intent(out)               :: cldbase(nx, ny)
     real, intent(out)               :: cldtop(nx, ny)
     real, intent(out)               :: ceiling(nx, ny)
     real, intent(out)               :: cldamt(nx, ny)

     ! locals
     integer                         :: i, j, k
     real                            :: icethresh
     real                            :: liqthresh
     real                            :: snowthresh

     ! intialize the arrays

     cldbase = 0.0
     cldtop = 0.0
     cldamt = 0.0
     ceiling = 1e+30  ! unlimited ceiling

     ! set thresholds for cloud determination based on grid
     ! resolution.  kind of hokey, but will get us by for now

     if (grid_spacing .le. 10000) then
        icethresh = 0.000005
        snowthresh = 0.000003
        liqthresh = 0.000003
     else
        icethresh = 0.000005
        snowthresh = 0.000025
        liqthresh = 0.000025
     end if

     ! loop through using these thresholds to determine cloudiness

     do j = 1, ny
        do i = 1, nx
           find_base: do k = 1, nz
              if ((cldliqmr(i, j, k) .ge. liqthresh) .or. &
                  (cldicemr(i, j, k) .ge. icethresh) .or. &
                  (snowmr(i, j, k) .ge. snowthresh)) then
                 cldbase(i, j) = height(i, j, k)
                 ceiling(i, j) = height(i, j, k) - topo(i, j)

                 ! for now all we can do is use cldamt as a yes/no
                 ! we should look at coming up with a fraction function

                 cldamt(i, j) = 1.0
                 exit find_base
              end if
           end do find_base
           find_top: do k = nz, 1, -1
              if ((cldliqmr(i, j, k) .ge. liqthresh) .or. &
                  (cldicemr(i, j, k) .ge. icethresh) .or. &
                  (snowmr(i, j, k) .ge. snowthresh)) then
                 cldtop(i, j) = height(i, j, k)
                 exit find_top
              end if
           end do find_top
        end do
     end do

     return
  end subroutine clouds

