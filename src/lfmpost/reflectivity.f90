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

  subroutine reflectivity(nx, ny, nz, rho, hgt, rainmr, icemr, snowmr, graupelmr, &
                          refl, max_refl, echo_tops)

     ! subroutine to compute estimated radar reflectivity from
     ! the precipitation mixing ratios.  will also return
     ! column max reflectivity and echo tops.  the estimation
     ! is done using formulae from kessler (1969) and
     ! rogers and yau (1989).
     !
     ! adapted from usaf weather agency routine.
     !
     !  brent shaw, noaa forecast system lab, dec 2000
     !

     implicit none

     integer, intent(in)             :: nx
     integer, intent(in)             :: ny
     integer, intent(in)             :: nz
     real, intent(in)                :: rho(nx, ny, nz)
     real, intent(in)                :: hgt(nx, ny, nz)
     real, intent(in)                :: rainmr(nx, ny, nz)
     real, intent(in)                :: icemr(nx, ny, nz)
     real, intent(in)                :: snowmr(nx, ny, nz)
     real, intent(in)                :: graupelmr(nx, ny, nz)
     real, intent(out)               :: refl(nx, ny, nz)
     real, intent(out)               :: max_refl(nx, ny)
     real, intent(out)               :: echo_tops(nx, ny)

     ! locals
     real, parameter                 :: svnfrth = 7.0/4.0
     real, parameter                 :: max_top_thresh = 5.0
     integer                         :: i, j, k
     real                            :: w
     max_refl = 0.0
     echo_tops = 1.0e37
     refl = 0.0

     do j = 1, ny
        do i = 1, nx
           do k = 1, nz

              ! compute the basic reflectivity
              !refl(i,j,k) =17300.0 * &
              !            (rho(i,j,k) * 1000.0 * &
              !             max(0.0,rainmr(i,j,k)))**svnfrth

              ! add the ice component
              !refl(i,j,k)=refl(i,j,k) + &
              !       38000.0*(rho(i,j,k) * 1000.0 * &
              !       max(0.0,icemr(i,j,k)+snowmr(i,j,k)+graupelmr(i,j,k)))**2.2

              ! convert to dbz
              !refl(i,j,k) = 10.*alog10(max(refl(i,j,k),1.0))

              ! test rams reflectivity algorithm
              w = 264083.11*(rainmr(i, j, k) + &
                             0.2*(icemr(i, j, k) + snowmr(i, j, k)) + &
                             2.0*graupelmr(i, j, k))
              w = max(1., w)
              refl(i, j, k) = 17.8*alog10(w)
              ! since we are going from the ground up, we can
              ! check threshold and set echo top.

              if (refl(i, j, k) .ge. max_top_thresh) echo_tops(i, j) = hgt(i, j, k)

           end do
           ! compute the max value in the column
           max_refl(i, j) = maxval(refl(i, j, :))
        end do
     end do
     return
  end subroutine reflectivity

