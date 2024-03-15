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

  subroutine the2t(thetae, p, tparcel)

     ! calculates temperature at any pressure level along a saturation
     ! adiabat by iteratively solving eq 2.76 from wallace and hobbs (1977).

     ! adpated from usaf routine, b. shaw, noaa/fsl

     use constants
     implicit none

     integer                          :: iter
     real, external                   :: mixsat
     real, intent(in)                 :: p
     real                             :: tcheck
     real                             :: theta
     real, intent(in)                 :: thetae
     real                             :: tovtheta
     real, intent(out)                :: tparcel
     logical                          :: converged

     converged = .false.
     tovtheta = (p/1000.0)**kappa
     tparcel = thetae/exp(lv*0.012/(cp*295.0))*tovtheta

     iterloop: do iter = 1, 50

        theta = thetae/exp(lv*mixsat(tparcel, p*100.)/(cp*tparcel))
        tcheck = theta*tovtheta

        if (abs(tparcel - tcheck) .lt. 0.05) then
           converged = .true.
           exit iterloop
        end if
        tparcel = tparcel + (tcheck - tparcel)*0.3

     end do iterloop
     if (.not. converged) then
        print *, 'warning!! thetae to temp calc did not converge!'
        print *, 'thetae and p:', thetae, p
     end if

     return
  end subroutine the2t
