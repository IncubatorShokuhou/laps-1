!dis    forecast systems laboratory
!dis    noaa/oar/erl/fsl
!dis    325 broadway
!dis    boulder, co     80303
!dis
!dis    forecast research division
!dis    local analysis and prediction branch
!dis    laps
!dis
!dis    this software and its documentation are in the public domain and
!dis    are furnished "as is."  the united states government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  they assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis    technical support to users.
!dis
!dis    permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  all modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  if significant modifications or enhancements
!dis    are made to this software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

program gsi_prep

!==========================================================
!  this program is a three-dimensional variational analysis
!  for laps. it can be either gsi or stmas depending on the
!  variational namelist parameter.
!
!  history:
!         creation: yuanfu xie        3-2006
!==========================================================

   implicit none

   ! laps configuration:
   call laps_conf

   ! laps prep:
   call laps_bkgd

   ! laps obs:
   call laps_obsv

   ! conversion to gsi:
   call gsi_bkobs

   ! deallocate dynamic memory:
   call laps_remv

   ! radar level iii bufr data:
   call gsi_radarbufr

   ! noaa polar satellite 1d bufr data: amsu-a, amsu-b, and hirs
   ! call gsi_noaa1d2bufr

end program gsi_prep
