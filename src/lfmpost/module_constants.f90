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

module constants

   ! constants used by mm5post

   real, parameter        :: cm2inch = 1.0/2.54
   real, parameter        :: cp = 1005.7
   real, parameter        :: grav = 9.81
   integer, parameter        :: imiss = -99999
   integer, parameter        :: ismth = 0
   integer, parameter        :: jsmth = 0
   real, parameter        :: lapse = 6.5e-03
   real, parameter        :: lv = 2.5e+06
   real, parameter        :: m2feet = 3.281
   real, parameter        :: mps2knts = 1.944
   real, parameter        :: mps2mph = 2.237
   real, parameter        :: pi = 3.1415927
   real, parameter        :: p0 = 100000.0
   real, parameter        :: r = 287.04
   real, parameter        :: rv = 461.5
   real, parameter        :: t0 = 273.15
   real, parameter        :: xmiss = -99999.9
   real, parameter        :: cpog = cp/grav
   real, parameter        :: cpor = cp/r
   real, parameter        :: deg2rad = pi/180.0
   real, parameter        :: e = r/rv
   real, parameter        :: gor = grav/r
   real, parameter        :: kappa = r/cp
   real, parameter        :: rad2deg = 180.0/pi
   real, parameter        :: rog = r/grav
   real, parameter        :: rvolv = rv/lv

end module constants

