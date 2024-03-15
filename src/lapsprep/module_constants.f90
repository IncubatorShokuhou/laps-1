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

module lapsprep_constants
   real, parameter :: pi = 3.14159265358
   real, parameter :: radius_of_earth = 6370. ! mm5/wrf
   real, parameter :: radians_per_degree = pi/180.
   real, parameter :: rdry = 287.1  ! dry air gas constant
   real, parameter :: g = 9.81      ! gravity

   ! microphysics constants for autoconversion of cloud liquid to
   ! rain and ice to snow, in kg/m**3
   real, parameter :: autoconv_lwc2rai = 0.0005
   real, parameter :: autoconv_ice2sno = 0.0005
   real, parameter :: lwc_min = 0.000001
   real, parameter :: ice_min = 0.000001
   real, parameter :: lcp_min = 0.0 ! simpler and more aggressive humidification
end module lapsprep_constants
