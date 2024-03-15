cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis 

       function zenith(lat,lon,slat,slon)

c      computes zenith angle give earth lat,lon and satellite subpoint
c      slat, slon.

c      returns zenith angle in degrees for use in ssec routines

c      input coordinates are all in radians.  e longitude =+

c      authored by dan birkenheuer   12/1/95


        real lat,lon,slat,slon
        real zenith
        real cos_of_zen


        cos_of_zen = sin (slat)*sin(lat)+cos(slat)*cos(lat)*
     1      cos(abs(slon-lon))

        zenith = acos(cos_of_zen)

        zenith = zenith + atan2 ( sin(zenith),6.6166-cos_of_zen )

c       convert to degrees

        zenith = zenith * 180./acos(-1.0)

        return

      end




