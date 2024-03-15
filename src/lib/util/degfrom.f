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
cdis
cdis   
cdis
        subroutine degfrom(xd,yd,deg_from)

        real deg_from,xd,yd,xx

             if(yd .eq. 0.0 .and. xd .eq. 0.0)then
             deg_from = 360.
             goto 1019
             endif

c   north of current storm.

             if(yd .eq. 0.0 .and. xd .gt. 0.0)then
             deg_from = 360.
             goto 1019
             endif

c   south of current storm.

             if(yd .eq. 0.0 .and. xd .lt. 0.0)then
             deg_from = 180.
             goto 1019
             endif

             xx = xd/yd
             deg_from = atand(xx)

c   left half.

             if(yd .lt. 0.0)then
             deg_from = 270. - deg_from
             endif

c   right half.

             if(yd .gt. 0.0)then
             deg_from = 90. - deg_from
             endif


1019     return
        end
