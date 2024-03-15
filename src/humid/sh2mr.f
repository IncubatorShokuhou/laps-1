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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis



        subroutine sh2mr (sh,mr,istatus)


c       $log: sh2mr.for,v $
c revision 1.1  1996/08/30  20:57:36  birk
c initial revision
c

c       given the specific humidity (g/g)
c       this routine computes the mixing ratio (g/g)

        real sh,mr
        real epsilon
        integer istatus

        epsilon = 0.622 !(ratio of molecular weight of water ~ 18,
c                       over dry air~29)

        istatus  = 1

        if (sh/epsilon .eq. 1.) then
                write(6,*) 'error in input to routine sh2mr'
                write(6,*) 'divide by zero has occurred'
                istatus = 0
        endif


        if (sh/epsilon .gt. 1.) then
                write(6,*) 'error in input to routine sh2mr'
                write(6,*) 'impossible condition has occurred'
                write(6,*) 'sh has value of ', sh
                istatus = 0
        endif

        mr = sh / (1.-sh/epsilon)

        return
        end





