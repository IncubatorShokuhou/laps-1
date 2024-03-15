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


        subroutine plotfield (data_anal,igrid,jgrid)


c       $log: plotfield.for,v $
c revision 1.1  1996/08/30  20:48:28  birk
c initial revision
c

c ported to the unix environment 1 oct 1993

        implicit none

        integer igrid,jgrid


        real data_anal(igrid,jgrid)

        real finc




        integer ioffp,ioffm
        real spval,epsval,cntmin,cntmax,cntint

        common /conre1/ioffp,spval,epsval,cntmin,cntmax,cntint,ioffm

        ioffp = 1   ! says that the missing value exists.
        spval = 1.0e+37  ! sets the missing value flag for plots
        ioffm = 0 !omit the conrec message on the plot -- doesn't work





c        call opngks


        ioffm = 1




        call conrec (data_anal,igrid,igrid,jgrid,0.,0.,finc,-1,0,0)
        call frame

c        call clsgks


        end
