cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
c
c
        subroutine mean_pres(num_sfc,p_s,pbar)
c
c*******************************************************************************
c
c       routine to calculate the mean surface or reduced pressure.
c
c       changes:
c               p. a. stamus    12-21-88        original
c
c       inputs/outputs:
c
c          variable     var type     i/o     description
c         ----------   ----------   -----   -------------
c          num_sfc         i          i      number of surface stations.
c          p_s             ra         i      surface (or reduced) pressure.
c          pbar            r          o      mean pressure of all the stations.
c
c       user notes:
c
c       1.  units are not changed in this routine.
c
c*******************************************************************************
c
        real p_s(num_sfc)
c
        sump = 0.
        cntp = 0.
c
        do 1 i=1,num_sfc
          if(p_s(i) .le. 0.) go to 1
          sump = sump + p_s(i)
          cntp = cntp + 1.
1       continue
c
        pbar = sump / cntp
c
        return
        end
