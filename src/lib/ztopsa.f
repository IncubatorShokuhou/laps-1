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

      real function ztopsa(z)

c*  this routine converts a height in meters into a pressure in a standard
c*  atmosphere in millibars.

      implicit none

      include 'constants.inc' ! psamslpa, gamma

      real t0,p0,p11,z11,c1,c2,z,flag,flg

      data flag,flg/1e-30,199998./
      data t0/288./
      data c1,c2/5.256,14600./
      data z11,p11/11000.,226.0971/

      p0 = psamslpa / 100.

      if (z.gt.flg) then
          ztopsa=flag
      else if (z.lt.z11) then
          ztopsa=p0*((t0-gamma*z)/t0)**c1
      else
          ztopsa=p11*10.**((z11-z)/c2)
      end if

      return
      end
