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
        subroutine razm_lat_lon_gm(r4_slat,r4_slon,r4_range,
     1          r4_azimuth,r4_tlat,r4_tlon,status)
        include 'trigd.inc'
c***given a range and azimuth from some site, calculate the latitude and
c   longitude. range is considered to be a great circle arc length.

c       j. wakefield    25 apr 86       original version.
c       s. albers       8  may 90       calculations done in double precision
c                                       this prevents 'singularity' at 0 az.
c       s. albers      24  aug 92       change status from .true. to 1

c argument      i/o     type                    description
c --------      ---     ----    -----------------------------------------------
c slat,slon      i      r*4     site location (deg).
c range          i      r*4     target range in km from the site.
c azimuth        i      r*4     target azimuth in degrees.
c tlat,tlon      o      r*4     location of target (deg).
c status         o      i*4     standard system status.


!       assumed earth radius is 6370000m
        real*8          km_per_deg
        parameter      (km_per_deg=111.17747d0) ! km per degree of latitude

c***parameter list variables
        real          r4_slat,r4_slon,r4_range,r4_azimuth,r4_tlat,r4_t
     1lon
        real*8          slat,slon,range,azimuth,tlat,tlon
        integer       status

c***local variables
        real*8          dist,cosdlon,dlon

c***library symbols
!       external ss$_normal, edf__ivarg
!       logical*4       ltest_diag_gg

c***begin razm_lat_lon_gm ------------------------------------------------------

c***check input arguments
        if(abs(r4_slat).gt.90. .or. abs(r4_slon).gt.180. .or.
     1   r4_range.lt..0   .or.
     1   r4_azimuth.lt..0 .or. r4_azimuth.gt.360.)then
         status=0
!        if(ltest_diag_gg())call output_diag_gg('razm_lat_lon_gm',status)
         return
        else
         status=1
        endif

        azimuth = r4_azimuth
        range = r4_range
        slat = r4_slat
        slon = r4_slon

c***do it.  note that calculations are done in degrees.
        dist=range/km_per_deg

        tlat=dasind(dcosd(azimuth)*dsind(dist)*dcosd(slat) 
     1           + dsind(slat)*dcosd(dist))

        cosdlon=(dcosd(dist) - dsind(slat)*dsind(tlat)) 
     1            / (dcosd(slat)*dcosd(tlat))
        if(abs(cosdlon).gt.1.)cosdlon=sign(1.d0,cosdlon)
        dlon=dacosd(cosdlon)
        if(azimuth.ge..0.and.azimuth.le.180.)then
         tlon=slon+dlon
        else
         tlon=slon-dlon
        endif

        r4_tlat = tlat
        r4_tlon = tlon

        return
        end
