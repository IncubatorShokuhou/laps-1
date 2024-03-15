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
        subroutine hum(t,td,rh,ni,nj,sat_t,sat_td)
c
c       g.s. stipanuk     1973            original version.
c       reference stipanuk paper entitled:
c            "algorithms for generating a skew-t, log p
c            diagram and computing selected meteorological
c            quantities."
c            atmospheric sciences laboratory
c            u.s. army electronics command
c            white sands missile range, new mexico 88002
c            33 pages
c       baker, schlatter  17-may-1982
c       p. stamus       02-08-89        changed to subroutine.
c                       06-14-89        pass in the sat_t/sat_td work arrays.
c                       08-21-89        add bad data check.
c                       09-20-89        add implicit none.
c                       04-13-90        drop scaling.
c                       11-15-90        bag bad data ck.  vectorize.
c
c   this function returns relative humidity (0.00 - 1.00) given the
c   temperature t and dew point td (k).  as calculated here,
c   relative humidity is the ratio of the actual vapor pressure to
c   the saturation vapor pressure.
c
c.....  note:  the temperature and dewpoint must be in degrees k.
c
        implicit none
        integer ni, nj, i, j
        real t(ni,nj), td(ni,nj), rh(ni,nj)
        real sat_t(ni,nj), sat_td(ni,nj)
c
c.....  calculate the actual and saturation vapor pressure, then the rh.
c
        call esat1(t,sat_t,ni,nj)
        call esat1(td,sat_td,ni,nj)
c
        do j=1,nj
        do i=1,ni
            rh(i,j) = (sat_td(i,j) / sat_t(i,j))
        enddo !i
        enddo !j
c
        return
        end
