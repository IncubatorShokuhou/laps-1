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
        subroutine esat1(t,es,ni,nj)
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
c       p. stamus         01-26-89      changed to subroutine.
c                         08-21-89      add bad data check.
c                         09-20-89      add implicit none.
c                         11-14-90      bad bad data ck. vectorize.
c
c   this function returns the saturation vapor pressure over
c   water (mb) given the temperature (k).
c   the algorithm is due to nordquist, w.s.,1973: "numerical approxima-
c   tions of selected meteorlolgical parameters for cloud physics prob-
c   lems," ecom-5475, atmospheric sciences laboratory, u.s. army
c   electronics command, white sands missile range, new mexico 88002.
c
c.....  note: t must be in degrees k coming into this routine.
c
        implicit none
        integer ni, nj, i, j
        real p1, p2, c1, term
        real t(ni,nj), es(ni,nj)
c
        do j=1,nj
        do i=1,ni
          p1 = 11.344 - 0.0303998 * t(i,j)
          p2 = 3.49149 - 1302.8844 / t(i,j)
          c1 = 23.832241 - 5.02808 * alog10(t(i,j))
          term = c1-1.3816e-7*10.**p1+8.1328e-3*10.**p2-2949.076/t(i,j)
          es(i,j) = 10. ** term
        enddo !i
        enddo !j
c
        return
        end
