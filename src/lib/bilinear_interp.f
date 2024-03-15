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
        subroutine bilinear_interp (n1,m1,g1,n2,m2,g2,dn,dm)

c       author: d. birkenheuer 18 oct 90
c       checked vectorizing      8 november 1990    (d. birkenheuer)

c       routine is designed to interpolate from a less dense grid (g1)
c       to a grid of higher density (g2)  where g1 is nested in g2.

c       for current purposes, g1 is lma array (8,8), n1=8, m1=8
c       g2 is laps array (57,57) n2=57,m2=57, and dn=dm=8 the spacing
c       between the nested gridpoints in the g2 array.

c       dn or dm can be computed by:
c       (n2-1)/(n1-1) = dn.....and similarly for dm

        implicit none

        integer n1,n2,m1,m2,dn,dm
        real g1(n1,m1),g2(n2,m2)

        real factori, factorj, fraci, fracj, sum
        integer inti,intj,i,j


        do j = 1,m2
        factorj = float(j-1)/float(dm) +1.
        intj = int(factorj)
        fracj = factorj-intj

        do i = 1,n2
        factori = float(i-1)/float(dn) +1.
        inti = int(factori)
        fraci = factori-inti

        sum = 0.

        g2(i,j) = g1(inti,intj) * (1.-fraci) * (1.-fracj)
        sum = sum + (1.-fraci) * (1.-fracj)

        if (inti+1 .le. n1) then
        g2(i,j) = g2(i,j) + g1(inti+1,intj) * (fraci) * (1.-fracj)
        sum = sum + (fraci) * (1.-fracj)
        endif


        if (intj+1 .le. m1) then
        g2(i,j) = g2(i,j) + g1(inti,intj+1) * (fracj) * (1.-fraci)
        sum = sum + (fracj) * (1.-fraci)
        endif

        if (inti+1 .le. n1  .and. intj+1  .le. m1)then
        g2(i,j) = g2(i,j) + g1(inti+1,intj+1) * (fraci) * (fracj)
        sum = sum + (fraci) * (fracj)
        endif

        g2(i,j) = g2(i,j)/sum

        enddo
        enddo


        return
        end

