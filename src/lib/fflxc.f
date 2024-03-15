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

        subroutine fflxc(ni,nj,m,scale,uh,vh,field,flxcnv,lat,lon 
     1  ,flu,flv,sigma,r_missing_data)

cdoc    calculates convergence given a wind field rotated to the map projection

!       original version s. albers                                      1986
!       made more general                                               1989
!       made even more general (adjustable dims)                        1991
!       added new map projections for sigma calc                        1997
!       input m is changed, hopefully balances changes in get_sigma     1998

        include 'trigd.inc'
        real m
        dimension field(ni,nj),flu(ni,nj),flv(ni,nj),sigma(ni,nj)
        real uh(ni,nj),vh(ni,nj),flxcnv(ni,nj),lat(ni,nj),lon(ni,nj)
        data  rpd/.0174532925/

c       write(6,18)
c18     format('        calculating flux convergence field')

c
c  calculate field flux and initialize arrays
c  units for field flux convergence are 1/scale * s**-1
c      coeff=-1/
c      (scale/(2*(delta*.0254mpin/m) meters/gridpoint) * .5147818mspkt)
!       coeff=-1./(scale/(2.*(delta*.0254/m))*.5147818)
        coeff=-1./(scale/(2./m))
c       write(6,*)' coeff = ',coeff
        nim1=ni-1
        nim2=ni-2
        njm1=nj-1
        njm2=nj-2
c
c       calculate sigma
c       this value represents  (m*rho(earth)*(1.+sin(phi0))/mpin)**2

        do 100 i=1,ni
        do 100 j=1,nj
!           sigma(i,j) = (1.+sind(phi0)) / (1.+ sind(lat(i,j))) 
!           sigma(i,j)=1.0

            call get_sigma(lat(i,j),lon(i,j),sigma(i,j),istatus)
            if(istatus .ne. 1)stop
c
            fl=field(i,j)/(coeff*sigma(i,j))
            flu(i,j)=fl* uh(i,j)
            flv(i,j)=fl* vh(i,j)
            flxcnv(i,j)=0.
100     continue

c
c
        dfludx=-3.*flu(   1,   1)+4.*flu(   2,   1)-flu(   3,   1)
        dflvdy=-3.*flv(   1,   1)+4.*flv(   1,   2)-flv(   1,   3)
 1000   flxcnv(1,1)=(dfludx+dflvdy)*sigma(1,1)**2
c
        do 2000 j=2,njm1
        dfludx=-3.*flu(   1,   j)+4.*flu(   2,   j)-flu(   3,   j)
        dflvdy=    flv(   1, j+1)  - flv(   1, j-1)
 2000   flxcnv(1,j)=(dfludx+dflvdy)*sigma(1,j)**2
c
        dfludx=-3.*flu(   1,  nj)+4.*flu(   2,  nj)-flu(   3,  nj)
        dflvdy= 3.*flv(   1,  nj)-4.*flv(   1,njm1)+flv(   1,njm2)
 3000   flxcnv(1,nj)=(dfludx+dflvdy)*sigma(1,nj)**2
c
        do 4000 i=2,nim1
        dfludx=    flu( i+1,   1)  - flu( i-1,   1)
        dflvdy=-3.*flv(   i,   1)+4.*flv(   i,   2)-flv(   i,   3)
 4000   flxcnv(i,1)=(dfludx+dflvdy)*sigma(i,1)**2
c
        do 5000 i=2,nim1
        ip1=i+1
        im1=i-1

        do 5000 j=2,njm1
        dfludx=    flu( ip1,   j)  - flu( im1,   j)
        dflvdy=    flv(   i, j+1)  - flv(   i, j-1)
        flxcnv(i,j)=(dfludx+dflvdy)*sigma(i,j)**2
c       call outptf(i,j,flu,flv,dfludx,dflvdy,flxcnv,uh,vh,sigma
c      .      ,nim2,njm2,field,scale,coeff)
 5000   continue
c
        do 6000 i=2,nim1
        dfludx=    flu( i+1,  nj)  - flu( i-1,  nj)
        dflvdy= 3.*flv(   i,  nj)-4.*flv(   i,njm1)+flv(   i,njm2)
 6000   flxcnv(i,nj)=(dfludx+dflvdy)*sigma(i,nj)**2
c
        dfludx= 3.*flu(  ni,   1)-4.*flu(nim1,   1)+flu(nim2,   1)
        dflvdy=-3.*flv(  ni,   1)+4.*flv(  ni,   2)-flv(  ni,   3)
 7000   flxcnv(ni,1)=(dfludx+dflvdy)*sigma(ni,1)**2
c
        do 8000 j=2,njm1
        dfludx= 3.*flu(  ni,   j)-4.*flu(nim1,   j)+flu(nim2,   j)
        dflvdy=    flv(  ni, j+1)  - flv(  ni, j-1)
 8000   flxcnv(ni,j)=(dfludx+dflvdy)*sigma(ni,j)**2
c
        dfludx= 3.*flu(  ni,  nj)-4.*flu(nim1,  nj)+flu(nim2,  nj)
        dflvdy= 3.*flv(  ni,  nj)-4.*flv(  ni,njm1)+flv(  ni,njm2)
 9000   flxcnv(ni,nj)=(dfludx+dflvdy)*sigma(ni,nj)**2
        return
        end

        subroutine get_sigma(rlat,rlon,sigma,istatus)

!       steve albers

cdoc    sigma (map factor) is defined to be one when we are located in the 
cdoc    intersection of the projection plane and the earth's surface. it varies
cdoc    from unity elsewhere.

!       equations from principles of meteorological analysis, walter saucier
!       pages 32,33

        include 'trigd.inc'
        real n
        character*6 c6_maproj

        integer init
        save init
        data init/0/

        save slat1,slat2,slon,c6_maproj

        if(init .eq. 0)then
            call get_standard_latitudes(slat1,slat2,istatus)
            if(istatus .ne. 1)then
                write(6,*)'get_sigma: bad istatus'
                return
            endif

            call get_standard_longitude(slon,istatus)
            if(istatus .ne. 1)then
                write(6,*)'get_sigma: bad istatus'
                return
            endif

            call get_c6_maproj(c6_maproj,istatus)
            if(istatus .ne. 1)then
                write(6,*)'get_sigma: bad istatus'
                return
            endif

            init = 1

        endif

        colat0 = 90. - slat1
        colat1 = 90. - slat1
        colat2 = 90. - slat2
        colat  = 90. - rlat

        if(c6_maproj .eq. 'plrstr')then ! polar stereo

!           rotate to arbitrary pole
            polat = slat2
            polon = slon
            rlat_ge = rlat
            rlon_ge = rlon
            call getops(rlat_ps,rlon_ps,rlat_ge,rlon_ge,polat,polon)       
!                          o       o        i      i      i     i

            call get_grid_spacing(grid_spacing_m,istatus)
            if(istatus .ne. 1)then
                write(6,*)'get_sigma: bad istatus'
                return
            endif

            call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                       ,grid_spacing_proj_m)

!           phi0 = slat1
            phi  = rlat_ps

!           check to see if you are at the opposite pole
            if(phi .eq. -90.)then      
                sigma = 0.
                istatus = 0
                write(6,*)'get_sigma: sigma undefined at opposite pole'       
                return
            endif

            sigma = (1. + sind(phi0)) / (1. + sind(phi))               ! eq. 13

        elseif(c6_maproj .eq. 'lambrt')then ! lambert

!           check to see if you are at a pole
            if(abs(rlat) .eq. 90.)then      
                sigma = 0.
                istatus = 0
                write(6,*)'get_sigma: sigma undefined at the pole'
                return
            endif

            if(slat1 .eq. slat2)then
                n = cosd(colat0)
                arg   =  tand(colat/2.) / tand(colat0/2.)
                sigma = (sind(colat0)   / sind(colat)   ) * arg**n     ! eq. 3

            else

                n = alog(sind(colat1)   /sind(colat2)   ) /            ! eq. 9
     1              alog(tand(colat1/2.)/tand(colat2/2.))

                arg   =  tand(colat/2.) / tand(colat1/2.)
                sigma = (sind(colat1)   / sind(colat)   ) * arg**n     ! eq. 10

            endif

        elseif(c6_maproj .eq. 'merctr')then ! mercator

!           check to see if you are at a pole
            if(abs(rlat) .eq. 90.)then      
                sigma = 0.
                istatus = 0
                write(6,*)'get_sigma: sigma undefined at the pole'
                return
            endif

            sigma = sind(colat1) / sind(colat)                         ! eq. 11

        else
            write(6,*)' invalid map projection in get_sigma: ',c6_maproj
            istatus = 0
            return

        endif

        istatus = 1
        return
        end
