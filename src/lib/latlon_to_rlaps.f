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

        subroutine latlon_to_rlapsgrid(rlat,rlon,lat,lon,ni,nj,ri,rj
     1                                ,istatus)

!       1991            steve albers
!       1994            steve albers - partially added lambert option
!       1997            steve albers - added both kinds of lambert projections
!                                    - as well as mercator projection
!       1997            steve albers - added local stereographic
!       2002            dan birkenheuer - added input test for valid lats

cdoc    this routine assumes a polar stereographic, lambert conformal,
cdoc    or mercator projection.
cdoc
cdoc    values returned are on the laps grid where i ranges from 1,ni and
cdoc    j ranges from 1,nj. istatus is set to 1 if a valid i,j was determined
cdoc    to be in the laps domain, except that a buffer of 0.5 grid points 
cdoc    around the perimenter is allowed. if the point is farther outside
cdoc    the domain, values of i,j will still be returned, while istatus is set
cdoc    to 0. i increased from grid west to grid east, and j increases from
cdoc    grid south to grid north.

        real rlat                         ! input lat
        real rlon                         ! input lon
        real lat(ni,nj),lon(ni,nj)        ! input (arrays of lat/lon)
        integer ni,nj                     ! input (laps dimensions)
        real ri,rj                        ! output (i,j on laps grid)
        integer istatus                   ! input / output

        save init,umin,umax,vmin,vmax,niprev,njprev
        data init/0/, niprev/0/, njprev/0/

        include 'grid_fname.cmn'

        if (abs(rlat) > 90.000) then
           if(istatus .ne. 100)then
              write(6,*) 'rejecting invalid latitude ',rlat
           endif
           istatus = -1
           return
        endif

        if(init.ne.nest .or. ni.ne.niprev .or. nj.ne.njprev)then
            call latlon_to_uv(lat(1,1),lon(1,1),umin,vmin,istatus)
            call latlon_to_uv(lat(ni,nj),lon(ni,nj),umax,vmax,istatus)

            write(6,101)umin,umax,vmin,vmax
101         format(1x,' initializing latlon_to_rlapsgrid',4f10.5)
            init = nest
            niprev = ni
            njprev = nj
        endif

        uscale = (umax - umin) / (float(ni) - 1.)
        vscale = (vmax - vmin) / (float(nj) - 1.)

        u0 = umin - uscale
        v0 = vmin - vscale

!       compute ulaps and vlaps

        call latlon_to_uv(rlat,rlon,ulaps,vlaps,istatus)

        if(umax .le. umin .or. vmax .le. vmin)then
            write(6,*)
     1      ' severe error: uvminmax are off in latlon_to_rlapsgrid'
            write(6,*)umin,umax,vmin,vmax
            stop
        endif

        if(uscale .eq. 0. .or. vscale .eq. 0.)then
            write(6,*)
     1      ' severe error: (u|v)scale = 0 in latlon_to_rlapsgrid'
            stop
        else
            ri = (ulaps - u0) / uscale
            rj = (vlaps - v0) / vscale
        endif

!       set status if location of point rounded off is on the laps grid
        if(nint(ri) .ge. 1 .and. nint(ri) .le. ni .and.
     1     nint(rj) .ge. 1 .and. nint(rj) .le. nj         )then
            istatus = 1
        else
            istatus = 0
        endif

        return
        end

        subroutine rlapsgrid_to_latlon(ri,rj,lat,lon,ni,nj,rlat,rlon
     1                                ,istatus)

!       1997            steve albers 

cdoc    this routine assumes a polar stereographic, lambert conformal,
cdoc    or mercator projection.

        real ri,rj                        ! input (i,j on laps grid)
        real lat(ni,nj),lon(ni,nj)        ! input (arrays of lat/lon)
        integer ni,nj                       ! input (laps dimensions)
        real rlat                         ! output lat
        real rlon                         ! output lon
        integer istatus                     ! output

        save init,umin,umax,vmin,vmax
        data init/0/

        include 'grid_fname.cmn'

        if(init .ne. nest)then
            call latlon_to_uv(lat(1,1),lon(1,1),umin,vmin,istatus)
            call latlon_to_uv(lat(ni,nj),lon(ni,nj),umax,vmax,istatus)

            write(6,101)umin,umax,vmin,vmax
101         format(1x,' initializing rlapsgrid_to_latlon',4f10.5)
            init = nest 
        endif

        uscale = (umax - umin) / (float(ni) - 1.)
        vscale = (vmax - vmin) / (float(nj) - 1.)

        u0 = umin - uscale
        v0 = vmin - vscale

!       compute lat,lon
        ulaps = u0 + uscale * ri
        vlaps = v0 + vscale * rj

        call uv_to_latlon(ulaps,vlaps,rlat,rlon,istatus)

        istatus = 1

        return
        end

        
        subroutine latlon_to_uv(rlat,rlon,u,v,istatus)

!       1997            steve albers 

cdoc    this routine assumes a polar stereographic, lambert conformal,
cdoc    or mercator projection.

        use mem_namelist, only: c6_maproj
     1                  ,slat1=>standard_latitude
     1                  ,slat2=>standard_latitude2
     1                  ,polat=>standard_latitude2
     1                  ,slon=>standard_longitude

        include 'trigd.inc'

        include 'grid_fname.cmn'

        if(c6_maproj .eq. 'plrstr')then ! polar stereo
            call latlon_to_uv_ps(rlat,rlon,slat1,polat,slon,u,v)

        elseif(c6_maproj .eq. 'lambrt')then ! lambert
            if(abs(slat2).eq.90. .and. slat2.ne.slat1)then
               print*,'error: lambert slat2 = 90.'
               istatus=1
               return
            endif
c           slat1 = standard_latitude
c           slat2 = standard_latitude2
c           slon = standard_longitude

            call latlon_to_uv_lc(rlat,rlon,slat1,slat2,slon,u,v)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif
c           slat1  = standard_latitude
c           cenlon = grid_cen_lon_cmn

            call latlon_to_uv_mc(rlat,rlon,slat1,cenlon,u,v)

        elseif(c6_maproj .eq. 'latlon')then ! latlon (cylindrical equidistant)
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif

            call latlon_to_uv_ll(rlat,rlon,cenlon,u,v)

        else
            write(6,*)'error in latlon_to_uv: unrecognized projection '
     1                ,c6_maproj       
            stop

        endif

        istatus = 1

        return
        end


        subroutine latlon_to_uv_ps(rlat_in,rlon_in,slat,polat,slon,u,v)

        include 'trigd.inc'

!       integer iwrite /0/
!       save iwrite

        if(abs(polat) .eq. 90.)then ! pole at n/s geographic pole
            if(.true.)then
                polon = slon
                call getops(rlat,rlon,rlat_in,rlon_in,polat,polon)
                rlon = rlon - 270.  ! counteract rotation done in getops

            else ! .false. (older simple method)
                rlat = rlat_in
                rlon = rlon_in

            endif

            b = rlon - slon         ! rotate relative to standard longitude
            s = polat / abs(polat)

        else                        ! local stereographic
            polon = slon
            call getops(rlat,rlon,rlat_in,rlon_in,polat,polon)
            b = rlon - 270.         ! rlon has zero angle pointing east
                                    ! b has zero angle pointing south
            s = 1.0

!           iwrite = iwrite + 1
!           if(rlon_in .eq. 53.99)then
!               write(6,1)rlat_in,rlon_in,b,rlat,rlon
!1              format('rlat_in | rlon_in | b | rlat | rlon',5f9.4)             
!           endif

        endif

        a=90.-rlat
        r = tand(a/2.)      ! consistent with haltiner & williams 1-21

!       b = angle measured counterclockwise from -v axis (zero angle south)
        u =   r * sind(b)
        v = (-r * cosd(b)) * s

!       if(rlon_in .eq. 53.99)then
!           write(6,2)rlat_in,rlon_in,b,rlat,rlon,u,v
!2          format('rlat_in | rlon_in | b | rlat | rlon | u | v'
!    1            ,5f10.4,2f10.6) 
!       endif

        return
        end

        subroutine latlon_to_uv_lc(rlat,rlon,slat1,slat2,slon,u,v)

        include 'trigd.inc'
        real n

!       difference between two angles, result is between -180. and +180.
        angdif(x,y)=mod(x-y+540.,360.)-180.

        call lambert_parms(slat1,slat2,n,s,rconst)

        r = (tand(45.-s*rlat/2.))**n
        u =    r*sind(n*angdif(rlon,slon))
        v = -s*r*cosd(n*angdif(rlon,slon))

        return
        end

        subroutine latlon_to_uv_mc(rlat,rlon,slat,cenlon,u,v)

        include 'trigd.inc'
        real pi, rpd

        parameter (pi=3.1415926535897932)
        parameter (rpd=pi/180.)

!       difference between two angles, result is between -180. and +180.
        angdif(x,y)=mod(x-y+540.,360.)-180.
 
        a = 90.-rlat
        b = cosd(slat)

        u = angdif(rlon,cenlon) * rpd * b
        v = alog(1./tand(a/2.))       * b

        return
        end

        subroutine latlon_to_uv_ll(rlat,rlon,cenlon,u,v)

        include 'trigd.inc'
        real pi, rpd

        parameter (pi=3.1415926535897932)
        parameter (rpd=pi/180.)

!       difference between two angles, result is between -180. and +180.
        angdif(x,y)=mod(x-y+540.,360.)-180.
 
        u = angdif(rlon,cenlon) * rpd 
        v = rlat                * rpd

        return
        end


        
        subroutine uv_to_latlon(u,v,rlat,rlon,istatus)

        use mem_namelist, only: c6_maproj
     1                  ,slat1=>standard_latitude
     1                  ,slat2=>standard_latitude2
     1                  ,polat=>standard_latitude2
     1                  ,slon=>standard_longitude

!       1997            steve albers 

!       this routine assumes a polar stereographic, lambert conformal,
!       or mercator projection.

        include 'trigd.inc'

        include 'grid_fname.cmn'

        if(c6_maproj .eq. 'plrstr')then ! polar stereo
c           slat1 = standard_latitude
c           polat = standard_latitude2
c           slon = standard_longitude

            call uv_to_latlon_ps(u,v,slat1,polat,slon,rlat,rlon)

        elseif(c6_maproj .eq. 'lambrt')then ! lambert
c           slat1 = standard_latitude
c           slat2 = standard_latitude2
c           slon = standard_longitude

            call uv_to_latlon_lc(u,v,slat1,slat2,slon,rlat,rlon)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif
c           slat1  = standard_latitude
c           cenlon = grid_cen_lon_cmn

            call uv_to_latlon_mc(u,v,slat1,cenlon,rlat,rlon)

        elseif(c6_maproj .eq. 'latlon')then ! latlon (cylindrical equidistant)
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif

            call uv_to_latlon_ll(u,v,cenlon,rlat,rlon)

        else
            write(6,*)'error in uv_to_latlon: unrecognized projection '
     1                ,c6_maproj
            stop

        endif

        istatus = 1

        return
        end

        subroutine uv_to_latlon_ps(u,v,slat,polat,slon
     1                                         ,rlat_out,rlon_out)

        include 'trigd.inc'

        if(abs(polat) .eq. 90.)then
            s = polat / abs(polat)
        else
            s = 1.0
        endif

        r=sqrt(u**2+v**2)

        if (r .eq. 0) then
            rlat=90.
            rlon=0.

        else                           
            a=2.* atand(r)               ! from haltiner & williams 1-21
            rlat=90.- a
            rlon = atan2d(s*v,u)
            rlon = rlon + 90.
        endif

        if(.true.)then ! rotate considering where the projection pole is
            polon = slon
!           this routine will rotate the longitude even if polat = +90.
            call pstoge(rlat,rlon,rlat_out,rlon_out,polat,polon)
        else
            rlat_out = rlat
            rlon_out = rlon + slon
        endif


        rlon_out = amod(rlon_out+540.,360.) - 180. ! convert to -180/+180 range


        return
        end

        subroutine uv_to_latlon_lc(u,v,slat1,slat2,slon,rlat,rlon)

        include 'trigd.inc'
        real n

        call lambert_parms(slat1,slat2,n,s,rconst)

!       rlon=slon + atand(-s*u/v) /n
!       rlat=(90.- 2.*atand((-  v/cosd(n*(rlon-slon)))**(1./n)))/s      

        angle  = atan2d(u,-s*v)
        rlat = (90.- 2.*atand((-s*v/cosd(angle))**(1./n))) / s      
        rlon = slon + angle / n

        rlon = mod(rlon+540.,360.) - 180.          ! convert to -180/+180 range

        return
        end

        subroutine uv_to_latlon_mc(u,v,slat,cenlon,rlat,rlon)

        include 'trigd.inc'
        parameter (pi=3.1415926535897932)
        parameter (rpd=pi/180.)

        b = cosd(slat)

        rlat_abs = 90. - atand(exp(-abs(v)/b)) * 2.

        if(v .gt. 0)then
            rlat =  rlat_abs
        else
            rlat = -rlat_abs
        endif

        rlon = u/b/rpd + cenlon
        rlon = mod(rlon+540.,360.) - 180.          ! convert to -180/+180 range

        return
        end

        subroutine uv_to_latlon_ll(u,v,cenlon,rlat,rlon)

        include 'trigd.inc'
        parameter (pi=3.1415926535897932)
        parameter (rpd=pi/180.)

        rlat = v/rpd

        rlon = u/rpd + cenlon
        rlon = mod(rlon+540.,360.) - 180.          ! convert to -180/+180 range

        return
        end


        function projrot_latlon(rlat,rlon,istatus)
        include 'trigd.inc'

!       projrot is the clockwise angle of a vector pointing true north relative
!       to a vector pointing to grid north. this is also the clockwise rotation
!       of a longitude line plotted on the model grid. units are degrees.

!       for azimuth or wind direction angles:

!       projrot = grid north - true north
!       true north = grid north - projrot
!       grid north = true north + projrot
       
!       added 8/4/2000 to make sure these are declared even if not passed in
        integer istatus
        real rlat, rlon

        real n

        save init
        data init/0/

        character*6  c6_maproj

        include 'grid_fname.cmn'

!       difference between two angles, result is between -180. and +180.
        angdif(x,y)=mod(x-y+540.,360.)-180.

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus.ne.1)then
           print*,'error returned from get_c6_maproj'
           return
        endif

        if(c6_maproj .eq. 'plrstr')then ! polar stereographic

            call get_standard_longitude(polon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_longitude'
               return
            endif
            call get_standard_latitudes(stdlat,polat,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_latitudes'
               return
            endif

            call get_grid_center(grid_cen_lat,grid_cen_lon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif


            if(polat .eq. +90.)then
                projrot_laps = polon - rlon

            elseif(polat .eq. -90.)then
                projrot_laps = rlon - polon 

            else ! abs(polat) .ne. 90.
!               if(grid_cen_lat .eq. polat .and. 
!    1             grid_cen_lon .eq. polon)then ! grid centered on proj pole
                if(.true.)then

                    if(init .eq. 0)then
                        write(6,*)
     1                   ' note: local stereographic projection.'
                        write(6,*)
     1                   ' using approximation for "projrot_latlon",'
     1                  ,' accurate calculation not yet in place.'
                        init = 1
                    endif

                    rn = cosd(90.-polat)
                    projrot_laps = rn * angdif(polon,rlon)      

                elseif(.true.)then
                    if(init .eq. 0)then
                        write(6,*)' error in projrot_latlon: '
                        write(6,*)' this type of local'
     1                  ,' stereographic projection not yet supported.'
                        write(6,*)' grid should be centered on'
     1                  ,' projection pole.'
                        init = 1
                    endif

                    projrot_laps = 0.
         
                else ! .false.
!                   find dx/lat and dy/lat, then determine projrot_laps

                endif

            endif ! polat

        elseif(c6_maproj .eq. 'lambrt')then ! lambert conformal

            call get_standard_latitudes(slat1,slat2,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_latitudes'
               return
            endif
            call get_standard_longitude(slon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_longitude'
               return
            endif

            call lambert_parms(slat1,slat2,n,s,rconst)

            call get_standard_longitude(stdlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_c6_maproj'
               return
            endif

            projrot_laps = n * s * angdif(stdlon,rlon)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            projrot_laps = 0.

        elseif(c6_maproj .eq. 'latlon')then ! latlon
            projrot_laps = 0.

        else
            write(6,*)
     1            'error in projrot_latlon: unrecognized projection '      
     1               ,c6_maproj       
            stop

        endif

        projrot_latlon = projrot_laps

        istatus = 1 
        return
        end

        subroutine projrot_latlon_2d(rlat,rlon,ni,nj,projrot_laps
     1                                              ,istatus)
        include 'trigd.inc'

cdoc    1997 steve albers    calculate map projection rotation, this is the
cdoc                         angle between the y-axis (grid north) and
cdoc                         true north. units are degrees.
!
!                            projrot_laps = (true north value of wind direction
!                                          - grid north value of wind direction)
       
!       added 8/4/2000 to make sure these are declared even if not passed in
        integer istatus
        real rlat(ni,nj), rlon(ni,nj), projrot_laps(ni,nj)

        real n

        save init
        data init/0/

        character*6  c6_maproj

        include 'grid_fname.cmn'

!       difference between two angles, result is between -180. and +180.
        angdif(x,y)=mod(x-y+540.,360.)-180.

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus.ne.1)then
           print*,'error returned from get_c6_maproj'
           return
        endif

        if(c6_maproj .eq. 'plrstr')then ! polar stereographic

            call get_standard_longitude(polon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_longitude'
               return
            endif
            call get_standard_latitudes(stdlat,polat,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_latitudes'
               return
            endif

            call get_grid_center(grid_cen_lat,grid_cen_lon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif


            if(polat .eq. +90.)then
                projrot_laps(:,:) = polon - rlon(:,:)

            elseif(polat .eq. -90.)then
                projrot_laps(:,:) = rlon(:,:) - polon 

            else ! abs(polat) .ne. 90.
!               if(grid_cen_lat .eq. polat .and. 
!    1             grid_cen_lon .eq. polon)then ! grid centered on proj pole
                if(.true.)then

                    if(init .eq. 0)then
                        write(6,*)
     1                   ' note: local stereographic projection.'
                        write(6,*)
     1                   ' using approximation for "projrot_latlon_2d",'
     1                  ,' accurate calculation not yet in place.'
                        init = 1
                    endif

                    rn = cosd(90.-polat)
                    do j = 1,nj
                    do i = 1,ni
                        projrot_laps(i,j) = rn * angdif(polon,rlon(i,j))
                    enddo ! i
                    enddo ! j

                elseif(.true.)then
                    if(init .eq. 0)then
                        write(6,*)' error in projrot_latlon_2d: '
                        write(6,*)' this type of local'
     1                  ,' stereographic projection not yet supported.'
                        write(6,*)' grid should be centered on'
     1                  ,' projection pole.'
                        init = 1
                    endif

                    projrot_laps = 0.
         
                else ! .false.
!                   find dx/lat and dy/lat, then determine projrot_laps

                endif

            endif ! polat

        elseif(c6_maproj .eq. 'lambrt')then ! lambert conformal

            call get_standard_latitudes(slat1,slat2,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_latitudes'
               return
            endif
            call get_standard_longitude(slon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_longitude'
               return
            endif

            call lambert_parms(slat1,slat2,n,s,rconst)

            call get_standard_longitude(stdlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_c6_maproj'
               return
            endif

            do j = 1,nj
            do i = 1,ni
                projrot_laps(i,j) = n * s * angdif(stdlon,rlon(i,j))
            enddo ! i
            enddo ! j

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            projrot_laps = 0.

        elseif(c6_maproj .eq. 'latlon')then ! latlon
            projrot_laps = 0.

        else
            write(6,*)
     1         'error in projrot_latlon_2d: unrecognized projection '
     1         ,c6_maproj       
            stop

        endif

        istatus = 1 
        return
        end



        function projrot_laps(rlon)
        include 'trigd.inc'

cdoc    this routine is being phased out. please try to use 'projrot_latlon'.

!       1997 steve albers    calculate map projection rotation, this is the
!                            angle between the y-axis (grid north) and
!                            true north. units are degrees.
!
!                            projrot_laps = (true north value of wind direction
!                                          - grid north value of wind direction)
       
!       added 8/4/2000 to make sure these are declared even if not passed in
        integer istatus
        real rlat, rlon

        real n

        save init
        data init/0/

        character*6  c6_maproj

        include 'grid_fname.cmn'

!       difference between two angles, result is between -180. and +180.
        angdif(x,y)=mod(x-y+540.,360.)-180.

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus.ne.1)then
           print*,'error returned from get_c6_maproj'
           return
        endif

        if(c6_maproj .eq. 'plrstr')then ! polar stereographic

            call get_standard_longitude(polon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_longitude'
               return
            endif
            call get_standard_latitudes(stdlat,polat,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_latitudes'
               return
            endif

            call get_grid_center(grid_cen_lat,grid_cen_lon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif


            if(polat .eq. +90.)then
                projrot_laps = polon - rlon

            elseif(polat .eq. -90.)then
                projrot_laps = rlon - polon 

            else ! abs(polat) .ne. 90.
!               if(grid_cen_lat .eq. polat .and. 
!    1             grid_cen_lon .eq. polon)then ! grid centered on proj pole
                if(.true.)then

                    if(init .eq. 0)then
                        write(6,*)
     1                   ' note: local stereographic projection.'
                        write(6,*)
     1                   ' using approximation for "projrot_laps",'
     1                  ,' accurate calculation not yet in place.'
                        init = 1
                    endif

                    rn = cosd(90.-polat)
                    projrot_laps = rn * angdif(polon,rlon)      

                elseif(.true.)then
                    if(init .eq. 0)then
                        write(6,*)' error in projrot_laps: '
                        write(6,*)' this type of local'
     1                  ,' stereographic projection not yet supported.'
                        write(6,*)' grid should be centered on'
     1                  ,' projection pole.'
                        init = 1
                    endif

                    projrot_laps = 0.
         
                else ! .false.
!                   find dx/lat and dy/lat, then determine projrot_laps

                endif

            endif ! polat

        elseif(c6_maproj .eq. 'lambrt')then ! lambert conformal

            call get_standard_latitudes(slat1,slat2,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_latitudes'
               return
            endif
            call get_standard_longitude(slon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_standard_longitude'
               return
            endif

            call lambert_parms(slat1,slat2,n,s,rconst)

            call get_standard_longitude(stdlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_c6_maproj'
               return
            endif

            projrot_laps = n * s * angdif(stdlon,rlon)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            projrot_laps = 0.

        elseif(c6_maproj .eq. 'latlon')then ! latlon
            projrot_laps = 0.

        else
            write(6,*)'error in projrot_laps: unrecognized projection '
     1               ,c6_maproj       
            stop

        endif

        projrot_latlon = projrot_laps

!       istatus = 1 ! not yet usable because of the entry
        return
        end


      subroutine check_domain(lat,lon,ni,nj,grid_spacing_m,intvl
     1                                                       ,istatus)

cdoc  this routine checks whether the lat/lon grid is consistent with
cdoc  map projection parameters as processed by latlon_to_rlapsgrid,
cdoc  and rlapsgrid_to_latlon. the grid size is also checked.
cdoc  this is a good sanity check of the netcdf static file, namelist files,
cdoc  as well as various grid conversion routines.

!     1997 steve albers

      include 'trigd.inc'
      real pi, rpd
      parameter (pi=3.1415926535897932)
      parameter (rpd=pi/180.)

      real lat(ni,nj),lon(ni,nj)

      character*6  c6_maproj

      istatus = 1
      tolerance_m = 1000.

      diff_grid_max = 0.

      write(6,*)
      write(6,*)' subroutine check_domain: checking latlon_to_rlapsgrid'

      call get_c6_maproj(c6_maproj,istatus)
      if(istatus.ne.1)then
         print*,'error from get_c6_maproj'
         return
      endif

      do i = 1,ni,intvl
      do j = 1,nj,intvl
          call latlon_to_rlapsgrid(lat(i,j),lon(i,j),lat,lon,ni,nj
     1                                              ,ri,rj,istat)

          if(istat .ne. 1)then
              write(6,*)' bad status from latlon_to_rlapsgrid'
              istatus = 0
              return
          endif

          diff_gridi = ri - float(i)
          diff_gridj = rj - float(j)
          diff_grid = sqrt(diff_gridi**2 + diff_gridj**2)

          if(diff_grid .gt. diff_grid_max)then
              diff_grid_max = diff_grid
              idmax = i
              jdmax = j
          endif

          diff_grid_max_m = diff_grid_max * grid_spacing_m

      enddo
      enddo

      write(6,*)' check_domain: max_diff (gridpoints) = ',diff_grid_max
     1         ,' at i/j',idmax,jdmax
      write(6,*)' check_domain: max_diff (approx m)   = '
     1                                                 ,diff_grid_max_m      

      if(diff_grid_max_m .gt. tolerance_m)then
          write(6,*)' warning: exceeded tolerance in check_domain'
     1               ,tolerance_m
          istatus = 0
      endif

!...........................................................................

      diff_ll_max = 0.

      write(6,*)' checking rlapsgrid_to_latlon'

      do i = 1,ni,intvl
      do j = 1,nj,intvl
          ri = i
          rj = j
          call rlapsgrid_to_latlon(ri,rj,lat,lon,ni,nj
     1                                  ,rlat,rlon,istat)

          if(istat .ne. 1)then
              write(6,*)' bad status from rlapsgrid_to_latlon'
              istatus = 0
              return
          endif

          diff_lli =  rlat - lat(i,j)
          diff_llj = (rlon - lon(i,j)) * cosd(lat(i,j))
          diff_ll = sqrt(diff_lli**2 + diff_llj**2)

          if(diff_ll .gt. diff_ll_max)then
              diff_ll_max = diff_ll
              idmax = i
              jdmax = j
          endif

          diff_ll_max_m = diff_ll_max * 110000. ! meters per degree

      enddo
      enddo

      write(6,*)' check_domain: max_diff (degrees) = ',diff_ll_max
     1         ,' at i/j',idmax,jdmax
      write(6,*)' check_domain: max_diff (approx m)   = ',diff_ll_max_m

      if(diff_ll_max_m .gt. tolerance_m)then
          write(6,*)' warning: exceeded tolerance in check_domain'
     1               ,tolerance_m
          istatus = 0
      endif

!...........................................................................

      call get_earth_radius(erad,istatus)
      if(istatus .ne. 1)then
          write(6,*)
     1    ' error calling get_earth_radius from check_domain'       
          return
      endif

      icen = ni/2+1
      jcen = nj/2+1

      diff_lat =  lat(icen,jcen+1) - lat(icen,jcen-1)
      diff_lon = (lon(icen,jcen+1) - lon(icen,jcen-1)) 
     1                        * cosd(lat(icen,jcen))

      dist = sqrt((diff_lat)**2 + (diff_lon)**2) / 2. * rpd * erad   

      if(abs(lat(icen,jcen)) .lt. 89.)then ! should be reasonably accurate
          write(6,*)
     1   ' measured grid spacing on earths surface at domain center is:'       
     1     ,dist       
      endif

!...........................................................................

      if(c6_maproj .ne. 'latlon')then
          call get_grid_spacing_actual(lat(icen,jcen),lon(icen,jcen)
     1                                      ,dist_calc,istatus)
          if(istatus .ne. 1)then
              write(6,*)
     1        ' error calling get_grid_spacing_actual from check_domain'       
              return
          endif

          write(6,*)
     1 ' calculated grid spacing on earths surface at domain center is:'       
     1     ,dist_calc

      endif

!...........................................................................

      call latlon_to_xy(lat(1,1),lon(1,1),erad,x1,y1)
      call latlon_to_xy(lat(1,2),lon(1,2),erad,x2,y2)

      dist = sqrt((x2-x1)**2 + (y2-y1)**2)

      write(6,*)
     1 ' grid spacing on projection plane using "latlon_to_xy" is:'      
     1 ,dist
      
!...........................................................................

      if(c6_maproj .eq. 'lambrt')then
          call get_standard_latitudes(std_lat1,std_lat2,istatus)
          if(istatus .ne. 1)then
              write(6,*)
     1        ' error calling get_standard_latitudes from check_domain'
              return
          endif

          call lambert_parms(std_lat1,std_lat2,constn,consts,constr)
          single_lat = (acosd(constn) - 90.)/(-consts)

          call get_grid_spacing_actual(single_lat,lon(icen,jcen)
     1                                      ,dist_calc,istatus)
          if(istatus .ne. 1)then
              write(6,*)
     1        ' error calling get_grid_spacing_actual from check_domain'       
              return
          endif

          write(6,*)' calculated grid spacing on '
     1      ,'earths surface at single lambert latitude:'       
     1     ,single_lat,dist_calc

      endif

!...........................................................................

      write(6,*)

      return
      end

      subroutine lambert_parms(slat1,slat2,n_out,s_out,rconst_out)

      include 'trigd.inc'
      real n,n_out

!     we only have to do the calculations once since the inputs are constants
      data init/0/
      save init,n,s,rconst 

      if(init .eq. 0)then ! calculate saved variables
          if(slat1 .ge. 0)then
              s = +1.
          else
              s = -1.
          endif

          colat1 = 90. - s * slat1
          colat2 = 90. - s * slat2

          if(slat1 .eq. slat2)then ! tangent lambert
              n = cosd(90.-s*slat1)
              rconst =       tand(colat1)    / tand(colat1/2.)**n

          else                     ! two standard latitudes
              n = alog(cosd(slat1)/cosd(slat2))/
     1            alog(tand(45.-s*slat1/2.)/tand(45.-s*slat2/2.))
              rconst =      (sind(colat1)/n) / tand(colat1/2.)**n

          endif

          init = 1

      endif

      n_out = n
      s_out = s
      rconst_out = rconst

      return
      end

      subroutine get_grid_spacing_actual(rlat,rlon
     1                                  ,grid_spacing_actual_m,istatus)

cdoc  calculate actual grid spacing at any given lat/lon location

      character*6 c6_maproj

      call get_standard_latitudes(slat1,slat2,istatus)
      if(istatus .ne. 1)then
          return
      endif

      call get_grid_spacing(grid_spacing_m,istatus)
      if(istatus .ne. 1)then
          write(6,*)
     1    ' error calling get_grid_spacing from get_grid_spacing_actual'       
          return
      endif

      call get_c6_maproj(c6_maproj,istatus)
      if(istatus .ne. 1)then
          return
      endif

      if(c6_maproj .eq. 'plrstr')then
          call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                 ,grid_spacing_proj_m)
      else
          grid_spacing_proj_m = grid_spacing_m
      endif

      call get_sigma(rlat,rlon,sigma,istatus)
      if(istatus .ne. 1)then
          write(6,*)
     1    ' error calling get_sigma from get_grid_spacing_actual'       
          return
      endif

      grid_spacing_actual_m = grid_spacing_proj_m / sigma

      return
      end

      subroutine get_grid_spacing_actual_xy(rlat,rlon
     1                        ,grid_spacing_actual_mx
     1                        ,grid_spacing_actual_my
     1                        ,istatus)

cdoc  calculate actual grid spacing (x,y directions) at any given lat/lon 
cdoc  location. this works for conformal or 'latlon' grids

      include 'trigd.inc'

      character*6 c6_maproj

      call get_standard_latitudes(slat1,slat2,istatus)
      if(istatus .ne. 1)then
          return
      endif

      call get_grid_spacing(grid_spacing_m,istatus)
      if(istatus .ne. 1)then
          write(6,*)
     1 ' error calling get_grid_spacing from get_grid_spacing_actual_xy'       
          return
      endif

      call get_c6_maproj(c6_maproj,istatus)
      if(istatus .ne. 1)then
          return
      endif

      if(c6_maproj .ne. 'latlon')then
          if(c6_maproj .eq. 'plrstr')then
              call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                     ,grid_spacing_proj_m)
          else
              grid_spacing_proj_m = grid_spacing_m
          endif

          call get_sigma(rlat,rlon,sigma,istatus)
          if(istatus .ne. 1)then
              write(6,*)
     1        ' error calling get_sigma from get_grid_spacing_actual'       
              return
          endif

          grid_spacing_actual_mx = grid_spacing_proj_m / sigma
          grid_spacing_actual_my = grid_spacing_proj_m / sigma

      else
          grid_spacing_actual_mx = grid_spacing_m * cosd(rlat)
          grid_spacing_actual_my = grid_spacing_m 

      endif

      return
      end


        subroutine get_grid_spacing_cen(grid_spacing_cen_m,istatus)

cdoc    calculate actual grid spacing at the center of the domain
cdoc    if we have a lat/lon domain the y direction spacing will be used

        call get_grid_center(grid_cen_lat,grid_cen_lon,istatus)
        if(istatus .ne. 1)return

        call get_grid_spacing_actual_xy(grid_cen_lat,grid_cen_lon
     1                                 ,grid_spacing_actual_mx
     1                                 ,grid_spacing_actual_my
     1                                 ,istatus)

        grid_spacing_cen_m = grid_spacing_actual_my
 
        return
        end


      subroutine get_ps_parms(slat1,slat2,grid_spacing_m                ! i
     1                       ,phi0,grid_spacing_proj_m)                 ! o

!     1998 steve albers

      include 'trigd.inc'

      logical l_secant
      data l_secant /.true./

!     secant projections are described in "principles of meteorological 
!     analysis", saucier, p. 33. 'phi_std' is the value of phi at the 
!     "standard latitude" as specified by input parameter. 'phi0' is the value 
!     of phi on the actual projection plane utilized for internal map 
!     projection calculations. 'phi0' represents the standard latitude in the
!     more generic sense.

!     we will eventually use the secant projection assumption (unless we run 
!     into software problems)

!     projection is tangent to earth's surface only if phi0 = 90.

      if(slat2 .eq. +90.)then     ! projection pole is at geographic north pole
          phi_std = slat1       

      elseif(slat2 .eq. -90.)then ! projection pole is at geographic south pole
          phi_std = -slat1

      else
          phi_std = +90.        ! we ignore standard lat for local sterographic
                                ! no need for this to be not equal to +90.
      endif

!     calculate grid_spacing_proj_m: grid spacing in the projection plane
      if(l_secant)then
          grid_spacing_proj_m = grid_spacing_m ! equal to parameter value 
                                               ! in namelist files
          phi0 = phi_std

      else ! tangent projection
          grid_spacing_proj_m = grid_spacing_m * 
     1                          2. / (1. + sind(phi_std))
          phi0 = 90.

      endif

      return
      end
