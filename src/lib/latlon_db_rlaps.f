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

        subroutine latlon_db_rlapsgrid(rlat,rlon,lat,lon,ni,nj,ri,rj
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
        implicit real*8 (a-h,o-z)

        real*8 rlat                         ! input lat
        real*8 rlon                         ! input lon
        real*8 lat(ni,nj),lon(ni,nj)        ! input (arrays of lat/lon)
        integer ni,nj                       ! input (laps dimensions)
        real*8 ri,rj                        ! output (i,j on laps grid)
        integer istatus                     ! input / output

        save init,umin,umax,vmin,vmax,niprev,njprev
        data init/0/, niprev/0/, njprev/0/

        include 'grid_fname.cmn'

        if (abs(rlat) > 90d0) then
           if(istatus .ne. 100)then
              write(6,*) 'rejecting invalid latitude ',rlat
           endif
           istatus = -1
           return
        endif

        if(init.ne.nest .or. ni.ne.niprev .or. nj.ne.njprev)then
            call latlon_db_uv(lat(1,1),lon(1,1),umin,vmin,istatus)
            call latlon_db_uv(lat(ni,nj),lon(ni,nj),umax,vmax,istatus)

            write(6,101)umin,umax,vmin,vmax
101         format(1x,' initializing latlon_db_rlapsgrid',4f10.5)
            init = nest
            niprev = ni
            njprev = nj
        endif

        uscale = (umax - umin) / (dble(ni) - 1d0)
        vscale = (vmax - vmin) / (dble(nj) - 1d0)

        u0 = umin - uscale
        v0 = vmin - vscale

!       compute ulaps and vlaps

        call latlon_db_uv(rlat,rlon,ulaps,vlaps,istatus)

        if(umax .le. umin .or. vmax .le. vmin)then
            write(6,*)
     1      ' severe error: uvminmax are off in latlon_db_rlapsgrid'
            write(6,*)umin,umax,vmin,vmax
            stop
        endif

        if(uscale .eq. 0d0 .or. vscale .eq. 0d0)then
            write(6,*)
     1        ' severe error: (u|v)scale = 0 in latlon_db_rlapsgrid'
            write(6,*)lat(1,1),lon(1,1),lat(1,1),lon(ni,nj)
            write(6,*)umin,umax,vmin,vmax
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

        subroutine rlapsgrid_db_latlon(ri,rj,lat,lon,ni,nj,rlat,rlon
     1                                ,istatus)

!       1997            steve albers 

cdoc    this routine assumes a polar stereographic, lambert conformal,
cdoc    or mercator projection.

        implicit real*8 (a-h,o-z)

        real*8 ri,rj                        ! input (i,j on laps grid)
        real*8 lat(ni,nj),lon(ni,nj)        ! input (arrays of lat/lon)
        integer ni,nj                       ! input (laps dimensions)
        real*8 rlat                         ! output lat
        real*8 rlon                         ! output lon
        integer istatus                     ! output

        save init,umin,umax,vmin,vmax
        data init/0/

        include 'grid_fname.cmn'

        if(init .ne. nest)then
            call latlon_db_uv(lat(1,1),lon(1,1),umin,vmin,istatus)
            call latlon_db_uv(lat(ni,nj),lon(ni,nj),umax,vmax,istatus)

            write(6,101)umin,umax,vmin,vmax
101         format(1x,' initializing rlapsgrid_db_latlon',4f10.5)
            init = nest 
        endif

        uscale = (umax - umin) / (dble(ni) - 1d0)
        vscale = (vmax - vmin) / (dble(nj) - 1d0)

        u0 = umin - uscale
        v0 = vmin - vscale

!       compute lat,lon
        ulaps = u0 + uscale * ri
        vlaps = v0 + vscale * rj

        call uv_db_latlon(ulaps,vlaps,rlat,rlon,istatus)

        istatus = 1

        return
        end

        
        subroutine latlon_db_uv(rlat,rlon,u,v,istatus)

!       1997            steve albers 

cdoc    this routine assumes a polar stereographic, lambert conformal,
cdoc    or mercator projection.

        use mem_namelist, only: c6_maproj
     1                  ,slat1_r4=>standard_latitude
     1                  ,slat2_r4=>standard_latitude2
     1                  ,polat_r4=>standard_latitude2
     1                  ,slon_r4=>standard_longitude

        include 'trigd.inc'

        implicit real*8 (a-h,o-z)

        include 'grid_fname.cmn'

        slat1 = slat1_r4
        slat2 = slat2_r4
        polat = polat_r4
        slon = slon_r4

        if(c6_maproj .eq. 'plrstr')then ! polar stereo
            call latlon_db_uv_ps(rlat,rlon,slat1,polat,slon,u,v)

        elseif(c6_maproj .eq. 'lambrt')then ! lambert
            if(abs(slat2).eq.90d0 .and. slat2.ne.slat1)then
               print*,'error: lambert slat2 = 90.'
               istatus=1
               return
            endif
c           slat1 = standard_latitude
c           slat2 = standard_latitude2
c           slon = standard_longitude

            call latlon_db_uv_lc(rlat,rlon,slat1,slat2,slon,u,v)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif
c           slat1  = standard_latitude
c           cenlon = grid_cen_lon_cmn

            call latlon_db_uv_mc(rlat,rlon,slat1,cenlon,u,v)

        elseif(c6_maproj .eq. 'latlon')then ! latlon (cylindrical equidistant)
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif

            call latlon_db_uv_ll(rlat,rlon,cenlon,u,v)

        else
            write(6,*)'error in latlon_db_uv: unrecognized projection '
     1                ,c6_maproj       
            stop

        endif

        istatus = 1

        return
        end


        subroutine latlon_db_uv_ps(dlat_in,dlon_in,dslat,dpolat,dslon
     1                            ,du,dv)

        include 'trigd.inc'

        real*8 dlat_in,dlon_in,dslat,dpolat,dslon,du,dv

!       convert to single precision input vars
        rlat_in = dlat_in
        rlon_in = dlon_in
        slat = dslat
        polat = dpolat
        slon = dslon

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

        endif

        a=90.-rlat
        r = tand(a/2.)      ! consistent with haltiner & williams 1-21

!       b = angle measured counterclockwise from -v axis (zero angle south)
        u =   r * sind(b)
        v = (-r * cosd(b)) * s

!       convert to double precision output vars
        du = u
        dv = v

        return
        end

        subroutine latlon_db_uv_lc(rlat,rlon,slat1,slat2,slon,u,v)

        include 'trigd.inc'

        implicit real*8 (a-h,o-z)

        real*8 n

!       difference between two angles, result is between -180. and +180.
        angdif(x,y)=mod(x-y+540d0,360d0)-180d0

        call lambert_parms_db(slat1,slat2,n,s,rconst)

        r = (tand(45d0-s*rlat/2.))**n
        u =    r*sind(n*angdif(rlon,slon))
        v = -s*r*cosd(n*angdif(rlon,slon))

        return
        end

        subroutine latlon_db_uv_mc(rlat,rlon,slat,cenlon,u,v)

        include 'trigd.inc'

        implicit real*8 (a-h,o-z)

        real*8 pi, rpd

        parameter (pi=3.1415926535897932)
        parameter (rpd=pi/180.)

!       difference between two angles, result is between -180. and +180.
        angdif(x,y)=mod(x-y+540d0,360d0)-180d0
 
        a = 90d0-rlat
        b = cosd(slat)

        u = angdif(rlon,cenlon) * rpd * b
        v = log(1d0/tand(a/2d0))      * b

        return
        end

        subroutine latlon_db_uv_ll(rlat,rlon,cenlon,u,v)

        include 'trigd.inc'

        implicit real*8 (a-h,o-z)

        real*8 pi, rpd

        parameter (pi=3.1415926535897932d0)
        parameter (rpd=pi/180d0)

!       difference between two angles, result is between -180. and +180.
        angdif(x,y)=mod(x-y+540d0,360d0)-180d0
 
        u = angdif(rlon,cenlon) * rpd 
        v = rlat                * rpd

        return
        end


        
        subroutine uv_db_latlon(u,v,rlat,rlon,istatus)

        use mem_namelist, only: c6_maproj
     1                  ,slat1_r4=>standard_latitude
     1                  ,slat2_r4=>standard_latitude2
     1                  ,polat_r4=>standard_latitude2
     1                  ,slon_r4=>standard_longitude

!       1997            steve albers 

!       this routine assumes a polar stereographic, lambert conformal,
!       or mercator projection.

        include 'trigd.inc'

        implicit real*8 (a-h,o-z)

        include 'grid_fname.cmn'

        slat1 = slat1_r4
        slat2 = slat2_r4
        polat = polat_r4
        slon = slon_r4

        if(c6_maproj .eq. 'plrstr')then ! polar stereo
c           slat1 = standard_latitude
c           polat = standard_latitude2
c           slon = standard_longitude

            call uv_db_latlon_ps(u,v,slat1,polat,slon,rlat,rlon)

        elseif(c6_maproj .eq. 'lambrt')then ! lambert
c           slat1 = standard_latitude
c           slat2 = standard_latitude2
c           slon = standard_longitude

            call uv_db_latlon_lc(u,v,slat1,slat2,slon,rlat,rlon)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif
c           slat1  = standard_latitude
c           cenlon = grid_cen_lon_cmn

            call uv_db_latlon_mc(u,v,slat1,cenlon,rlat,rlon)

        elseif(c6_maproj .eq. 'latlon')then ! latlon (cylindrical equidistant)
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'error returned from get_grid_center'
               return
            endif

            call uv_db_latlon_ll(u,v,cenlon,rlat,rlon)

        else
            write(6,*)'error in uv_db_latlon: unrecognized projection '
     1                ,c6_maproj
            stop

        endif

        istatus = 1

        return
        end

        subroutine uv_db_latlon_ps(u,v,slat,polat,slon
     1                                         ,rlat_out,rlon_out)

        include 'trigd.inc'

        implicit real*8 (a-h,o-z)

        real flat,flon,flat_out,flon_out,fpolat,fpolon

        if(abs(polat) .eq. 90.)then
            s = polat / abs(polat)
        else
            s = 1.0
        endif

        r=sqrt(u**2+v**2)

        if (r .eq. 0) then
            rlat=90d0
            rlon=0d0

        else                           
            a=2d0 * atand(r)               ! from haltiner & williams 1-21
            rlat=90d0 - a
            rlon = atan2d(s*v,u)
            rlon = rlon + 90d0
        endif

        if(.true.)then ! rotate considering where the projection pole is
            polon = slon
!           this routine will rotate the longitude even if polat = +90.
            flat = rlat
            flon = rlon
            fpolat = polat
            fpolon = polon
            call pstoge(flat,flon,flat_out,flon_out,fpolat,fpolon)
            rlat_out = flat_out
            rlon_out = flon_out
        else
            rlat_out = rlat
            rlon_out = rlon + slon
        endif


        rlon_out = mod(rlon_out+540d0,360d0) - 180d0 ! convert to -180/+180 range


        return
        end

        subroutine uv_db_latlon_lc(u,v,slat1,slat2,slon,rlat,rlon)

        include 'trigd.inc'

        implicit real*8 (a-h,o-z)

        real*8 n

        call lambert_parms_db(slat1,slat2,n,s,rconst)

!       rlon=slon + atand(-s*u/v) /n
!       rlat=(90.- 2.*atand((-  v/cosd(n*(rlon-slon)))**(1./n)))/s      

        angle  = atan2d(u,-s*v)
        rlat = (90d0- 2d0*atand((-s*v/cosd(angle))**(1d0/n))) / s      
        rlon = slon + angle / n

        rlon = mod(rlon+540d0,360d0) - 180d0       ! convert to -180/+180 range

        return
        end

        subroutine uv_db_latlon_mc(u,v,slat,cenlon,rlat,rlon)

        include 'trigd.inc'

        implicit real*8 (a-h,o-z)

        parameter (pi=3.1415926535897932d0)
        parameter (rpd=pi/180d0)

        b = cosd(slat)

        rlat_abs = 90d0 - atand(exp(-abs(v)/b)) * 2d0

        if(v .gt. 0d0)then
            rlat =  rlat_abs
        else
            rlat = -rlat_abs
        endif

        rlon = u/b/rpd + cenlon
        rlon = mod(rlon+540d0,360d0) - 180d0       ! convert to -180/+180 range

        return
        end

        subroutine uv_db_latlon_ll(u,v,cenlon,rlat,rlon)

        include 'trigd.inc'

        implicit real*8 (a-h,o-z)

        parameter (pi=3.1415926535897932d0)
        parameter (rpd=pi/180d0)

        rlat = v/rpd

        rlon = u/rpd + cenlon
        rlon = mod(rlon+540d0,360d0) - 180d0       ! convert to -180/+180 range

        return
        end

      subroutine lambert_parms_db(slat1,slat2,n_out,s_out,rconst_out)

      include 'trigd.inc'

      implicit real*8 (a-h,o-z)

      real*8 n,n_out

!     we only have to do the calculations once since the inputs are constants
      data init/0/
      save init,n,s,rconst 

      if(init .eq. 0)then ! calculate saved variables
          if(slat1 .ge. 0)then
              s = +1d0
          else
              s = -1d0
          endif

          colat1 = 90d0 - s * slat1
          colat2 = 90d0 - s * slat2

          if(slat1 .eq. slat2)then ! tangent lambert
              n = cosd(90d0-s*slat1)
              rconst =       tand(colat1)    / tand(colat1/2.)**n

          else                     ! two standard latitudes
              n = log(cosd(slat1)/cosd(slat2))/
     1            log(tand(45d0-s*slat1/2.)/tand(45d0-s*slat2/2d0))
              rconst =      (sind(colat1)/n) / tand(colat1/2d0)**n

          endif

          init = 1

      endif

      n_out = n
      s_out = s
      rconst_out = rconst

      return
      end
