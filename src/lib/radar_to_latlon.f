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


        subroutine radar_to_latlon_old(lat_grid,lon_grid,height_grid
     1                  ,azimuth,slant_range,elev
     1                  ,rlat_radar,rlon_radar,rheight_radar)

        include 'trigd.inc'
        implicit real (a-z)


        if(rlat_radar .eq. 0.0)then
            write(6,*)' warning, radar coords not initialized'
        endif

        rpd = 3.141592653589/180.
        mpd = 111194.
        radius_earth = 6371.e3
        radius_earth_8_thirds = 6371.e3 * 2.6666666

        hor_dist = slant_range * cosd(elev)

        curvature = hor_dist **2 / radius_earth_8_thirds
        height_grid =
     1  slant_range * sind(elev) + curvature + rheight_radar

        height_factor =
     1   (radius_earth + 0.5 * (rheight_radar + height_grid))
     1  /radius_earth

        delta_x = sind(azimuth) * hor_dist / height_factor
        delta_y = cosd(azimuth) * hor_dist / height_factor

        delta_lat = delta_y / mpd

        lat_grid = rlat_radar + delta_lat

        cos_factor =  cosd( 0.5 * (rlat_radar + lat_grid ) )

        delta_lon = delta_x / mpd / cos_factor

        lon_grid = rlon_radar + delta_lon

        return

        end


        subroutine radar_to_latlon(lat_grid,lon_grid,height_grid
     1                  ,azimuth,slant_range,elev
     1                  ,rlat_radar,rlon_radar,rheight_radar)

cdoc    calculate radar echo location given radar location and az/ran/elev 
cdoc    of radar echo

        include 'trigd.inc'
        implicit real (a-z)

        real rlat_radar             ! i    (degrees)
        real rlon_radar             ! i    (degrees)
        real rheight_radar          ! i    (meters)
        real azimuth                ! i    (degrees)
        real slant_range            ! i    (meters)
        real elev                   ! i    (degrees)
        real lat_grid               ! o    (degrees)
        real lon_grid               ! o    (degrees)
        real height_grid            ! o    (meters)

        integer i_status

        if(rlat_radar .eq. 0.0)then
            write(6,*)' warning, radar coords not initialized'
        endif

        rpd = 3.141592653589/180.
        mpd = 111194.
        radius_earth = 6371.e3
        radius_earth_8_thirds = 6371.e3 * 2.6666666

        hor_dist = slant_range * cosd(elev)

        curvature = hor_dist **2 / radius_earth_8_thirds
        height_grid =
     1  slant_range * sind(elev) + curvature + rheight_radar

        height_factor =
     1   (radius_earth + 0.5 * (rheight_radar + height_grid))
     1  /radius_earth

        r_range = hor_dist / height_factor * .001

!       write(12,*)r_range,azimuth

        call razm_lat_lon_gm(
     1          rlat_radar,                                           ! i
     1          rlon_radar,                                           ! i
     1          r_range,                                              ! i
     1          azimuth,                                              ! i
     1          lat_grid,                                             ! o
     1          lon_grid,                                             ! o
     1          i_status )                                            ! o

        if(i_status .ne. 1)then
            write(6,*)
     1         ' error: status check failed after razm_lat_lon_gm call'

            difflat = lat_grid - rlat_radar
            difflon = lon_grid - rlon_radar
            if(r_range .gt. 10000. .and. abs(difflat) .lt. .01
     1                             .and. abs(difflon) .lt. .01)then
                write(6,*)
     1           ' error: qc check failed after razm_lat_lon_gm call'
                write(6,*)'difflat,difflon,r_range'
     1                    ,difflat,difflon,r_range
                write(6,*)'azimuth,rlat_radar,rlon_radar'
     1                    ,azimuth,rlat_radar,rlon_radar
            endif

            stop
        endif

        return

        end

