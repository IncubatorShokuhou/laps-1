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

!       miscellaneous conversion routines

!       steve albers           1990

        function omega_to_w(omega,pressure_pa)

cdoc    convert omega to w

!       omega units are pascals / second
!       pressure units are pascals
!       w units are meters / second        

        real omega_to_w

        real scale_height
        parameter (scale_height = 8000.) ! meters

        omega_to_w = - (omega / pressure_pa) * scale_height

        return
        end

        function w_to_omega(w,pressure_pa)

cdoc    convert w to omega

!       omega units are pascals / second
!       pressure units are pascals
!       w units are meters / second        

        real w_to_omega

        real scale_height
        parameter (scale_height = 8000.) ! meters

        w_to_omega = - (w * pressure_pa) / scale_height

        return
        end

        subroutine latlon_to_radar(lat_grid,lon_grid,height_grid
     1                            ,azimuth,slant_range,elev
     1                            ,rlat_radar,rlon_radar,rheight_radar)       

cdoc    convert lat/lon/elev to radar azimuth / slant range / elevation angle

        include 'trigd.inc'

        implicit real (a-z)

        if(rlat_radar .eq. 0.0)then
            write(6,*)' warning, radar coords not initialized'
        endif

        rpd = 3.141592653589/180.
        mpd = 111194.
        radius_earth = 6371.e3
        radius_earth_8_thirds = 6371.e3 * 2.6666666
        diff = 0.

        delta_lat = lat_grid - rlat_radar
        delta_lon = lon_grid - rlon_radar


        cos_factor =  cosd( 0.5 * (rlat_radar + lat_grid ) )


        delta_y = delta_lat * mpd

        delta_x = delta_lon * mpd * cos_factor

        height_factor =
     1   (radius_earth + 0.5 * (rheight_radar + height_grid))
     1   /radius_earth

        hor_dist = sqrt(delta_x**2 + delta_y**2) * height_factor

        delta_z = height_grid - rheight_radar

        slant_range = sqrt(hor_dist**2 + delta_z**2)

        if(hor_dist .gt. 0.0)then
            azimuth = atan3(delta_x,delta_y)/rpd
        else
            azimuth = 0.0
        endif

        curvature = hor_dist **2 / radius_earth_8_thirds
        height_diff = delta_z - curvature


        if(slant_range .gt. 0.0)then
            elev = atan2 ( height_diff , hor_dist ) / rpd
        else
            elev = 0.0
        endif


c       write(6,1)slant_range,delta_x,delta_y,azimuth
c       1                       ,hor_dist,curvature
1       format(1x,7f10.0)

        return

        end



        function height_to_zcoord(height_m,istatus) ! using standard atmosphere

!       steve albers fsl
cdoc    convert height to fractional z-coordinate
cdoc    note that this routine works with the standard atmosphere.
cdoc    when the vertical grid is pressure, the height is converted to
cdoc    pressure, then the interpolation to the vertical grid is performed.
cdoc    thus if the height is midway between two laps levels in height space,
cdoc    the value of height_to_zcoord will not have a fraction of 0.5.
cdoc    if the pressure is midway between two laps levels, then the
cdoc    value of height_to_zcoord will have a fraction of 0.5.

        implicit real (a-z)

        integer istatus

        logical ltest_vertical_grid

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        if(ltest_vertical_grid('height'))then
           print*, 'call is obsolete, please report this message to '
           print*, 'and how it occured to laps-bugs@fsl.noaa.gov'
!           height_to_zcoord = height_m / height_interval

        elseif(ltest_vertical_grid('pressure'))then
            pressure_pa = ztopsa(height_m) * 100.
            if(pressure_pa .eq. 9999900.)pressure_pa = 1e-2
            height_to_zcoord = zcoord_of_pressure(pressure_pa)

            if(height_to_zcoord .eq. r_missing_data)then
                istatus = 0
                height_to_zcoord = 0.
            endif

        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' this routine supports pressure or height'
            istatus = 0
            return

        endif

        istatus = 1
        return
        end


        function height_to_zcoord2(height_m,heights_3d
     1                            ,ni,nj,nk,i,j,istatus)

!       1994 steve albers fsl (original)
!       1998 steve albers fsl (overhaul)

cdoc    convert from height to fractional z-coordinate. 3-d heights 
cdoc    are used as a reference.
cdoc    note that this routine works with the real atmosphere.
cdoc    when the vertical grid is pressure (e.g.), the height is converted to
cdoc    pressure, then the interpolation to the vertical grid is performed.
cdoc    thus if the height is midway between two laps levels in height space,
cdoc    the value of 'height_to_zcoord2' will not have a fraction of 0.5.
cdoc    if the pressure is midway between two laps levels, then the
cdoc    value of 'height_to_zcoord2' will have a fraction of 0.5.

        implicit real (a-z)

        integer i,j,k,ni,nj,nk,k_ref,istatus

        real heights_3d(ni,nj,nk)

        logical ltest_vertical_grid,l_valid_grid

        data k_ref /1/
        save k_ref

        data init /0/
        save init,l_valid_grid

        if(init .eq. 0)then ! do this just one time for efficiency
           l_valid_grid = .false.
           if(ltest_vertical_grid('height'))then
              print*, 'height grid not supported in height_to_zcoord2'
!             height_to_zcoord2 = height_m / height_interval
              istatus = 0
              return
           elseif(ltest_vertical_grid('pressure'))then
              l_valid_grid = .true.
           elseif(ltest_vertical_grid('sigma_p'))then
              l_valid_grid = .true.
           else
           endif
           init = 1
        endif

        if(l_valid_grid)then ! valid (pressure) grid with 3-d heights supplied
            height_to_zcoord2 = nk+1 ! default value is off the grid

            k = k_ref

            if(height_m .gt. heights_3d(i,j,nk))then
                height_to_zcoord2 = nk+1 
!               write(6,101)k_ref,height_m,heights_3d(i,j,nk)
!101            format('  note: above domain in height_to_zcoord2,'       
!    1                ,' k_ref,h,h(nk)',i3,2e11.4)
                istatus = 0
                return

            elseif(height_m .lt. heights_3d(i,j,1))then
                height_to_zcoord2 = 0
                write(6,102)k_ref,height_m,heights_3d(i,j,1)
102             format('  warning: below domain in height_to_zcoord2,'
     1                ,' k_ref,h,h(1)',i3,2e11.4)
                istatus = 0
                return

            endif ! input height is outside domain

            do iter = 1,nk
                if(heights_3d(i,j,k+1) .ge. height_m .and.
     1             heights_3d(i,j,k)   .le. height_m         )then
                    thickness = heights_3d(i,j,k+1) - heights_3d(i,j,k)
                    fraction = (height_m - heights_3d(i,j,k))/thickness
                    pressure_low  = zcoord_of_level(k)
                    pressure_high = zcoord_of_level(k+1)
                    diff_log_space = log(pressure_high/pressure_low)
                    pressure = pressure_low 
     1                         * exp(diff_log_space*fraction)
                    height_to_zcoord2 = zcoord_of_pressure(pressure)

                    goto999

                elseif(height_m .gt. heights_3d(i,j,k+1))then
                    k = min(k+1,nk-1)

                elseif(height_m .lt. heights_3d(i,j,k))then
                    k = max(k-1,1)

                endif

            enddo ! iter

            height_to_zcoord2 = 0
            write(6,*)' error, iteration limit in height_to_zcoord2'
            istatus = 0
            return

        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' this routine supports pressure or height'
            istatus = 0
            return

        endif

999     k_ref = k       ! successful return
        istatus = 1
        return

        end


        function height_to_zcoord3(height_m,heights_3d,zcoords_1d
     1                            ,ni,nj,nk,i,j,istatus)

!       2000 steve albers

cdoc    convert from height to fractional z-coordinate. 3-d heights and
cdoc    1-d z-coordinate values are used as a reference.
cdoc    note that this routine works with the real atmosphere.
cdoc    the type of interpolation is similar to that in 'height_to_zcoord2'.
cdoc    when the vertical grid is pressure (e.g.), the height is converted to
cdoc    pressure, then the interpolation to the vertical grid is performed.
cdoc    thus if the height is midway between two laps levels in height space,
cdoc    the value of 'height_to_zcoord3' will not have a fraction of 0.5.
cdoc    if the pressure is midway between two laps levels, then the
cdoc    value of 'height_to_zcoord3' will have a fraction of 0.5.

        implicit real (a-z)

        integer i,j,k,ni,nj,nk,k_ref,istatus

        real heights_3d(ni,nj,nk)
        real zcoords_1d(nk)

        logical ltest_vertical_grid

        data k_ref /1/
        save k_ref

        if(ltest_vertical_grid('height'))then
           print*, 'call is obsolete, please report this message to '
           print*, 'and how it occured to laps-bugs@fsl.noaa.gov'
!          height_to_zcoord3 = height_m / height_interval

        elseif(ltest_vertical_grid('pressure'))then
            height_to_zcoord3 = nk+1 ! default value is off the grid

            k = k_ref

            if(height_m .gt. heights_3d(i,j,nk))then
                height_to_zcoord3 = nk+1 
!               write(6,101)k_ref,height_m,heights_3d(i,j,nk)
!101            format('  note: above domain in height_to_zcoord3,'       
!    1                ,' k_ref,h,h(nk)',i3,2e11.4)
                istatus = 0
                return

            elseif(height_m .lt. heights_3d(i,j,1))then
                height_to_zcoord3 = 0
                write(6,102)k_ref,height_m,heights_3d(i,j,1)
102             format('  warning: below domain in height_to_zcoord3,'
     1                ,' k_ref,h,h(1)',i3,2e11.4)
                istatus = 0
                return

            endif ! input height is outside domain

            do iter = 1,nk
                if(heights_3d(i,j,k+1) .ge. height_m .and.
     1             heights_3d(i,j,k)   .le. height_m         )then
                    thickness = heights_3d(i,j,k+1) - heights_3d(i,j,k)
                    fraction = (height_m - heights_3d(i,j,k))/thickness
                    pressure_low  = zcoords_1d(k)
                    pressure_high = zcoords_1d(k+1)
                    diff_log_space = log(pressure_high/pressure_low)
                    pressure = pressure_low 
     1                         * exp(diff_log_space*fraction)
                    height_to_zcoord3 = zcoord_of_pressure(pressure)

                    goto999

                elseif(height_m .gt. heights_3d(i,j,k+1))then
                    k = min(k+1,nk-1)

                elseif(height_m .lt. heights_3d(i,j,k))then
                    k = max(k-1,1)

                endif

            enddo ! iter

            height_to_zcoord3 = 0
            write(6,*)' error, iteration limit in height_to_zcoord3'
            istatus = 0
            return

        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' this routine supports pressure or height'
            istatus = 0
            return

        endif

999     k_ref = k       ! successful return
        istatus = 1
        return

        end


        function height_to_pressure(height_m,heights_3d
     1                             ,pressures_1d,ni,nj,nk,i,j)

cdoc    convert height to pressure, using 3-d heights and 1-d pressures as
cdoc    a reference.

        implicit real (a-z)
        integer i,j,k,ni,nj,nk

        logical ltest_vertical_grid

        real heights_3d(ni,nj,nk)
        real pressures_1d(nk)

        if(ltest_vertical_grid('height'))then
           print*, 'call is obsolete, please report this message to '
           print*, 'and how it occured to laps-bugs@fsl.noaa.gov'

        elseif(ltest_vertical_grid('pressure'))then
            height_to_pressure = -999. ! default value is off the grid
            do k = 1,nk-1
                if(heights_3d(i,j,k+1) .ge. height_m .and.
     1           heights_3d(i,j,k)   .le. height_m         )then
                    thickness = heights_3d(i,j,k+1) - heights_3d(i,j,k)
                    fraction = (height_m - heights_3d(i,j,k))/thickness
                    pressure_low = pressures_1d(k)
                    pressure_high = pressures_1d(k+1)
                    diff_log_space = log(pressure_high/pressure_low)
                    height_to_pressure
     1          = pressure_low * exp(diff_log_space*fraction)
                endif

            enddo ! k

        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' this routine supports pressure or height'
            istatus = 0
            return

        endif

        return
        end


        subroutine pressure_to_height(pres_pa,heights_3d
     1                               ,ni,nj,nk,i,j,height_out,istatus)       

cdoc    convert pressure to height, using a 3-d height field for reference

        real heights_3d(ni,nj,nk)

        logical ltest_vertical_grid

        if(ltest_vertical_grid('height'))then
            istatus = 0
            return

        elseif(ltest_vertical_grid('pressure'))then
            if(pres_pa .lt. pressure_of_level(nk) .or.
     1         pres_pa .gt. pressure_of_level(1)       )then
                istatus = 0
                return
            endif


            rk = zcoord_of_logpressure(pres_pa)
            if(rk .lt. 1. .or. rk .gt. float(nk))then
                istatus = 0
                return

            else
                if(rk .eq. float(nk))then
                    klow = nk - 1
                    khigh = nk
                else
                    klow = int(rk)
                    khigh = klow + 1
                endif

                frac = rk - klow

                height_out = heights_3d(i,j,klow)  * (1. - frac)
     1                     + heights_3d(i,j,khigh) * (frac)

                istatus = 1
                return

            endif ! rk is within domain

        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' this routine supports pressure or height'
            istatus = 0
            return

        endif

        istatus = 0
        return
        end

        subroutine pres_to_ht(pres_pa,pres_3d,heights_3d
     1                       ,ni,nj,nk,i,j,height_out,istatus)       

cdoc    convert pressure to height, using 3-d pres & ht fields for reference
cdoc    this currently does a linear interpolation in the vertical. we can 
cdoc    change this to a log interpolation later if needed.

        real pres_3d(ni,nj,nk)
        real heights_3d(ni,nj,nk)

!       rk = zcoord_of_logpressure(pres_pa)

        rk = rlevel_of_field(pres_pa,pres_3d,ni,nj,nk,i,j,istatus)
        if(istatus .ne. 1)return

        if(rk .lt. 1. .or. rk .gt. float(nk))then
            istatus = 0
            return

        else
            if(rk .eq. float(nk))then
                klow = nk - 1
                khigh = nk
            else
                klow = int(rk)
                khigh = klow + 1
            endif

            frac = rk - klow

            height_out = heights_3d(i,j,klow)  * (1. - frac)
     1                 + heights_3d(i,j,khigh) * (frac)

            istatus = 1
            return

        endif ! rk is within domain

        istatus = 0
        return
        end

        subroutine pres_to_ht_2d(pres_pa_2d,pres_pa_3d,heights_3d
     1                          ,temp_3d,ni,nj,nk,height_out_2d,istatus)       

        real pres_pa_3d(ni,nj,nk)
        real heights_3d(ni,nj,nk)
        real temp_3d(ni,nj,nk)
        real pres_pa_2d(ni,nj)
        real height_out_2d(ni,nj)

        include 'constants.inc'
        real c1
        parameter (c1 = ep_1) 

        real c2
        parameter (c2 = r_d / grav)  

!       interpolate to find heights on the input pressure surface
!       simplified for now until we take the layer mean temp

        do i = 1,ni
        do j = 1,nj

!         find closest pressure level
          delta_p_min = 1e10
          do k = 1,nk
              delta_p_abs = abs(pres_pa_2d(i,j) - pres_pa_3d(i,j,k))
              if(delta_p_abs .lt. delta_p_min)then
                  delta_p_min = delta_p_abs
                  kref = k
              endif
          enddo ! k

          delta_p = pres_pa_2d(i,j) - pres_pa_3d(i,j,kref)
          delta_h_simple = -delta_p_min * 0.1 ! simplified hypsometric equation

          alog_term = alog(pres_pa_2d(i,j) / pres_pa_3d(i,j,kref))
          t_k = temp_3d(i,j,kref)
          delta_h = -c2 * t_k * alog_term

          height_out_2d(i,j) = heights_3d(i,j,kref) + delta_h

        enddo ! j
        enddo ! i

        return
        end


        function height_of_level(level)

cdoc    calculate the height of a given pressure level, using standard atmos.
cdoc    works only for constant pressure levels

        implicit real (a-z)

        integer level,istatus

        logical ltest_vertical_grid

        if(ltest_vertical_grid('height'))then
           print*, 'call is obsolete, please report this message to '
           print*, 'and how it occured to laps-bugs@fsl.noaa.gov'
!            height_of_level = height_interval * level

        elseif(ltest_vertical_grid('pressure'))then
            height_of_level = psatoz(pressure_of_level(level) * .01)

        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' this routine supports pressure or height'
            istatus = 0
            return

        endif

        istatus = 1
        return
        end

        function zcoord_of_level(level)

cdoc    calculate zcoord (e.g. pressure) of a given level. 
cdoc    works only for constant pressure levels
 
        use mem_namelist, only: nk_laps

        real, allocatable, dimension(:) :: sigma_1d_out

        logical ltest_vertical_grid_lc

        if(ltest_vertical_grid_lc('height'))then
           print*, 'call is obsolete, please report this message to '
           print*, 'and how it occured to laps-bugs@fsl.noaa.gov'
!           zcoord_of_level = height_interval * level

        elseif(ltest_vertical_grid_lc('pressure'))then
            zcoord_of_level = pressure_of_level(level)

        elseif(ltest_vertical_grid_lc('sigma_ht'))then
            allocate(sigma_1d_out(nk_laps))
            call get_sigma_1d(nk,sigma_1d_out,istatus)
            zcoord_of_level = sigma_1d_out(level)
            deallocate (sigma_1d_out)

        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' this routine supports pressure or height'
            istatus = 0
            return

        endif

        istatus = 1
        return
        end


        function pressure_of_level(level)

cdoc    calculate pressure of a given integer level. 
cdoc    works only for constant pressure levels

        real, allocatable, dimension(:) :: pres_1d

        integer level, istatus, istat_alloc

        call get_laps_dimensions(nk,istatus)
        if(istatus .ne. 1)stop

        allocate(pres_1d(nk), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate pres_1d'
            stop
        endif

        call get_pres_1d(i4time,nk,pres_1d,istatus)
        if(istatus .ne. 1)stop

        pressure_of_level = pres_1d(level)

        deallocate(pres_1d)

        istatus = 1
        return
        end


        function pressure_of_rlevel(rlevel)

cdoc    obtain pressure of a given real (fractional) level. 
cdoc    works only for constant pressure levels

        call get_laps_dimensions(nk,istatus)
        if(istatus .ne. 1)stop

        l1 = int(rlevel)
        if(l1 .eq. nk)l1 = l1-1
        l2 = int(rlevel)+1
        frac = rlevel - l1

        pres1 = pressure_of_level(l1)
        pres2 = pressure_of_level(l2)

        pressure_of_rlevel = pres1 * (1. - frac) + pres2 * frac

        istatus = 1
        return
        end


        function rlevel_of_field(value,field_3d,ni,nj,nk,i,j,istatus)       

cdoc    find z coordinate given a field value, i, j, and the whole 3-d field

        implicit real (a-z)

        integer i,j,k,ni,nj,nk,k_ref,istatus,isign

        real field_3d(ni,nj,nk)

        logical ltest_vertical_grid

        data k_ref /1/
        save k_ref

        if(.true.)then                            
            if(field_3d(i,j,nk) .gt. field_3d(i,j,1))then
                rsign = 1.0
                isign = 1
            else
                rsign = -1.0
                isign = -1
            endif

            rlevel_of_field = nk+1 ! default value is off the grid

            k = k_ref

            if((value - field_3d(i,j,nk)) * rsign .gt. 0.)then
                rlevel_of_field = nk+1 
!               write(6,101)k_ref,value,field_3d(i,j,nk)
!101            format('  note: above domain in rlevel_of_field,'       
!    1                ,' k_ref,h,h(nk)',i3,2e11.4)
                istatus = 0
                return

            elseif((value - field_3d(i,j,1)) * rsign .lt. 0.)then
                rlevel_of_field = 0
                write(6,102)isign,value,field_3d(i,j,1),field_3d(i,j,nk)
102             format('  warning: below domain in rlevel_of_field,'
     1                ,' isign,h,h(1),h(nk)',i3,3e11.4)
                write(6,*)i,j,field_3d(i,j,:)
                istatus = 0
                return

            endif ! input height is outside domain

            do iter = 1,nk
                if( (field_3d(i,j,k+1) - value) * rsign .ge. 0. .and.
     1              (field_3d(i,j,k)   - value) * rsign .le. 0.  )then
                    thickness = field_3d(i,j,k+1) - field_3d(i,j,k)
                    fraction = (value - field_3d(i,j,k))/thickness

                    rlevel_of_field = k + fraction

                    goto999

                elseif((value - field_3d(i,j,k+1)) * rsign .gt. 0.)then
                    k = min(k+1,nk-1)

                elseif((value - field_3d(i,j,k))   * rsign .lt. 0.)then
                    k = max(k-1,1)

                endif

            enddo ! iter

            rlevel_of_field = 0
            write(6,*)' error, iteration limit in rlevel_of_field'
            istatus = 0
            return

        endif

999     k_ref = k       ! successful return
        istatus = 1

        return
        end


        function rlevel_of_logfield(value,field_3d,ni,nj,nk,i,j,istatus)       

cdoc    find z coordinate given a field value, i, j, and the whole 3-d field
cdoc    log vertical interpolation is used. 

        implicit real (a-z)

        integer i,j,k,ni,nj,nk,k_ref,istatus,isign

        real field_3d(ni,nj,nk)

        logical ltest_vertical_grid

        data k_ref /1/
        save k_ref

        if(ltest_vertical_grid('height'))then
            print*, 'call is obsolete, please report this message to '       
            print*, 'and how it occured to laps-bugs@fsl.noaa.gov'
!           rlevel_of_logfield = value / height_interval

        elseif(ltest_vertical_grid('pressure'))then
            if(field_3d(i,j,nk) .gt. field_3d(i,j,1))then
                rsign = 1.0
                isign = 1
            else
                rsign = -1.0
                isign = -1
            endif

            rlevel_of_logfield = nk+1 ! default value is off the grid

            k = k_ref

            if((value - field_3d(i,j,nk)) * rsign .gt. 0.)then
                rlevel_of_logfield = nk+1 
!               write(6,101)k_ref,value,field_3d(i,j,nk)
!101            format('  note: above domain in rlevel_of_logfield,'       
!    1                ,' k_ref,h,h(nk)',i3,2e11.4)
                istatus = 0
                return

            elseif((value - field_3d(i,j,1)) * rsign .lt. 0.)then
                rlevel_of_logfield = 0
                write(6,102)k_ref,value,field_3d(i,j,1)
102             format('  warning: below domain in rlevel_of_logfield,'
     1                ,' k_ref,h,h(1)',i3,2e11.4)
                istatus = 0
                return

            endif ! input height is outside domain

            do iter = 1,nk
                if( (field_3d(i,j,k+1) - value) * rsign .ge. 0. .and.
     1              (field_3d(i,j,k)   - value) * rsign .le. 0.  )then
                    thickness = field_3d(i,j,k+1) - field_3d(i,j,k)

                    fraction = log(value/field_3d(i,j,k)) / 
     1                         log(field_3d(i,j,k+1)/field_3d(i,j,k))

                    rlevel_of_logfield = k + fraction

                    goto999

                elseif((value - field_3d(i,j,k+1)) * rsign .gt. 0.)then
                    k = min(k+1,nk-1)

                elseif((value - field_3d(i,j,k))   * rsign .lt. 0.)then
                    k = max(k-1,1)

                endif

            enddo ! iter

            rlevel_of_logfield = 0
            write(6,*)' error, iteration limit in rlevel_of_logfield'
            istatus = 0
            return

        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' this routine supports pressure or height'
            istatus = 0
            return

        endif

999     k_ref = k       ! successful return
        istatus = 1

        return
        end


        function zcoord_of_pressure(pres_pa)

cdoc    convert pressure to a real (fractional) level. 
cdoc    works only for constant pressure levels.

        real, allocatable, dimension(:) :: pres_1d

        logical ltest_vertical_grid,l_valid_grid

        data init /0/
        save init,l_valid_grid

        if(init .eq. 0)then ! do this just one time for efficiency
           l_valid_grid = .false.
           if(ltest_vertical_grid('height'))then
              print*, 'height grid not supported in height_to_zcoord2'
              istatus = 0
              return
           elseif(ltest_vertical_grid('pressure'))then
              l_valid_grid = .true.
           endif
           init = 1
        endif

        if(l_valid_grid)then
            call get_laps_dimensions(nk,istatus)
            if(istatus .ne. 1)stop

            allocate(pres_1d(nk), stat=istat_alloc )
            if(istat_alloc .ne. 0)then
                write(6,*)' error: could not allocate pres_1d'
                stop
            endif

            call get_pres_1d(i4time,nk,pres_1d,istatus)
            if(istatus .ne. 1)stop

            arg = rlevel_of_field(pres_pa,pres_1d,1,1,nk,1,1,istatus)       

            if(istatus .ne. 1)then
                call get_r_missing_data(r_missing_data,istatus)
                if(istatus .ne. 1)stop
                zcoord_of_pressure = r_missing_data

            else
                zcoord_of_pressure = arg

            endif    

            deallocate(pres_1d)

        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' zcoord_of_pressure supports pressure or height'  
            istatus = 0
            return

        endif

        istatus = 1
        return
        end


        function zcoord_of_logpressure(pres_pa)

cdoc    convert pressure to a real (fractional) level in log space. 
cdoc    works only for constant pressure levels.

        logical ltest_vertical_grid

        if(ltest_vertical_grid('height'))then

        elseif(ltest_vertical_grid('pressure'))then
            call get_r_missing_data(r_missing_data,istatus)

            rz = zcoord_of_pressure(pres_pa)

            if(rz .eq. r_missing_data)then
                zcoord_of_logpressure = r_missing_data
                return
            endif

            call get_laps_dimensions(nk,istatus)

            if(rz .eq. float(nk))then
                zcoord_of_logpressure = rz
                return
            endif

            pres_low  = pressure_of_level(int(rz))
            pres_high = pressure_of_level(int(rz)+1)

            frac = log(pres_pa/pres_low) / log(pres_high/pres_low)

            zcoord_of_logpressure = int(rz) + frac

        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' this routine supports pressure or height'
            istatus = 0
            return

        endif

        istatus = 1
        return
        end


        subroutine   uvgrid_to_radar(u_grid,     ! i
     1                               v_grid,     ! i
     1                               t_radar,    ! o
     1                               r_radar,    ! o
     1                               azimuth,    ! i
     1                               latitude,   ! i
     1                               longitude)  ! i

cdoc    convert u and v (grid north) to tangential and radial velocity,
cdoc    given the radar azimuth and geographic longitude.

        real longitude

        call   uvgrid_to_disptrue(u_grid,
     1                            v_grid,
     1                            di_true,
     1                            speed,
     1                            latitude,
     1                            longitude)

        call disptrue_to_radar(di_true,
     1                         speed,
     1                         t_radar,
     1                         r_radar,
     1                         azimuth)

        return
        end


        subroutine   uvtrue_to_radar(u_true,  ! i
     1                               v_true,  ! i
     1                               t_radar, ! o
     1                               r_radar, ! o
     1                               azimuth) ! i

cdoc    convert u and v (true north) to tangential and radial velocity,
cdoc    given the radar azimuth.

        call   uv_to_disp(u_true,
     1                    v_true,
     1                    di_true,
     1                    speed)

        call disptrue_to_radar(di_true,
     1                         speed,
     1                         t_radar,
     1                         r_radar,
     1                         azimuth)

        return
        end


        subroutine   uvgrid_to_disptrue(u_grid,    ! i
     1                                  v_grid,    ! i
     1                                  di_true,   ! o
     1                                  speed,     ! o
     1                                  latitude,  ! i
     1                                  longitude) ! i

cdoc    convert u and v (grid north) to dir and speed (true north),
cdoc    given the longitude.

        real latitude, longitude

        call get_config(istatus)

        if(istatus .ne. 1)then
            write(6,*)' error, get_laps_config not successfully called'       
            stop
        endif

        speed = sqrt( u_grid **2 + v_grid **2 )

        if(speed .gt. 0)then
            di_grid = atan3d(-u_grid,-v_grid)
!           switched the sign on this (in 2016)
            di_true = di_grid - projrot_latlon(latitude,longitude
     1                                        ,istatus)
            di_true = mod(di_true+360.,360.)
        else
            di_true = 0.
        endif

        return
        end



        subroutine   disptrue_to_radar(di_true, ! i
     1                                 speed,   ! i
     1                                 t_radar, ! o
     1                                 r_radar, ! o
     1                                 azimuth) ! i

cdoc    convert dir and speed (true north) to tangential and radial velocity,
cdoc    given the radar azimuth.

        include 'trigd.inc'

!       real longitude

        angle_diff = di_true - azimuth
        t_radar = speed * (-sind(angle_diff))
        r_radar = speed * (-cosd(angle_diff))

        return
        end


        subroutine   radar_to_uvgrid(t_radar,   ! i
     1                               r_radar,   ! i
     1                               u_grid,    ! o
     1                               v_grid,    ! o
     1                               azimuth,   ! i
     1                               latitude,  ! i
     1                               longitude) ! i

cdoc    convert tangential and radial velocity to u and v (grid north),
cdoc    given the radar azimuth and geographic longitude.

        real latitude, longitude

        call radar_to_disptrue(di_true,
     1                         speed,
     1                         t_radar,
     1                         r_radar,
     1                         azimuth)

        call   disptrue_to_uvgrid(di_true,
     1                            speed,
     1                            u_grid,
     1                            v_grid,
     1                            longitude)


        return
        end

        subroutine   radar_to_uvtrue(t_radar,  ! i
     1                               r_radar,  ! i
     1                               u_true,   ! o
     1                               v_true,   ! o
     1                               azimuth)  ! i

cdoc    convert tangential and radial velocity to u and v (true north),
cdoc    given the radar azimuth.

        call radar_to_disptrue(di_true,
     1                         speed,
     1                         t_radar,
     1                         r_radar,
     1                         azimuth)

        call   disp_to_uv(di_true,
     1                    speed,
     1                    u_true,
     1                    v_true)


        return
        end

        subroutine   radar_to_disptrue(di_true,  ! o
     1                                 speed,    ! o
     1                                 t_radar,  ! i
     1                                 r_radar,  ! i
     1                                 azimuth)  ! i

cdoc    convert tangential and radial velocity to dir and speed (true north),
cdoc    given the radar azimuth.

!       real longitude

        speed = sqrt(t_radar**2 + r_radar**2)

        if(speed .gt. 0)then
            di_radar = atan3d(-t_radar,-r_radar)
            di_true = mod(di_radar + azimuth + 360.,360.)
        else
            di_true = 0.
        endif

        return
        end

        subroutine   disptrue_to_uvgrid(di_true,     ! i
     1                                  speed,       ! i
     1                                  u_grid,      ! o
     1                                  v_grid,      ! o
     1                                  longitude)   ! i

cdoc    convert dir and speed (true north) to u and v (grid north)

        real latitude, longitude

        call get_config(istatus)

        if(istatus .ne. 1)then
            write(6,*)' error, get_laps_config not successfully called'       
            stop
        endif

        latitude = -999. ! since lat is not yet passed in

!       switched the sign on this (in 2016)
        di_grid = di_true + projrot_latlon(latitude,longitude
     1                                             ,istatus)

        call         disp_to_uv(di_grid,
     1                  speed,
     1                  u_grid,
     1                  v_grid)

        return
        end

        subroutine   uvtrue_to_uvgrid(u_true,    ! i
     1                                v_true,    ! i
     1                                u_grid,    ! o
     1                                v_grid,    ! o
     1                                longitude) ! i

cdoc    convert wind vector from true north to grid north, given the longitude.

        real latitude, longitude

        call get_config(istatus)

        if(istatus .ne. 1)then
            write(6,*)' error, get_laps_config not successfully called'       
            stop
        endif

        latitude = -999. ! since lat is not yet passed in

        angle = projrot_latlon(latitude,longitude,istatus)

        call         rotate_vec(u_true,
     1                  v_true,
     1                  u_grid,
     1                  v_grid,
     1                  angle)

        return
        end

        subroutine   uvtrue_to_uvgrid_2d(u_true,    ! i
     1                                   v_true,    ! i
     1                                   u_grid,    ! o
     1                                   v_grid,    ! o
     1                                   longitude, ! i
     1                                   ni,        ! i
     1                                   nj)        ! i

cdoc    convert wind vector from true north to grid north, given the longitude.

        real u_true(ni,nj), v_true(ni,nj)
        real u_grid(ni,nj), v_grid(ni,nj)
        real latitude(ni,nj), longitude(ni,nj)

        real angle(ni,nj)

        call get_config(istatus)

        if(istatus .ne. 1)then
            write(6,*)' error, get_laps_config not successfully called'       
            stop
        endif

        latitude = -999. ! since lat is not yet passed in
 
        call projrot_latlon_2d(latitude,longitude,ni,nj,angle,istatus)

        call         rotate_vec_2d(u_true,
     1                  v_true,
     1                  u_grid,
     1                  v_grid,
     1                  angle,
     1                  ni,nj)

        return
        end

        subroutine   uvgrid_to_uvtrue(u_grid,    ! i
     1                                v_grid,    ! i
     1                                u_true,    ! o
     1                                v_true,    ! o
     1                                longitude) ! i

cdoc    convert wind vector from grid north to true north, given the longitude

        real latitude, longitude

        call get_config(istatus)

        if(istatus .ne. 1)then
            write(6,*)' error, get_laps_config not successfully called'       
            stop
        endif

        latitude = -999. ! since lat is not yet passed in

        angle = -projrot_latlon(latitude,longitude,istatus)


        call         rotate_vec(u_grid,
     1                          v_grid,
     1                          u_true,
     1                          v_true,
     1                          angle)

        return
        end



      subroutine rotate_vec(u1,v1,u2,v2,angle)

cdoc  rotate vector (u1,v1) through a clockwise angle to obtain vector (u2,v2)
cdoc  u points east and v points north

      include 'trigd.inc'
      u2 =  u1 * cosd(angle) + v1 * sind(angle)
      v2 = -u1 * sind(angle) + v1 * cosd(angle)

      return
      end

      subroutine rotate_vec_2d(u1,v1,u2,v2,angle,ni,nj)

cdoc  rotate vector (u1,v1) through a clockwise angle to obtain vector (u2,v2)
cdoc  u points east and v points north

      parameter (pi = 3.1415926535897932)
      parameter (rpd = pi / 180.)

      real u1(ni,nj),v1(ni,nj),u2(ni,nj),v2(ni,nj),angle(ni,nj)

      u2(:,:) =  u1(:,:) * cos(angle(:,:)*rpd) 
     1         + v1(:,:) * sin(angle(:,:)*rpd)
      v2(:,:) = -u1(:,:) * sin(angle(:,:)*rpd) 
     1         + v1(:,:) * cos(angle(:,:)*rpd)

      return
      end


        subroutine   disp_to_uv(dir,
     1                          speed,
     1                          u,
     1                          v)
cdoc    convert dir and speed to u and v

        include 'trigd.inc'
        u  = - sind(dir) * speed
        v  = - cosd(dir) * speed

        return
        end


        subroutine   uv_to_disp(u,
     1                          v,
     1                          dir,
     1                          speed)

cdoc    convert u and v to dir and speed

        speed = sqrt( u**2 + v**2 )

        if(speed .gt. 0)then
            dir = atan3d(-u,-v)
!           dir = mod(dir+360.,360.)
        else
            dir = 0.
        endif

        return
        end

        function k_to_f(x)
cdoc    convert kelvin to fahrenheit
        real k_to_f
        k_to_f = (x - 273.15) * 1.8 + 32.
        return
        end

        function f_to_k(x)
cdoc    convert fahrenheit to kelvin
        real f_to_k
        f_to_k = (x - 32.) / 1.8 + 273.15
        return
        end

        function k_to_c(x)
cdoc    convert kelvin to celsius
        real k_to_c
        k_to_c = (x - 273.15)
        return
        end

        function c_to_k(x)
cdoc    convert celsius to kelvin
        real c_to_k
        c_to_k = (x + 273.15)
        return
        end

        function f_to_c(x)
cdoc    convert fahrenheit to celsius
        real f_to_c
        f_to_c = (x - 32.) / 1.8 
        return
        end
c
        function c_to_f(t_c)
cdoc    convert celsius to fahrenheit
        c_to_f = (t_c * 9./5.) + 32.             ! c to f
        return
        end
c
c

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine conv_meters_to_inches(data_in, numarr,prodno,
     1  imax, jmax)
c   subroutine to convert accumulations meters to inches
        real data_in(imax, jmax, numarr)
        integer prodno
c   begin
        do j=1,jmax
        do i=1,imax
          data_in(i,j,prodno) = data_in(i,j,prodno) / 0.0254
        enddo !i
        enddo !j
c
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

