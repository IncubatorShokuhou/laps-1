cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis
        subroutine plot_obs(k_level,l_ask_questions,asc9_tim
     1    ,i_radar_start,i_radar_end,namelist_parms,plot_parms
     1    ,imax,jmax,kmax,n_plotted,grid_ra_ref,grid_ra_vel,lat,lon
     1    ,topo,grid_spacing_m,mode)

!       steve a         nov  1989       original version
!       steve a         nov  1991       adjustable dimensions

        include 'trigd.inc'
        include 'lapsplot.inc'

        common /zoom/ zoom

        real grid_ra_ref(imax,jmax,kmax)
        real grid_ra_vel(imax,jmax,kmax)
        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)

        real aspect_a(imax,jmax)          ! local

        integer n_plotted(imax,jmax)

        real riob_a(1000000)
        real rjob_a(1000000)
        logical l_plot, l_plot_thisob

        character*150 directory
        character*31 ext,c3_obsext
        character*10  units_2d
        character*125 comment_2d
        character*13 filename13
        character*3 var_2d
        character*1 c1_plottype
        character*8 c8_obstype

        character*13 fileradar
        character*24 atime_24
        character*9 asc9_tim,asc9_tim_radar
        logical l_ask_questions

        data init/0/

        real
     1  azimuth_deg,
     1  x,
     1  y,
     1  spd_kt,
     1  speed_ms,
     1  dir,
     1  du,
     1  projrot,
     1  range_km,
     1  pix_per_km,
     1  mspkt

        data
     1  projrot/0./,
     1  mspkt/.518/

        character
     1  c_obs_type*1,c_map*1,c_mode*1,c_anl*1,c_radial*1,c_pin*1,
     1  c_sfcpro*1,string,wx*25

        common /plotobs/ c_obs_type,c_map,c_mode,c_anl,c_radial

        integer               vert_rad_pirep,
     1                          vert_rad_sao,
     1                          vert_rad_meso,
     1                          vert_rad_prof

!       *** set up constants

        call get_r_missing_data(r_missing_data, istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in plot_obs'
        endif

        call get_vert_rads (    vert_rad_pirep,
     1                          vert_rad_sao,
     1                          vert_rad_meso,
     1                          vert_rad_prof,
     1                          istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in plot_obs'
        endif

!       mode (1: call from 'wd'/'ob' options, 2: call from 'rd'/'vi')

        grid_scale = 3.1
        pix_per_km = 1.65

        x0 = 170.
        y0 = 245.
        pix_per_km = 1.37
        dubase = 2.4

        i_low_offset  = -nint(pix_per_km * 1.51 * grid_scale)
        i_high_offset = i_low_offset 
     1                + nint(pix_per_km * 3.02 * grid_scale)
        dusmall = dubase * pix_per_km/1.65

        jmax_ref = 209
        size_factor = (float(jmax_ref-1) / 300.) 

        if(imax .gt. jmax)then ! horizontally oriented box
            size_factor = size_factor / (float(imax)/float(jmax))       
        endif

        size_factor = size_factor * plot_parms%obs_size

        write(6,*)' subroutine plot_obs: size_factor = ',size_factor

        size_prof = 3.   * size_factor / zoom
        size_pirep = 3.0 * size_factor / zoom
        size_maps = 2.   * size_factor / zoom
        size_vad = 2.    * size_factor / zoom
        size_anl = 0.7   * size_factor / zoom
        size_suw = 1.    * size_factor / zoom
        size_radar = 1.  * size_factor
        size_meso = 2.   * size_factor / zoom

        aspect = 1.0 ! initialize to default value

        if(namelist_parms%l_sphere)then
            do i = 1,imax
            do j = 1,jmax
                arg = cosd(lat(i,j))
                ratio_log = nint(log(arg) / log(0.5))
                projfrac = 0.5**ratio_log

                interval_i = nint(float(interval) / projfrac)

                if(arg .gt. 0.)then
                    aspect_a(i,j) = 1.0 / arg
                else
                    aspect_a(i,j) = 1.0
                endif

            enddo ! j
            enddo ! i

        else
            aspect_a = 1.0

        endif

        lun_in = 5

        call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)

        write(6,*)' subroutine plot_obs at ',asc9_tim,': level '
     1           ,k_level

        call setusv_dum(2hin,201)

!       plot radar obs   *********************************************************************

        write(6,95)
95      format(' want radial velocities     [y,n,a,<ret>=y] ? ',$)
        if(l_ask_questions)read(lun_in,211)c_radial

        if(c_radial(1:1) .eq. 'n')goto205

        write(6,*)' skipping call to ht_of_level to get retheight'       
!       retheight = ht_of_level(k_level)

        c1_plottype = 'y'

        do i_radar = i_radar_start,i_radar_end

          write(6,*)' mode / radar = ',mode,i_radar

          if(mode .eq. 1)then
            call cv_asc_i4time(asc9_tim,i4time_needed)

            if(i_radar .le. 9)then
                write(ext,151)i_radar
151             format('v0',i1)
            else
                write(ext,152)i_radar
152             format('v',i2)
            endif

            var_2d = 'vel'

            write(6,*)' calling get_laps_3dgrid'

            call get_laps_3dgrid(i4time_needed,1200,i4time_found
     1          ,imax,jmax,kmax,ext,var_2d
     1          ,units_2d,comment_2d,grid_ra_vel,istatus)

            if(istatus.ne.1)goto221

            write(6,*)comment_2d

          endif

          do k_grid = 1,kmax

            if(k_grid .gt. k_level .and. k_level .gt. 0)goto221

            do j_grid = 1,jmax
            do i_grid = 1,imax

              if(k_level .eq. 0)then
                k_sfc = 
     1              nint(height_to_zcoord(topo(i_grid,j_grid),istatus))     
              else
                k_sfc = -99
              endif

              if(k_grid  .eq. k_level                   .or.
     1           k_level .eq. 0 .and. k_sfc .eq. k_grid      )then

                if(grid_ra_vel(i_grid,j_grid,k_grid)
     1                          .ne. r_missing_data)then

                  n_plotted(i_grid,j_grid) = n_plotted(i_grid,j_grid)+1      

                  call plot_vr(i_grid,j_grid
     1              ,grid_ra_vel(i_grid,j_grid,k_grid),imax,jmax
     1              ,c1_plottype,n_plotted(i_grid,j_grid))

                endif

              endif ! at the right level to plot

            enddo ! i_grid
            enddo ! j_grid

          enddo ! k_grid

221       continue

        enddo ! i_radar

        if(mode .eq. 2)return

!******  derived radar obs ********************************

205     call cv_asc_i4time(asc9_tim,i4time)

        write(6,210)
210     format(' derived radar obs  [r, ret = none]',30x,'? ',$)    
        if(l_ask_questions)read(lun_in,211)c_obs_type
211     format(a1)

!       if(c_obs_type(1:1) .ne. ' ')then
        if(c_obs_type(1:1) .eq. 'r')then

        write(6,*)
        write(6,*)' derived radar obs'

        call setusv_dum(2hin,5) ! orange

        lun = 61

        do i_radar = i_radar_start,i_radar_end

          if(i_radar .le. 9)then
            write(ext,251)i_radar
 251        format('d0',i1)
          else
            write(ext,252)i_radar
 252        format('d',i2)
          endif

          call get_directory(ext,directory,len_dir)
          open(lun
     1        ,file=directory(1:len_dir)//filename13(i4time,ext(1:3))      
     1        ,status='old',err=1300)

1211      read(61,*,end=1300)ri,rj,k,dir,speed_ms
!421      format(1x,f6.3,f8.3,i2,2f6.1)
          ri = ri + 1.
          rj = rj + 1.
          k_ob = k + 1.

          k_sfc = 
     1    nint(height_to_zcoord(topo(nint(ri),nint(rj)),istatus))

          if(k_ob .eq. k_level .or.
     1       k_level .eq. 0 .and. k_sfc .eq. k_ob)then

!               write(6,421)alat,alon,k,dir,speed_ms
1322            format(' old ',2f8.3,2i5)
1323            format(' new ',2f8.3,10x,2i5)

                spd_kt = speed_ms / mspkt

                if(nint(ri) .ge. 1 .and. nint(ri) .le. imax .and.
     1             nint(rj) .ge. 1 .and. nint(rj) .le. jmax      )then

                    aspect = aspect_a(nint(ri),nint(rj))

                else
                    aspect = 1.0

                endif

                write(6,*)nint(ri),nint(rj),dir,spd_kt,aspect

!               note the 'false' passed as dxx winds are grid north
                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_radar,aspect,'false')

          elseif(k_ob .gt. k_level)then
                goto1300

          endif ! k_ob .eq. k_level

          goto1211

1300      continue

        enddo ! i_radar

        endif


!       plot suw obs   ****************************************************************************************************************

        if(c_obs_type(1:1) .eq. 's')then

222         call setusv_dum(2hin,201)

            lun = 32
            ext = 'suw'
            call get_directory(ext,directory,len_dir)
            open(lun
     1          ,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1          ,status='old',err=240)

230         read(32,*,end=240)azimuth_deg,range_km,elev_deg,dir,speed_ms

            range_m = range_km * 1000.

            range_km = range_km * cosd(elev_deg)

            call radar_to_latlon(alat,alon,retheight
     1                      ,azimuth_deg,range_m,elev_deg
     1                      ,rlat_radar,rlon_radar,rheight_radar)

            k = nint(height_to_zcoord(retheight,istatus))

            if ( k .eq. k_level) then ! plot even range_km rings

              if ( range_km .gt. 20. )  then

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = speed_ms  / mspkt

c               write(6,112)elev_deg,k,range_km,azimuth_deg,dir,spd_kt
112             format(1x,f6.1,i4,2f7.0,4x,2f7.0,i4)

                if(nint(ri) .ge. 1 .and. nint(ri) .le. imax .and.
     1             nint(rj) .ge. 1 .and. nint(rj) .le. jmax      )then

                    aspect = aspect_a(nint(ri),nint(rj))

                else
                    aspect = 1.0

                endif

                call plot_windob(dir,spd_kt,ri,rj,lat,lon
     1                          ,imax,jmax,size_suw,aspect,'true')
              endif

            end if ! valid wind

        goto230
240     continue

        endif

        dularge = dusmall * 2.

!       plot vad wind  *************************************************
        write(6,*)
        write(6,*)' plotting vad wind'

        call setusv_dum(2hin,229)

        lun = 21
        ext = 'vad'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=30)

        do while (.true.)
            read(21,*,end=30)el,dir,speed_ms

            retheight = el * 1000. + rheight_radar

            k = nint(height_to_zcoord(retheight,istatus))

            write(6,*)el,retheight,k

            if(speed_ms .lt. 90. .and. k .eq. k_level)then

                alat = rlat_radar
                alon = rlon_radar

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = speed_ms  / mspkt

                write(6,111)alat,alon,max(dir,-99.),spd_kt
111             format(1x,2f8.1,4x,f7.0,f7.0,i4,2f8.3)

                if(nint(ri) .ge. 1 .and. nint(ri) .le. imax .and.
     1             nint(rj) .ge. 1 .and. nint(rj) .le. jmax      )then

                    aspect = aspect_a(nint(ri),nint(rj))

                else
                    aspect = 1.0

                endif

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_vad,aspect,'true')


            endif
        enddo

 30     continue

!       plot sfc/metar winds  ***********************************************
        write(6,97)
97      format(' plot sfc / profiler winds?     [y,n,<ret>=y] ? ',$)
        if(l_ask_questions)read(lun_in,211)c_sfcpro

        if(c_sfcpro(1:1) .ne. 'n')then

50          write(6,*)' sfc/metar data'

            call setusv_dum(2hin,14) ! royalblue

            lun = 32
            ext = 'sag'
            call get_directory(ext,directory,len_dir)
            open(lun,file=directory(1:len_dir)//
     1           filename13(i4time,ext(1:3)),status='old',err=811)

55          read(32,*,end=60)ri,rj,rk,dir,speed_ms
            k = nint(rk)

            if(k_level .gt. 0)then
                wt_vert = weight_vertical(retheight,k_level,istatus)
            else
                wt_vert = 0.
            endif

            if(abs(k - k_level) .le. vert_rad_sao .or. k_level .eq. 0
     1                                                           )then

                if(nint(ri) .ge. 1 .and. nint(ri) .le. imax .and.
     1             nint(rj) .ge. 1 .and. nint(rj) .le. jmax      )then

                    aspect = aspect_a(nint(ri),nint(rj))

                else
                    aspect = 1.0

                endif

                write(6,111)ri,rj,max(dir,-99.),speed_ms
     1                     ,k,wt_vert,aspect      

                if(k .eq. k_level)then
                    call setusv_dum(2hin,12) ! aqua
                else
                    call setusv_dum(2hin,15) ! slate blue
                endif

                spd_kt = speed_ms  / mspkt

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_meso,aspect,'true')

            endif ! k .eq. k_level

            goto55

60          continue

            close(32)

            call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)

!           plot profile winds  ***********************************************
811         write(6,*)
            write(6,*)
     1          ' profile winds (e.g. profiler, tower, raob, dropsonde)'       

            lun = 32
            ext = 'prg'
            call get_directory(ext,directory,len_dir)
            open(lun,file=directory(1:len_dir)//
     1           filename13(i4time,ext(1:3)),status='old',err=911)

            do while (.true.)
                read(32,*,end=41,err=39)ri,rj,rk,dir,speed_ms,c8_obstype       
 33             format(1x,5f10.0,1x,a8)

                k = nint(rk)
!               write(6,*)k,alat,alon,retheight,dir,speed_ms
!               k_sfc = nint(height_to_zcoord(topo(i,j),istatus))
                k_sfc = 2

                if(k .eq. k_level
     1    .or. k_level .eq. 0 .and. k .eq. k_sfc
     1   .and. dir .ne. r_missing_data
     1   .and. speed_ms .ne. r_missing_data
     1                                                  )then

!                   call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!                   call latlon_ram_laps(alat,alon,x,y,init,'p')

                    spd_kt = speed_ms  / mspkt

                    if(nint(ri) .ge. 1 .and. nint(ri) .le. imax .and.
     1                 nint(rj) .ge. 1 .and. nint(rj) .le. jmax  )then

                        aspect = aspect_a(nint(ri),nint(rj))

                    else
                        aspect = 1.0

                    endif

                    if(c8_obstype(1:4) .eq. 'raob')then
                        call setusv_dum(2hin,3)  ! red
                    else
                        call setusv_dum(2hin,17) ! lavender
                    endif

                    call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                              ,size_prof,aspect,'true')
                    write(6,311,err=38)ri,rj,dir,spd_kt,aspect
     1                                ,c8_obstype
311                 format(1x,2f8.1,4x,f7.0,f7.0,f8.3,1x,a8)
38                  continue

                endif ! k .eq. k_level

39              continue

            enddo

41          continue

            close(32)

        endif ! plot sfc and profiler winds

!       plot pirep winds  ***********************************************
911     write(6,96)
 96     format(' plot acars / pireps winds?     [y,w,n,<ret>=y] ? ',$)
        if(l_ask_questions)read(lun_in,211)c_pin
        write(6,*)
        if(c_pin(1:1) .eq. 'n')then
            write(6,*)' plotting wisdom / cloud drift winds'
            write(6,*)' acars / pireps winds - not plotted'
        elseif(c_pin(1:1) .eq. 'w')then
            write(6,*)' plotting wisdom balloon winds'
            write(6,*)
     1           ' acars / pireps / cloud drift winds - not plotted'      
        else
            write(6,*)' acars / pireps / wisdom / cloud drift winds'
        endif

        lun = 32
        ext = 'pig'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=1011)

        iob = 0 ! initialize for obs density routine
        dist_plot = namelist_parms%dist_plot_ua ! threshold in grid points

        do while (.true.)
            read(32,*,end=51)ri,rj,rk,dir,speed_ms,c3_obsext

            k = nint(rk)

            if(abs(rk - float(k_level)) .le. float(vert_rad_pirep))then

                if(k .eq. k_level)then
                    if(c3_obsext .eq. 'pin')then
                        call setusv_dum(2hin,12) ! aqua
                    elseif(c3_obsext .eq. 'wis')then
                        call setusv_dum(2hin,5)  ! orange
                    else
                        call setusv_dum(2hin,8)  ! green-yellow
                    endif
                else
                    if(c3_obsext .eq. 'pin')then
                        call setusv_dum(2hin,15) ! slate blue
                    elseif(c3_obsext .eq. 'wis')then
                        call setusv_dum(2hin,26) ! pale orange
                    else
                        call setusv_dum(2hin,34) ! grey
                    endif
                endif

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = speed_ms  / mspkt

                if(nint(ri) .ge. 1 .and. nint(ri) .le. imax .and.
     1             nint(rj) .ge. 1 .and. nint(rj) .le. jmax      )then

                    aspect = aspect_a(nint(ri),nint(rj))

                else
                    aspect = 1.0

                endif

!               check if acars/pireps are included
                if(c3_obsext .ne. 'pin' .or. c_pin(1:1) .eq. 'y')then 
                    l_plot_thisob = .true.
                else
                    l_plot_thisob = .false.
                endif

!               option for plotting wisdom only
                if(c_pin(1:1) .eq. 'w')then
                    if(c3_obsext .eq. 'wis')then 
                        l_plot_thisob = .true.
                    else
                        l_plot_thisob = .false.
                    endif
                endif                        

                if(l_plot_thisob)then 
                  call check_ob_density(riob_a,rjob_a,1000000
     1                                 ,iob,ri,rj,dist_plot,l_plot)

                  if(l_plot)then ! passes density criteria
                    iob = iob + 1
                    riob_a(iob) = ri
                    rjob_a(iob) = rj


                    write(6,921)ri,rj,rk,max(dir,-99.),spd_kt,c3_obsext
921                 format(1x,3f8.1,4x,f7.0,f7.0,2x,a3)
                    call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_pirep,aspect,'true')
                  else
                    write(6,922)ri,rj,rk,max(dir,-99.),spd_kt,c3_obsext
922                 format(1x,3f8.1,4x,f7.0,f7.0,2x,a3,' not plotted')
                  endif
 
                endif ! acars/pireps are included

            endif ! k .eq. k_level

        enddo

51      continue

        close(32)



!       plot model winds  ***********************************************
930     write(6,*)
        write(6,*)' model winds'

        call setusv_dum(2hin,203)


        lun = 32
        ext = 'mpg'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=1011)

        do while (.true.)

            read(32,*,end=932)
     1          k_grid,dum,dum,i_grid,j_grid,u_grid,v_grid

            alat = lat(i_grid,j_grid)
            alon = lon(i_grid,j_grid)

            if(k_grid .eq. k_level)then

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')
                call uvgrid_to_disptrue(u_grid,
     1                          v_grid,
     1                          dir,
     1                          speed_ms,
     1                          alat,
     1                          alon)

                spd_kt = speed_ms / mspkt

                if(nint(ri) .ge. 1 .and. nint(ri) .le. imax .and.
     1             nint(rj) .ge. 1 .and. nint(rj) .le. jmax      )then

                    aspect = aspect_a(nint(ri),nint(rj))

                else
                    aspect = 1.0

                endif

!               write(6,111)alat,alon,max(dir,-99.),spd_kt,aspect
                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_maps,aspect,'true')

            endif ! k .eq. k_level

        enddo ! i

932     continue

        close(32)

1011    continue

        return
        end


        subroutine plot_temp_obs(k_level,i4time,imax,jmax,kmax
     1                          ,r_missing_data,lat,lon,topo,zoom
     1                          ,plot_parms)

        include 'lapsplot.inc'

        character*1 c_pin
        character*3 ext
        character*8 c8_obstype
        character*150 directory
        character*150 filename
        character*13 filename13

        logical l_found_file

!       plot temperature obs  ***********************************************

!       size_temp = 8. * float(max(imax,jmax)) / 300.
        size_temp = 1.1 ! 3.33

        write(6,*)
        write(6,*)' plot temperature obs, size_temp = ',size_temp

        write(6,96)
 96     format(' plot acars temps?     [y,n,<ret>=y] ? ',$)
        read(5,211)c_pin
 211    format(a1)

        lun = 32
        ext = 'tmg'
        call get_directory(ext,directory,len_dir)

        filename = directory(1:len_dir)//filename13(i4time,ext(1:3))
        call s_len(filename,len_file)

        open(lun,file=filename(1:len_file),status='old',err=31)
        l_found_file = .true.
        go to 32
       
 31     write(6,*)' could not open ',filename(1:len_file)
        l_found_file = .false.
        go to 42

 32     td = r_missing_data
        p = r_missing_data
        dir = r_missing_data
        spd_kt = r_missing_data
        gust = r_missing_data

        nobs_temp = 0

        do while (.true.) ! count the temperature obs
            read(32,*,end=41,err=36)ri,rj,rk,t_k,c8_obstype
 36         continue

            k = nint(rk)
            k_sfc = 2

            if( k .eq. k_level  
     1                      .or.
     1         (k_level .eq. 0 .and. k .eq. k_sfc)  )then

              if(t_k .ne. r_missing_data)then

                nobs_temp = nobs_temp + 1
 
              endif ! t_k .ne. r_missing_data

            endif ! k .eq. k_level

        enddo

41      continue

        rewind(32)

42      write(6,*)' number of temperature obs = ',nobs_temp

        do while (l_found_file) ! plot the temp obs
            read(32,*,end=141,err=150)ri,rj,rk,t_k,c8_obstype
150         continue

            k = nint(rk)
            k_sfc = 2

            if( k .eq. k_level  
     1                      .or.
     1         (k_level .eq. 0 .and. k .eq. k_sfc)  )then

              if(t_k .ne. r_missing_data)then

                t_c = t_k - 273.15

!               spd_kt = speed_ms  / mspkt

                iflag = 3

                if(c8_obstype(1:3) .eq. 'ras')then      ! rass
                    icol_in = 12 ! aqua
                elseif(c8_obstype(1:3) .eq. 'wis')then  ! wisdom balloons
                    icol_in = 12 ! aqua
                elseif(c8_obstype(1:3) .eq. 'rao')then  ! raob
                    icol_in = 7  ! yellow
                elseif(c8_obstype(1:3) .eq. 'rad')then  ! radiometer
                    icol_in = 5  ! orange
                elseif(c8_obstype(1:3) .eq. 'dro')then  ! dropsonde
                    icol_in = 14 ! royal blue
                elseif(c8_obstype(1:3) .eq. 'tow')then  ! tower
                    icol_in = 14 ! royal blue
                elseif(c8_obstype(1:2) .eq. 'sa')then   ! satsnd
                    icol_in = 17 ! lavender
                elseif(c8_obstype(1:4) .eq. 'goes')then ! goes satellite
                    icol_in = 17 ! lavender
                elseif(c8_obstype(1:4) .eq. 'poes')then ! poes satellite
                    icol_in = 16 ! dark violet
                elseif(c8_obstype(1:5) .eq. 'metar')then ! metar
                    icol_in = 11 ! green
                else                                    ! acars
                    icol_in = 3  ! red
                endif

                call setusv_dum(2hin,icol_in)

                iflag_cv = 0

                if(icol_in .ne. 3 .or. c_pin .eq. 'y')then ! filter out acars?
                  call plot_mesoob(dir,spd_kt,gust,t_c,td,p,ri,rj
     1                            ,lat,lon,imax,jmax,size_temp
     1                            ,zoom,nobs_temp
     1                            ,icol_in,du_loc,wx
     1                            ,iflag,iflag_cv,namelist_parms
     1                            ,plot_parms)


                  write(6,111,err=121)ri,rj,t_c,c8_obstype
111               format(1x,3f8.1,1x,a8)
121               continue

                else
                  write(6,*)' access restricted for ',c8_obstype

                endif

              endif ! t_k .ne. r_missing_data

            endif ! k .eq. k_level

        enddo

141     continue

        close(32)

        return
        end

        subroutine plot_td_obs(k_level,i4time,imax,jmax,kmax
     1                        ,r_missing_data,lat,lon,topo,zoom
     1                        ,namelist_parms,plot_parms
     1                        ,k_mb,mode,field_2d,i_overlay)

        include 'lapsplot.inc'
        include 'icolors.inc'

        character*3 ext, t1
        character*5 pw            
        character*8 c8_obstype
        character*150 directory
        character*150 filename
        character*13 filename13
        character*100 c_label 
        character*10 c_staname

        real field_2d(imax,jmax)

        logical l_found_file

!       plot dewpoint obs  ***********************************************

!       size_temp = 8. * float(max(imax,jmax)) / 300.
        size_temp = 1.1 ! 3.33

        write(6,*)

        if(mode .eq. 1)then
            write(6,*)' plot dewpoint obs, size_td = ',size_temp
            call get_border(imax,jmax,x_1,x_2,y_1,y_2)
            call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax),1)
        elseif(mode .eq. 2)then
            write(6,*)' plot sh obs (from dewpoint), size = ',size_temp
            call get_border(imax,jmax,x_1,x_2,y_1,y_2)
            call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax),1)
        elseif(mode .ge. 3)then
            write(6,*)' plot pw (gps) obs, mode/size = ',mode,size_temp
!           c_label = 'integrated vapor obs (cm)'
            c_label = 'iwv obs, grid-obs (cm)'
            i_overlay = i_overlay + 1
            call setusv_dum(2hin,icolors(i_overlay))
            call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
            call make_fnam_lp(i4time,asc_tim_9,istatus)
            call write_label_lplot(imax,jmax,c_label,asc_tim_9
     1                            ,plot_parms,namelist_parms
     1                            ,i_overlay,'hsect')
            call get_border(imax,jmax,x_1,x_2,y_1,y_2)
            call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax),1)
        endif                   

        lun = 32
        ext = 'hmg'
        call get_directory(ext,directory,len_dir)

        filename = directory(1:len_dir)//filename13(i4time,ext(1:3))
        call s_len(filename,len_file)

        open(lun,file=filename(1:len_file),status='old',err=31)
        l_found_file = .true.
        go to 32
       
 31     write(6,*)' could not open ',filename(1:len_file)
        l_found_file = .false.
        go to 42

 32     continue
        p = r_missing_data
        dir = r_missing_data
        spd_kt = r_missing_data
        gust = r_missing_data

        nobs_temp = 0

        do while (.true.) ! count the dewpoint obs
            read(32,*,end=41,err=36)ri,rj,rk,t_k,c8_obstype,c_staname
 36         continue

            k = nint(rk)
            k_sfc = 2

            if( k .eq. k_level  
     1                      .or.
     1         (k_level .eq. 0 .and. k .eq. k_sfc)  )then

              if(t_k .ne. r_missing_data)then

                nobs_temp = nobs_temp + 1
 
              endif ! t_k .ne. r_missing_data

            endif ! k .eq. k_level

        enddo

41      continue

        rewind(32)

42      write(6,*)' number of dewpoint obs = ',nobs_temp

        do while (l_found_file) ! plot the dewpoint obs
            if(mode .lt. 3)then
              read(32,*,end=141,err=150)ri,rj,rk,t_k,c8_obstype
     1                                              ,c_staname       
            else
              read(32,*,end=141,err=150)ri,rj,rk,t_k,c8_obstype
            endif
150         continue

            k = nint(rk)
            k_sfc = 2

            if( k .eq. k_level  
     1                      .or.
     1         (k_level .eq. 0 .and. k .eq. k_sfc)  )then

              if(t_k .ne. r_missing_data)then

                if(mode .eq. 1)then
                  t_c = t_k - 273.15
                elseif(mode .eq. 2)then
                  t_c = t_k
                elseif(mode .eq. 3)then                  
                  if(t_k .lt. 0.)goto 140
                  t_c = t_k
                  c_staname = ' '

                  if(trim(c8_obstype) .ne. 'gpstpw')then
                    write(6,*)' skipping: ',ri,rj,rk,c8_obstype
                    goto 140         
                  else
                    write(6,*)' keeping:  ',ri,rj,rk,c8_obstype
                  endif
                elseif(mode .eq. 4)then                  
                  if(t_k .lt. 0.)goto 140
                  field_value = field_2d(nint(ri),nint(rj)) * 100.
                  t_c = field_value - t_k 
                  c_staname = ' '

                  if(trim(c8_obstype) .ne. 'gpstpw')then
                    write(6,*)' skipping: ',ri,rj,rk,c8_obstype
                    goto 140         
                  else
                    write(6,*)' keeping:  ',ri,rj,rk,c8_obstype,t_k
     1                                     ,field_value,t_c
                  endif
                endif

!               spd_kt = speed_ms  / mspkt

                iflag = 3

                if(c8_obstype(1:3) .eq. 'ras')then      ! rass
                    icol_in = 12 ! aqua
                elseif(c8_obstype(1:3) .eq. 'rao')then  ! raob
                    icol_in = 7  ! yellow
                elseif(c8_obstype(1:3) .eq. 'rad')then  ! radiometer
                    icol_in = 5  ! orange
                elseif(c8_obstype(1:3) .eq. 'dro')then  ! dropsonde
                    icol_in = 14 ! royal blue
                elseif(c8_obstype(1:2) .eq. 'sa')then   ! satsnd
                    icol_in = 17 ! lavender
                elseif(c8_obstype(1:4) .eq. 'goes')then ! goes satellite
                    icol_in = 17 ! lavender
                elseif(c8_obstype(1:4) .eq. 'poes')then ! poes satellite
                    icol_in = 16 ! dark violet
                elseif(c8_obstype(1:5) .eq. 'metar')then ! metar
                    icol_in = 11 ! green
                else                                    ! acars / gps
                    icol_in = 3  ! red
                endif

                if(mode .lt. 3)then
                    call setusv_dum(2hin,icol_in)
                endif

                if(mode .eq. 2)then ! convert from td to q
                    svp = es(t_c)
                    t_c = (svp / float(k_mb)) * 1000.
                endif

                obs_size = plot_parms%contour_line_width

                zoom_max = 1.5
                if(zoom .lt. zoom_max)then
                    if(nobs_temp .gt. 30)then
                        zoom_eff = 1.0                 ! smaller obs 
                    else
                        zoom_eff = zoom / zoom_max     ! larger obs
                    endif
                else
                    zoom_eff = zoom / zoom_max         ! larger obs
                endif

                zoom_eff = zoom_eff / obs_size

                du = float(jmax) / 252.
                du2 = (du / zoom_eff) * 0.6
                xsta = ri
                ysta = rj
                angd = 0.
                cntr = 0.
                charsize = .0040 / zoom_eff

!               plot station name if it exists
                if(c_staname .ne. 'noname')then
                    call s_len(c_staname,len_name)
                    call pchiqu(xsta, ysta-du2*3.5, 
     1                          c_staname(1:len_name),
     1                          charsize,angd,cntr)      
                endif

!               plot ob location (with a plus sign)
                if(mode .ge. 3 .or.
     1             c8_obstype(1:3) .eq. 'rao')then 
                    call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
                    call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)
                endif

!               plot ob
                if(mode .lt. 3)then
                  write(t1,100,err=101) nint(t_c)
 100              format(i3)
 101              call left_justify(t1)
                  call s_len(t1,len_t1)
                  if(c8_obstype(1:3) .eq. 'rao')then 
                      xloc = xsta + du2*4.0
                  else
                      xloc = xsta
                  endif
                  call pcmequ(xloc,ysta,t1(1:len_t1),charsize,angd,cntr)
                  write(6,102,err=103)ri,rj,t_c,t1,c8_obstype,c_staname
 102              format(1x,3f8.1,1x,a8,1x,a10)
 103              continue
                else ! gps pw
                  write(pw,110,err=111) t_k             
 110              format(f5.2)
 111              call left_justify(pw)
                  call pchiqu(xsta+du2*7.0,ysta+du2*3.0,trim(pw)
     1                       ,charsize,angd,cntr)
 
                  if(mode .eq. 3)then
                    write(6,121,err=122)ri,rj,t_c,pw,c8_obstype
     1                                              ,c_staname
 121                format(1x,3f8.1,1x,a8,1x,a10)
 122                continue

                  elseif(mode .eq. 4)then ! plot difference ob
                    write(pw,110,err=131) t_c             
 131                call left_justify(pw)
                    call pchiqu(xsta+du2*7.0,ysta-du2*3.0,trim(pw)
     1                       ,charsize,angd,cntr)
                    write(6,132,err=133)ri,rj,t_k,t_c,c8_obstype
     1                                               ,c_staname
 132                format(1x,4f8.1,1x,a8,1x,a10)
 133                continue                  

                  endif
                endif

!               write(6,111,err=121)ri,rj,t_c,t1,c8_obstype,c_staname
!11             format(1x,3f8.1,1x,a8,1x,a10)
!21             continue

              endif ! t_k .ne. r_missing_data

            endif ! k .eq. k_level

140     enddo

141     continue

        close(32)

        return
        end

        subroutine plot_td_sndobs(k_level,i4time,imax,jmax,kmax
     1                        ,r_missing_data,lat,lon,topo,zoom
     1                        ,plot_parms)

        include 'lapsplot.inc'

        integer max_snd,max_snd_levels
        parameter (max_snd = 30000)
        parameter (max_snd_levels = 200)

        real heights_3d(imax,jmax,kmax)
        real pres_3d(imax,jmax,kmax)
        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)
         
        real dum1_3d(imax,jmax,kmax) ! dummy variable

!       arrays returned from read_tdsnd
        real lat_tdsnd(max_snd)
        real lon_tdsnd(max_snd)
        real bias_htlow(max_snd)
        real ob_pr_td(max_snd,kmax) ! vertically interpolated tsnd temp
        real inst_err_tdsnd(max_snd)
        character*5 c5_name(max_snd) 
        character*8 c8_sndtype(max_snd) 

        write(6,*)' subroutine plot_td_sndobs'

        i4_window_raob_file = 3600
        ilaps_cycle_time = 3600

        call get_pres_3d(i4time,imax,jmax,kmax,pres_3d,istatus)

!       set up standard atmosphere for heights
        do k = 1,kmax
        do i = 1,imax
        do j = 1,jmax
            heights_3d(i,j,k) = psatoz(pres_3d(i,j,k)*.01) 
        enddo ! j
        enddo ! i
        enddo ! k
       
!       read dewpoint obs from hmg file (generating this if needed)
        call  read_tdsnd(i4time,heights_3d,dum1_3d,               ! input
     1                   pres_3d,                                 ! input
     1                   lat_tdsnd,lon_tdsnd,                     ! output
     1                   lat,lon,                                 ! input
     1                   max_snd,max_snd_levels,                  ! input
     1                   ob_pr_td,inst_err_tdsnd,                 ! output
     1                   c5_name,c8_sndtype,                      ! output
     1                   .true.,.false.,                          ! input
     1                   i4_window_raob_file,                     ! input
     1                   bias_htlow,                              ! output
     1                   n_tdsnd,                                 ! output
     1                   ilaps_cycle_time,                        ! input
     1                   imax,jmax,kmax,                          ! input
     1                   r_missing_data,                          ! input
     1                   istatus)                                 ! output

        write(6,*)' back from read_tdsnd, # soundings = ',n_tdsnd

        return
        end

        subroutine check_ob_density(riob_a,rjob_a,maxobs
     1                             ,iob,ri,rj,dist_plot,l_plot)

!       determine whether to plot an ob based on its distance to the previously
!       plotted obs

        real riob_a(maxobs) ! i/o
        real rjob_a(maxobs) ! i/o
        real dist_plot    ! i    number of gridpoints allowed between obs
        integer iob       ! i/o
        logical l_plot    ! o

        l_plot = .false.

        dist_plot_sq = dist_plot**2

        if(iob .eq. 0)then
            l_plot = .true.
        else
            l_plot = .true.
            do i = 1,iob
                dist_sq = (ri - riob_a(i))**2 + (rj - rjob_a(i))**2
                if(dist_sq .lt. dist_plot_sq)then
                    l_plot = .false.
                endif
            enddo ! i
        endif        

        if(l_plot)then ! add this ob into the arrays of plotted obs
            iob = iob + 1
            riob_a(iob) = ri
            rjob_a(iob) = rj
        endif

        return
        end
