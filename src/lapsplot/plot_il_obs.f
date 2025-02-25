
        subroutine plot_il_obs(k_level,i4time,imax,jmax,kmax
     1                        ,r_missing_data,lat,lon,topo,zoom
     1                        ,i_overlay,namelist_parms,plot_parms)

        include 'lapsplot.inc'
        include 'icolors.inc'

        integer max_snd,max_snd_levels
        parameter (max_snd = 10000)
        parameter (max_snd_levels = 200)

        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)
         
!       profile stuff
        real lat_pr(max_snd)
        real lon_pr(max_snd)
        real elev_pr(max_snd)

        integer i4time_ob_pr(max_snd)
        integer nlevels_obs_pr(max_snd)

        character*8 obstype(max_snd)
        character*3 ext, t1
        character*100 c_label
        character*9 asc_tim_9

!       arrays returned from read_snd_metadata
        real lat_tdsnd(max_snd)
        real lon_tdsnd(max_snd)
        real cloud_base_temp(max_snd)
        real cloud_liquid(max_snd)
        character*5 c5_name(max_snd)
        character*8 c8_sndtype(max_snd) 

        write(6,*)' subroutine plot_il_obs'

        i4_window_raob_file = 3600
        ilaps_cycle_time = 3600
        ext = 'snd'

        write(6,*)' calling read_snd_metadata'

        call read_snd_metadata(lun,i4time,ext                         ! i
     1                        ,max_snd,max_snd_levels                 ! i
     1                        ,lat,lon,imax,jmax                      ! i
     1                        ,n_profiles                             ! o
     1                        ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr   ! o
     1                        ,c5_name,i4time_ob_pr,obstype           ! o
     1                        ,cloud_base_temp,cloud_liquid           ! o
     1                        ,istatus)                               ! o

        write(6,*)' back from read_snd_metadata, # soundings = '
     1            ,n_profiles

!       plot label
        c_label = 'integrated liquid obs (mm)'
        i_overlay = i_overlay + 1
        call setusv_dum(2hin,icolors(i_overlay))
        call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
        call make_fnam_lp(i4time,asc_tim_9,istatus)
        call write_label_lplot(imax,jmax,c_label,asc_tim_9
     1                            ,plot_parms,namelist_parms
     1                            ,i_overlay,'hsect')       

!       plot obs
        k_level = 0
        k_mb = 0
        mode = 3

        call get_border(imax,jmax,x_1,x_2,y_1,y_2)
        call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax),1)

        write(6,*)' plotting integrated liquid obs from snd file'
        do i = 1,n_profiles
            if(cloud_liquid(i) .ne. r_missing_data)then
                write(6,*)' i,c5_name,cloud_liquid= ',i,c5_name(i)
     1                                               ,cloud_liquid(i)

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

                call latlon_to_rlapsgrid(lat_pr(i),lon_pr(i),lat,lon
     1                                  ,imax,jmax,xsta,ysta,istatus)

                du = float(jmax) / 252.
                du2 = (du / zoom_eff) * 0.6
                angd = 0.
                cntr = 0.
                charsize = .0040 / zoom_eff

!               plot station name & location
                call s_len(c5_name(i),len_name)
                call pchiqu(xsta, ysta-du2*3.5, 
     1                      c5_name(i)(1:len_name),
     1                      charsize,angd,cntr)      
                call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
                call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)

!               plot observation
                rplot = cloud_liquid(i)*1000.
                write(t1,100,err=20) rplot
 20             call left_justify(t1)
                call s_len(t1,len_t1)
                call pcmequ(xsta+du2*4.0,ysta,t1(1:len_t1),charsize
     1                                                    ,angd,cntr)       
 100            format(f3.1)

                write(6,111,err=121)xsta,ysta,rplot
     1                             ,t1,obstype(i),c5_name(i)
111             format(1x,2f8.1,f8.2,1x,a8,1x,a10)
121             continue

            endif
        enddo ! i

        return
        end

