
        subroutine plot_sounding(i4time_ref,lun,NX_L,NY_L,NZ_L
     1                          ,r_missing_data,laps_cycle_time,maxstns
     1                          ,i_overlay,plot_parms,namelist_parms)       

        include 'lapsplot.inc'

        real pres_3d(NX_L,NY_L,NZ_L)
        real field_3d(NX_L,NY_L,NZ_L)
!       real rh_3d(NX_L,NY_L,NZ_L)

        real pres_2d(NX_L,NY_L)
        real t_2d(NX_L,NY_L)
        real td_2d(NX_L,NY_L)
        real pw_2d(NX_L,NY_L)
        real lat(NX_L,NY_L)
        real lon(NX_L,NY_L)

        real temp_vert(NZ_L)
        real rh_vert(NZ_L)
        real sh_vert(NZ_L)
        real td_vert(NZ_L)

        real pres_1d(NZ_L)
        real logp_1d(NZ_L), logp_bottom, logp_top, logp, logp_sfc
  
        real k_to_c, make_td

        character*1 c_prodtype
        character*3 var_2d
        character*150  directory
        character*31  ext
        character*10  units_2d
        character*125 comment_2d
        character*9 a9time
        character*4 fcst_hhmm
        character*3 c3_string
        character*4 c4_string
        character*33 c33_label
        character*16 c16_latlon
        character*11 c_pw

        integer i_overlay

        include 'icolors.inc'

        common /image/ n_image

        skewt(t_c,logp) = t_c - (logp - logp_bottom) * 32.

        n_image = 0

        write(6,*)' Input x grid point for sounding...'

        read(5,*)xsound

        write(6,*)' Input y grid point for sounding...'

        read(5,*)ysound

        tlow_c = -30.
        thigh_c = +50.

!       Get 3-D pressure field
        call get_pres_3d(i4_valid,NX_L,NY_L,NZ_L,pres_3d,istatus)
        if(istatus .ne. 1)go to 900

        isound = nint(xsound)
        jsound = nint(ysound)

        ext = 'static'

!       Get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

        var_2d='LAT'
        call read_static_grid(NX_L,NY_L,var_2d,lat,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-lat'
            return
        endif

        var_2d='LON'
        call read_static_grid(NX_L,NY_L,var_2d,lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-lon'
            return
        endif

        rlat = lat(isound,jsound)
        rlon = lon(isound,jsound)

        write(c16_latlon,101)rlat,rlon
 101    format(2f8.2)

        do iz = 1,NZ_L
            pres_1d(iz) = pres_3d(isound,jsound,iz)
            logp_1d(iz) = log(pres_1d(iz))
        enddo ! iz

!       Read appropriate 3-D fields
50      call input_product_info(i4time_ref              ! I
     1                         ,laps_cycle_time         ! I
     1                         ,3                       ! I
     1                         ,c_prodtype              ! O
     1                         ,ext                     ! O
     1                         ,directory               ! O
     1                         ,a9time                  ! O
     1                         ,fcst_hhmm               ! O
     1                         ,i4_initial              ! O
     1                         ,i4_valid                ! O
     1                         ,istatus)                ! O

!       Read Temperature
        if(c_prodtype .eq. 'A')then
            iflag_temp = 1 ! Returns Ambient Temp (K)

            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,NX_L,NY_L,NZ_L,field_3d,istatus)
            if(istatus .ne. 1)goto900

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            c33_label = 'LAPS Analysis Sounding'

        elseif(c_prodtype .eq. 'N')then
            call get_directory('balance',directory,len_dir)
            ext = 'lt1'
            directory = directory(1:len_dir)//ext(1:3)

            var_2d = 'T3'

            call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_nearest       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,NX_L,NY_L,NZ_L,field_3d,istatus)       

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            c33_label = 'LAPS Balanced Sounding'

        elseif(c_prodtype .eq. 'B' .or. c_prodtype .eq. 'F')then ! Bkg or Fcst
            var_2d = 'T3'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istatus)

            call make_fnam_lp(i4_valid,a9time,istatus)

            if(c_prodtype .eq. 'B')then
                c33_label = 'LAPS Background Sounding '//fcst_hhmm

            elseif(c_prodtype .eq. 'F')then

!               Note that is would be possible to get the c_model from the
!               directory variable with some type of basename extraction
!               subroutine

                c33_label = 'LAPS Forecast Sounding '//fcst_hhmm

            endif

        else
            write(6,*)' Unknown choice, will quit'
            go to 900

        endif

        i_overlay = i_overlay + 1

        call interp_3d(field_3d,temp_vert,xsound,xsound,ysound,ysound,       
     1                 NX_L,NY_L,NZ_L,1,NZ_L,r_missing_data)

        istat_td = 0

!       Read RH/SH 
        if(c_prodtype .eq. 'A')then ! Read RH
            var_2d = 'RHL'
            ext = 'lh3'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_rh)
            if(istat_rh .ne. 1)goto100

        elseif(c_prodtype .eq. 'N')then ! Read RH
            call get_directory('balance',directory,len_dir)
            ext = 'lh3'
            directory = directory(1:len_dir)//ext(1:3)

            var_2d = 'RHL'

            call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_nearest       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,NX_L,NY_L,NZ_L,field_3d,istat_rh)       
            if(istat_rh .ne. 1)goto100

        elseif(c_prodtype .eq. 'B' .or. c_prodtype .eq. 'F')then ! Bkg or Fcst
            var_2d = 'SH'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_sh)
            if(istat_sh .ne. 1)goto100

        else
            write(6,*)' Sorry, RH/SH not yet supported for prodtype: '
     1               ,c_prodtype
            istat_rh = 0
            goto100

        endif

        if(c_prodtype .eq. 'A' .or. c_prodtype .eq. 'N')then
            call interp_3d(field_3d,rh_vert,xsound,xsound,ysound,ysound,       
     1                     NX_L,NY_L,NZ_L,1,NZ_L,r_missing_data)

            do iz = 1,NZ_L
                t_c         = k_to_c(temp_vert(iz)) 
                td_vert(iz) = DWPT(t_c,rh_vert(iz))
            enddo ! iz

            istat_td = 1

        else
            call interp_3d(field_3d,sh_vert,xsound,xsound,ysound,ysound,       
     1                     NX_L,NY_L,NZ_L,1,NZ_L,r_missing_data)

            t_ref = -199.
            do iz = 1,NZ_L
                t_c         = k_to_c(temp_vert(iz)) 
                p_mb        = pres_1d(iz)/100.
                q_gkg       = sh_vert(iz) * 1000.
                td_vert(iz) = make_td(p_mb,t_c,q_gkg,t_ref)
            enddo ! iz

            istat_td = 1

        endif

!       Read in sfc data (pressure, temp, dewpoint, tpw)
        if(c_prodtype .eq. 'A')then ! Read LSX
            ext = 'lsx'

            var_2d = 'PS'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,pres_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto100

            var_2d = 'T'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,t_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto100

            var_2d = 'TD'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,td_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto100

            ext = 'lh4'
            var_2d = 'TPW'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,pw_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto100

        elseif(c_prodtype .eq. 'B' .or. c_prodtype .eq. 'F')then ! Bkg or Fcst
            write(6,*)' Look for Bkg/Fcst sfc fields'

            call s_len(directory,len_dir)

            write(6,*)' Orig directory = ',directory
            write(6,*)' ext = ',ext

            if(c_prodtype .eq. 'B')then
                directory = directory(1:len_dir)//'../lgb'
            else ! Fcst
                ext = 'fsf'
                call s_len(directory,len_dir)
                do i = 1,len_dir-2
                    if(directory(i:i+2) .eq. 'fua')then
                        directory(i:i+2) = ext(1:3)
                        write(6,*)'Substituted directory string'
                    endif
                enddo ! i
            endif

            write(6,*)' New directory = ',directory
            write(6,*)' ext = ',ext

            var_2d = 'PSF'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,NX_L,NY_L
     1                              ,pres_2d
     1                              ,istat_sfc)
            if(istat_sfc .ne. 1)goto100

            var_2d = 'TSF'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,NX_L,NY_L
     1                              ,t_2d
     1                              ,istat_sfc)
            if(istat_sfc .ne. 1)goto100


            var_2d = 'DSF'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,NX_L,NY_L
     1                              ,td_2d
     1                              ,istat_sfc)
            if(istat_sfc .ne. 1)goto100

        else
            istat_sfc = 0
            go to 100

        endif

!       Convert sfc variables
        p_sfc_pa = pres_2d(isound,jsound)
        logp_sfc = log(p_sfc_pa)
        t_sfc_k  = t_2d(isound,jsound)
        td_sfc_k = td_2d(isound,jsound)
        pw_sfc   = pw_2d(isound,jsound)

        write(6,*)' Sfc P = ', p_sfc_pa
        write(6,*)' TPW = ', pw_sfc

!       Read Wind (a la xsect)

 100    continue

!       Set up scaling for plots
        pbottom = 110000.
        ptop = 10000.

        logp_bottom = log(pbottom)
        logp_top = log(ptop)

        width = thigh_c - tlow_c
        height = logp_top - logp_bottom

        box_low = 0.1
        box_high = 0.9
        box_diff = box_high - box_low
        rlabel_low  = box_low  - box_diff / 8.
        rlabel_high = box_high + box_diff/ 8.

!       Scale to allow plotting within box
        call set(box_low   , box_high    , box_low , box_high,
     1           tlow_c, thigh_c, logp_bottom, logp_top, 1)

!       call setusv_dum(2hIN,7) ! Yellow
        call setusv_dum(2hIN,32) ! Yellow

!       Plot temp scale
        do it = -80, nint(thigh_c), 10
            y1 = logp_bottom
            y2 = logp_top
            x1 = skewt(float(it),y1)
            x2 = skewt(float(it),y2)
            call line(x1,y1,x2,y2) 
        enddo ! it

!       Plot theta scale
!       do it = -80, nint(thigh_c), 10
!           y1 = logp_bottom
!           y2 = logp_top
!           x1 = skewt(float(it),y1)
!           x2 = skewt(float(it),y2)
!           call line(x1,y1,x2,y2) ! or plot a curve?
!       enddo ! it

!       Plot Horizontal lines at 100mb intervals
        do ip = 100,1100,100
            p_pa = float(ip) * 100.
            x1 = tlow_c
            x2 = thigh_c
            y1 = log(p_pa)
            y2 = log(p_pa)
            call line(x1,y1,x2,y2) 
        enddo ! ip

!       Plot Rest of Box
        x1 = tlow_c
        y1 = logp_bottom
        x2 = thigh_c
        y2 = logp_bottom
        call line(x1,y1,x2,y2) 

        x1 = tlow_c
        y1 = logp_top
        x2 = thigh_c
        y2 = logp_top
        call line(x1,y1,x2,y2) 

        x1 = tlow_c
        y1 = logp_bottom
        x2 = tlow_c
        y2 = logp_top
        call line(x1,y1,x2,y2) 

        x1 = thigh_c
        y1 = logp_bottom
        x2 = thigh_c
        y2 = logp_top
        call line(x1,y1,x2,y2) 

!       Scale to allow plotting outside box
        call set(rlabel_low, rlabel_high , rlabel_low, rlabel_high,
     1           tlow_c-width/8., thigh_c+width/8., 
     1           logp_bottom-height/8., logp_top+height/8.   , 1)

        call setusv_dum(2hIN,32) ! Yellow

!       Plot Pressure labels at 100mb intervals
        do ip = 100,1100,100
            p_pa = float(ip) * 100.
            y = log(p_pa) 

!           Pressure
            x = tlow_c - 3. * .8/box_diff
            write(c4_string,2014)ip
            call pwrity (x, y, c4_string, 4, 1, 0, 0)
2014        format(i4)
        enddo ! ip

!       Plot Temperature Labels
        do it = nint(tlow_c), nint(thigh_c), 10
            t_c = float(it)
            write(c3_string,2013)it
2013        format(i3)
            x = t_c
            y = logp_bottom + .04 * .8/box_diff
            call pwrity (x, y, c3_string, 3, 1, 0, 0)
        enddo

!       Plot lat/lon info
        x = thigh_c - 10.
        y = logp_top + height*.05
        call pwrity (x, y, c16_latlon, 16, 1, 0, 0)

!       Do sounding plots on non-skew T, log P scale
        call set(box_low   , box_high    , box_low , box_high,
     1           tlow_c, thigh_c, logp_bottom, logp_top, 1)

        call setusv_dum(2hIN,icolors(i_overlay))

!       Plot temp and dewpoint sounding
        do iz = 2,NZ_L

            if(istat_sfc .eq. 0      .OR.
     1         pres_1d(iz-1) .le. p_sfc_pa          )then  ! above sfc or no 

                call setusv_dum(2hIN,icolors(i_overlay))   ! sfc data   
                                                                    
                y1 = logp_1d(iz-1)
                y2 = logp_1d(iz)

                x1 = skewt(k_to_c(temp_vert(iz-1)),y1)
                x2 = skewt(k_to_c(temp_vert(iz)),y2)
                call line(x1,y1,x2,y2) 

                if(istat_td .eq. 1)then
                    x1 = skewt(td_vert(iz-1),y1)
                    x2 = skewt(td_vert(iz),y2)
                    CALL GSLN (3)
                    call line(x1,y1,x2,y2) 
                    CALL GSLN (1)
                endif

            elseif(pres_1d(iz) .gt. p_sfc_pa)then          ! below sfc
                call setusv_dum(2hIN,34)

                y1 = logp_1d(iz-1)
                y2 = logp_1d(iz)

                x1 = skewt(k_to_c(temp_vert(iz-1)),y1)
                x2 = skewt(k_to_c(temp_vert(iz)),y2)
                call line(x1,y1,x2,y2) 

                if(istat_td .eq. 1)then
                    x1 = skewt(td_vert(iz-1),y1)
                    x2 = skewt(td_vert(iz),y2)
                    CALL GSLN (3)
                    call line(x1,y1,x2,y2) 
                    CALL GSLN (1)
                endif

            else                                           ! straddles the sfc
!               Plot line below the sfc              
                call setusv_dum(2hIN,34)

                y1 = logp_1d(iz-1)
                y2 = logp_sfc

                x1 = skewt(k_to_c(temp_vert(iz-1)),y1)
                x2 = skewt(k_to_c(t_sfc_k),y2)
                call line(x1,y1,x2,y2) 

                if(istat_td .eq. 1)then
                    x1 = skewt(td_vert(iz-1),y1)
                    x2 = skewt(k_to_c(td_sfc_k),y2)
                    CALL GSLN (3)
                    call line(x1,y1,x2,y2) 
                    CALL GSLN (1)
                endif

!               Plot line above the sfc              
                call setusv_dum(2hIN,icolors(i_overlay))

                y1 = logp_sfc
                y2 = logp_1d(iz)

                x1 = skewt(k_to_c(t_sfc_k),y1)
                x2 = skewt(k_to_c(temp_vert(iz)),y2)
                call line(x1,y1,x2,y2) 

                if(istat_td .eq. 1)then
                    x1 = skewt(k_to_c(td_sfc_k),y1)
                    x2 = skewt(td_vert(iz),y2)
                    CALL GSLN (3)
                    call line(x1,y1,x2,y2) 
                    CALL GSLN (1)
                endif

            endif

        enddo ! iz


!       Plot time/source label
        if(box_low .eq. 0.1)then
            call set(0., 1., 0., 1.
     1              ,0., 1., 0., 1., 1)
            call write_label_lplot(100,94,c33_label,a9time
     1                            ,plot_parms,namelist_parms
     1                            ,i_overlay,'sound')
        elseif(box_low .eq. 0.2)then
            call set(0., 1., 0., 1.
     1              ,0., 1., 0., 1., 1)
            call write_label_lplot(100,94,c33_label,a9time
     1                            ,plot_parms,namelist_parms
     1                            ,i_overlay,'sound')
        endif

!       Plot TPW value
        ix = 800
        iy = 180
        rsize = .010
        write(c_pw,890)pw_sfc*100. ! convert M to CM
 890    format('IWV = ',f5.2)
        CALL PCHIQU (cpux(ix),cpux(iy),c_pw,rsize,0,-1.0)

        write(6,*)' Sounding has been plotted...'

        go to 50

 900    continue

        call sflush

        return
        end
