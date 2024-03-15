
        subroutine plot_sounding(i4time_ref,lun,nx_l,ny_l,nz_l
     1                          ,r_missing_data,laps_cycle_time,maxstns
     1                          ,i_overlay,plot_parms,namelist_parms
     1                          ,l_plotobs)       

        use mem_namelist, only: max_snd_grid, max_snd_levels

        include 'lapsplot.inc'

        real pres_3d(nx_l,ny_l,nz_l)
        real field_3d(nx_l,ny_l,nz_l)
!       real rh_3d(nx_l,ny_l,nz_l)

        real pres_2d(nx_l,ny_l)
        real t_2d(nx_l,ny_l)
        real td_2d(nx_l,ny_l)
        real u_2d(nx_l,ny_l)
        real v_2d(nx_l,ny_l)
        real pw_2d(nx_l,ny_l)
        real cape_2d(nx_l,ny_l)
        real lil_2d(nx_l,ny_l)
        real lic_2d(nx_l,ny_l)
        real lat(nx_l,ny_l)
        real lon(nx_l,ny_l)
        real topo(nx_l,ny_l)

        real temp_vert(max_snd_levels)
        real ht_vert(max_snd_levels)
        real u_vert(max_snd_levels)
        real v_vert(max_snd_levels)
        real rh_vert(max_snd_levels)
        real sh_vert(max_snd_levels) ! dimensionless
        real td_vert(max_snd_levels)
        real lwc_vert(max_snd_levels)
        real ice_vert(max_snd_levels)
        real rain_vert(max_snd_levels)
        real snow_vert(max_snd_levels)
        real pice_vert(max_snd_levels)

        real pres_1d(max_snd_levels)
        real logp_1d(max_snd_levels), logp_bottom, logp_top, logp
     1                              , logp_sfc
        real lil_sfc, lic_sfc, lil_cpt, lic_cpt
  
        real k_to_c, make_td, make_ssh

        character*1 c_prodtype, c_plotobs
        character*3 var_2d
        character*150  directory
        character*31  ext
        character*10  units_2d
        character*125 comment_2d
        character*40 c_model
        character*9 a9time
        character*5 fcst_hhmm
        character*3 c3_string
        character*4 c4_string
        character*40 c_label
        character*16 c16_latlon
        character*11 c_pw
        character*11 c_cape
        character*20 c20_x, c20_y
        character*255 new_dataroot
        logical l_latlon, l_parse, l_plotobs

        integer i_overlay

        include 'icolors.inc'

!       sounding observation declarations
        real ob_pr_ht_obs(max_snd_grid,max_snd_levels)
        real ob_pr_pr_obs(max_snd_grid,max_snd_levels)
        real ob_pr_u_obs(max_snd_grid,max_snd_levels) 
        real ob_pr_v_obs(max_snd_grid,max_snd_levels)
        real ob_pr_t_obs(max_snd_grid,max_snd_levels)
        real ob_pr_td_obs(max_snd_grid,max_snd_levels)
        real lat_pr(max_snd_grid)
        real lon_pr(max_snd_grid)
        real elev_pr(max_snd_grid)
        integer nlevels_obs_pr(max_snd_grid)
        integer i4time_ob_pr(max_snd_grid)
        character*5 c5_name, c5_name_a(max_snd_grid), c5_name_min
        character*8 obstype(max_snd_grid)

        include 'lapsparms.for'
        integer mxl
        parameter (mxl=max_lvls+1) ! number of 3d levels plus the sfc level

        common/indx/ p(mxl),t(mxl),td(mxl),ht(mxl),pbecr(20,4)
     1  ,tdfcr(20,2),vel(20)
     1  ,temdif(mxl),partem(mxl),pbe(mxl)
     #  ,dd85,ff85,dd50,ff50
        real lcl,li,k_index

        common /image/ n_image

        skewt(t_c,logp) = t_c - (logp - logp_bottom) * 32.

        pbe = 0. ! initialize

        nsmooth = plot_parms%obs_size
        if(nsmooth .ne. 3)then
            nsmooth = 1
        endif

        write(6,*)
        write(6,*)' subroutine plot_sounding: nsmooth is ',nsmooth

        itd = 2 ! dashed dewpoint lines

        n_image = 0

        ext = 'static'

!       get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

        var_2d='lat'
        call read_static_grid(nx_l,ny_l,var_2d,lat,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading laps static-lat'
            return
        endif

        var_2d='lon'
        call read_static_grid(nx_l,ny_l,var_2d,lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading laps static-lon'
            return
        endif

        var_2d='avg'
        call read_static_grid(nx_l,ny_l,var_2d,topo,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading laps static-topo'
            return
        endif

 80     write(6,*)
        write(6,*)' input x grid point (or latitude) for sounding...'
        read(5,*)c20_x

        call s_len(c20_x,lenx)
        l_latlon = l_parse(c20_x(1:lenx),'l')

        if(l_latlon)then ! x value was flagged as latitude with "l" at the end 
            read(c20_x(1:lenx-1),*)soundlat

            write(6,*)' input longitude for sounding...'       
            read(5,*)c20_y
            call s_len(c20_y,leny)
            if(l_parse(c20_y(1:leny),'l'))then
                read(c20_y(1:leny-1),*)soundlon
            else
                read(c20_y(1:leny),*)soundlon
            endif

            call latlon_to_rlapsgrid(soundlat,soundlon,lat,lon
     1                              ,nx_l,ny_l,xsound,ysound,istatus)       

            if(istatus .ne. 1)then
                write(6,*)' station is outside domain - try again...'
                return
            endif

            if(xsound .lt. 1. .or. xsound .gt. float(nx_l) .or.
     1         ysound .lt. 1. .or. ysound .gt. float(ny_l)    )then
                write(6,*)' station is outside domain - try again...'
                return
            endif

        else
            read(c20_x(1:lenx),*)xsound
            write(6,*)' input y grid point for sounding...'
            read(5,*)ysound
        endif

 40     write(6,*)' enter c_plotobs'
        read(5,*)c_plotobs
        if(c_plotobs .eq. '1')then
            l_plotobs = .true.
        elseif(c_plotobs .eq. '0')then
            l_plotobs = .false.
        else
            write(6,*)' unknown c_plotobs, will quit ',c_plotobs
            go to 900
        endif

        if(.true.)then ! force config with new dataroot
            write(6,*)' enter new dataroot:'
            read(5,17)new_dataroot
 17         format(a)
            call s_len(new_dataroot,lenroot)

            if(new_dataroot(1:1) .eq. 'q')then
                write(6,*)' unknown dataroot, will quit'
                go to 900
            else
                write(6,*)' new dataroot is ',new_dataroot(1:lenroot) 
            endif

            call force_get_laps_config(new_dataroot(1:lenroot),istatus)
            if(istatus .ne. 1)then
                write(6,*)' bad status returned from force_laps_config'
                return
            else
                write(6,*)' forced config to ',new_dataroot(1:lenroot)
            endif
        endif

        tlow_c = -30.
        thigh_c = +50.

!       get 3-d pressure field
        call get_pres_3d(i4_valid,nx_l,ny_l,nz_l,pres_3d,istatus)
        if(istatus .ne. 1)go to 900

        isound = nint(xsound)
        jsound = nint(ysound)

        rlat = lat(isound,jsound)
        rlon = lon(isound,jsound)

        write(c16_latlon,101)rlat,rlon
 101    format(2f8.2)

        do iz = 1,nz_l
            pres_1d(iz) = pres_3d(isound,jsound,iz)
            logp_1d(iz) = log(pres_1d(iz))
        enddo ! iz

        if(l_plotobs .eqv. .false.)then

          n_lvls_snd = nz_l

!         read appropriate 3-d fields
50        call input_product_info(i4time_ref            ! i
     1                         ,laps_cycle_time         ! i
     1                         ,3                       ! i
     1                         ,c_prodtype              ! o
     1                         ,ext                     ! o
     1                         ,directory               ! o
     1                         ,a9time                  ! o
     1                         ,fcst_hhmm               ! o
     1                         ,i4_initial              ! o
     1                         ,i4_valid                ! o
     1                         ,istatus)                ! o

!         read temperature
          if(c_prodtype .eq. 'a')then
            iflag_temp = 1 ! returns ambient temp (k)

            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,nx_l,ny_l,nz_l,field_3d,istatus)
            if(istatus .ne. 1)goto900

            call make_fnam_lp(i4time_nearest,a9time,istatus)

            if(nsmooth .eq. 3)then
                c_label = 'analysis sounding (3x3 smoothing)'
            else
                c_label = 'analysis sounding'
            endif

          elseif(c_prodtype .eq. 'n')then
            call get_directory('balance',directory,len_dir)
            ext = 'lt1'
            directory = directory(1:len_dir)//ext(1:3)

            var_2d = 't3'

            call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_nearest       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,nx_l,ny_l,nz_l,field_3d,istatus)       

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            c_label = 'balanced sounding'

          elseif(c_prodtype .eq. 'b' .or. c_prodtype .eq. 'f')then ! bkg or fcst
            var_2d = 't3'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istatus)
            if(istatus .ne. 1)goto900

            call make_fnam_lp(i4_valid,a9time,istatus)

            if(c_prodtype .eq. 'b')then
                c_label = 'background sounding '//fcst_hhmm

            elseif(c_prodtype .eq. 'f')then

                call directory_to_cmodel(directory,c_model)

                c_label = 'forecast sounding '//fcst_hhmm//' '
     1                                               //trim(c_model)
            endif

          else
            write(6,*)' unknown choice, will quit'
            go to 900

          endif

          i_overlay = i_overlay + 1

          if(nsmooth .gt. 1)then
              call smooth_box_3d(field_3d,nx_l,ny_l,nz_l,nsmooth)
          endif

          call interp_3d(field_3d,temp_vert,xsound,xsound,ysound,ysound,       
     1                   nx_l,ny_l,nz_l,1,nz_l,r_missing_data)

        else ! plotobs is true
          iflag_temp = 2 ! returns height

          call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                    ,nx_l,ny_l,nz_l,field_3d,istatus)
          if(istatus .ne. 1)goto900

          i4time_snd = (i4time_ref/laps_cycle_time) * laps_cycle_time
          lun_snd = 12
          mode = 3 ! key levels off of any data
          ext = 'snd'
          write(6,*)' call read_snd_data for i4time: ',i4time_snd
          n_profiles = 0

          call read_snd_data(lun_snd,i4time_snd,ext                     ! i
     1                         ,max_snd_grid,max_snd_levels             ! i
     1                         ,lat,lon,topo,nx_l,ny_l,nz_l             ! i
     1                         ,field_3d,.true.                         ! i
     1                         ,mode                                    ! i
     1                         ,n_profiles                              ! i/o
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr    ! o
     1                         ,c5_name_a,i4time_ob_pr,obstype          ! o
     1                         ,ob_pr_ht_obs,ob_pr_pr_obs               ! o
     1                         ,ob_pr_u_obs,ob_pr_v_obs                 ! o
     1                         ,ob_pr_t_obs,ob_pr_td_obs                ! o
     1                         ,istatus)                                ! o

          distmin = 999999.
          write(6,*)' n_profiles is ',n_profiles
          do is = 1,n_profiles
            call latlon_to_rlapsgrid(lat_pr(is),lon_pr(is),lat,lon
     1                              ,nx_l,ny_l,xob,yob,istatus)       
            dist = sqrt((xob-xsound)**2 + (yob-ysound)**2)
            write(6,*)is,lat_pr(is),lon_pr(is),c5_name_a(is),dist
            if(dist .lt. distmin)then
              distmin = dist
              ismin = is
              c5_name_min = c5_name_a(is)
            endif
          enddo ! is 
     
          if(n_profiles .gt. 0)then
              write(6,*)' minimum distance sounding ob is ',ismin
     1                           ,nlevels_obs_pr(ismin),c5_name_min
          endif

          i_overlay = i_overlay + 1
          write(6,*)' sounding ob data, nlevels is: '
     1                                ,nlevels_obs_pr(ismin)
          ll = 0

          std_corr = 1.0

          write(6,*)
     1'   il   ll    htob     prob      prstd     prlvl        t     td'

          do il = 1,nlevels_obs_pr(ismin)
              if(ob_pr_ht_obs(ismin,il) .ne. r_missing_data
     1      .or. ob_pr_pr_obs(ismin,il) .ne. r_missing_data)then

                if(ob_pr_ht_obs(ismin,il) .ne. r_missing_data     ! only ht
     1       .and. ob_pr_pr_obs(ismin,il) .eq. r_missing_data)then        
                  pres_std_pa = ztopsa(ob_pr_ht_obs(ismin,il))*100.
                  pres_lvl_pa = pres_std_pa * std_corr              

                elseif(ob_pr_ht_obs(ismin,il) .ne. r_missing_data ! only pres
     1       .and.     ob_pr_pr_obs(ismin,il) .eq. r_missing_data)then        
                  pres_lvl_pa = ob_pr_pr_obs(ismin,il)*100.

                else                                              ! both present
!                 recalculate correction to standard atmosphere
                  pres_lvl_pa = ob_pr_pr_obs(ismin,il)*100.
                  pres_std_pa = ztopsa(ob_pr_ht_obs(ismin,il))*100.
                  std_corr = pres_lvl_pa / pres_std_pa

                endif
              
                if(ob_pr_t_obs(ismin,il) .ne. r_missing_data)then
                  ll = ll + 1
                  pres_1d(ll) = pres_lvl_pa
                  logp_1d(ll) = log(pres_lvl_pa)
                  temp_vert(ll) = c_to_k(ob_pr_t_obs(ismin,il))
                  td_vert(ll) = ob_pr_td_obs(ismin,il)
                  istat_td = 1
                  write(6,55)il,ll,ob_pr_ht_obs(ismin,il)               
     1                     ,ob_pr_pr_obs(ismin,il) 
     1                     ,pres_std_pa
     1                     ,pres_lvl_pa
     1                     ,ob_pr_t_obs(ismin,il) 
     1                     ,ob_pr_td_obs(ismin,il) 
55                format(1x,2i5,6f10.3)
                else
                  write(6,*)' temp info missing for level ',il
                endif

              else
                write(6,*)' height info missing for level ',il
              endif

          enddo ! il

          n_lvls_snd = ll ! nlevels_obs_pr(ismin)

!         c_label = 'sounding ob '//c5_name_min
!         write(c_label,60)c5_name_min,lat_pr(ismin),lon_pr(ismin)
          write(c_label,60)trim(obstype(ismin)),trim(c5_name_min)
     1                    ,lat_pr(ismin),lon_pr(ismin)
 60       format(a,' ob ',a5,2f9.2)

          p_sfc_pa = pres_1d(1) ! 110000.
          call make_fnam_lp(i4time_ob_pr(ismin),a9time,istatus)

          pw_sfc = r_missing_data
          cape_sfc = r_missing_data

        endif ! l_plotobs

        if(l_plotobs .eqv. .false.)then

!         read height
          if(c_prodtype .eq. 'a')then
            iflag_temp = 2 ! returns height

            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,nx_l,ny_l,nz_l,field_3d,istatus)
            if(istatus .ne. 1)goto300

          elseif(c_prodtype .eq. 'n')then
            call get_directory('balance',directory,len_dir)
            ext = 'lt1'
            directory = directory(1:len_dir)//ext(1:3)

            var_2d = 'ht'

            call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_nearest       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,nx_l,ny_l,nz_l,field_3d,istatus)       

          elseif(c_prodtype .eq. 'b' .or. c_prodtype .eq. 'f')then ! bkg or fcst
            var_2d = 'ht'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istatus)
            if(istatus .ne. 1)goto300

          else
            write(6,*)' unknown choice, will quit'
            go to 900

          endif

          call interp_3d(field_3d,ht_vert,xsound,xsound,ysound,ysound,       
     1                   nx_l,ny_l,nz_l,1,nz_l,r_missing_data)

!         read rh/sh 
 300      istat_td = 0
          if(c_prodtype .eq. 'a')then ! read rh
            var_2d = 'rhl'
            ext = 'lh3'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,nx_l,ny_l,nz_l       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_rh)
            if(istat_rh .ne. 1)goto1000

          elseif(c_prodtype .eq. 'n')then ! read rh
            call get_directory('balance',directory,len_dir)
            ext = 'lh3'
            directory = directory(1:len_dir)//ext(1:3)

            var_2d = 'rhl'

            call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_nearest       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,nx_l,ny_l,nz_l,field_3d,istat_rh)       
            if(istat_rh .ne. 1)goto1000

          elseif(c_prodtype .eq. 'b' .or. c_prodtype .eq. 'f')then ! bkg or fcst
            var_2d = 'sh'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_sh)
            if(istat_sh .ne. 1)goto1000

          else
            write(6,*)' sorry, rh/sh not yet supported for prodtype: '
     1               ,c_prodtype
            istat_rh = 0
            goto1000

          endif

          if(c_prodtype .eq. 'a' .or. c_prodtype .eq. 'n')then
            if(nsmooth .gt. 1)then
                call smooth_box_3d(field_3d,nx_l,ny_l,nz_l,nsmooth)
            endif
            call interp_3d(field_3d,rh_vert,xsound,xsound,ysound,ysound,       
     1                     nx_l,ny_l,nz_l,1,nz_l,r_missing_data)

            do iz = 1,nz_l
                t_c         = k_to_c(temp_vert(iz)) 
                td_vert(iz) = dwpt(t_c,rh_vert(iz))
                p_mb        = pres_1d(iz)/100.
                sh_vert(iz) = make_ssh(p_mb,t_c,rh_vert(iz),-132.)
            enddo ! iz

            istat_td = 1

          else
            if(nsmooth .gt. 1)then
                call smooth_box_3d(field_3d,nx_l,ny_l,nz_l,nsmooth)
            endif
            call interp_3d(field_3d,sh_vert,xsound,xsound,ysound,ysound,       
     1                     nx_l,ny_l,nz_l,1,nz_l,r_missing_data)

            t_ref = -199.
            do iz = 1,nz_l
                t_c         = k_to_c(temp_vert(iz)) 
                p_mb        = pres_1d(iz)/100.
                q_gkg       = sh_vert(iz) * 1000.
                td_vert(iz) = make_td(p_mb,t_c,q_gkg,t_ref)
                rh_vert(iz) = humidity(t_c,td_vert(iz))
            enddo ! iz

            istat_td = 1

          endif

!         read cloud liquid
          istat_lwc = 0
          if(c_prodtype .eq. 'a')then ! read cloud liquid
            var_2d = 'lwc'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,nx_l,ny_l,nz_l       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_lwc)
          elseif(c_prodtype .eq. 'f')then 
            var_2d = 'lwc'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_lwc)
!           if(istat_lwc .ne. 1)goto1000
          endif

          if(istat_lwc .eq. 1)then
            call interp_3d(field_3d,lwc_vert,xsound,xsound
     1                    ,ysound,ysound,nx_l,ny_l,nz_l,1,nz_l
     1                    ,r_missing_data)
          else
            lwc_vert = -999.
          endif

!       read cloud ice
          istat_ice = 0
          if(c_prodtype .eq. 'a')then ! read cloud ice
            var_2d = 'ice'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,nx_l,ny_l,nz_l       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_ice)
          elseif(c_prodtype .eq. 'f')then 
            var_2d = 'ice'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_ice)
!           if(istat_ice .ne. 1)goto1000
          endif

          if(istat_ice .eq. 1)then
            call interp_3d(field_3d,ice_vert,xsound,xsound
     1                    ,ysound,ysound,nx_l,ny_l,nz_l,1,nz_l
     1                    ,r_missing_data)
          else
            ice_vert = -999.
          endif

!         read precipitating rain
          istat_rain = 0
          if(c_prodtype .eq. 'a')then ! read precipitating rain
            var_2d = 'rai'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,nx_l,ny_l,nz_l       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_rain)
          elseif(c_prodtype .eq. 'f')then 
            var_2d = 'rai'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_rain)
!           if(istat_rain .ne. 1)goto1000
          endif

          if(istat_rain .eq. 1)then
            call interp_3d(field_3d,rain_vert,xsound,xsound
     1                    ,ysound,ysound,nx_l,ny_l,nz_l,1,nz_l
     1                    ,r_missing_data)
          else
            rain_vert = -999.
          endif

!         read precipitating snow
          istat_snow = 0
          if(c_prodtype .eq. 'a')then ! read precipitating snow
            var_2d = 'sno'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,nx_l,ny_l,nz_l       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_snow)
          elseif(c_prodtype .eq. 'f')then 
            var_2d = 'sno'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_snow)
!           if(istat_snow .ne. 1)goto1000
          endif

          if(istat_snow .eq. 1)then
            call interp_3d(field_3d,snow_vert,xsound,xsound
     1                    ,ysound,ysound,nx_l,ny_l,nz_l,1,nz_l
     1                    ,r_missing_data)
          else
            snow_vert = -999.
          endif

!         read precipitating ice
          istat_pice = 0
          if(c_prodtype .eq. 'a')then ! read precipitating ice
            var_2d = 'pic'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,nx_l,ny_l,nz_l       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_pice)
          elseif(c_prodtype .eq. 'f')then 
            var_2d = 'pic'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_pice)
!           if(istat_pice .ne. 1)goto1000
          endif

          if(istat_pice .eq. 1)then
            call interp_3d(field_3d,pice_vert,xsound,xsound
     1                    ,ysound,ysound,nx_l,ny_l,nz_l,1,nz_l
     1                    ,r_missing_data)
          else
            pice_vert = -999.
          endif

!         read 3-d u wind component
          istat_u = 0
          if(c_prodtype .eq. 'a')then ! read 3-d u wind
            var_2d = 'u3'
            ext = 'lw3'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,nx_l,ny_l,nz_l       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_u)
          elseif(c_prodtype .eq. 'f')then
            var_2d = 'u3'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_u)
!           if(istat_u .ne. 1)goto1000
          endif

          if(istat_u .eq. 1)then
            call interp_3d(field_3d,u_vert,xsound,xsound
     1                    ,ysound,ysound,nx_l,ny_l,nz_l,1,nz_l
     1                    ,r_missing_data)
          else
            u_vert = -999.
          endif

!       read 3-d v wind component
          istat_v = 0
          if(c_prodtype .eq. 'a')then ! read 3-d v wind
            var_2d = 'v3'
            ext = 'lw3'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,nx_l,ny_l,nz_l       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_v)
          elseif(c_prodtype .eq. 'f')then
            var_2d = 'v3'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_v)
!           if(istat_v .ne. 1)goto1000
          endif

          if(istat_v .eq. 1)then
            call interp_3d(field_3d,v_vert,xsound,xsound
     1                    ,ysound,ysound,nx_l,ny_l,nz_l,1,nz_l
     1                    ,r_missing_data)
          else
            u_vert = -999.
          endif

!       read in sfc data (pressure, temp, dewpoint, u, v, tpw, cape)
          if(c_prodtype .eq. 'a')then ! read lsx
            ext = 'lsx'

            var_2d = 'ps'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,pres_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 't'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,t_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'td'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,td_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'u'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,u_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'v'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,v_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            ext = 'lh4'
            var_2d = 'tpw'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,pw_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            ext = 'lst'
            var_2d = 'pbe'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,cape_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            ext = 'lil'
            var_2d = 'lil'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,lil_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'lic'
            call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,lic_2d,0,istat_lic)
            if(istat_lic .ne. 1)lic_2d = r_missing_data

          elseif(c_prodtype .eq. 'b' .or. c_prodtype .eq. 'f')then ! bkg or fcst
            write(6,*)' look for bkg/fcst sfc fields'

            call s_len(directory,len_dir)

            write(6,*)' orig directory = ',directory
            write(6,*)' ext = ',ext

            if(c_prodtype .eq. 'b')then
                directory = directory(1:len_dir)//'../lgb'
            else ! fcst
                ext = 'fsf'
                call s_len(directory,len_dir)
                do i = 1,len_dir-2
                    if(directory(i:i+2) .eq. 'fua')then
                        directory(i:i+2) = ext(1:3)
                        write(6,*)'substituted directory string'
                    endif
                enddo ! i
            endif

            write(6,*)' new directory = ',directory
            write(6,*)' ext = ',ext

            var_2d = 'psf'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,nx_l,ny_l
     1                              ,pres_2d
     1                              ,istat_sfc)
!           if(istat_sfc .ne. 1)goto1000

            var_2d = 'tsf'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,nx_l,ny_l
     1                              ,t_2d
     1                              ,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'dsf'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,nx_l,ny_l
     1                              ,td_2d
     1                              ,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'tpw'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,nx_l,ny_l
     1                              ,pw_2d
     1                              ,istat_tpw)
            if(istat_tpw .ne. 1)pw_2d = r_missing_data

            var_2d = 'pbe'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,nx_l,ny_l
     1                              ,cape_2d
     1                              ,istat_pbe)
            if(istat_pbe .ne. 1)cape_2d = r_missing_data

            var_2d = 'lil'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,nx_l,ny_l
     1                              ,lil_2d
     1                              ,istat_lil)
            if(istat_lil .ne. 1)lil_2d = r_missing_data

            var_2d = 'lic'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,nx_l,ny_l
     1                              ,lic_2d
     1                              ,istat_lic)
            if(istat_lic .ne. 1)lic_2d = r_missing_data

          else
            istat_sfc = 0
            go to 1000

          endif

!         convert sfc variables
          p_sfc_pa = pres_2d(isound,jsound)
          logp_sfc = log(p_sfc_pa)
          t_sfc_k  = t_2d(isound,jsound)
          td_sfc_k = td_2d(isound,jsound)
          pw_sfc   = pw_2d(isound,jsound)
          cape_sfc = cape_2d(isound,jsound)
          lil_sfc  = lil_2d(isound,jsound)
          lic_sfc  = lic_2d(isound,jsound)
          topo_sfc = topo(isound,jsound)

        endif ! l_plotobs is false

        write(6,*)
        write(6,102)c16_latlon,isound,jsound,a9time
 102    format(' ascii sounding data at lat,lon =',a16
     1        ,'   i,j = ',2i5,5x,' time: ',a9)     
        write(6,*)
        write(6,*)' lvl    p(mb)     ht(m)      t(c)'//
     1            '      td(c)      rh(%)    cld liq  cld ice'//   
     1            '    rain      snow      pcpice   u (m/s)   v (m/s)'       
        write(6,*)'                                                '
     1          //'g/m**3    g/m**3    g/m**3    g/m**3    g/m**3'

        do iz = 1,n_lvls_snd
            if(l_plotobs .eqv. .false.)then
              iz_test = min((iz+1),n_lvls_snd)
            else
              iz_test = iz
            endif
            if(pres_1d(iz_test)/100. .le. p_sfc_pa/100.)then
                write(6,1)iz,
     1                pres_1d(iz)/100.,
     1                ht_vert(iz),
     1                k_to_c(temp_vert(iz)), 
     1                td_vert(iz),
     1                rh_vert(iz),
     1                lwc_vert(iz)*1000.,
     1                ice_vert(iz)*1000.,
     1                rain_vert(iz)*1000.,
     1                snow_vert(iz)*1000.,
     1                pice_vert(iz)*1000.,
     1                u_vert(iz),
     1                v_vert(iz)
 1              format(i4,5f10.2,5f10.3,2f10.1)
            endif
        enddo ! iz
        write(6,*)

!       compute integrated cloud liquid/ice
        call integrate_slwc(lwc_vert*1000.,ht_vert,1,1,nz_l,lil_cpt)
        call integrate_slwc(ice_vert*1000.,ht_vert,1,1,nz_l,lic_cpt)

        td_sfc_c = k_to_c(td_sfc_k)
        p_sfc_mb = p_sfc_pa/100.
        sh_sfc = make_ssh(p_sfc_mb,td_sfc_c,1.0,-132.) / 1000.

        call integrate_tpw(sh_vert,sh_sfc,pres_1d,p_sfc_pa
     1                    ,1,1,nz_l,pw_cpt)

!       compute cape and other indices
        twet_snow = +1.3 ! also in deriv.nl
        p(1) = p_sfc_mb                                  
        t(1) = k_to_c(t_sfc_k)
        td(1) = td_sfc_c
        ht(1) = topo_sfc 

        do k = 1,nz_l
            if(pres_1d(k) .lt. p_sfc_pa)then ! first level above sfc
                n_first_level = k
                goto100
            endif
        enddo
100     write(6,*)' n_first_level = ',n_first_level

        n = 1
        do k = n_first_level,nz_l
            n = n + 1
            p(n) = pres_1d(k) / 100.       ! pa to mb
            t(n) = k_to_c(temp_vert(k))
            td(n)= td_vert(k)
            ht(n)= ht_vert(k)     
        enddo ! k

        i4_elapsed = ishow_timer()

        io = 2
        idebug = 2
        nlevel = n
        call sindx(nlevel,li,si,bli,tt,sweat
     1                ,twet_snow                                        ! i
     1                ,hwb0,hwb_snow                                    ! o
     1                ,plcl,lcl,ccl
     1                ,tconv,io
     1                ,icp,ict,k_index,tmax,pbeneg,pbepos,t500,pbli
     1                ,velneg,water,ihour,idebug,istatus)
        if(istatus .ne. 1)then
            write(6,*)' warning: bad istatus returned from sindx'
!           return
        endif

        i4_elapsed = ishow_timer()

        write(6,*)' sfc p (mb) = ', p_sfc_mb     
     1           ,' terrain height (m)= ',topo_sfc 
        write(6,*)' sfc t (c) = ', k_to_c(t_sfc_k)
        write(6,*)' sfc td (c) = ', td_sfc_c
        write(6,*)' tpw (cm) [read in/computed] = ', 
     1                            pw_sfc*100.,pw_cpt*100.
        write(6,*)' sbcape (j/kg) = ', cape_sfc,pbepos
        write(6,*)' integrated cloud water (mm) [read in/computed] = ',
     1                            lil_sfc*1000.,lil_cpt*1000.
        write(6,*)' integrated cloud ice (mm) [read in/computed] = ',
     1                            lic_sfc*1000.,lic_cpt*1000.
        write(6,*)

!       read wind (a la xsect)

 1000   continue

!       set up scaling for plots
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

!       scale to allow plotting within box
        call set(box_low   , box_high    , box_low , box_high,
     1           tlow_c, thigh_c, logp_bottom, logp_top, 1)

!       call setusv_dum(2hin,7) ! yellow
        call setusv_dum(2hin,32) ! yellow

!       plot temp scale
        do it = -80, nint(thigh_c), 10
            y1 = logp_bottom
            y2 = logp_top
            x1 = skewt(float(it),y1)
            x2 = skewt(float(it),y2)
            call line(x1,y1,x2,y2) 
        enddo ! it

!       plot mixing ratio scale
        call setusv_dum(2hin,23) ! medium gray
        do im = -1,5                      
          w_gkg = 2.**im                  
          do ip = 1000,100,-100
            p1 = ip + 100
            p2 = ip
            y1 = log(p1*100.)    
            y2 = log(p2*100.)  
            tdsat1 = tmr(w_gkg,p1) 
            tdsat2 = tmr(w_gkg,p2)
            x1 = skewt(tdsat1,y1)
            x2 = skewt(tdsat2,y2)
            call line(x1,y1,x2,y2) ! plot a curve
          enddo ! ip

!         add label for this mixing ratio value
          x = tmr(w_gkg,1060.) - 0.4 
          y = log(106000.) ! + .04 * .8/box_diff   
          if(im .ge. 0)then
              write(c3_string,2013)nint(w_gkg)
          else
              write(c3_string,2012)w_gkg
2012          format(f3.1)
          endif
          call pwrity (x, y, c3_string, 3, 0, 0, 0)
        enddo ! im

!       plot theta(e) scale
        call setusv_dum(2hin,32) ! yellow
!       do it = -80, nint(thigh_c), 10
!           y1 = logp_bottom
!           y2 = logp_top
!           x1 = skewt(float(it),y1)
!           x2 = skewt(float(it),y2)
!           call line(x1,y1,x2,y2) ! or plot a curve?
!       enddo ! it

!       plot horizontal lines at 100mb intervals
        do ip = 100,1100,100
            p_pa = float(ip) * 100.
            x1 = tlow_c
            x2 = thigh_c
            y1 = log(p_pa)
            y2 = log(p_pa)
            call line(x1,y1,x2,y2) 
        enddo ! ip

!       plot rest of box
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

!       scale to allow plotting outside box
        call set(rlabel_low, rlabel_high , rlabel_low, rlabel_high,
     1           tlow_c-width/8., thigh_c+width/8., 
     1           logp_bottom-height/8., logp_top+height/8.   , 1)

        call setusv_dum(2hin,32) ! yellow

!       plot pressure labels at 100mb intervals
        do ip = 100,1100,100
            p_pa = float(ip) * 100.
            y = log(p_pa) 

!           pressure
            x = tlow_c - 3. * .8/box_diff
            write(c4_string,2014)ip
            call pwrity (x, y, c4_string, 4, 1, 0, 0)
2014        format(i4)
        enddo ! ip

!       plot temperature labels
        do it = nint(tlow_c), nint(thigh_c), 10
            t_c = float(it)
            write(c3_string,2013)it
2013        format(i3)
            x = t_c
            y = logp_bottom + .04 * .8/box_diff
            call pwrity (x, y, c3_string, 3, 1, 0, 0)
        enddo

!       plot lat/lon info
        x = thigh_c - 10.
        y = logp_top + height*.05
        call pwrity (x, y, c16_latlon, 16, 1, 0, 0)

!       do sounding plots on non-skew t, log p scale
        call set(box_low   , box_high    , box_low , box_high,
     1           tlow_c, thigh_c, logp_bottom, logp_top, 1)

        call setusv_dum(2hin,icolors(i_overlay))

        write(6,*)' sounding overlay is ',i_overlay

!       plot temp and dewpoint sounding
        do iz = 2,n_lvls_snd

            if(istat_sfc .eq. 0      .or.
     1         pres_1d(iz-1) .le. p_sfc_pa          )then  ! above sfc or no 

                call setusv_dum(2hin,icolors(i_overlay))   ! sfc data   
                                                                    
                y1 = logp_1d(iz-1)
                y2 = logp_1d(iz)

                x1 = skewt(k_to_c(temp_vert(iz-1)),y1)
                x2 = skewt(k_to_c(temp_vert(iz)),y2)
                call line(x1,y1,x2,y2) 

                if(istat_td .eq. 1)then
                    x1 = skewt(td_vert(iz-1),y1)
                    x2 = skewt(td_vert(iz),y2)
                    call gsln (itd)
                    call line(x1,y1,x2,y2) 
                    call gsln (1)
                endif

            elseif(pres_1d(iz) .gt. p_sfc_pa)then          ! below sfc
                call setusv_dum(2hin,34)

                y1 = logp_1d(iz-1)
                y2 = logp_1d(iz)

                x1 = skewt(k_to_c(temp_vert(iz-1)),y1)
                x2 = skewt(k_to_c(temp_vert(iz)),y2)
                call line(x1,y1,x2,y2) 

                if(istat_td .eq. 1)then
                    x1 = skewt(td_vert(iz-1),y1)
                    x2 = skewt(td_vert(iz),y2)
                    call gsln (itd)
                    call line(x1,y1,x2,y2) 
                    call gsln (1)
                endif

            else                                           ! straddles the sfc
!               plot line below the sfc              
                call setusv_dum(2hin,34)

                y1 = logp_1d(iz-1)
                y2 = logp_sfc

                x1 = skewt(k_to_c(temp_vert(iz-1)),y1)
                x2 = skewt(k_to_c(t_sfc_k),y2)
                call line(x1,y1,x2,y2) 

                if(istat_td .eq. 1)then
                    x1 = skewt(td_vert(iz-1),y1)
                    x2 = skewt(k_to_c(td_sfc_k),y2)
                    call gsln (itd)
                    call line(x1,y1,x2,y2) 
                    call gsln (1)
                endif

!               plot line above the sfc              
                call setusv_dum(2hin,icolors(i_overlay))

                y1 = logp_sfc
                y2 = logp_1d(iz)

                x1 = skewt(k_to_c(t_sfc_k),y1)
                x2 = skewt(k_to_c(temp_vert(iz)),y2)
                call line(x1,y1,x2,y2) 

                if(istat_td .eq. 1)then
                    x1 = skewt(k_to_c(td_sfc_k),y1)
                    x2 = skewt(td_vert(iz),y2)
                    call gsln (itd)
                    call line(x1,y1,x2,y2) 
                    call gsln (1)
                endif

            endif

        enddo ! iz


!       plot time/source label
        if(box_low .eq. 0.1)then
            call set(0., 1., 0., 1.
     1              ,0., 1., 0., 1., 1)
            call write_label_lplot(100,94,c_label,a9time
     1                            ,plot_parms,namelist_parms
     1                            ,i_overlay,'sound')
        elseif(box_low .eq. 0.2)then
            call set(0., 1., 0., 1.
     1              ,0., 1., 0., 1., 1)
            call write_label_lplot(100,94,c_label,a9time
     1                            ,plot_parms,namelist_parms
     1                            ,i_overlay,'sound')
        endif

!       plot tpw value
        if(pw_sfc .ne. r_missing_data)then
            ix = 795
            iy = 240
            rsize = .010
            write(c_pw,890)pw_sfc*100. ! convert m to cm
 890        format('iwv = ',f5.2)
            call pchiqu (cpux(ix),cpux(iy),c_pw,rsize,0,-1.0)
        endif

!       plot cape value
        if(cape_sfc .ne. r_missing_data)then
            ix = 795
            iy = 220
            rsize = .010
            write(c_cape,891)nint(cape_sfc)
 891        format('cape = ',i4)
            call pchiqu (cpux(ix),cpux(iy),c_cape,rsize,0,-1.0)
        endif

        write(6,*)' sounding has been plotted...'

        go to 40

 900    continue

        call sflush

        return
        end
