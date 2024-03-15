
        subroutine get_vis(i4time,sol_alt,l_use_vis,l_use_vis_add        ! i
     1                    ,l_use_vis_partial,lat,lon,idb,jdb             ! i
     1                    ,i4_sat_window,i4_sat_window_offset            ! i
     1                    ,rlaps_land_frac,topo                          ! i
     1                    ,cvr_snow,tgd_sfc_k                            ! i
     1                    ,offset_vis_i,offset_vis_j                     ! i
     1                    ,di_dh_vis,dj_dh_vis                           ! o
     1                    ,cloud_frac_vis_a,sat_albedo,reflectance       ! o
     1                    ,mode_refl,ihist_alb,static_albedo,sfc_albedo  ! o
     1                    ,vis_snow_max                                  ! o
     1                    ,subpoint_lat_clo,subpoint_lon_clo             ! o 
     1                    ,comment                                       ! o
     1                    ,ni,nj,nk,r_missing_data                       ! i
     1                    ,istat_vis_a                                   ! o
     1                    ,istatus)                                      ! o

!       steve albers 1997 through 2018

        use mem_namelist, only: grid_spacing_m, solalt_thr_vis

        integer ihist_alb(-10:20)
        integer ihist_alb_sfc(-10:20)
        integer ihist_frac_sat(-10:20)
        integer istat_vis_a(ni,nj)     ! cloud mask based on vis (image space)

        real lat(ni,nj), lon(ni,nj)
        real sfc_albedo(ni,nj), sfc_albedo_lwrb(ni,nj)
        real static_albedo(ni,nj)   ! static albedo database
        real sat_data_in(ni,nj)
        real sat_albedo(ni,nj) ! cloud albedo from reflectance/sfc alb
                               ! image space
        real reflectance(ni,nj)     
        real cvr_snow(ni,nj)       ! gridpoint space
        real vis_snow_max(ni,nj)   ! gridpoint space
        real tgd_sfc_k(ni,nj)
        real rlaps_land_frac(ni,nj)
        real topo(ni,nj)
        real subpoint_lat_clo(ni,nj)
        real subpoint_lon_clo(ni,nj)
        real offset_vis_i(ni,nj)    ! sat i minus actual i
        real offset_vis_j(ni,nj)    ! sat j minus actual j
        real di_dh_vis(ni,nj)                      
        real dj_dh_vis(ni,nj)                      
        real i_fill_seams(ni,nj)                      
        real dum2d(ni,nj)                      

!       this stuff is for reading vis data from lvd file
        real sol_alt(ni,nj)
        real sol_alt_sat(ni,nj)
        real cloud_frac_vis_a(ni,nj) ! cloud albedo with clear alb subtracted
        integer mxstn
        parameter (mxstn = 100)       ! max number of "stations" in data file
        character*9 filename
        character*31 ext

        character*3 lvd_ext
        data lvd_ext /'lvd'/

        character var*3,comment*125,units*10
        logical l_use_vis, l_use_vis_add, l_use_vis_partial

!       l_use_vis_add = .false.
        icount_vis_add_potl = 0
        visthr = 0.0  ! for adding visible
        sol_alt_sat(:,:) = sol_alt(:,:)

        if(grid_spacing_m .ne. 3000.)then ! consider namelist 'mode_refl' parm
            mode_refl = 1 ! 0,1 to use reflectance instead of albedo from lvd       
        else
            mode_refl = 1
        endif

!       initialize histograms
        do i = -10,20
            ihist_alb(i) = 0
            ihist_alb_sfc(i) = 0
            ihist_frac_sat(i) = 0
        enddo ! i

!       initialize
        do i = 1,ni
        do j = 1,nj
            sat_albedo(i,j) = r_missing_data
            sat_data_in(i,j) = r_missing_data
            reflectance(i,j) = r_missing_data
            cloud_frac_vis_a(i,j) = r_missing_data
            istat_vis_a(i,j) = 0
        enddo ! j
        enddo ! i

        call get_sfc_albedo(ni,nj,lat,r_missing_data,i4time              ! i
     1                     ,rlaps_land_frac,topo                         ! i
     1                     ,cvr_snow,tgd_sfc_k                           ! i
     1                     ,sfc_albedo,sfc_albedo_lwrb                   ! o
     1                     ,static_albedo,istat_sfc_alb)                 ! o

!       determine whether to use vis / albedo data
        if(.not. l_use_vis)then
            write(6,*)' note: l_use_vis set to not use vis data'
            istatus = 0
            return
        endif

!       read satellite reflectance / albedo data
        ntrys=2
        itry = 1
        istatus = 0
        do while (itry .le. ntrys .and. 
     1            istatus .ne. 1 .and. istatus .ne. -1)
            i4time_try = i4time - ((itry-1)*900)
            write(6,*)' getting the vis data from lvd file - try = '
     1                                                      ,itry      
            call make_fnam_lp(i4time_try,filename,istatus)

            ext = lvd_ext
            if(mode_refl .eq. 0)then
                var = 'alb'
            else
                var = 'svs'
            endif
            ilevel = 0
            call get_laps_2dvar(i4time_try+i4_sat_window_offset
     1                     ,i4_sat_window
     1                     ,i4time_nearest,lat,lon
     1                     ,subpoint_lat_clo,subpoint_lon_clo      ! o 
     1                     ,ext,var,units
     1                     ,comment,ni,nj,sat_data_in,ilevel,istatus)
            write(6,*)' istatus from sat data = ',var,istatus
            write(6,*)' comment from get_vis = ',comment
            itry = itry + 1
        enddo ! itry

        if(istatus .ne. 1 .and. istatus .ne. -1)then
            write(6,*)' no vis / albedo available'
            write(6,*)' return from get_vis'
            istatus = 0
            return
        endif

!       compute parallax info
        write(6,*)' call get_parallax_info in get_vis (vis)'
        call get_parallax_info(ni,nj,i4time                              ! i
     1                        ,lat,lon                                   ! i
     1                        ,subpoint_lat_clo,subpoint_lon_clo         ! i
     1                        ,di_dh_vis,dj_dh_vis,i_fill_seams)         ! o

        ridbg = float(idb) - di_dh_vis(idb,jdb) * topo(idb,jdb)
        rjdbg = float(jdb) - dj_dh_vis(idb,jdb) * topo(idb,jdb)

        idbg = min(max(nint(ridbg),1),ni)
        jdbg = min(max(nint(rjdbg),1),nj)

        write(6,81)idb,jdb,idbg,jdbg,cvr_snow(idbg,jdbg),
     1             tgd_sfc_k(idb,jdb),sfc_albedo(idb,jdb),
     1             sfc_albedo_lwrb(idb,jdb),
     1             static_albedo(idb,jdb),topo(idbg,jdbg)
81      format(' sncv/tgd/sfalb/sfalwb/stal/topo',4i5,5f8.3,f8.1
     1                                           ,' ctr (get_sfc_alb)')

!       possibly 'sfc_albedo_lwrb' should be corrected for parallax        
        if(.false.)then
            call shift_parallax(di_dh_vis,dj_dh_vis,ni,nj                ! i
     1                         ,sfc_albedo_lwrb,topo                     ! i
     1                         ,dum2d)                                   ! o
            sfc_albedo_lwrb(:,:) = dum2d(:,:)
        endif

!       initial test for missing albedo (and partial data coverage)
!       loop in satellite i,j (image space)
        write(6,*)' solalt_thr_vis = ',solalt_thr_vis
        do i = 1,ni
        do j = 1,nj
            ioff = min(max(i+nint(offset_vis_i(i,j)),1),ni)
            joff = min(max(j+nint(offset_vis_j(i,j)),1),nj)

            if(mode_refl .eq. 0)then
                sat_albedo(i,j) = sat_data_in(ioff,joff)
            else
!               convert to reflectance
!               reflectance = (sat_data_in(ioff,joff) / 256.) * 1.2
                if(sat_data_in(ioff,joff) .ne. r_missing_data)then
                    reflectance(i,j) = sat_data_in(ioff,joff)
                endif

!               should this be done elsewhere according to how other
!               routines use sfc_albedo?                
                if(sol_alt(i,j) .gt. solalt_thr_vis .and.
     1             reflectance(i,j) .ne. r_missing_data)then
                    if(i .eq. idb .and. j .eq. jdb)then
                        iverbose = 1
                    else
                        iverbose = 0
                    endif
                    call refl_to_albedo2(reflectance(i,j)              ! i
     1                                  ,sol_alt_sat(i,j)              ! i
     1                                  ,sfc_albedo(i,j),iverbose      ! i
     1                                  ,vis_snow_max(i,j)             ! o
     1                                  ,cloud_albedo)                 ! o (cla)
                else
                    cloud_albedo = r_missing_data
                endif

                sat_albedo(i,j) = cloud_albedo
                comment = 'satellite cloud albedo'
            endif

            if(sat_albedo(i,j) .eq. r_missing_data .and. 
     1                          (.not. l_use_vis_partial)      )then
                write(6,*)
     1              ' no vis / albedo available (missing data found)'          
                write(6,*)' return from get_vis'
                sat_albedo = r_missing_data
                istatus = 0
                return
            endif
        enddo ! j
        enddo ! i

        n_missing_albedo = 0

!       loop in satellite i,j (image space)
        do i = 1,ni
        do j = 1,nj

          rig = float(i) - di_dh_vis(i,j) * topo(i,j)
          rjg = float(j) - dj_dh_vis(i,j) * topo(i,j)

          ig = min(max(nint(rig),1),ni)
          jg = min(max(nint(rjg),1),nj)

!         we will now only use the vis data if the solar alt > 'vis_alt_thr'
!         0 degrees now used to allow 30 min slack in data availability
          if(sol_alt(i,j) .lt. 0.0)then
              if(sat_albedo(i,j) .ne. r_missing_data)then
                  write(6,'(" error -  sat_albedo not missing:"
     1                    ,2f8.3)')sol_alt(i,j),sat_albedo(i,j)
                  stop
              endif
!             sat_albedo(i,j) = r_missing_data
          endif

          if(sat_albedo(i,j) .ne. r_missing_data)then

!           translate the sat_albedo into cloud fraction

!           store histogram information for satellite data
            iscr_alb  = nint(sat_albedo(i,j)*10.)
            iscr_alb  = min(max(iscr_alb,-10),20)
            ihist_alb(iscr_alb) = ihist_alb(iscr_alb) + 1

            if(sfc_albedo_lwrb(i,j) .ne. r_missing_data)then
                iscr_alb_sfc  = nint(sfc_albedo_lwrb(i,j)*10.)
                iscr_alb_sfc  = min(max(iscr_alb_sfc,-10),20)
                ihist_alb_sfc(iscr_alb_sfc) = 
     1          ihist_alb_sfc(iscr_alb_sfc) + 1
            endif

!           'cloud_frac_vis' is what's used in 'insert_vis'
            clear_albedo = sfc_albedo_lwrb(i,j)
            if(mode_refl .eq. 0)then
              cloud_frac_vis = albedo_to_cloudfrac2(clear_albedo
     1                                           ,sat_albedo(i,j))
!    1                                           ,sat_albedo(ig,jg))
            else ! reflectance mode
!             cloud_frac_vis = albedo_to_cloudfrac2(clear_albedo
!    1                                           ,sat_albedo(i,j)*0.6)
!    1                                           ,sat_albedo(ig,jg))

!             should a rayleigh correction be included here?
              if(.false.)then
                  cloud_frac_vis = (sat_albedo(i,j) - clear_albedo)
     1                           / (1.0             - clear_albedo)
              else
                  cloud_frac_vis = sat_albedo(i,j)
              endif
              cloud_frac_vis = min(max(cloud_frac_vis,0.),1.)              
            endif

            iscr_frac_sat = nint(cloud_frac_vis*10.)
            iscr_frac_sat = min(max(iscr_frac_sat,-10),20)
            ihist_frac_sat(iscr_frac_sat) = 
     1      ihist_frac_sat(iscr_frac_sat) + 1

!           make sure satellite cloud fraction is between 0 and 1
            if(cloud_frac_vis .le. 0.0)cloud_frac_vis = 0.0
            if(cloud_frac_vis .ge. 1.0)cloud_frac_vis = 1.0

            cloud_frac_vis_a(i,j) = cloud_frac_vis

!           is there enough of a signal from the vis to say a cloud is present?
!           consider doing this comparison with parallax info (ig,jg)
            if(       cloud_frac_vis_a(i,j) .gt. visthr
     1          .and. sfc_albedo(i,j) .ne. r_missing_data
!    1          .and. sfc_albedo(i,j) .le. 0.3 ! test now done in 'cloud_top'
     1          .and. cvr_snow(ig,jg) .lt. 0.01
     1          .and. l_use_vis_add                         )then
                istat_vis_a(i,j) = 1
                icount_vis_add_potl = icount_vis_add_potl + 1
            endif

            if(i .eq. idb .and. j .eq. jdb)then
              if(mode_refl .eq. 0)then
                write(6,91)clear_albedo,sat_albedo(ig,jg),cloud_frac_vis
91              format(' clralb/satalb/cf_vis ',3f9.3,' ctr')
              elseif(mode_refl .eq. 1)then
                write(6,92)reflectance(ig,jg),sol_alt(ig,jg)
     1                    ,sat_albedo(ig,jg),clear_albedo,cloud_frac_vis
     1                    ,istat_vis_a(i,j)
92              format(' refl/salt/salb/clralb/cf_vis/visadd  ',5f8.3
     1                                                       ,i2,' ctr')
              endif
            endif

          else
            n_missing_albedo =  n_missing_albedo + 1
            cloud_frac_vis_a(i,j) = r_missing_data

          endif

        enddo ! i
        enddo ! j

        write(6,*)
        write(6,*)' n_missing_albedo = ',n_missing_albedo
        write(6,*)

        write(6,*)' number of potential visible clear/add = '
     1            ,ni*nj-n_missing_albedo,icount_vis_add_potl

        if(n_missing_albedo .eq. ni*nj)then ! return with status = 0
            write(6,*)' all albedos were missing - return from get_vis'
            istatus = 0
            return
        endif

        write(6,*)
        write(6,*)'              histograms'
        write(6,*)' i       ',
     1  ' sat albedo  sfc albedo  cld frac sat'
        do i = -5,15
            write(6,11)i,ihist_alb(i),ihist_alb_sfc(i),ihist_frac_sat(i)       
11          format(i4,i12,i12,i12)
        enddo ! i

        write(6,21)minval(cloud_frac_vis_a),maxval(cloud_frac_vis_a)
     1            ,cloud_frac_vis_a(ni/2,nj/2)
21      format(' end of get_vis: cloud_frac_vis_a range/center '
     1        ,2f7.3)

        write(6,*)
        istatus = 1

        return
        end

        subroutine get_sfc_albedo(ni,nj,lat,r_missing_data,i4time    ! i
     1                           ,rlaps_land_frac,topo               ! i
     1                           ,cvr_snow,tgd_sfc_k                 ! i
     1                           ,sfc_albedo,sfc_albedo_lwrb         ! o
     1                           ,static_albedo                      ! o
     1                           ,istat_sfc_alb)                     ! o

!       this returns the surface albedo. this is from the static database
!       to yield a less confident "lower bound". if we are confident that
!       snow/ice is not present, then the 'sfc_albedo' array is also 
!       populated and can be used more at face value.

        character*3 var
        real lat(ni,nj)
        real sfc_albedo(ni,nj)      ! populated with "reliable" values that
                                    ! may include land/sea snow/ice

        real sfc_albedo_lwrb(ni,nj) ! populated with lower bound (i.e. from
                                      ! static database)

        real static_albedo(ni,nj)   ! static albedo database

        real rlaps_land_frac(ni,nj)
        real topo(ni,nj)
        real cvr_snow(ni,nj)
        real tgd_sfc_k(ni,nj)

        write(6,*)' subroutine get_sfc_albedo...'

        istat_sfc_alb = 0
        static_albedo = r_missing_data

        call get_static_field_interp('albedo',i4time,ni,nj
     1                               ,static_albedo,istat_sfc_alb)       
        if(istat_sfc_alb .ne. 1)then ! read sfc albedo from fixed database
            write(6,*)' monthly albedo data n/a, look for fixed data'
            var = 'alb'
            call read_static_grid(ni,nj,var,static_albedo,istat_sfc_alb)
        endif

        write(6,*)' static albedo range: ',minval(static_albedo)
     1                                    ,maxval(static_albedo)

        icount_albedo = 0
        icount_albedo_lwrb = 0

        do i = 1,ni
        do j = 1,nj
            if(istat_sfc_alb .eq. 1 .and. 
     1         static_albedo(i,j) .ne. r_missing_data)then ! static data avalbl
                sfc_albedo_lwrb(i,j) = static_albedo(i,j)

!               'sfc_albedo' is set to missing if we aren't confident enough
!               in its value for use with visible satellite. over water a
!               warm sfc/gnd temperature is used to imply no ice cover and
!               thus a confident value.

                if(rlaps_land_frac(i,j) .le. 0.25)then         ! over water
!                 if((abs(lat(i,j)) .le. 40. .and. topo(i,j) .le.  100.)
                  if((tgd_sfc_k(i,j) .gt. 278.15 .and. 
     1                                             topo(i,j) .le.  100.)
     1                                     .or.
     1               (abs(lat(i,j)) .le. 20. .and. topo(i,j) .le. 3000.)
     1                                          )then          ! it's reliable
                      sfc_albedo(i,j) = static_albedo(i,j)
                  else                               ! sea/lake ice possible?
                      sfc_albedo(i,j) = r_missing_data
                  endif

                else                                           ! over land
                  if(cvr_snow(i,j) .ne. r_missing_data)then
                      sfc_albedo(i,j) = 
     1                max(static_albedo(i,j),cvr_snow(i,j))
                  else
                      sfc_albedo(i,j)      = r_missing_data
                  endif
                endif

            else ! static database not available
                sfc_albedo_lwrb(i,j) = 0.2097063 
                sfc_albedo(i,j)      = r_missing_data

            endif

            if(sfc_albedo_lwrb(i,j) .ne. r_missing_data)then
                icount_albedo_lwrb = icount_albedo_lwrb + 1
            endif

            if(sfc_albedo(i,j) .ne. r_missing_data)then
                icount_albedo = icount_albedo + 1
            endif

        enddo ! j
        enddo ! i

        write(6,*)' number of sfc albedo and lower bound points = '       
     1           ,icount_albedo,icount_albedo_lwrb

        return
        end

        function albedo_to_cloudfrac2(clear_albedo,albedo)

!       cloud_albedo_ref = .4485300
        cloud_albedo_ref = .38       

        arg = albedo

        call stretch2(clear_albedo,cloud_albedo_ref,0.,1.,arg)

        albedo_to_cloudfrac2 = arg

        return
        end

        subroutine shift_parallax(di_dh,dj_dh,ni,nj,field,ht
     1                           ,field_shift)

        real di_dh(ni,nj)                      
        real dj_dh(ni,nj)                      
        real ht(ni,nj)
        real field(ni,nj)
        real field_shift(ni,nj)

        field_shift(:,:) = field(:,:)

        do i = 1,ni
        do j = 1,nj
            is = i + di_dh(i,j) * ht(i,j)
            js = j + dj_dh(i,j) * ht(i,j)
            if(is .ge. 1 .and. is .le. ni .and.
     1         js .le. 1 .and. js .le. nj      )then
                field_shift(is,js) = field(i,j)
            endif
        enddo ! j
        enddo ! i

        return
        end

       
        subroutine refl_to_albedo2(reflectance,solalt,sfc_albedo      ! i
     1                            ,iverbose                           ! i
     1                            ,vis_snow_max                       ! o
     1                            ,cloud_albedo)                      ! o

        include 'trigd.inc'        

        real land_refl,land_refl_dir,land_refl_dif,land_albedo

!       convert reflectance to cloud (+land) albedo
!       note that two solutions may be possible with low sun
!       sfc_albedo can be accounted for?
!       this is presently under development and is being called.
!       'cloud_albedo' assumes 'sfc_albedo' that may include snow
!       'snow_albedo' assumes no snow cover and is an upper bound

        pf_land = 2.0 ! approximate for now for low phase angle
        land_refl_dir = sind(solalt_eff) * 0.9 * sfc_albedo * pf_land ! direct
        land_refl_dif =                    0.1 * sfc_albedo           ! diffuse
        land_refl = land_refl_dir + land_refl_dif

        alb_thn = reflectance           ! modify by phase angle?
        solalt_eff = max(solalt,6.)
!       alb_thk = reflectance / sind(solalt_eff)

        refl_thk = sind(solalt_eff)
        if(reflectance .gt. land_refl .and. refl_thk .gt. land_refl)then
            alb_thk = (reflectance - land_refl)
     1              / (refl_thk    - land_refl)
        else
            alb_thk = 0.
        endif
        alb_thk = max(alb_thk,0.)

        frac_thk = 1.0 ! cloud_albedo
        frac_thn = 1.0 - frac_thk

        cloud_albedo = alb_thn * frac_thn + alb_thk * frac_thk
        air_refl = 0.05 ! approximate rayleigh for 600nm, near nadir

        if(iverbose .eq. 1)then
            write(6,1)reflectance,solalt_eff,sfc_albedo
     1               ,land_refl,air_refl,alb_thk,cloud_albedo
1           format(' refl/salt/salb/lrfl/arfl/athk/cldalb ',7f8.3,
     1             ' ctr (refl_2_alb_c)')
        endif

!       consider snow over dry land
        land_albedo = 0.15
        pf_land = 2.0 ! approximate for now for low phase angle
     
        snow_albedo = 0.7
        pf_snow = 1.0 ! approximate for now for low phase angle

        land_refl_dir = sind(solalt_eff) * 0.9 * land_albedo * pf_land ! direct
        land_refl_dif =                    0.1 * land_albedo           ! diffuse
        land_refl = land_refl_dir + land_refl_dif

        snow_refl = sind(solalt_eff) * snow_albedo * pf_snow
        if(reflectance .ge. land_refl .and. snow_refl .gt. land_refl
     1                                .and. solalt_eff .gt. 10.    )then     
            vis_snow_max = (reflectance - land_refl)
     1                   / (snow_refl   - land_refl)
            vis_snow_max = min(vis_snow_max,1.)            
        elseif(reflectance .lt. land_refl)then
            vis_snow_max = 0.
        else
            vis_snow_max = 1.
        endif        

        if(iverbose .eq. 1)then
            write(6,2)reflectance,solalt_eff,snow_albedo
     1               ,land_refl,snow_refl,vis_snow_max
2           format(' refl/salt/salb/lrfl/srfl/vissnwmx',4x,6f8.3,8x,
     1             ' ctr (refl_2_alb_s)')
        endif

        return
        end
        
        
