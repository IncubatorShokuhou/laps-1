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

        subroutine laps_deriv_sub(i4time,
     1                        nx_l,ny_l,
     1                        nz_l,
     1                        n_pirep,
     1                        maxstns,
     1                        max_snd_grid,max_snd_levels,
     1                        n_prods,
     1                        iprod_number,
     1                        temp_3d,
     1                        heights_3d,
     1                        rh_3d_pct,
     1                        pres_sfc_pa,
     1                        t_sfc_k,
     1                        dbz_max_2d,istatus_lps,    ! o
     1                        twet_snow,l_cloud_only,    ! o
     1                        j_status,istatus)

        use cloud_rad ! cloud radiation and microphysics parameters
        use constants_laps, only: r

        esl(x) = 6.1121*exp(17.67*x/(x+243.5))

        integer       ss_normal,sys_bad_prod,sys_no_data,
     1                  sys_abort_prod

        parameter (ss_normal      =1, ! success
     1             sys_bad_prod   =2, ! inappropriate data, insufficient data
     1             sys_no_data    =3, ! no data
     1             sys_abort_prod =4) ! failed to make a prod

!       1991     steve albers - original version
!       1993 mar steve albers - add lmt product
!       1995 jul steve albers - fix coverage threshold on cloud base
!                               added cloud ceiling
!                               added sfc or 2-d cloud type.
!       1995 nov 1  s. albers - add diagnostic output of cloud ceiling
!                               in precip type comparisons
!       1995 nov 2  s. albers - use sao's to add drizzle in sfc "thresholded"
!                               precip type calculation (lct-ptt)
!       1995 nov 10 s. albers - better handling of cloud tops at top of domain
!                               when "three-dimensionalizing" radar data
!       1995 nov 29 s. albers - improve use of sao's to add snow in ptt field.
!                               cloud ceiling threshold replaced with
!                               thresholds on cloud cover and sfc dewpoint
!                               depression.
!       1995 dec 4  s. albers - use sao's to add rain in sfc "thresholded"
!                               precip type calculation (lct-ptt)
!       1995 dec 13 s. albers - now calls get_radar_ref
!       1996 aug 22 s. albers - now calls read_radar_3dref
!       1996 oct 10 s. albers - max sao cloud cover is now 1.00 + some other
!                               misc cleanup.
!       1997 jul 31 k. dritz  - removed include of lapsparms.for.
!       1997 jul 31 k. dritz  - added call to get_i_perimeter.
!       1997 jul 31 k. dritz  - removed parameter statements for ix_low,
!                               ix_high, iy_low, and iy_high, and instead
!                               compute them dynamically (they are not used
!                               as array bounds, only passed in a call).
!       1997 jul 31 k. dritz  - added nx_l, ny_l as dummy arguments.
!       1997 jul 31 k. dritz  - added call to get_r_missing_data.
!       1997 jul 31 k. dritz  - removed parameter statements for default_base,
!                               default_top, and default_ceiling, and instead
!                               compute them dynamically.
!       1997 jul 31 k. dritz  - removed parameter statement for nhor.  now
!                               initialize c1_name_array dynamically instead
!                               of with a data statement.
!       1997 jul 31 k. dritz  - added nz_l as dummy argument.
!       1997 jul 31 k. dritz  - added call to get_ref_base.
!       1997 aug 01 k. dritz  - added maxstns, ix_low, ix_high, iy_low, and
!                               iy_high as arguments in call to insert_sao.
!       1997 aug 01 k. dritz  - also now pass r_missing_data to barnes_r5.
!       1997 aug 01 k. dritz  - pass r_missing_data to insert_sat.
!       1997 aug 01 k. dritz  - pass ref_base to rfill_evap.

!       prevents clearing out using satellite (hence letting saos dominate)
!       below this altitude (m agl)
        real surface_sao_buffer
        parameter (surface_sao_buffer = 800.)

        real thresh_cvr,default_top,default_base,default_clear_cover
     1                   ,default_ceiling

        parameter       (thresh_cvr = 0.65) ! used to "binaryize" cloud cover

        parameter       (default_clear_cover = .01)

        real thresh_cvr_base,thresh_cvr_top,thresh_cvr_ceiling
        parameter (thresh_cvr_base = 0.1)
        parameter (thresh_cvr_top  = 0.1)
        parameter (thresh_cvr_ceiling = thresh_cvr)

        real thresh_thin_lwc_ice     ! threshold cover for thin cloud lwc/ice
        parameter (thresh_thin_lwc_ice = 0.050) ! .020

        real vis_radar_thresh_cvr,vis_radar_thresh_dbz
        parameter (vis_radar_thresh_cvr = 0.2)  ! 0.2, 0.0
        parameter (vis_radar_thresh_dbz = 10.)  ! 5. , -99.

        real lat(nx_l,ny_l),lon(nx_l,ny_l)
        real topo(nx_l,ny_l)
        real rlaps_land_frac(nx_l,ny_l)
        real solar_alt(nx_l,ny_l)
        real solar_ha(nx_l,ny_l)

        real k_to_c

        logical l_cloud_only
        logical l_packed_output
        logical l_evap_radar

        data l_packed_output /.false./
        data l_evap_radar /.false./

        logical l_fill
        logical l_flag_mvd
        logical l_flag_cloud_type
        logical l_flag_icing_index
        logical l_flag_bogus_w, l_bogus_radar_w, l_deep_vv
        logical l_flag_pcp_type        
        logical l_parse

        logical l_sao_lso
        data l_sao_lso /.true./ ! do things the new way?

        logical l_perimeter
        data l_perimeter /.true./ ! use saos just outside domain?

        include 'laps_cloud.inc'

!       nominal cloud heights. actual ones used are fitted to the terrain.
        real cld_hts_new(kcloud)

        data cld_hts_new/1200.,1300.,1400.,1500.,1600.,1700.,1800.,
     11900.,2000.,2100.,2200.,2400.,2600.,2800.,3000.,3200.,
     23400.,3600.,3800.,4000.,4500.,5000.,5500.,6000.,6500.,
     37000.,7500.,8000.,8500.,9000.,9500.,10000.,11000.,12000.,
     413000.,14000.,15000.,16000.,17000.,18000.,19000.,20000./

        integer iarg

        equivalence (cld_hts,cld_hts_new)

        real clouds_3d(nx_l,ny_l,kcloud)

        integer ihist_alb(-10:20)

        real cloud_top(nx_l,ny_l)
        real cloud_base(nx_l,ny_l)
        real cloud_ceiling(nx_l,ny_l)

        real cldtop_m(nx_l,ny_l)
        real cldtop_m_co2(nx_l,ny_l)
        real cldtop_m_tb8(nx_l,ny_l)

        real cld_pres_1d(kcloud)
        real pres_3d(nx_l,ny_l,nz_l)
        real clouds_3d_pres(nx_l,ny_l,nz_l)

        real cvhz(nx_l,ny_l)
        real cvhz1(nx_l,ny_l),cvew1(nx_l,kcloud)
        real cvr_max(nx_l,ny_l),cvew2(nx_l,kcloud)
        real cvr_sao_max(nx_l,ny_l)
        real cvr_snow_cycle(nx_l,ny_l)
        real cvr_water_temp(nx_l,ny_l)
        real cvr_snow(nx_l,ny_l)
        real band8_mask(nx_l,ny_l)

        character*4 radar_name
        character*31 radarext_3d_cloud
        real radar_ref_3d(nx_l,ny_l,nz_l)
        real heights_3d(nx_l,ny_l,nz_l)

        real vv_to_height_ratio_cu
        real vv_to_height_ratio_sc
        real vv_for_st

        real mvd_3d(nx_l,ny_l,nz_l)
!       real lwc_res_3d(nx_l,ny_l,nz_l)
        real w_3d(nx_l,ny_l,nz_l)

        integer icing_index_3d(nx_l,ny_l,nz_l)

        integer cldpcp_type_3d(nx_l,ny_l,nz_l) ! also contains 3d precip type

!       output array declarations
        real out_array_3d(nx_l,ny_l,nz_l)

        real, allocatable, dimension(:,:,:) :: slwc
        real, allocatable, dimension(:,:,:) :: cice
        real slwc_int(nx_l,ny_l)
        real cice_int(nx_l,ny_l)
        real rain_int(nx_l,ny_l)
        real snow_int(nx_l,ny_l)
        real pice_int(nx_l,ny_l)

        real pcpcnc(nx_l,ny_l,nz_l)
        real raicnc(nx_l,ny_l,nz_l)
        real snocnc(nx_l,ny_l,nz_l)
        real piccnc(nx_l,ny_l,nz_l)

        real cldamt(nx_l,ny_l)
        real cldalb_in(nx_l,ny_l)
        real cldalb_out(nx_l,ny_l)
        real cldod_out(nx_l,ny_l)
        real cldod_out_l(nx_l,ny_l)
        real cldod_out_i(nx_l,ny_l)
        real odint_l(nx_l,ny_l)
        real odint_i(nx_l,ny_l)
        real odint(nx_l,ny_l)
        real visibility(nx_l,ny_l)
        real simvis(nx_l,ny_l)
        real static_albedo(nx_l,ny_l)
        real sfc_albedo(nx_l,ny_l)

!       real snow_2d(nx_l,ny_l)

        character*2 c2_precip_types(0:10)

        character*20 c_z2m

        data c2_precip_types
     1  /'  ','rn','sn','zr','sl','ha','l ','zl','  ','  ','  '/

        character*3 c3_pt_flag
        character*1 c1_r,c1_s
        character*8 c8_project

        integer i2_pcp_type_2d(nx_l,ny_l)
        real r_pcp_type_2d(nx_l,ny_l)

        real dum1_array(nx_l,ny_l)
        real dum2_array(nx_l,ny_l)
        real dum3_array(nx_l,ny_l)
        real dum4_array(nx_l,ny_l)

      ! used for "potential" precip type
        logical l_mask_pcptype(nx_l,ny_l)
        integer ibase_array(nx_l,ny_l)
        integer itop_array(nx_l,ny_l)

        logical l_unresolved(nx_l,ny_l)

        character*1 c1_name_array(nx_l,ny_l)
        character*9 filename

        character*35 time
        character*13 filename13

        integer max_fields
        parameter (max_fields = 10)

        character*255 c_filespec
        character var*3,comment*125,directory*150,ext*31,units*10
        character*3 exts(20)
        character*3 var_a(max_fields)
        character*125 comment_a(max_fields)
        character*10  units_a(max_fields)

!       arrays used to read in satellite data
        real tb8_k(nx_l,ny_l)
        real tb8_cold_k(nx_l,ny_l)
        real albedo(nx_l,ny_l)
        real cloud_frac_vis_a(nx_l,ny_l)
        real cloud_frac_co2_a(nx_l,ny_l)

        real temp_3d(nx_l,ny_l,nz_l)
        real rh_3d_pct(nx_l,ny_l,nz_l)
        real model_3d(nx_l,ny_l,nz_l)

        real t_sfc_k(nx_l,ny_l)
        real t_gnd_k(nx_l,ny_l)
        real sst_k(nx_l,ny_l)
        real td_sfc_k(nx_l,ny_l)
        real pres_sfc_pa(nx_l,ny_l)

!       declarations for lso file stuff
        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real cover_s(maxstns)
        real t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real dd_s(maxstns), ff_s(maxstns)
        real ddg_s(maxstns), ffg_s(maxstns)
        real vis_s(maxstns)
        real pstn_s(maxstns),pmsl_s(maxstns),alt_s(maxstns)
        real store_hgt(maxstns,5),ceil(maxstns),lowcld(maxstns)
        real cover_a(maxstns),rad_s(maxstns)
        integer obstime(maxstns),kloud(maxstns),idp3(maxstns)
        character store_emv(maxstns,5)*1,store_amt(maxstns,5)*4
        character wx_s(maxstns)*8, obstype(maxstns)*8
        character atime*24, infile*256

        integer station_name_len
        parameter (station_name_len = 3)                   
        character c_stations(maxstns)*(station_name_len)    

        character asc_tim_9*9

        real ri_s(maxstns), rj_s(maxstns)

!       product # notification declarations
        integer j_status(20),iprod_number(20)

!       stuff for 2d fields
        real ref_mt_2d(nx_l,ny_l)
        real dbz_low_2d(nx_l,ny_l)
        real dbz_max_2d(nx_l,ny_l)

!       sfc precip and cloud type (lct file)
        real r_pcp_type_thresh_2d(nx_l,ny_l)
        real r_cld_type_2d(nx_l,ny_l)

        character*40 c_vars_req
        character*180 c_values_req

        character*3 lso_ext
        data lso_ext /'lso'/

        istat = init_timer()

        write(6,*)' welcome to the laps_deriv_sub (derived cloud prods)'

        idb = (nx_l/2) + 1
        jdb = (ny_l/2) + 1

        call get_static_field_interp('albedo',i4time,nx_l,ny_l 
     1                              ,static_albedo,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in obtaining static albedo'
            return
        endif
        sfc_albedo = static_albedo

        allocate( slwc(nx_l,ny_l,nz_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate slwc'
            stop
        endif

        allocate( cice(nx_l,ny_l,nz_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate cice'
            stop
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error calling get_r_missing_data'
           stop
        endif

        call get_deriv_parms(mode_evap,l_bogus_radar_w,              ! o
     1                       l_deep_vv,l_cloud_only,                 ! o
     1                       vv_to_height_ratio_cu,                  ! o
     1                       vv_to_height_ratio_sc,                  ! o
     1                       vv_for_st,                              ! o
     1                       c_z2m,                                  ! o
     1                       thresh_cvr_cty_vv,thresh_cvr_lwc,       ! o
     1                       twet_snow,                              ! o
     1                       hydrometeor_scale_cldliq,               ! o
     1                       hydrometeor_scale_cldice,               ! o
     1                       hydrometeor_scale_pcp,                  ! o
     1                       istatus)                                ! o
        if (istatus .ne. 1) then
           write (6,*) 'error calling get_deriv_parms'
           stop
        endif

        if(mode_evap .gt. 0)l_evap_radar = .true.

        default_base     = r_missing_data
        default_top      = r_missing_data
        default_ceiling  = r_missing_data

        do j = 1,ny_l
           do i = 1,nx_l
              c1_name_array(i,j) = ' '
           enddo
        enddo

        call get_ref_base(ref_base,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting ref_base'
           stop
        endif

c determine the source of the radar data
        c_vars_req = 'radarext_3d'

        call get_static_info(c_vars_req,c_values_req,1,istatus)

        if(istatus .eq. 1)then
            radarext_3d_cloud = c_values_req(1:3)
        else
            write(6,*)' error getting radarext_3d'
            goto 9999
        endif

!       radarext_3d_cloud = radarext_3d
!       radarext_3d_cloud = 'v02'

        write(6,*)' radarext_3d_cloud = ',radarext_3d_cloud

c read in laps lat/lon and topo
        call get_laps_domain_95(nx_l,ny_l,lat,lon,topo
     1           ,rlaps_land_frac,grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting laps domain'
            goto 9999
        endif

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' error getting laps_cycle_time'
            goto 9999
        endif

        call get_pres_3d(i4time,nx_l,ny_l,nz_l,pres_3d,istatus)
        if(istatus .ne. 1)goto 9999

        n_lc3 = 1
        n_lps = 2
        n_lcb = 3
        n_lcv = 4
        n_lcp = 5
        n_lwc = 6
        n_lil = 7
        n_lct = 8
        n_lmd = 9
        n_lco = 10
        n_lrp = 11
        n_lty = 12
        n_lmt = 13

        exts(n_lc3) = 'lc3'
        exts(n_lcp) = 'lcp'
        exts(n_lwc) = 'lwc'
        exts(n_lil) = 'lil'
        exts(n_lcb) = 'lcb'
        exts(n_lct) = 'lct'
        exts(n_lcv) = 'lcv'
        exts(n_lmd) = 'lmd'
        exts(n_lco) = 'lco'
        exts(n_lps) = 'lps'
        exts(n_lrp) = 'lrp'
        exts(n_lty) = 'lty'
        exts(n_lmt) = 'lmt'

        do i = 1,13
            j_status(i) = sys_no_data
        enddo

        n_prods = 9
        iprod_start = 5
        iprod_end = 13

        if(.true.)then                    ! read data, then calc derived fields
            i4_elapsed = ishow_timer()

            write(6,*)
            write(6,*)'reading lc3,lt1,lps,lsx,lcv,lcb '
     1               ,'to calculate derived fields'

!           read in data (lc3 - clouds_3d)
            ext = 'lc3'
            call get_clouds_3dgrid(i4time,i4time_lc3
     1                         ,nx_l,ny_l,kcloud,ext
     1                         ,clouds_3d,cld_hts,cld_pres_1d,istatus)
            if(istatus .ne. 1 .or. i4time .ne. i4time_lc3)then
                write(6,*)' error reading 3d clouds'
                goto 999
            endif

!           read in data (lps - radar_ref_3d)
            var = 'ref'
            ext = 'lps'
            call get_laps_3d(i4time,nx_l,ny_l,nz_l
     1       ,ext,var,units,comment,radar_ref_3d,istatus_lps)

            if(istatus_lps .ne. 1)then
                write(6,*)' warning: could not read lps 3d ref, filling'
     1                   ,' array with r_missing_data'
                call constant_3d(radar_ref_3d,r_missing_data
     1                          ,nx_l,ny_l,nz_l)           
                istat_radar_2dref = 0
                istat_radar_3dref = 0
                istat_radar_3dref_orig = 0

            else  ! istatus_lps = 1
                read(comment,510)istat_radar_2dref,istat_radar_3dref
     1                          ,istat_radar_3dref_orig
 510            format(23x,3i3)

!               obtain column max ref
                call get_max_reflect(radar_ref_3d,nx_l,ny_l,nz_l
     1                              ,ref_base,dbz_max_2d)

            endif ! istatus_lps

            write(6,*)' istatus_lps = ',istatus_lps

            var = 'lcv'
            ext = 'lcv'
            call get_laps_2d(i4time,ext,var,units,comment
     1                  ,nx_l,ny_l,cvr_max,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading cvr_max analysis - abort'
                goto999
            endif

            var = 'cla'
            ext = 'lcv'
            call get_laps_2d(i4time,ext,var,units,comment
     1                  ,nx_l,ny_l,cldalb_in,istatus)

            if(istatus .ne. 1 .and. istatus .ne. -1)then
                write(6,*)' error reading cloud albedo analysis - abort'
                goto999
            endif

            var = 'cce'
            ext = 'lcb'
            call get_laps_2d(i4time,ext,var,units,comment
     1                   ,nx_l,ny_l,cloud_ceiling,istatus)
            if(istatus .ne. 1 .and. istatus .ne. -1)then
                write(6,*)' error reading cld_ceiling analysis - abort'       
                goto999
            endif

            var = 'lcb'
            ext = 'lcb'
            call get_laps_2d(i4time,ext,var,units,comment
     1                   ,nx_l,ny_l,cloud_base,istatus)
            if(istatus .ne. 1 .and. istatus .ne. -1)then
                write(6,*)' error reading cld_base analysis - abort'       
                goto999
            endif

!           access sao data from lso files
            ext = 'lso'
            call get_directory(ext,directory,len_dir) ! returns directory
            infile = directory(1:len_dir)//filename13(i4time,ext(1:3))
            call read_surface_old(infile,maxstns,atime,n_meso_g,
     1           n_meso_pos,
     1           n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,
     1           n_obs_pos_g,n_obs_b,
     1           n_obs_pos_b,                   ! we use this as an obs count
     1           c_stations,obstype,        
     1           lat_s,lon_s,elev_s,wx_s,t_s,td_s,dd_s,ff_s,ddg_s,
     1           ffg_s,pstn_s,pmsl_s,alt_s,kloud,ceil,lowcld,cover_a,
     1           rad_s,idp3,store_emv,
     1           store_amt,store_hgt,vis_s,obstime,istat_sfc)

            i4_elapsed = ishow_timer()

            call get_c8_project(c8_project,istatus)
            if(istatus .ne. 1)goto 999

            do i = 1,n_obs_pos_b 

!               does this station not report precip?
                if(obstype(i)(1:4) .eq. 'meso'     .or.
     1             obstype(i)(1:4) .eq. 'cdot'     .or.
     1             obstype(i)(7:8) .eq. '1a'       .or.
     1             wx_s(i)(1:7)    .eq. 'unknown'  .or.       ! already in lso
     1             l_parse(c8_project,'afgwc')     .or.
     1             l_parse(c8_project,'afwa')             )then
                    wx_s(i) = 'unknown'

                endif

            enddo ! i

        endif

!       calculate derived fields
        write(6,*)
        write(6,*)' calculating derived fields'

!       write out cloud grid in pressure coordinates
        write(6,*)' writing out grid in pressure coordinates'

        do k = kcloud,1,-1
          write(6,601)k,clouds_3d(idb,jdb,k)
601       format('  clouds_3d      ctr (lc3)',i4,f9.3)      
        enddo ! k

        call interp_height_pres_fast(nx_l,ny_l,nz_l,kcloud
     1  ,clouds_3d_pres,clouds_3d,heights_3d,cld_hts,istatus)

        do k = nz_l,1,-1
          write(6,621)k,clouds_3d_pres(idb,jdb,k),heights_3d(idb,jdb,k)
621       format('  clouds_3d_pres ctr (lcp)',i4,f9.3,f10.1)      
        enddo ! k

        var = 'lcp'
        ext = 'lcp'
        units = 'fractional'
        comment = 'laps cloud cover'
        call put_laps_3d(i4time,ext,var,units,comment,clouds_3d_pres
     1                                          ,nx_l,ny_l,nz_l)
        j_status(n_lcp) = ss_normal
        i4_elapsed = ishow_timer()


!       calculate and write out lwc, mvd, and icing index
        write(6,*)
        write(6,*)' calling lwc etc. routine (get_cloud_deriv)'
        iflag_slwc = 13 ! new smf lwc
        l_flag_mvd = .true.
        l_flag_cloud_type = .true.
        l_flag_icing_index = .true.
        l_flag_bogus_w = .true.
        l_flag_pcp_type = .true.

!       if(.not. l_flag_cloud_type)then ! read in instead
!       endif

        call get_cloud_deriv(
     1                nx_l,ny_l,nz_l,clouds_3d,cld_hts,
     1                temp_3d,rh_3d_pct,heights_3d,pres_3d,
     1                istat_radar_3dref,radar_ref_3d,grid_spacing_cen_m,       
     1                l_mask_pcptype,ibase_array,itop_array,
     1                iflag_slwc,slwc,cice,
     1                thresh_cvr_cty_vv,thresh_cvr_lwc,
     1                l_flag_cloud_type,cldpcp_type_3d,
     1                l_flag_mvd,mvd_3d,
     1                l_flag_icing_index,icing_index_3d,
     1                vv_to_height_ratio_cu,                               ! i
     1                vv_to_height_ratio_sc,                               ! i
     1                vv_for_st,                                           ! i
     1                l_flag_bogus_w,w_3d,l_bogus_radar_w,                 ! i
     1                l_deep_vv,                                           ! i
     1                twet_snow,                                           ! i
     1                l_flag_pcp_type,                                     ! i
     1                istatus)                                             ! o
        if(istatus .ne. 1)then
            write(6,*)' bad status return from get_cloud_deriv'
            goto 999
        else
            write(6,*)' good status return from get_cloud_deriv'
        endif

        i4_elapsed = ishow_timer()

        if(l_flag_cloud_type)then

!           write 3d cloud type
!           4 most significant bits are precip type, other 4 are cloud type
            do k = 1,nz_l
            do j = 1,ny_l
            do i = 1,nx_l
                iarg = cldpcp_type_3d(i,j,k)
                out_array_3d(i,j,k) = iarg - iarg/16*16         ! 'cty'
            enddo
            enddo
            enddo

            ext = 'cty'
            var = 'cty'
            units = 'none'
            comment = 
     1         'cloud type: (1-10) - st,sc,cu,ns,ac,as,cs,ci,cc,cb'
            call put_laps_3d(i4time,ext,var,units
     1                      ,comment,out_array_3d             
     1                      ,nx_l,ny_l,nz_l)

            i4_elapsed = ishow_timer()

        endif ! l_flag_cloud_type

!       calculate "sfc" cloud type
!       now this is simply the cloud type of the lowest significant
!       (> 0.65 cvr) layer present. if a cb is present, it is used instead.
        write(6,*)' calculate sfc cloud type'
        do i = 1,nx_l
        do j = 1,ny_l
            r_cld_type_2d(i,j) = 0

!           pick lowest "significant" layer
            do k = nz_l,1,-1
                cld_type_3d = out_array_3d(i,j,k)
                if(cld_type_3d .gt. 0)then
                    r_cld_type_2d(i,j) = cld_type_3d
                endif
            enddo ! k

!           pick a cb if present
            do k = nz_l,1,-1
                cld_type_3d = out_array_3d(i,j,k)
                if(cld_type_3d .eq. 10)then
                    r_cld_type_2d(i,j) = cld_type_3d
                endif
            enddo ! k

        enddo ! j
        enddo ! i


!       convert slwc and cice by applying scale factor to parcel method values
        if(hydrometeor_scale_cldliq .ge. 0.)then
            ratio_cldliq =  hydrometeor_scale_cldliq
        else
            ratio_cldliq = -hydrometeor_scale_cldliq / 
     1                  (grid_spacing_cen_m/1000.)
        endif

        if(hydrometeor_scale_cldice .ge. 0.)then
            ratio_cldice =  hydrometeor_scale_cldice
        else
            ratio_cldice = -hydrometeor_scale_cldice / 
     1                  (grid_spacing_cen_m/1000.)
        endif

        cld_ice_ub_gpm3 = .03

        do k = 1,nz_l
        do j = 1,ny_l
        do i = 1,nx_l
            if(slwc(i,j,k) .ne. r_missing_data)then
                slwc(i,j,k) = (slwc(i,j,k) * ratio_cldliq)
            endif
            if(cice(i,j,k) .ne. r_missing_data)then
!               cice(i,j,k) = (cice(i,j,k) * ratio_cldice)
                cice(i,j,k) = min(cice(i,j,k),cld_ice_ub_gpm3)                
        
            endif
        enddo 
        enddo
        enddo

        write(6,*)
        write(6,*)' inserting thin clouds into lwc/ice fields'
        call insert_thin_lwc_ice(clouds_3d,clouds_3d_pres,heights_3d
     1       ,temp_3d,cldalb_in,cld_hts,nx_l,ny_l,nz_l,kcloud,idb,jdb
     1       ,thresh_thin_lwc_ice       
     1       ,pres_3d,slwc,cice,istatus)

        if(istatus .ne. 1)then
            write(6,*)' bad status return from insert_thin_lwc_ice'
            goto 999
        endif

!       note these bounds can introduce artifacts related to metar influence.
!       it would be better to design 'insert_thin_lwc_ice' to put more of        
!       the hydrometeors at lower levels upstream so limits set here would
!       have less effect.
        if(.false.)then

!         apply temperature dependent upper bound to cloud ice
          do k = 1,nz_l
          do j = 1,ny_l
          do i = 1,nx_l
!           temperature dependent threshold, we can consider basing this
!           from saturation vapor density (svd * .025) = (es/rt) * .10
!           where .10 means 10% of water saturation
            es_pa = esl(temp_3d(i,j,k) - 273.15) * 100.
            cld_ice_ub_gpm3 = 1e3 * .10 * es_pa / (r * temp_3d(i,j,k))
!           if(temp_3d(i,j,k) .gt. 243.15)then
!               cld_ice_ub_gpm3 = 0.1
!           else
!               cld_ice_ub_gpm3 = 0.03
!           endif
            if(cice(i,j,k) .ne. r_missing_data .and.
     1         cice(i,j,k) .gt. cld_ice_ub_gpm3      )then
               cice(i,j,k) = cld_ice_ub_gpm3
            endif         
            if(i .eq. idb .and. j .eq. jdb)then
                write(6,701)k,temp_3d(i,j,k),es_pa,cld_ice_ub_gpm3
     1                       ,cice(i,j,k),slwc(i,j,k)
701             format('k/t/e/cld_ice_uprb_gpm3/cice/slwc',i4,f8.2,f8.1
     1                       ,3f9.5)
            endif
          enddo
          enddo
          enddo

          i4_elapsed = ishow_timer()

!         apply upper bound of 0.35 g/m**3
          slwc = min(slwc,0.35)
          cice = min(cice,0.35)

        endif ! apply bounds

!       calculate and write integrated lwc
!       note slwc/cice is here in g/m**3
        write(6,*)
        write(6,*)' calculating integrated lwc and cice'

!       this routine can also return od using 'reff_clwc_f' and 'reff_cice_f'
        if(.false.)then
          call integrate_slwc(slwc,heights_3d,nx_l,ny_l,nz_l,slwc_int)
          call integrate_slwc(cice,heights_3d,nx_l,ny_l,nz_l,cice_int)

        else ! integrate slwc and check od
          ihtype = 1
          call integrate_slwc_od(slwc,heights_3d,temp_3d      ! i
     1                               ,nx_l,ny_l,nz_l,ihtype   ! i
     1                               ,slwc_int,odint_l)       ! o
          ihtype = 2
          call integrate_slwc_od(cice,heights_3d,temp_3d      ! i
     1                               ,nx_l,ny_l,nz_l,ihtype   ! i
     1                               ,cice_int,odint_i)       ! o
          odint = odint_l + odint_i
          write(6,801)odint(idb,jdb),cldalb_in(idb,jdb)
801       format(' ctr odint/cldalb_in:',2f9.4)

        endif

!       1e6 factor converts from metric tons/m**2 to g/m**2 (lwp)
        write(6,*)' integrated slwc range (g/m**2) is '
     1                       ,minval(slwc_int)*1e6,maxval(slwc_int)*1e6

        write(6,*)' integrated cice range (g/m**2) is '
     1                       ,minval(cice_int)*1e6,maxval(cice_int)*1e6

!       calculate cloud optical depth and cloud albedo
!       constants are mass extinction efficiency (mee) and mass backscatter
!       efficiency (mbe).
!       note the 1e3 term converts 'slwc_int' or 'clwc_int' units from mg/m**2
!       to kg/m**2 (si units)

        const_lwp = (1.5 * 1e3) / (rholiq     * reff_clwc)
        const_lwp_bks = const_lwp * bksct_eff_clwc

        const_iwp = (1.5 * 1e3) / (rholiq     * reff_cice)
        const_iwp_bks = const_iwp * bksct_eff_cice

        const_rwp = (1.5 * 1e3) / (rholiq     * reff_rain)
        const_rwp_bks = const_rwp * bksct_eff_rain
  
        const_swp = (1.5 * 1e3) / (rhosnow    * reff_snow)
        const_swp_bks = const_swp * bksct_eff_snow

        const_gwp = (1.5 * 1e3) / (rhograupel * reff_graupel)
        const_gwp_bks = const_gwp * bksct_eff_graupel

!       note that tau = od = mee times lwp

!       cloud amount is opacity of cloud liquid and cloud ice hydrometeors
        do j = 1,ny_l
        do i = 1,nx_l
            cldamt(i,j) = 1. - (exp( -(const_lwp * slwc_int(i,j) 
     1                               + const_iwp * cice_int(i,j))) ) 
        enddo ! i
        enddo ! j

        i4_elapsed = ishow_timer()

!       derived radar/precip stuff
        if(istat_radar_3dref .eq. 1)then ! lmt

            if(l_evap_radar)then 

                write(6,*)' calling rfill_evap: mode_evap = ',mode_evap       

!               use lps reflectivity field

!               use cloud base field

                i4_elapsed = ishow_timer()

!               apply evaporation subroutine
                call rfill_evap(radar_ref_3d,nx_l,ny_l,nz_l
     1          ,cloud_base,lat,lon,topo,mode_evap
     1          ,temp_3d,rh_3d_pct,cldpcp_type_3d,heights_3d,istatus
     1          ,ref_base)

                i4_elapsed = ishow_timer()

            endif ! l_evap_radar, etc.

          ! do max tops

            i4_elapsed = ishow_timer()

            write(6,*)' getting max tops'

            call get_maxtops(radar_ref_3d,heights_3d,nx_l,ny_l,nz_l
     1                      ,ref_mt_2d)

!           get laps reflectivities at the surface (or immediately above it)
            write(6,*)' getting low level reflectivity'
            call get_low_ref(radar_ref_3d,pres_sfc_pa,nx_l,ny_l,nz_l
     1                      ,dbz_low_2d)

            istat_pty = 0

!           do sfc precip type
            var = 'td'
            ext = 'lsx'
            call get_laps_2d(i4time,ext,var,units,comment
     1                  ,nx_l,ny_l,td_sfc_k,istatus)
            if(istatus .ne. 1)then
                write(6,*)
     1          ' error reading sfc td - sfc precip type not computed'       
                goto700
            endif

            i4_elapsed = ishow_timer()

!           note that pres_sfc_pa was already read in above
            call get_sfc_preciptype(pres_sfc_pa,t_sfc_k,td_sfc_k
     1                             ,cldpcp_type_3d,twet_snow
     1                             ,dbz_low_2d,i2_pcp_type_2d
     1                             ,nx_l,ny_l,nz_l)

!           compute thresholded precip type
            do i = 1,nx_l
            do j = 1,ny_l
                iarg = i2_pcp_type_2d(i,j)
                r_pcp_type_2d(i,j) = iarg/16

!               apply a threshold to the sfc precip type that depends
!               on both the low level reflectivity and the "potential"
!               sfc precip type.

                if(r_pcp_type_2d(i,j) .eq. 1.0        ! rain
     1        .or. r_pcp_type_2d(i,j) .eq. 3.0        ! zr
     1        .or. r_pcp_type_2d(i,j) .eq. 4.0        ! ip
     1        .or. r_pcp_type_2d(i,j) .eq. 5.0        ! hail
     1                                                      )then

                    r_pcp_type_thresh_2d(i,j) = r_pcp_type_2d(i,j)

                    if(dbz_low_2d(i,j) .lt. 13.0)then
                        r_pcp_type_thresh_2d(i,j) = 0. ! no precip
                    endif
                else ! apply dry threshold to snow

                    call nowrad_virga_correction(r_pcp_type_2d(i,j), ! i
     1                                    r_pcp_type_thresh_2d(i,j), ! o
     1                                    t_sfc_k(i,j),              ! i
     1                                    td_sfc_k(i,j),             ! i
     1                                    istat_radar_3dref_orig)    ! i

                    if(.true.)then ! testcode
                        if(r_pcp_type_2d(i,j) .eq. 2.0 .and.
     1                      r_pcp_type_thresh_2d(i,j) .eq. 0.0)then
                            dbz_low_2d(i,j) = ref_base
                        endif
                    endif

                endif

            enddo
            enddo

            i4_elapsed = ishow_timer()

!           add sao drizzle to the thresh sfc precip type (if applicable out
!           of radar range)
            call sao_drizzle_correction(r_pcp_type_thresh_2d,nx_l,ny_l
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k
     1              ,cloud_ceiling,r_missing_data)

            i4_elapsed = ishow_timer()

!           add sao rain to the thresh sfc precip type (if applicable out
!           of radar range)
            call sao_rain_correction(r_pcp_type_thresh_2d,nx_l,ny_l
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k,twet_snow
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)

            i4_elapsed = ishow_timer()

!           add sao snow to the thresh sfc precip type (if applicable out
!           of radar range)
            call sao_snow_correction(r_pcp_type_thresh_2d,nx_l,ny_l
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k,twet_snow
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)

            i4_elapsed = ishow_timer()

!           increase radar reflectivity threshold for precip in those areas
!           where the radar has echo, but the sao says no precip.
            call sao_precip_correction(r_pcp_type_thresh_2d,nx_l,ny_l
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)

            i4_elapsed = ishow_timer()

!           write lmt/llr
!           note that these arrays start off with 1 as the first index
            var_a(1) = 'lmt'
            var_a(2) = 'llr'
            ext = 'lmt'
            units_a(1) = 'm'
            units_a(2) = 'dbz'
            comment_a(1) = 'laps maximum tops'
            comment_a(2) = 'laps low level reflectivity'

            call move(ref_mt_2d, out_array_3d(1,1,1),nx_l,ny_l)
            call move(dbz_low_2d,out_array_3d(1,1,2),nx_l,ny_l)

            call put_laps_multi_2d(i4time,ext,var_a,units_a,
     1          comment_a,out_array_3d,nx_l,ny_l,2,istatus)

            if(istatus .eq. 1)then
                j_status(n_lmt) = ss_normal
                write(6,*)' success in writing out lmt'
            else
                write(6,*)' error detected writing out lmt'
            endif


!           compare precip type to the obs
            iarg = 0
            if(istat_radar_3dref .eq. 1 .and. istat_sfc .eq. 1)then
              call make_fnam_lp(i4time,asc_tim_9,istatus)
              write(6,*)' ',asc_tim_9,' __'
              write(6,*)' comparing precip type to the obs:',n_obs_pos_b
              write(6,*)' sta pty ptt dbz   t(anl)  td(anl)'//
     1                  '  t(ob)  td(ob)  pty(ob)'//
     1                  '        tw(anl) tw(ob) elev(ob) celg'
              do i = 1,n_obs_pos_b
                call latlon_to_rlapsgrid(lat_s(i),lon_s(i),lat,lon
     1                     ,nx_l,ny_l,ri,rj,istatus)

                i_i = nint(ri)
                i_j = nint(rj)

                if(i_i .ge. 1 .and. i_i .le. nx_l .and.
     1             i_j .ge. 1 .and. i_j .le. ny_l            )then
                    i_pty = i2_pcp_type_2d(i_i,i_j) / 16
                    i_ptt = r_pcp_type_thresh_2d(i_i,i_j)

!                   does laps analyze or does the sao report wx here
                    if(  (i_pty .gt. 0 .or. wx_s(i) .ne. '        ')

!                        and is this an sao station that reports precip?
!    1                  .and.  obstype(i)(1:4) .ne. 'meso'
!    1                  .and.  obstype(i)(1:4) .ne. 'cdot'
!    1                  .and.  obstype(i)(7:8) .ne. '1a'
     1                  .and.  wx_s(i)(1:7)    .ne. 'unknown'
     1                                                         )then
                        c3_pt_flag = ' __'
                    else
                        c3_pt_flag = '   '
                    endif

                    t_sfc_c  = k_to_c(t_sfc_k(i_i,i_j))
                    td_sfc_c = k_to_c(td_sfc_k(i_i,i_j))
                    p_sfc_mb = pres_sfc_pa(i_i,i_j) / 100.

                    tw_sfc_c = tw(t_sfc_c,td_sfc_c,p_sfc_mb)

                    t_s_c  = f_to_c(t_s(i))
                    td_s_c = f_to_c(td_s(i))

                    if(t_s_c .gt. -50. .and. td_s_c .gt. -50.)then
                        tw_s_c = tw(t_s_c,td_s_c,p_sfc_mb)
                    else
                        tw_s_c = -99.
                    endif

                    if(cloud_ceiling(i_i,i_j) .ne. r_missing_data)
     1                                                     then
                        iceil = nint(cloud_ceiling(i_i,i_j))
                    else
                        iceil = 99999
                    endif

!                   snow
                    call parse_wx_pcp(wx_s(i),'s',ipresent,istatus)       
                    if(ipresent .eq. 1)then
                        c1_s = 's'
                    else
                        c1_s = ' '
                    endif

!                   rain
                    call parse_wx_pcp(wx_s(i),'r',ipresent,istatus)       
                    if(ipresent .eq. 1)then
                        c1_r = 'r'
                    else
                        c1_r = ' '
                    endif

                    call filter_string(wx_s(i))

                    write(6,1101,err=1102)c_stations(i)(1:3)
     1                          ,c2_precip_types(i_pty)
     1                          ,c2_precip_types(i_ptt)
     1                          ,int(dbz_low_2d(i_i,i_j))
     1                          ,t_sfc_c
     1                          ,td_sfc_c
     1                          ,t_s_c
     1                          ,td_s_c
     1                          ,wx_s(i)
     1                          ,c3_pt_flag
     1                          ,tw_sfc_c
     1                          ,tw_s_c
     1                          ,elev_s(i)
     1                          ,iceil
     1                          ,cvr_max(i_i,i_j)
     1                          ,c1_r,c1_s,obstype(i)(7:8)
1101                format(1x,a3,2x,a2,2x,a2,i4,4f8.1,3x,a8,2x,a3
     1                                ,3f8.1,i7,f5.2,1x,2a1,1x,a2
     1                                ,' ptvrf')
1102            endif ! ob is in domain
              enddo ! i
            endif

            istat_pty = 1

        else ! set sfc preciptype to missing
            write(6,*)
     1      ' no sfc preciptype calculated due to lack of radar data'       
            do i = 1,nx_l
            do j = 1,ny_l
                r_pcp_type_thresh_2d(i,j) = r_missing_data
                r_pcp_type_2d(i,j) = r_missing_data
            enddo ! j
            enddo ! i

        endif ! istat_radar_3dref

 700    continue

        if(istat_radar_3dref .eq. 1)then
            if(l_flag_bogus_w .and. l_bogus_radar_w) then    ! adan add
!             re-calculate cloud bogus omega within radar echo area
!             add by adan
              call get_radar_deriv(nx_l,ny_l,nz_l,grid_spacing_cen_m,       
     1                           r_missing_data,
     1                           radar_ref_3d,clouds_3d,cld_hts,
     1                           temp_3d,heights_3d,pres_3d,
     1                           ibase_array,itop_array,thresh_cvr,
     1                           vv_to_height_ratio_cu,                 ! i
     1                           cldpcp_type_3d,w_3d,istat_radar_deriv)       
              if(istat_radar_deriv .ne. 1)then
                write(6,*)' bad status return from get_radar_deriv'
              endif
            endif                           ! l_flag_bogus_w (adan add)

            i4_elapsed = ishow_timer()

            write(6,*)' computing precip concentration'

!           calculate 3d precip concentration in kg/m**3
            call cpt_pcp_cnc(radar_ref_3d,temp_3d           ! input
     1                                  ,rh_3d_pct          ! input
     1                                  ,cldpcp_type_3d     ! input
     1                                  ,nx_l,ny_l,nz_l     ! input
     1                                  ,c_z2m              ! input
     1                                  ,pres_3d            ! input
     1                                  ,pcpcnc             ! output
     1                                  ,raicnc             ! output
     1                                  ,snocnc             ! output
     1                                  ,piccnc)            ! output

!           calculate integrated rainwater
            write(6,*)
            write(6,*)' calculating integrated rainwater'
            call integrate_slwc(raicnc,heights_3d,nx_l,ny_l,nz_l
     1                                                     ,rain_int)
            write(6,*)' integrated rain range is ',minval(rain_int)
     1                                            ,maxval(rain_int)

            i4_elapsed = ishow_timer()

!           calculate integrated snow
            write(6,*)
            write(6,*)' calculating integrated snow'
            call integrate_slwc(snocnc,heights_3d,nx_l,ny_l,nz_l
     1                                                     ,snow_int)
            write(6,*)' integrated snow range is ',minval(snow_int)
     1                                            ,maxval(snow_int)

            i4_elapsed = ishow_timer()

!           calculate integrated precipitating ice
            write(6,*)
            write(6,*)' calculating integrated precip. ice'
            call integrate_slwc(piccnc,heights_3d,nx_l,ny_l,nz_l
     1                                                     ,pice_int)
            write(6,*)' integrated pice range is ',minval(pice_int)
     1                                            ,maxval(pice_int)

            i4_elapsed = ishow_timer()

            nf = 6

        else
            nf = 2

        endif

!       rain, snow, and graupel are added for the cldalb & simvis computation
        do j = 1,ny_l
        do i = 1,nx_l
!           albedo = (b * tau) / (1. + b * tau)
            btau = (const_lwp_bks * slwc_int(i,j) + 
     1              const_iwp_bks * cice_int(i,j) + 
     1              const_rwp_bks * rain_int(i,j) + 
     1              const_swp_bks * snow_int(i,j) + 
     1              const_gwp_bks * pice_int(i,j)   ) 
            cldalb_out(i,j) = btau / (1. + btau) 

            cldod_out(i,j) = const_lwp * slwc_int(i,j)
     1                     + const_iwp * cice_int(i,j)  
     1                     + const_rwp * rain_int(i,j)  
     1                     + const_swp * snow_int(i,j)  
     1                     + const_gwp * pice_int(i,j)  

            cldod_out_l(i,j) = const_lwp * slwc_int(i,j)
     1                       + const_rwp * rain_int(i,j)  

            cldod_out_i(i,j) = const_iwp * cice_int(i,j)  
     1                       + const_swp * snow_int(i,j)  
     1                       + const_gwp * pice_int(i,j)  

            simvis(i,j) = cldalb_out(i,j) + (1.-cldalb_out(i,j))**2 
     1          * (sfc_albedo(i,j)/(1.-cldalb_out(i,j)*sfc_albedo(i,j)))

            if(i .eq. idb .and. j .eq. jdb)then
                cvrmax = maxval(clouds_3d(i,j,:))
                write(6,1201)cldod_out(i,j),cldod_out_l(i,j)
     1                      ,cldod_out_i(i,j)
 1201           format(' ctr cloud od/l/i',3f9.3)
                write(6,1202)cvrmax,cldod_out(i,j),cldalb_in(i,j)
     1                      ,cldalb_out(i,j),sfc_albedo(i,j),simvis(i,j)
 1202           format(' ctr cloud cvr/od/albi-o/sfcalb/smv:',6f9.3)
                write(6,1203)slwc_int(i,j),cice_int(i,j)
 1203           format(' ctr slwc_int/cice_int (mg/m**2)',2f10.6)
            endif

        enddo ! i
        enddo ! j

!       convert slwc and cice from g/m**3 to kg/m**3
        do k = 1,nz_l
        do j = 1,ny_l
        do i = 1,nx_l
            if(slwc(i,j,k) .ne. r_missing_data)then
                slwc(i,j,k) = (slwc(i,j,k) * 1.0) / 1e3
            endif
            if(cice(i,j,k) .ne. r_missing_data)then
                cice(i,j,k) = (cice(i,j,k) * 1.0) / 1e3
            endif
        enddo 
        enddo
        enddo

!       calculate low level hydrometeor fields and visibility
!       these constants can be derived from microphysical constants
        clwc2alpha = 1.5 / (rholiq  * reff_clwc)
        cice2alpha = 1.5 / (rholiq  * reff_cice)
        rain2alpha = 1.5 / (rholiq  * reff_rain)
        snow2alpha = 1.5 / (rhosnow * reff_snow)
        pice2alpha = 1.5 / (rhograupel * reff_graupel)

        a1 = clwc2alpha
        b1 = cice2alpha
        c1 = rain2alpha
        d1 = snow2alpha
        e1 = pice2alpha
      
        iwrite = 0

        do j = 1,ny_l
        do i = 1,nx_l
            k_topo = int(zcoord_of_pressure(pres_sfc_pa(i,j)))
            slwc_low = slwc(i,j,k_topo+1)
            cice_low = cice(i,j,k_topo+1)
            rain_low = raicnc(i,j,k_topo+1)
            snow_low = snocnc(i,j,k_topo+1)
            pice_low = piccnc(i,j,k_topo+1)
            alpha = a1 * slwc_low + b1 * cice_low + c1 * rain_low 
     1            + d1 * snow_low + e1 * pice_low
            visibility(i,j) = 2.8 / max(alpha,.00004) 
            if(alpha .gt. 0.0002 .and. iwrite .le. 5 .or.
     1             i .eq. idb .and. j .eq. jdb      )then
                write(6,*)'visibility terms at',i,j,k_topo
                write(6,*)slwc_low,cice_low,rain_low,snow_low,pice_low
     1                   ,dbz_low_2d(i,j),alpha,visibility(i,j)
                iwrite = iwrite + 1
            endif
        enddo ! i
        enddo ! j
      
!       write lil file
!       note that these arrays start off with 1 as the first index
        ext = 'lil'
        var_a(1) = 'lil'
        var_a(2) = 'lic'
        var_a(3) = 'cod'
        var_a(4) = 'cla'
        var_a(5) = 'vis'
        var_a(6) = 'smv'
        units_a(1) = 'm'
        units_a(2) = 'm'
        units_a(3) = ' '
        units_a(4) = ' '
        units_a(5) = 'm'
        units_a(6) = ' '
        comment_a(1) = 'integrated cloud liquid'
        comment_a(2) = 'integrated cloud ice'
        comment_a(3) = 'cloud optical depth'
        comment_a(4) = 'cloud albedo'
        comment_a(5) = 'visibility'
        comment_a(6) = 'visible albedo'

        call move(slwc_int,  out_array_3d(1,1,1),nx_l,ny_l)
        call move(cice_int,  out_array_3d(1,1,2),nx_l,ny_l)
        call move(cldod_out, out_array_3d(1,1,3),nx_l,ny_l)
        call move(cldalb_out,out_array_3d(1,1,4),nx_l,ny_l)
        call move(visibility,out_array_3d(1,1,5),nx_l,ny_l)
        call move(simvis    ,out_array_3d(1,1,6),nx_l,ny_l)

        call put_laps_multi_2d(i4time,ext,var_a,units_a,
     1      comment_a,out_array_3d,nx_l,ny_l,6,istatus)

        if(istatus .eq. 1)then
            j_status(n_lil) = ss_normal
            write(6,*)' success in writing out lil'
        else
            write(6,*)' error detected writing out lil'
        endif

!       apply hydrometeor scale to precip
        if(hydrometeor_scale_pcp .ge. 0.)then
            ratio_pcp =  hydrometeor_scale_pcp
        else
            ratio_pcp = -hydrometeor_scale_pcp / 
     1                  (grid_spacing_cen_m/1000.)
        endif

        pcpcnc = pcpcnc * ratio_pcp
        raicnc = raicnc * ratio_pcp
        snocnc = snocnc * ratio_pcp
        piccnc = piccnc * ratio_pcp

!       write out cloud liquid water, cloud ice and precip content fields
        var_a(1) = 'lwc'
        var_a(2) = 'ice'
        var_a(3) = 'pcn'
        var_a(4) = 'rai'
        var_a(5) = 'sno'
        var_a(6) = 'pic'

        ext = 'lwc'

        units_a(1) = 'kg/m**3'
        units_a(2) = 'kg/m**3'
        units_a(3) = 'kg/m**3'
        units_a(4) = 'kg/m**3'
        units_a(5) = 'kg/m**3'
        units_a(6) = 'kg/m**3'

        comment_a(1) = 'cloud liquid water content - laps smith feddes'       
        comment_a(2) = 'cloud ice content - laps smith feddes'
        comment_a(3) = 'precipitate content'
        comment_a(4) = 'rain content'
        comment_a(5) = 'snow content'
        comment_a(6) = 'precipitating ice content'

        slwc_max = maxval(slwc)
        cice_max = maxval(cice)
        pcpcnc_max = maxval(pcpcnc)
        raicnc_max = maxval(raicnc)
        snocnc_max = maxval(snocnc)
        piccnc_max = maxval(piccnc)

        write(6,*)' max range (g/m**3) slwc = ',slwc_max*1e3
        write(6,*)' max range (g/m**3) cice = ',cice_max*1e3
        write(6,*)' max range (g/m**3) pcpcnc = ',pcpcnc_max*1e3
        write(6,*)' max range (g/m**3) raicnc = ',raicnc_max*1e3
        write(6,*)' max range (g/m**3) snocnc = ',snocnc_max*1e3
        write(6,*)' max range (g/m**3) piccnc = ',piccnc_max*1e3

        if(slwc_max   .le. 1e6 .and. cice_max   .le. 1e6 .and.
     1     pcpcnc_max .le. 1e6 .and. raicnc_max .le. 1e6 .and.
     1     snocnc_max .le. 1e6 .and. piccnc_max .le. 1e6       )then

            call put_laps_3d_multi(i4time,ext,var_a,units_a,comment_a
     1                        ,slwc,cice
     1                        ,pcpcnc,raicnc
     1                        ,snocnc,piccnc
     1                        ,nx_l,ny_l,nz_l
     1                        ,nx_l,ny_l,nz_l
     1                        ,nx_l,ny_l,nz_l
     1                        ,nx_l,ny_l,nz_l
     1                        ,nx_l,ny_l,nz_l
     1                        ,nx_l,ny_l,nz_l
     1                        ,nf,istatus)
            if(istatus .eq. 1)j_status(n_lwc) = ss_normal
        else
            write(6,*)' error: large hydrometeor values lwc not written'
        endif

!       write out mean volume diameter field (potentially compressible)
        ext = 'lmd'
        var = 'lmd'
        units = 'm'
        comment = 'mean volume diameter of cloud droplets'
        call put_laps_3d(i4time,ext,var,units,comment,mvd_3d
     1                                           ,nx_l,ny_l,nz_l)
        j_status(n_lmd) = ss_normal

!       write out icing index field  (potentially compressible)
        ext = 'lrp'
        do k = 1,nz_l
        do j = 1,ny_l
        do i = 1,nx_l
            out_array_3d(i,j,k) = icing_index_3d(i,j,k)
        enddo
        enddo
        enddo

        var = 'lrp'
        units = 'none'
        comment = 'icing severity index '//
     1  '1-ltcnt,2-mdcnt,3-hvcnt,4-ltint,5-mdint,6-hvint'
        call put_laps_3d(i4time,ext,var,units,comment,out_array_3d
     1                                  ,nx_l,ny_l,nz_l)

        j_status(n_lrp) = ss_normal

!       write 3d precip type
!       4 most significant bits are precip type, other 4 are cloud type
        do k = 1,nz_l
        do j = 1,ny_l
        do i = 1,nx_l
            iarg = cldpcp_type_3d(i,j,k)
            out_array_3d(i,j,k) = iarg/16                   ! 'pty'
        enddo
        enddo
        enddo

        ext = 'pty'
        var = 'pty'
        units = 'none'
        comment = 'precip type: 0-none,1-rain,2-snow,3-zr,4-ip,5-hail'
        call put_laps_3d(i4time,ext,var,units
     1                  ,comment,out_array_3d(1,1,1)
     1                  ,nx_l,ny_l,nz_l)

        i4_elapsed = ishow_timer()

!       write sfc precip and cloud type
        var_a(1) = 'pty'
        var_a(2) = 'ptt'
        var_a(3) = 'sct'
        ext = 'lct'
        units_a(1) = 'undim'
        units_a(2) = 'undim'
        units_a(3) = 'undim'
        comment_a(1) = 'laps precip type (unthresholded): '//
     1                 '0:np 1:rn 2:sn 3:zr 4:sl 5:ha 6:l  7:zl'
        comment_a(2) = 'laps precip type (refl threshold): '//
     1                 '0:np 1:rn 2:sn 3:zr 4:sl 5:ha 6:l  7:zl'
        comment_a(3) = 'laps cloud type '//
     1                 '0:clr 1:st 2:sc 3:cu 4:ns 5:ac '//
     1                 '6:as 7:cs 8:ci 9:cc 10: cb'

        call move(r_pcp_type_2d,       out_array_3d(1,1,1),nx_l,ny_l)
        call move(r_pcp_type_thresh_2d,out_array_3d(1,1,2),nx_l,ny_l)
        call move(r_cld_type_2d,       out_array_3d(1,1,3),nx_l,ny_l)        

        call put_laps_multi_2d(i4time,ext,var_a,units_a,
     1                 comment_a,out_array_3d,nx_l,ny_l,3,istatus)

        if(istatus .eq. 1)j_status(n_lct) = ss_normal

!       write out cloud derived omega field
        var = 'com'
        ext = 'lco'
        units = 'pa/s'
        comment = 'laps cloud derived omega'
        call put_laps_3d(i4time,ext,var,units,comment,w_3d
     1                  ,nx_l,ny_l,nz_l)
        j_status(n_lco) = ss_normal

        i4_elapsed = ishow_timer()

        istatus = 1

!       read sounding metadata to get integrated cloud liquid obs
!       lun = 89
!       ext = 'snd'
!       call read_snd_metadata(lun,i4time,ext                         ! i
!    1                        ,max_snd_grid,max_snd_levels            ! i
!    1                        ,lat,lon,imax,jmax                      ! i
!    1                        ,n_profiles                             ! o
!    1                        ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr   ! o
!    1                        ,c5_name,i4time_ob_pr,obstype           ! o
!    1                        ,cloud_base_temp,cloud_liquid           ! o
!    1                        ,istatus)                               ! o

!       write(6,*)' back from read_snd_metadata, # soundings = '
!    1            ,n_profiles

999     continue

        write(6,*)' notifications'
        do i = iprod_start,iprod_end
            write(6,*)' ',exts(i),' ',j_status(i),' ',i
        enddo ! i

9999    deallocate( slwc )
        deallocate( cice )

        write(6,*)' end of laps_deriv_sub'

        return
        end


