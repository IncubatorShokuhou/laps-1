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

        subroutine lapswind_plot(c_display,i4time_ref,lun,nx_l,ny_l,
     1                           nz_l, max_radars, l_radars,
     1                           r_missing_data,
     1                           laps_cycle_time,zoom,density,
     1                           dyn_low,dyn_high,dx,dy,
     1                           plot_parms,namelist_parms,ifield_found)

!       1995        steve albers         original version
!       1995 dec 8  steve albers         automated pressure range
!       97-aug-14     ken dritz     added nx_l, ny_l, nz_l as dummy args
!       97-aug-14     ken dritz     added max_radars as dummy arg
!       97-aug-14     ken dritz     added r_missing_data as dummy arg
!       97-aug-14     ken dritz     added laps_cycle_time as dummy arg
!       97-aug-14     ken dritz     removed include of lapsparms.for
!       97-aug-14     ken dritz     pass nx_l, ny_l, r_missing_data, and
!                                   laps_cycle_time to plot_cont
!       97-aug-14     ken dritz     pass nx_l, ny_l (a second time) and
!                                   r_missing_data, laps_cycle_time to
!                                   plot_barbs
!       97-aug-14     ken dritz     pass nx_l, ny_l, and laps_cycle_time to
!                                   plot_grid
!       97-aug-14     ken dritz     pass nx_l, ny_l, and laps_cycle_time to
!                                   plot_cldpcp_type
!       97-aug-14     ken dritz     pass nx_l, ny_l, and laps_cycle_time to
!                                   plot_stations
!       97-aug-17     ken dritz     pass r_missing_data to divergence
!       97-aug-25     steve albers  removed equivalence for uv_2d.
!                                   removed equivalence for slwc_int.
!                                   removed equivalence for slwc_2d.
!                                   removed /lapsplot_cmn1/ and /lapsplot_cmn2/
!       97-sep-24     john smart    added display funtionality for
!                                   polar orbiter (lrs).
!       98-mar-23        "          added lvd subdirectory flexibility.

        use mem_namelist, only: max_snd_grid,max_snd_levels
     1                         ,model_fcst_intvl,precip_cycle_time

        include 'trigd.inc'

        include 'constants.inc'

        include 'lapsplot.inc'

        real lat(nx_l,ny_l),lon(nx_l,ny_l),topo(nx_l,ny_l)

        real,  allocatable  :: static_grid(:,:)

        character*1 c_display, qtype, tunits, c_prodtype, clvl_soil
        character*1 cansw
        character*10 c_sat_plot
        character*13 filename
        character*14 a14_time
        character*3 c3_site
        character*4 c4_string
        character*5 c5_string,fcst_hhmm
        character*4 c4_log
        character*7 c7_string
        character*9 c9_string,a9_start,a9_end
        character*20 colortable, btemp_table
        character infile*255
        character*40 c_model
        character*200 new_dataroot
        character*40 vert_grid
        character i4_to_byte

        real clow,chigh,cint_ref
        data clow/-200./,chigh/+400/,cint_ref/10./

        integer idum1_array(nx_l,ny_l)
        integer contable(0:1,0:1)

        real dum1_array(nx_l,ny_l)

      ! used for "potential" precip type
        logical iflag_mvd,iflag_icing_index,iflag_cloud_type
     1         ,iflag_bogus_w
        logical iflag_snow_potential, l_plot_image, l_image
        logical l_low_fill, l_high_fill

        logical lmask_rqc_3d(nx_l,ny_l,1)
        real rqc(nx_l,ny_l,1)

        integer ibase_array(nx_l,ny_l)
        integer itop_array(nx_l,ny_l)

        character*3 c_field
        character*2 c_metacode
        character*4 c_type, c_type_i
        character*3 cstatic
        character*3 c_bkg
        character c19_label*19,c30_label*30,c_label*100

!       integer ity,ily,istatus
!       data ity/35/,ily/1010/

!       stuff to read in wind file
        integer kwnd
        parameter (kwnd = 3)
        real u_2d(nx_l,ny_l) ! wrt true north
        real v_2d(nx_l,ny_l) ! wrt true north
        real u_2d1(nx_l,ny_l,1) ! wrt true north
        real v_2d1(nx_l,ny_l,1) ! wrt true north
        real w_2d(nx_l,ny_l)
        real liw(nx_l,ny_l)
!       real helicity(nx_l,ny_l)
        real vas(nx_l,ny_l)
        real cint
        real uv_2d(nx_l,ny_l,2)

        real dir(nx_l,ny_l)
        real spds(nx_l,ny_l)

        real sndr_po(19,nx_l,ny_l)

        character*4 var_2d, var_2d_in
        character*150  directory
        character*31  ext
        character*20  units_2d
        character*4   lvl_coord_2d
        character*125 comment_2d
        character*9 comment_a,comment_b

!       for reading in radar data
        real dummy_array(nx_l,ny_l)
        real radar_array(nx_l,ny_l)
!       real radar_array_adv(nx_l,ny_l)

        real v_nyquist_in_a(max_radars)
        real rlat_radar_a(max_radars), rlon_radar_a(max_radars) 
        real rheight_radar_a(max_radars)
        integer i4time_radar_a(max_radars)
        integer n_vel_grids_a(max_radars)
        integer       ioffset(max_radars)
        integer       joffset(max_radars)

        logical l_offset
        parameter (l_offset = .false.)

        character*4 radar_name,radar_name_a(max_radars)
        character*31 ext_radar,ext_radar_a(max_radars)

!       real omega_3d(nx_l,ny_l,nz_l)
        real grid_ra_ref(nx_l,ny_l,nz_l,l_radars)

!       real grid_ra_vel(nx_l,ny_l,nz_l,max_radars)
!       real grid_ra_nyq(nx_l,ny_l,nz_l,max_radars)
        real, allocatable, dimension(:,:,:,:) :: grid_ra_vel
        real, allocatable, dimension(:,:,:,:) :: grid_ra_nyq

        integer idx_radar(max_radars)

        real grid_ra_ref_dum(1,1,1,1)
        real grid_ra_vel_dum(1,1,1,1)
        real field_3d(nx_l,ny_l,nz_l)
        real pres_3d(nx_l,ny_l,nz_l)

!       real lifted(nx_l,ny_l)
        real height_2d(nx_l,ny_l)
        real temp_2d(nx_l,ny_l)
        real tw_sfc_k(nx_l,ny_l)
        real td_2d(nx_l,ny_l)
        real pres_2d(nx_l,ny_l)
        real temp_3d(nx_l,ny_l,nz_l)
        real pressures_mb(nz_l)

!       real slwc_int(nx_l,ny_l)
        real column_max(nx_l,ny_l)
        integer i_array(nx_l,ny_l)

        real field2_2d(nx_l,ny_l)
        real cice_2d(nx_l,ny_l)
        real field_2d(nx_l,ny_l)
        real field_2d_buf(nx_l,ny_l)
        real field_2d_sum(nx_l,ny_l)
        real field_2d_diff(nx_l,ny_l)

        real snow_2d(nx_l,ny_l)
        real snow_2d_buf(nx_l,ny_l)
        real precip_2d(nx_l,ny_l)
        real precip_2d_buf(nx_l,ny_l)
        real accum_2d(nx_l,ny_l)

        real dx(nx_l,ny_l)
        real dy(nx_l,ny_l)

!       local variables used in
        logical l_mask(nx_l,ny_l)
        integer ipcp_1d(nz_l)

        integer iarg

        real cloud_cvr(nx_l,ny_l)
!       real cloud_2d(nx_l,ny_l)

        character*255 c_filespec_ra
        character*255 c_filespec_src
        data c_filespec_src/'*.src'/

        character*255 c_filespec
        character*255 cfname
        character*2   cchan

        include 'satellite_dims_lvd.inc'

        character*15  clvdvars(maxchannel)
        data clvdvars/'[svs, svn, alb]',
     1                '[s3a, s3c,    ]',
     1                '[s4a, s4c,    ]',
     1                '[s8a, s8w, s8c]',
     1                '[sca, scc,    ]',
     1                '[sca, scc,    ]'/

        logical lfndtyp

        logical lapsplot_pregen,l_precip_pregen,l_pregen,l_radar_read
        data lapsplot_pregen /.true./

!       real heights_3d(nx_l,ny_l,nz_l)

        real p_1d_pa(nz_l)
        real rh_2d(nx_l,ny_l)
        real sh_2d(nx_l,ny_l)

        real k_to_f, k_to_c, c_to_k
        real make_rh

        include 'laps_cloud.inc'
        include 'bgdata.inc'

        real clouds_3d(nx_l,ny_l,kcloud)
        real cld_pres(kcloud)

        common /supmp1/ dummy,part

!       i_image: whether this particular plot is an image
        common /image/ n_image, i_image 

        common /plot_field_cmn/ i_plotted_field

!       common /conre1/ioffp,spval,epsval,cntmin,cntmax,cntint,ioffm

        character asc9_tim*9, asc_tim_24*24
        character asc9_tim_r*9, a9time*9
        character asc9_tim_n*9
        character c9_radarage*9

        character*9   c_fdda_mdl_src(maxbgmodels)
        character*10  cmds
        character*10  c10_grid_fname
        character*200 c_dataroot

        character*6 c_vnt_units
        character*7 c_units_type

c       include 'satellite_dims_lvd.inc'
        include 'satellite_common_lvd.inc'

        data mode_lwc/2/

        integer i_first
        save i_first
        data i_first /1/

        ifield_found = 0
        i_plotted_field = 0

        call find_domain_name(c_dataroot,c10_grid_fname,istatus)

        call get_vertical_grid(vert_grid,istatus)

        ndim_read = 2

        ialloc_vel = 0
        i4time_temp = 0

        icen = nx_l/2+1
        jcen = ny_l/2+1

        c_vnt_units = namelist_parms%c_vnt_units
        c_units_type = namelist_parms%c_units_type

        plot_parms%iraster = namelist_parms%iraster

!       surface temperature ranges
        if(namelist_parms%l_discrete)then
            sfctf_h = 120.
            sfctc_h = 50.
        else
            sfctf_h = 125.
            sfctc_h = 50.
        endif

        sfctf_l = -60.
        sfctc_l = -50.

        sfctdf_h = 120.
        sfctdf_l = -20.

        sfctdc_h =  50.
        sfctdc_l = -30.

        btemp_l = -80.
        btemp_h = +40.
        btemp_table = 'linear'

!       surface wind range
        chigh_sfcwind = namelist_parms%chigh_sfcwind

!       cape range
        chigh_cape = namelist_parms%chigh_cape

!       upslope moisture flux
        umf_l =  -40.
        umf_h = +120.

        hel_l = -800.
        hel_h = +800.

        cloud_albedo_f = 1.00 
        cloud_albedo_a = 1.00

        call get_pres_3d(i4time_ref,nx_l,ny_l,nz_l,pres_3d,istatus) 

        call get_max_radar_files(max_radar_files,istatus)      

        i_overlay = 0
        n_image = 0
        jdot = 1   ! 1 = dotted county boundaries, 0 = solid county boundaries
        part = 0.9 ! for plotting routines
        igrid = 0

        lagt = 10800

        ioffm = 1 ! don't plot label stuff in conrec

!       get fdda_model_source from parameter file
        call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models,istatus)

        ext = 'static'

!       get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

        call s_len(c10_grid_fname,lf)
        if(c10_grid_fname(1:lf).eq.'nest7grid')then
           var_2d='lat'
        else
           var_2d='lac'   !wrfsi c-stagger
        endif
        call read_static_grid(nx_l,ny_l,var_2d,lat,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading laps static-lat'
            return
        endif

        if(c10_grid_fname(1:lf).eq.'nest7grid')then
           var_2d='lon'
        else
           var_2d='loc'   !wrfsi c-stagger
        endif 
        call read_static_grid(nx_l,ny_l,var_2d,lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading laps static-lon'
            return
        endif
        if(c10_grid_fname(1:lf).eq.'nest7grid')then
           var_2d='avg'
        else
           var_2d='avc'   !wrfsi c-stagger
        endif

        call read_static_grid(nx_l,ny_l,var_2d,topo,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading domain static-topo'
            return
        endif

        i_balance = 0

1200    continue

        plot_parms%l_discrete = namelist_parms%l_discrete

        write(6,11)
11      format(//'  select field: (append with "i" for image)',       
     1      /'     [wd,wb,wr,wf,bw] wind'
     1      ,' (lw3/lwm, lga/lgb, fua/fsf, laps-bkg, bal), '
     1      /'     [wo,co,bo,lo,fo] '      
     1      ,'anlyz/cloud/balance/bkg/fcst omega'      
     1      /'     radar: [ra] intermediate vrc, [rf] anal/fcst fields'       
     1      /'            radar intermediate vxx - ref [rv], vel [rd]'     
     1      /
     1      /'     sfc: [p,pm,ps,al,pp,tf,tc,df,dc,ws,wp,vv,hu,ta'          
     1      ,',th,te,vo,mr,mc,dv,ha,ma]'
     1      /'          [sp,cs,vs,tw,fw,hi,gf]'
     1      /'          [of,oc,ov,os,op,og,qf,qc,qv,qs,qp,qg] obs plots'    
     1      /'          [st,mw] obs/mesowx locations'    
     1      /'          [bs] sfc background/forecast, '
     1                ,'[by,bn] balance toggle'
     1      /10x,'[li,lw,he,pe,ne,um] li, li*w, helcty, cape, cin, umf'
     1      /10x,'[s] other stability indices, [sm] soil moisture'
     1      /
     1      /'     temp: [t,tb,tr,to,bt] (laps,lga,fua,obs,bal)'       
     1      ,', [pt,pb] theta, bal theta'
     1      /'     hgts: [ht,hb,hr,bh] (laps,lga,fua,bal),'
     1      /'           [hh,bl,lf] ht of temp sfc, pbl, fire wx')

        write(6,12)
 12     format(
     1       /'     humidity: [br,fr,lq,rb] (lga;fua;lq3;bal)'       
     1       /'               [ho,qo,wv](tdobs;qobs;pwobs)'            
     1       /'               [pw] precipitable h2o'            
     1       /
     1       /'     clouds/precip: [ci] cloud ice,'
     1       ,' [ls] cloud lwc'
     1       /'         [is] integrated cloud lwc  '
     1       /'         [mv] mean volume drop diam,   [ic] icing index,'       
     1       /'         [cc] cld ceiling (agl),'
     1       ,' [cb,ct] cld base/top (msl)'      
     1       /'         [cv/cg] cloud cover (2-d)'
     1       ,' [cy,py] cloud/precip type'
     1       /'         [pc,rn,sn,pi] pcp conc, [sa/pa] snow/pcp accum,'      
     1       ,' [sc/csc] snow cvr'
     1      //'     static info: [gg] '
     1       /'     [lv(d),lr(lsr),v3,v5,po,lc] lvd; lsr; vcf; '             
     1       ,'tsfc-11u; polar orbiter, lcv'
     1      //'     difference field: [di] '
     1      //' ',52x,'[q] quit/display ? ',$)

 15     format(a3)

        read(lun,16)c_type
 16     format(a4)

        call s_len(c_type,len_type)
        if(len_type .eq. 0)goto1200

        if(c_type(len_type:len_type) .eq. 'i' 
     1                  .and. c_type .ne. 'di'
     1                  .and. c_type .ne. 'li'
     1                  .and. c_type .ne. 'cwi'
     1                  .and. c_type .ne. 'ci')then
            i_image = 1
            c_type_i = c_type(1:len_type-1)
        else
            i_image = 0
            c_type_i = c_type
        endif

        if(c_type .eq. 'm')then
            write(6,*)' plot just map background'
            ifield_found = 1
            call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                         ,namelist_parms,plot_parms)
            goto 1200
        endif

        write(6,*)' c_type is: ',c_type

        if(c_type(1:2) .eq. 'by')then
            i_balance = 1
            goto 1200
        elseif(c_type(1:2) .eq. 'bn')then
            i_balance = 0
            goto 1200
        endif

        plot_parms%color_power = 1.

        if(c_type(1:2) .eq. 'fc')then ! force config with new dataroot
            write(6,*)' time is: ',asc9_tim
            write(6,*)' enter new dataroot:'
            read(lun,17)new_dataroot
 17         format(a)
            call s_len(new_dataroot,lenroot)
            call force_get_laps_config(new_dataroot(1:lenroot),istatus)
            if(istatus .ne. 1)then
                write(6,*)' bad status returned from force_laps_config'
                return
            endif
            call get_lapsplot_parms(namelist_parms,istatus)
            if(c_type(1:3) .eq. 'fcf')then
                call frame
!               ifield_found = 0 ! latest experiment
                i_overlay = 0
                n_image = 0
            endif
            goto 1200
        endif

        if(c_type(1:2) .eq. 'di' .or. c_type(1:2) .eq. 'ml' .or. 
     1     c_type(1:2) .eq. 'dt')then

            if(c_type(1:2) .eq. 'di')then
              write(6,*)' plotting difference field of last two entries' 
     1                 ,asc9_tim
              call diff_miss(field_2d,field_2d_buf,field_2d_diff
     1                                            ,nx_l,ny_l)       

              c_label = 'difference field (b - a)'
              colortable = 'hues'

            elseif(c_type(1:2) .eq. 'dt')then
              if(c_type(3:3) .eq. '4')then
                  thresh = 40.
              elseif(c_type(3:3) .eq. '3')then
                  thresh = 30.
              else
                  thresh = 20.
              endif

              write(6,*)
     1         ' plotting contingency table of last two entries', thresh 

              lmask_rqc_3d = .true.

              ext = 'lcv'
              call get_laps_2d(i4time_ref,ext,'rqc',units_2d
     1                        ,comment_2d,nx_l,ny_l,rqc,istatus)
              if(istatus .ne. 1)then
                  write(6,*)' error reading 2d rqc analysis'
              else
                  write(6,*)' apply rqc to mask (under construction)'
                  where(rqc .ne. 3.0)lmask_rqc_3d = .false.
              endif

              call calc_contable_3d(
     1             field_2d_buf,field_2d,thresh,nx_l,ny_l,1,    ! i
     1             lmask_rqc_3d,r_missing_data,                 ! i
     1             field_2d_diff)                               ! o

              plot_parms%iraster = +1
              plot_parms%l_discrete = .true.
              colortable = 'cont' ! 'ref'                 
              where(field_2d_diff .eq. 0.0)field_2d_diff = +1.4 ! hit      
              where(field_2d_diff .eq. 1.0)field_2d_diff = +0.7 ! miss      
              where(field_2d_diff .eq. 2.0)field_2d_diff = +2.2 ! false pos
              where(field_2d_diff .eq. 3.0)field_2d_diff = +0.0 ! correct neg
              where(field_2d_diff .eq. r_missing_data)field_2d_diff=.100 ! outside mask
              dyn_low = 0.0
              dyn_high = 3.0

              lun_out = 6
              ilow = 1
              ihigh = nx_l
              jlow = 1
              jhigh = ny_l
              call contingency_table(field_2d_buf,field_2d       ! i
     1                              ,nx_l,ny_l,1                 ! i
     1                              ,thresh,thresh               ! i
     1                              ,lun_out                     ! i
     1                              ,ilow,ihigh,jlow,jhigh       ! i
     1                              ,lmask_rqc_3d                ! i
     1                              ,contable)                   ! o

              call skill_scores(contable,lun_out                   ! i
     1              ,frac_coverage                                 ! o
     1              ,frac_obs                                      ! o
     1              ,frac_fcst                                     ! o
     1              ,bias                                          ! o
     1              ,ets)                                          ! o

              write(c_label,41)nint(thresh),bias,ets
 41           format(i2,'dbz contingency table (b-a) bias =',f5.2
     1                                             ,' ets =',f6.3)

            else
              write(6,*)' plotting product field of last two entries' 
              call multar_miss(field_2d,field_2d_buf,field_2d_diff
     1                                            ,nx_l,ny_l)       

              c_label = 'product field (b - a)'
              colortable = 'hues'
            endif

!           use scale from the most recent plot?
!           scale = 1.
            write(6,*)' scale/cint: ',scale,cint   

            if(.false.)then ! experimental
                colortable = 'hues'
                call plot_field_2d(i4time_3dw,c_type,field_2d_diff,scale       
     1                            ,namelist_parms,plot_parms
     1                            ,clow,chigh,cint,c_label
     1                            ,i_overlay,c_display,lat,lon,jdot
     1                            ,nx_l,ny_l,r_missing_data,colortable)       

            elseif(c_type(3:3) .ne. 'i' .and. c_type(4:4) .ne. 'i')then
              ! contour plot
                call contour_settings(field_2d_diff,nx_l,ny_l
     1               ,clow,chigh,cint,zoom,density,scale)      

                call plot_cont(field_2d_diff,scale,clow,chigh,cint,
     1               asc9_tim,namelist_parms,plot_parms,       
     1               c_label,i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)

            else ! image plot
!               plot_parms%iraster = -1

                if(dyn_low  .eq. r_missing_data .or. 
     1             dyn_high .eq. r_missing_data)then ! initialize range
                    call array_range(field_2d_diff,nx_l,ny_l,rmin,rmax
     1                              ,r_missing_data)

                    rmin = rmin/scale
                    rmax = rmax/scale

                    rscale = max(abs(rmin),abs(rmax))
                    rmin = -rscale
                    rmax = +rscale
                    dyn_low  = rmin ! save for subsequent times
                    dyn_high = rmax ! save for subsequent times

                else ! use previously initialized values
                    rmin = dyn_low
                    rmax = dyn_high

                endif

                write(6,*)' ccpfil for diff plot range = ',rmin,rmax
     1                                                    ,scale

                call ccpfil(field_2d_diff,nx_l,ny_l,rmin,rmax,colortable
     1                     ,n_image,scale,'hsect',plot_parms
     1                     ,namelist_parms)    
                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum('in',7)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)     

            endif

        elseif(c_type(1:2) .eq. 'mn')then
            write(6,*)' plotting sum field of last two entries'       
            call diff_miss(field_2d,field_2d_buf,field_2d_diff
     1                                          ,nx_l,ny_l)       

            c_label = 'sum field (b + a)'

!           use scale from the most recent plot?
!           scale = 1.
            write(6,*)' scale/cint: ',scale,cint   

            if(.false.)then ! experimental
                colortable = 'hues'
                call plot_field_2d(i4time_3dw,c_type,field_2d_diff,scale       
     1                            ,namelist_parms,plot_parms
     1                            ,clow,chigh,cint,c_label
     1                            ,i_overlay,c_display,lat,lon,jdot
     1                            ,nx_l,ny_l,r_missing_data,colortable)       

            elseif(c_type(3:3) .ne. 'i')then ! contour plot
                call contour_settings(field_2d_diff,nx_l,ny_l
     1               ,clow,chigh,cint,zoom,density,scale)      

                call plot_cont(field_2d_diff,scale,clow,chigh,cint,
     1               asc9_tim,namelist_parms,plot_parms,       
     1               c_label,i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)

            else ! image plot
!               plot_parms%iraster = -1

                if(dyn_low  .eq. r_missing_data .or. 
     1             dyn_high .eq. r_missing_data)then ! initialize range
                    call array_range(field_2d_diff,nx_l,ny_l,rmin,rmax
     1                              ,r_missing_data)

                    rmin = rmin/scale
                    rmax = rmax/scale

                    rscale = max(abs(rmin),abs(rmax))
                    rmin = -rscale
                    rmax = +rscale
                    dyn_low  = rmin ! save for subsequent times
                    dyn_high = rmax ! save for subsequent times

                else ! use previously initialized values
                    rmin = dyn_low
                    rmax = dyn_high

                endif

                write(6,*)' ccpfil for diff plot range = ',rmin,rmax
     1                                                    ,scale

                call ccpfil(field_2d_diff,nx_l,ny_l,rmin,rmax,'hues'
     1                     ,n_image,scale,'hsect',plot_parms
     1                     ,namelist_parms)    
                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum('in',7)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)     

            endif

        else ! non-difference plot
            if(igrid .eq. 1)then
                write(6,*)' copying field_2d to field_buf for diff optn'       
                call move(field_2d,field_2d_buf,nx_l,ny_l)       
            endif

            scale = 1.0 ! initialize to default value

        endif

        igrid = 1

        if(    c_type_i      .eq. 'wd' .or. c_type_i      .eq. 'wb'  ! wind fields
     1    .or. c_type_i(1:2) .eq. 'co' .or. c_type_i      .eq. 'wr'
     1    .or. c_type_i      .eq. 'wf' .or. c_type_i(1:2) .eq. 'bw'
     1    .or. c_type_i(1:2) .eq. 'bo' .or. c_type_i      .eq. 'lo'
     1    .or. c_type_i(1:2) .eq. 'fo' .or. c_type_i(1:2) .eq. 'wo')then       

            if(c_type_i .eq. 'wd')then
                ext = 'lwm'

            elseif(c_type_i .eq. 'wb')then
                call make_fnam_lp(i4time_ref,asc9_tim,istatus)
                ext = 'lga'

            elseif(c_type_i .eq. 'wr')then
                ext = 'fua'

            elseif(c_type_i(1:2) .eq. 'co')then
                ext = 'lco'

            elseif(c_type_i(1:2) .eq. 'wo')then
                ext = 'lw3'

            elseif(c_type_i .eq. 'lo')then
                ext = 'lga'

            elseif(c_type_i .eq. 'fo')then
                ext = 'fua'

            elseif(c_type_i .eq. 'wf')then
                ext = 'lw3'
            elseif(c_type_i .eq. 'bw'.or.c_type_i.eq.'bo')then
                ext = 'balance'
            endif


            if(c_type_i .eq. 'wd')then
                write(6,13)
13              format(
     1    '     enter level in mb, -1 = mean, 0 = sfc',24x,'? ',$)
            else
                write(6,14)
14              format(
     1    '     enter level in mb, -1 = pbl mean, 0 = sfc',20x,'? ',$)
            endif

            call input_level(lun,k_level,k_mb,pres_3d,nx_l,ny_l,nz_l)       

!!! move and add by huiling yuan and steve albers to avoid 'lwm' empty direcotry, may 25, 2010
            if(c_type_i.eq.'bw'.or.c_type_i.eq.'bo')then
               call get_filespec(ext,1,c_filespec,istatus)
               ext='lw3'
               call s_len(c_filespec,ilen)
               c_filespec=c_filespec(1:ilen)//ext(1:3)//'/*.'//ext
            else
               if(k_level .ge. 1 .and. ext(1:3) .eq. 'lwm')then
                 ext='lw3'
               endif
               call get_filespec(ext,2,c_filespec,istatus)
            endif
!!!  end move block to avoid 'lwm' empty directory 

            if(c_type_i.ne.'lo' .and. c_type_i .ne. 'fo' 
     1                          .and. c_type_i .ne. 'wr'
     1                          .and. c_type_i .ne. 'wb')then
               write(6,*)
               write(6,*)'    looking for laps wind data: ',ext(1:3)
               call get_file_time(c_filespec,i4time_ref,i4time_3dw)

            else 
               i4time_3dw = i4time_ref

            endif

            call make_fnam_lp(i4time_3dw,asc9_tim,istatus)

            if(c_type_i.eq.'bw' .or. c_type_i.eq.'bo')ext='balance'

            if(c_type_i(1:2) .eq.'co' .or. c_type_i(1:2).eq.'bo' .or.
     1         c_type_i(1:2) .eq.'lo' .or. c_type_i(1:2).eq.'fo' .or. 
     1         c_type_i(1:2) .eq.'wo'                         )then
                c_field = 'w'
                goto115
            endif

            if(k_level .gt. -1)then

                if(k_level .eq. 0)then ! sfc winds
                    write(6,102)
102                 format(/
     1                  '  field [di,sp,u,v,dv,vc (barbs),ob (obs)]'
     1                  ,27x,'? ',$)
                    read(lun,15)c_field

                    if(c_type_i .eq. 'wd')then
                        ext = 'lwm'
                    elseif(c_type_i .eq. 'wb')then
                        ext = 'lgb'
                    elseif(c_type_i .eq. 'wr')then
                        ext = 'fsf'
                    endif

                    if(ext(1:3) .eq. 'lgb' .or. ext(1:3) .eq. 'fsf')then       
                        call input_background_info(
     1                              ext                     ! i
     1                             ,directory,c_model       ! o
     1                             ,i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,asc9_tim            ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o
                        if(istatus.ne.1)goto1200

                        level=0

                        var_2d = 'usf'

                        write(6,*)' reading sfc wind data from: '
     1                            ,ext(1:3),' ',var_2d

                        call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,level,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 u_2d,istatus)

                        if(istatus.ne.1)goto1200

                        var_2d = 'vsf'

                        write(6,*)' reading sfc wind data from: '
     1                            ,ext(1:3),' ',var_2d

                        call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,level,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 v_2d,istatus)

                        i4time_3dw = i4_valid
                        call make_fnam_lp(i4time_3dw,asc9_tim
     1                                   ,istatus)     
                        write(6,*)' valid time = ',asc9_tim

                    elseif(c_type_i.eq.'bw')then
                        write(6,*)' balanced lsx was chosen'

                        call get_directory('balance',directory,lend)
                        directory=directory(1:lend)//'lsx/'

                        ext = 'lsx'

                        var_2d = 'u'
                        call get_2dgrid_dname(directory,i4time_3dw
     1                    ,laps_cycle_time*100,i4time_heights,ext,var_2d      
     1                    ,units_2d,comment_2d,nx_l,ny_l,u_2d,0
     1                    ,istatus)   
                        if(istatus .ne. 1)goto 1200   

                        var_2d = 'v'
                        call get_2dgrid_dname(directory,i4time_3dw
     1                    ,laps_cycle_time*100,i4time_heights,ext,var_2d      
     1                    ,units_2d,comment_2d,nx_l,ny_l,v_2d,0
     1                    ,istatus)      
                        if(istatus .ne. 1)goto 1200   

                    else ! lwm (unbalanced)
!                       call get_directory(ext,directory,len_dir)
!                       c_filespec = directory(1:len_dir)//'*.'//ext(1:3)      

                        call get_filespec(ext,2,c_filespec,istatus)

                        var_2d = 'su'
                        call get_laps_2d(i4time_3dw,ext,var_2d
     1                      ,units_2d,comment_2d,nx_l,ny_l,u_2d,istatus)

                        var_2d = 'sv'
                        call get_laps_2d(i4time_3dw,ext,var_2d
     1                      ,units_2d,comment_2d,nx_l,ny_l,v_2d,istatus)

                    endif

                else if(k_level .gt. 0)then
                    write(6,103)
103                 format(/
     1                       '  field [di,sp,u,v,om,dv,vo,pv,va,vc'      
     1                      ,' (barbs), ob (obs))]'   
     1                                          ,14x,'? ',$)
                    read(lun,15)c_field

                    if(ext(1:3) .eq. 'lwm')then
                      ext = 'lw3'
                    endif

                    write(6,*)' ext = ',ext

                    call make_fnam_lp(i4time_3dw,asc9_tim,istatus)      
                    write(6,*)' valid time = ',asc9_tim

                    if(c_field .ne. 'w ' .and. c_field .ne. 'ob')then
                      write(6,*)' calling get_uv_2d for ',ext
                      call get_uv_2d(i4time_3dw,k_level,uv_2d,ext
     1                             ,nx_l,ny_l,fcst_hhmm,c_model,istatus)      
                      call make_fnam_lp(i4time_3dw,asc9_tim,istatus)

                      if(c_type_i .eq. 'wf')then

!                       calculate wind difference vector (lw3 - model first guess)
                        var_2d = 'u3'
                        call get_modelfg_3d(i4time_3dw,var_2d
     1                           ,nx_l,ny_l,nz_l,field_3d,istatus) 
                        call multcon(field_3d(1,1,k_level),-1.
     1                        ,nx_l,ny_l)      
                        call add(field_3d(1,1,k_level),uv_2d(1,1,1),u_2d
     1                        ,nx_l,ny_l)      

                        var_2d = 'v3'
                        call get_modelfg_3d(i4time_3dw,var_2d
     1                           ,nx_l,ny_l,nz_l,field_3d,istatus) 
                        call multcon(field_3d(1,1,k_level),-1.
     1                                       ,nx_l,ny_l)      
                        call add(field_3d(1,1,k_level),uv_2d(1,1,2),v_2d       
     1                                   ,nx_l,ny_l)      

                      else ! c_type_i .ne. 'wf'
                        call move(uv_2d(1,1,1),u_2d,nx_l,ny_l)
                        call move(uv_2d(1,1,2),v_2d,nx_l,ny_l)

                      endif ! c_type_i .eq. 'wf'
 
                    endif ! c_field = 'w'

                    if(c_field .eq. 'pv')then ! read 3d temperature field
                      if(c_type_i .eq. 'wd')then
                        ext = 'lt1'
                        call get_temp_3d(i4time_ref,i4time_nearest
     1                                  ,1,nx_l,ny_l,nz_l
     1                                  ,temp_3d,istatus)
                        if(istatus .ne. 1)goto1200
                      elseif(c_type_i.eq.'bw')then
                        ext = 'lt1'
                        call get_temp_3d(i4time_ref,i4time_nearest
     1                                  ,4,nx_l,ny_l,nz_l
     1                                  ,temp_3d,istatus)
                        if(istatus .ne. 1)goto1200
                      elseif(c_type_i .eq. 'wb')then
                        ext = 'lga'
                        var_2d = 't3'
                        call input_background_info(
     1                              ext                     ! i
     1                             ,directory,c_model       ! o
     1                             ,i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,asc9_tim                ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o
                        if(istatus.ne.1)goto1200 
                        call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,temp_3d
     1                              ,istatus)
                        if(istatus .ne. 1)then
                          write(6,*)' error reading grid ',var_2d,' '
     1                                                    ,ext,istatus       
                          goto1200
                        endif
                      elseif(c_type_i .eq. 'wr')then
                        ext = 'fua'
                        var_2d = 't3'
                        call input_background_info(
     1                              ext                     ! i
     1                             ,directory,c_model       ! o
     1                             ,i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,asc9_tim                ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o
                        if(istatus.ne.1)goto1200
                        call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,temp_3d
     1                              ,istatus)
                        if(istatus .ne. 1)then
                          write(6,*)' error reading grid ',var_2d,' '
     1                                                    ,ext,istatus       
                          goto1200
                        endif
                      endif

                    endif ! pv field temperature read

                endif ! k_level > 0

            elseif(k_level .eq. -1)then ! read mean winds from 2d grids

                if(c_type_i .eq. 'wd')then
                    ext = 'lwm'
                elseif(c_type_i .eq. 'wb')then
                    ext = 'lgb'
                elseif(c_type_i .eq. 'wr')then
                    ext = 'fsf'
                endif

                if(ext(1:3) .eq. 'lgb' .or. ext(1:3) .eq. 'fsf')then       
                    call input_background_info(
     1                              ext                     ! i
     1                             ,directory,c_model       ! o
     1                             ,i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,asc9_tim                ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o
                    if(istatus.ne.1)goto1200

                    level=0

                    var_2d = 'upb'

                    write(6,*)' reading pbl wind data from: '
     1                            ,ext(1:3),' ',var_2d

                    call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,level,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 u_2d,istatus)

                    if(istatus.ne.1)goto1200

                    var_2d = 'vpb'

                    write(6,*)' reading pbl wind data from: '
     1                            ,ext(1:3),' ',var_2d

                    call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,level,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 v_2d,istatus)

                    i4time_3dw = i4_valid
                    write(6,*)' valid time = ',asc9_tim

                else ! ext = lwm
                    write(6,*)' getting lwm mean wind file'

                    ext = 'lwm'
                    call get_directory(ext,directory,len_dir)
                    var_2d = 'mu'
                    call get_laps_2d(i4time_3dw,ext,var_2d
     1              ,units_2d,comment_2d,nx_l,ny_l,u_2d,istatus)
                    var_2d = 'mv'
                    call get_laps_2d(i4time_3dw,ext,var_2d
     1              ,units_2d,comment_2d,nx_l,ny_l,v_2d,istatus)

                endif ! ext

                write(6,104)
104             format(/'  field [di,sp,u,v,vc (barbs)]   ',25x,'? ',$)
                read(lun,15)c_field

            endif ! k_level

!  ***      display wind data  ******************************************************

115         if(c_field(1:2) .eq. 'di' .or. c_field(1:2) .eq. 'sp')then
                do i = 1,nx_l
                do j = 1,ny_l
                    if(u_2d(i,j) .eq. r_missing_data
     1            .or. v_2d(i,j) .eq. r_missing_data)then
                        dir(i,j)  = r_missing_data
                        spds(i,j) = r_missing_data
                    else
                        call uvgrid_to_disptrue(u_2d(i,j),
     1                                  v_2d(i,j),
     1                                  dir(i,j),
     1                                  spds(i,j),
     1                                  lat(i,j),
     1                                  lon(i,j)     )
                        spds(i,j) = spds(i,j) / mspkt
                    endif
                enddo ! j
                enddo ! i

                field_2d = spds ! support for diff option

            endif

            if(c_field(1:2) .eq. 'di')then
                c19_label = ' isogons   (deg)   '
                call mklabel(k_mb,c19_label,c_label)

                call plot_cont(dir,1e0,clow,chigh,30.,asc9_tim,
     1              namelist_parms,plot_parms,c_label,i_overlay,
     1              c_display,lat,lon,jdot,       
     1              nx_l,ny_l,r_missing_data,laps_cycle_time)

            else if(c_field(1:2) .eq. 'sp')then
                if(c_type_i .eq. 'wb' .or. c_type_i .eq. 'wr')then
                    call mk_fcst_hlabel(k_mb,'isotachs',fcst_hhmm       
     1                                 ,ext(1:3),'kt'
     1                                 ,c_model,c_label)

                else
                    c19_label = ' isotachs (kt) '//ext(1:3)
                    if(c_type_i .eq. 'bw')then
                        c19_label = c19_label(1:18)//'b'
                    endif

                    call mklabel(k_mb,c19_label,c_label)
     
                endif

                if(k_level .gt. 0 .and. k_mb .lt. 500)then
                    cint = 10.
                    clow = 0.
                    chigh = 200.
                elseif(k_level .gt. 0)then
                    cint = 5.
                    clow = 0.
                    chigh = 100.
                else ! k_level = 0 (sfc)
                    cint = 5.
                    clow = 0.
                    chigh = chigh_sfcwind
                endif

                scale = 1.0

                call plot_field_2d(i4time_3dw,c_field,spds,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'spectral')

            else if(c_field(1:1) .eq. 'u' .or. 
     1              c_field(1:1) .eq. 'u')then
                if(c_type_i .eq. 'wf')then
                    c19_label = ' u  diff        m/s'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(c_type_i .eq. 'wb')then
                    call mk_fcst_hlabel(k_mb,'u',fcst_hhmm
     1                                 ,ext(1:3),'m/s'
     1                                 ,c_model,c_label)
                elseif(c_type_i .eq. 'wr')then
                    call mk_fcst_hlabel(k_mb,'u',fcst_hhmm
     1                                 ,ext(1:3),'m/s'
     1                                 ,c_model,c_label)
                elseif(c_type_i .eq. 'bw')then
                    c19_label = ' u  (bal)       m/s'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(ext(1:3) .eq. 'lwm')then
                    c19_label = ' u  (lwm)       m/s'
                    call mklabel(k_mb,c19_label,c_label)
                else
                    c19_label = ' u  (anal)      m/s'
                    call mklabel(k_mb,c19_label,c_label)
                endif

                if(k_level .gt. 0 .and. k_mb .lt. 500)then
                    cint = 10.
                    clow = -100.
                    chigh = 100.
                else
                    cint = 5.
                    clow =  -25.
                    chigh = +25.
                endif

                scale = 1.0
!               call contour_settings(u_2d,nx_l,ny_l,clow,chigh,cint
!    1                               ,zoom,density,scale)

                call plot_field_2d(i4time_3dw,c_field,u_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

                call move(u_2d,field_2d,nx_l,ny_l)

            else if(c_field .eq. 'v' .or. c_field .eq. 'v' .or. 
     1                                    c_field .eq. 'vi')then
                if(c_type_i .eq. 'wf')then
                    c19_label = ' v  diff        m/s'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(c_type_i .eq. 'wb')then
                    call mk_fcst_hlabel(k_mb,'v',fcst_hhmm
     1                                 ,ext(1:3),'m/s'
     1                                 ,c_model,c_label)
                elseif(c_type_i .eq. 'wr')then
                    call mk_fcst_hlabel(k_mb,'v',fcst_hhmm
     1                                 ,ext(1:3),'m/s'
     1                                 ,c_model,c_label)
                elseif(c_type_i .eq. 'bw')then
                    c19_label = ' v  (bal)       m/s'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(ext(1:3) .eq. 'lwm')then
                    c19_label = ' v  (lwm)       m/s'
                    call mklabel(k_mb,c19_label,c_label)
                else
                    c19_label = ' v  (anal)      m/s'
                    call mklabel(k_mb,c19_label,c_label)
                endif

                if(k_level .gt. 0 .and. k_mb .lt. 500)then
                    cint = 10.
                    clow = -100.
                    chigh = 100.
                else
                    cint = 5.
                    clow =  -25.
                    chigh = +25.
                endif

                scale = 1.0
!               call contour_settings(v_2d,nx_l,ny_l,clow,chigh,cint
!    1                               ,zoom,density,scale)

                call plot_field_2d(i4time_3dw,c_field,v_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

                call move(v_2d,field_2d,nx_l,ny_l)

            else if(c_field .eq. 'vc' .or. c_field .eq. 'ob')then
                if(c_type_i .eq. 'wf')then
!                   c19_label = ' wind diff (kt)    '
!                   call mklabel(k_mb,c19_label,c_label)

                    call mk_fcst_hlabel(k_mb,'wind diff',fcst_hhmm
     1                                 ,ext(1:3),'kt'
     1                                 ,c_model,c_label)

                elseif(c_type_i.eq.'wb')then
!                   c19_label = ' wind '//ext(1:3)//' '
!    1                                      //fcst_hhmm//'   kt'
!                   call mklabel(k_mb,c19_label,c_label)

                    call mk_fcst_hlabel(k_mb,'wind',fcst_hhmm
     1                                 ,ext(1:3),'kt'
     1                                 ,c_model,c_label)

                elseif(c_type_i.eq.'wr')then
                    if(k_level .eq. -1)then
                        call mk_fcst_hlabel(k_mb,'pbl mean wind'
     1                                 ,fcst_hhmm     
     1                                 ,ext(1:3),'kt'
     1                                 ,c_model,c_label)
                    elseif(k_level .eq. 0)then
                        call mk_fcst_hlabel(k_mb,'sfc wind',fcst_hhmm
     1                                 ,ext(1:3),'kt'
     1                                 ,c_model,c_label)
                    else
                        call mk_fcst_hlabel(k_mb,'wind',fcst_hhmm
     1                                 ,ext(1:3),'kt'
     1                                 ,c_model,c_label)
                    endif

                elseif(c_type_i.eq.'bw'                       )then
                    c19_label = ' wind  (bal)    kt'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(c_field .eq. 'ob')then
                    c19_label = ' wind  (obs)    kt'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(ext(1:3) .eq. 'lwm')then
                    c19_label = ' wind  (lwm)    kt'
                    call mklabel(k_mb,c19_label,c_label)
                else
                    c19_label = ' wind  (anl)    kt'
                    call mklabel(k_mb,c19_label,c_label)
                endif

                if(k_level .ne. 0)then ! upper air
                    nxz = (float(nx_l) /1.2) / sqrt(zoom)
                    nyz = float(ny_l) / sqrt(zoom)
                    ngz = max(nxz,nyz)
                    interval = (ngz / 50) + 1
                    size = float(interval) * .11
                else                   ! sfc
                    nxz = float(nx_l) / zoom
                    nyz = float(ny_l) / zoom
                    interval = int(max(nxz,nyz) / 85.) + 1
                    size = float(interval) * .15

                endif

                if(c_field .eq. 'ob')then
                    igrid = 0
                    if(i4time_3dw .eq. 0)then
                        i4time_3dw = i4time_ref
                        call make_fnam_lp(i4time_3dw,asc9_tim
     1                                   ,istatus)
                    endif
                    write(6,*)' wind barb ob size = ',size
                endif

                call plot_barbs(u_2d,v_2d,lat,lon,topo,size,zoom
     1               ,interval,asc9_tim,namelist_parms,plot_parms      
     1               ,c_label,c_field,k_level,i_overlay,c_display       
     1               ,nx_l,ny_l,nz_l,max_radars
!    1               ,grid_ra_ref_dum,grid_ra_vel_dum       
     1               ,nx_l,ny_l,r_missing_data,laps_cycle_time,jdot)

            else if(c_field .eq. 'w' .or. c_field .eq. 'om')then ! omega 
                if(c_type_i(1:2) .eq. 'co')then
                    write(6,*)' reading cloud omega'
                    var_2d = 'com'
                    ext = 'lco'
                    call get_laps_2dgrid(i4time_3dw,0,i4time_nearest
     1                                  ,ext,var_2d,units_2d,comment_2d
     1                                  ,nx_l,ny_l,w_2d,k_mb,istatus)
                    call mklabel(k_mb
     1                     ,' cloud omega 0.1ubar/s',c_label)       

                    i4_valid = i4time_nearest

                elseif(c_type_i(1:2) .eq. 'wo')then
                    write(6,*)' reading lw3 omega'
                    var_2d = 'om'
                    ext = 'lw3'
                    call get_laps_2dgrid(i4time_3dw,0,i4time_nearest,
     1              ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                          ,w_2d,k_mb,istatus)
                    call mklabel(k_mb
     1                     ,' anlyz omega 0.1ubar/s',c_label)       

                    i4_valid = i4time_nearest

                else if(c_type_i(1:2) .eq. 'bo')then
                    write(6,*)' reading balanced omega'
                    var_2d = 'om'
                    ext = 'lw3'
                    call get_directory('balance',directory,lend)
                    directory=directory(1:lend)//'lw3/'
                    call get_2dgrid_dname(directory,i4time_3dw
     1              ,laps_cycle_time*100,i4time_heights,ext,var_2d
     1              ,units_2d,comment_2d,nx_l,ny_l,w_2d,k_mb,istatus)       

                    call mklabel(k_mb
     1                     ,' balnc omega 0.1ubar/s',c_label)       

                    i4_valid = i4time_3dw

                else if(c_type_i(1:2) .eq. 'lo' .or. 
     1                  c_type_i(1:2) .eq. 'fo')then
                    call input_background_info(
     1                              ext                     ! i
     1                             ,directory,c_model       ! o
     1                             ,i4time_3dw              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,asc9_tim                ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o
                   if(istatus.ne.1)goto1200

                   write(6,211)ext(1:3)
                   var_2d = 'om'

                   call read_laps(i4_initial,i4_valid,directory,ext
     1             ,nx_l,ny_l,1,1,var_2d,k_mb,lvl_coord_2d,units_2d
     1             ,comment_2d,w_2d,istatus)

                   call s_len(c_model,len_model)
                   if(ext(1:3) .eq. 'fua' .and. len_model .gt. 0)then
!                      call mklabel(k_mb,' '//fcst_hhmm
!    1                         //' '//c_model(1:len_model)
!    1                         //' om ubar/s',c_label)

                       call mk_fcst_hlabel(k_mb,'omega',fcst_hhmm
     1                                 ,ext(1:3),'0.1ubar/s'
     1                                 ,c_model,c_label)
                   else
                       call mklabel(k_mb,' '//fcst_hhmm
     1                    //' '//ext(1:3)//' om 0.1ubar/s',c_label)
                   endif
                endif

                do j = 1,ny_l
                do i = 1,nx_l
                    if(w_2d(i,j) .eq. r_missing_data)then
                        w_2d(i,j) = 0.
                    endif
                enddo ! i
                enddo ! j

                scale = 1e-2
                if(i_image .eq. 0)then
                    cint = -1. * 2. ** (-density)
                else
                    cint = -1.0
                endif

!               if(plot_parms%iraster .lt. 1)then
!                   plot_parms%l_discrete = .true.
!               endif

                call get_grid_spacing_cen(grid_spacing_m,istatus)
                if(grid_spacing_m .ge. 5500.)then
                    plot_parms%iraster = 1
                    chigh = 40.
                    clow = -40.
                else
                    plot_parms%iraster = 1
                    chigh = 400.
                    clow = -400.
                endif

                call plot_field_2d(i4_valid,c_type_i,w_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'omega')       

                call move(w_2d,field_2d,nx_l,ny_l)

            elseif(c_field .eq. 'dv')then ! display divergence field
                call divergence(u_2d,v_2d,field_2d,lat,lon,nx_l,ny_l
     1                         ,.true.,r_missing_data)

                c19_label = ' dvrg (cptd) 1e-5/s'

                if(c_type_i(1:2) .eq. 'wf')then
                    c19_label = ' div  (diff) 1e-5/s'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(c_type_i(1:2) .eq. 'wb')then
                    c19_label = ' div  ('//ext(1:3)//')  1e-5/s'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(c_type_i(1:2) .eq. 'wd')then
                    c19_label = ' div  ('//ext(1:3)//')  1e-5/s'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(c_type_i(1:2) .eq. 'wr')then
                    call mk_fcst_hlabel(k_mb,'divergence',fcst_hhmm       
     1                                 ,ext(1:3),'1e-5/s'
     1                                 ,c_model,c_label)
                elseif(c_type_i(1:2) .eq. 'bw')then
                    c30_label = ' divergence (balanced)  1e-5/s'
                    call mklabel(k_mb,c30_label,c_label)
                else
                    c19_label = ' div  (anal) 1e-5/s'
                    call mklabel(k_mb,c19_label,c_label)
                endif


                scale = 1e-5
!               call contour_settings(field_2d,nx_l,ny_l,clow,chigh,cint       
!    1                               ,zoom,density,scale)

                clow = -200.
                chigh = +200.
                cint = 20.

                call plot_field_2d(i4time_3dw,c_type_i,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

            elseif(c_field .eq. 'vo')then ! display vorticity field
                call vorticity_abs(u_2d,v_2d,field_2d,lat,lon,nx_l,ny_l       
     1                            ,dx,dy,.true.,r_missing_data)

                write(6,*)' u range ',minval(u_2d),maxval(u_2d)    
                write(6,*)' v range ',minval(v_2d),maxval(v_2d)    
                write(6,*)' dx range ',minval(dx),maxval(dx)    
                write(6,*)' dy range ',minval(dy),maxval(dy)    
                write(6,*)' vorticity range ',minval(field_2d)    
     1                                       ,maxval(field_2d)

                if(c_type_i(1:2) .eq. 'wf')then
                    c19_label = ' vort (diff) 1e-5/s'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(c_type_i(1:2) .eq. 'wb')then
                    c19_label = ' vort (lga)  1e-5/s'
                    call mklabel(k_mb,c19_label,c_label)
                elseif(c_type_i(1:2) .eq. 'wr')then
                    c19_label = ' vort (fua)  1e-5/s'

!                   note that c_model is blank in this case
                    call mk_fcst_hlabel(k_mb,'vort',fcst_hhmm
     1                                 ,ext(1:3),'1e-5/s'
     1                                 ,c_model,c_label)

                elseif(c_type_i(1:2) .eq. 'bw')then
                    c19_label = ' vort (bal)  1e-5/s'
                    call mklabel(k_mb,c19_label,c_label)
                else
                    c19_label = ' vort (anal) 1e-5/s'
                    call mklabel(k_mb,c19_label,c_label)
                endif

                scale = 1e-5
             
                clow = -20.
                chigh = +80.
                cint = 5.

                call plot_field_2d(i4time_3dw,c_type_i,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'spectral')

            elseif(c_field .eq. 'pv')then ! display potential vorticity field

!               this should also read in the 3d temperature field
                u_2d1(:,:,1) = u_2d(:,:)
                v_2d1(:,:,1) = v_2d(:,:)
                call calc_potvort(i4time,u_2d1,v_2d1,temp_3d,field_3d       
     1                           ,lat,lon,nx_l,ny_l,nz_l,1,k_level
     1                           ,.true.,dx,dy,r_missing_data,istatus)       

                field_2d(:,:) = field_3d(:,:,1)

                if(c_type_i(1:2) .eq. 'wf')then
                    c19_label = ' pvort (diff)  pvu '
                    call mklabel(k_mb,c19_label,c_label)
                elseif(c_type_i(1:2) .eq. 'wb')then
                    c19_label = ' pvort (lga)   pvu '
                    call mklabel(k_mb,c19_label,c_label)
                elseif(c_type_i(1:2) .eq. 'wr')then
                    c19_label = ' pvort (fua)   pvu '

!                   note that c_model is blank in this case
                    call mk_fcst_hlabel(k_mb,'pvort',fcst_hhmm
     1                                 ,ext(1:3),'pvu'
     1                                 ,c_model,c_label)

                elseif(c_type_i(1:2) .eq. 'bw')then
                    c19_label = ' pvort (bal)   pvu '
                    call mklabel(k_mb,c19_label,c_label)
                else
                    c19_label = ' pvort (anal)  pvu '
                    call mklabel(k_mb,c19_label,c_label)
                endif

                scale = 1e-6
                
                clow = -5.
                chigh = +5.
                cint = 0.5

                call plot_field_2d(i4time_3dw,c_type_i,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'spectral')

            elseif(c_field .eq. 'va')then ! vorticity advection field
                call vorticity_abs(u_2d,v_2d,field_2d,lat,lon,nx_l,ny_l       
     1                            ,dx,dy,.true.,r_missing_data)
                call cpt_advection(field_2d,u_2d,v_2d,dx,dy,nx_l,ny_l
     1                            ,field2_2d)
                call mklabel(k_mb,' vort adv 1e-9/s**2',c_label)

                scale = 1e-9
                call contour_settings(field2_2d,nx_l,ny_l
     1                           ,clow,chigh,cint,zoom,density,scale)      

                call plot_cont(field2_2d,scale,clow,chigh,cint,
     1               asc9_tim,namelist_parms,plot_parms,
     1               c_label,i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)

            endif

        elseif(c_type .eq. 'lw')then ! read in liw field from 3d grids

            write(6,*)
            write(6,*)'    looking for laps li*w data:'

            ext = 'liw'
            call get_directory(ext,directory,len_dir)
            c_filespec = directory(1:len_dir)//'*.'//ext(1:3)
            call get_file_time(c_filespec,i4time_ref,i4time_3dw)
            call make_fnam_lp(i4time_3dw,asc9_tim,istatus)

            write(6,*)' getting pregenerated li * omega file'

            var_2d = 'liw'
            call get_laps_2d(i4time_3dw,ext,var_2d
     1          ,units_2d,comment_2d,nx_l,ny_l,liw,istatus) ! k-pa/s


            chigh = 50.

            call plot_cont(liw,1e0,0.0,50.0,-0.5,asc9_tim,
     1              namelist_parms,plot_parms,
     1              'sfc li x 600mb omega  pa-k/s     ',i_overlay
     1              ,c_display,lat,lon,jdot
     1              ,nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif(c_type_i .eq. 'um')then 
            write(6,*)
            write(6,*)'    looking for upslope moisture flux:'

            ext = 'liw'
            call get_directory(ext,directory,len_dir)
            c_filespec = directory(1:len_dir)//'*.'//ext(1:3)
            call get_file_time(c_filespec,i4time_ref,i4time_3dw)
            call make_fnam_lp(i4time_3dw,asc9_tim,istatus)

            write(6,*)' getting pregenerated li * omega file'

            var_2d = 'umf'
            call get_laps_2d(i4time_3dw,ext,var_2d
     1          ,units_2d,comment_2d,nx_l,ny_l,field_2d,istatus) ! k-pa/s

            plot_parms%iraster = 1

            c_label = 'upslope moisture flux (cm-m/s)'
            call plot_field_2d(i4time_3dw,c_type,field_2d,.01
     1                        ,namelist_parms,plot_parms
     1                        ,umf_l,umf_h,10.,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'upflux')

        elseif(c_type_i .eq. 'li')then ! read in li 'field_2d' from 3d grids
            if(lapsplot_pregen)then
                write(6,*)' getting li from lst'
!               read in li data
                var_2d = 'li'
                ext = 'lst'
                call get_laps_2dgrid(i4time_ref,lagt,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                          ,field_2d,0,istatus)

            else
                call get_laps_2dgrid(i4time_ref,lagt,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                          ,field_2d,0,istatus)

            endif

            if(istatus .ne. 1)then
                write(6,*)' error reading lifted index data'
                goto1200
            endif

            call make_fnam_lp(i4time_nearest,asc9_tim_n,istatus)

            c_label = '    sfc lifted index     (k)     '

            call plot_field_2d(i4time_nearest,c_type,field_2d,1.0
     1                        ,namelist_parms,plot_parms
     1                        ,+10.,-10.,2.,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'spectral')


        elseif(c_type_i .eq. 's')then ! read in lst data generically
            ext = 'lst'

            write(6,825)
 825        format(/'  select field (var_2d):  '
     1       /
     1       /'     lst (stability): [li,pbe,nbe,si,tt,k,lcl,wb0] ? ',$)       

            read(lun,824)var_2d
 824        format(a)
            call upcase(var_2d,var_2d)

            level=0
            call get_laps_2dgrid(i4time_ref,lagt
     1                          ,i4time_nearest
     1                          ,ext,var_2d,units_2d,comment_2d
     1                          ,nx_l,ny_l
     1                          ,field_2d,0,istatus)

            if(istatus .ne. 1 .and. istatus .ne. -1)then
                write(6,*)' error reading stability analysis ',var_2d
                goto1200
            endif

            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

            call s_len(comment_2d,len_comment)
            call s_len(units_2d,len_units)

            if(len_units .gt. 0)then
                if(units_2d(1:len_units) .eq. 'm')then
                    c_label = comment_2d(1:len_comment)
     1                   //'   ('//units_2d(1:len_units)//'-msl)    '
                else
                    c_label = comment_2d(1:len_comment)
     1                      //'   ('//units_2d(1:len_units)//')    '
                endif
            else
                c_label = comment_2d(1:len_comment)
            endif

            scale = 1.
            call contour_settings(field_2d,nx_l,ny_l,clow,chigh,cint
     1                           ,zoom,density,scale)       

!           call plot_cont(field_2d,scale,clow,chigh,cint
!    1                    ,asc9_tim,namelist_parms,plot_parms
!    1                    ,c_label,i_overlay,c_display,lat,lon,jdot
!    1                    ,nx_l,ny_l,r_missing_data,laps_cycle_time)

            call plot_field_2d(i4time_nearest,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

        elseif(c_type .eq. 'tw')then
	    if (laps_cycle_time.eq.0)then
	        i4time_temp=i4time_ref
	    else
                i4time_temp = i4time_ref / laps_cycle_time 
     1                                   * laps_cycle_time
            endif
!           read in surface temp data
            var_2d = 't'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time,i4time_temp       
     1                          ,ext,var_2d,units_2d,comment_2d
     1                          ,nx_l,ny_l,temp_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' laps sfc temp not available'
                return
            endif

!           read in surface dewpoint data
            var_2d = 'td'
            ext = 'lsx'
            call get_laps_2d(i4time_temp,
     1                       ext,var_2d,units_2d,comment_2d,
     1                       nx_l,ny_l,td_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' laps sfc dewpoint not available'
                return
            endif

!           read in surface pressure data
            var_2d = 'ps'
            ext = 'lsx'
            call get_laps_2d(i4time_temp,
     1                       ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                       ,pres_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' laps sfc pressure not available'
                return
            endif

            call get_tw_approx_2d(temp_2d,td_2d,pres_2d,nx_l,ny_l
     1                           ,tw_sfc_k)

            call make_fnam_lp(i4time_temp,asc9_tim_n,istatus)

            zero_c = 273.15

            do i = 1,nx_l
            do j = 1,ny_l
                field_2d(i,j) = tw_sfc_k(i,j) - zero_c
            enddo ! j
            enddo ! i


            call plot_cont(field_2d,1e-0,-30.,+30.,2.,asc9_tim_n,
     1                    namelist_parms,plot_parms,
     1                    '    sfc wetbulb (approx) (c)     '
     1                    ,i_overlay,c_display,lat,lon,jdot,
     1                     nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'ms' .or. c_type .eq. 'ob'
     1                          .or. c_type .eq. 'st' ! station locations  
     1                          .or. c_type .eq. 'of' ! air t,td in f
     1                          .or. c_type .eq. 'oc' ! air t,td in c
     1                          .or. c_type .eq. 'os' ! stations  
     1                          .or. c_type .eq. 'ov' ! sky cover, visibility
     1                          .or. c_type(1:2) .eq. 'op' ! precip 
     1                          .or. c_type .eq. 'og' ! soil/water t + solar rad
     1                          .or. c_type .eq. 'or' ! solar radiation only
     1                          .or. c_type .eq. 'oh' ! humidity (gps-pw)
     1                          .or. c_type .eq. 'ow' ! wind only
     1                          .or. c_type .eq. 'qf' ! qc air t,td in f
     1                          .or. c_type .eq. 'qc' ! qc air t,td in c
     1                          .or. c_type .eq. 'qs' ! qc stations
     1                          .or. c_type .eq. 'qv' ! qc sky cover, visib
     1                          .or. c_type .eq. 'qp' ! qc precip
     1                          .or. c_type .eq. 'qg' ! qc soil/water t
     1                          .or. c_type .eq. 'mw' ! mesowx points
     1                                                )then

            i4time_plot = i4time_ref ! / laps_cycle_time * laps_cycle_time
            call get_filespec('lso',2,c_filespec,istatus)
            call get_file_time(c_filespec,i4time_ref,i4time_plot)
            call make_fnam_lp(i4time_plot,a9time,istatus)

            if(c_type(2:2) .eq. 's' )iflag = 0 ! station locations
            if(c_type      .eq. 'ms')iflag = 1
            if(c_type(2:2) .ne. 's' )iflag = 2
            if(c_type(1:2) .eq. 'mw')iflag = 0 ! mesowx points locations
            if(c_type(1:2) .eq. 'st')iflag = 0 ! station locations

            c_label = '                                 '

            igrid = 0

            call plot_stations(a9time,c_label,c_type,i_overlay
     1                        ,namelist_parms,plot_parms
     1                        ,max_snd_grid,max_snd_levels
     1                        ,c_display,lat,lon,topo,c_file,iflag
     1                        ,nx_l,ny_l,nz_l,laps_cycle_time,zoom)

        elseif(c_type(1:2) .eq. 'he')then
            write(6,*)
!           write(6,*)'    looking for 3d laps wind data:'
!           call get_file_time(c_filespec,i4time_ref,i4time_3dw)

            write(6,*)' getting pregenerated helicity file'
            var_2d = 'lhe'
            ext = 'lhe'

            call get_laps_2dgrid(i4time_ref,lagt,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                  ,field_2d,0,istatus)

            i4time_3dw = i4time_nearest
            call make_fnam_lp(i4time_3dw,asc9_tim,istatus)
            
            abs_max = 0
            do i = 1,nx_l
            do j = 1,ny_l
                abs_max = max(abs_max,abs(field_2d(i,j)))
            enddo ! j
            enddo ! i

            write(6,*)' max helicity magnitude = ',abs_max

            c_field = 'he'
            kwind = 0
            clow = hel_l
            chigh = hel_h

            if(abs_max .gt. 1.)then ! new way
                scale = 1.
                c_label = 'storm rel helicity m**2/s**2     '           
                cint = 40.
            else                    ! old way
                scale = 1e-4
                c_label = 'helicity  sfc-500 1e-4m/s**2     '          
                cint = 5.
            endif

            scale = 1.
            call plot_field_2d(i4time_3dw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'spectral')

        elseif(c_type(1:2) .eq. 'bl' .or. c_type(1:2) .eq. 'lf')then       
            write(6,*)

            write(6,826)
 826        format(/'  select field (var_2d):  '
     1       /
     1       /' pbl | firewx: [ptp,pdm | vnt,ham,hah,fwi,cwi] ? ',$)       

            read(lun,824)var_2d
!824        format(a)
            call upcase(var_2d,var_2d)

            if(c_type(1:2) .eq. 'bl')then
                ext = 'pbl'
            else
                ext = 'lfr'
            endif

            write(6,*)' getting pregenerated file ',ext(1:3),' ',var_2d       

            call get_laps_2dgrid(i4time_ref,lagt,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                  ,field_2d,0,istatus)
            if(istatus .eq. 0)then
                write(6,*)' cant find ',var_2d,' analysis ',istatus
                goto1200
            else
                write(6,*)' comment is ',trim(comment_2d)
            endif

            i4time_3dw = i4time_nearest
            call make_fnam_lp(i4time_3dw,asc9_tim,istatus)

            if(units_2d(1:2) .eq. 'pa')then
                units_2d = 'hpa'
                scale = 100.
            elseif(units_2d(1:4) .eq. 'none')then
                units_2d = '          '
                scale = 1.
            else ! default value
                scale = 1.
            endif

            colortable = 'spectral'
            
            if(var_2d(1:3) .eq. 'ptp')then
                clow = 500.
                chigh = 1100.
                cint = 0.
            elseif(var_2d(1:3) .eq. 'pdm')then
                if(namelist_parms%c_pbl_depth_units .eq. 'english')then
                    clow = 0.
                    chigh = 8000.
                    cint = 1000.
                    units_2d = 'ft'
                    scale = 1. / ft_per_m
                else ! 'metric'
                    clow = 0.
                    chigh = 2400.
                    cint = 400.
                    units_2d = 'm'
                    scale = 1. 
                endif
            elseif(var_2d(1:3) .eq. 'vnt')then
                if(c_vnt_units .eq. 'kt-ft')then
                    clow = 150.
                    scale = 1000. / (ft_per_m / mspkt) ! convert from m**2/s 
                                                       ! to ft-kt (inverse)
                    units_2d = 'kt-ft x1000'
                else
                    clow = 5000.
                    units_2d = c_vnt_units
                endif

                chigh = 0.
                cint = 0.

                colortable = 'vnt'

                if(i_image .eq. 1)then
                    plot_parms%icol_barbs = +1 ! keep future barbs plots bright
                endif

            elseif(var_2d(1:2) .eq. 'ha')then
                clow = 2.
                chigh = 6.
                cint = 1.0

                colortable = 'haines'

            elseif(var_2d(1:3) .eq. 'fwi')then
                clow = 0.
                chigh = 40.
                cint = 0.

            elseif(var_2d(1:3) .eq. 'cwi')then
                clow =  0.
                chigh = 1.
                cint = 1.

                colortable = 'cwi'

            else ! default values
                clow = 0.
                chigh = 0.
                cint = 0.
            endif

            call s_len2(units_2d,len_units)
            if(len_units .gt. 0)then
                c_label = comment_2d(1:26)//' ('//units_2d(1:len_units)
     1                                          //')'       
            else
                c_label = comment_2d(1:26)
            endif

            call plot_field_2d(i4time_3dw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,colortable)

        elseif(c_type .eq. 'lv'.or.c_type .eq. 'lr' )then

         if(c_type .eq. 'lv' )then

          ext = 'lvd'

          call config_satellite_lvd(istatus)
          if(istatus.ne.1)then
             write(6,*)'warning config_satellite_lvd status = ',istatus
!            return
          endif


          write(6,*)' available satellites are:'
          do k=1,maxsat
           if(isats(k).eq.1)then
            write(6,*)c_sat_id(k)
           endif
          enddo

          write(6,105)
 105      format(/' enter one satellite you would like to plot.'
     1       /' to choose among multiple satellites enter "multi" : '      
     1       ,$)       

          read(lun,*)c_sat_plot
          call filter_string(c_sat_plot)

          write(6,*)' entered satellite is ',c_sat_plot

          isl = 1
          ish = maxsat

!         determine which satellite was desired
          if(c_sat_plot .ne. 'multi')then 
              do k=1,maxsat
                  if(isats(k).eq.1)then
                      if(c_sat_id(k) .eq. c_sat_plot)then
                          write(6,*)' you have selected satellite # ',k       
                          isl = k
                          ish = k
                      endif
                  endif
              enddo
          endif

          write(6,*)

          do k=isl,ish
           if(isats(k).eq.1)then

            if(isl .ne. ish)then
                write(6,114)c_sat_id(k)
114             format(5x,'plot the data for ',a6,34x,' [y/n] ? ',$)
                read(lun,*)cansw
            else
                cansw = 'y'
            endif

            if(cansw.eq.'y'.or.cansw.eq.'y')then

             call get_directory(ext,directory,len_dir)
             directory=directory(1:len_dir)//trim(c_sat_id(k))//'/'
c
c determine which channels have been processed for this satellite
c
             j=0
             lfndtyp=.false.
             do while(.not.lfndtyp.and.j.le.maxtype)
              j=j+1
              if(itypes(j,k).eq.1)then
               ist=j
               lfndtyp=.true.
               write(6,*)' found type for ',j,k
              endif
             enddo

             write(6,118)
118          format(5x,'select field',1x,'(vis, 3.9, 6.7, 11.2, 12, 13)'
     1              /32x,' [enter 1, 2, 3, 4, 5, 6; neg for img] ? ',$)
             read(lun,*)ilvd

             if(ilvd .lt. 0)then
                 l_plot_image = .true.
                 ilvd = -ilvd
             else
                 l_plot_image = .false.
             endif
             
             write(6,*)
             write(6,*)'    looking for laps lvd data:'
c
             if(ichannels(ilvd,ist,k).eq.1)then
              write(6,121)clvdvars(ilvd)
121           format(5x,'select 2d var name:',3x,a15,10x, $)
              read(lun,*)var_2d
              call upcase(var_2d,var_2d)
             else
              print*,'this channel was not processed ',ilvd,ist,k
              print*,'ichannels: ',ichannels
              goto 119
             endif

             call get_2dgrid_dname(directory
     1               ,i4time_ref,100000,i4time_nearest,ext,var_2d
     1               ,units_2d,comment_2d,nx_l,ny_l,vas,0,istatus)

             if(istatus .eq. 0)then
              write(6,*)' cant find ',var_2d,' analysis ',istatus
              goto1200
             endif

             colortable = 'linear'
             if(ilvd.gt.1)then
              c_label='b-temps (c): '//c_sat_id(k)//'/'//var_2d
              vasmx=-255.
              vasmn=255.
              do i = 1,nx_l
              do j = 1,ny_l
                 if(vas(i,j).ne.r_missing_data)then
                    vas(i,j) = vas(i,j) - 273.15
                    vasmx=int(max(vas(i,j),vasmx))
                    vasmn=int(min(vas(i,j),vasmn))
                 endif
              enddo
              enddo
              clow = -80.
              chigh = +40.
              cint = (vasmx-vasmn)/10.
              scale = 1e0
              scale_l = btemp_h       ! for image plots
              scale_h = btemp_l       ! for image plots
              colortable = btemp_table
             elseif(var_2d.eq.'alb')then
              c_label='albedo '//c_sat_id(k)
              scale_l = 0.00          ! for image plots
              scale_h = cloud_albedo_a! for image plots
             elseif(var_2d.eq.'svs')then
              c_label='vis / svs reflectance -   '//c_sat_id(k)
              scale_l = 0.0           ! for image plots
              scale_h = 1.2           ! for image plots
             else
              c_label='vis counts (normalized) - '//c_sat_id(k)
              scale_l = 20.           ! for image plots
              scale_h = 110.          ! for image plots
             endif

             call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

             if(l_plot_image)then
                 call ccpfil(vas,nx_l,ny_l,scale_l,scale_h,colortable
     1                ,n_image,scale,'hsect',plot_parms,namelist_parms)
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                 call setusv_dum('in',7)
                 call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                 call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                              ,namelist_parms,plot_parms)

             else ! contours
                 if(ilvd .eq. 1)then ! visible data
                     if(var_2d.eq.'alb' .or. var_2d .eq. 'svs')then
                         clow = 0.0
                         chigh = 1.
                         cint = 0.1
                         scale = 1e0
                     else
                         clow = 0.0
                         chigh = 256.
                         cint = 05.
                         scale = 1e0
                     endif
                 else                ! ir data
                     call contour_settings(vas,nx_l,ny_l
     1                                   ,clow,chigh,cint
     1                                   ,zoom,density,scale)
                 endif

                 call plot_cont(vas,scale,clow,chigh,cint,asc9_tim,
     1             namelist_parms,plot_parms,c_label,
     1             i_overlay,c_display,lat,lon,jdot,
     1             nx_l,ny_l,r_missing_data,laps_cycle_time)

             endif
             write(6,*)' copying vas array to field_2d to support diff'       
             call move(vas,field_2d,nx_l,ny_l) ! supports the diff option

            endif !(cansw)
           endif  !(isats)
119       enddo   !(maxsat)

         else     !(c_type='lr'?)

          print*,' lsr plotting currently not available'

         endif    !(c_type='lv'?)

        elseif(c_type .eq. 'v3')then

         var_2d = 'alb'
         call get_2dgrid_dname(directory
     1        ,i4time_ref,100000,i4time_nearest,ext,var_2d
     1        ,units_2d,comment_2d,nx_l,ny_l,vas,0,istatus)

         if(istatus .eq. 0)then
            write(6,*)' cant find alb analysis'
            goto1200
         endif
         c_label = 'vis cld frac (tenths) '//c_sat_id(k)
         clow  =  -6.
         chigh = +16.
         cint = 2.0
         scale = 1e-1

         do i = 1,nx_l
         do j = 1,ny_l
          if(vas(i,j) .ne. r_missing_data)then
             vas(i,j) = albedo_to_cloudfrac(vas(i,j))
          endif
         enddo
         enddo
         call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

         call plot_cont(vas,scale,clow,chigh,cint,asc9_tim,
     1        namelist_parms,plot_parms,c_label,i_overlay,
     1        c_display,lat,lon,jdot,       
     1        nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'v5')then

         var_2d = 's8a'
         call get_2dgrid_dname(directory
     1       ,i4time_ref,100000,i4time_nearest,ext,var_2d
     1       ,units_2d,comment_2d,nx_l,ny_l,vas,0,istatus)
         if(istatus .eq. 0)then
            write(6,*)' cant find vas/s8a analysis'
            goto1200
         endif
         c_label = 'sfc t - band 8  (k)  -'//c_sat_id(k)
         clow = -8.0
         chigh = 20.
         cint = 4.
         scale = 1e0

!  get sfc t to take the difference...
         ext = 'lsx'
         var_2d = 't'
         call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1       ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l,dum1_array,0
     1       ,istatus)
         if(istatus .ne. 1)then
            write(6,*)' cant find vas/s8a analysis'
            goto1200
         endif

         do i = 1,nx_l
         do j = 1,ny_l
           vas(i,j) = dum1_array(i,j) - vas(i,j)
         enddo
         enddo

         call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

         call plot_cont(vas,scale,clow,chigh,cint,asc9_tim,
     1       namelist_parms,plot_parms,c_label,i_overlay,c_display,
     1       lat,lon,jdot,nx_l,ny_l,r_missing_data,laps_cycle_time)


        elseif( c_type .eq. 'po' )then

          call make_fnam_lp(i4time_ref,asc9_tim,istatus)
          ext = 'lsr'
          call get_directory(ext,directory,len_dir)
          cfname=directory(1:len_dir)//asc9_tim//'_12.lsr'
          open(14,file=cfname,form='unformatted',status='old',err=18)
          goto 27
18        cfname=directory(1:len_dir)//asc9_tim//'_14.lsr'
          open(14,file=cfname,form='unformatted',status='old',err=19)

          n=index(cfname,' ')-1
          write(6,*)'reading ',cfname(1:n)
          goto 28
27        write(6,*)'reading ',cfname(1:n)
28        read(14)sndr_po
          write(6,*)'enter the channel [1-19]'
          read(5,25)ichan
25        format(i2)
          do j=1,ny_l
          do i=1,nx_l
             if(sndr_po(ichan,i,j).ne.r_missing_data)then
                vas(i,j)=sndr_po(ichan,i,j)-273.15
             endif
          enddo
          enddo

          clow = -80.
          chigh = +40.
          cint = 5.
          write(cchan,111)ichan
 111      format(i2)
          if(cchan(1:1).eq.' ')cchan(1:1)='0'
          if(cchan(2:2).eq.' ')cchan(2:2)='0'

          c_label = 'polar orbiter channel '//cchan//' deg c'

          call plot_cont(vas,1e0,clow,chigh,cint,asc9_tim,
     1       namelist_parms,plot_parms,c_label,i_overlay,c_display,
     1       lat,lon,jdot,nx_l,ny_l,r_missing_data,laps_cycle_time)

          goto 21
19        write(6,*)'not able to open an lsr file ', asc9_tim
21        continue

        elseif(c_type(1:2) .eq. 'lc')then ! satellite/cloud data in lcv
            write(6,122)
 122        format(5x,'select field name [s8a, s3a, d39, alb] :    ',$)       
            read(lun,*)var_2d_in
            call upcase(var_2d_in,var_2d_in)

            ext = 'lcv'

            if(var_2d_in(1:3) .ne. 'd39')then
                var_2d = var_2d_in
                call get_laps_2dgrid(i4time_ref,3*laps_cycle_time
     1                              ,i4time_nearest,ext,var_2d
     1                              ,units_2d,comment_2d,nx_l,ny_l
     1                              ,vas,0,istatus)
                if(istatus .eq. 0)then
                    write(6,*)' cant find ',var_2d,' analysis ',istatus
                    goto1200
                endif

                colortable = 'linear'

            else ! 'd39' field which is the difference of 's3a' - 's8a'
                var_2d = 's8a'
                call get_laps_2dgrid(i4time_ref,3*laps_cycle_time
     1                              ,i4time_nearest,ext,var_2d
     1                              ,units_2d,comment_2d,nx_l,ny_l
     1                              ,vas,0,istatus)
                if(istatus .eq. 0)then
                    write(6,*)' cant find ',var_2d,' analysis ',istatus
                    goto1200
                endif

                var_2d = 's3a'
                call get_laps_2dgrid(i4time_nearest,0
     1                              ,i4time_nearest,ext,var_2d
     1                              ,units_2d,comment_2d,nx_l,ny_l
     1                              ,field_2d,0,istatus)
                if(istatus .eq. 0)then
                    write(6,*)' cant find ',var_2d,' analysis ',istatus
                    goto1200
                endif

!               subtract the two satellite fields
                call diff(field_2d,vas,vas,nx_l,ny_l)

                colortable = 'hues'

            endif

            c_label='lcv '//comment_2d(1:24)

            if(var_2d_in .eq. 's8a' .or. var_2d_in .eq. 's3a')then
              c_label(30:38) = '    deg c'
              do i = 1,nx_l
              do j = 1,ny_l
                 if(vas(i,j).ne.r_missing_data)then
                    vas(i,j) = vas(i,j) - 273.15
                 endif
              enddo
              enddo

              scale = 1e0
              scale_l = +40.          ! for image plots
              scale_h = btemp_l       ! for image plots
            elseif(var_2d_in.eq.'alb')then
!             c_label='albedo '//c_sat_id(k)
              scale_l = 0.00          ! for image plots
              scale_h = cloud_albedo_a! for image plots
            elseif(var_2d_in.eq.'cla')then
              scale_l = 0.00          ! for image plots
              scale_h = cloud_albedo_a! for image plots
            elseif(var_2d_in.eq.'cwt')then
              scale_l = 0.00          ! for image plots
              scale_h = cloud_albedo_a! for image plots
            elseif(var_2d_in.eq.'svs')then
              c_label='vis reflectance - '//c_sat_id(k)
              scale_l = 30.           ! for image plots
              scale_h = 100.          ! for image plots
            elseif(var_2d_in.eq.'d39')then
              c_label='3.9u - 11u difference  (k)'
              scale_l = -10           ! for image plots
              scale_h = +10.          ! for image plots
            elseif(var_2d_in.eq.'swi')then
              c_label='downward short wave (w/m**2)'              
              scale_l = 0.            ! for image plots
              scale_h = 1000.         ! for image plots
              colortable = 'spectral'
            else
!             c_label='vis counts (normalized) - '//c_sat_id(k)
              scale_l = 30.           ! for image plots
              scale_h = 230.          ! for image plots
            endif

            scale = 1.0 

            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

            if(c_type(3:3) .eq. 'i')then
                call ccpfil(vas,nx_l,ny_l,scale_l,scale_h,colortable
     1                     ,n_image,scale,'hsect',plot_parms
     1                     ,namelist_parms)
                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum('in',7)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)

            else ! contours
                if(var_2d_in.eq.'alb')then
                    clow = 0.0
                    chigh = 1.
                    cint = 0.1
                    scale = 1e0
                else
!                   clow = -10.0
!                   chigh = 40.
!                   cint = 10.
                    scale = 1e0
                    call contour_settings(vas,nx_l,ny_l,clow,chigh,cint
     1                                   ,zoom,density,scale)
                endif

                call plot_cont(vas,scale,clow,chigh,cint,asc9_tim,
     1             namelist_parms,plot_parms,c_label,i_overlay,
     1             c_display,lat,lon,jdot,
     1             nx_l,ny_l,r_missing_data,laps_cycle_time)

             endif
             write(6,*)' copying vas array to field_2d to support diff'       
             call move(vas,field_2d,nx_l,ny_l) ! supports the diff option

        elseif( c_type .eq. 'ra' .or. c_type .eq. 'cl'
     1    .or.  c_type .eq. 'rr' .or. c_type .eq. 'rf'
     1    .or.  c_type .eq. 'rd' .or. c_type .eq. 'rv')then

            ndim_read = 2 ! default

            if(c_type .eq. 'ra')mode = 1
            if(c_type .eq. 'cl')mode = 2
cabdel	    
            if (laps_cycle_time.eq.0)then
	        i4time_tmp1= i4time_file
                i4time_tmp2=i4time_ref-2400
            else
                i4time_tmp1 = (i4time_ref)/laps_cycle_time 
     1                      * laps_cycle_time
                i4time_tmp2 = (i4time_ref-2400)/laps_cycle_time 
     1                      * laps_cycle_time
        endif
cabdel	
            if(c_type .eq. 'rr')then
                if(i4time_ref .ne. i4time_tmp1)then
                    i4time_get = i4time_tmp2
                else
                    i4time_get = i4time_ref
                endif
            else
                i4time_get = i4time_ref
            endif

            if(c_type .ne. 'rf')then ! radar intermediate data files

!             if(.not. l_radar_read)then
              if(.true.)then

                if(c_type .eq. 'ra')then ! read data from vrc files

                  if(.false.)then
!                     obtain height field
                      ext = 'lt1'
                      var_2d = 'ht'
                      call get_laps_3dgrid(
     1                   i4time_get,10000000,i4time_ht,
     1                   nx_l,ny_l,nz_l,ext,var_2d
     1                  ,units_2d,comment_2d,field_3d,istatus)
                      if(istatus .ne. 1)then
                        write(6,*)' error locating height field'
                        return
                      endif

                      call get_radar_ref(i4time_get,100000,i4time_radar
     1                ,mode
     1                ,.true.,nx_l,ny_l,nz_l,lat,lon,topo,.true.,.true.       
     1                ,field_3d
     1                ,grid_ra_ref,n_ref
     1                ,rlat_radar,rlon_radar,rheight_radar,istat_2dref
     1                ,istat_3dref)
                  endif

                  ext = 'vrc'
                  call get_filespec(ext,2,c_filespec,istatus)
                  call get_file_time(c_filespec,i4time_get,i4time_radar)       

                  write(6,*)' calling read_radar_2dref...'
                  call read_radar_2dref(i4time_radar,radar_name
     1                                 ,nx_l,ny_l,field_2d,istat_2dref)       

                  write(6,*)' istat_2dref = ',istat_2dref

                  call get_ref_base(ref_base,istatus)

                  write(6,*)' leaving missing data points alone'
                  i_miss = 0
                  do i = 1,nx_l
                  do j = 1,ny_l
                      if(field_2d(i,j) .eq. r_missing_data)then
!                         field_2d(i,j) = ref_base
                          i_miss = i_miss + 1
                      endif
                  enddo ! j
                  enddo ! i

                  pct_miss = float(i_miss) / float(nx_l*ny_l) * 100.

                  write(6,*)' percent coverage = ', 100. - pct_miss

                elseif(c_type .eq. 'rd')then ! read data from v01, v02, etc.

                  write(6,*)' reading velocity data from the radars'

                  if(ialloc_vel .eq. 0)then ! allocate 4d radar arrays
                      allocate(grid_ra_vel(nx_l,ny_l,nz_l,max_radars)       
     1                        ,stat=istat_alloc)
                      if(istat_alloc .ne. 0)then
                          write(6,*)
     1                        ' error: could not allocate grid_ra_vel'
                          stop
                      endif

                      allocate(grid_ra_nyq(nx_l,ny_l,nz_l,max_radars)       
     1                        ,stat=istat_alloc)
                      if(istat_alloc .ne. 0)then
                          write(6,*)
     1                        ' error: could not allocate grid_ra_nyq'
                          stop
                      endif

                      ialloc_vel = 1

                  endif ! ialloc_vel

                  nx_r = nx_l
                  ny_r = ny_l

                  call get_multiradar_vel(
     1                i4time_get,100000000,i4time_radar_a
     1               ,max_radars,n_radars,ext_radar_a,r_missing_data   
     1               ,nx_l,ny_l,nz_l,lat,lon                            ! i
     1               ,nx_r,ny_r,igrid_r                                 ! i
     1               ,grid_ra_vel,grid_ra_nyq,idx_radar,v_nyquist_in_a       
     1               ,ioffset,joffset                                   ! o
     1               ,l_offset                                          ! i
     1               ,n_vel_grids_a
     1               ,rlat_radar_a,rlon_radar_a,rheight_radar_a
     1               ,radar_name_a      
     1               ,istat_radar_vel,istat_radar_nyq)

                  if(istat_radar_vel .eq. 1)then
                    write(6,*)' radar 3d vel data successfully read in'       
     1                       ,(n_vel_grids_a(i),i=1,n_radars)
                  else
                    write(6,*)' radar 3d vel data not successfully read'     
     1                       ,(n_vel_grids_a(i),i=1,n_radars)
                    write(6,*)' istat_radar_vel = ',istat_radar_vel
                    return
                  endif

                  i4time_radar = i4time_radar_a(1)

                elseif(c_type .eq. 'rv')then
!                   ask which radar number (extension)
                    write(6,*)
                    write(6,2026)
2026                format('  enter radar extension (for reflectivity)'      
     1                                                     ,28x,'? ',$)       
                    read(lun,*)ext_radar

                    write(6,*)' reading reflectivity data from radar '
     1                                          ,ext_radar

!                   obtain height field
                    ext = 'lt1'
                    var_2d = 'ht'
                    call get_laps_3dgrid(
     1                   i4time_radar,1000000,i4time_ht,
     1                   nx_l,ny_l,nz_l,ext,var_2d
     1                  ,units_2d,comment_2d,field_3d,istatus)
                    if(istatus .ne. 1)then
                        write(6,*)' error locating height field'
                        call get_pres_3d(i4time_radar
     1                                ,nx_l,ny_l,nz_l,pres_3d,istatus)
                        if(istatus .ne. 1)then
                            write(6,*)' error getting pressure field'      
                            goto 1200
                        else
                            write(6,*)' convert pres to ht'
     1                              ,' - using standard atmosphere'     
                            do k = 1,nz_l
                            do j = 1,ny_l
                            do i = 1,nx_l
                                field_3d(i,j,k) = 
     1                                  psatoz(pres_3d(i,j,k)*.01) 
                            enddo ! i
                            enddo ! j
                            enddo ! k
                        endif
                    endif ! height field may be necessary

                    write(6,*)

                    call get_ref_base(ref_base,istatus)

                    call get_filespec(ext_radar,2,c_filespec,istatus)
                    call get_file_time(c_filespec,i4time_get
     1                                ,i4time_radar)    

                    if(ext_radar .ne. 'vrz')then
                        l_low_fill = namelist_parms%l_low_fill
                        l_high_fill = namelist_parms%l_high_fill
                    else
                        l_low_fill = .false.                        
                        l_high_fill = .false.                        
                    endif

                    i4_tol = 1200

                    call read_radar_3dref(i4time_radar,
     1               i4_tol,i4_ret,                                   ! i/o
!    1               .true.,ref_base,
     1               .true.,r_missing_data,
     1               nx_l,ny_l,nz_l,ext_radar,
     1               lat,lon,topo,l_low_fill,l_high_fill,
     1               field_3d,
     1               grid_ra_ref,
     1               rlat_radar,rlon_radar,rheight_radar,radar_name,
     1               n_ref_grids,istat_radar_2dref,istat_radar_3dref)

!                   i4time_radar = i4time_get

                endif ! c_type

              endif ! l_radar_read

              l_radar_read = .true.

              call make_fnam_lp(i4time_radar,asc9_tim_r,istatus)

            endif

            if(c_type .eq. 'rv')then                     ! vxx reflectivity
                write(6,2021)
2021            format(/'  [rf-i] reflectivity data, '
     1                 /'  [mr-i] column max ref'       
     1                 /' ',61x,' [q] quit ? ',$)

            elseif(c_type .eq. 'rd')then                 ! vxx velocity
                write(6,2022)
2022            format(/'  [ve] velocity contours, '  
     1                 ,' [vi] velocity image '
!    1                 ,'[f1] 1 hr fcst max reflectivity,'
     1                 /' ',61x,' [q] quit ? ',$)

            else                                         ! non vxx reflectivity
                write(6,2023)
2023            format(/'  [rf-i] reflectivity data, '
     1                 /'  [mr-i] column max ref, [mt] max tops,'       
     1                 /'  [lr-i] low lvl refl, '
     1                 ,'[f1] 1 hr fcst max reflectivity,'
     1                 /' ',61x,' [q] quit ? ',$)

            endif

            read(lun,824)c_field

            if(c_field(3:3) .eq. 'i')then
                i_image = 1
            else
                i_image = 0
            endif               

            if(  c_field(1:2) .eq. 'rf' .or. c_field(1:2) .eq. 'rv' 
     1      .or. c_field(1:2) .eq. 'vi' .or. c_field(1:2) .eq. 've')then       
                write(6,2025)
2025            format('         enter level in mb ',45x,'? ',$)
                call input_level(lun,k_level,k_mb,pres_3d
     1                          ,nx_l,ny_l,nz_l)       
            endif

            if(      c_type .eq. 'rf' 
     1         .and. c_field(1:2) .ne. 'mr'
     1         .and. c_field(1:2) .ne. 'lr')then

              call input_product_info(i4time_get            ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

              var_2d = 'ref'

              if(c_prodtype .eq. 'a')then
!               obtain lps reflectivity field
                write(6,*)' reading lps radar volume reflectivity'
                ext = 'lps'
                call get_laps_3dgrid(
     1                   i4time_get,1000000,i4time_radar,
     1                   nx_l,ny_l,nz_l,ext,var_2d
     1                  ,units_2d,comment_2d,grid_ra_ref,istatus)

                call make_fnam_lp(i4time_radar,asc9_tim_r,istatus)

                ndim_read = 3

              elseif(c_prodtype .eq. 'f')then
                call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,k_mb,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 field_2d,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' could not read forecast field'       
                    goto1200
                endif

!               i4time_radar = i4_valid

                call directory_to_cmodel(directory,c_model)

              endif

            endif

            if(c_field(1:2) .eq. 'mr')then ! column max reflectivity data
                plot_parms%icol_barbs = +1 ! keep future barbs plots bright
                if(c_type .eq. 'rf')then
                    write(6,*)' getting analyzed lmr file'
                    var_2d = 'r'
                    ext = 'lmr'
!                   i4time_hour = (i4time_radar+laps_cycle_time/2)
!    1                          /laps_cycle_time * laps_cycle_time
                    call get_laps_2dgrid(i4time_ref,10000,i4time_radar
     1                                  ,ext,var_2d,units_2d,comment_2d       
     1                                  ,nx_l,ny_l,radar_array,0
     1                                  ,istatus)

                    c_label = 'composite reflectivity (analysis) '    

                elseif(c_type .eq. 'ra')then
                    radar_array = field_2d
                    c_label = 'composite ref (interim-vrc) '    

                else
                    write(6,*)' calling get_max_ref'
                    call get_max_reflect(grid_ra_ref,nx_l,ny_l,nz_l
     1                                  ,r_missing_data,radar_array)
                    c_label = 'composite ref (intermediate) '    

                endif

                call make_fnam_lp(i4time_radar,asc9_tim_r,istatus)

                write(6,*)' c_field = ',c_field,' ',c_field(3:3)

                field_2d = radar_array ! diff support

!               display r field
                if(c_field(3:3) .ne. 'i')then
                    call plot_cont(radar_array,1e0,0.,chigh,cint_ref,
     1               asc9_tim_r,namelist_parms,plot_parms,c_label,
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)
                else
                    call ccpfil(radar_array,nx_l,ny_l,-10.,70.,'ref'
     1                         ,n_image,1e0,'hsect',plot_parms
     1                         ,namelist_parms)
                    call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                    call setusv_dum('in',7)
                    call write_label_lplot(nx_l,ny_l,c_label
     1                                    ,asc9_tim_r,plot_parms
     1                                    ,namelist_parms,i_overlay
     1                                    ,'hsect')       
                    call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                                 ,namelist_parms,plot_parms)

                endif

            elseif(c_field(1:2) .eq. 'lr')then ! low lvl reflectivity data
                if(c_type .eq. 'rf')then
                    write(6,*)' getting analyzed lmt file'
                     
c abdel       
                    if (laps_cycle_time.eq.0)then
	                i4time_hour= (i4time_radar+laps_cycle_time/2)
	            else
                        i4time_hour = (i4time_radar+laps_cycle_time/2)      
     1                              /laps_cycle_time * laps_cycle_time
	            endif
c abdel	
                    var_2d = 'llr'
                    ext = 'lmt'
                    call get_laps_2dgrid(i4time_ref,laps_cycle_time*100       
     1                              ,i4time_lr,ext,var_2d,units_2d
     1                              ,comment_2d,nx_l,ny_l
     1                              ,radar_array,0,istatus)

                    call make_fnam_lp(i4time_lr,asc9_tim_r,istatus)

                    c_label = 'low lvl refl (analyzed/dbz)  '

                elseif(c_type .eq. 'ra')then
                    radar_array = field_2d
                    c_label = 'low lvl refl  (interim-vrc)  '    

                else
                    write(6,*)' error: unknown c_type'
                    stop

                endif

!               display r field
                if(c_field(3:3) .ne. 'i')then
                    call plot_cont(radar_array,1e0,0.,chigh,cint_ref
     1                        ,asc9_tim_r
     1                        ,namelist_parms,plot_parms,c_label
     1                        ,i_overlay,c_display,lat,lon
     1                        ,jdot,nx_l,ny_l,r_missing_data
     1                        ,laps_cycle_time)

                else
                    plot_parms%icol_barbs = +1 ! keep future barbs plots bright
                    call ccpfil(radar_array,nx_l,ny_l,-10.,70.,'ref'
     1                         ,n_image,1e0,'hsect',plot_parms
     1                         ,namelist_parms)
                    call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                    call setusv_dum('in',7)
                    call write_label_lplot(nx_l,ny_l,c_label
     1                                    ,asc9_tim_r,plot_parms
     1                                    ,namelist_parms
     1                                    ,i_overlay,'hsect')       
                    call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                                 ,namelist_parms,plot_parms)

                endif

            elseif(c_field(1:2) .eq. 'rf')then
                plot_parms%icol_barbs = +1 ! keep future barbs plots bright
                if(c_type .eq. 'rv')then
!                   c19_label = 'reflectivity '//ext_radar(1:3)//' '
                    call filter_string(radar_name)
                    c19_label = ' ref (dbz) '//radar_name(1:4)//' '
     1                                       //ext_radar(1:3)
                else
                    c19_label = '  reflectivity anl '
                endif

                call mklabel(k_mb,c19_label,c_label)
        
                scale = 1e0
                clow = 0.

                if(c_prodtype .eq. 'f' .and. ndim_read .eq. 2)then
                    call mk_fcst_hlabel(k_mb,comment_2d,fcst_hhmm
     1                                 ,ext(1:3),units_2d
     1                                 ,c_model,c_label)

                    scale = 1e0
                    call plot_field_2d(i4_valid,c_field
     1                        ,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,-10.,70.,cint_ref,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'ref')

                elseif(c_field(3:3) .ne. 'i')then
                    call plot_cont(grid_ra_ref(1,1,k_level,1)
     1                        ,scale,clow,chigh
     1                        ,cint_ref,asc9_tim_r,namelist_parms
     1                        ,plot_parms,c_label,i_overlay
     1                        ,c_display,lat,lon,jdot,nx_l        
     1                        ,ny_l,r_missing_data,laps_cycle_time)

                    write(6,*)' cp ref to field_2d to support diff'       
                    call move(grid_ra_ref(1,1,k_level,1),
     1                        field_2d,nx_l,ny_l) ! supports the diff option

                else
                    plot_parms%icol_barbs = +1 ! keep future barbs plots bright

                    call ccpfil(grid_ra_ref(1,1,k_level,1)
     1                         ,nx_l,ny_l,-10.,70.,'ref',n_image,1e0
     1                         ,'hsect',plot_parms,namelist_parms)
                    call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                    call setusv_dum('in',7)
                    call write_label_lplot(nx_l,ny_l,c_label
     1                                    ,asc9_tim_r,plot_parms
     1                                    ,namelist_parms
     1                                    ,i_overlay,'hsect')       
                    call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                                 ,namelist_parms,plot_parms)

                    write(6,*)' cp ref to field_2d to support diff'       
                    call move(grid_ra_ref(1,1,k_level,1),
     1                        field_2d,nx_l,ny_l) ! supports the diff option

                endif

            elseif(c_field .eq. 've')then
                call mklabel(k_mb,'  radial vel  (kt) ',c_label)

                write(6,2031)
2031            format('         enter radar # (of ones available)  '
     1                                                    ,28x,'? ',$)
                read(lun,*)i_radar

                call make_fnam_lp(i4time_radar_a(i_radar),asc9_tim_r
     1                                                      ,istatus)

                call plot_cont(grid_ra_vel(1,1,k_level,i_radar),.518
     1                        ,clow,chigh,5.,asc9_tim_r,namelist_parms       
     1                        ,plot_parms,c_label
     1                        ,i_overlay,c_display,lat,lon        
     1                        ,jdot,nx_l,ny_l,r_missing_data
     1                        ,laps_cycle_time)                          

            elseif(c_field .eq. 'vi')then
                call mklabel(k_mb,'  radial vel  (kt) ',c_label)

                write(6,*)'               # of radars available = '
     1                    ,n_radars

                write(6,2032)
2032            format('         enter radar # (of ones available - '
     1                ,'use 0 for multi-radar plot) ','? ',$)
                read(lun,*)i_radar  

                idum1_array = 0

                if(i_radar .eq. 0)then ! multi radars
                    do i_radar = 1,n_radars
                        write(6,*)'         plotting radar # ',i_radar
                        call plot_obs(k_level,.false.,asc9_tim
     1                      ,i_radar,i_radar
     1                      ,namelist_parms,plot_parms
     1                      ,nx_l,ny_l,nz_l,idum1_array,grid_ra_ref
     1                      ,grid_ra_vel(1,1,1,i_radar)    
     1                      ,lat,lon,topo,grid_spacing_m,2)
                    enddo ! i_radar

                else ! single radar
                    call make_fnam_lp(i4time_radar_a(i_radar)
     1                               ,asc9_tim_r,istatus)

                    write(6,*)' plotting radar # ',i_radar
     1                       ,' at ',asc9_tim_r

                    call plot_obs(k_level,.false.,asc9_tim_r
     1                      ,i_radar,i_radar
     1                      ,namelist_parms,plot_parms
     1                      ,nx_l,ny_l,nz_l,idum1_array,grid_ra_ref
     1                      ,grid_ra_vel(1,1,1,i_radar)    
     1                      ,lat,lon,topo,grid_spacing_m,2)

                endif

                n_image = n_image + 1

                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum('in',7)
                call write_label_lplot(nx_l,ny_l,c_label
     1                                ,asc9_tim_r,plot_parms
     1                                ,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)

            elseif(c_field(1:2) .eq. 'mt')then ! do max tops
	        if (laps_cycle_time.eq.0)then
	           i4time_hour = (i4time_radar+laps_cycle_time/2)
	        else
                   i4time_hour = (i4time_radar+laps_cycle_time/2)
     1                           /laps_cycle_time * laps_cycle_time
                endif
                var_2d = 'lmt'
                ext = 'lmt'
                call get_laps_2dgrid(i4time_hour,laps_cycle_time*100
     1                              ,i4time_lr,ext,var_2d,units_2d
     1                              ,comment_2d,nx_l,ny_l
     1                              ,radar_array,0,istatus)

                highest_top_m = 0.

                do i = 1,nx_l
                do j = 1,ny_l
                    highest_top_m = max(radar_array(i,j),highest_top_m)
                enddo ! j
                enddo ! i

!               generate contour range and interval
                cont_high = int(highest_top_m/1000.)
                cont_low = int(max(cont_high/2.,4.))

                if(cont_high - cont_low .gt. 4)then
                    cint = 2.
                else
                    cint = 1.
                endif

!               create a floor to the array for better contouring
                if(i_image .eq. 0)then
                  do i = 1,nx_l
                  do j = 1,ny_l
                    rfloor = cont_low - 0.1
                    radar_array(i,j) = radar_array(i,j) / 1000.
                    radar_array(i,j) =
     1                  max(radar_array(i,j),rfloor)
                  enddo ! j
                  enddo ! i
                  scale = 1e0
                else
                  scale = 1e3
                endif

!               display max tops
                write(6,*)' displaying max tops, cint = '
     1                          ,cont_low,cont_high,cint

!               call plot_cont(radar_array,1e0,cont_low
!    1                        ,cont_high,cint,asc9_tim_r,namelist_parms       
!    1                        ,plot_parms
!    1                        ,'max echo tops    (km msl)        '
!    1                        ,i_overlay,c_display,lat,lon
!    1                        ,jdot,nx_l,ny_l,r_missing_data
!    1                        ,laps_cycle_time)

                plot_parms%iraster = 1

                c_label = 'max echo tops    (km msl)'
                colortable = 'spectral'
                clow = 0.
                chigh = 20.
                cint = 1.
                call plot_field_2d(i4_valid,c_type,radar_array,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,colortable)

            elseif(c_field .eq. 'f1')then ! fcst max reflectivity

                if(lapsplot_pregen)then
                    write(6,*)' getting pregenerated radar fcst file'

                    var_2d = 'r06'
                    ext = 'lmr'
c abdel       
                    if (laps_cycle_time.eq.0)then	
	    	        i4time_hour = (i4time_radar+laps_cycle_time/2)
		    else  
                        i4time_hour = (i4time_radar+laps_cycle_time/2)
     1                              /laps_cycle_time * laps_cycle_time
                    endif
cabdel
                    call make_fnam_lp(i4time_hour,asc9_tim_r,istatus)
                    call get_laps_2d(i4time_hour,ext,var_2d,units_2d
     1                    ,comment_2d,nx_l,ny_l,radar_array_adv,istatus)       

                endif ! pregenerated file


!               display advected reflectivity field
                call plot_cont(radar_array_adv,1e0,0.,chigh,cint_ref,
     1          asc9_tim_r,namelist_parms,plot_parms,
     1          'max reflectivity  1 hr fcst    ',
     1          i_overlay,c_display,lat,lon,jdot,
     1          nx_l,ny_l,r_missing_data,laps_cycle_time)

!           elseif(c_field .eq. 'nt')then
!               l_radar_read = .false.
!               goto2010

            endif ! c_field

        elseif( c_type .eq. 'sr')then
            mode = 1

            call make_fnam_lp(i4time_radar,asc9_tim_r,istatus)

!           read in surface temp data
            var_2d = 't'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(
     1          i4time_radar,laps_cycle_time,i4time_temp
     1         ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1         ,temp_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' laps sfc temp not available'
                return
            endif

!           read in surface dewpoint data
            var_2d = 'td'
            ext = 'lsx'
            call get_laps_2d(i4time_temp,
     1                       ext,var_2d,units_2d,comment_2d,
     1                       nx_l,ny_l,td_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' laps sfc dewpoint not available'
                return
            endif

!           read in surface pressure data
            var_2d = 'ps'
            ext = 'lsx'
            call get_laps_2d(i4time_temp,
     1                       ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                       ,pres_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' laps sfc pressure not available'
                return
            endif

            call get_tw_approx_2d(temp_2d,td_2d,pres_2d,nx_l,ny_l,tw_sfc
     1_k)

            call zs(radar_array,temp_2d,td_2d,pres_2d,tw_sfc_k,nx_l,ny_l
     1                                                  ,snow_2d)

            c_label = 'snowfall rate         in/hr      '
            scale = 1. / (10. * (100./2.54) * 3600.) ! denom = (in snow/hr) / (m/s)

            call plot_cont(snow_2d,scale,0.,0.,-0.1,asc9_tim_r,
     1          namelist_parms,plot_parms,c_label,  
     1          i_overlay,c_display,lat,lon,jdot,
     1          nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif( c_type(1:2) .eq. 'sa' .or. c_type(1:2) .eq. 'pa' 
     1     .or. c_type(1:2) .eq. 'fa' .or. c_type(1:2) .eq. 's4')then       
            ipcp_cycle_time = precip_cycle_time
            if(c_type(1:2) .eq. 'sa')then
                write(6,1321)
1321            format('     ','enter # of hours of snow accumulation,',
     1          ' [-99 for storm total]     ','? ',$)
                var_2d = 'sto'
                ext = 'l1s'
            elseif(c_type(1:2) .eq. 'pa')then
                write(6,1322)
1322            format('     ','enter # of hours of precip accumulation,
     1',
     1          ' [-99 for storm total]   ','? ',$)
                var_2d = 'rto'
                ext = 'l1s'
            elseif(c_type(1:2) .eq. 's4')then
                write(6,1323)
1323            format('     '
     1          ,'enter # of hours of stageiv precip accum',
     1          ' [-99 for storm total]   ','? ',$)
                var_2d = 'ppt'
                ext = 'st4'
            elseif(c_type(1:2) .eq. 'fa')then
                write(6,*)' enter model for fsf'
                read(lun,*)ext
                write(6,*)'ext = ',ext
                write(6,*)' enter model initialization time'
                read(lun,*)asc9_tim
                write(6,*)'asc9_tim = ',asc9_tim
                call i4time_fname_lp(asc9_tim,i4_dum,istatus)
                write(6,1324)
1324            format('     '
     1          ,'enter # of hours of fcst precip accum',
     1          ' [-99 for storm total]   ','? ',$)
                var_2d = 'r01'
                if(trim(ext) .eq. 'nam-nh')then
                    ipcp_cycle_time = 21600
                else
                    ipcp_cycle_time = model_fcst_intvl
                endif
            endif

            read(lun,*)r_hours
            write(6,*)r_hours

            call get_directory(ext,directory,len_dir)

!           cycle over at :28 after (if input time is not on the hour)?
            if(i4time_ref .ne. (i4time_ref / 3600) * 3600)then
                i4time_ref1 = (i4time_ref-1680)/laps_cycle_time
     1                                        * laps_cycle_time
            else
                i4time_ref1 = i4time_ref
            endif

            if(r_hours .eq. -99.)then ! storm total
                write(6,*)' getting storm total accum from file'
                c9_string = 'storm tot'
                call get_laps_2dgrid(i4time_ref,10000000,i4time_stm_tot
     1                  ,ext,var_2d
     1                  ,units_2d,comment_2d,nx_l,ny_l,accum_2d,0
     1                                                  ,istatus)
                write(6,*)' storm total was reset at ',comment_2d(1:9)
                call i4time_fname_lp(comment_2d(1:9),i4time_reset,istatu
     1s)
                istatus = 1
                num_hr_accum = (i4time_stm_tot - i4time_reset) / 3600
                num_mn_accum = ((i4time_stm_tot - i4time_reset) - 
     1                           num_hr_accum*3600) / 60
                i4time_accum = i4time_stm_tot
                i4time_end = i4time_stm_tot
                i4time_start = i4time_reset

!               encode(7,2017,c7_string)min(num_hr_accum,999)
                if(num_mn_accum .eq. 0)then
                    write(c7_string,2017)min(num_hr_accum,999)
2017                format(i4,' hr')
                else
                    write(c7_string,2018)min(num_hr_accum,999)
     1                                  ,num_mn_accum
2018                format(i2,'h',i2,'m ')
                endif
     
                if(c_type(1:2) .eq. 'sa')then
                    if(c_units_type .eq. 'english')then
                        c_label = 'stm tot snow acc (in)'//c7_string
                    else ! metric
                        c_label = 'stm tot snow acc (mm)'//c7_string
                    endif
                else
                    if(c_units_type .eq. 'english')then
                        c_label = 'stm tot prcp acc (in)'//c7_string
                    else ! metric
                        c_label = 'stm tot prcp acc (mm)'//c7_string
                    endif
                endif

            elseif(.true.)then ! precip interval (via r_hours)
!               near realtime - look for snow accumulation files
!               if(i4time_now_gg() - i4time_ref1 .lt. 300)then ! real time radar
                if(.false.)then
                   !find latest time of radar data
                    if(.true.)then ! read mhr packed data
                        c_filespec = c_filespec_ra
                    endif

                    call get_file_time(c_filespec,i4time_ref1
     1                                ,i4time_accum)
                else
                    i4time_accum = i4time_ref1
                endif

                i4time_end = i4time_accum
                i4time_interval = nint(r_hours * 3600.)
                i4time_start = i4time_end - i4time_interval

                nf = 1

                call make_fnam_lp(i4time_start,a9_start,istatus)
                call make_fnam_lp(i4time_end,a9_end,istatus)

                write(6,*)' range of precip interval is ',a9_start,' '
     1                                                   ,a9_end
 
                write(6,*)' calling get_interval_precip',ipcp_cycle_time
                call get_interval_precip(var_2d(1:1),ext(1:3)
     1                                  ,i4time_start,i4time_end
     1                                  ,ipcp_cycle_time,i4_dum
     1                                  ,r_missing_data       
     1                                  ,nx_l,ny_l,nf,accum_2d,istatus)

                write(c9_string,2029)r_hours
2029            format(f5.1,' hr ')

            endif

            if(istatus .ne. 1)goto1200

            call make_fnam_lp(i4time_accum,asc9_tim_r,istatus)

            if(c_units_type .eq. 'english')then
                scale = 1. / ((100./2.54)) ! denominator = (in/m)
                if(r_hours .ne. -99.)then ! already have label for storm total
                    if(c_type(1:2) .eq. 'sa')then
                        c_label = c9_string//' snow accum  (in)'
                    elseif(c_type(1:2) .eq. 'pa')then
                        c_label = c9_string//' prcp accum  (in)'
                    elseif(c_type(1:2) .eq. 's4')then
                        c_label = c9_string//' stage iv prcp accum (in)'
                    elseif(c_type(1:2) .eq. 'fa')then
                        c_label = c9_string//' '//trim(ext)
     1                                     //' prcp accum (in)'
                    endif
                endif

            else ! metric
                scale = .001               ! numerator   = (m/mm)
                if(r_hours .ne. -99.)then ! already have label for storm total
                    if(c_type(1:2) .eq. 'sa')then
                        c_label = c9_string//' snow accum  (mm)'
                    else
                        c_label = c9_string//' prcp accum  (mm)'
                    endif
                endif
            endif ! units type

            clow = 0.

            if(c_type(1:2) .eq. 'pa' .or. c_type(1:2) .eq. 's4')then
                if(abs(r_hours) .gt. 1.0)then
                    cint = -0.01 ! -0.05
                else
                    cint = -0.01
                endif
                chigh = 10.
                if(r_hours .eq. -99.)then 
                    colortable = 'acc_sto'
                else
                    colortable = 'acc_inc'
                endif
            else ! 'sa'
                if(abs(r_hours) .gt. 1.0)then
                    cint = -0.5
                else
                    cint = -0.5
                endif
                chigh = 20.
                if(r_hours .eq. -99.)then 
                    colortable = 'sno_sto'
                else
                    colortable = 'sno_inc'
                endif
            endif

            call condition_precip(nx_l,ny_l,c_type,accum_2d,scale
     1                           ,abs(cint))

            call plot_field_2d(i4time_accum,c_type,accum_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,colortable)

        elseif( c_type .eq. 'rx')then
            write(6,1311)
1311        format('     ','enter # of hours of radar data,',
     1          ' [-99 for storm total]     ','? ',$)
            read(lun,*)r_hours

            write(6,*)

            i4time_accum = (i4time_ref+60) / 120 * 120 ! rounded off time
            i4time_end = i4time_accum
            i4time_interval = nint(r_hours * 3600.)
            i4time_start = i4time_end - i4time_interval

            if(.true.)then ! get entire time span from radar etc. data
                 write(6,*)
     1           ' getting entire time span of accumulation from radar '       
     1          ,'data, etc.'
                 call get_radar_max_pd(i4time_start,i4time_end
     1                ,nx_l,ny_l,nz_l,max_radar_files
     1                ,lat,lon,topo,grid_ra_ref
     1                ,dummy_array,radar_array,frac_sum,istatus)

            endif

!           encode(9,2029,c9_string)r_hours
            write(c9_string,2029)r_hours
!2029        format(f5.1,' hr ')

            c_label = c9_string//' reflctvty history      '

            if(istatus .ne. 1)goto1200

            call make_fnam_lp(i4time_accum,asc9_tim_r,istatus)

            scale = 1.

            if(.false.)then
                write(6,*)' writing lrx field in current directory'
                directory = '[]'
                ext = 'lrx'
                var_2d = 'lrx'
                call put_laps_2d(i4time_accum,ext,var_2d
     1          ,units_2d,comment_2d,nx_l,ny_l,radar_array,istatus)
            endif

            call plot_cont(radar_array,scale,
     1          0.,80.,cint_ref,asc9_tim_r,namelist_parms,plot_parms,
     1          c_label,i_overlay,c_display,lat,lon,jdot,
     1          nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif( c_type_i .eq. 't'  .or. c_type_i .eq. 'pt'
     1     .or. c_type_i .eq. 'bt' .or. c_type_i .eq. 'pb')then

            write(6,1513)
1513        format('     enter level in mb',48x,'? ',$)
            call input_level(lun,k_level,k_mb,pres_3d,nx_l,ny_l,nz_l)     

            if(c_type_i .eq. 'pt')then
                iflag_temp = 0 ! returns potential temperature
                call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                          ,nx_l,ny_l,nz_l,temp_3d,istatus)

                call mklabel(k_mb,'  potential temp  k'
     1                        ,c_label)      

                do i = 1,nx_l
                do j = 1,ny_l
                    field_2d(i,j) = temp_3d(i,j,k_level)
                enddo ! j
                enddo ! i

                scale = 1.
                call contour_settings(field_2d,nx_l,ny_l
     1                           ,clow,chigh,cint
     1                           ,zoom,density,scale)

            elseif(c_type_i .eq. 'pb')then
                iflag_temp = 3 ! returns balanced potential temperature
                call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                          ,nx_l,ny_l,nz_l,temp_3d,istatus)

                call mklabel(k_mb,'  balanced theta  k'
     1                                 ,c_label)   

                do i = 1,nx_l
                do j = 1,ny_l
                    field_2d(i,j) = temp_3d(i,j,k_level)
                enddo ! j
                enddo ! i

                scale = 1.
                call contour_settings(field_2d,nx_l,ny_l
     1                           ,clow,chigh,cint
     1                           ,zoom,density,scale)

            elseif(c_type_i .eq. 't ')then
                call get_temp_2d(i4time_ref,lagt,i4time_nearest
     1                          ,k_mb,nx_l,ny_l,field_2d,istatus)

                call mklabel(k_mb,' temperature      c'
     1                        ,c_label)       

            elseif(c_type_i .eq. 'bt')then
                var_2d = 't3'
                ext='lt1'

                call get_directory('balance',directory,lend)
                directory=directory(1:lend)//'lt1/'
                call get_2dgrid_dname(directory
     1           ,i4time_ref,laps_cycle_time*10000,i4time_nearest
     1           ,ext,var_2d,units_2d,comment_2d
     1           ,nx_l,ny_l,field_2d,k_mb,istatus)       


                call mklabel(k_mb,' temp (bal)       c'
     1                        ,c_label)

            endif

            if(c_type_i .eq. 't ' .or. c_type_i .eq. 'bt')then

                do i = 1,nx_l
                do j = 1,ny_l
                    field_2d(i,j) = k_to_c(field_2d(i,j))
                enddo ! j
                enddo ! i

                write(6,*)' k_mb for t plot = ',k_mb
                if(k_mb .eq. 300)then 
                    cint = 5.
                    clow = -60.
                    chigh = -30.
                elseif(k_mb .eq. 500)then
                    cint = 5.
                    clow = -40.
                    chigh = 0.
                elseif(k_mb .eq. 700)then
                    cint = 5.
                    clow = -30.
                    chigh = +30.
                elseif(k_mb .eq. 850)then
                    cint = 5.
                    clow = -10.
                    chigh = +30.
                elseif(k_mb .eq. 1000)then
                    cint = 5.
                    clow = 0.
                    chigh = +40.
                else
                    scale = 1.
                    call contour_settings(field_2d,nx_l,ny_l
     1                           ,clow,chigh,cint
     1                           ,zoom,density,scale)
                endif

            endif

            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

!           call get_pres_3d(i4time_nearest,nx_l,ny_l,nz_l,pres_3d
!    1                                     ,istatus)

!           if(pres_3d(icen,jcen,k_level) .le. 80000.)then
!               clow =  0.
!               chigh = 0.
!               cint = 2.
!           else
!               clow =  0.
!               chigh = 0.
!               cint = 5.
!           endif

            call plot_field_2d(i4time_nearest,c_type
     1                        ,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

            i4time_temp = i4time_nearest

        elseif(c_type .eq. 'hh')then
            write(6,1515)
1515        format('     enter temperature surface to display '
     1                      ,'height of (deg c)',11x,'? ',$)
            read(lun,*)temp_in_c
            temp_in_k = temp_in_c + 273.15

            iflag_temp = 1 ! returns ambient temperature

            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                          ,nx_l,ny_l,nz_l,temp_3d,istatus)
            if(istatus .ne. 1)goto1200

!           obtain height field
            ext = 'lt1'
            var_2d = 'ht'
            call get_laps_3dgrid(i4time_ref,10000000,i4time_ht
     1                          ,nx_l,ny_l,nz_l,ext,var_2d
     1                          ,units_2d,comment_2d,field_3d,istatus)
            if(istatus .ne. 1)then
                write(6,*)' error locating height field'
                return
            endif

            write(c_label,1516)nint(temp_in_c)
1516        format('height of ',i3,'c lvl (hft msl)     ')

            do j = 1,ny_l
            do i = 1,nx_l
                height_2d(i,j) = 0.

                do k = 1,nz_l-1
                    if(temp_3d(i,j,k  ) .gt. temp_in_k  .and.
     1                 temp_3d(i,j,k+1) .le. temp_in_k       )then

!                       desired temperature occurs in this layer
                        frac_k = (     temp_in_k - temp_3d(i,j,k))/
     1                         (temp_3d(i,j,k+1) - temp_3d(i,j,k))

                        height_2d(i,j) = field_3d(i,j,k) * (1.0-frac_k)       
     1                 +                 field_3d(i,j,k+1) * frac_k

                        height_2d(i,j) = height_2d(i,j) * 3.281

                    endif
                enddo ! k
            enddo ! i
            enddo ! j

            clow = 0.
            chigh = 0.
            cint = 5.

            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

            call plot_cont(height_2d,1e2,clow,chigh,cint,asc9_tim,
     1       namelist_parms,plot_parms,c_label,i_overlay,c_display,
     1       lat,lon,jdot,nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif(c_type .  eq. 'la' .or. c_type   .eq. 'lj' .or.
     1         c_type   .eq. 'sj' .or. c_type_i .eq. 'ls' .or.
     1         c_type_i .eq. 'ci' .or. c_type_i .eq. 'pc' .or.
     1         c_type_i .eq. 'rn' .or. c_type_i .eq. 'sn' .or.
     1         c_type_i .eq. 'pi' .or. c_type_i .eq. 'cn'     )then
            write(6,1514)
1514        format('     enter level in mb; or [-1] for max in column'
     1                          ,21x,'? ',$)

            call input_level(lun,k_level,k_mb,pres_3d,nx_l,ny_l,nz_l)       
c abdel       
            if (laps_cycle_time.eq.0)then
	        i4time_lwc = i4time_ref
	    else
                i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time
	    endif
            if(c_type .eq. 'la')then ! returns adiabatic lwc
                if(k_level .gt. 0)then
                    call mklabel(k_mb,
     1          ' adiabt lwc  g/m^3 ',c_label)
                else
                    c_label = 'maximum adiabatic lwc g/m^3      '
                endif

            elseif(c_type .eq. 'lj')then ! returns adjusted lwc
                if(k_level .gt. 0)then
                    call mklabel(k_mb,
     1          ' adjstd lwc  g/m^3 ',c_label)
                else
                    c_label = 'maximum adjusted  lwc g/m^3      '
                endif

            elseif(c_type .eq. 'sj')then ! returns adjusted slwc
                if(k_level .gt. 0)then
                    call mklabel(k_mb,
     1          ' adjstd slwc g/m^3 ',c_label)
                else
                    c_label = 'maximum adjusted slwc g/m^3      '
                endif

            elseif(c_type_i .eq. 'ls')then ! returns new smith - feddes lwc
                if(k_level .gt. 0)then
                    call mklabel(k_mb,
     1                             ' cloud lwc g/m^3   ',c_label)
                else
                    c_label = 'column max lwc        g/m^3      '
                endif

            elseif(c_type_i .eq. 'ci')then ! returns cloud ice
                if(k_level .gt. 0)then
                    call mklabel(k_mb,
     1                             ' cloud ice g/m^3   ',c_label)
                else
                    c_label = 'max smith-feddes  ice g/m^3      '
                endif

            elseif(c_type_i .eq. 'cn')then ! returns cloud condensate
                if(k_level .gt. 0)then
                    call mklabel(k_mb,
     1                     ' cloud condensate  g/m^3 ',c_label)
                else
                    c_label = 'max smith-feddes  ice g/m^3      '
                endif

            elseif(c_type_i .eq. 'pc')then
                if(k_level .gt. 0)then
                    call mklabel(k_mb,
     1                             ' precip conc g/m^3 ',c_label)
                endif

            elseif(c_type_i .eq. 'rn')then
                if(k_level .gt. 0)then
                    call mklabel(k_mb,
     1                             ' precip rain g/m^3 ',c_label)
                endif

            elseif(c_type_i .eq. 'sn')then
                if(k_level .gt. 0)then
                    call mklabel(k_mb,
     1                             ' precip snow g/m^3 ',c_label)
                endif

            elseif(c_type_i .eq. 'pi')then
                if(k_level .gt. 0)then
                    call mklabel(k_mb,
     1                             ' precip ice g/m^3 ',c_label)
                endif

            endif

            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            if(c_type_i .eq. 'pc')then
                var_2d = 'pcn'
                plot_parms%color_power = 0.3
            elseif(c_type_i .eq. 'rn')then
                var_2d = 'rai'
                plot_parms%color_power = 0.3
            elseif(c_type_i .eq. 'sn')then
                var_2d = 'sno'
                plot_parms%color_power = 0.3
            elseif(c_type_i .eq. 'ls')then
                var_2d = 'lwc'
                plot_parms%color_power = 0.3
            elseif(c_type_i .eq. 'ci')then
                var_2d = 'ice'
                plot_parms%color_power = 0.3
            elseif(c_type_i .eq. 'pi')then
                var_2d = 'pic'
                plot_parms%color_power = 0.3
            else
                var_2d = 'lwc'
                plot_parms%color_power = 0.3
            endif

            if(c_prodtype .eq. 'a')then
                write(6,*)' getting pregenerated lwc file'
                ext = 'lwc'
                call get_directory(ext,directory,len_dir)

                if(k_mb .eq. -1)then ! get 3d grid
                    call get_laps_3dgrid(i4time_ref,86400,i4time_cloud,       
     1                                   nx_l,ny_l,nz_l,ext,var_2d,
     1                                   units_2d,comment_2d,field_3d,
     1                                   istatus) 

                else ! get 2d horizontal slice from 3d grid
                  call get_laps_2dgrid(i4time_ref,86400,i4time_cloud,
     1                                   ext,var_2d,units_2d,
     1                                   comment_2d,nx_l,ny_l,field_2d,       
     1                                   k_mb,istatus)

                  if(c_type_i .eq. 'cn')then
                      var_2d = 'ice'
                      call get_laps_2dgrid(i4time_ref,86400,
     1                                   i4time_cloud,      
     1                                   ext,var_2d,units_2d,
     1                                   comment_2d,nx_l,ny_l,
     1                                   field_2d_buf,       
     1                                   k_mb,istatus)
                      field_2d = field_2d + field_2d_buf
                      write(6,*)' adding ice to get condensate'
                      comment_2d = 'cloud condensate'
                  endif

                endif

                i4time_lwc = i4time_cloud

                cint = -0.002

            elseif(c_prodtype .eq. 'f')then
                if(k_mb .ne. -1)then ! get 2d grid
                    call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,k_mb,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 field_2d,istatus)

                else ! get 3d grid
                    call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istatus)
                endif

                if(istatus .ne. 1)then
                    write(6,*)' could not read forecast field'       
                    goto1200
                endif

!               c_label(11:29) = ' fua '//var_2d(1:4)
!    1                             //fcst_hhmm//' g/m^3'

                call directory_to_cmodel(directory,c_model)

                units_2d = 'g/m**3'

                call mk_fcst_hlabel(k_mb,comment_2d,fcst_hhmm
     1                                 ,ext(1:3),units_2d
     1                                 ,c_model,c_label)

                i4time_lwc = i4_valid
                cint = -0.002

            else
                goto1200

            endif

            call make_fnam_lp(i4time_lwc,asc9_tim,istatus)

            write(6,*)' ascii valid time = ',asc9_tim

            clow = 0.
            chigh = 1.0

            if(k_level .gt. 0)then ! plot slwc on const pressure sfc
               if(c_prodtype .eq. 'f')then
                   call plot_field_2d(i4time_lwc,c_type
     1                        ,field_2d,1e-3
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'linear')

               elseif(c_type_i .ne. 'ci')then
                   call plot_field_2d(i4time_lwc,c_type
     1                        ,field_2d,1e-3
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'linear')

               else ! c_type .eq. 'ci'
                   call plot_field_2d(i4time_lwc,c_type
     1                        ,field_2d,1e-3
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'linear')

               endif ! c_type

            else ! find maximum value in column
               do j = 1,ny_l
               do i = 1,nx_l
                   column_max(i,j) = 0.
                   if(c_type_i .ne. 'ci')then
                     do k = 1,nz_l
                       column_max(i,j) = 
     1                 max(column_max(i,j),field_3d(i,j,k)) ! slwc_3d
                     enddo ! k
                   else
                     do k = 1,nz_l
                       column_max(i,j) = 
     1                 max(column_max(i,j),field_3d(i,j,k)) ! cice_3d
                     enddo ! k
                   endif
               enddo ! i
               enddo ! j

               call subcon(column_max,1e-30,field_2d,nx_l,ny_l)

!              call plot_cont(field_2d,1e-3,
!    1                        clow,chigh,cint,asc9_tim,
!    1                        namelist_parms,plot_parms,c_label,
!    1                        i_overlay,c_display,lat,lon,
!    1                        jdot,nx_l,ny_l,r_missing_data,
!    1                        laps_cycle_time)

               call plot_field_2d(i4time_lwc,c_type
     1                        ,field_2d,1e-3
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'linear')
            endif

        elseif(c_type .eq. 'mv' .or. c_type .eq. 'ic')then
            write(6,1514)

            call input_level(lun,k_level,k_mb,pres_3d,nx_l,ny_l,nz_l) 
c abdel       
            if (laps_cycle_time.eq.0)then	
	        i4time_lwc = i4time_ref
	    else
                i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time
	    endif
c abdel
            if(c_type .eq. 'mv')then
                if(k_level .gt. 0)then
                    call mklabel(k_mb
     1                            ,'     mvd     m^-6  ',c_label)     
                else
                    c_label = 'mean volume diameter  m^-6       '
                endif

                write(6,*)' getting pregenerated lmd file'
                var_2d = 'lmd'
                ext = 'lmd'

            elseif(c_type .eq. 'ic')then
                if(k_level .gt. 0)then
                    call mklabel(k_mb,'   icing index     '
     1                            ,c_label)
                else
                    c_label = '        laps icing index         '
                endif

                write(6,*)' getting pregenerated lrp file'
                var_2d = 'lrp'
                ext = 'lrp'

            endif ! c_type .eq. 'ic'

            if(k_mb .eq. -1)then ! get 3d grid
                call get_laps_3dgrid(i4time_ref,10000000
     1                                  ,i4time_cloud
     1                                  ,nx_l,ny_l,nz_l,ext,var_2d
     1                                  ,units_2d,comment_2d,field_3d ! slwc_3d
     1                                  ,istatus)

            else ! get 2d horizontal slice from 3d grid
                call get_laps_2dgrid(i4time_ref,10000000
     1                                  ,i4time_cloud
     1                                  ,ext,var_2d
     1                                  ,units_2d,comment_2d,nx_l
     1                                  ,ny_l,field_2d,k_mb,istatus)

            endif


            call make_fnam_lp(i4time_cloud,asc9_tim,istatus)

            if(c_type .eq. 'mv')then
                clow = 10.
                chigh = 26.
                cint = 2.

                if(k_level .gt. 0)then ! plot mvd on const pressure sfc
                   if(.true.)then
                       call subcon(mvd_2d,1e-30,field_2d,nx_l,ny_l)
                       call plot_cont(field_2d,0.9999e-6,
     1                   clow,chigh,cint,asc9_tim,namelist_parms,
     1                   plot_parms,c_label,
     1                   i_overlay,c_display,lat,lon,jdot,
     1                   nx_l,ny_l,r_missing_data,laps_cycle_time)
                   else
                       call plot_cont(field_3d(1,1,k_level),0.9999e-6, ! mvd_3d
     1                   clow,chigh,cint,asc9_tim,
     1                   namelist_parms,plot_parms,c_label,
     1                   i_overlay,c_display,lat,lon,jdot,
     1                   nx_l,ny_l,r_missing_data,laps_cycle_time)
                   endif

                else ! find maximum value in column
                   do j = 1,ny_l
                   do i = 1,nx_l
                       column_max(i,j) = -1e-30
                       do k = 1,nz_l
                           column_max(i,j) = max(column_max(i,j)
     1                                          ,field_3d(i,j,k))      ! mvd_3d
                       enddo ! k
                   enddo ! i
                   enddo ! j

                   call plot_cont(column_max,0.9999e-6,
     1               clow,chigh,cint,asc9_tim,namelist_parms,
     1               plot_parms,c_label,
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)

                endif

            elseif(c_type .eq. 'ic')then
                clow = 0.
                chigh = 10.
                cint = 1.0

                if(k_level .gt. 0)then ! plot on const pressure sfc
                   if(.true.)then
                       call plot_cont(field_2d,1e0,clow,chigh,cint
     1                      ,asc9_tim
     1                      ,namelist_parms,plot_parms,c_label
     1                      ,i_overlay,c_display     
     1                      ,lat,lon,jdot
     1                      ,nx_l,ny_l,r_missing_data,laps_cycle_time)

                   endif

                else ! find maximum value in column
                   if(.true.)then
                       do j = 1,ny_l
                       do i = 1,nx_l
                           column_max(i,j) = -1e-30
                           do k = 1,nz_l
                            if(field_3d(i,j,k) .gt. 0.)column_max(i,j) ! slwc_3d
     1                     = max(column_max(i,j),field_3d(i,j,k)+.01)       
                           enddo ! k
                       enddo ! i
                       enddo ! j
                   endif

                   call plot_cont(column_max,1e0,
     1                  clow,chigh,cint,asc9_tim,namelist_parms,
     1                  plot_parms,c_label,
     1                  i_overlay,c_display,lat,lon,jdot,
     1                  nx_l,ny_l,r_missing_data,laps_cycle_time)

                endif ! k_level

            endif ! c_type

        elseif(c_type .eq. 'cy')then
1524        write(6,1517)
1517        format('     enter lvl (mb); or [0] 2d cldtyp'
!    1          ,' [-1] low cloud,'
!    1          ,' [-2] high cloud'
     1          ,' ? ',$)

1525        read(lun,*)k_level
            k_mb = k_level

            if(k_level .lt. 0)then
                write(6,*)' try again'
                goto1524
            endif

            if(.true.)then ! read 2d cloud type field

                if(k_level .gt. 0)then ! read from 3-d cloud type
                    pressure = float(k_level*100)
                    k_level = nint(rlevel_of_field(pressure,pres_3d
     1                            ,nx_l,ny_l,nz_l,icen,jcen,istatus))
                    k_mb    = nint(pres_3d(icen,jcen,k_level) / 100.)
                    ext = 'cty'
                    var_2d = 'cty'
                    call mklabel
     1                    (k_mb,'     cloud type    ',c_label)

                else                   ! read from 2-d cloud type
                    ext = 'lct'
                    var_2d = 'sct'
                    c_label = '      laps    2-d cloud type     '

                endif


                call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1            ,i4time_cloud,ext,var_2d
     1            ,units_2d,comment_2d,nx_l,ny_l,field_2d,k_mb,istatus)

                if(istatus .ne. 1 .and. istatus .ne. -1)then
                    write(6,*)' error reading cloud type'
                    goto 1200
                endif

                call make_fnam_lp(i4time_cloud,asc9_tim,istatus)

!               convert from real to integer
                do i = 1,nx_l
                do j = 1,ny_l
                    i_array(i,j) = int(field_2d(i,j))
                enddo ! i
                enddo ! j

                call plot_cldpcp_type(i_array
     1                ,asc9_tim,namelist_parms,plot_parms
     1                ,c_label,c_type,k,i_overlay,c_display
     1                ,lat,lon,idum1_array
     1                ,nx_l,ny_l,laps_cycle_time,jdot)

            else ! old archaic code

            endif ! k_level .eq. 0


        elseif(c_type .eq. 'tp' .or. c_type .eq. 'py')then ! precip type
1624        write(6,1617)
1617        format('     enter level in mb; [0] for surface,'
     1          ,' or [-1] for sfc thresholded: ','? ',$)

1625        read(lun,*)k_level
            k_mb = k_level

            if(k_level .gt. 0)then
                call mklabel(k_mb,'    precip type    ',c_label)
                ndim=3
            elseif(k_level .eq.  0)then
                c_label = 'sfc precip type   (nothresh)     '
                ndim=2
            elseif(k_level .eq. -1)then
                c_label = 'sfc precip type   (thresh)       '
                ndim=2
            endif

            if(k_level .eq. -1)then
                var_2d = 'ptt'
                k_level = 0
            else
                var_2d = 'pty'
            endif

            if(k_level .gt. 0)then
               pressure = float(k_level*100)
               k_level = nint(rlevel_of_field(pressure,pres_3d
     1                       ,nx_l,ny_l,nz_l,icen,jcen,istatus))
               k_mb    = nint(pres_3d(icen,jcen,k_level) / 100.)
            endif
c abdel       
            if(laps_cycle_time.eq.0)then 
               i4time_pcp = i4time_ref
	    else     
               i4time_pcp = i4time_ref/laps_cycle_time * laps_cycle_time
	    endif
            l_precip_pregen = .true.

            call input_product_info(i4time_pcp              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,ndim                    ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,asc9_tim                ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            if(k_level .gt. 0)then ! plot precip type on const pressure sfc
                if(c_prodtype .eq. 'a')then
                    write(6,*)' reading pregenerated precip type field'
                    ext = 'pty'
                    call get_laps_2dgrid(i4time_pcp,laps_cycle_time
     1                    ,i4time_nearest,ext,var_2d
     1                    ,units_2d,comment_2d,nx_l,ny_l
     1                    ,field_2d,k_mb,istatus)
                    call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

                elseif(c_prodtype .eq. 'f')then
                    call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,k_mb,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 field_2d,istatus)
                    if(istatus .ne. 1)then
                        write(6,*)' could not read forecast field'       
                        goto1200
                    endif
                    c_label(11:33) = ' fua     '//var_2d(1:4)
     1                                 //fcst_hhmm//'      '

                else
                    write(6,*)' not yet supported'
                    go to 1200

                endif

            elseif(k_level .eq. 0)then ! extract surface precip type field

                if(c_prodtype .eq. 'a')then

                  ! read sfc precip type from lty field
                    write(6,*)
     1              ' reading pregenerated sfc precip type field '
     1                  ,var_2d      

!                   var_2d was defined earlier in the if block
                    ext = 'lct'
                    call get_laps_2dgrid(i4time_pcp,laps_cycle_time
     1                                  ,i4time_temp,
     1                      ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                          ,field_2d,0,istatus)
                    if(istatus .ne. 1)goto1200
                    call make_fnam_lp(i4time_temp,asc9_tim,istatus)

                elseif(c_prodtype .eq. 'f')then
                    write(6,*)' experimental model sfc precip type...'       

                    level = 0
                    var_2d = 'spt'

                    call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,level,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 field_2d,istatus)
                    if(istatus .ne. 1)then
                        write(6,*)' could not read forecast field'       
                        goto1200
                    endif
                    c_label(11:33) = ' fua     '//var_2d(1:4)
     1                                 //fcst_hhmm//'      '

                else
                    write(6,*)' not yet supported'
                    go to 1200

                endif ! c_prodtype

            endif ! k_level

!           convert from real to integer
            do i = 1,nx_l
            do j = 1,ny_l
                i_array(i,j) = field_2d(i,j)
            enddo ! i
            enddo ! j

            call plot_cldpcp_type(i_array
     1             ,asc9_tim,namelist_parms,plot_parms
     1             ,c_label,c_type,k,i_overlay,c_display  
     1             ,lat,lon,idum1_array
     1             ,nx_l,ny_l,laps_cycle_time,jdot)

        elseif(c_type_i .eq. 'ia' .or. c_type_i .eq. 'ij'
     1    .or. c_type_i .eq. 'ie' .or. c_type_i .eq. 'is'
     1    .or. c_type_i .eq. 'in' .or. c_type_i .eq. 'od'
     1    .or. c_type_i .eq. 'ca' .or. c_type_i .eq. 'sv')then       

          call input_product_info(    i4time_ref            ! i
     1                             ,laps_cycle_time         ! i
     1                             ,2                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o


          if(c_type_i .eq. 'od')then
              var_2d = 'cod'
          elseif(c_type_i .eq. 'ca')then
              var_2d = 'cla'
          elseif(c_type_i .eq. 'sv')then
              var_2d = 'smv'
          elseif(c_type_i .ne. 'ie')then
              var_2d = 'lil'
          else
              var_2d = 'lic'
          endif
          level = 0
          if(c_prodtype .eq. 'a')then
              ext = 'lil'
              call get_laps_2dgrid(i4time_ref,86400,i4time_cloud,
     1                             ext,var_2d,units_2d,comment_2d,
     1                             nx_l,ny_l,field_2d,0,istatus)

              if(c_type_i .eq. 'od')then
                  c_label = 'cloud optical depth (from hydrometeors) '
              elseif(c_type_i .eq. 'ca')then
                  c_label = 'cloud albedo (from hydrometeors) '
              elseif(c_type_i .eq. 'sv')then
                  c_label = 'simulated visible albedo '
              elseif(c_type_i .ne. 'ie')then
                  c_label = 'integrated cloud liquid (mm) '
              else
                  c_label = 'integrated cloud ice (mm)    '
              endif

!             add liquid and ice to get total condensate
              if(c_type_i .eq. 'in')then
                  field_2d_sum(:,:) = field_2d(:,:)
                  var_2d = 'lic'
                  call get_laps_2dgrid(i4time_ref,86400,i4time_cloud,     
     1                             ext,var_2d,units_2d,comment_2d,
     1                             nx_l,ny_l,field_2d,0,istatus)
                  field_2d(:,:) = field_2d_sum(:,:) + field_2d(:,:)
                  c_label = 'integrated cloud condensate (mm)    '
              endif

          elseif(c_prodtype .eq. 'f')then
              call read_laps(i4_initial,i4_valid,directory,
     1                       ext,nx_l,ny_l,1,1,       
     1                       var_2d,level,lvl_coord_2d,
     1                       units_2d,comment_2d,
     1                       field_2d,istatus)

              if(istatus .ne. 1)then
                  write(6,*)' could not read forecast field'       
                  goto1200
              endif

              call directory_to_cmodel(directory,c_model)

              call mk_fcst_hlabel(level,comment_2d,fcst_hhmm
     1                                 ,ext(1:3),units_2d
     1                                 ,c_model,c_label)

!             add liquid and ice to get total condensate
              if(c_type_i .eq. 'in')then
                  field_2d_sum(:,:) = field_2d(:,:)
                  var_2d = 'lic'
                  call read_laps(i4_initial,i4_valid,directory,
     1                       ext,nx_l,ny_l,1,1,       
     1                       var_2d,level,lvl_coord_2d,
     1                       units_2d,comment_2d,
     1                       field_2d,istatus)

                  if(istatus .ne. 1)then
                      write(6,*)' could not read forecast field'       
                      goto1200
                  endif

                  field_2d(:,:) = field_2d_sum(:,:) + field_2d(:,:)

                  comment_2d = 'integrated cloud condensate'
                  units_2d = 'mm'
                  call mk_fcst_hlabel(level,comment_2d,fcst_hhmm
     1                                 ,ext(1:3),units_2d
     1                                 ,c_model,c_label)
              endif

              i4time_cloud = i4_valid

          endif

          call make_fnam_lp(i4time_cloud,asc9_tim,istatus)

          clow = 0.
          cint = -0.1

          if(c_type_i .eq. 'od')then
              scale = 1.
              chigh = +160.
              cint = -0.2
              colortable = 'linear'
              plot_parms%color_power = 0.50
          elseif(c_type_i .eq. 'ca')then
              scale = 1.
              chigh = +1.
              cint = 0.2
              colortable = 'linear'
              plot_parms%color_power = 1.00
          elseif(c_type_i .eq. 'sv')then
              scale = 1.
              chigh = +1.
              cint = 0.2
              colortable = 'linear'
              plot_parms%color_power = 1.00
          else
              scale = 1e-3 ! data are in m, plot is in mm
              chigh = +2.
              colortable = 'tpw'
              plot_parms%color_power = 0.33
          endif

          call plot_field_2d(i4time_cloud,c_type_i,field_2d
     1                        ,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,colortable)


        elseif(c_type(1:2) .eq. 'pe' .or. c_type(1:2) .eq. 'ne')then        
          ext = 'lst'

          if(c_type(1:2) .eq. 'pe')then
              var_2d = 'pbe'

              call get_laps_2dgrid(i4time_ref,10000000,i4time_temp,
     1        ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                          ,field_2d,0,istatus)

          else
              var_2d = 'nbe'

              call get_laps_2dgrid(i4time_ref,10000000,i4time_temp,
     1        ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                          ,field_2d,0,istatus)

          endif

          call make_fnam_lp(i4time_temp,asc9_tim,istatus)


          scale = 1.
          call make_fnam_lp(i4time_temp,asc9_tim,istatus)

          if(c_type(1:2) .eq. 'pe')then
              call condition_cape(nx_l,ny_l,c_type,r_missing_data
     1                           ,field_2d)

              c_label = 'sbcape              (j/kg)       '
              clow = 0.
              chigh = chigh_cape
              cint = +500.

          elseif(c_type(1:2) .eq. 'ne')then
!             change flag value (for now)
              do i = 1,nx_l
              do j = 1,ny_l

!                 if(field_2d(i,j) .eq. -1e6)then
!                 if(abs(field_2d(i,j)) .ge. +1e6)then
!                     field_2d(i,j) = r_missing_data
!                 endif

                  if(field_2d(i,j) .ge. -2.)then
                      if(field_2d(i,j) .eq. r_missing_data)then
                          field_2d(i,j) = +999. ! r_missing_data
                      else
                          field_2d(i,j) = +0.1
                      endif

                  elseif(field_2d(i,j) .lt. -500.)then
                      field_2d(i,j) = -500.                      

                  endif

              enddo ! j
              enddo ! i

              c_label = 'cin                 (j/kg)       '
              clow = -500  !   0.
              chigh = 50.  !  50.
              cint =  50.  ! -10.

          endif


          if(c_type(3:3) .ne. 'i')then ! contour plot
!             call contour_settings(field_2d,nx_l,ny_l
!    1                               ,clow,chigh,cint,zoom,density,1.)       

              call plot_cont(field_2d,scale,clow,chigh,cint,asc9_tim       
     1           ,namelist_parms,plot_parms,c_label,i_overlay
     1           ,c_display,lat,lon,jdot
     1           ,nx_l,ny_l,r_missing_data,laps_cycle_time)

          else ! image plot
              where(field_2d .eq. r_missing_data)field_2d = 0.
              call ccpfil(field_2d,nx_l,ny_l,clow,chigh,'cpe'
     1                   ,n_image,scale,'hsect',plot_parms
     1                   ,namelist_parms)    
              call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
              call setusv_dum('in',7)
              call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                              ,plot_parms,namelist_parms,i_overlay       
     1                              ,'hsect')
              call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                           ,namelist_parms,plot_parms)

          endif ! image plot

c
c j. smart - 4/19/99. updated moisture plotting. in addition, added two
c                     more switches for lga/fua plotting.
c
        elseif(c_type(1:2) .eq. 'lq' .or. c_type(1:2).eq.'rb')then
c
c j. smart - 4/19/99. lq is laps-lq3 either sh or rh
c    "       4/25/02  rb is laps-balance lq3  "
c
            print*,'you selected plotting of 3d humidity data '

            write(6,1513)
            call input_level(lun,k_level,k_mb,pres_3d,nx_l,ny_l,nz_l)       

            write(6,1615)
1615        format(10x,'plot rh, q, or td [r/q/d]  ? ',$)
            read(5,*)qtype

            if(qtype .eq. 'q' .or. qtype .eq. 'd')then ! q or td

              var_2d = 'sh '
              ext = 'lq3'

              plot_parms%iraster = 1

              if(c_type(1:2).eq.'rb')then ! balanced
                 call get_directory('balance',directory,lend)
                 directory=directory(1:lend)//'lq3/'
                 call get_2dgrid_dname(directory
     1             ,i4time_ref,laps_cycle_time*100,i4time_heights
     1             ,ext,var_2d,units_2d,comment_2d
     1             ,nx_l,ny_l,field_2d,k_mb,istatus)

                 if(qtype .eq. 'q')then
                   call mklabel(k_mb,' spec hum  (bal) g/kg',c_label)
                 else
                   call mklabel(k_mb,' dewpoint  (bal) ',c_label)
                 endif

              else ! analysis
                 call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                             ,i4time_heights
     1                             ,ext,var_2d,units_2d,comment_2d
     1                             ,nx_l,ny_l,field_2d,k_mb,istatus)

                 if(qtype .eq. 'q')then
                   call mklabel(k_mb,' spec hum (anal) g/kg',c_label)
                 else
                   call mklabel(k_mb,' dewpoint (anal) ',c_label)
                 endif

              endif   

              if(istatus.ne.1 .and. istatus.ne.-1)then
                 print*,'no plotting for the requested time period'

              else

                if(qtype .eq. 'q')then
                 clow = 0.
                 chigh = +25.
                 cint = 1.0

                 cint = 0.0 ! -0.1

                 scale = 1e-3

                 plot_parms%color_power = 0.7

                 call plot_field_2d(i4time_heights,c_type_i,field_2d
     1                        ,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'tpw')

                else ! dewpoint
                  write(6,*)' calculate analysis tdew - level ',k_level
                  do i = 1,nx_l
                  do j = 1,ny_l
!                     field_2d(i,j) = make_td(
!    1                                pres_3d(i,j,k_level)*.01,99.
!    1                               ,field_2d(i,j)*.001,-132.)
                      field_2d(i,j) = tdew(pres_3d(i,j,k_level)*.01
     1                                    ,field_2d(i,j))
                  enddo ! j
                  enddo ! i                       

                  call contour_settings(field_2d,nx_l,ny_l
     1                            ,clow,chigh,cint,zoom,density,scale)      
                  clow = -100.
                  chigh = +30.
!                 cint = 10.

                  call plot_field_2d(i4time_heights,c_type_i,field_2d
     1                        ,1.
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'moist')

                endif 

              endif

            elseif(qtype .eq. 'r')then ! rh
              ext = 'lh3'
              if(c_type(1:2).eq.'rb')then
                 var_2d = 'rhl'
                 write(6,*)' reading rhl / ',var_2d
                 call mklabel(k_mb,' balanced rh (liq) %'
     1                                 ,c_label)
                 call get_directory('balance',directory,lend)
                 directory=directory(1:lend)//'lh3/'
                 call get_2dgrid_dname(directory
     1             ,i4time_ref,laps_cycle_time*100,i4time_heights
     1             ,ext,var_2d,units_2d,comment_2d
     1             ,nx_l,ny_l,field_2d,k_mb,istatus)

              else
                 write(6,1616)
1616             format(10x,'plot rh3 or rhl [3/l]? ',$)
                 read(5,*)qtype


                 if(qtype .eq. '3')then
                  var_2d = 'rh3'
                  write(6,*)' reading rh3 / ',var_2d
                  call mklabel(k_mb,' laps rh     (rh3) %'
     1                                 ,c_label)     
                 elseif(qtype .eq. 'l')then
                  var_2d = 'rhl'
                  write(6,*)' reading rhl / ',var_2d
                  call mklabel(k_mb,' laps rh     (liq) %'
     1                                 ,c_label)     
                 endif
                 call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                             ,i4time_heights
     1                             ,ext,var_2d,units_2d,comment_2d
     1                             ,nx_l,ny_l,field_2d,k_mb,istatus)

                 if(qtype .eq. '3')then
                    call mklabel(k_mb,' laps rh     (rh3) %'
     1                                 ,c_label)
                 elseif(qtype .eq. 'l')then
                    call mklabel(k_mb,' laps rh     (liq) %'
     1                                 ,c_label)
                 endif
              endif
 
              if(istatus.ne. 1)then
                 print*,'no plotting for the requested time period'
              else

                clow = 0.
                chigh = +100.
                cint = 10.

                call make_fnam_lp(i4time_heights,asc9_tim,istatus)

                scale = 1e0
                call plot_field_2d(i4time_heights,c_type
     1                        ,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'moist')

              endif
            endif

        elseif(c_type_i .eq. 'br' .or. c_type_i .eq. 'fr')then
c
c j. smart - 4/19/99. br is laps-lga either sh or rh (sh is converted to rh).
c
            ext='lga'
            if(c_type_i.eq.'fr')ext='fua'

            print*,'      plotting ',ext(1:3),' humidity data'

            call input_background_info(
     1                              ext                     ! i
     1                             ,directory,c_model       ! o
     1                             ,i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,asc9_tim                ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o
            if(istatus.ne.1)goto1200

            write(6,1513)
            call input_level(lun,k_level,k_mb,pres_3d,nx_l,ny_l,nz_l)       

            write(6,1615)
            read(5,*)qtype

            var_2d = 'sh '
            call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,k_mb,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 sh_2d,istat_sh)

            var_2d = 'rh3'
            call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,k_mb,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 rh_2d,istat_rh)

!           test for valid (non-missing) rh
            if(istat_rh .eq. 1)then
                istat_rh = 0
                do i=1,nx_l
                do j=1,ny_l
                    if(rh_2d(i,j) .ne. r_missing_data)then
                        istat_rh = 1
                    endif
                enddo ! j
                enddo ! i
            endif

            if(istat_rh .ne. 1 .and. istat_sh .ne. 1)then
                print*,' rh/sh not obtained from ',ext(1:3)
                print*,'no plotting of data for requested time period'
                goto1200
            endif


            if(.true.)then

                if(qtype.eq.'q' .and. istat_sh .eq. 1)then
                    write(6,*)' plotting q directly, range is: ',
     1                        minval(sh_2d),maxval(sh_2d)

!                   call mklabel(k_mb,' '//fcst_hhmm
!    1                         //' '//ext(1:3)//' q  (x1e3)',c_label)
                    comment_2d = 'q'     ! for label
                    units_2d = '(x1e3)'  ! for label
                    call mk_fcst_hlabel(k_mb,comment_2d,fcst_hhmm
     1                                 ,ext(1:3),units_2d
     1                                 ,c_model,c_label)

                    clow = 0.
                    chigh = +25.
                    cint = 1.0
c                   cint = -1.

                    scale = 1e-3
                    plot_parms%color_power = 0.7

                    call plot_field_2d(i4_valid,c_type_i,sh_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'tpw')

                    call move(sh_2d,field_2d,nx_l,ny_l) ! supports diff option

                elseif(qtype .eq. 'd' .and. istat_sh .eq. 1)then
                    write(6,*)' calculate model tdew from p and sh'
                    write(6,*)' q range is: '
     1                       ,minval(sh_2d),maxval(sh_2d)

                    do i = 1,nx_l
                    do j = 1,ny_l
!                       field_2d(i,j) = make_td(
!    1                                  pres_3d(i,j,k_level)*.01,99.
!    1                                 ,sh_2d(i,j)*.001,-132.)
                        field_2d(i,j) = tdew(pres_3d(i,j,k_level)*.01
     1                                      ,sh_2d(i,j))
                    enddo ! j
                    enddo ! i                       

                    clow = -100.
                    chigh = +30.
                    cint = 10.

                    comment_2d = 'dewpoint         (c)'

                    call mk_fcst_hlabel(k_mb,comment_2d,fcst_hhmm
     1                                 ,ext(1:3),units_2d
     1                                 ,c_model,c_label)

                    call plot_field_2d(i4_valid,c_type_i,field_2d       
     1                                ,1.
     1                                ,namelist_parms,plot_parms
     1                                ,clow,chigh,cint,c_label
     1                                ,i_overlay,c_display,lat,lon,jdot
     1                                ,nx_l,ny_l,r_missing_data,'moist')

                elseif(qtype .eq. 'r')then
                    if(istat_sh .eq. 1 .and. istat_rh .ne. 1)then

                        write(6,1635)
1635                    format(10x
     1                       ,'input t_ref for rh calc [deg c] ? ',$)
                        read(5,*)t_ref

                        var_2d = 't3 '
                        call read_laps(i4_initial,i4_valid,directory,
     1                                 ext,nx_l,ny_l,1,1,       
     1                                 var_2d,k_mb,lvl_coord_2d,
     1                                 units_2d,comment_2d,
     1                                 temp_2d,istatus)

                        if(istatus.ne.1)then
                            print*,var_2d, ' not obtained from '
     1                            ,ext(1:3)
                        endif

!                       call mklabel(k_mb,' '//fcst_hhmm
!    1                         //' '//ext(1:3)//' rh %cptd ',c_label)

                        comment_2d = 'rh'   ! for label
                        units_2d = '%cptd'  ! for label

                        call make_fnam_lp(i4_valid,asc9_tim,istatus)

                        do i = 1,nx_l
                        do j = 1,ny_l
                            rh_2d(i,j)=make_rh(float(k_mb)
     1                         ,temp_2d(i,j)-273.15
     1                         ,sh_2d(i,j)*1000.,t_ref)*100. ! q_3d
                        enddo ! j
                        enddo ! i

                    elseif(istat_rh .eq. 1)then
                        write(6,1636)
1636                    format(10x,'ok to plot rh as read in ? ',$)
                        read(5,*)directory   

                        if(directory(1:1) .eq. 'n' 
     1                .or. directory(1:1) .eq. 'n')then
                            goto1200
                        endif

                        units_2d = '%' ! for label

                    else
                        write(6,*)' rh/sh not obtained...'
                        goto1200

                    endif ! istat_rh / istat_sh

                    call mk_fcst_hlabel(k_mb,comment_2d,fcst_hhmm
     1                                 ,ext(1:3),units_2d
     1                                 ,c_model,c_label)

                    clow = 0.
                    chigh = +100.
                    cint = 10.

                    field_2d = rh_2d ! supports diff option better

                    call plot_field_2d(i4_valid,c_type_i,field_2d,1e0       
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'moist')

                endif ! plot rh
            endif ! true

        elseif(c_type_i .eq. 'hb' .or. c_type_i .eq. 'tb' .or.
     1         c_type_i .eq. 'hr' .or. c_type_i .eq. 'tr'     )then
            
            if(c_type_i(1:1) .eq. 'h')then
                var_2d = 'ht'
            else
                var_2d = 't3'
            endif

            if(c_type_i(2:2) .eq. 'b')then
                ext = 'lga'
            else
                ext = 'fua'
            endif

            call input_background_info(
     1                              ext                     ! i
     1                             ,directory,c_model       ! o
     1                             ,i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,asc9_tim                ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o
            if(istatus.ne.1)goto1200

!           call get_pres_3d(i4_valid,nx_l,ny_l,nz_l,field_3d,istatus)       

            write(6,1513)
            call input_level(lun,k_level,k_mb,pres_3d,nx_l,ny_l,nz_l)       

            call read_laps(i4_initial,i4_valid,directory,ext,
     1          nx_l,ny_l,1,1,       
     1          var_2d,k_mb,lvl_coord_2d,units_2d,comment_2d,
     1          field_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading grid ',var_2d,' ',ext,istatus
                goto1200
            endif


            if(c_type_i(1:1) .eq. 'h')then
                scale = 10.

!               call mklabel(k_mb,' laps '//ext(1:3)//' height dm'
!    1                        ,c_label)

!               call mklabel(k_level,ext(1:3)//' '
!    1                         //fcst_hhmm//' fcst ht dm',c_label)

!               call mklabel(k_mb,' '//fcst_hhmm
!    1                         //' '//ext(1:3)//' height dm',c_label)

                call mk_fcst_hlabel(k_mb,'height',fcst_hhmm
     1                                 ,ext(1:3),'dm'
     1                                 ,c_model,c_label)

                clow = 0.
                chigh = 0.
                call array_range(field_2d,nx_l,ny_l,rmin,rmax
     1                          ,r_missing_data)

                range = (rmax-rmin) / scale
                range_eff = range / density
 
                if(range_eff .gt. 40)then
                    cint = 6.
                elseif(range_eff .gt. 20)then
                    cint = 3.
                elseif(range_eff .gt. 8)then
                    cint = 2.
                else ! range_eff < 8
                    cint = 1.
                endif

            else  
!               call mklabel(k_mb,' '//fcst_hhmm
!    1                         //' '//ext(1:3)//' temp    c',c_label)

                call mk_fcst_hlabel(k_mb,'temperature',fcst_hhmm
     1                                 ,ext(1:3),'deg c'
     1                                 ,c_model,c_label)
                scale = 1.

                do i = 1,nx_l
                do j = 1,ny_l
                    field_2d(i,j) = field_2d(i,j) - 273.15
                enddo ! j
                enddo ! i

                write(6,*)' k_mb for t fua/lga plot = ',k_mb
                if(vert_grid .eq. 'sigma_ht')then
                    call contour_settings(field_2d,nx_l,ny_l
     1                                   ,clow,chigh,cint
     1                                   ,zoom,density,scale)       
                elseif(k_mb .eq. 300)then 
                    cint = 5.
                    clow = -60.
                    chigh = -30.
                elseif(k_mb .eq. 500)then
                    cint = 5.
                    clow = -40.
                    chigh = 0.
                elseif(k_mb .eq. 700)then
                    cint = 5.
                    clow = -30.
                    chigh = +30.
                elseif(k_mb .eq. 850)then
                    cint = 5.
                    clow = -10.
                    chigh = +30.
                elseif(k_mb .eq. 1000)then
                    cint = 5.
                    clow = 0.
                    chigh = +40.
                else
                    call contour_settings(field_2d,nx_l,ny_l
     1                                   ,clow,chigh,cint
     1                                   ,zoom,density,scale)       
                endif


            endif

            call make_fnam_lp(i4_valid,asc9_tim,istatus)

            call plot_field_2d(i4_valid,c_type_i,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

            i4time_temp = i4_valid

        elseif(c_type .eq. 'to')then
            write(6,1513)
            call input_level(lun,k_level,k_mb,pres_3d,nx_l,ny_l,nz_l)       

            if(laps_cycle_time .eq. 0)then
	    i4time_temp = i4time_ref 
	    else
                i4time_temp = (i4time_ref / laps_cycle_time) 
     1                        * laps_cycle_time
            endif

!           plot temperature obs from the tmg file
            call plot_temp_obs(k_level,i4time_temp,nx_l,ny_l,nz_l
     1                        ,r_missing_data,lat,lon,topo,zoom
     1                        ,plot_parms)

        elseif(c_type .eq. 'ho' .or. c_type .eq. 'qo')then
            write(6,1513)
            call input_level(lun,k_level,k_mb,pres_3d,nx_l,ny_l,nz_l)       

            if(laps_cycle_time.eq. 0)then
	    i4time_temp = i4time_ref
	    else
                i4time_temp = (i4time_ref / laps_cycle_time) 
     1                        * laps_cycle_time
            endif

!           read dewpoint obs from the snd file (and make hmg file)
            write(6,*)' reading snd file and making hmg file'
            call plot_td_sndobs(k_level,i4time_temp,nx_l,ny_l,nz_l
     1                         ,r_missing_data,lat,lon,topo,zoom
     1                         ,plot_parms)

            if(c_type .eq. 'ho')then
                mode = 1 ! dewpoint obs
            else
                mode = 2 ! sh obs (from dewpoint)
            endif

!           plot obs from the hmg file 
            write(6,*)' plotting dewpoint from hmg file'
            call plot_td_obs(k_level,i4time_temp,nx_l,ny_l,nz_l
     1                      ,r_missing_data,lat,lon,topo,zoom
     1                      ,namelist_parms,plot_parms
     1                      ,k_mb,mode,field_2d,i_overlay)

        elseif(c_type .eq. 'wv')then
            k_level = 0
            k_mb = 0

            if(laps_cycle_time .eq. 0)then
	    i4time_temp = i4time_ref
	    else
                i4time_temp = (i4time_ref / laps_cycle_time) 
     1                        * laps_cycle_time
            endif

!           read dewpoint obs from the snd file (and make hmg file)
            write(6,*)' reading snd file and gps to make hmg file'
            call plot_td_sndobs(k_level,i4time_temp,nx_l,ny_l,nz_l
     1                         ,r_missing_data,lat,lon,topo,zoom
     1                         ,plot_parms)

!           mode = 3 ! pw obs (from gps)
            mode = 4 ! pw obs difference from field_2d (from gps)

!           use this if the hmg file becomes available
            write(6,*)' plotting gps wv obs from hmg file'
            call plot_td_obs(k_level,i4time_temp,nx_l,ny_l,nz_l
     1                      ,r_missing_data,lat,lon,topo,zoom
     1                      ,namelist_parms,plot_parms
     1                      ,k_mb,mode,field_2d,i_overlay)

        elseif(c_type .eq. 'il')then
            k_level = 0
            k_mb = 0

            if(laps_cycle_time .eq. 0)then
	    i4time_temp = i4time_ref
	    else
                i4time_temp = (i4time_ref / laps_cycle_time) 
     1                        * laps_cycle_time
            endif

            write(6,*)' plotting integrated liquid obs from snd file'
            call plot_il_obs(k_level,i4time_temp,nx_l,ny_l,nz_l
     1                      ,r_missing_data,lat,lon,topo,zoom
     1                      ,i_overlay,namelist_parms,plot_parms)

        elseif(c_type(1:2) .eq. 'ht'.or. c_type(1:2) .eq. 'bh')then
            write(6,1513)
            call input_level(lun,k_level,k_mb,pres_3d,nx_l,ny_l,nz_l)       

            var_2d = 'ht'

            ext='lt1'

            if(c_type(1:2) .eq. 'bh' )then
               call get_directory('balance',directory,lend)
               directory=directory(1:lend)//'lt1/'
               call get_2dgrid_dname(directory
     1             ,i4time_ref,laps_cycle_time*100,i4time_heights
     1             ,ext,var_2d,units_2d,comment_2d
     1             ,nx_l,ny_l,field_2d,k_mb,istatus)

               call mklabel(k_mb,' height  (bal)   dm',c_label)       

            else ! 'ht'
               call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                             ,i4time_heights
     1                             ,ext,var_2d,units_2d,comment_2d
     1                             ,nx_l,ny_l,field_2d,k_mb,istatus)       
               call mklabel(k_mb,' height          dm',c_label)       

            endif

            if(istatus .ne. 1)then
                write(6,*)' error reading laps height analysis'
                goto1200
            endif

            scale = 10.

            clow =  0.
            chigh = 0.
            call array_range(field_2d,nx_l,ny_l,rmin,rmax
     1                      ,r_missing_data)

            range = (rmax-rmin) / scale 
            range_eff = range / density

            if(range_eff .gt. 40)then
                cint = 6.
            elseif(range_eff .gt. 20)then
                cint = 3.
            elseif(range_eff .gt. 8)then
                cint = 2.
            else ! range_eff < 8
                cint = 1.
            endif

            call make_fnam_lp(i4time_heights,asc9_tim,istatus)

            call plot_field_2d(i4time_heights,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

        elseif(c_type(1:2) .eq. 'pw')then
            var_2d = 'tpw'
            ext = 'lh4'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,
     1              ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                     ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading precipitable water'
                goto1200
            endif

            c_label = 'total precipitable water  cm     '

!           cint = 0.25
            scale = 1e-2 ! data in m, plot in cm

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call contour_settings(field_2d,nx_l,ny_l,clow,chigh,cint
     1                           ,zoom,density,scale)       

            clow = 0.00
            chigh = namelist_parms%chigh_tpw
            plot_parms%color_power = namelist_parms%power_tpw

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'tpw')

        elseif(c_type(1:2) .eq. 'tt' 
     1    .or. c_type(1:2) .eq. 'tf' .or. c_type(1:2) .eq. 'gf' 
     1    .or. c_type(1:2) .eq. 'gc' .or. c_type(1:2) .eq. 'tc')then

            if(c_type(1:1) .eq. 't')then
                var_2d = 't'
            elseif(c_type(1:1) .eq. 'g')then
                var_2d = 'tgd'
            endif

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)

            if(c_type(1:1) .eq. 't')then
                c_label = 'sfc temperature   '
            elseif(c_type(1:1) .eq. 'g')then
                c_label = 'ground temperature'
            endif

            call s_len2(c_label,len_c_label)

            if(c_units_type .eq. 'english')then
                do i = 1,nx_l
                do j = 1,ny_l
                    field_2d(i,j) = k_to_f(field_2d(i,j))
                enddo ! j
                enddo ! i

                c_label = c_label(1:len_c_label)//' (f)'
                sfct_l = sfctf_l
                sfct_h = sfctf_h

            elseif(c_units_type .ne. 'english')then
                do i = 1,nx_l
                do j = 1,ny_l
                    field_2d(i,j) = k_to_c(field_2d(i,j))
                enddo ! j
                enddo ! i

                c_label = c_label(1:len_c_label)//' (c)'
                sfct_l = sfctc_l
                sfct_h = sfctc_h

            endif

!           add balance to label
            call s_len2(c_label,len_c_label)
            if(i_balance .eq. 1)then
                c_label = c_label(1:len_c_label)//' (bal)'
            endif

            if(istatus .ne. 1)then
                write(6,*)' error reading surface temps'
                goto1200
            endif

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            if(c_type(3:3) .ne. 'i')then ! contour plot
                call contour_settings(field_2d,nx_l,ny_l
     1                               ,clow,chigh,cint,zoom,density,1.)       

                call plot_cont(field_2d,1e-0,clow,chigh,cint,asc9_tim       
     1           ,namelist_parms,plot_parms,c_label,i_overlay
     1           ,c_display,lat,lon,jdot
     1           ,nx_l,ny_l,r_missing_data,laps_cycle_time)

            else ! image plot
                call ccpfil(field_2d,nx_l,ny_l,sfct_l,sfct_h,'temp'       
     1                     ,n_image,1e-0,'hsect',plot_parms
     1                     ,namelist_parms)    
                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum('in',7)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)

            endif

        elseif(c_type(1:2) .eq. 'td' .or. c_type(1:2) .eq. 'df'
     1                               .or. c_type(1:2) .eq. 'dc')then
            var_2d = 'td'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1             ,i4time_pw
     1             ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                     ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface td'
                goto1200
            endif

            if(c_units_type .eq. 'english')then
                do i = 1,nx_l
                do j = 1,ny_l
                    field_2d(i,j) = k_to_f(field_2d(i,j))
                enddo ! j
                enddo ! i

                c_label = 'sfc dew point       (f)          '
                sfctd_l = sfctdf_l
                sfctd_h = sfctdf_h

            elseif(c_units_type .ne. 'english')then
                do i = 1,nx_l
                do j = 1,ny_l
                    field_2d(i,j) = k_to_c(field_2d(i,j))
                enddo ! j
                enddo ! i
                c_label = 'sfc dew point       (c)          '
                sfctd_l = sfctdc_l
                sfctd_h = sfctdc_h

            endif

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

!           add balance to label
            call s_len2(c_label,len_c_label)
            if(i_balance .eq. 1)then
                c_label = c_label(1:len_c_label)//' (bal)'
            endif

            if(c_type(3:3) .ne. 'i')then ! contour plot
                call contour_settings(field_2d,nx_l,ny_l
     1                               ,clow,chigh,cint,zoom,density,1.)       

                call plot_cont(field_2d,scale,clow,chigh,cint
     1               ,asc9_tim,namelist_parms,plot_parms
     1               ,c_label,i_overlay,c_display
     1               ,lat,lon,jdot
     1               ,nx_l,ny_l,r_missing_data,laps_cycle_time)

            else ! image plot
                call ccpfil(field_2d,nx_l,ny_l,sfctd_h,sfctd_l,'hues'       
     1                     ,n_image,scale,'hsect',plot_parms
     1                     ,namelist_parms)    
                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum('in',7)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)

            endif

        elseif(c_type(1:2) .eq. 'hi')then
            var_2d = 'hi'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,
     1              ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                     ,field_2d,0,istatus)

!           k to f
            do i = 1,nx_l
            do j = 1,ny_l
                if(field_2d(i,j) .ne. r_missing_data)then
                    field_2d(i,j) = k_to_f(field_2d(i,j))
                endif
            enddo ! j
            enddo ! i

            if(istatus .ne. 1 .and. istatus .ne. -1)then
                write(6,*)' error reading heat index'
                goto1200
            endif

            c_label = 'heat index          (f)          '

            scale = 1.

            clow = sfctf_l
            chigh = sfctf_h
            cint = 0.

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'temp')

        elseif(c_type(1:2) .eq. 'mc')then
            var_2d = 'mrc'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1             ,i4time_pw,
     1              ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                     ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface moisture convergence'
                goto1200
            endif

            c_label = 'sfc mstr flux conv  (x 1e-4 s-1) '

            clow = -40.
            chigh = +40.
            cint = 10.

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            scale = 1e-4

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,chigh,clow,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'spectralr')

        elseif(c_type .eq. 'ws')then ! surface wind
            var_2d = 'u'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,u_2d,0,istatus)      
            var_2d = 'v'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,v_2d,0,istatus)      

            if(istatus .ne. 1)then
                write(6,*)' error reading surface wind'
                goto1200
            endif

            c_label = 'surface wind            (kt)     '

            nxz = float(nx_l) / zoom
            nyz = float(ny_l) / zoom

            interval = int(max(nxz,nyz) / 85.) + 1
            size = float(interval) * .15

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_barbs(u_2d,v_2d,lat,lon,topo,size,zoom,interval       
     1                     ,asc9_tim,namelist_parms,plot_parms
     1                     ,c_label,c_field,k_level
     1                     ,i_overlay,c_display
     1                     ,nx_l,ny_l,nz_l,max_radars
!    1                     ,grid_ra_ref_dum,grid_ra_vel_dum
     1                     ,nx_l,ny_l,r_missing_data,laps_cycle_time
     1                     ,jdot)

        elseif(c_type .eq. 'wp')then ! planetary boundary layer mean wind
            ext = 'lfr'
            var_2d = 'upb'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,u_2d,0,istatus)      
            var_2d = 'vpb'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,v_2d,0,istatus)      

            if(istatus .ne. 1)then
                write(6,*)' error reading pbl wind'
                goto1200
            endif

            c_label = 'pbl mean wind            (kt)    '

            nxz = float(nx_l) / zoom
            nyz = float(ny_l) / zoom

            interval = int(max(nxz,nyz) / 65.) + 1
            size = float(interval) * .15

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_barbs(u_2d,v_2d,lat,lon,topo,size,zoom,interval       
     1                     ,asc9_tim,namelist_parms,plot_parms
     1                     ,c_label,c_field,k_level
     1                     ,i_overlay,c_display
     1                     ,nx_l,ny_l,nz_l,max_radars
!    1                     ,grid_ra_ref_dum,grid_ra_vel_dum
     1                     ,nx_l,ny_l,r_missing_data,laps_cycle_time
     1                     ,jdot)

        elseif(c_type .eq. 'w6')then ! 0-6km shear
            ext = 'lhe'
            var_2d = 'shu'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,u_2d,0,istatus)      
            var_2d = 'shv'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,v_2d,0,istatus)      

            if(istatus .ne. 1)then
                write(6,*)' error reading 0-6km shear vector'
                goto1200
            endif

            c_label = '0-6km agl wind shear     (kt)    '

            nxz = float(nx_l) / zoom
            nyz = float(ny_l) / zoom

            interval = int(max(nxz,nyz) / 65.) + 1
            size = float(interval) * .15

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_barbs(u_2d,v_2d,lat,lon,topo,size,zoom,interval       
     1                     ,asc9_tim,namelist_parms,plot_parms
     1                     ,c_label,c_field,k_level
     1                     ,i_overlay,c_display
     1                     ,nx_l,ny_l,nz_l,max_radars
!    1                     ,grid_ra_ref_dum,grid_ra_vel_dum
     1                     ,nx_l,ny_l,r_missing_data,laps_cycle_time
     1                     ,jdot)

        elseif(c_type_i .eq. 's6')then ! 0-6km shear vector magnitude
            ext = 'lhe'
            var_2d = 'shu'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,u_2d,0,istatus)      
            var_2d = 'shv'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,v_2d,0,istatus)      

            if(istatus .ne. 1)then
                write(6,*)' error reading 0-6km shear vector'
                goto1200
            endif

            c_label = '0-6km agl wind shear     (kt)    '

            nxz = float(nx_l) / zoom
            nyz = float(ny_l) / zoom

            interval = int(max(nxz,nyz) / 65.) + 1
            size = float(interval) * .15

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            do i = 1,nx_l
            do j = 1,ny_l
                if(u_2d(i,j) .eq. r_missing_data
     1            .or. v_2d(i,j) .eq. r_missing_data)then
                    dir(i,j)  = r_missing_data
                    spds(i,j) = r_missing_data
                else
                    call uvgrid_to_disptrue(u_2d(i,j),
     1                              v_2d(i,j),
     1                              dir(i,j),
     1                              spds(i,j),
     1                              lat(i,j),
     1                              lon(i,j)     )
                    spds(i,j) = spds(i,j) / mspkt
                endif
            enddo ! j
            enddo ! i

            field_2d = spds ! support for diff option

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            c_label = '0-6km agl shear vector magnitude (kt)'
            scale = 1.
            clow = 0.
            chigh = 100.
            cint = 0.
            plot_parms%color_power = 1.0
            colortable = 'spectral'

            call plot_field_2d(i4time_pw,c_type,spds,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,colortable)

            if(i_image .eq. 1 .and. (.not. namelist_parms%l_sphere)
     1                                                         )then      
                plot_parms%icol_barbs = +1 ! keep future sfc barbs plots bright
            endif

        elseif(c_type .eq. 'bs')then ! surface backgrounds 
            write(6,711)
 711        format('   background extension [lgb,fsf]',5x,'? ',$)
            read(lun,712)ext
 712        format(a3)

            call input_background_info(
     1                              ext                     ! i
     1                             ,directory,c_model       ! o
     1                             ,i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,asc9_tim              ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o
            if(istatus.ne.1)goto1200

            if(ext.eq.'lgb')then
               write(6,723)
 723           format(/'  select field (var_2d-i):  '
     1          /
     1          /'     sfc: [usf,vsf,psf,tsf,dsf,fsf,slp,p] ? ',$)

            else
               write(6,725)
 725           format(/'  select field (var_2d-i):  '
     1          /
     1          /'  sfc: [usf,vsf,psf,tsf,dsf,rh,slp,th,the'       
     1                 ,',pbe,nbe,lhe,llr,lmr,lcv,s01,sto,'
     1                 /10x,'ptp,pdm,vnt,hah,ham,fwi,lwo,swo,tpw,li] ? '
     1                 ,$)       

            endif

            read(lun,724)var_2d
 724        format(a)
            call upcase(var_2d,var_2d)

            write(6,*)' variable selected is ',var_2d

            call s_len(var_2d,len_var)
            if(var_2d(len_var:len_var) .eq. 'i' .and. 
     1         var_2d(1:len_var) .ne. 'fwi'     .and.
     1         var_2d(1:len_var) .ne. 'li'                  )then
                l_image = .true.
                var_2d = var_2d(1:len_var-1)
            else
                l_image = .false.
            endif

!           if(var_2d .eq. 'lcv')then
!               l_image = .true.
!           endif

            level=0

            write(6,*)' var_2d is ',var_2d

            if(var_2d .ne. 'msf')then
                call read_laps(i4_initial,i4_valid,directory,ext
     1             ,nx_l,ny_l,1,1,var_2d,level,lvl_coord_2d
     1             ,units_2d,comment_2d,field_2d,istatus)

                if(istatus .ne. 1)then
                    write(6,*)' error reading grid ',var_2d,' ',ext
     1                                              ,istatus       
                    goto1200
                endif
            else
                call read_laps(i4_initial,i4_valid,directory,ext
     1             ,nx_l,ny_l,1,1,'dsf',level,lvl_coord_2d
     1             ,units_2d,comment_2d,td_2d,istatus)

                if(istatus .ne. 1)then
                    write(6,*)' error reading grid dsf ',ext,istatus       
                    goto1200
                endif

                call read_laps(i4_initial,i4_valid,directory,ext
     1             ,nx_l,ny_l,1,1,'psf',level,lvl_coord_2d
     1             ,units_2d,comment_2d,field_2d,istatus)

                if(istatus .ne. 1)then
                    write(6,*)' error reading grid psf ',ext,istatus  
                    goto1200
                endif

!               calculate mixing ratio from dewpoint and pressure
                do i = 1,nx_l
                do j = 1,ny_l
                    field_2d(i,j) = w( k_to_c(td_2d(i,j))
     1                                ,field_2d(i,j)/100. ) * 1e-3
                enddo ! j
                enddo ! i

            endif

            scale = 1. ! default value

            if(var_2d .eq. 'tsf' .or.
     1         var_2d .eq. 'dsf' .or.
     1         var_2d .eq. 'tgd' .or.
     1         var_2d .eq. 't'   .or.
     1         var_2d .eq. 'td'       )then

                write(6,726)
 726            format(10x,'plot fahrenheit or celsius [f/c]  ? ',$)
                read(5,*)tunits

                write(6,*)' converting sfc data to f/c'

!               kelvin conversion to f or c
                do i = 1,nx_l
                do j = 1,ny_l
                    if(tunits .ne. 'c' .or. l_image)then
                        field_2d(i,j) = k_to_f(field_2d(i,j))
                        units_2d = 'deg f'
                    else
                        field_2d(i,j) = k_to_c(field_2d(i,j))
                        units_2d = 'deg c'
                    endif
                enddo ! j
                enddo ! i

            elseif(var_2d .eq. 'ps'  .or. var_2d .eq. 'psf'
     1        .or. var_2d .eq. 'msl' .or. var_2d .eq. 'slp'
     1        .or. var_2d .eq. 'p'                         )then
                scale = 100.
                units_2d = 'hpa'

            elseif(var_2d(2:3) .eq. '01' .or. var_2d(2:3) .eq. 'to')then       
                if(c_units_type .eq. 'english')then
                    scale = 1. / ((100./2.54)) ! denom = (in/m)
                    units_2d = 'in'
                else ! metric
                    scale = .001               ! numer = (m/mm)
                    units_2d = 'mm'
                endif

            elseif(var_2d .eq. 'tpw')then       
                scale = 1e-2    ! data are in m, display is cm
                units_2d = 'cm'

            elseif(var_2d .eq. 'umf')then       
                scale = 1e-2    
                units_2d = 'cm-m/s'
                plot_parms%iraster = 1

            elseif(var_2d .eq. 'vis')then       
                scale = 1000.    
                units_2d = 'km'
                plot_parms%iraster = 1

            elseif(var_2d .eq. 'msf')then       
                scale = 1e-3  ! calculated field 
                units_2d = 'g/kg'
                comment_2d = 'sfc mixing ratio'

            elseif(var_2d .eq. 'rsf')then       
                scale = 1e-3                        
                units_2d = 'g/kg'
                comment_2d = 'sfc spec humidity'

            elseif(var_2d .eq. 'pdm')then 
                if(namelist_parms%c_pbl_depth_units .eq. 'english')then       
                    units_2d = 'ft'
                    scale = 1. / ft_per_m
                endif

            elseif(var_2d .eq. 'vnt')then 
                if(c_vnt_units .eq. 'kt-ft')then
                    scale = 1000. / (ft_per_m / mspkt) ! convert from m**2/s 
                                                       ! to kt-ft (inverse)
                    units_2d = 'kt-ft x1000'
                endif

            elseif(var_2d .eq. 'lwo')then ! temporary fix for fsf inconsistency
!               comment_2d = 'brightness temperature'
!               units_2d = 'deg k'

                comment_2d = 'olr'
                units_2d = 'w/m**2'

!               convert olr radiance to brightness temperature
!               write(6,*)' converting radiance to brightness temp'
!               do i = 1,nx_l
!               do j = 1,ny_l
!                   if(field_2d(i,j) .ne. r_missing_data)then
!                       field_2d(i,j) = rad_to_temp(field_2d(i,j))       
!                   endif
!               enddo ! j
!               enddo ! i
!               where(field_2d(:,:) .ne. r_missing_data)
!    !                field_2d(:,:) = rad_to_temp(field_2d(:,:))

            elseif(units_2d(1:4) .eq. 'none' .or.
     1             units_2d(1:4) .eq. 'none'      )then
                units_2d = '          '

            endif

            call s_len2(comment_2d,len_fcst)
            call s_len2(units_2d,len_units)
            call s_len(c_model,len_model)

            call upcase(c_model,c_model)

            if(ext(1:3) .eq. 'fsf')then
                len_fcst = min(len_fcst,30)

                call mk_fcst_hlabel(0,comment_2d(1:len_fcst),fcst_hhmm       
     1                                 ,ext(1:3),units_2d(1:len_units)
     1                                 ,c_model,c_label)

            else ! lgb
                len_fcst = 25 
                call mk_fcst_hlabel(0,comment_2d(1:len_fcst),fcst_hhmm       
     1                                 ,ext(1:3),units_2d(1:len_units)
     1                                 ,c_model,c_label)

            endif

            write(6,*)' c_label = ',c_label

            call make_fnam_lp(i4_valid,asc9_tim,istatus)

            write(6,*)' l_image = ',l_image

            if(.not. l_image)then ! surface background contours
                if(var_2d .eq. 's01')then
                    clow = 0.0
                    chigh = 100.0
                    cint = -.1
                else
                    call contour_settings(field_2d,nx_l,ny_l
     1                            ,clow,chigh,cint,zoom,density,scale)      
                endif

                call plot_cont(field_2d,scale,clow,chigh,cint
     1                        ,asc9_tim,namelist_parms,plot_parms
     1                        ,c_label,i_overlay,c_display
     1                        ,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,laps_cycle_time)       
            
            else                  ! surface background / forecast images
                if(var_2d .eq. 'llr' .or. var_2d .eq. 'lmr')then
                    call ccpfil(field_2d,nx_l,ny_l,-10.0,70.0,'ref'
     1                         ,n_image,scale,'hsect',plot_parms
     1                         ,namelist_parms)        
                elseif(var_2d .eq. 'tsf' .or. var_2d .eq. 'tgd')then
                    call ccpfil(field_2d,nx_l,ny_l,sfctf_l,sfctf_h      
     1                         ,'temp',n_image,scale,'hsect'
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'dsf')then
                    call ccpfil(field_2d,nx_l,ny_l,sfctdf_h,sfctdf_l
     1                         ,'hues',n_image,scale,'hsect'
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'slp')then
                    call ccpfil(field_2d,nx_l,ny_l,1040.,960.
     1                         ,'spectral',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 's8a')then
                    call ccpfil(field_2d,nx_l,ny_l
     1                         ,c_to_k(btemp_h),c_to_k(btemp_l)
     1                         ,btemp_table,n_image,scale,'hsect'
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'lwo')then
!                   call ccpfil(field_2d,nx_l,ny_l,313.15,223.15
                    call ccpfil(field_2d,nx_l,ny_l,313.15,100.15
     1                         ,'linear',n_image,scale,'hsect'
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'swo')then
                    call ccpfil(field_2d,nx_l,ny_l,0.0,500.
     1                         ,'linear',n_image,scale,'hsect'
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'swi')then
                    call ccpfil(field_2d,nx_l,ny_l,0.0,1000.
     1                         ,'spectral',n_image,scale,'hsect'
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'vis')then
                    call ccpfil(field_2d,nx_l,ny_l,20.,0.
     1                         ,'linear',n_image,scale,'hsect'
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'pbe')then
!                   call condition_cape(nx_l,ny_l,'pei',r_missing_data
!    1                                 ,field_2d)
                    call ccpfil(field_2d,nx_l,ny_l,0.0,chigh_cape
     1                         ,'cpe',n_image,scale,'hsect'
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'nbe')then
                    call ccpfil(field_2d,nx_l,ny_l,-500.,+50.
     1                         ,'cpe',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'lhe')then
                    call ccpfil(field_2d,nx_l,ny_l,hel_l,hel_h
     1                         ,'spectral',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'the')then
                    call ccpfil(field_2d,nx_l,ny_l,250.,370.
     1                         ,'hues',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'rh')then
                    call ccpfil(field_2d,nx_l,ny_l,0.,100.
     1                         ,'moist',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'tpw')then
                    plot_parms%color_power = namelist_parms%power_tpw
                    call ccpfil(field_2d,nx_l,ny_l,0.
     1                         ,namelist_parms%chigh_tpw
     1                         ,'tpw',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'msf' .or. var_2d .eq. 'rsf')then
                    plot_parms%color_power = 0.7
                    call ccpfil(field_2d,nx_l,ny_l,0.,25.0
     1                         ,'tpw',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'umf')then
                    call ccpfil(field_2d,nx_l,ny_l,umf_l,umf_h
     1                         ,'upflux',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'vnt')then
                    if(c_vnt_units .eq. 'kt-ft')then
                        call ccpfil(field_2d,nx_l,ny_l,150.,0.
     1                             ,'vnt',n_image,scale,'hsect' 
     1                             ,plot_parms,namelist_parms) 
                    else
                        call ccpfil(field_2d,nx_l,ny_l,5000.,0.
     1                             ,'vnt',n_image,scale,'hsect' 
     1                             ,plot_parms,namelist_parms) 
                    endif
                    plot_parms%icol_barbs = +1 ! keep future barbs plots bright
                elseif(var_2d(1:2) .eq. 'ha')then
                    call ccpfil(field_2d,nx_l,ny_l,2.,6.
     1                         ,'haines',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'fwi')then
                    call ccpfil(field_2d,nx_l,ny_l,0.,40.
     1                         ,'spectralr',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'fwx')then
                    call ccpfil(field_2d,nx_l,ny_l,0.,20.
     1                         ,'spectralr',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'ptp')then
                    call ccpfil(field_2d,nx_l,ny_l,50000.,110000.
     1                         ,'spectral',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d .eq. 'pdm' .or. var_2d .eq. 'blh')then       
                    if(c_units_type .eq. 'metric')then
                        chigh = 2400.
                    else ! english
                        chigh = 8000.
                    endif

                    call ccpfil(field_2d,nx_l,ny_l,0.,chigh
     1                         ,'spectral',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                elseif(var_2d(2:3) .eq. '01' .or.           ! precip
     1                 var_2d(2:3) .eq. 'to')then
                    if(var_2d(1:1) .eq. 'r')then
                        chigh = +10.
                        if(var_2d(2:3) .eq. '01')then
                            colortable = 'acc_inc'
                        else
                            colortable = 'acc_sto'
                        endif
                    else
                        chigh = +20.
                        if(var_2d(2:3) .eq. '01')then
                            colortable = 'sno_inc'
                        else
                            colortable = 'sno_sto'
                        endif
                    endif

                    call condition_precip(nx_l,ny_l,'pai',field_2d
     1                                   ,scale,.01)      

                    call ccpfil(field_2d,nx_l,ny_l,0.,chigh ! *scale
     1                         ,colortable,n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 

                elseif(var_2d .eq. 'lcv' .or. var_2d .eq. 'cla' .or. 
     1                 var_2d .eq. 'smv')then
                    call ccpfil(field_2d,nx_l,ny_l,0.0,cloud_albedo_f
     1                         ,'linear',n_image,scale,'hsect'
     1                         ,plot_parms,namelist_parms)        
!               elseif(var_2d .eq. 'lil')then
!                   call ccpfil(field_2d,nx_l,ny_l,0.0,2.0,'moist'
!    1                         ,n_image,scale,'hsect',plot_parms
!    1                         ,namelist_parms)        
                elseif(var_2d .eq. 'lhf')then
                    call ccpfil(field_2d,nx_l,ny_l,-100.,+500.
     1                         ,'spectral',n_image,scale,'hsect' 
     1                         ,plot_parms,namelist_parms) 
                else
!                   call array_range(field_2d,nx_l,ny_l,rmin,rmax
!    1                              ,r_missing_data)
!                   call ccpfil(field_2d,nx_l,ny_l
!    1                         ,rmin/scale,rmax/scale,'spectral'       
!    1                         ,n_image,scale,'hsect',plot_parms
!    1                         ,namelist_parms)        

                    call contour_settings(field_2d,nx_l,ny_l
     1                           ,clow,chigh,cint
     1                           ,zoom,density,scale)       

                    call ccpfil(field_2d,nx_l,ny_l
     1                         ,clow,chigh,'spectral'       
     1                         ,n_image,scale,'hsect',plot_parms
     1                         ,namelist_parms)        
                endif

                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum('in',7)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)

            endif

        elseif(c_type_i .eq. 'p' 
     1    .or. c_type_i .eq. 'pm' 
     1    .or. c_type_i .eq. 'pp' )then ! reduced/msl/pert pres

            if(c_type_i .eq. 'p')then
                var_2d = 'p'
            elseif(c_type_i .eq. 'pm')then
                var_2d = 'msl'
            elseif(c_type_i .eq. 'pp')then
                var_2d = 'pp'
            endif

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,field_2d,0
     1                          ,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            if(c_type_i .eq. 'p')then
                c_label = comment_2d(1:23)//' (mb)     '       
            elseif(c_type_i .eq. 'pm')then
                c_label = 'msl pressure            (mb)     '
            elseif(c_type_i .eq. 'pp')then
                c_label = 'perturbation pressure   (mb)     '
            endif

!           add balance to label
            call s_len2(c_label,len_c_label)
            if(i_balance .eq. 1)then
                c_label = c_label(1:len_c_label)//' (bal)'
            endif

            clow = 0.
            chigh = 0.
!           cint = 1. ! 3.

            scale = 100.

            call array_range(field_2d,nx_l,ny_l,pres_low_pa,pres_high_pa
     1                      ,r_missing_data)

            if(c_type_i .eq. 'pp' .and. i_image .eq. 1)then ! p_prime image
                clow = +5.0
                chigh = -5.0
                cint = 0.1

                plot_parms%ncols = nint(abs(chigh-clow)/cint)
 
            else ! default pressure plot
                pres_high_mb = pres_high_pa / scale
                pres_low_mb  = pres_low_pa / scale
                range = ((pres_high_pa - pres_low_pa) / scale) 
     1                                    / (sqrt(zoom) * density)

                if(range .gt. 35.)then
                    cint = 4.
                elseif(range .gt. 8.)then
                    cint = 2.
                else ! range < 8
                    cint = 1.
                endif

                icint   = cint
                clow = int(pres_low_mb) / 4 * 4
                chigh = (int(pres_high_mb) / icint) * icint + icint
            endif

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

        elseif(c_type_i .eq. 'ps' .or. c_type_i .eq. 'al')then ! surface pres
            var_2d = 'ps'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            if(c_type_i .eq. 'ps')then
                c_label = 'surface pressure        (mb)     '
            else ! convert to altimeter setting
                c_label = 'altimeter setting       (mb)     '
                do i = 1,nx_l
                do j = 1,ny_l
!                   find standard atmosphere p at this elevation
                    pstd_pa = ztopsa(topo(i,j)) * 100.
                    pstn_pa = field_2d(i,j)
                    alt_pa  = psamslpa * (pstn_pa / pstd_pa)
                    field_2d(i,j) = alt_pa
                enddo ! j
                enddo ! i
            endif

!           add balance to label
            call s_len2(c_label,len_c_label)
            if(i_balance .eq. 1)then
                c_label = c_label(1:len_c_label)//' (bal)'
            endif

            scale = 100.

            call contour_settings(field_2d,nx_l,ny_l,clow,chigh,cint
     1                           ,zoom,density,scale)

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

!           call plot_cont(field_2d,scale,clow,chigh,cint
!    1             ,asc9_tim,namelist_parms,plot_parms
!    1             ,c_label,i_overlay,c_display
!    1             ,lat,lon,jdot
!    1             ,nx_l,ny_l,r_missing_data,laps_cycle_time)

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'spectral')

        elseif(c_type .eq. 'vv')then
            var_2d = 'vv'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1          ,i4time_pw,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                     ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'sfc vert velocity     (cm/s)     '

            clow = -200.
            chigh = +200.
            cint = 10.

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_cont(field_2d,1e-2,clow,chigh,cint
     1                    ,asc9_tim,namelist_parms,plot_parms
     1                    ,c_label,i_overlay,c_display,lat,lon,jdot
     1                    ,nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif(c_type(1:2) .eq. 'hu')then
            var_2d = 'rh'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,
     1              ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                     ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'sfc rel humidity   (%)     '

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            clow = 0.
            chigh = 100.
            cint = 10.


            call plot_field_2d(i4time_pw,c_type,field_2d,1e0
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'moist')

        elseif(c_type(1:2) .eq. 'sm')then
            var_2d = 'lsm'
            ext = 'lm1'

!           input level of soil moisture
            write(6,735)
 735        format(35x,'input soil level           [1,2,3] ? ',$)
            read(5,*)lvl_soil

            write(clvl_soil,736)lvl_soil
 736        format(i1)

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,-lvl_soil,istatus)

            if(istatus .ne. 1 .and. istatus .ne. -1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'soil moisture level '//clvl_soil
     1                                      //' (volumetric) '

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            clow = 0.
            chigh = 0.5
            cint = 0.05

            plot_parms%ncols = nint(abs(chigh-clow)/cint)
            plot_parms%iraster = 1

            call plot_field_2d(i4time_pw,c_type,field_2d,1e0
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'moist')

        elseif(c_type .eq. 'ta')then
            var_2d = 'tad'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,field_2d
     1                          ,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'sfc temp adv (x 1e-5 dg k/s)     '

            clow = -100.
            chigh = +100.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_cont(field_2d,1e-5,clow,chigh,cint,
     1          asc9_tim,namelist_parms,plot_parms,
     1          c_label,i_overlay,c_display,lat,lon,jdot,
     1          nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif(c_type(1:2) .eq. 'th')then
            var_2d = 'th'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1             ,i4time_pw,
     1              ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                     ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'sfc potential temp   (deg k)     '

            clow = +240.
            chigh = +330.
            cint = 2.

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

!           call plot_cont(field_2d,1e-0,clow,chigh,cint,asc9_tim
!    1                    ,namelist_parms,plot_parms,c_label
!    1                    ,i_overlay,c_display,lat,lon,jdot
!    1                    ,nx_l,ny_l,r_missing_data,laps_cycle_time)

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

        elseif(c_type(1:2) .eq. 'te')then
            var_2d = 'the'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'sfc equiv potl temp  (deg k)     '

            clow = +250.
            chigh = +370.
            cint = 2.
            scale = 1e0

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

        elseif(c_type(1:2) .eq. 'vo')then
            var_2d = 'vor'

            if(i_balance .eq. 1)then
                ext = 'balance'
                c_label = 'sfc rel vorticity (balanced)  (1e-5/s)'
            else
                ext = 'lsx'
                c_label = 'sfc rel vorticity   (1e-5/s)'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw
     1                          ,ext,var_2d,units_2d,comment_2d
     1                          ,nx_l,ny_l,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            if(grid_spacing_m .le. 1500)then
                clow = -160.
                chigh = +160.
                cint = 10.
            else
                clow = -80.
                chigh = +80.
                cint = 5.
            endif

            scale = 1e-5

            plot_parms%iraster = 1

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

        elseif(c_type(1:2) .eq. 'mr')then
            var_2d = 'mr'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            field_2d = field_2d * .001

            c_label = 'sfc mixing ratio      (g/kg)     '

            clow = 0.
            chigh = +25.
            cint = 1.

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

!           call plot_cont(field_2d,1e-0,clow,chigh,cint,
!    1        asc9_tim,namelist_parms,plot_parms,
!    1        c_label,i_overlay,c_display,lat,lon,jdot,
!    1        nx_l,ny_l,r_missing_data,laps_cycle_time)

            plot_parms%color_power = 0.7

            call plot_field_2d(i4time_pw,c_type,field_2d,1e-3
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'tpw')

        elseif(c_type(1:2) .eq. 'dv')then
            var_2d = 'div'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'sfc divergence  (x 1e-5 s-1)     '

            if(grid_spacing_m .le. 1500.)then
                chigh = +80.
                clow = -80.
                cint = 20.
            else
                chigh = +40.
                clow = -40.
                cint = 10.
            endif

            scale = 1e-5
!           call contour_settings(field_2d,nx_l,ny_l,clow,chigh,cint       
!    1                           ,zoom,density,scale)
            plot_parms%iraster = 1

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,chigh,clow,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'spectralr')

        elseif(c_type .eq. 'ha')then ! theta advection
            var_2d = 'tha'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                                     ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'sfc theta adv   (x 1e-5 k/s)     '

            clow = -100.
            chigh = +100.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_cont(field_2d,1e-5,clow,chigh,cint,
     1        asc9_tim,namelist_parms,plot_parms,
     1        c_label,i_overlay,c_display,lat,lon,jdot,
     1        nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'ma')then
            var_2d = 'mra'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                                     ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'sfc mstr adv (x 1e-5 g/kg/s)     '

            clow = -100.
            chigh = +100.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_cont(field_2d,1e-5,clow,chigh,cint,
     1        asc9_tim,namelist_parms,plot_parms,
     1        c_label,i_overlay,c_display,lat,lon,jdot,
     1        nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif(c_type(1:2) .eq. 'sp')then

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            var_2d = 'u'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,u_2d,0,istatus)
            var_2d = 'v'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,v_2d,0,istatus)       


            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            plot_parms%l_discrete = .false.

            do i = 1,nx_l
            do j = 1,ny_l
                    if(u_2d(i,j) .eq. r_missing_data
     1            .or. v_2d(i,j) .eq. r_missing_data)then
                        dir(i,j)  = r_missing_data
                        spds(i,j) = r_missing_data
                    else
                        call uvgrid_to_disptrue(u_2d(i,j),
     1                                  v_2d(i,j),
     1                                  dir(i,j),
     1                                  spds(i,j),
     1                                  lat(i,j),
     1                                  lon(i,j)     )
                        spds(i,j) = spds(i,j) / mspkt
                        if(c_type(1:3) .eq. 'sp3')then
                            spds(i,j) = spds(i,j)**3
                        endif
                    endif
            enddo ! j
            enddo ! i

            field_2d = spds ! support for diff option

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            if(c_type(1:3) .eq. 'sp3')then
                c_label = 'sfc wind power       (1000kt**3) '
                scale = 1000.
                clow = 0.
                chigh = (chigh_sfcwind**3 / scale) * 0.2
                cint = 0.
                plot_parms%color_power = 0.4  
                colortable = 'power' 
            else
                c_label = 'sfc wind speed          (kt)     '
                scale = 1.
                clow = 0.
                chigh = chigh_sfcwind
                cint = 0.
                plot_parms%color_power = 1.0
                colortable = 'spectral'
            endif

            call plot_field_2d(i4time_pw,c_type,spds,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,colortable)

            if(i_image .eq. 1 .and. (.not. namelist_parms%l_sphere)
     1                                                         )then      
                plot_parms%icol_barbs = +1 ! keep future sfc barbs plots bright
            endif

        elseif(c_type_i .eq. 'u' .or. c_type_i .eq. 'v')then
            call upcase(c_type,c_type)

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            var_2d = c_type(1:1)
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l,field_2d,0
     1                          ,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'sfc wind '//c_type(1:1)//'       (m/s)        '       

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            if(i_image .eq. 1)then
                clow = -20.
                chigh = +20.
                cint = 5.
            else
                call contour_settings(field_2d,nx_l,ny_l
     1                           ,clow,chigh,cint,zoom,density,1.)       
            endif

            scale = 1.
            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

        elseif(c_type .eq. 'cs')then
            var_2d = 'css'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'colorado severe storm index     '

            clow = 0.
            chigh = +100.
            cint = 10.

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_cont(field_2d,1e-0,clow,chigh,cint,
     1        asc9_tim,namelist_parms,plot_parms,
     1        c_label,i_overlay,c_display,lat,lon,jdot,
     1        nx_l,ny_l,r_missing_data,laps_cycle_time)

        elseif(c_type_i .eq. 'vs')then
            var_2d = 'vis'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lil'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)
            if(istatus .ne. 1)then
                write(6,*)' error reading lil visibility, try lsx '
     1                                                      ,var_2d
                ext = 'lsx'
                call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)
                if(istatus .ne. 1)then
                    goto1200
                else
                    c_label = 'sfc visibility  (miles - lsx file) '
                endif
            else
                c_label = 'sfc visibility  (miles - lil file) '
            endif

            clow = 30.
            chigh = +0.
            cint = -0.1
            scale = 1600.
            plot_parms%color_power = 2.0 

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'linear')

        elseif(c_type(1:2) .eq. 'fw')then
            var_2d = 'fwx'

            if(i_balance .eq. 1)then
                ext = 'balance'
            else
                ext = 'lsx'
            endif

            call get_laps_2dgrid(
     1               i4time_ref,laps_cycle_time*100,i4time_pw,
     1               ext,var_2d,units_2d,comment_2d,nx_l,ny_l,
     1               field_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading surface ',var_2d
                goto1200
            endif

            c_label = 'fire weather          (0-20)     '

            call make_fnam_lp(i4time_pw,asc9_tim,istatus)

            if(c_type(3:3) .ne. 'i')then
                clow = 0.
                chigh = +20.
                cint = 2.0

                call plot_cont(field_2d,1.,clow,chigh,cint,
     1           asc9_tim,namelist_parms,plot_parms,
     1           c_label,i_overlay,c_display,lat,lon,jdot,       
     1           nx_l,ny_l,r_missing_data,laps_cycle_time)
            else
                call ccpfil(field_2d,nx_l,ny_l,0.0,20.0,'spectral'
     1                     ,n_image,1.,'hsect',plot_parms
     1                     ,namelist_parms)      
                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum('in',7)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)

            endif

        elseif(c_type(1:2) .eq. 'sc' .or. c_type(1:3) .eq. 'csc')then       
            if(c_type(1:2) .eq. 'sc')then
                var_2d = 'sc'
                ext = 'lm2'
                c_label = 'snow cover analysis    (%) '
            else
                var_2d = 'csc'
                ext = 'lcv'
                c_label = 'satobs (csc) snow cover (%)'
            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)

            if(istatus .ne. 1 .and. istatus .ne. -1)then
                write(6,*)' error reading snow cover ',ext(1:3)
     1                   ,' ',var_2d
                goto1200
            endif

            clow = 0.
            chigh = +100.
            cint = 25.
            scale=0.01

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'linear')

        elseif(c_type(1:2) .eq. 'sd')then       
            var_2d = 'sc'
            ext = 'lm2'
            c_label = 'snow depth analysis    (in) '

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)

            if(istatus .ne. 1 .and. istatus .ne. -1)then
                write(6,*)' error reading snow depth ',ext(1:3)
     1                   ,' ',var_2d
                goto1200
            endif

            clow = 0.
            chigh = +20.
            cint = 1.

            scale = 1. / (100./2.54) ! denominator = (in/m)

            call plot_field_2d(i4time_pw,c_type,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'linear')

        elseif(c_type(1:2) .eq. 'cb' .or. c_type(1:2) .eq. 'ct' 
     1                               .or. c_type(1:2) .eq. 'cc')then
            ext = 'lcb'

            if(c_type(1:2) .eq. 'cb')then ! cloud base
                var_2d = 'lcb'
                c_label = 'cloud base         m   msl       '
!               c_label = 'cloud prlx top  (smoothed)       '
                clow = 0.
                chigh = 10000.
                chigh_img = 8000. ! 8000. or 16000.
                cint = 1000.

            elseif(c_type(1:2) .eq. 'ct')then
                var_2d = 'lct'
                c_label = 'cloud top          m   msl       '
                clow = 2000.
                chigh = 20000.
                chigh_img = 16000.
                cint = 1000.

            elseif(c_type(1:2) .eq. 'cc')then ! cloud ceiling
                var_2d = 'cce'
!               c_label = 'cloud ceiling      m   agl       '
                c_label = 'cloud prlx top (unsmoothed)      '
                clow = 0.
                chigh = 0.
                chigh_img = 16000. ! 8000.
                cint = -100.

            endif

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_nearest,ext,var_2d
     1                          ,units_2d,comment_2d,nx_l,ny_l
     1                          ,field_2d,0,istatus)

            if(istatus .ne. 1 .and. istatus .ne. -1)then
                write(6,*)' error reading ',ext,var_2d
                goto1200
            endif


            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

            if(c_type(3:3) .ne. 'i')then
                call plot_cont(field_2d,1e0,clow,chigh,cint,
     1               asc9_tim,namelist_parms,plot_parms,c_label,       
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)

            else
                where(field_2d .eq. r_missing_data)field_2d = 0.
                call ccpfil(field_2d,nx_l,ny_l,clow,chigh_img,'linear'
     1                     ,n_image,1e0,'hsect',plot_parms
     1                     ,namelist_parms)       
                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum('in',7)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)

            endif

        elseif(c_type .eq. 'cv' .or. c_type .eq. 'cg')then
            write(6,2514)
2514        format('     enter level (1-42); [-bbbb] for mb; '
     1                          ,'or [0] for max in column',5x,'? ',$)
            read(lun,*)k_level

            if(k_level .gt. 0)then
                var_2d = 'lc3'
                ext = 'lc3'

                call get_laps_2dgrid(i4time_ref,86400,i4time_nearest,
     1           ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                  ,cloud_cvr,k_level,istatus)

            elseif(k_level .lt. 0)then ! k_level is -pressure in mb
               var_2d = 'lcp'
               ext = 'lcp'

               call get_laps_2dgrid(i4time_ref,86400,i4time_nearest,
     1           ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                  ,cloud_cvr,-k_level,istatus)

            else ! k_level .eq. 0
               var_2d = 'lcv'
               ext = 'lcv'

               call get_laps_2dgrid(i4time_ref,86400,i4time_nearest,
     1           ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                  ,cloud_cvr,-k_level,istatus)

            endif

            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

!           get cloud cover
            if(k_level .gt. 0)then
                read(comment_2d,3515)cloud_height
3515            format(e20.8)

                write(c_label,3516)nint(cloud_height)
3516            format(i5,'  m msl   cloud cover       ')

                write(6,*)' lvl_cld = ',lvl_cld

            elseif(k_level .eq. 0)then
                c_label = 'cloud cover (fraction)           '

            else ! k_level .lt. 0
                write(c_label,3517)-k_level
3517            format(i5,'  mb    cloud cover         ')

            endif

            if(c_type .eq. 'cv')then
                clow = 0.1
                chigh = 0.9
                cint = 0.1
                call plot_cont(cloud_cvr,1e0,clow,chigh,cint,
     1               asc9_tim,namelist_parms,plot_parms,c_label,        
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)

            else ! 'cg'
                write(6,*)' calling solid fill cloud plot'
                
                if(nx_l*ny_l .gt. 1000000)then
                    colortable = 'linear_reduced'
                else
                    colortable = 'linear'
                endif

                plot_parms%icol_barbs = +1 ! keep future barbs plots bright

                call ccpfil(cloud_cvr,nx_l,ny_l,0.0,1.0,colortable     
     1                     ,n_image,1e0,'hsect',plot_parms
     1                     ,namelist_parms)     
                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum('in',7)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)

            endif

            field_2d = cloud_cvr ! supports diff option

        elseif(c_type(1:2) .eq. 'gg')then

           if(.not. allocated (static_grid))then
              allocate (static_grid(nx_l,ny_l))
           endif

           call make_fnam_lp(i4time_ref,asc9_tim,istatus)

           write(6,219)
219        format(5x,'select static field:'
     1/' [tn-i,lf-i,  gr,  lu, al-i,ts-i,gn-i,sn-i,sl-i]'
     1/'  ter/lndfrac/grid/use/alb/temp/green/slp ? ',$)
           read(lun,*)cstatic

           if(cstatic(1:2).eq.'tn')then
              clow = -400.
              chigh = +5000.
              cint = +200.
              c_label = 'static terrain (m)'

              topomax = -9999.
              do i = 1,nx_l
              do j = 1,ny_l
                  topomax = max(topomax,topo(i,j))
              enddo ! j
              enddo ! i

              field_2d(:,:) = topo(:,:) ! support for difference plots

              if(cstatic .eq. 'tni')then
                write(6,*)' calling solid fill plot - max terrain = '       
     1                   ,topomax
                scale = 4000.
                call ccpfil(topo,nx_l,ny_l,0.0,scale,'linear',n_image
     1                     ,1e0,'hsect',plot_parms,namelist_parms)     
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)
              else
                scale = 1.
                call contour_settings(topo,nx_l,ny_l
     1               ,clow,chigh,cint,zoom,density,scale)      

                call plot_cont(topo,scale,
     1               clow,chigh,cint,asc9_tim,
     1               namelist_parms,plot_parms,c_label,      
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)
              endif
              i4time_topo = 0

           elseif(cstatic(1:2) .eq. 'gr')then
               call plot_grid(i_overlay,c_display,lat,lon,
     1                     nx_l,ny_l,laps_cycle_time)

           elseif(cstatic(1:2) .eq. 'lf')then
              var_2d='ldf'
              call read_static_grid(nx_l,ny_l,var_2d,static_grid
     1                             ,istatus)
              if(istatus .ne. 1)then
                 print*,' warning: could not read laps static-slope-lat'
              endif

              clow = .01
              chigh = .99
              cint = .49
              c_label = 'static land fraction         '
              if(cstatic .eq. 'lfi')then
                write(6,*)' calling solid fill plot'
                scale = 1.
                call ccpfil(static_grid,nx_l,ny_l,0.0,scale,'linear'
     1                     ,n_image,1e0,'hsect',plot_parms
     1                     ,namelist_parms)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)
              else

                call plot_cont(static_grid,1e0,
     1               clow,chigh,cint,asc9_tim,
     1               namelist_parms,plot_parms,c_label,      
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)

                i4time_topo = 0
              endif
           elseif(cstatic(1:2) .eq. 'lu')then
              var_2d='use'
              call read_static_grid(nx_l,ny_l,var_2d,static_grid
     1,istatus)
              if(istatus .ne. 1)then
                 print*,' warning: could not read static-landuse'
                 return
              endif
              clow = 0.
              chigh = 24.
              cint = 1.
              c_label = 'static land use           '
!             call plot_cont(static_grid,1e0,
!    1               clow,chigh,cint,asc9_tim,
!    1               namelist_parms,plot_parms,c_label,      
!    1               i_overlay,c_display,lat,lon,jdot,
!    1               nx_l,ny_l,r_missing_data,laps_cycle_time)
              i4time_topo = 0

              plot_parms%iraster = 1

              call plot_field_2d(i4time_topo,c_field,static_grid,1e0       
     1                        ,namelist_parms,plot_parms
     1                        ,clow,chigh,cint,c_label
     1                        ,i_overlay,c_display,lat,lon,jdot
     1                        ,nx_l,ny_l,r_missing_data,'hues')

           elseif(cstatic(1:2) .eq. 'al')then

              call get_static_field_interp('albedo',i4time_ref
     1,nx_l,ny_l,static_grid,istatus)
              if(istatus .ne. 1)then
                 print*,' warning: could not read static-albedo'
                 return
              endif

              call get_mxmn_2d(nx_l,ny_l,static_grid,chigh
     1                        ,clow,imx,jmx,imn,jmn)
              cint = (chigh-clow)/10.
              if(cint.eq.0.0)cint=.05

c             clow = 0.
c             chigh = 1.0
c             cint = .05

              c_label = 'interpolated albedo     '
              asc9_tim=asc9_tim(1:5)//'1800'

              if(cstatic .eq. 'ali' .or. i_image .eq. 1)then
                write(6,*)' calling solid fill plot'
                call ccpfil(static_grid,nx_l,ny_l,clow,chigh,'linear'
     1                     ,n_image,1e0,'hsect',plot_parms
     1                     ,namelist_parms)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)
              else
                call plot_cont(static_grid,1e0,
     1               clow,chigh,cint,asc9_tim,
     1               namelist_parms,plot_parms,c_label,
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)
              endif
              call move(static_grid,field_2d,nx_l,ny_l) ! supports diff option
              i4time_topo = 0

           elseif(cstatic(1:2) .eq. 'gn')then

              call get_static_field_interp('green',i4time_ref
     1,nx_l,ny_l,static_grid,istatus)
              if(istatus .ne. 1)then
                 print*,' warning: could not read static-albedo'
                 return
              endif

              static_grid=static_grid/100.

c             call get_mxmn_2d(nx_l,ny_l,static_grid,chigh
c    1                        ,clow,imx,jmx,imn,jmn)
c             cint = (chigh-clow)/10.
c             if(cint.eq.0.0)cint=0.1

              clow = 0.0
              chigh = 1.0

              c_label = 'interpolated green fraction     '
              call make_fnam_lp(i4time_ref,asc9_tim,istatus)
              asc9_tim=asc9_tim(1:5)//'1800'

              if(cstatic .eq. 'gni')then

                write(6,*)' calling solid fill plot'
                call ccpfil(static_grid,nx_l,ny_l,clow,chigh,'green'
     1                     ,n_image,1e0,'hsect',plot_parms
     1                     ,namelist_parms)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)
              else
                cint = .10
                call plot_cont(static_grid,1e0,
     1               clow,chigh,cint,asc9_tim,
     1               namelist_parms,plot_parms,c_label,
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)
              endif

           elseif(cstatic(1:2) .eq. 'ts')then

              var_2d='tmp'
              call read_static_grid(nx_l,ny_l,var_2d,static_grid
     1,istatus)

              if(istatus .ne. 1)then
                 print*,' warning: could not read static: soil temp'
                 return
              endif

              call get_mxmn_2d(nx_l,ny_l,static_grid,chigh
     1                        ,clow,imx,jmx,imn,jmn)

              cint = (chigh-clow)/10.
              c_label = 'mean annual soil temp (k)  '

              if(cstatic .eq. 'tsi')then
                write(6,*)' calling solid fill plot'
                call ccpfil(static_grid,nx_l,ny_l,clow,chigh
     1               ,'hues',n_image,1e0,'hsect',plot_parms
     1               ,namelist_parms)
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)
              else
                call plot_cont(static_grid,1e0,
     1               clow,chigh,cint,asc9_tim,
     1               namelist_parms,plot_parms,c_label,
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)
              endif
              i4time_topo = 0

           elseif(cstatic(1:2) .eq. 'sn')then

              var_2d='sln'
              call read_static_grid(nx_l,ny_l,var_2d,static_grid
     1,istatus)
              if(istatus .ne. 1)then
                 print*,' warning: could not read static-slope-lon'
                 return
              endif
              clow = -1.
              chigh = 1.0
              cint = .025
              c_label = 'longitude component terrain slope '

              if(cstatic .eq. 'sni')then
                write(6,*)' calling solid fill plot'
                call get_mxmn_2d(nx_l,ny_l,static_grid,rmx2d
     1                          ,rmn2d,imx,jmx,imn,jmn)
                call ccpfil(static_grid,nx_l,ny_l,rmn2d,rmx2d
     1                     ,'linear',n_image,1e0,'hsect',plot_parms
     1                     ,namelist_parms)       
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)
              else
                call plot_cont(static_grid,1e0,
     1               clow,chigh,cint,asc9_tim,
     1               namelist_parms,plot_parms,c_label,
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)
              endif
              i4time_topo = 0

           elseif(cstatic(1:2) .eq. 'sl')then

              var_2d='slt'
              call read_static_grid(nx_l,ny_l,var_2d,static_grid
     1,istatus)
              if(istatus .ne. 1)then
                 print*,' warning: could not read static-slope-lon'
                 return
              endif
              clow = -1.
              chigh = 1.0
              cint = .025
              c_label = 'latitude component terrain slope  '

              if(cstatic .eq. 'sli')then
                write(6,*)' calling solid fill plot'
                call get_mxmn_2d(nx_l,ny_l,static_grid,rmx2d
     1                          ,rmn2d,imx,jmx,imn,jmn)
                call ccpfil(static_grid,nx_l,ny_l,rmn2d,rmx2d
     1                     ,'linear',n_image,1e0,'hsect',plot_parms
     1                     ,namelist_parms)       
                call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                                ,plot_parms,namelist_parms
     1                                ,i_overlay,'hsect')       
                call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                             ,namelist_parms,plot_parms)
              else
                call plot_cont(static_grid,1e0,
     1               clow,chigh,cint,asc9_tim,
     1               namelist_parms,plot_parms,c_label,
     1               i_overlay,c_display,lat,lon,jdot,
     1               nx_l,ny_l,r_missing_data,laps_cycle_time)
              endif
              i4time_topo = 0

              deallocate (static_grid)

           endif !cstatic

        elseif(c_type .eq. 'cf')then
           call frame
           close(8)

        elseif(c_type .eq. 'q ')then
           goto9000


        endif ! c_field

!       set ifield_found except for intermediate query values
        if(i_plotted_field .eq. 1)then
            write(6,*)' setting ifield_found to 1 in hsect '
            ifield_found = 1
        else
            write(6,*)' not setting ifield_found in hsect'
     1                ,i_plotted_field
        endif

        goto1200

9000    if (.false.) then ! special frame calls
            if(c_display .eq. 'm' .or. c_display .eq. 'p')then
                call frame
            else
                if(c_display .eq. 't')then
                elseif(c_display .eq. 'r')then
                    call frame2(c_display)
                else
                    call frame
                endif
            endif
        endif

 211    format(/' enter yydddhhmmhhmm or hhmm for file: ',a3,$)
 221    format(a13)

        return
        end


        subroutine plot_cont(array,scale,clow,chigh,cint,
     1    asc_tim_9,namelist_parms,plot_parms,c_label,
     1    i_overlay,c_display,lat,lon,jdot,
     1    nx_l,ny_l,r_missing_data,laps_cycle_time)

!       97-aug-14     ken dritz     added nx_l, ny_l as dummy arguments
!       97-aug-14     ken dritz     added r_missing_data, laps_cycle_time
!                                   as dummy arguments
!       97-aug-14     ken dritz     removed include of lapsparms.for

        include 'lapsplot.inc'

        common /mcolor/mini,maxi

        common /plot_field_cmn/ i_plotted_field

        real lat(nx_l,ny_l),lon(nx_l,ny_l)

        character c_label*(*),asc_tim_9*9,c_metacode*2,asc_tim_24*24
        character*1 c_display

        real array(nx_l,ny_l)
        real array_plot(nx_l,ny_l)

!       integer ity,ily,istatus
!       data ity/35/,ily/1010/

        include 'icolors.inc'

        write(6,*)' subroutine plot_cont...',clow,chigh,cint,scale

        i_plotted_field = 1

        y_spacing = 3

        write(6,1505)c_label,scale,asc_tim_9
1505    format(7x,a,4x,'units = ',1pe9.0,6x,a9)

        if(asc_tim_9 .ne. '         ')then
            call i4time_fname_lp(asc_tim_9,i4time_file,istatus)
            call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
!           asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '
        else
            asc_tim_24 = '                        '
        endif

        n_missing = 0
        vmax = -1e30
        vmin = 1e30

        do i = 1,nx_l
        do j = 1,ny_l
            if(array(i,j) .ne. r_missing_data)then
                array_plot(i,j) = array(i,j) / scale 
                vmax = max(vmax,array_plot(i,j))
                vmin = min(vmin,array_plot(i,j))
            else
                array_plot(i,j) = array(i,j) 
                n_missing = n_missing + 1
            endif
        enddo ! i
        enddo ! j

        x_spacing = nx_l / 28

        do j=ny_l,1,-y_spacing
            write(6,500)
     1  (nint(min(max(array_plot(i,j),-99.),999.)),i=1,nx_l,x_spacing)
500         format(1x,42i3)
        enddo ! j

!       set map background stuff
        if(c_display .eq. 'r' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'l' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 't' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'e' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'n' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'o' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'c' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'v')then
            goto990
        elseif(c_display .eq. 'p')then ! generate a map background only
            c_metacode = 'm'
            call lapsplot(array_plot,nx_l,ny_l,clow,chigh,cint
     1                   ,plot_parms,namelist_parms,lat,lon
     1                   ,c_metacode,jdot)
            goto990
        else
            c_metacode = 'c'
        endif

        if(c_metacode .eq. 'c ')then
            i_overlay = i_overlay + 1

        else if(c_metacode .eq. 'm ')then
            write(6,*)' c_metacode,i_overlay = ',c_metacode,i_overlay

            if(c_display .eq. 'r')then
                call lapsplot(array_plot,nx_l,ny_l
     1                       ,clow,chigh,cint,plot_parms
     1                       ,namelist_parms,lat,lon
     1                       ,c_metacode,jdot)
            endif

            c_metacode = 'c '
            i_overlay = 1

!           i4time_plot = i4time_file/laps_cycle_time*laps_cycle_time
!!      1                                            -laps_cycle_time
            call setusv_dum('in',34)

            iflag = 0

            call get_maxstns(maxstns,istatus)
            if (istatus .ne. 1) then
               write (6,*) 'error getting value of maxstns'
               stop
            endif

            write(6,*)' not calling plot_station_locations: 1'
        endif

        write(6,*)' plotting: c_metacode,i_overlay = ',
     1                        c_metacode,i_overlay

        if(clow .eq. 0. .and. chigh .eq. 0. .and. cint .gt. 0.)then
            clow =  (nint(vmin/cint)-1) * cint
            chigh = (nint(vmax/cint)+1) * cint
        endif

        write(6,*)' clow,high,cint ',clow,chigh,cint
        write(6,*)' max/min/n_missing = ',vmax,vmin,n_missing

        call setusv_dum('in',icolors(i_overlay))
        write(6,*)' plotting field, color = '
     1                       ,icolors(i_overlay)

        if(c_metacode .ne. 'n ')then
            if(c_metacode .eq. 'c ')then
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                 call write_label_lplot(nx_l,ny_l,c_label,asc_tim_9
     1                                 ,plot_parms,namelist_parms
     1                                 ,i_overlay,'hsect')      
            endif

            if(c_display .ne. 't')then
                mini = icolors(i_overlay)
                maxi = icolors(i_overlay)
            else
                mini = icolors_tt(i_overlay)
                maxi = icolors_tt(i_overlay)
            endif

            call lapsplot(array_plot,nx_l,ny_l,clow,chigh,cint
     1                   ,plot_parms,namelist_parms
     1                   ,lat,lon,c_metacode,jdot)
        endif

990     return
        end



        subroutine plot_barbs(u,v,lat,lon,topo,size,zoom,
     1  interval,asc_tim_9,namelist_parms,plot_parms,
     1  c_label,
     1  c_field,k_level,i_overlay,c_display,imax,jmax,kmax,max_radars,       
!    1  grid_ra_ref_dum,grid_ra_vel_dum,
     1  nx_l,ny_l,r_missing_data,
     1  laps_cycle_time,jdot)      

        include 'lapsplot.inc'

        common /plot_field_cmn/ i_plotted_field

        character c_label*(*),asc_tim_9*9,c_metacode*2,asc_tim_24*24
        character c_field*2,c_display*1

        real u(nx_l,ny_l)
        real v(nx_l,ny_l)
        real lat(nx_l,ny_l)
        real lon(nx_l,ny_l)
        real topo(nx_l,ny_l)

        real n_plotted(nx_l,ny_l)
        real grid_ra_ref(imax,jmax,kmax)
        real grid_ra_vel(imax,jmax,kmax)

!       integer ity,ily,istatus
!       data ity/35/,ily/1010/

        include 'icolors.inc'

!       i_image: whether this particular plot is an image
        common /image/ n_image, i_image 

        logical l_obs

        call get_grid_spacing_cen(grid_spacing_m,istatus)

        write(6,*)' subroutine plot_barbs: size = ',size

        i_plotted_field = 1

        write(6,*)' i_plotted_field = ',i_plotted_field

        write(6,1505)c_label,asc_tim_9
1505    format(2x,a,2x,a9)

        call i4time_fname_lp(asc_tim_9,i4time_file,istatus)
        call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
!       asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

!       set map background stuff
        if(c_display .eq. 'r' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'l' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 't' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'e' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'n' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'o' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'c' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'v')then
            goto990
        else
            c_metacode = 'c'
        endif

        if(c_metacode .eq. 'c ')then
            i_overlay = i_overlay + 1
        else if(c_metacode .eq. 'm ')then
            write(6,*)' c_metacode,i_overlay = ',c_metacode,i_overlay

            if(c_display .eq. 'r')then
                call lapsplot(array_plot,nx_l,ny_l,clow,chigh,cint
     1                       ,plot_parms,namelist_parms
     1                       ,lat,lon,c_metacode,jdot)
            endif

            c_metacode = 'c '
            i_overlay = 1
c abdel	    
            if (laps_cycle_time.eq.0)then
	      i4time_plot = i4time_file
            else
              i4time_plot = i4time_file/laps_cycle_time*laps_cycle_time
            endif
c abdel	  
            call setusv_dum('in',34)

            iflag = 0

            call get_maxstns(maxstns,istatus)
            if (istatus .ne. 1) then
               write (6,*) 'error getting value of maxstns'
               stop
            endif

            write(6,*)' not calling plot_station_locations: 2'
        endif

        write(6,*)' plotting: c_metacode,i_overlay = ',
     1                          c_metacode,i_overlay

        call setusv_dum('in',icolors(i_overlay))

        if(c_metacode .ne. 'n ')then
            if(c_metacode .eq. 'y ' .or. c_metacode .eq. 'c ')then
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                 call write_label_lplot(nx_l,ny_l,c_label,asc_tim_9
     1                                 ,plot_parms,namelist_parms
     1                                 ,i_overlay,'hsect')      
            endif


            if(c_metacode .eq. 'y ' .or. c_metacode .eq. 'c ')then
                if(c_field .eq. 'ob')then
                    call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltyp
     1e)
                    write(6,2041)
2041                format('         enter radar # (0=all)',40x,'? ',$)       
                    read(5,*)i_radar

                    n_plotted = 0

                    if(i_radar .gt. 0)then         ! single radar plot
                        i_radar_start = i_radar
                        i_radar_end = i_radar
                    else                           ! multi-radar plot
                        i_radar_start = 1
                        i_radar_end = max_radars
                    endif

                    write(6,*)' start/end radars = '
     1                       ,i_radar_start,i_radar_end

!                   call plot_obs(k_level,.true.,asc_tim_9(1:7)//'00'
                    call plot_obs(k_level,.true.,asc_tim_9(1:9)
     1                  ,i_radar_start,i_radar_end
     1                  ,namelist_parms,plot_parms
     1                  ,imax,jmax,kmax,n_plotted
     1                  ,grid_ra_ref,grid_ra_vel,lat,lon,topo
     1                  ,grid_spacing_m,1)
                    return
                endif

                if(n_image .ge. 1 .and. 
     1             plot_parms%icol_barbs .lt. +1)then
                    call setusv_dum('in',22)
                    write(6,*)' plotting quasi black wind barbs'
     1                       ,n_image,plot_parms%icol_barbs
                else
                    call setusv_dum('in',icolors(i_overlay))
                    write(6,*)' plotting wind barbs, color = '
     1                       ,icolors(i_overlay)
                endif

                call get_border(nx_l,ny_l,x_1,x_2,y_1,y_2)

                call set(x_1,x_2,y_1,y_2,1.,float(nx_l),1.,float(ny_l)
     1                                                             ,1)       

                call plot_winds_2d(u,v,interval,size,zoom
     1          ,nx_l,ny_l,lat,lon,r_missing_data,namelist_parms)
!               call frame

            endif

        endif

990     return
        end


        subroutine plot_grid(i_overlay,c_display,lat,lon,
     1                       nx_l,ny_l,laps_cycle_time)

        include 'lapsplot.inc'

!       97-aug-14     ken dritz     added nx_l, ny_l as dummy arguments
!       97-aug-14     ken dritz     added laps_cycle_time as dummy argument
!       97-aug-14     ken dritz     changed laps_domain_file to hardwire
!       97-aug-14     ken dritz     removed include of lapsparms.for

        character c_label*33,asc_tim_9*9,c_metacode*2,asc_tim_24*24
        character c_field*2,c_display*1

        real lat(nx_l,ny_l)
        real lon(nx_l,ny_l)

!       integer ity,ily,istatus
!       data ity/35/,ily/1010/

        include 'icolors.inc'

        logical l_obs

!       set map background stuff
        if(c_display .eq. 'r' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'l' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 't' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'e' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'n' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'o' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'c' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'v')then
            goto990
        else
            c_metacode = 'c'
        endif

        if(c_metacode .eq. 'c ')then
            i_overlay = i_overlay + 1
        else if(c_metacode .eq. 'm ')then
            write(6,*)' c_metacode,i_overlay = ',c_metacode,i_overlay

            if(c_display .eq. 'r')then
                call lapsplot(array_plot,nx_l,ny_l,clow,chigh,cint
     1                       ,plot_parms,namelist_parms
     1                       ,lat,lon,c_metacode,jdot)
            endif

            c_metacode = 'c '
            i_overlay = 1

!           i4time_plot = i4time_file/laps_cycle_time*laps_cycle_time
!       1                                            -laps_cycle_time
            call setusv_dum('in',34)

        endif

        write(6,*)' plotting: c_metacode,i_overlay = ',
     1                          c_metacode,i_overlay

        call setusv_dum('in',icolors(i_overlay))

        if(c_metacode .ne. 'n ')then
            if(c_metacode .eq. 'y '
     1    .or. c_metacode .eq. 'c ')then
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
!                call pwrity(cpux(320),cpux(ity),c_label,33,2,0,0)      
!                call pwrity
!    1                (cpux(800),cpux(ity),asc_tim_24(1:17),17,2,0,0)
                 call write_label_lplot(nx_l,ny_l,c_label,asc_tim_9
     1                                 ,plot_parms,namelist_parms
     1                                 ,i_overlay,'hsect')      
            endif


            if(c_metacode .eq. 'y ' .or. c_metacode .eq. 'c ')then

                call setusv_dum('in',icolors(i_overlay))

!               call get_border(nx_l,ny_l,x_1,x_2,y_1,y_2)

!               call set(x_1,x_2,y_1,y_2,1.,float(nx_l),1.,float(ny_l))

                call plot_grid_2d(nx_l,ny_l,lat,lon)

            endif

        endif

990     return
        end

        subroutine plot_cldpcp_type(icldpcp_type_2d
     1     ,asc_tim_9,namelist_parms,plot_parms
     1     ,c_label,c_field,k_level,i_overlay,c_display
     1     ,lat,lon,ifield_2d
     1     ,nx_l,ny_l,laps_cycle_time,jdot)

!       97-aug-14     ken dritz     added nx_l, ny_l, laps_cycle_time as
!                                   dummy arguments
!       97-aug-14     ken dritz     removed include of lapsparms.for

        include 'lapsplot.inc'

        character c_label*33,asc_tim_9*9,c_metacode*2,asc_tim_24*24
        character c_field*2,c_display*1

        integer icldpcp_type_2d(nx_l,ny_l)
        real lat(nx_l,ny_l)
        real lon(nx_l,ny_l)
        integer ifield_2d(nx_l,ny_l)

        integer iarg

!       integer ity,ily,istatus
!       data ity/35/,ily/1010/

        include 'icolors.inc'

        logical l_obs

        write(6,1505)c_label,asc_tim_9
1505    format(2x,a33,2x,a9)

        call i4time_fname_lp(asc_tim_9,i4time_file,istatus)
        call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
!       asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

!       set map background stuff
        if(c_display .eq. 'r' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'l' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 't' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'e' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'n' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'o' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'c' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'v')then

            interval = 2
            size = 8.0

            call plot_types_2d(icldpcp_type_2d,interval,size,c_field
     1                        ,.false.,plot_parms
     1                        ,nx_l,ny_l,lat,lon,ifield_2d)
            goto990
        else
            c_metacode = 'c'
        endif

        if(c_metacode .eq. 'c ')then
            i_overlay = i_overlay + 1
        else if(c_metacode .eq. 'm ')then
            write(6,*)' c_metacode,i_overlay = ',c_metacode,i_overlay

            call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                         ,namelist_parms,plot_parms)       

            c_metacode = 'c '
            i_overlay = 1

c abdel       
            if (laps_cycle_time.eq.0)then
               i4time_plot = i4time_file
	    else
               i4time_plot = i4time_file/laps_cycle_time*laps_cycle_time
!       1                                               -laps_cycle_time
            endif
            call setusv_dum('in',34)

            iflag = 0

            call get_maxstns(maxstns,istatus)
            if (istatus .ne. 1) then
               write (6,*) 'error getting value of maxstns'
               stop
            endif

            write(6,*)' not calling plot_station_locations: 3'
        endif

        write(6,*)' plotting: c_metacode,i_overlay = ',
     1                          c_metacode,i_overlay

        call i4time_fname_lp(asc_tim_9,i4time_file,istatus)
        call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
!       asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

        call setusv_dum('in',icolors(i_overlay))

        if(c_metacode .ne. 'n ')then
            if(c_metacode .eq. 'y '
     1    .or. c_metacode .eq. 'c ')then
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                 call write_label_lplot(nx_l,ny_l,c_label,asc_tim_9
     1                                 ,plot_parms,namelist_parms
     1                                 ,i_overlay,'hsect')      
            endif


            if(c_metacode .eq. 'y ' .or. c_metacode .eq. 'c ')then
                call setusv_dum('in',icolors(i_overlay))

                if(max(nx_l,ny_l) .gt. 61)then
                    interval = 4
                else
                    interval = 2
                endif

                size = 8.0

                call plot_types_2d(icldpcp_type_2d,interval,size,c_field      
     1                            ,.true.,plot_parms,nx_l,ny_l,lat,lon
     1                            ,ifield_2d)
!               call frame

            endif

        endif

990     return
        end

        subroutine plot_stations(asc_tim_9,c_label,c_field,i_overlay
     1   ,namelist_parms,plot_parms,max_snd_grid,max_snd_levels
     1   ,c_display,lat,lon,topo,c_file,iflag
     1   ,nx_l,ny_l,nz_l,laps_cycle_time,zoom)

!       97-aug-14     ken dritz     added nx_l, ny_l, laps_cycle_time as
!                                   dummy arguments
!       97-aug-14     ken dritz     removed include of lapsparms.for

        character c_label*(*),asc_tim_9*9,c_metacode*2,asc_tim_24*24
        character c_field*4,c_display*1
        character*(*) c_file

        real lat(nx_l,ny_l)
        real lon(nx_l,ny_l)
        real topo(nx_l,ny_l)

        integer iarg

!       integer ity,ily,istatus
!       data ity/35/,ily/1010/

        include 'icolors.inc'
        include 'lapsplot.inc'

        logical l_obs

        write(6,1505)c_label(1:33),asc_tim_9
1505    format(2x,a33,2x,a9)

        call i4time_fname_lp(asc_tim_9,i4time_file,istatus)
        call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
!       asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

!       set map background stuff
        if(c_display .eq. 'r' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'v')then
            goto990
        else
            c_metacode = 'c'
        endif

        if(c_metacode .eq. 'c ')then
            i_overlay = i_overlay + 1
        else if(c_metacode .eq. 'm ')then

            if(iflag .eq. 1)then
                jdot = 0 ! solid boundaries
            else
                jdot = 1 ! dotted boundaries
            endif

            write(6,*)' c_metacode,i_overlay = ',c_metacode,i_overlay
            write(6,*)' iflag,jdot = ',iflag,jdot

            call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                         ,namelist_parms,plot_parms)

            c_metacode = 'c '
            i_overlay = 1
c abdel       
            if (laps_cycle_time.eq.0)then
	        i4time_plot = i4time_file
	    else
                i4time_plot = i4time_file ! /laps_cycle_time*laps_cycle_time
!       1                                            -laps_cycle_time
	    endif

            call setusv_dum('in',34) ! grey
!           call setusv_dum(2hin,11)

            call get_maxstns(maxstns,istatus)
            if (istatus .ne. 1) then
               write (6,*) 'error getting value of maxstns'
               stop
            endif

            write(6,*)' not calling plot_station_locations: 4 '
     1               ,iflag,c_metacode

        endif

        write(6,*)' plotting: c_metacode,i_overlay = ',
     1                          c_metacode,i_overlay

        call i4time_fname_lp(asc_tim_9,i4time_file,istatus)
        call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
!       asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

        call setusv_dum('in',icolors(i_overlay))

        if(c_metacode .ne. 'n ')then
            if(c_metacode .eq. 'y '
     1  .or. c_metacode .eq. 'c ')then
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
!                call pwrity
!    1                (cpux(320),cpux(ity),c_label,33,2,0,0)
!                call pwrity
!    1                (cpux(800),cpux(ity),asc_tim_24(1:17),17,2,0,0)

!                if(iflag .ge. 1)then
!                    c_label = 'sfc obs'
!                endif

!                call write_label_lplot(nx_l,ny_l,c_label,asc_tim_9
!    1                                 ,plot_parms,namelist_parms
!    1                                 ,i_overlay,'hsect')      
            endif

            if(c_metacode .eq. 'y ' .or. c_metacode .eq. 'c ')then

                if(namelist_parms%c_ob_color .eq. 'white')then
                    call setusv_dum('in',1)                     ! white
                elseif(iflag .eq. 2)then ! obs with station locations?
                    call setusv_dum('in',icolors(i_overlay))
                else                 ! station locations by themselves
                    call setusv_dum('in',34)                    ! grey
                endif

                if(max(nx_l,ny_l) .gt. 61)then
                    interval = 4
                else
                    interval = 2
                endif

                size = 1.0

                call get_maxstns(maxstns,istatus)
                if(istatus .ne. 1) then
                   write (6,*) 'error getting value of maxstns'
                   stop
                endif

                write(6,*)' calling plot_station_locations: 5 '
     1                   ,iflag,c_metacode       
                call plot_station_locations(i4time_file,lat,lon,topo
     1                      ,nx_l,ny_l,nz_l,iflag,maxstns,c_field,zoom
     1                      ,namelist_parms,plot_parms
     1                      ,max_snd_grid,max_snd_levels                    ! i
     1                      ,asc_tim_24,c_label,i_overlay)       
            endif

        endif

990     return
        end


        subroutine mklabel(k_level,c19_label,c_label)

!       97-aug-17     ken dritz     removed include of lapsparms.for

        character c19_label*(*),c_label*(*)
        character*40 vert_grid

        call get_vertical_grid(vert_grid,istatus)

!       initialize with blanks
        len_label = len(c_label)
        do i = 1,len_label
            c_label(i:i) = ' '
        enddo ! i

        if(k_level .gt. 0)then
             if(vert_grid .eq. 'sigma_ht')then
                 write(c_label,101)k_level,c19_label
101              format(i5,' m ',a,5x)

             elseif(vert_grid .eq. 'pressure')then
                if(k_level .gt. 50)then ! k_level is given in pressure (mb)
                    ipres = k_level

                else                    ! k_level is level number
                    ipres = k_level
!                   ipres = nint(zcoord_of_level(k_level)/100.)

                endif

                write(c_label,102)ipres,c19_label
102             format(i4,' hpa',a,5x)

             endif
        else if(k_level .eq. 0)then
            write(c_label,103)c19_label
103         format('surface',a,6x)

        else if(k_level .eq. -1)then
            write(c_label,104)
104         format('steering winds                   ')

        endif

        return
        end

        subroutine get_border(ni,nj,x_1,x_2,y_1,y_2)

        if(ni .eq. nj)then
            x_1 = .05
            x_2 = .95
            y_1 = .05
            y_2 = .95
        elseif(ni .lt. nj)then
            ratio = float(ni-1) / float(nj-1)
            x_1 = .50 - .45 * ratio
            x_2 = .50 + .45 * ratio
            y_1 = .05
            y_2 = .95
        elseif(ni .gt. nj)then
            ratio = float(nj-1) / float(ni-1)
            x_1 = .05
            x_2 = .95
            y_1 = .50 - .45 * ratio
            y_2 = .50 + .45 * ratio
        endif

        return
        end



        subroutine setusv_dum(c2_dum,icol_in)

        character*2 c2_dum

        common /icol_index/ icol

!       icol = min(icol_in,35)
        icol = icol_in

!       write(6,*)' color # ',icol,icol_in

        call gstxci(icol)            
        call gsplci(icol)          
        call gspmci(icol)           
        call gsfaci(icol)                 

        call pcseti('cc',icol)
        call pcseti('oc',icol)
        call pcseti('sc',icol)

        return
        end


        subroutine write_label_lplot(ni,nj,c_label,a9time
     1                              ,plot_parms,namelist_parms
     1                              ,i_overlay,c5_sect)       

        include 'lapsplot.inc'

        character*(*) c_label
        character*24 asc_tim_24,asc_tim_24_in
        character*4 c4_grid
        character*5 c5_sect,c5_grid
        character*9 a9time
        character*20 c_ul 

        common /image/ n_image
        common /icol_index/ icol_common

!       namelist_parms%time_zone = 0.
!       namelist_parms%c3_time_zone = 'utc'

!       call upcase(c_label,c_label)

        write(6,*)' write_label_lplot:',c_label

        call s_len2(c_label,len_label)
        do i = 1,len_label
            if(c_label(i:i) .eq. ':')then
                c_label(i:i) = ' '
            endif
        enddo ! i

        if(a9time .ne. '         ')then
            call i4time_fname_lp(a9time,i4time_lbl,istatus)       
            i4time_lbl = i4time_lbl+nint(namelist_parms%time_zone*3600)      
            call cv_i4tim_asc_lp(i4time_lbl,asc_tim_24_in,istatus)      
        else
            asc_tim_24_in = '                        '
        endif

        asc_tim_24 = asc_tim_24_in(1:14)//asc_tim_24_in(16:17)//' '      

        i_label = i_overlay + n_image

        call get_border(ni,nj,x_1,x_2,y_1,y_2)

        jsize_t = 2 

!       top label
        y_2 = y_2 + .0225 ! .025

        icol_save = icol_common

        call setusv_dum('  ',7)

!       add grid resolution for top label (short format)
        call get_grid_spacing_cen(grid_spacing_m,istatus)

        if(grid_spacing_m .ge. 999.5)then ! 1-km or greater
            igrid_spacing = nint(grid_spacing_m/1000.)
            if(igrid_spacing .le. 99 .and. igrid_spacing .ge. 1
     1                               .and. c5_sect .ne. 'xsectx'  )then       
                write(c4_grid,1)igrid_spacing
 1              format(i2,'km')
            else
                c4_grid = '    '
            endif
        else                              ! < 1-km
            igrid_spacing = nint(grid_spacing_m)
            if(igrid_spacing .le. 999 .and. igrid_spacing .ge. 1
     1                                .and. c5_sect .ne. 'xsectx' )then       
                write(c4_grid,2)igrid_spacing
 2              format(i3,'m')
            else
                c4_grid = '    '
            endif
        endif

!       grid resolution (long format)
        if(grid_spacing_m .ge. 999.5 .and. 
     1     grid_spacing_m .le. 9950.       )then             
            ihundreds = nint(grid_spacing_m / 100.)
            if(((ihundreds/10) * 10) .ne. ihundreds)then
                write(c5_grid,3)nint(grid_spacing_m / 100.) / 10.0
 3              format(f3.1,'km')
            else
                c5_grid = c4_grid
            endif
        else
            c5_grid = c4_grid
        endif

!       set for zoom
        zfrac = 1.0 / plot_parms%zoom_wdw

!       note that "square" case works for aspect ratio up to 1.192
        frame_factx = 1.0  ! / 0.75
        frame_facty = 1.0  ! / 0.8

        zxcen = (0.5 + ((plot_parms%xcen - 0.5) * frame_factx)) * 1023.
        zycen = (0.5 + ((plot_parms%ycen - 0.5) * frame_facty)) * 1023.

!       top left label               

        iy = y_2 * 1024
!       ix = 170 

        if(c5_sect .eq. 'hsect')then
            ix = 70
            rsize = .013
        else
            ix = 70
            rsize = .011
        endif

        ix = zxcen + (float(ix-512) / plot_parms%zoom_wdw)
        iy = zycen + (float(iy-512) / plot_parms%zoom_wdw)

        call s_len2(namelist_parms%c_institution,len_inst)

        if(c5_sect .eq. 'xsect')then
            len_inst = min(len_inst,14)

            if(len_inst .le. 9)then
                c_ul = namelist_parms%c_institution(1:len_inst)
     1                 //' laps '//c4_grid
            else
                c_ul = namelist_parms%c_institution(1:len_inst)//' '
     1                 //c4_grid
            endif

        else ! hsect

            if(len_inst .le. 9 .and. .false.)then
                c_ul = namelist_parms%c_institution(1:len_inst)
     1                 //' laps '//c5_grid
            else
                c_ul = namelist_parms%c_institution(1:len_inst)//' '
     1                 //c5_grid
            endif

        endif

        rsize_zoom = rsize/plot_parms%zoom_wdw

        call pchiqu (cpux(ix),cpux(iy),c_ul,rsize_zoom,0,-1.0)   

!       if(c5_sect .eq. 'sound')then
!           call pwrity(cpux(ix),cpux(iy),'noaa/fsl',8,jsize_t,0,0)
!       endif

        call setusv_dum('  ',icol_save)

!       bottom label
        jsize_b = 1 ! [valid range is 0-2]
        rsize_b = jsize_b + 2.

        if(c5_sect .eq. 'hsect')then
            vspace=.0065
            v1=.0040
            if(i_label .ge. 3)i_label = 3 - i_label ! add labels on top
        else
            vspace=.0075
            v1=.0045
        endif

        if(jsize_b .eq. 2)then
            y_1 = y_1 - .025 - .035 * float(i_label-1)
        else
            y_1 = y_1 - rsize_b*v1 - rsize_b * vspace
     1                              * float(i_label-1)
        endif

        rsize = .010 ! .011

!       field on bottom left
        if(c5_sect .eq. 'xsect')then
            ix = 130
        else
            ix = 100 ! 130
        endif

        iy = y_1 * 1024

!       if zoom is one, coordinates are original ix,iy
!       as zoom goes towards infinity coordinates go zxcen,zycen
!       ix = nint( (float(ix) * zfrac) + zxcen * (1.0 - zfrac))
        ix = zxcen + (float(ix-512) / plot_parms%zoom_wdw)
!       iy = nint( (float(iy) * zfrac) + zycen * (1.0 - zfrac))
!       iy = nint(                       zycen) ! for testing
        iy = zycen + (float(iy-512) / plot_parms%zoom_wdw)

        if(plot_parms%zoom_wdw .gt. 0.0)then
            write(6,*)'frame_factxy,zxcen,zycen,ixy'
     1               ,frame_factx,frame_facty,zxcen,zycen,ix,iy
        endif

        rsize_zoom = rsize/plot_parms%zoom_wdw

        if(len_label .gt. 0)then
            call pchiqu (cpux(ix),cpux(iy),c_label(1:len_label)
     1                  ,rsize_zoom,0,-1.0)      
        else
            write(6,*)' note that label has zero length...'
        endif

        if(.false.)then ! another test
            do ry = 0.0,1.0,.10
                iy = nint(ry * 1023.)
                write(c_label,11)iy
11              format(i4)
                len_label = 4
                write(6,*)' test label ',c_label(1:4),ix,iy
                call pchiqu (cpux(ix),cpux(iy),c_label(1:len_label)
     1                  ,rsize_zoom,0,-1.0)      
            enddo
        endif

!       time on bottom right
        if(c5_sect .eq. 'hsect')then
            ix = 672
        else
            ix = 622
        endif
        iy = y_1 * 1024

!       set for zoom
!       ix = nint( (float(ix) * zfrac) + zxcen * (1.0 - zfrac))
        ix = zxcen + (float(ix-512) / plot_parms%zoom_wdw)
!       iy = nint( (float(iy) * zfrac) + zycen * (1.0 - zfrac))
        iy = zycen + (float(iy-512) / plot_parms%zoom_wdw)

        call downcase(asc_tim_24(5:10),asc_tim_24(5:10))
        call pchiqu (cpux(ix),cpux(iy),'vt '//asc_tim_24(1:17)
     1                                //namelist_parms%c3_time_zone
     1                                ,rsize_zoom,0,-1.0)

!       resolution on bottom center
!       ix = 520
!       iy = y_1 * 1024
!       call pchiqu (cpux(ix),cpux(iy),c4_grid,rsize,0,-1.0)

        return
        end

        subroutine plot_field_2d(i4time,c_type_in,field_2d,scale
     1                        ,namelist_parms,plot_parms
     1                        ,clow_in,chigh_in,cint_in
     1                        ,c_label,i_overlay,c_display,lat,lon
     1                        ,jdot,nx_l,ny_l,r_missing_data,colortable)

        include 'lapsplot.inc'

        character*(*) c_type_in, c_label, c_display, colortable
        character*10 c_type
        character*9 asc9_tim

        real field_2d(nx_l,ny_l)
        real field_2d_plot(nx_l,ny_l)

        real lat(nx_l,ny_l)
        real lon(nx_l,ny_l)

!       i_image: whether this particular plot is an image
        common /image/ n_image, i_image 
        common /zoom/  zoom, density

        write(6,*)' subroutine plot_field_2d: i4time = ',i4time

        write(6,*)' inputs clow/chigh/cint = ',clow_in,chigh_in,cint_in

        c_type = c_type_in

        call make_fnam_lp(i4time,asc9_tim,istatus)

        if(clow_in .gt. chigh_in)then
            chigh = clow_in
            clow = chigh_in
        else
            chigh = chigh_in
            clow = clow_in
        endif

        if(cint_in .lt. 0.)then
!           clow_img = abs(cint_in)
            clow_img = clow_in
            chigh_img = chigh_in
        else
            clow_img = clow_in
            chigh_img = chigh_in
        endif

        call s_len(c_type,len_type)

        call downcase(c_type,c_type)

        if(plot_parms%l_hinterp_zoom)then
            field_2d_plot = field_2d ! interpolate based on zoom
        else
            field_2d_plot = field_2d
        endif

!       if( (c_type(len_type:len_type)  .ne. 'i' .or. c_type  .eq. 'hi')       
!    1                          .and. 
        if(                   i_image .eq. 0                       )then       
            write(6,*)' plot_field_2d - contour plot ',c_type
            if(cint_in .eq. 0.)then
                call contour_settings(field_2d_plot,nx_l,ny_l
     1                               ,clow,chigh,cint       
     1                               ,zoom,density,scale)

            elseif(clow_in .eq. 0. .and. chigh_in .eq. 0.)then
                call contour_settings(field_2d_plot,nx_l,ny_l
     1                               ,clow,chigh,cint       
     1                               ,zoom,density,scale)

                cint = cint_in

            else
                cint = cint_in
                chigh = chigh + 5. * cint ! add 5 extra contours

            endif

            call plot_cont(field_2d_plot,scale,clow,chigh,cint
     1                        ,asc9_tim,namelist_parms,plot_parms
     1                        ,c_label,i_overlay
     1                        ,c_display,lat,lon,jdot,nx_l,ny_l
     1                        ,r_missing_data,laps_cycle_time)

        else ! image plot
            write(6,*)' plot_field_2d - image plot ',c_type
     1               ,clow_img,chigh_img,plot_parms%color_power

            if(clow_img .eq. 0. .and. chigh_img .eq. 0.)then
                call contour_settings(field_2d_plot,nx_l,ny_l
     1                               ,clow_img,chigh_img,cint       
     1                               ,zoom,density,scale)
            endif

            if(colortable .eq. 'linear')then
                    plot_parms%icol_barbs = +1 ! keep future barbs plots bright
            endif

            call ccpfil(field_2d_plot,nx_l,ny_l
     1                 ,clow_img,chigh_img ! *scale      
     1                 ,colortable,n_image,scale,'hsect',plot_parms
     1                 ,namelist_parms)       
            call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
            call setusv_dum('in',7)
            call write_label_lplot(nx_l,ny_l,c_label,asc9_tim
     1                            ,plot_parms,namelist_parms
     1                            ,i_overlay,'hsect')       
            call lapsplot_setup(nx_l,ny_l,lat,lon,jdot
     1                         ,namelist_parms,plot_parms)

        endif ! image plot

        return
        end

        subroutine mk_fcst_hlabel(k_mb,comment_2d,fcst_hhmm_in
     1                                ,ext
     1                                ,units_2d
     1                                ,c_model
     1                                ,c_label)

        use mem_namelist, only: model_fcst_intvl

        character*(*) comment_2d,ext,units_2d,c_model,c_label

        character*5 fcst_hhmm_in,fcst_hhmm

        character*40 vert_grid

        ! yuanfu added this declaration of logical function l_parse:
        logical l_parse

        call get_vertical_grid(vert_grid,istatus)

        len_model_max = 7

        c_label = ' '

        call s_len2(comment_2d,len_fcst)
        write(6,*)'comment_2d = ',comment_2d(1:len_fcst)

        call s_len2(units_2d,len_units)
        write(6,*)'units_2d = ',units_2d(1:len_units)

        ! yuanfu removed ".eqv. .true." as it is unnecessary and problem for gfortran:
        if(l_parse(units_2d,'percent'))then
            units_2d = '%'
            len_units = 1
        endif

        call s_len2(c_model,len_model)

        call s_len(fcst_hhmm_in,length_fcst_in)

        write(6,*)' mk_fcst_hlabel: len_model = ',len_model

        if(ext .eq. 'lga' .and. len_model .eq. 0)then
            write(6,*)' c_model has zero length, using lga in label'
            c_model = 'lga'
        endif

        call s_len2(c_model,len_model)
        call upcase(c_model,c_model)

        if(k_mb .gt. 0)then      ! 3d field
            if(vert_grid .eq. 'pressure')then
                write(c_label,102)k_mb
102             format(i4,' hpa ')
            else
                write(c_label,103)k_mb
103             format(i5,' m  ')
            endif
            ic = 10  ! position where comment info should begin
        elseif(k_mb .eq. -1)then ! column max                
            c_label(1:8) = 'col max '
            ic = 9
        else                     ! sfc field
            ic = 1   ! position where comment info should begin
        endif

        if(fcst_hhmm_in(length_fcst_in-1:length_fcst_in) .eq. '00'
     1                 .and.            model_fcst_intvl .ge. 3600 )then       
            fcst_hhmm = fcst_hhmm_in(1:length_fcst_in-2)//'hr '
        else
            fcst_hhmm = fcst_hhmm_in
        endif

        ist = 36 ! position where forecast time should begin

        if(len_units .gt. 0)then
            c_label(ic:ic+len_fcst+len_units+2) = comment_2d(1:len_fcst)
     1                         //' ('//units_2d(1:len_units)//')'
        else
            c_label(ic:ic+len_fcst+len_units) = comment_2d(1:len_fcst)
        endif

!       fcst time info
        c_label(ist:ist+length_fcst_in+1) = 
     1                fcst_hhmm(1:length_fcst_in)//' '

        if(c_model(4:4) .eq. '-' .and. 
     1     len_model .gt. len_model_max)then 
            len_model_label = len_model - 4
            len_model_label = min(len_model_label,len_model_max)
            ims = ist + length_fcst_in + 1
            ime = ims + len_model_label + 5 
            c_label(ims:ime) = c_model(5:len_model)//' fcst'       
        elseif(len_model .gt. 0)then
            len_model = min(len_model,len_model_max)
            ims = ist + length_fcst_in + 1
            ime = ims + len_model + 4 
            c_label(ims:ime) = c_model(1:len_model)//' fcst'       
        else
            c_label(ist+length_fcst_in+1:ist+length_fcst_in+4) = 'fcst'       
        endif

        write(6,*)'c_label = ',c_label

        return
        end


        subroutine directory_to_cmodel(directory,c_model)

        character*(*) directory,c_model
        character*255 dir_local

        call s_len(directory,len_full)
        if(directory(len_full:len_full) .eq. '/')then
            dir_local = directory(1:len_full-1) ! trim slash at the end
        else
            dir_local = directory(1:len_full) 
        endif        

        call get_directory_length(dir_local,len_path)

        call s_len(dir_local,len_path_and_model)

        write(6,*)' directory = ',directory
        write(6,*)' dir_local = ',dir_local

        istart = len_path+1
        if(istart .le. len_path_and_model)then        
            c_model = dir_local(istart:len_path_and_model)
        else
            write(6,*)' could not decode model name '
        endif

        write(6,*)' c_model = ',c_model

        return
        end
