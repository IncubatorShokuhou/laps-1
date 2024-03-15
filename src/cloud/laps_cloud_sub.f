
        subroutine laps_cloud_sub(i4time,
     1                        nx_l,ny_l,
     1                        nz_l,
     1                        n_pirep,
     1                        maxstns,
     1                        max_cld_snd,
     1                        i_diag,
     1                        n_prods,
     1                        iprod_number,
     1                        isplit,
     1                        j_status)

!       obtain cloud parameters
        use mem_namelist, only: l_use_vis,l_use_vis_add                
     1                         ,l_use_vis_partial                      
     1                         ,l_use_39,latency_co2                   
     1                         ,pct_req_lvd_s8a                        
     1                         ,cld_weight_modelfg
     1                         ,i4_sat_window,i4_sat_window_offset     
     1                         ,l_use_metars,l_use_radar,iwrite_output
     1                         ,l_use_pireps,i_varadj,l_corr_parallax

        include 'trigd.inc'
        include 'cloud.inc'

        integer       ss_normal,sys_bad_prod,sys_no_data,
     1                sys_abort_prod

        parameter (ss_normal      =1, ! success
     1             sys_bad_prod   =2, ! inappropriate data, insufficient data
     1             sys_no_data    =3, ! no data
     1             sys_abort_prod =4) ! failed to make a prod

!       prevents clearing out using satellite (hence letting saos dominate)
!       below this altitude (m agl)
        real surface_sao_buffer
        parameter (surface_sao_buffer = 800.)

        real thresh_cvr_ceiling,default_top,default_base
     1                           ,default_clear_cover,default_ceiling       

        parameter       (thresh_cvr_ceiling = 0.65) ! used to "binaryize" cloud cover

        parameter       (default_clear_cover = .001) 

        real thresh_cvr_base,thresh_cvr_top
        parameter (thresh_cvr_base = 0.1)
        parameter (thresh_cvr_top  = 0.1)

        real thresh_thin_lwc_ice     ! threshold cover for thin cloud lwc/ice
        parameter (thresh_thin_lwc_ice = 0.1)

        real vis_radar_thresh_cvr,vis_radar_thresh_dbz
        parameter (vis_radar_thresh_cvr = 0.2)  ! 0.2, 0.0
        parameter (vis_radar_thresh_dbz = 30.)  ! 10., 5. , -99.

        real lat(nx_l,ny_l),lon(nx_l,ny_l)
        real topo(nx_l,ny_l)
        real rlaps_land_frac(nx_l,ny_l)
        real solar_alt(nx_l,ny_l)
        real solar_az(nx_l,ny_l)
        real solar_ha(nx_l,ny_l)

        logical l_packed_output 
        logical l_use_co2_mode1, l_use_co2_mode2
        logical l_evap_radar, l_get_cloudtype
        logical l_trust_narrowband ! should narrowband be trusted (i.e. not 
                                   ! qc'd out) in the absence of other data?

        logical lstat_co2_a(nx_l,ny_l)

        data l_packed_output /.false./
        data l_evap_radar /.false./
        data l_trust_narrowband /.false./ ! overriden to .true. for nowrad

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

        equivalence (cld_hts,cld_hts_new)

        real, allocatable, dimension(:,:,:) :: cldcv1
        real, allocatable, dimension(:,:,:) :: cf_modelfg
        real, allocatable, dimension(:,:,:) :: t_modelfg
        real, allocatable, dimension(:,:,:) :: sh_modelfg
        real, allocatable, dimension(:,:,:) :: ref_modelfg
!       real cldcv1(nx_l,ny_l,kcloud)
!       real cf_modelfg(nx_l,ny_l,kcloud)
!       real t_modelfg(nx_l,ny_l,kcloud)

        real clouds_3d(nx_l,ny_l,kcloud)

        integer ista_snd(max_cld_snd)

        real, allocatable, dimension(:,:) :: cld_snd
        real, allocatable, dimension(:,:) :: wt_snd
!       real cld_snd(max_cld_snd,kcloud)
!       real wt_snd(max_cld_snd,kcloud)

        real cvr_snd(max_cld_snd)
        integer i_snd(max_cld_snd)
        integer j_snd(max_cld_snd)

        integer ihist_alb(-10:20)

        real cloud_top(nx_l,ny_l)
        real cloud_base(nx_l,ny_l)
        real cloud_ceiling(nx_l,ny_l)

        real cldtop_m(nx_l,ny_l)
        real cldtop_co2_m(nx_l,ny_l)
        real cldtop_tb8_m(nx_l,ny_l)
        real cldtop_co2_pa_a(nx_l,ny_l)
        real ht_sao_top(nx_l,ny_l)

        real cldcv_sao(nx_l,ny_l,kcloud)
        real cld_pres_1d(kcloud)
        real pressures_pa(nz_l)
        real pres_3d(nx_l,ny_l,nz_l)
        real wtcldcv(nx_l,ny_l,kcloud)

        real cvhz(nx_l,ny_l)
        real cvhz1(nx_l,ny_l),cvew1(nx_l,kcloud)
        real cvr_max(nx_l,ny_l),cvew2(nx_l,kcloud)
        real cvr_sao_max(nx_l,ny_l)
        real cvr_snow_cycle(nx_l,ny_l)
        real cvr_water_temp(nx_l,ny_l)
        real cvr_snow(nx_l,ny_l)
        real vis_snow_max(nx_l,ny_l)
        real plot_mask(nx_l,ny_l)
        real plot_maskr(nx_l,ny_l)

        character*4 radar_name
        character*31 radarext_3d_cloud
        real radar_ref_3d(nx_l,ny_l,nz_l)
        real closest_radar(nx_l,ny_l)
        integer istat_radar_2dref_a(nx_l,ny_l)
        integer istat_radar_3dref_a(nx_l,ny_l)
        logical lstat_radar_3dref_orig_a(nx_l,ny_l)
       
        real heights_3d(nx_l,ny_l,nz_l)

!       output array declarations
        real out_array_3d(nx_l,ny_l,9)

!       real snow_2d(nx_l,ny_l)

        character*2 c2_precip_types(0:10)

        data c2_precip_types
     1  /'  ','rn','sn','zr','sl','ha','l ','zl','  ','  ','  '/

        integer i2_pcp_type_2d(nx_l,ny_l)
        real r_pcp_type_2d(nx_l,ny_l)

        real dum1_array(nx_l,ny_l)
        real dum2_array(nx_l,ny_l)
        real dum3_array(nx_l,ny_l)
        real dum4_array(nx_l,ny_l)

!       arguments for calling get_cloud_deriv
        logical l_mask_pcptype(nx_l,ny_l)
        integer ibase_array(nx_l,ny_l)
        integer itop_array(nx_l,ny_l)
        logical   l_flag_cloud_type,l_flag_mvd,l_flag_icing_index
     1           ,l_flag_bogus_w,l_bogus_radar_w,l_deep_vv
     1           ,l_flag_pcp_type

        logical l_unresolved(nx_l,ny_l)

        character*1 c1_name_array(nx_l,ny_l,kcloud)

        integer max_fields
        parameter (max_fields = 10)

        character*255 c_filespec
        character var*3,comment*125,directory*150,ext*31,units*10
        character*125 comment_tb8,comment_t39,comment_sst,comment_alb       
        character*3 exts(20)
        character*3 var_a(max_fields)
        character*125 comment_a(max_fields)
        character*10  units_a(max_fields)

!       arrays used to read in satellite data
        real tb8_k(nx_l,ny_l)
        real tb8_k_offset(nx_l,ny_l)
        real t39_k(nx_l,ny_l)
        real tb8_cold_k(nx_l,ny_l)
        real sat_albedo(nx_l,ny_l) ! cloud albedo from reflectance/sfc alb
        real sat_refl(nx_l,ny_l)   ! satellite reflectance
        real static_albedo(nx_l,ny_l)              ! static albedo database
        real sfc_albedo(nx_l,ny_l)           
        real cloud_frac_vis_a(nx_l,ny_l) ! cloud albedo with clear alb subtracted
        real cloud_frac_co2_a(nx_l,ny_l)
        real subpoint_lat_clo_vis(nx_l,ny_l)
        real subpoint_lon_clo_vis(nx_l,ny_l)
        real subpoint_lat_clo_s8a(nx_l,ny_l)
        real subpoint_lon_clo_s8a(nx_l,ny_l)
        real di_dh_ir(nx_l,ny_l)                      
        real dj_dh_ir(nx_l,ny_l)                      
        real di_dh_vis(nx_l,ny_l)                      
        real dj_dh_vis(nx_l,ny_l)                      
        real offset_ir_i(nx_l,ny_l)     ! sat i minus actual i
        real offset_ir_j(nx_l,ny_l)     ! sat j minus actual j
        real offset_vis_i(nx_l,ny_l)    ! sat i minus actual i
        real offset_vis_j(nx_l,ny_l)    ! sat j minus actual j
        real buff(nx_l,ny_l)
        real cldht_prlx_top(nx_l,ny_l)
        real cldht_prlx_unsm(nx_l,ny_l)
        integer i_fill_seams(nx_l,ny_l) ! differentiate ir/vis?

        integer istat_39_a(nx_l,ny_l)
        integer istat_39_add_a(nx_l,ny_l)
        integer istat_vis_potl_a(nx_l,ny_l)  ! image space
        integer istat_vis_added_a(nx_l,ny_l) ! gridpoint space?

        real cloud_albedo(nx_l,ny_l)    ! cloud albedo (corrected parallax)
        real cloud_od(nx_l,ny_l)        ! cloud optical depth
        real cloud_op(nx_l,ny_l)        ! cloud opacity

        real temp_3d(nx_l,ny_l,nz_l)

        real t_sfc_k(nx_l,ny_l)
        real t_gnd_k(nx_l,ny_l)
        real sst_k(nx_l,ny_l)
        real td_sfc_k(nx_l,ny_l)
        real tgd_sfc_k(nx_l,ny_l)
        real pres_sfc_pa(nx_l,ny_l)

!       declarations for lso file stuff
        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real cover_s(maxstns)
        real t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns)
        real ffg_s(maxstns), vis_s(maxstns)
        real pstn_s(maxstns),pmsl_s(maxstns),alt_s(maxstns)
        real store_hgt(maxstns,5),ceil(maxstns),lowcld(maxstns)
        real cover_a(maxstns),rad_s(maxstns),solar_ea(maxstns)
        integer obstime(maxstns),kloud(maxstns),idp3(maxstns)
        character store_emv(maxstns,5)*1,store_amt(maxstns,5)*4
        character wx_s(maxstns)*8
        character atype(maxstns)*6
        character reptype(maxstns)*6

        integer station_name_len
        parameter (station_name_len = 3)                   
        character c_stations(maxstns)*(station_name_len)    

        real ri_s(maxstns), rj_s(maxstns)

!       product # notification declarations
        integer j_status(20),iprod_number(20)

!       stuff for 2d fields
        real ref_mt_modelfg(nx_l,ny_l)
        real dbz_low_2d(nx_l,ny_l)
        real dbz_max_2d(nx_l,ny_l)
        real rqc_2d(nx_l,ny_l)
        real swi_2d(nx_l,ny_l)

!       sfc precip and cloud type (lct file)
        real r_pcp_type_thresh_2d(nx_l,ny_l)
        real r_cld_type_2d(nx_l,ny_l)

        character*40 c_vars_req
        character*180 c_values_req

        character*3 lso_ext        
        data lso_ext /'lso'/

        allocate( cldcv1(nx_l,ny_l,kcloud), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate cldcv1'
        endif

        allocate( cf_modelfg(nx_l,ny_l,kcloud), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate cf_modelfg'
        endif

        allocate( t_modelfg(nx_l,ny_l,kcloud), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate t_modelfg'
        endif

        allocate( sh_modelfg(nx_l,ny_l,nz_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate sh_modelfg'
        endif

        allocate( ref_modelfg(nx_l,ny_l,nz_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate ref_modelfg'
        endif

        allocate( cld_snd(max_cld_snd,kcloud), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate cld_snd'
        endif

        allocate( wt_snd(max_cld_snd,kcloud), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate wt_snd'
        endif

        istat = init_timer()

        write(6,*)' welcome to the laps gridded cloud analysis'

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error calling get_r_missing_data'
           stop
        endif

        lstat_radar_3dref_orig_a = .false.
        cloud_albedo = r_missing_data ! initialize

        call get_i_perimeter(i_perimeter,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error calling get_i_perimeter'
           stop
        endif

        ix_low  = 1    - i_perimeter
        ix_high = nx_l + i_perimeter
        iy_low  = 1    - i_perimeter
        iy_high = ny_l + i_perimeter

        nx_dim_lut = nx_l + i_perimeter - 1
        ny_dim_lut = ny_l + i_perimeter - 1

        default_base     = r_missing_data
        default_top      = r_missing_data
        default_ceiling  = r_missing_data

        do j = 1,ny_l
           do i = 1,nx_l
              c1_name_array(i,j,:) = ' '
           enddo
        enddo

        call get_ref_base(ref_base,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting ref_base'
           stop
        endif

!     this should work out to be slightly larger than needed
      n_fnorm =
     1      1.6 * ( (nx_dim_lut*nx_dim_lut) + (ny_dim_lut*ny_dim_lut) )


c determine the source of the radar data
        c_vars_req = 'radarext_3d'

        call get_static_info(c_vars_req,c_values_req,1,istatus)

        if(istatus .eq. 1)then
            radarext_3d_cloud = c_values_req(1:3)
        else
            write(6,*)' error getting radarext_3d'
            return
        endif

        write(6,*)' radarext_3d_cloud = ',radarext_3d_cloud

c read in laps lat/lon and topo
        call get_laps_domain_95(nx_l,ny_l,lat,lon,topo
     1           ,rlaps_land_frac,grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting laps domain'
            return
        endif

        write(6,*)' actual grid spacing in domain center = '
     1                              ,grid_spacing_cen_m

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' error getting laps_cycle_time'
            return
        endif


!       fill namelist_parms data structure
        namelist_parms%l_use_metars = l_use_metars
        namelist_parms%l_use_radar = l_use_radar

        if(latency_co2 .ge. 0)then
            l_use_co2_mode1 = .true.
        else
            l_use_co2_mode1 = .false.
        endif

        l_use_co2_mode2 = .false.

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

        n_prods = 4
        iprod_start = 1
        iprod_end = n_prods

        ext = 'lc3'

c dynamically adjust the heights of the cloud levels to the terrain
        write(6,*)' fitting cloud heights to the terrain'
        topo_min = 1e10
        do j = 1,ny_l
        do i = 1,nx_l
            topo_min = min(topo_min,topo(i,j))
        enddo
        enddo

        topo_min_buff = topo_min - .005 ! put in a slight buffer so that
                                        ! cloud height grid extends below topo

        write(6,*)'      old ht     new ht      lowest terrain point = '       
     1           ,topo_min

        range_orig = cld_hts(kcloud) - cld_hts(1)
        range_new =  cld_hts(kcloud) - topo_min_buff
        range_ratio = range_new / range_orig

        do k = 1,kcloud
            cld_hts_orig = cld_hts(k)
            cld_hts(k) = cld_hts(kcloud)
     1                - (cld_hts(kcloud)-cld_hts_orig) * range_ratio
            write(6,21)k,cld_hts_orig,cld_hts(k)
21          format(1x,i3,2f10.3)
        enddo ! k

        if(cld_hts(1) .ge. topo_min)then
            write(6,*)' error: topo extends at or below edge of '
     1               ,'cloud height grid',topo_min,cld_hts(1)
            return
        endif

        write(6,*)' initializing fields'
c initialize the fields
        if(l_perimeter)then
          do k = 1,kcloud
          do n = 1,max_cld_snd
              cld_snd(n,k) = r_missing_data
          enddo
          enddo
        else
          do k = 1,kcloud
          do j = 1,ny_l
          do i = 1,nx_l
            cldcv1(i,j,k) = r_missing_data
          enddo
          enddo
          enddo
        endif

c read in laps temps
        var = 't3'
        ext = 'lt1'
        call get_laps_3d(i4time,nx_l,ny_l,nz_l
     1  ,ext,var,units,comment,temp_3d,istatus)

        if(istatus .ne. 1)then
            write(6,*)' error reading 3d temps - get_modelfg not called'
            goto999
        endif

c read in laps heights
        var = 'ht'
        ext = 'lt1'
        call get_laps_3d(i4time,nx_l,ny_l,nz_l
     1  ,ext,var,units,comment,heights_3d,istatus)

        if(istatus .ne. 1)then
            write(6,*)' error reading 3d heights - get_modelfg not calle
     1d'
            goto999
        endif

c obtain model first guess cloud cover field (along with reflectivity, water vapor)
        call get_modelfg(cf_modelfg,t_modelfg,sh_modelfg,ref_modelfg    ! o
     1           ,default_clear_cover,r_missing_data                    ! i
     1           ,temp_3d,heights_3d,cld_hts                            ! i
     1              ,i4time,ilaps_cycle_time                            ! i
     1                  ,nx_l,ny_l,nz_l,kcloud                          ! i
     1                  ,istatus)                                       ! o

c read in radar data
!       get time of radar file of the indicated appropriate extension
        call get_filespec(radarext_3d_cloud(1:3),2,c_filespec,istatus)
        call get_file_time(c_filespec,i4time,i4time_radar)

        if(namelist_parms%l_use_radar)then
            i4_tol = 1200
        else
            write(6,*)' withholding radar data from cloud analysis'
            i4_tol = -1
        endif

        call read_multiradar_3dref(i4time,                               ! i
     1                 i4_tol,i4_ret,                                    ! i
     1                 .true.,ref_base,                                  ! i
     1                 nx_l,ny_l,nz_l,radarext_3d_cloud,                 ! i
     1                 lat,lon,topo,.true.,.true.,                       ! i
     1                 heights_3d,                                       ! i
     1                 radar_ref_3d,                                     ! o
     1                 rlat_radar,rlon_radar,rheight_radar,radar_name,   ! o
     1                 iqc_2dref,closest_radar,                          ! o
     1                 n_ref_grids,n_2dref,n_3dref,istat_radar_2dref_a,  ! o  
     1                 istat_radar_3dref_a)                              ! o

        if(iqc_2dref .eq. 1)then
            write(6,*)' good quality 2d radar (e.g. nowrad)'
            l_trust_narrowband = .true.
        else
            write(6,*)' lesser quality 2d radar (e.g. narrowband)'      
        endif

        rqc_2d = 0.
        where (istat_radar_2dref_a .eq. 1)rqc_2d = 2.
        where (istat_radar_3dref_a .eq. 1)rqc_2d = 3.

!       else
!           write(6,*)'radar data outside time window'
!           n_ref_grids = 0
!           call constant_i(istat_radar_2dref_a,0,nx_l,ny_l)       
!           call constant_i(istat_radar_3dref_a,0,nx_l,ny_l)

!       endif

c blend in first guess radar
        n_fg_radar = 0
        n_fg_echoes = 0

        call get_maxtops(ref_modelfg,heights_3d
     1                  ,nx_l,ny_l,nz_l,ref_mt_modelfg)

        write(6,*)' modelfg maxtops range is: '
     1      ,minval(ref_mt_modelfg),maxval(ref_mt_modelfg)
        write(6,*)' closest radar range is: '
     1      ,minval(closest_radar),maxval(closest_radar)

        range_thresh_lo = 200000.
        range_thresh_hi = 300000.

        radius_earth_8_thirds = 6371.e3 * 2.6666666
        aterm = 1. / radius_earth_8_thirds

        elev_ang_thr = 0.5 ! angle of radar horizon 
        bterm = tand(elev_ang_thr)

        do i = 1,nx_l
        do j = 1,ny_l

!           determine dynamic cutoff distance between radar data and first guess
            if(ref_mt_modelfg(i,j) .gt. 0.)then ! find distance from radar given
                cterm = ref_mt_modelfg(i,j)     ! echo height and elev angle
                hor_dist = (sqrt(4.*aterm*cterm + bterm**2.) - bterm)       
     1                                / (2.*aterm)
                range_thresh = 
     1              min(max(hor_dist,range_thresh_lo),range_thresh_hi)
            else
                range_thresh = range_thresh_hi
            endif

            if(rqc_2d(i,j) .eq. 0. .or. 
     1         closest_radar(i,j) .gt. range_thresh)then
                do k = 1,nz_l
                    if(ref_modelfg(i,j,k) .ne. r_missing_data)then
                        radar_ref_3d(i,j,k) = ref_modelfg(i,j,k)
                        rqc_2d(i,j) = 1.
                        n_fg_radar = n_fg_radar + 1
                        if(ref_modelfg(i,j,k) .gt. ref_base)then
                            n_fg_echoes = n_fg_echoes + 1
                        endif
                    endif ! first guess is present
                enddo ! k
            endif ! radar obs data are absent
        enddo ! j
        enddo ! i

        frac_fg_radar = float(n_fg_radar) / float(nx_l*ny_l*nz_l)

        write(6,*)'first guess radar used over ',frac_fg_radar*100.
     1           ,'% of domain'

        write(6,*)'number of first guess echoes is ',n_fg_echoes

c read in and insert sao data as cloud soundings
!       read in surface pressure
        var = 'ps'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var,units,comment
     1                  ,nx_l,ny_l,pres_sfc_pa,istatus)
        if(istatus .ne. 1)then
            write(6,*)
     1  ' error reading surface pres analyses - abort cloud analysis'
            goto999
        endif

        var = 't'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var,units,comment
     1                  ,nx_l,ny_l,t_sfc_k,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading sfc temps - abort cloud analysis'
            goto999
        endif

        var = 'td'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var,units,comment
     1                  ,nx_l,ny_l,td_sfc_k,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading sfc td - abort cloud analysis'
            goto999
        endif

        var = 'tgd'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var,units,comment
     1                  ,nx_l,ny_l,tgd_sfc_k,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading sfc tgd - abort cloud analysis'
            goto999
        endif

        write(6,*)
        write(6,*)' call ingest/insert sao routines'
        n_cld_snd = 0
        call insert_sao(i4time,cldcv1,cf_modelfg,t_modelfg             ! i
     1  ,cld_hts,default_clear_cover,namelist_parms                    ! i
     1  ,lat,lon,topo,t_sfc_k,wtcldcv                                  ! i
     1  ,c1_name_array,l_perimeter,ista_snd
     1  ,cvr_snd,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd
     1  ,nx_l,ny_l,kcloud                                              ! i
     1  ,n_obs_pos_b,lat_s,lon_s,c_stations    ! returned for precip type comp
     1  ,wx_s,t_s,td_s                         !    "      "    "     "
     1  ,elev_s                                ! ret for comparisons
     1  ,rad_s,solar_ea                        !  "   "       "   
     1  ,istat_sfc,maxstns,ix_low,ix_high,iy_low,iy_high)
   
        
        if(istat_sfc .ne. 1)then
            write(6,*)' no sao data inserted: aborting cloud analysis'
            goto999
        endif


c read in and insert pirep data as cloud soundings
        if(l_use_pireps .eqv. .true.)then
          write(6,*)' using pireps'

          call insert_pireps(i4time,cld_hts
     1      ,default_clear_cover
     1      ,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd
     1      ,lat,lon,nx_l,ny_l,kcloud,ix_low,ix_high,iy_low,iy_high
     1      ,n_pirep,istatus)
          if(istatus .ne. 1)then
            write(6,*)' error: bad status from insert_pireps,'
     1               ,' aborting cloud analysis'
            goto999
          endif
          i4_elapsed = ishow_timer()
        else
          write(6,*)' withholding pireps'
        endif

c read in and insert co2 slicing data as cloud soundings
        if(l_use_co2_mode1)then
            call insert_co2ctp(i4time,cld_hts,heights_3d                  ! i
     1            ,nx_l,ny_l,nz_l,kcloud,r_missing_data                   ! i
     1            ,l_use_co2_mode1,latency_co2                            ! i
     1            ,default_clear_cover                                    ! i
     1            ,lat,lon,ix_low,ix_high,iy_low,iy_high                  ! i
     1            ,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd       ! i/o
     1            ,istatus)                                               ! o
        endif

c do analysis to horizontally spread sao, pirep, and optionally co2 data
        write(6,*)
        if(l_use_co2_mode1)then
            write(6,*)' analyzing sfc obs, pirep, and co2-slicing data'       
        else
            write(6,*)' analyzing sfc obs and pirep data'
        endif

        max_obs = n_cld_snd * kcloud

!       set weight for using model background clouds beyond a certain effective
!       radius of influence from the sfc obs/pireps
!       weight_modelfg = 0.    ! model wt inactive, obs used to infinite radius
!       weight_modelfg = 1.    ! model used beyond ~100km from nearest obs
!       weight_modelfg = .01   ! model used beyond ~250km from nearest obs
!       weight_modelfg = .0001 ! model used beyond ~630km from nearest obs
 
        weight_modelfg = cld_weight_modelfg 

        call barnes_r5(clouds_3d,nx_l,ny_l,kcloud,cldcv1,wtcldcv
     1     ,cf_modelfg,l_perimeter,cld_snd,wt_snd,r_missing_data
     1     ,grid_spacing_cen_m,i_snd,j_snd,n_cld_snd,max_cld_snd
     1     ,max_obs,weight_modelfg,nx_dim_lut,ny_dim_lut
     1     ,ix_low,ix_high,iy_low,iy_high,n_fnorm,istatus) 
        if(istatus .ne. 1)then
            write(6,*)
     1      ' error: bad status from barnes_r5, aborting cloud analysis'
            goto999
        endif

!       use model first guess in upper part of domain (or use zero)
        ktop = 22 ! restrict where analysis occurs (helps with ir/vis)

!       no restriction where analysis occurs (mainly for ir)
!       ktop = kcloud 

        if(ktop .lt. kcloud)then
            write(6,*)' use model first guess w/o barnes above k = '
     1               ,ktop
            clouds_3d(:,:,ktop+1:kcloud) = cf_modelfg(:,:,ktop+1:kcloud)
!           write(6,*)' use zero first guess above k =  ',ktop
!           clouds_3d(:,:,ktop+1:kcloud) = 0.          
        endif

!       cloud cover qc check
        call qc_clouds_3d(clouds_3d,nx_l,ny_l,kcloud)

!       hold current cloud cover from sao/pireps in a buffer array
        do k = 1,kcloud
        do j = 1,ny_l
        do i = 1,nx_l
            cldcv_sao(i,j,k) = clouds_3d(i,j,k)
        enddo
        enddo
        enddo

        var = 'sc'
        ext = 'lm2'
        call get_laps_2d(i4time-ilaps_cycle_time,ext,var,units,comment
     1                  ,nx_l,ny_l,cvr_snow,istat_cvr_snow)

        if(istat_cvr_snow .ne. 1)then               ! try using snow cover field
!       if(istat_cvr_snow .ne. 1 .or. .true.)then   ! don't use snow cover
            write(6,*)' error reading snow cover'
            do i = 1,nx_l
            do j = 1,ny_l
                cvr_snow(i,j) = r_missing_data
            enddo
            enddo
        endif

        i4_elapsed = ishow_timer()

c read in satellite data
        call get_sat_data(i4time,i4_sat_window,i4_sat_window_offset,     ! i
     1                    nx_l,ny_l,r_missing_data,                      ! i
     1                    l_use_39,l_use_co2_mode2,latency_co2,          ! i
     1                    lat,lon,                                       ! i
     1                    subpoint_lat_clo_s8a,subpoint_lon_clo_s8a,     ! o
     1                    tb8_k,istat_tb8,comment_tb8,                   ! o
     1                    t39_k,istat_t39,comment_t39,                   ! o
     1                    sst_k,istat_sst,comment_sst,                   ! o
     1                    cldtop_co2_pa_a,cloud_frac_co2_a,              ! o
     1                    istat_co2,lstat_co2_a)                         ! o

        i4_elapsed = ishow_timer()

!       default (e.g. for goes-16, coms)
        offset_ir_i(:,:) = 0.
        offset_ir_j(:,:) = 0.

!       positive i offset will shift the image west        
!       positive j offset will shift the image south

        where(abs(subpoint_lon_clo_s8a(:,:)-( -75.)) .lt.  8.) ! goes 16
            offset_ir_i(:,:)  = +6500./grid_spacing_cen_m
            offset_ir_j(:,:)  = -5500./grid_spacing_cen_m
!           offset_ir_i(:,:)  =  -500./grid_spacing_cen_m
!           offset_ir_j(:,:)  = +2500./grid_spacing_cen_m
            offset_vis_i(:,:) = -1500./grid_spacing_cen_m
            offset_vis_j(:,:) = +3000./grid_spacing_cen_m
            offset_vis_i(:,:) = -3000./grid_spacing_cen_m
            offset_vis_j(:,:) = +4000./grid_spacing_cen_m
            offset_vis_i(:,:) = +2000./grid_spacing_cen_m
            offset_vis_j(:,:) = -2000./grid_spacing_cen_m
        end where
        where(abs(subpoint_lon_clo_s8a(:,:)-(-135.)) .lt. 20.) ! goes w
            offset_ir_i(:,:)  = -3000./grid_spacing_cen_m
            offset_ir_j(:,:)  =     0.
            offset_vis_i(:,:) =     0./grid_spacing_cen_m
            offset_vis_j(:,:) =     0./grid_spacing_cen_m
        end where
!       where(abs(subpoint_lon_clo_s8a(:,:)-( -75.)) .lt.  8.) ! goes e
!           offset_ir_i(:,:)  = +4700./grid_spacing_cen_m
!           offset_ir_j(:,:)  = -3200./grid_spacing_cen_m
!           offset_vis_i(:,:) =     0./grid_spacing_cen_m
!           offset_vis_j(:,:) =     0./grid_spacing_cen_m
!       end where

!       apply offsets to 11u ir data
        do i = 1,nx_l
        do j = 1,ny_l
            ioff = min(max(i+nint(offset_ir_i(i,j)),1),nx_l)
            joff = min(max(j+nint(offset_ir_j(i,j)),1),ny_l)
            buff(i,j) = tb8_k(ioff,joff)
        enddo ! j
        enddo ! i
        tb8_k(:,:) = buff(:,:)

!       apply offsets to 3.9u data        
        do i = 1,nx_l
        do j = 1,ny_l
            ioff = min(max(i+nint(offset_ir_i(i,j)),1),nx_l)
            joff = min(max(j+nint(offset_ir_j(i,j)),1),ny_l)
            buff(i,j) = t39_k(ioff,joff)
        enddo ! j
        enddo ! i
        t39_k(:,:) = buff(:,:)

!       calculate solar altitude
        do j = 1,ny_l
        do i = 1,nx_l
            call solar_position(lat(i,j),lon(i,j),i4time,solar_alt(i,j)
     1                                     ,solar_dec,solar_ha(i,j))
            call equ_to_altaz_d(solar_dec,solar_ha(i,j),lat(i,j)         ! i
     1                                   ,altdum,solar_az(i,j))          ! o
        enddo
        enddo

!       yescloud
        idb = 515
        jdb = 185

!       nocloud
        idb = 516
        jdb = 188

!       nocloud
        idb = 516
        jdb = 188

!       snow clouds (mtns)
        idb = 440
        jdb = 323

!       domain center
        idb = (nx_l/2) + 1
        jdb = (ny_l/2) + 1

!       snow clouds (plains)
        idb = 850 ! 591
        jdb =  30 ! 511

        icen = idb
        jcen = jdb

        write(6,*)' ctr idb/jdb = ',idb,jdb

        write(6,*)'solar dec/ha ',solar_dec,solar_ha(icen,jcen),
     1            ' at lat/lon ',lat(icen,jcen),lon(icen,jcen)
        write(6,*)'solar altitude = ',solar_alt(icen,jcen)

!       cloud cover qc check
        call qc_clouds_3d(clouds_3d,nx_l,ny_l,kcloud)

        write(6,*)' call get_parallax_info before insert_sat (ir)'
        call get_parallax_info(nx_l,ny_l,i4time                          ! i
     1                        ,lat,lon                                   ! i
     1                        ,subpoint_lat_clo_s8a,subpoint_lon_clo_s8a ! i
     1                        ,di_dh_ir,dj_dh_ir,i_fill_seams)           ! o

!       consider passing out parallax info, used with vis & cvr_snow, to be
!       later used in other ways such as ir
!       returns 'cloud_frac_vis_a' (cloud albedo with clear alb subtracted)
!       returns 'sat_albedo' (cloud albedo from reflectance/sfc alb, though)
!       sfc alb isn't yet subtracted). 'sat_albedo' isn't yet used in analysis

        call get_vis(i4time,solar_alt,l_use_vis,l_use_vis_add            ! i
     1              ,l_use_vis_partial,lat,lon,idb,jdb                   ! i
     1              ,i4_sat_window,i4_sat_window_offset                  ! i
     1              ,rlaps_land_frac,topo                                ! i
     1              ,cvr_snow,tgd_sfc_k                                  ! i
     1              ,offset_vis_i,offset_vis_j                           ! i
     1              ,di_dh_vis,dj_dh_vis                                 ! o
     1              ,cloud_frac_vis_a,sat_albedo,sat_refl,mode_refl      ! o
     1              ,ihist_alb,static_albedo,sfc_albedo                  ! o
     1              ,vis_snow_max                                        ! o
     1              ,subpoint_lat_clo_vis,subpoint_lon_clo_vis           ! o 
     1              ,comment_alb                                         ! o
     1              ,nx_l,ny_l,kcloud,r_missing_data                     ! o
     1              ,istat_vis_potl_a,istat_vis)                         ! o

        call get_istat_39(t39_k,tb8_k,solar_alt,r_missing_data           ! i
     1                   ,rlaps_land_frac,nx_l,ny_l                      ! i
     1                   ,static_albedo                                  ! i
     1                   ,istat_39_a)                                    ! o

!       cloud cover qc check
        call qc_clouds_3d(clouds_3d,nx_l,ny_l,kcloud)

        call get_pres_3d(i4time,nx_l,ny_l,nz_l,pres_3d,istatus)

        if(i_varadj .eq. 1)then
            write(6,*)' call cloud_var before insert_sat'
            t_gnd_k = t_sfc_k ! initialize with preliminary value
            call cloud_var(i4time,lat,lon
     1                    ,nx_l,ny_l,nz_l,kcloud,heights_3d,temp_3d
     1                    ,t_gnd_k,clouds_3d,cld_hts,tb8_k
     1                    ,cloud_frac_vis_a
     1                    ,subpoint_lat_clo_s8a,subpoint_lon_clo_s8a  ! i 
     1                    ,r_missing_data                             ! i
     1                    ,di_dh_ir,dj_dh_ir)                         ! i 
        endif

        if(l_corr_parallax .eqv. .false.)then
            write(6,*)' turn off parallax correction for remainder'
            di_dh_ir = 0. ! for testing
            dj_dh_ir = 0. ! for testing
            di_dh_vis = 0. ! for testing
            dj_dh_vis = 0. ! for testing
        endif

        i4_elapsed = ishow_timer()

        write(6,*)' offset_ir',subpoint_lon_clo_s8a(idb,jdb)
     1                        ,offset_ir_i(idb,jdb)
     1                        ,offset_ir_j(idb,jdb)
        write(6,*)' offset_vis',subpoint_lon_clo_s8a(idb,jdb)
     1                        ,offset_vis_i(idb,jdb)
     1                        ,offset_vis_j(idb,jdb)
        write(6,*)' dij_dh_ir',di_dh_ir(idb,jdb)
     1                        ,dj_dh_ir(idb,jdb)
        write(6,*)' dij_dh_vis',di_dh_vis(idb,jdb)
     1                         ,dj_dh_vis(idb,jdb)

        call insert_sat(i4time,clouds_3d,cldcv_sao,cld_hts,lat,lon,
     1       pct_req_lvd_s8a,default_clear_cover,                       ! i
     1       tb8_k,istat_tb8,                                           ! i
     1       sst_k,istat_sst,                                           ! i
     1       istat_39_a, l_use_39,                                      ! i
     1       di_dh_ir,dj_dh_ir,                                         ! i
     1       di_dh_vis,dj_dh_vis,                                       ! i
     1       i_fill_seams,idb,jdb,                                      ! i
     1       offset_ir_i,offset_ir_j,                                   ! i
     1       istat_39_add_a,                                            ! o
     1       tb8_cold_k,                                                ! o
     1       grid_spacing_cen_m,surface_sao_buffer,                     ! i
     1       cloud_frac_vis_a,istat_vis_potl_a,                         ! i
     1       istat_vis_added_a,                                         ! o
     1       solar_alt,solar_ha,solar_dec,                              ! i
     1       lstat_co2_a, cloud_frac_co2_a, cldtop_co2_pa_a,            ! i
     1       rlaps_land_frac,                                           ! i
     1       topo,heights_3d,temp_3d,t_sfc_k,td_sfc_k,pres_sfc_pa,      ! i
     1       t_modelfg,sh_modelfg,pres_3d,                              ! i
     1       cvr_snow,nx_l,ny_l,kcloud,nz_l,r_missing_data,sfc_albedo,  ! i
     1       t_gnd_k,                                                   ! o
     1       cldtop_co2_m,cldtop_tb8_m,cldtop_m,ht_sao_top,             ! o
     1       cldht_prlx_top,cldht_prlx_unsm,                            ! o
     1       istatus)                                                   ! o

        if(istatus .ne. 1)then
            write(6,*)' error: bad status returned from insert_sat'
            goto999
        endif

        write(6,291)(heights_3d(idb,jdb,k)
     1              ,clouds_3d(idb,jdb,k),k=kcloud,1,-1)
291     format(' cldcv section 1b (after insert_sat):'
     1                  /'    ht      cvr',50(/f8.1,f8.3,'   ctr1b'))

!       cloud cover qc check
        call qc_clouds_3d(clouds_3d,nx_l,ny_l,kcloud)

        write(6,*)' cloud top (band 8 vs. ht_sao_top)'
        scale = .0001
        write(6,301)
301     format('  cloud top (km)             band 8                ',
     1                 23x,'               ht_sao_top')
        call array_plot(cldtop_tb8_m,ht_sao_top,nx_l,ny_l,'horz cv'
     1                 ,c1_name_array(:,:,1),kcloud,cld_hts,scale)

        write(6,*)' cloud top (band 8 vs. satellite analysis)'
        scale = .0001
        write(6,302)
302     format('  cloud top (km)             band 8                ',
     1                 23x,'          satellite analysis')
        call array_plot(cldtop_tb8_m,cldtop_m,nx_l,ny_l,'horz cv'
     1                 ,c1_name_array(:,:,1),kcloud,cld_hts,scale)

        i4_elapsed = ishow_timer()

        n_radar_2dref=0
        n_radar_3dref=0
        n_radar_3dref_orig=0

c       three dimensionalize radar data if necessary (e.g. nowrad)
        write(6,*)' three dimensionalizing radar data'

!       clear out radar echo above the highest cloud top
        k_ref_def = nint(zcoord_of_pressure(float(700*100)))

        do j = 1,ny_l
        do i = 1,nx_l
            
            if(istat_radar_2dref_a(i,j) .eq. 1 .and.
     1         istat_radar_3dref_a(i,j) .eq. 0       )then

                if(radar_ref_3d(i,j,1) .gt. ref_base)then

                    k_topo = int(zcoord_of_pressure(pres_sfc_pa(i,j)))
                    k_ref = max(k_ref_def,k_topo+2)

!                   k_ref = k_ref_def

                    cloud_top_m = default_top

!                   test for cloud top
                    do k = kcloud-1,1,-1
                        if(clouds_3d(i,j,k  ) .gt. thresh_cvr_ceiling 
     1                                       .and.       
     1                     clouds_3d(i,j,k+1) .le. thresh_cvr_ceiling 
     1                                                            )then 
                            cloud_top_m = 0.5 * (cld_hts(k) 
     1                                         + cld_hts(k+1))
                            if(cloud_top_m .gt. heights_3d(i,j,nz_l)
     1                                                            )then
                                cloud_top_m = heights_3d(i,j,nz_l)
                            endif
                           goto150
                        endif
                    enddo ! k

!                   add third dimension to radar echo
150                 if(abs(cloud_top_m) .le. 1e6)then ! valid cloudtop
                        k_cloud_top = nint(height_to_zcoord2(cloud_top_m
     1                          ,heights_3d,nx_l,ny_l,nz_l,i,j,istatus))      
                        if(istatus .ne. 1)then
                            write(6,*)' error: bad status returned'
     1                          ,' from height_to_zcoord2'
     1                          ,cloud_top_m,heights_3d(i,j,nz_l),i,j       
                            goto999
                        endif

                        k_cloud_top = max(k_ref,k_cloud_top)
                    else
                        if(l_trust_narrowband)then
                            k_cloud_top = k_ref
                        else
                            k_cloud_top = 1 ! will set radar to no echo
                        endif
                    endif

                    do k = nz_l,k_cloud_top,-1
                        radar_ref_3d(i,j,k) = ref_base
                    enddo ! k

                endif ! radar echo at this grid point

!               if(k_cloud_top .gt. 1)then
!               we have three-dimensionalized this grid point
                n_radar_2dref = n_radar_2dref + 1
                n_radar_3dref = n_radar_3dref + 1
                istat_radar_3dref_a(i,j) = 1
!               endif

            elseif(istat_radar_2dref_a(i,j) .eq. 1 .and.
     1             istat_radar_3dref_a(i,j) .eq. 1       )then

!               grid point is already fully three dimensional
                n_radar_2dref = n_radar_2dref + 1
                n_radar_3dref = n_radar_3dref + 1
                n_radar_3dref_orig = n_radar_3dref_orig + 1
                lstat_radar_3dref_orig_a(i,j) = .true.

            endif ! is this grid point 2-d or 3-d?

        enddo ! i
        enddo ! j

        if(n_radar_2dref .ge. 1)then
            istat_radar_2dref = 1
        else
            istat_radar_2dref = 0
        endif

        if(n_radar_3dref .ge. 1)then
            istat_radar_3dref = 1
        else
            istat_radar_3dref = 0
        endif

!       note, if orig data is a mixture of 2d and 3d, the istat here gets set
!       to 3d.
        if(n_radar_3dref_orig .ge. 1)then
            istat_radar_3dref_orig = 1
        else
            istat_radar_3dref_orig = 0
        endif

        write(6,*)' n_radar 2dref/3dref_orig/3dref = ' 
     1           ,n_radar_2dref,n_radar_3dref_orig,n_radar_3dref

        write(6,*)' istat_radar 2dref/3dref_orig/3dref = ' 
     1           ,istat_radar_2dref,istat_radar_3dref_orig
     1           ,istat_radar_3dref

        i4_elapsed = ishow_timer()

!       generate field of radar coverage (2d/3d)
        do i = 1,nx_l
        do j = 1,ny_l
            if(rqc_2d(i,j) .eq. 1.)then
                plot_maskr(i,j) = 10.0 ! 3d first guess data
            elseif(lstat_radar_3dref_orig_a(i,j))then
                plot_maskr(i,j) = 30.0 ! 3d original data
            elseif(istat_radar_2dref_a(i,j) .eq. 1)then
                plot_maskr(i,j) = 20.0 ! 2d original data
            else                    
                plot_maskr(i,j) = 0.0  ! no data
            endif
        enddo ! j
        enddo ! i

!       write(6,*)' cldcv section 1c:',clouds_3d(idb,jdb,:)

c insert radar data
        if(n_radar_3dref .gt. 0)then
            call get_max_reflect(radar_ref_3d,nx_l,ny_l,nz_l,ref_base       
     1                          ,dbz_max_2d)
 
!           generate ascii prelim reflectivity plot
            write(6,1201)
1201        format('                  radar coverage',50x
     1                              ,'prelim radar max reflectivity')       
            scale = 0.01
            call array_plot(plot_maskr,dbz_max_2d,nx_l,ny_l,'horz cv'
     1                     ,c1_name_array(:,:,1),kcloud,cld_hts,scale) ! radar 

            call insert_radar(i4time,clouds_3d,cld_hts
     1          ,temp_3d,t_sfc_k,td_sfc_k                            ! i
     1          ,grid_spacing_cen_m,nx_l,ny_l,nz_l                   ! i
     1          ,kcloud,cloud_base,ref_base                          ! i
     1          ,topo,solar_alt,r_missing_data                       ! i
     1          ,radar_ref_3d,dbz_max_2d                             ! i/o
     1          ,vis_radar_thresh_dbz                                ! i
     1          ,l_unresolved                                        ! o
     1          ,heights_3d                                          ! i
     1          ,istatus)                                            ! o

            if(istatus .ne. 1)then
                write(6,*)
     1          ' error: bad status returned from insert_radar'      
                goto999
            endif

        endif

        write(6,391)(heights_3d(idb,jdb,k)
     1              ,clouds_3d(idb,jdb,k),k=kcloud,1,-1)
391     format(' cldcv section 2 (after insert_radar):'
     1                  /'    ht      cvr',50(/f8.1,f8.3,'   ctr2'))

        if(i_varadj .eq. 1)then
            write(6,*)' call cloud_var before insert vis'
            call cloud_var(i4time,lat,lon
     1                    ,nx_l,ny_l,nz_l,kcloud,heights_3d,temp_3d
     1                    ,t_gnd_k,clouds_3d,cld_hts,tb8_k
     1                    ,cloud_frac_vis_a
     1                    ,subpoint_lat_clo_vis,subpoint_lon_clo_vis  ! i 
     1                    ,r_missing_data                             ! i
     1                    ,di_dh_vis,dj_dh_vis)                       ! i
        endif

        i4_elapsed = ishow_timer()

c       insert visible / 3.9u satellite in clearing step
        if(istat_vis .eq. 1 .or. (istat_t39 .eq. 1 .and. l_use_39) )then
!       if(.false.)then

!           this routine will set 'cloud_albedo' from 'cloud_frac_vis_a'
!           'sat_albedo isn't yet used but could be used for 'mode_refl=1'
!           once it is refined with 'sfc_albedo' adjustment. 'cloud_albedo' 
!           is now parallax corrected.
           

            call insert_vis(i4time,clouds_3d,cld_hts
     1        ,topo,cloud_frac_vis_a,sat_albedo,mode_refl,ihist_alb   ! i
     1        ,istat_39_a,l_use_39,idb,jdb                            ! i
     1        ,nx_l,ny_l,kcloud,r_missing_data                        ! i
     1        ,vis_radar_thresh_cvr,vis_radar_thresh_dbz              ! i
     1        ,istat_radar_3dref,radar_ref_3d,nz_l,ref_base
     1        ,solar_alt,solar_az                                     ! i
     1        ,di_dh_vis,dj_dh_vis,i_fill_seams                       ! i
     1        ,cldtop_tb8_m,cldht_prlx_top                            ! i
     1        ,cloud_albedo,cloud_od,cloud_op                         ! o
     1        ,dbz_max_2d,surface_sao_buffer,istatus)

        else
            write(6,'(" ctr3 skipping call to insert_vis")')
            write(6,'(" ctr3 cloud albedo",g12.4)')cloud_albedo(idb,jdb)
        endif

        write(6,491)(heights_3d(idb,jdb,k)
     1              ,clouds_3d(idb,jdb,k),k=kcloud,1,-1)
491     format(' cldcv section 3 (after insert_vis):'
     1        /'    ht      cvr',50(/f8.1,f8.3,'   ctr3'))

c clear out stuff below ground
        do k = 1,kcloud
        do j = 1,ny_l
        do i = 1,nx_l
            if(cld_hts(k) .lt. topo(i,j))
     1                   clouds_3d(i,j,k) = default_clear_cover
        enddo
        enddo
        enddo

c ascii plots in horizontal and vertical slices

c       horizontal slices

        do k=4,kcloud,2
            call slice(cf_modelfg,nx_l,ny_l,kcloud,cvhz1
     1                ,nx_l,ny_l,1,0,0,k,0,0)
            call slice(cldcv_sao,nx_l,ny_l,kcloud,cvr_max
     1                ,nx_l,ny_l,1,0,0,k,0,0)
            write(6,401)k,cld_hts(k)
401         format(4x,'lvl',i4,f8.0,' m     model first guess only',
     1            20x,'              with point data added')
            scale = 1.
            call array_plot(cvhz1,cvr_max,nx_l,ny_l,'horz cv'
     1                     ,c1_name_array(:,:,k),kcloud,cld_hts,scale)
        enddo ! k

        do k=4,kcloud,4
            call slice(cldcv_sao,nx_l,ny_l,kcloud,cvhz1
     1                ,nx_l,ny_l,1,0,0,k,0,0)
            call slice(clouds_3d,nx_l,ny_l,kcloud,cvr_max
     1                ,nx_l,ny_l,1,0,0,k,0,0)
            write(6,402)k,cld_hts(k)
402         format(4x,'lvl',i4,f8.0,' m     before satellite/radar',
     1            20x,'              after satellite/radar')
            scale = 1.
            call array_plot(cvhz1,cvr_max,nx_l,ny_l,'horz cv'
     1                     ,c1_name_array(:,:,k),kcloud,cld_hts,scale)
        enddo ! k

c       ew slices
        write(6,*)
        do j=10,ny_l,20
            call slice(cf_modelfg,nx_l,ny_l,kcloud,cvew1
     1                ,nx_l,kcloud,0,1,0,0,j,0)
            call slice(cldcv_sao,nx_l,ny_l,kcloud,cvew2
     1                ,nx_l,kcloud,0,1,0,0,j,0)
            write(6,501)j
501         format(5x,'  j =',i4,10x,'  model first guess only       ',
     1          21x,'           with point data added')

            scale = 1.
            call array_plot(cvew1,cvew2,nx_l,kcloud,'vert cv'
     1                     ,c1_name_array,kcloud,cld_hts,scale)
        enddo ! j

c       ew slices
        write(6,*)
        do j=10,ny_l,20
            call slice(cldcv_sao,nx_l,ny_l,kcloud,cvew1
     1                ,nx_l,kcloud,0,1,0,0,j,0)
            call slice(clouds_3d,nx_l,ny_l,kcloud,cvew2
     1                ,nx_l,kcloud,0,1,0,0,j,0)
            write(6,502)j
502         format(5x,'  j =',i4,10x,'  before satellite/radar       ',
     1          21x,'           after satellite/radar')

            scale = 1.
            call array_plot(cvew1,cvew2,nx_l,kcloud,'vert cv'
     1                     ,c1_name_array,kcloud,cld_hts,scale)
        enddo ! j

!       get max cloud cover
        do j = 1,ny_l
        do i = 1,nx_l
            cvr_sao_max(i,j) = 0.
            cvr_max(i,j) = 0.
            do k = 1,kcloud
                cvr_sao_max(i,j) = max(cvr_sao_max(i,j)
     1                                ,cldcv_sao(i,j,k))   
                cvr_max(i,j) = max(cvr_max(i,j),clouds_3d(i,j,k))
            enddo ! k
        enddo ! i
        enddo ! j

!       log info on cloud holes with cover at least .10 less than all neighbors
        do i = 2,nx_l-1
        do j = 2,ny_l-1
            isign_hole = 0
            do ii = i-1,i+1
            do jj = j-1,j+1
                if(ii .ne. i .or. jj .ne. j)then
                    diff = cvr_max(i,j) - cvr_max(ii,jj)
                    isign_hole = isign_hole + int(sign(1.0,diff+.10))

                endif ! not at the center point

            enddo ! jj
            enddo ! ii

            if(isign_hole .eq. -8)then
                write(6,*)' hole detected ',i,j,cvr_max(i,j)

                do jj = j+1,j-1,-1
                    write(6,511,err=512)
     1                    nint(cvr_max(i-1,jj)*100)
     1                   ,nint(cvr_max(i  ,jj)*100)
     1                   ,nint(cvr_max(i+1,jj)*100)
     1                   ,nint2(tb8_k(i-1,jj),1),nint2(tb8_k(i,jj),1)      
     1                   ,nint2(tb8_k(i+1,jj),1)
     1                   ,nint(t_gnd_k(i-1,jj)),nint(t_gnd_k(i,jj))        
     1                   ,nint(t_gnd_k(i+1,jj))
     1                   ,nint2(tb8_k(i-1,jj)-t_gnd_k(i-1,jj),1)
     1                   ,nint2(tb8_k(i  ,jj)-t_gnd_k(i  ,jj),1)
     1                   ,nint2(tb8_k(i+1,jj)-t_gnd_k(i+1,jj),1)
     1                   ,nint(topo(i-1,jj)),nint(topo(i,jj))        
     1                   ,nint(topo(i+1,jj))
     1                   ,nint(cvr_sao_max(i-1,jj)*100)
     1                   ,nint(cvr_sao_max(i  ,jj)*100)
     1                   ,nint(cvr_sao_max(i+1,jj)*100)
     1                   ,nint2(cloud_frac_vis_a(i-1,jj),100)
     1                   ,nint2(cloud_frac_vis_a(i  ,jj),100)
     1                   ,nint2(cloud_frac_vis_a(i+1,jj),100)
 511                format(1x,3i3,4x,3i4,4x,3i4,4x,3i4,4x,3i5,4x,3i3
     1                    ,4x,3i3)
 512            enddo ! jj

            endif ! cloud hole detected

        enddo ! j
        enddo ! i

        if(istat_radar_3dref .eq. 1)then
            call get_max_reflect(radar_ref_3d,nx_l,ny_l,nz_l
     1                          ,r_missing_data,dbz_max_2d)

            call compare_cloud_radar(radar_ref_3d,dbz_max_2d,cvr_max
     1         ,ref_base,cloud_frac_vis_a
     1         ,vis_radar_thresh_cvr,vis_radar_thresh_dbz,r_missing_data       
     1         ,nx_l,ny_l,nz_l)
        endif

        write(6,601)
601     format('  max cloud cover              sao/pirep           ',
     1            20x,'              final analysis')
        scale = 1.
        call array_plot(cvr_sao_max,cvr_max,nx_l,ny_l,'horz cv'
     1                 ,c1_name_array(:,:,1),kcloud,cld_hts,scale)

        write(6,701)
701     format('  max cloud cover           visible satellite     ',
     1            20x,'              final analysis')
        scale = 1.
        call array_plot(cloud_frac_vis_a,cvr_max,nx_l,ny_l,'horz cv'
     1                 ,c1_name_array(:,:,1),kcloud,cld_hts,scale)

        i4_elapsed = ishow_timer()

        do k = 1,nz_l
            write(6,101)k,heights_3d(idb,jdb,k)
101         format(1x,i3,f10.2)
        enddo ! k

        i4_elapsed = ishow_timer()

!       write out cloud field
        do k=1,kcloud
            rlevel = height_to_zcoord2(cld_hts(k),heights_3d
     1                   ,nx_l,ny_l,nz_l,idb,jdb,istatus)   

            if(rlevel .le. float(nz_l))then
                cld_pres_1d(k) = pressure_of_rlevel(rlevel)
            else
                cld_pres_1d(k) = 0. ! cloud height level above pressure grid
            endif

            write(6,102)k,cld_hts(k),cld_pres_1d(k),istatus
102         format(1x,i3,2f10.1,i2)

        enddo ! k

        ext = 'lc3'

        write(6,*)' calling put_clouds_3d'
        call put_clouds_3d(i4time,ext,clouds_3d,cld_hts,cld_pres_1d
     1                                       ,nx_l,ny_l,kcloud,istatus)

        j_status(n_lc3) = ss_normal
        i4_elapsed = ishow_timer()

!       get cloud bases and tops
2000    write(6,*)' calculating cloud base and top'
        do j = 1,ny_l
        do i = 1,nx_l
            cloud_base(i,j) = default_base
            cloud_ceiling(i,j) = default_ceiling
            cloud_top(i,j)  = default_top

            if(cvr_max(i,j) .ge. thresh_cvr_base)then

!             test for cloud base (msl)
              do k = kcloud-1,1,-1
                if(clouds_3d(i,j,k  ) .lt. thresh_cvr_base .and.
     1             clouds_3d(i,j,k+1) .ge. thresh_cvr_base  )then
                    cloud_base(i,j) = 0.5 * (cld_hts(k) + cld_hts(k+1))
                    cloud_base(i,j) = max(cloud_base(i,j),topo(i,j))
                endif
              enddo ! k

            endif ! clouds exist in this column


            if(cvr_max(i,j) .ge. thresh_cvr_top)then

!             test for cloud top (msl)
              do k = 1,kcloud-1
                if(clouds_3d(i,j,k  ) .gt. thresh_cvr_top .and.
     1             clouds_3d(i,j,k+1) .le. thresh_cvr_top  )then
                    cloud_top(i,j) = 0.5 * (cld_hts(k) + cld_hts(k+1))
                endif
              enddo ! k

            endif ! clouds exist in this column

            if(cvr_max(i,j) .ge. thresh_cvr_ceiling)then

!             test for cloud ceiling (agl)
              do k = kcloud-1,1,-1
                if(clouds_3d(i,j,k  ) .lt. thresh_cvr_base .and.
     1             clouds_3d(i,j,k+1) .ge. thresh_cvr_base  )then
                    cloud_ceiling(i,j) = 0.5 * 
     1                    (cld_hts(k) + cld_hts(k+1))
     1                                                - topo(i,j)
                    cloud_ceiling(i,j) = max(cloud_ceiling(i,j),0.)
                endif
              enddo ! k

            endif ! clouds exist in this column

        enddo ! i
        enddo ! j

        i4_elapsed = ishow_timer()

!       calculate cloud analysis implied snow cover
        call cloud_snow_cvr(cvr_max,cloud_frac_vis_a,cldtop_tb8_m
     1          ,tb8_k,topo,di_dh_vis,dj_dh_vis,i_fill_seams          ! i
     1          ,nx_l,ny_l,idb,jdb                                    ! i
     1          ,rlaps_land_frac,tgd_sfc_k,vis_snow_max               ! i
     1          ,grid_spacing_cen_m,r_missing_data,cvr_snow_cycle)    ! i/o

        i4_elapsed = ishow_timer()

!       more ascii plots
        write(6,801)
801     format('                            visible satellite     ',
     1            20x,'      csc  (cycle)  snow cover')
        scale = 1.
        call array_plot(cloud_frac_vis_a,cvr_snow_cycle,nx_l,ny_l
     1                 ,'horz cv',c1_name_array(:,:,1),kcloud,cld_hts
     1                 ,scale)       

        write(6,901)
901     format('                     lm2 (overall) snow cover      ',
     1            20x,'      csc  (cycle)  snow cover')
        scale = 1.
        call array_plot(cvr_snow,cvr_snow_cycle,nx_l,ny_l,'horz cv'
     1                  ,c1_name_array(:,:,1),kcloud,cld_hts,scale)

        do i = 1,nx_l
        do j = 1,ny_l
            if(cldtop_m(i,j)  .ne. r_missing_data .and.
     1         cvr_max(i,j)   .ge. 0.1                            )then       
                plot_mask(i,j) = cloud_top(i,j)            ! set cloud top mask
            else
                plot_mask(i,j) = r_missing_data
            endif
        enddo ! j
        enddo ! i

        write(6,*)' cloud top (band 8 vs. final analysis)'
        scale = .0001
        write(6,1001)
1001    format('  cloud top (km)             band 8                ',
     1            20x,'              final analysis')
        call array_plot(cldtop_tb8_m,plot_mask,nx_l,ny_l,'horz cv'
     1                  ,c1_name_array(:,:,1),kcloud,cld_hts,scale) ! plot band 8 mask


        write(6,1101)
1101    format('  max cloud cover   3.9u     (nocld:cld:added=3:7:9)',      
     1       18x,'              final analysis')

!       set 3.9 micron plot mask
        do i = 1,nx_l
        do j = 1,ny_l
            if(istat_39_a(i,j) .eq. -1)then
                plot_mask(i,j) = 0.3 
            elseif(istat_39_a(i,j) .eq. 0)then
                plot_mask(i,j) = 0.0 
            elseif(istat_39_a(i,j) .eq. 1)then
                if(istat_39_add_a(i,j) .eq. 1)then
                    plot_mask(i,j) = 0.9
                else
                    plot_mask(i,j) = 0.7
                endif
            endif
        enddo ! j
        enddo ! i

        scale = 1.
        call array_plot(plot_mask,cvr_max,nx_l,ny_l,'horz cv'
     1                 ,c1_name_array(:,:,1),kcloud,cld_hts,scale) ! 3.9u mask



!       write out lcb file (cloud base, top, and ceiling fields)
        if(iwrite_output .ge. 1)then
            ext = 'lcb'
            call get_directory(ext,directory,len_dir)

!           call move(cloud_base    ,out_array_3d(1,1,1),nx_l,ny_l)
            call move(cldht_prlx_top,out_array_3d(1,1,1),nx_l,ny_l)

            call move(cloud_top     ,out_array_3d(1,1,2),nx_l,ny_l)

!           call move(cloud_ceiling ,out_array_3d(1,1,3),nx_l,ny_l)
            call move(cldht_prlx_unsm,out_array_3d(1,1,3),nx_l,ny_l)

            call put_clouds_2d(i4time,directory,ext,nx_l,ny_l
     1                                  ,out_array_3d,istatus)
            if(istatus .eq. 1)j_status(n_lcb) = ss_normal
        endif ! iwrite_output

!       this is where we will eventually split the routines, additional data
!       is necessary for more derived fields

        if(n_radar_3dref .gt. 0)then ! write out data (lps - radar_ref_3d)
            write(6,*)' writing out 3d radar reflectivity field'

            var = 'ref'
            ext = 'lps'
            units = 'dbz'
            write(comment,490)istat_radar_2dref,istat_radar_3dref
     1                       ,istat_radar_3dref_orig
 490        format('laps radar reflectivity',3i3)
            call put_laps_3d(i4time,ext,var,units,comment,radar_ref_3d       
     1                                                ,nx_l,ny_l,nz_l)
            j_status(n_lps) = ss_normal

!           generate ascii final reflectivity plot
            write(6,1202)
1202        format('                  radar coverage',50x
     1                              ,'final radar max reflectivity')
            scale = 0.01
            call array_plot(plot_maskr,dbz_max_2d,nx_l,ny_l,'horz cv'
     1                     ,c1_name_array(:,:,1),kcloud,cld_hts,scale) ! radar 

        endif ! n_radar_3dref

        call compare_analysis_to_saos(nx_l,ny_l,cvr_sao_max
     1  ,cloud_frac_vis_a,tb8_k,t_gnd_k,t_sfc_k,cvr_max,r_missing_data
     1  ,dbz_max_2d,cld_snd,ista_snd,max_cld_snd,cld_hts,kcloud
     1  ,n_cld_snd,c_stations,lat_s,lon_s,elev_s,maxstns)

!       reread solar data from latest lso file
        call read_surface_sa(i4time,maxstns,                       ! i
     1   n_obs_b,c_stations,reptype,atype,                         ! o
     1   lat_s,lon_s,elev_s,wx_s,t_s,td_s,                         ! o
     1   kloud,store_amt,store_hgt,                                ! o
     1   rad_s,solar_ea,obstime,istatus)                           ! o

        call compare_analysis_to_rad(i4time,nx_l,ny_l,cvr_sao_max  ! i
     1  ,solar_alt,cvr_snow,cloud_albedo,idb,jdb                   ! i
     1  ,cloud_frac_vis_a,tb8_k,t_gnd_k,td_sfc_k,cvr_max,r_missing_data
     1  ,dbz_max_2d,cld_snd,ista_snd,max_cld_snd,cld_hts,kcloud
     1  ,rad_s,n_cld_snd,c_stations,lat_s,lon_s,elev_s
     1  ,maxstns,n_obs_b,swi_2d)

!       write lcv file
        if(iwrite_output .ge. 1)then
            do i = 1,nx_l
            do j = 1,ny_l
                cvr_water_temp(i,j) = r_missing_data
                ioff = min(max(i + nint(offset_ir_i(i,j)),1),nx_l)
                joff = min(max(j + nint(offset_ir_j(i,j)),1),ny_l)
                tb8_k_offset(i,j) = tb8_k(ioff,joff)
                if(i .eq. idb .and. j .eq. jdb)then
                    write(6,*)'tb8_offset',i,j,ioff,joff
                    write(6,'(" ctr4 cld_frac_vis_a",g12.4)')
     1                               cloud_frac_vis_a(i,j)
                    write(6,'(" ctr4 cloud albedo",g12.4)')
     1                               cloud_albedo(i,j)
                    write(6,'(" ctr4 sat albedo",g12.4)')sat_albedo(i,j)
                endif
            enddo ! j
            enddo ! i

            ext = 'lcv'
            var_a(1) = 'lcv'
            var_a(2) = 'csc'
            var_a(3) = 'cwt'
            var_a(4) = 's8a'
            var_a(5) = 's3a'
            var_a(6) = 'alb'
            var_a(7) = 'cla'
            var_a(8) = 'rqc'
            var_a(9) = 'swi'
            units_a(1) = 'undim'
            units_a(2) = 'undim'
            units_a(3) = 'k'
            units_a(4) = 'k'
            units_a(5) = 'k'
            units_a(6) = ' '
            units_a(7) = ' '
            units_a(8) = ' '
            units_a(9) = 'w/m**2'
            comment_a(1) = 'laps cloud cover'
            comment_a(2) = 'laps cloud analysis implied snow cover'
            comment_a(3) = 'satellite reflectance'
            comment_a(4) = comment_tb8
            comment_a(5) = comment_t39
            comment_a(6) = 'cloud-frac-vis-a'      ! alb
            comment_a(7) = 'cloud-albedo    '      ! cla
            comment_a(8) = 'laps radar quality'
            comment_a(9) = 'downward solar radiation'

            call move(cvr_max       ,out_array_3d(1,1,1),nx_l,ny_l)
            call move(cvr_snow_cycle,out_array_3d(1,1,2),nx_l,ny_l)
            call move(sat_refl      ,out_array_3d(1,1,3),nx_l,ny_l)

!           offset for navigation though not for parallax
            call move(tb8_k           ,out_array_3d(1,1,4),nx_l,ny_l)
            call move(t39_k           ,out_array_3d(1,1,5),nx_l,ny_l) ! s3a
            call move(cloud_frac_vis_a,out_array_3d(1,1,6),nx_l,ny_l) ! alb

            call move(cloud_albedo    ,out_array_3d(1,1,7),nx_l,ny_l) ! cla
            call move(rqc_2d          ,out_array_3d(1,1,8),nx_l,ny_l)
            call move(swi_2d          ,out_array_3d(1,1,9),nx_l,ny_l)

            call put_laps_multi_2d(i4time,ext,var_a,units_a,
     1              comment_a,out_array_3d,nx_l,ny_l,9,istatus)

            if(istatus .eq. 1)j_status(n_lcv) = ss_normal
        endif ! iwrite_output

500     continue

        l_get_cloudtype = .false.

        if(l_get_cloudtype)then ! compute just the cloud type
            l_flag_cloud_type = .true.
            l_flag_mvd = .false.
            l_flag_icing_index = .false.
            l_flag_bogus_w = .false.
            l_deep_vv = .false.
!           call get_cloud_deriv(nx_l,ny_l,nz_l,clouds_3d,cld_hts         ! i
!!   1                          temp_3d,rh_3d_dum,heights_3d,pres_3d,     ! i
!!   1                          istat_radar,radar_3d,grid_spacing_cen_m,  ! i  
!    1                          l_mask_pcptype,                           ! o
!    1                          ibase_array,itop_array,                   ! o
!!   1                          iflag_slwc,slwc_3d,cice_3d,
!!   1                          thresh_cvr_cty_vv,thresh_cvr_lwc,
!!   1                          l_flag_cloud_type,cldpcp_type_3d,         ! i/o
!!   1                          l_flag_mvd,mvd_3d,
!!   1                          l_flag_icing_index,icing_index_3d,
!!   1                          vv_to_height_ratio_cu,                    ! i
!!   1                          vv_to_height_ratio_sc,                    ! i
!!   1                          vv_for_st,                                ! i
!!   1                          l_flag_bogus_w,omega_3d,l_bogus_radar_w,
!    1                          l_deep_vv,                                ! i
!    1                          twet_snow_dum,                            ! i
!!   1                          l_flag_pcp_type,                          ! i
!!   1                          istatus)                                  ! o

!           write cty field          

        endif

!       if needed 't_modelfg' and 'sh_modelfg' are available.
!       'rh_modelfg' is sometimes calculated in 'get_model_fg'

        if(i_varadj .eq. 1)then
            write(6,*)' call cloud_var at end of analysis'
            call cloud_var(i4time,lat,lon
     1                    ,nx_l,ny_l,nz_l,kcloud,heights_3d,temp_3d
     1                    ,t_gnd_k,clouds_3d,cld_hts,tb8_k
     1                    ,cloud_frac_vis_a
     1                    ,subpoint_lat_clo_vis,subpoint_lon_clo_vis  ! i 
     1                    ,r_missing_data                             ! i
     1                    ,di_dh_vis,dj_dh_vis)                       ! i
        endif

999     continue

        write(6,*)' notifications'
        do i = iprod_start,iprod_end
            write(6,*)' ',exts(i),' ',j_status(i),' ',i
        enddo ! i

        write(6,*)' end of cloud analysis package'

        deallocate(cldcv1)
        deallocate(cf_modelfg)
        deallocate(t_modelfg)
        deallocate(sh_modelfg)
        deallocate(ref_modelfg)
        deallocate(cld_snd)
        deallocate(wt_snd)

        return
        end



        subroutine put_clouds_2d(i4time,directory,ext,imax,jmax
     1                          ,field_2dcloud,istatus)

        integer  nfields
        parameter (nfields = 3)

        character*(*) directory
        character*31 ext

        character*125 comment_2d(nfields)
        character*10 units_2d(nfields)
        character*3 var_2d(nfields)
        integer lvl,lvl_2d(nfields)
        character*4 lvl_coord_2d(nfields)

        real field_2dcloud(imax,jmax,nfields)

        write(6,11)directory,ext(1:5)
11      format(' writing 2d clouds ',a50,1x,a5,1x,a3)

        lvl = 0

        lvl_coord_2d(1) = 'msl'
        lvl_coord_2d(2) = 'msl'
        lvl_coord_2d(3) = 'agl'

        var_2d(1) = 'lcb'
        var_2d(2) = 'lct'
        var_2d(3) = 'cce'

        comment_2d(1) = 'laps cloud base'
        comment_2d(2) = 'laps cloud top'
        comment_2d(3) = 'laps cloud ceiling'

        do k = 1,nfields
            lvl_2d(k) = lvl
            units_2d(k) = 'm'
        enddo

        call write_laps_data(i4time,directory,ext,imax,jmax,
     1  nfields,nfields,var_2d,lvl_2d,lvl_coord_2d,units_2d,
     1                     comment_2d,field_2dcloud,istatus)

        return
        end

        subroutine put_clouds_3d(i4time,ext,clouds_3d,cld_hts
     1                          ,cld_pres_1d,ni,nj,nk,istatus)

!       1997 jul 31 k. dritz  - removed include of lapsparms.for, which was
!                               not actually needed for anything.

        character*150 directory
        character*31 ext

        integer nz_cloud_max
        parameter (nz_cloud_max = 42)

        character*125 comment_3d(nz_cloud_max)
        character*10 units_3d(nz_cloud_max)
        character*3 var_3d(nz_cloud_max),var_2d
        integer lvl_3d(nz_cloud_max)
        character*4 lvl_coord_3d(nz_cloud_max)

        real clouds_3d(ni,nj,nk)
        real cld_hts(nk)
        real cld_pres_1d(nk)

        call get_directory(ext,directory,len_dir)

        var_2d = 'lc3'

        write(6,11)directory,ext(1:5),var_2d
11      format(' writing 3d ',a50,1x,a5,1x,a3)

        do k = 1,nk
            units_3d(k)   = 'fractional'
            lvl_3d(k) = k
            lvl_coord_3d(k) = 'msl'

            var_3d(k) = var_2d

            write(comment_3d(k),1)cld_hts(k),cld_pres_1d(k)
1           format(2e20.8,' height msl, pressure')

        enddo ! k

        call write_laps_data(i4time,directory,ext,ni,nj,
     1  nk,nk,var_3d,lvl_3d,lvl_coord_3d,units_3d,
     1                     comment_3d,clouds_3d,istatus)

        return
        end


        subroutine cloud_snow_cvr(cvr_max,cloud_frac_vis_a,cldtop_tb8_m
     1             ,tb8_k,topo,di_dh,dj_dh,i_fill_seams                ! i
     1             ,ni,nj,idb,jdb                                      ! i
     1             ,rlaps_land_frac,tgd_sfc_k                          ! i
     1             ,vis_snow_max                                       ! i
     1             ,grid_spacing_cen_m,r_missing_data,cvr_snow_cycle)

        real cvr_max(ni,nj)          ! input
        real cloud_frac_vis_a(ni,nj) ! input
        real cldtop_tb8_m(ni,nj)     ! input
        real tb8_k(ni,nj)            ! input
        real tgd_sfc_k(ni,nj)        ! input
        real topo(ni,nj)             ! input
        real rlaps_land_frac(ni,nj)  ! input
        real di_dh(ni,nj)            ! input           
        real dj_dh(ni,nj)            ! input
        integer i_fill_seams(ni,nj)  ! input
        real vis_snow_max(ni,nj)     ! input
        real cvr_snow_cycle(ni,nj)   ! output

        logical l_cvr_max            ! local

        n_csc_pts = 0
        n_no_csc_pts = 0
        n_clear_pts = 0
        n_cld_pts = 0

        n_snow = 0
        n_nosnow = 0
        n_missing = 0

        cvr_snow_cycle = r_missing_data ! initialize

!       buffer is needed, particularly if parallax isn't fully
!       considered for ir data and for cvr_max value
        ip = max(nint(8000./grid_spacing_cen_m),1)

        write(6,*)' snow cover perimeter is ',ip

!       loop in gridpoint space
        do i = 1,ni
        do j = 1,nj

            l_cvr_max = .false.

!           make sure all neighbors are clear in addition to the grid point
            il = max(i-ip,1)
            ih = min(i+ip,ni)
            do ii = il,ih
                jl = max(j-ip,1)
                jh = min(j+ip,nj)
                do jj = jl,jh
                    if(cvr_max(ii,jj) .gt. 0.1)l_cvr_max = .true.
                enddo ! jj
            enddo ! ii

            if(                cvr_max(i,j) .le. 0.1        )then! no cld cover
                n_clear_pts = n_clear_pts + 1
            else
                n_cld_pts = n_cld_pts + 1
            endif

!           add parallax correction to snow cover
            ig = i + nint(di_dh(i,j) * topo(i,j))
            jg = j + nint(dj_dh(i,j) * topo(i,j))
            ig = max(min(ig,ni),1)
            jg = max(min(jg,nj),1)

            it = i - nint(di_dh(i,j) * topo(i,j))
            jt = j - nint(dj_dh(i,j) * topo(i,j))
            it = max(min(it,ni),1)
            jt = max(min(jt,nj),1)
            if(i_fill_seams(i,j) .ne. 0)then
                itn = min(i,it)
                itx = max(i,it)
            else
                itn = it
                itx = it
            endif
!           itn = i
!           itx = i
!           jt = j

!           consider visible satellite
            if(.not. l_cvr_max)then                             ! no cld cover
                cvr_snow_cycle(i,j) = vis_snow_max(ig,jg)
            endif

            if(vis_snow_max(ig,jg) .lt. 0.2)then
                cvr_snow_cycle(i,j) = 0.
            endif

            if(rlaps_land_frac(i,j) .le. 0.25 .and. 
     1         topo(i,j) .le. 10.                   )then
                iocean = 1 
            else
                iocean = 0 
            endif
            tb8_snow_thr = 281.15 - float(iocean) * 6.

!           if ground/ocean is warm then no snow is present
            if(tb8_k(i,j) .ne. r_missing_data
     1       .and.        tb8_k(i,j) .gt. tb8_snow_thr
     1       .and.        cvr_max(i,j) .le. 0.1     )then
                cvr_snow_cycle(i,j) = 0.
            endif

!           if no ir set snow cover based on ocean "ground" temp
            if(tb8_k(i,j) .eq. r_missing_data .and.         
     1         iocean .eq. 1                  .and. 
     1         tgd_sfc_k(i,j) .ne. r_missing_data   )then
              if(tgd_sfc_k(i,j) .gt. 274.15)then      
                cvr_snow_cycle(i,j) = 0.
              elseif(tgd_sfc_k(i,j) .lt. 272.15)then
                cvr_snow_cycle(i,j) = 1.
              else
                cvr_snow_cycle(i,j) = (274.15-tgd_sfc_k(i,j))/2.
              endif
            endif

            if(i .eq. idb .and. j .eq. jdb)then
                write(6,11)idb,jdb,l_cvr_max,vis_snow_max(i,j)
     1                    ,tb8_k(i,j),cvr_snow_cycle(i,j)
11              format(2i5,' vssnmx/tb8/cvsncyc',l2,f8.3,f8.1,f8.3
     1                     ' ctr (cvr_snw_cyc)')
            endif

        enddo ! j
        enddo ! i

        do i = 1,ni
        do j = 1,nj

!           count various categories of csc
            if(cvr_snow_cycle(i,j) .ne. r_missing_data)then
                if(cvr_snow_cycle(i,j) .gt. 0.1)then
                    n_snow = n_snow + 1
                else
                    n_nosnow = n_nosnow + 1
                endif
            else
                n_missing = n_missing + 1
            endif

        enddo ! j
        enddo ! i

        write(6,*)' # csc/nocsc/clr/cld pts = ',n_csc_pts,n_no_csc_pts
     1                                       ,n_clear_pts,n_cld_pts
        write(6,1)n_snow,n_nosnow,n_missing
 1      format('  # snow/nosnow/missing ',3i7)

        return
        end


        subroutine compare_analysis_to_saos(ni,nj,cvr_sao_max
     1  ,cloud_frac_vis_a,tb8_k,t_gnd_k,t_sfc_k,cvr_max,r_missing_data
     1  ,dbz_max_2d,cld_snd,ista_snd,max_cld_snd,cld_hts,kcloud
     1  ,n_cld_snd,c_stations,lat_s,lon_s,elev_s,maxstns)

        real cloud_frac_vis_a(ni,nj),tb8_k(ni,nj),t_gnd_k(ni,nj)
     1        ,t_sfc_k(ni,nj),cvr_max(ni,nj),cvr_sao_max(ni,nj)
     1        ,dbz_max_2d(ni,nj)

        real cld_snd(max_cld_snd,kcloud)
        integer ista_snd(max_cld_snd)
        real cld_hts(kcloud)

        character c_stations(maxstns)*(*)
        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)

        character*3 c3_discrep
        character*1 c1_c

        do j = 1,nj
        do i = 1,ni
            cvr_max(i,j) = min(cvr_max(i,j),1.00)
        enddo ! i
        enddo ! j

        if(.true.)then
            write(6,*)
     1      ' comparing cloud/sat/sfc data at sao/pirep locations'
            iwrite = 0

            n_ovc = 0
            n_tovc = 0
            n_sct = 0
            n_bkn = 0
            n_tsct = 0
            n_tbkn = 0
            ovc_sum = 0.
            tovc_sum = 0.
            sct_sum = 0.
            bkn_sum = 0.
            tsct_sum = 0.
            tbkn_sum = 0.


            do isnd = 1,n_cld_snd
              ista = ista_snd(isnd)
              if(ista .ne. 0)then
                call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista),lat,lon
     1                          ,ni,nj,ri,rj,istatus)

                i_i = nint(ri)
                i_j = nint(rj)

                if(i_i .ge. 3 .and. i_i .le. ni-2 .and.
     1             i_j .ge. 3 .and. i_j .le. nj-2            )then

                    if(iwrite .eq. iwrite/20*20)then
                        write(6,*)
                        write(6,*)'sta   i    j   vis frac tb8_k  '
     1                  //'t_gnd_k t_sfc_k  cldsnd cv_sa_mx cvr_mx '
     1                  //'snd-ht  dbz        9pt    25pt'
                    endif

!                   calculate 9pt cover
                    cvr_9pt = 0.
                    do ii = -1,1
                    do jj = -1,1
                        cvr_9pt = cvr_9pt + cvr_max(i_i+ii,i_j+jj)
                    enddo ! jj
                    enddo ! ii
                    cvr_9pt = cvr_9pt / 9.

!                   calculate 25pt cover
                    cvr_25pt = 0.
                    do ii = -2,2
                    do jj = -2,2
                        cvr_25pt = cvr_25pt + cvr_max(i_i+ii,i_j+jj)
                    enddo ! jj
                    enddo ! ii
                    cvr_25pt = cvr_25pt / 25.

                    iwrite = iwrite + 1

                    cld_snd_max = -.999
                    height_max = 0.
                    do k = 1,kcloud
                        if(cld_snd(isnd,k) .ne. r_missing_data)then
                            cld_snd_max = max(cld_snd_max
     1                                       ,cld_snd(isnd,k))
                            height_max = cld_hts(k)
                        endif
                    enddo ! k

!                   cld_snd_max = cvr_snd(isnd)

!                   flag significant discrepancies between saos and analysis
                    ht_12000 = elev_s(ista) + 12000./3.281 + 1.
                    c3_discrep = '   '

                    if(cvr_max(i_i,i_j) - cld_snd_max .le. -0.25)then
                        if(cvr_max(i_i,i_j) .le. 0.1)then
                            c3_discrep = ' **'
                        else
                            c3_discrep = ' . '
                        endif
                    endif

                    if(cvr_max(i_i,i_j) - cld_snd_max .ge. +0.25)then
                        if(height_max .gt. ht_12000)then
                            if(cld_snd_max .le. 0.1)then
                               c3_discrep = ' **'
                            else
                               c3_discrep = ' . '
                            endif
                        endif
                    endif

                    if(cvr_max(i_i,i_j) - cld_snd_max .ge. +0.50)then
                        if(height_max .gt. ht_12000)c3_discrep = ' **'
                    endif

                    if(cvr_max(i_i,i_j) - cld_snd_max .le. -0.50)then
                        c3_discrep = ' **'
                    endif

                    if(cvr_sao_max(i_i,i_j) .lt. cld_snd_max - 0.1)then       
                        c3_discrep = ' ss'
                    endif

                    icat1 = 1
                    if(cvr_25pt .ge. 0.10)icat1 = 2
                    if(cvr_25pt .ge. 0.50)icat1 = 3
                    if(cvr_25pt .ge. 0.90)icat1 = 4

                    icat2 = 1
                    if(cld_snd_max .ge. 0.10)icat2 = 2
                    if(cld_snd_max .ge. 0.50)icat2 = 3
                    if(cld_snd_max .ge. 0.90)icat2 = 4

                    if(icat1 .ne. icat2
     1            .and. abs(cld_snd_max - cvr_25pt) .gt. .10
     1                              .and.
     1           (height_max .gt. ht_12000 .or. cld_snd_max .eq. 1.00)
     1                                                             )then       
                        c1_c = 'c'
                    else
                        c1_c = ' '
                    endif

                    if(abs(cld_snd_max - 1.00) .lt. .01)then
                        n_ovc = n_ovc + 1
                        ovc_sum = ovc_sum + cvr_25pt
                    endif

                    if(height_max .gt. ht_12000)then
                        if(abs(cld_snd_max - .60) .lt. .01)then
                            n_tovc = n_tovc + 1
                            tovc_sum = tovc_sum + cvr_25pt
                        endif
                        if(abs(cld_snd_max - .25) .lt. .01)then
                            n_sct = n_sct + 1
                            sct_sum = sct_sum + cvr_25pt
                        endif
                        if(abs(cld_snd_max - .70) .lt. .01)then
                            n_bkn = n_bkn + 1
                            bkn_sum = bkn_sum + cvr_25pt
                        endif
                        if(abs(cld_snd_max - .15) .lt. .01)then
                            n_tsct = n_tsct + 1
                            tsct_sum = tsct_sum + cvr_25pt
                        endif
                        if(abs(cld_snd_max - .40) .lt. .01)then
                            n_tbkn = n_tbkn + 1
                            tbkn_sum = tbkn_sum + cvr_25pt
                        endif
                    endif ! > 12000 agl type station

                    write(6,1111,err=1112)c_stations(ista)(1:3)
     1                           ,i_i,i_j
     1                           ,cloud_frac_vis_a(i_i,i_j)
     1                           ,tb8_k(i_i,i_j)
     1                           ,t_gnd_k(i_i,i_j)
     1                           ,t_sfc_k(i_i,i_j)
     1                           ,cld_snd_max
     1                           ,cvr_sao_max(i_i,i_j)
     1                           ,cvr_max(i_i,i_j)
     1                           ,height_max
     1                           ,dbz_max_2d(i_i,i_j)
     1                           ,c3_discrep
     1                           ,cvr_9pt
     1                           ,cvr_25pt
     1                           ,c1_c
1111                format(1x,a3,2i5,f8.2,3f8.1,3f8.2,f8.0,f5.0       
     1                    ,a3,f8.2,f7.2,1x,a1)

1112            endif ! ob is in domain
              endif ! ista .ne. 0 (valid value)
            enddo ! isnd

            if(n_sct .gt. 0)write(6,11)n_sct,sct_sum/float(n_sct)
11          format(' mean analysis value for  sct (.25) sao  = ',i3,f7.2
     1)

            if(n_bkn .gt. 0)write(6,12)n_bkn,bkn_sum/float(n_bkn)
12          format(' mean analysis value for  bkn (.70) sao  = ',i3,f7.2
     1)

            if(n_tsct .gt. 0)write(6,13)n_tsct,tsct_sum/float(n_tsct)
13          format(' mean analysis value for -sct (.20) sao  = ',i3,f7.2
     1)

            if(n_tbkn .gt. 0)write(6,14)n_tbkn,tbkn_sum/float(n_tbkn)
14          format(' mean analysis value for -bkn (.40) sao  = ',i3,f7.2
     1)

            if(n_ovc .gt. 0)write(6,15)n_ovc,ovc_sum/float(n_ovc)
15          format(' mean analysis value for  ovc (1.00) sao = ',i3,f7.2
     1)

            if(n_tovc .gt. 0)write(6,16)n_tovc,tovc_sum/float(n_tovc)
16          format(' mean analysis value for -ovc  (.60) sao = ',i3,f7.2
     1)

            write(6,*)

        endif ! do sao comparison to analysis

        return
        end


        subroutine compare_cloud_radar(radar_ref_3d,dbz_max_2d,cvr_max
     1          ,ref_base,cloud_frac_vis_a
     1          ,vis_radar_thresh_cvr,vis_radar_thresh_dbz
     1          ,r_missing_data,ni,nj,nk)

        real radar_ref_3d(ni,nj,nk)                   ! i
        real dbz_max_2d(ni,nj)                        ! i
        real cvr_max(ni,nj)                           ! i
        real cloud_frac_vis_a(ni,nj)                  ! i

!       this routine compares the cloud and radar fields and flags
!       remaining differences that weren't caught in earlier processing

        write(6,*)'comparing clouds and radar (_rdr)'
        write(6,*)'vis_radar_thresh_cvr = ',vis_radar_thresh_cvr
        write(6,*)'vis_radar_thresh_dbz = ',vis_radar_thresh_dbz

        do i = 1,ni
        do j = 1,nj
            if(       cvr_max(i,j)    .lt. 0.2
     1          .and. dbz_max_2d(i,j) .gt. ref_base
     1          .and. dbz_max_2d(i,j) .ne. r_missing_data
     1                                                            )then

!             we have a discrepancy between the vis and radar

              iblank_radar = 0

              if(cvr_max(i,j) .eq. cloud_frac_vis_a(i,j))then ! cvr = vis

                if(dbz_max_2d(i,j) .lt. vis_radar_thresh_dbz)then
                  write(6,1)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                     ,cloud_frac_vis_a(i,j)
1                 format(' vis_rdr: cvr/dbz/vis <',2i4,f8.2,f8.1,f8.2)

                else
                  write(6,11)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                      ,cloud_frac_vis_a(i,j)
11                format(' vis_rdr: cvr/dbz/vis >',2i4,f8.2,f8.1,f8.2)

                endif


              elseif(cvr_max(i,j) .lt. cloud_frac_vis_a(i,j))then

!                 don't blame the vis            cvr < vis

                  if(dbz_max_2d(i,j) .lt. vis_radar_thresh_dbz)then

!                     we will block out the radar
                      write(6,2)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
2                     format(' cld_rdr: cvr/dbz/vis <',2i4,f8.2
     1                      ,f8.1,f8.2,' blank out radar')

                      iblank_radar = 1

                  else ! radar is too strong to block out

                      write(6,3)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
3                     format(' cld_rdr: cvr/dbz/vis >',2i4,f8.2
     1                                                 ,f8.1,f8.2)

                  endif ! ref < vis

              elseif(cvr_max(i,j) .gt. cloud_frac_vis_a(i,j))then

!                 don't know if vis lowered cloud cover        cvr > vis

                  if(dbz_max_2d(i,j) .lt. vis_radar_thresh_dbz)then

!                     we will block out the radar
                      write(6,4)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
4                     format(' ???_rdr: cvr/dbz/vis <',2i4,f8.2
     1                      ,f8.1,f8.2,' blank out radar')

                      iblank_radar = 1

                  else ! radar is too strong to block out

                      write(6,5)i,j,cvr_max(i,j),dbz_max_2d(i,j)
     1                       ,cloud_frac_vis_a(i,j)
5                     format(' ???_rdr: cvr/dbz/vis >',2i4,f8.2
     1                                                 ,f8.1,f8.2)

                  endif ! ref < vis

              endif ! cover compare to vis

              if(iblank_radar .eq. 1)then        ! take action
                  dbz_max_2d(i,j) = ref_base
                  do kl = 1,nk
                      radar_ref_3d(i,j,kl) = ref_base
                  enddo ! kl
              endif

            endif ! radar echo with low cloud cover
        enddo ! j
        enddo ! i

        return
        end

        function nint2(x,ifactor)

!       this routine helps scale the arguments for ascii debug printing

        call get_r_missing_data(r_missing_data,istatus)

        if(x .ne. r_missing_data)then
            nint2 = nint(x*float(ifactor))
        else
            nint2 = 999999
        endif

        return
        end

        subroutine qc_clouds_3d(clouds_3d,nx_l,ny_l,kcloud)

        real clouds_3d(nx_l,ny_l,kcloud)
        logical l_poss_extrap ! used to allow for edge effects from 'barnes_r5'

        nskip_max = 4 ! 'see barnes_r5'

        i4_elapsed = ishow_timer()

        write(6,*)' subroutine qc_clouds_3d...'

        do i = 1,nx_l
        do j = 1,ny_l
            if(nx_l-i .le. nskip_max .or. ny_l-j .le. nskip_max)then
                l_poss_extrap = .true. ! extrapolation edge effects possible
            else
                l_poss_extrap = .false.
            endif

            do k = 1,kcloud
!               call qc_clouds_0d(i,j,k,clouds_3d(i,j,k)
!    1                           ,nx_l,ny_l,l_poss_extrap)

!               subroutine code reproduced in calling routine for efficiency
                clouds_0d = clouds_3d(i,j,k)

                if(clouds_0d .gt. 1.0)then 
                    if(.not. l_poss_extrap)then
                        if(clouds_0d .gt. 1.001)then
                            write(6,*)
     1                          ' error, clouds_0d > 1',i,j,k,clouds_0d       
                            stop
                        else ! just over 1.0 with no edge effect
                            write(6,*)
     1                      ' warning, clouds_0d > 1 - reset'
     1                      ,i,j,k,clouds_0d
                            clouds_0d = 1.0
                        endif
                    else
                        write(6,*)
     1                 ' warning, clouds_0d > 1 - reset for edge effect'       
     1                  ,i,j,k,clouds_0d
                        clouds_0d = 1.0
                    endif

                elseif(clouds_0d .lt. 0.0)then 
                    if(l_poss_extrap)then
                        qc_thr = -200.0
                    else
                        qc_thr = -0.0005
                    endif
                    if(clouds_0d .lt. qc_thr)then
                        write(6,*)' error, clouds_0d << 0',i,j,k
     1                                    ,clouds_0d   
                        stop
                    else 
                        write(6,*)
     1                 ' warning, clouds_0d < 0 - reset for edge effect'       
     1                  ,i,j,k,clouds_0d
                        clouds_0d = 0.
                    endif
                endif

                clouds_3d(i,j,k) = clouds_0d

            enddo ! k
        enddo ! j
        enddo ! i

        i4_elapsed = ishow_timer()

        return
        end

        subroutine qc_clouds_0d(i,j,k,clouds_3d
     1                         ,nx_l,ny_l,l_poss_extrap)

        real clouds_3d
        logical l_poss_extrap ! used to allow for edge effects from 'barnes_r5'

        if(clouds_3d .gt. 1.0)then 
            if(.not. l_poss_extrap)then
                if(clouds_3d .gt. 1.001)then
                    write(6,*)' error, clouds_3d > 1',i,j,k,clouds_3d
                    stop
                else ! just over 1.0 with no edge effect
                    write(6,*)
     1              ' warning, clouds_3d > 1 - reset'
     1              ,i,j,k,clouds_3d
                    clouds_3d = 1.0
                endif
            else
                write(6,*)
     1          ' warning, clouds_3d > 1 - reset for edge effect'
     1          ,i,j,k,clouds_3d
                clouds_3d = 1.0
            endif
        endif

        if(clouds_3d .lt. 0.0)then 
            write(6,*)' error: clouds_3d < 0 - reset',i,j,k,clouds_3d       
            clouds_3d = 0.0
        endif

        return
        end

        subroutine get_parallax_info(ni,nj,i4time                      ! i
     1                              ,lat,lon                           ! i
     1                              ,subpoint_lat_clo,subpoint_lon_clo ! i
     1                              ,di_dh,dj_dh,i_fill_seams)         ! o

        include 'trigd.inc'

        real subpoint_lat_clo(ni,nj)           ! i
        real subpoint_lon_clo(ni,nj)           ! i
        real lat(ni,nj)                        ! i
        real lon(ni,nj)                        ! i
        real alt(ni,nj)                        ! l (emission angle)
        real azi(ni,nj)                        ! l
        real phase(ni,nj)                      ! l
        real spec(ni,nj)                       ! l
        real dx(ni,nj)                         ! l
        real dy(ni,nj)                         ! l
        real projrot_laps(ni,nj)               ! l

        integer i_fill_seams(ni,nj)            ! o
        real di_dh(ni,nj)                      ! o
        real dj_dh(ni,nj)                      ! o

        i4_elapsed = ishow_timer()

!       satellite geometry for parallax offset
        write(6,*)' calling satgeom...'
        range_m = 42155680.00
        call satgeom(i4time,lat,lon,ni,nj 
     1      ,subpoint_lat_clo,subpoint_lon_clo
     1      ,range_m,r_missing_data,phase,spec,alt,azi,istatus)

        call get_grid_spacing_array(lat,lon,ni,nj,dx,dy)
        call projrot_latlon_2d(lat,lon,ni,nj,projrot_laps,istatus)

!       calculate parallax offset (sat/lvd grid index minus analysis grid
!       index)
        do i = 1,ni
        do j = 1,nj
            if(alt(i,j) .gt. 0.)then
                ds_dh = tand(90. - alt(i,j))
                azi_grid = azi(i,j) - projrot_laps(i,j)
                di_dh(i,j) = (ds_dh / dx(i,j)) * (-sind(azi_grid))
                dj_dh(i,j) = (ds_dh / dy(i,j)) * (-cosd(azi_grid))
            else
                di_dh(i,j) = 0.
                dj_dh(i,j) = 0.
            endif
        enddo ! i
        enddo ! j

        do j = 1,nj
        do i = 1,ni
            i_fill_seams(i,j) = 0
            im1 = max(i-1,1)
            ip1 = min(i+1,ni)
            if(subpoint_lon_clo(ip1,j) .gt. subpoint_lon_clo(i,j))then       
                i_fill_seams(i,j) = -1 ! w side of seam
            endif
            if(subpoint_lon_clo(im1,j) .lt. subpoint_lon_clo(i,j))then
                i_fill_seams(i,j) = +1 ! e side of seam
            endif
        enddo ! i
        enddo ! j

        i4_elapsed = ishow_timer()

        return
        end
