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
        program laps_deriv_main

        integer j_status(20)

        character*9 a9_time

        call get_systime(i4time,a9_time,istatus)
        if(istatus .ne. 1)go to 999

        write(6,*)' systime = ',a9_time

        call get_grid_dim_xy(nx_l,ny_l,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting horizontal domain dimensions'
           go to 999
        endif

        call get_laps_dimensions(nz_l,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting vertical domain dimension'
           go to 999
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting r_missing_data'
           go to 999
        endif
          
        call laps_deriv(i4time,
     1                  nx_l,ny_l,
     1                  nz_l,
     1                  r_missing_data,
     1                  j_status)

999     continue

        end
          
        subroutine laps_deriv(i4time,
     1                  nx_l,ny_l,
     1                  nz_l,
     1                  r_missing_data,
     1                  j_status)

        use mem_namelist, only: max_snd_grid, max_snd_levels

        integer j_status(20),iprod_number(20)

        real temp_3d(nx_l,ny_l,nz_l)
        real rh_3d_pct(nx_l,ny_l,nz_l)
        real td_3d_k(nx_l,ny_l,nz_l)
        real heights_3d(nx_l,ny_l,nz_l)
        real u_3d(nx_l,ny_l,nz_l)
        real v_3d(nx_l,ny_l,nz_l)

        real temp_sfc_k(nx_l,ny_l)
        real pres_sfc_pa(nx_l,ny_l)
        real rh_sfc_pct(nx_l,ny_l)
        real tpw_2d(nx_l,ny_l)         ! units are m
        real u_sfc_ms(nx_l,ny_l)
        real v_sfc_ms(nx_l,ny_l)

        real dbz_max_2d(nx_l,ny_l)

        real lat(nx_l,ny_l)
        real lon(nx_l,ny_l)
        real topo(nx_l,ny_l)
        real ldf(nx_l,ny_l)

        character*31 ext

        character*10  units_2d
        character*125 comment_2d
        character*3 var_2d

        logical l_cloud_only

!       get parameters for laps_deriv_sub call
        call get_meso_sao_pirep(n_meso,n_sao,n_pirep,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting n_pirep'
           go to 999
        endif

        call get_maxstns(maxstns,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting maxstns'
           go to 999
        endif

        max_cld_snd = maxstns + n_pirep

!       read data for laps_deriv_sub and put_stability calls

!       read lt1 - temp_3d
        var_2d = 't3'
        ext = 'lt1'
        call get_laps_3d(i4time,nx_l,ny_l,nz_l
     1      ,ext,var_2d,units_2d,comment_2d,temp_3d,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading 3d temp'
            return
        endif
        call qc_field_3d('t3',temp_3d,nx_l,ny_l,nz_l,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading 3d temp'
            return
        endif

!       read lt1 - heights_3d
        var_2d = 'ht'
        ext = 'lt1'
        call get_laps_3d(i4time,nx_l,ny_l,nz_l
     1      ,ext,var_2d,units_2d,comment_2d,heights_3d,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading 3d heights'
            return
        endif
        call qc_field_3d('ht',heights_3d,nx_l,ny_l,nz_l,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading 3d heights'
            return
        endif

!       read rh
        var_2d = 'rhl'
        ext = 'lh3'
        call get_laps_3d(i4time,nx_l,ny_l,nz_l
     1      ,ext,var_2d,units_2d,comment_2d,rh_3d_pct,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading 3d rh (lh3/rhl)'
            return
        endif
        call qc_field_3d('rhl',rh_3d_pct,nx_l,ny_l,nz_l,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading 3d rh (lh3/rhl)'
            return
        endif

!       read in surface temp data
        var_2d = 't'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,nx_l,ny_l,temp_sfc_k,istatus)

        if(istatus .ne. 1)then
            write(6,*)' laps sfc temp not available'
            go to 999
        endif

!       read in surface pressure data
        var_2d = 'ps'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,nx_l,ny_l,pres_sfc_pa,istatus)

        if(istatus .ne. 1)then
            write(6,*)' laps sfc pres not available'
            go to 999
        endif

!       read in surface rh data
        var_2d = 'rh'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,nx_l,ny_l,rh_sfc_pct,istatus)

        if(istatus .ne. 1)then
            write(6,*)' laps sfc rh not available'
            go to 999
        endif

!       read in surface u component
        var_2d = 'u'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,nx_l,ny_l,u_sfc_ms,istatus)

        if(istatus .ne. 1)then
            write(6,*)' laps sfc u not available'
            go to 999
        endif

!       read in surface v component
        var_2d = 'v'
        ext = 'lsx'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,nx_l,ny_l,v_sfc_ms,istatus)

        if(istatus .ne. 1)then
            write(6,*)' laps sfc v not available'
            go to 999
        endif

!       read in tpw data
        var_2d = 'tpw'
        ext = 'lh4'
        call get_laps_2d(i4time,ext,var_2d,units_2d,comment_2d
     1                  ,nx_l,ny_l,tpw_2d,istatus)

        if(istatus .ne. 1)then
            write(6,*)' laps tpw not available'
            go to 999
        endif

        write(6,*)
        write(6,*)' calling laps_deriv_sub'
        call laps_deriv_sub(i4time,                  ! i
     1                  nx_l,ny_l,                   ! i
     1                  nz_l,                        ! i
     1                  n_pirep,                     ! i
     1                  maxstns,                     ! i
     1                  max_snd_grid,max_snd_levels, ! i
     1                  n_prods,
     1                  iprod_number,
     1                  temp_3d,                     ! i
     1                  heights_3d,                  ! i
     1                  rh_3d_pct,                   ! i
     1                  pres_sfc_pa,                 ! i
     1                  temp_sfc_k,                  ! i
     1                  dbz_max_2d,istat_lps,        ! o
     1                  twet_snow,l_cloud_only,      ! o
     1                  j_status,                    ! o
     1                  istatus1)                    ! o

        if(l_cloud_only)then
            write(6,*)' skipping stability and fire sections'
            go to 999
        endif

!       istat_lps = 0 ! temporary until this is wired in'
        istat_twet_snow = 1 ! if we get this far the read was successful

        call get_laps_domain_95(nx_l,ny_l,lat,lon,topo
     1                         ,ldf,grid_spacing_m,istatus)       
        if(istatus .ne. 1)then
            write(6,*)' error getting laps domain'
            go to 999
        endif

        if(istat_twet_snow .eq. 1)then
            call get_laps_cycle_time(laps_cycle_time,istatus)
            if (istatus .ne. 1) then
                write (6,*) 'error getting laps cycle time'
                go to 999
            endif

            i4_elapsed = ishow_timer()

            write(6,*)
            write(6,*)' calling put_stability'
            call put_stability(
     1           i4time                          ! i
     1          ,nx_l,ny_l,nz_l                  ! i
     1          ,heights_3d                      ! i
     1          ,lat,lon,topo                    ! i
     1          ,laps_cycle_time                 ! i
     1          ,temp_3d                         ! i
     1          ,rh_3d_pct                       ! i
     1          ,temp_sfc_k                      ! i
     1          ,pres_sfc_pa                     ! i
     1          ,twet_snow                       ! i
     1          ,td_3d_k                         ! o
     1          ,istat_lst)                      ! o
        else
            write(6,*)' put_stability not called for lst file'

        endif

        i4_elapsed = ishow_timer()

!       if we need space we can deallocate rh_3d_pct here
!       if we need space we can allocate u_3d, v_3d here

        write(6,*)
        write(6,*)' calling put_derived_wind_prods'
        call put_derived_wind_prods(nx_l,ny_l,nz_l           ! input
     1          ,nx_l,ny_l,nz_l                              ! input (sic)
     1          ,max_radars_dum,r_missing_data               ! input
     1          ,i4time                                      ! input
     1          ,dbz_max_2d,istat_lps                        ! input
     1          ,lat,lon,topo,ldf                            ! input
     1          ,tpw_2d                                      ! input
     1          ,heights_3d                                  ! input
     1          ,u_3d,v_3d)                                  ! output

        write(6,*)
        if(istat_lst .eq. 1)then

            write(6,*)' calling fire_fields'
            call fire_fields(nx_l,ny_l,nz_l,temp_3d,td_3d_k          ! i
     1                      ,u_3d,v_3d                               ! i
     1                      ,temp_sfc_k,pres_sfc_pa                  ! i
     1                      ,rh_sfc_pct,u_sfc_ms,v_sfc_ms            ! i
     1                      ,r_missing_data,i4time                   ! i
     1                      ,istatus)                                ! o

        else
            write(6,*)' skipping call to fire_fields'

        endif

 999    write(6,*)' end of subroutine laps_deriv'

        return

        end

 
       subroutine get_deriv_parms(mode_evap,l_bogus_radar_w,              ! o
     1                            l_deep_vv,l_cloud_only,                 ! o
     1                            vv_to_height_ratio_cu,                  ! o
     1                            vv_to_height_ratio_sc,                  ! o
     1                            vv_for_st,                              ! o
     1                            c_z2m,                                  ! o
     1                            thresh_cvr_cty_vv,thresh_cvr_lwc,       ! o
     1                            twet_snow,                              ! o
     1                            hydrometeor_scale_cldliq,               ! o
     1                            hydrometeor_scale_cldice,               ! o
     1                            hydrometeor_scale_pcp,                  ! o
     1                            istatus)                                ! o

       real vv_to_height_ratio_cu
       real vv_to_height_ratio_sc
       real vv_for_st
       real thresh_cvr_cty_vv,thresh_cvr_lwc,twet_snow
       real hydrometeor_scale_cldliq,hydrometeor_scale_cldice
       real hydrometeor_scale_pcp

       character*20 c_z2m

       logical l_bogus_radar_w, l_deep_vv, l_cloud_only

       namelist /deriv_nl/ mode_evap, l_bogus_radar_w, l_deep_vv,
     1                     l_cloud_only,
     1                     vv_to_height_ratio_cu,
     1                     vv_to_height_ratio_sc,
     1                     vv_for_st,
     1                     c_z2m,
     1                     thresh_cvr_cty_vv,thresh_cvr_lwc,
     1                     twet_snow,
     1                     hydrometeor_scale_cldliq,
     1                     hydrometeor_scale_cldice,
     1                     hydrometeor_scale_pcp
 
       character*150 static_dir,filename

!      set default value(s)
       c_z2m = 'albers'
       hydrometeor_scale_cldliq = 0.2
       hydrometeor_scale_cldice = 0.2
       hydrometeor_scale_pcp    = 1.0
       l_cloud_only = .false.
 
       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/deriv.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,deriv_nl,err=901)
       close(1)

       print*,'success reading deriv_nl in ',filename
       write(*,deriv_nl)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading deriv_nl in ',filename
       write(*,deriv_nl)
       istatus = 0
       return

       end
