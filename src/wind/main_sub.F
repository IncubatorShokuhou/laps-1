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


        subroutine lapswind_anal  (i4time_sys,
     1              ntmin,ntmax,                                     ! i
     1              n_meso,n_sao,                                    ! i
     1              heights_3d,                                      ! i
     1              istat_wind)                                      ! o

!       1990         steve albers  (original version)
!       1992 apr 30  steve albers  (mod for multi-radar doppler velocities)
!       1993 feb 20  steve albers  subroutinize the profiler processing
!       1993 apr     steve albers  accept rams background
!       1993 nov     steve albers  partial oe interface, 1st guess subroutine
!                                  eliminate arrays.f from other subroutines
!                                  eliminate common blocks
!       1994 oct     steve albers  add i4time_radar_a array
!                                  move usfc, vsfc from lw3 to lwm file
!       1994 nov 01  steve albers  pass additional arguments into laps_anl
!       1995 aug     steve albers  some of the u and v arrays are now
!                                  equivalenced to make passing them around
!                                  into a generic barnes routine easier.
!       1995 nov     k. brewster   changes for nyquist velocity array
!                                  and sounding profiles
!       1995 nov 28  steve albers  iradar_cycle_time = 600 (from 60)
!       1995 dec  8  steve albers  calls to separate ref and vel radar access
!                                  routines (in getradar.f)
!       1996 aug     steve albers  mods to allow option to rotate to grid
!                                  north.
!       1997 feb     steve albers  cleaned up some unneeded arrays/variables
!       1997 mar     steve albers  removed equiv statements.
!       1997 jun     ken dritz     removed nz_l_max from the argument list
!                                  in the call to get_fg_wind (it is not
!                                  actually used in get_fg_wind).
!       1997 jun     ken dritz     changed return statement to stop statement
!                                  in handling of istatus .ne. 1 after return
!                                  from call to get_domain_laps.
!       1997 jun     ken dritz     added r_missing_data as dummy argument.
!       1997 jun     ken dritz     removed parameter declaration of npts
!                                  as nx_l * ny_l (because it is unused).
!       1997 jun     ken dritz     removed parameter declaration of k_3d
!                                  as 2 * nz_l (because it is unused).
!       1997 jun     ken dritz     pass n_pirep in call to rdpirep.
!       1997 jun     ken dritz     pass nx_l and ny_l a second time to
!                                  vert_wind (once for dummy arguments ni,nj
!                                  and once for newly added dummy arguments
!                                  nx_l_max and ny_l_max).
!       1997 jun     ken dritz     pass r_missing_data to vert_wind.
!       1998 feb     steve albers  call rdpirep with cloud drift winds also

        use mem_namelist ! for these global parameters:
                         ! nx_l,ny_l,nk_laps,max_radars
                         ! ,r_missing_data,i2_missing_data
                         ! ,laps_cycle_time
                         ! 
                         ! plus a set of wind parameters:

        use mem_grid, only: lat,lon,topo

        use mem_wind

#ifdef single_exec
        use mem_background, only: bgfields%u3, bgfields%v3
#endif

        real, pointer, dimension(:,:,:) :: uanl, vanl, wanl 
        real, pointer, dimension(:,:) :: uanl_sfcitrp,vanl_sfcitrp 

        include 'barnesob.inc'
        type (barnesob) :: obs_point(max_wind_obs)                           

        include 'windparms.inc'

!       laps analysis grids
        real grid_laps_wt(nx_l,ny_l,nk_laps)
        real grid_laps_u(nx_l,ny_l,nk_laps)
        real grid_laps_v(nx_l,ny_l,nk_laps)

        include 'main_sub.inc'

!       housekeeping
        integer       ss_normal,rtsys_no_data
        parameter(ss_normal        =1, ! success
     1            rtsys_no_data    =3)

        character*3 exts(20),ext_in

        real rlat_radar(max_radars),rlon_radar(max_radars)
     1                               ,rheight_radar(max_radars)
        character*4 radar_name(max_radars)
        real v_nyquist_in(max_radars)
        integer n_vel_grids(max_radars),i4time_radar_a(max_radars)
        integer       ioffset(max_radars)
        integer       joffset(max_radars)
        logical       l_offset_radar 

!       stuff for call to laps analysis
!       real upass1(nx_l,ny_l,nk_laps),vpass1(nx_l,ny_l,nk_laps)

!        real uanl(nx_l,ny_l,nk_laps),vanl(nx_l,ny_l,nk_laps) ! wrt true north
!        real, allocatable, dimension(:,:,:) :: wanl

#ifdef multi_exec
        real, allocatable, dimension(:,:,:,:) ::
     +                     u_mdl_bkg_4d, v_mdl_bkg_4d
#else
        real, pointer,dimension(:,:,:,:) ::
     +                     u_mdl_bkg_4d, v_mdl_bkg_4d 
#endif

        dimension u_laps_fg(nx_l,ny_l,nk_laps),
     1            v_laps_fg(nx_l,ny_l,nk_laps)

        integer  n_3d_fields
        parameter (n_3d_fields = 3)

!       real outarray_4d(nx_l,ny_l,nk_laps,n_3d_fields)

        character*125 comment_2d,comment_a(2)
        character*10 units_2d,units_a(2)

        character*31 ext
        character*200 fname

!       character*31 ext_fg

!       stuff for sfc winds
!        real uanl_sfcitrp(nx_l,ny_l),vanl_sfcitrp(nx_l,ny_l)

        real, allocatable, dimension(:,:,:,:) :: grid_ra_vel
        real, allocatable, dimension(:,:,:,:) :: grid_ra_nyq
!       real grid_ra_vel(nx_l,ny_l,nk_laps,max_radars)
!       real grid_ra_nyq(nx_l,ny_l,nk_laps,max_radars)

        integer idx_radar_a(max_radars)
        character*31 ext_radar(max_radars)

        real heights_3d(nx_l,ny_l,nk_laps)
        real heights_1d(nk_laps)

!       local for laps_anl
        integer n_obs_lvl(nk_laps)

        integer i4_loop_total
        save i4_loop_total
        data i4_loop_total/0/

        istat_wind = 0

        nz_l = nk_laps
        
        ! set up pointers to mem_wind arrays so names do not need to be changed
        !    can not do this with the use rename list
        uanl => wind%uanl
        vanl => wind%vanl
        uanl_sfcitrp => wind%uanl_sfcitrp
        vanl_sfcitrp => wind%vanl_sfcitrp

        call get_directory('log',fname,len)
        open(15,file=fname(1:len)//'wind_stats.log'
     1      ,status='unknown',err=999)

        write(6,*)' welcome to the laps wind analysis'

        if(max_radars .eq. 0 .and. l_use_radial_vel)then
            write(6,*)' warning: max_radars = 0, '
     1               ,'setting l_use_radial_vel to .false.'
            l_use_radial_vel = .false.
        endif

!       get actual grid spacing valid at the gridpoint nearest the center
        icen = nx_l/2 + 1
        jcen = ny_l/2 + 1
        call get_grid_spacing_actual_xy(lat(icen,jcen),lon(icen,jcen)       
     1                        ,grid_spacing_actual_mx
     1                        ,grid_spacing_actual_my
     1                        ,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error return from get_grid_spacing_actual_xy'       
            stop
        endif

        grid_spacing_m = grid_spacing_actual_my

!       write(15,*)' i4time of data = ',i4time_sys
10      continue

        istat = init_timer()

        n_prods_out = 1 ! (any number > 0)

!       housekeeping
        n_prods = 7

        n_msg = 1
        n_pig = 2
        n_prg = 3
        n_sag = 4
        n_d00 = 5
        n_lw3 = 6
        n_lwm = 7

        exts(n_msg) = 'msg'
        exts(n_pig) = 'pig'
        exts(n_prg) = 'prg'
        exts(n_sag) = 'sag'
        exts(n_d00) = 'd00'
        exts(n_lw3) = 'lw3'
        exts(n_lwm) = 'lwm'

        write(6,*)' swlat ',lat(   1,   1)
        write(6,*)' nwlat ',lat(   1,ny_l)
        write(6,*)' selat ',lat(nx_l,   1)
        write(6,*)' nelat ',lat(nx_l,ny_l)

        write(6,*)' swlon ',lon(   1,   1)
        write(6,*)' nwlon ',lon(   1,ny_l)
        write(6,*)' selon ',lon(nx_l,   1)
        write(6,*)' nelon ',lon(nx_l,ny_l)

!       initialize grids with the missing data value
        do k = 1,nz_l
        do j = 1,ny_l
        do i = 1,nx_l
            grid_laps_u(i,j,k) = r_missing_data
            grid_laps_v(i,j,k) = r_missing_data
            grid_laps_wt(i,j,k) = r_missing_data
        enddo
        enddo
        enddo

!  ***  remap radial velocities to laps grid  *****************************
        i4_elapsed = ishow_timer()

!       i4_tol = max(laps_cycle_time / 2, iradar_cycle_time / 2)
        i4_tol = 900 ! seconds

        call get_l_offset_radar(nx_l,ny_l,grid_spacing_m,         ! i
     1                          nx_r,ny_r,igrid_r,l_offset_radar) ! o

        allocate(grid_ra_vel(nx_r,ny_r,nz_l,max_radars)
     1           ,stat=istat_alloc)             
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate grid_ra_vel'      
                        stop
        else
            write(6,*)' allocated grid_ra_vel ',nx_r,ny_r
        endif

        allocate(grid_ra_nyq(nx_r,ny_r,nz_l,max_radars)
     1           ,stat=istat_alloc)             
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate grid_ra_nyq'      
                        stop
        else
            write(6,*)' allocated grid_ra_nyq ',nx_r,ny_r
        endif

        if(l_use_radial_vel)then
            write(6,*)
            write(6,*)' reading radial velocity data'
            call get_multiradar_vel(i4time_sys,i4_tol                ! i
     1       ,i4time_radar_a                                         ! o
     1       ,max_radars                                             ! i
     1       ,n_radars                                               ! o
     1       ,ext_radar                                              ! o
     1       ,r_missing_data,nx_l,ny_l,nz_l,lat,lon                  ! i
     1       ,nx_r,ny_r,igrid_r                                      ! i
     1       ,grid_ra_vel,grid_ra_nyq,idx_radar_a,v_nyquist_in       ! o
     1       ,ioffset,joffset                                        ! o
     1       ,l_offset_radar                                         ! i
     1       ,n_vel_grids                                            ! o
     1       ,rlat_radar,rlon_radar,rheight_radar,radar_name         ! o
     1       ,istat_radar_vel,istat_radar_nyq)                       ! o

            if(n_radars .gt. 0)then
                write(6,*)
                write(6,5545)(radar_name(i),i=1,n_radars)
5545            format(' retrieved radar names:',30(1x,a4))
            endif

            if(istat_radar_vel .eq. 1)then
                write(6,*)' radar 3d vel data successfully read in'
     1                      ,(n_vel_grids(i),i=1,n_radars)
            else
                write(6,*)' radar 3d vel data not successfully read in'
     1                      ,(n_vel_grids(i),i=1,n_radars)
            endif

        else
            n_radars = 0
            istat_radar_vel = 0
            write(6,*)' not using radar 3d vel, l_use_radial_vel = '
     1               ,l_use_radial_vel

        endif ! l_use_radial_vel

        write(6,*)

!  ***  access model first guess*********************************************
        i4_elapsed = ishow_timer()

#ifdef multi_exec
        allocate (  u_mdl_bkg_4d(nx_l,ny_l,nz_l,ntmin:ntmax) )
        allocate (  v_mdl_bkg_4d(nx_l,ny_l,nz_l,ntmin:ntmax) )
         
        call get_fg_wind_new(
     1           i4time_sys,laps_cycle_time                    ! input
     1          ,nx_l,ny_l,nz_l                                ! input
     1          ,ntmin,ntmax                                   ! input
     1          ,u_mdl_bkg_4d,v_mdl_bkg_4d                     ! output
     1          ,u_laps_fg,v_laps_fg                           ! output
     1          ,istatus)                                      ! output 
        if(istatus .ne. 1)then
            write(6,*)' abort laps wind analysis - no first guess info'
            goto 999
        endif
#else
        u_mdl_bkg_4d => bgfields%u3
        v_mdl_bkg_4d => bgfields%v3
        u_laps_fg = bgfields%u3(:,:,:,0)
        v_laps_fg = bgfields%v3(:,:,:,0)
#endif

        write(6,*)' column of first guess winds'
        write(6,*)'j  v_laps_fg(23,j,13)'
        do j = 1,ny_l
            write(6,5555)j,v_laps_fg(23,j,13)
 5555       format(i4,4f7.2)
        enddo ! j

!  ***  read in wind obs **************************************************

        i4_elapsed = ishow_timer()

!       initialize observation data structure
        do i = 1,max_wind_obs
            obs_point(i)%l_withhold = .false.
        enddo ! i

        write(6,*)' calling get_wind_3d_obs'

        nobs_point = 0

        call get_wind_3d_obs(
     1            nx_l,ny_l,nz_l,                                 ! i
     1            r_missing_data,i2_missing_data,                 ! i
     1            i4time_sys,heights_3d,heights_1d,               ! i
     1            max_pr,max_pr_levels,weight_prof,l_use_raob,    ! i
     1            l_use_cdw,                                      ! i
     1            n_sao,n_pirep,                                  ! i
     1            lat,lon,topo,                                   ! i
     1            ntmin,ntmax,                                    ! i
     1            u_mdl_bkg_4d, v_mdl_bkg_4d,                     ! i
     1            grid_laps_u,grid_laps_v,grid_laps_wt,           ! o
     1            max_wind_obs,obs_point,nobs_point,              ! i/o
     1            rlat_radar(1),rlon_radar(1),rheight_radar(1),   ! i
     1            istat_radar_vel,n_vel_grids(1),                 ! i
     1            grid_ra_vel(1,1,1,1),                           ! i
     1            istatus_remap_pro,                              ! o
     1            istatus                )                        ! o

        if(istatus .ne. 1)then
            write(6,*)' abort laps wind analysis'
            goto 999
        endif

        if(istatus_remap_pro .ne. 1)then
            goto 999
        endif

        num_wind_obs = nobs_point

#ifdef multi_exec
        deallocate (u_mdl_bkg_4d, v_mdl_bkg_4d)
#endif

        i4_elapsed = ishow_timer()

        if((.not. l_grid_north_bkg) .and. l_grid_north_anal)then       
                write(6,*)' rotating first guess to grid north'

                do k = 1, nz_l
                do j = 1, ny_l
                do i = 1, nx_l
                    call uvtrue_to_uvgrid(
     1                           u_laps_fg(i,j,k),v_laps_fg(i,j,k)
     1                          ,u_grid   ,v_grid
     1                          ,lon(i,j)           )
                    u_laps_fg(i,j,k) = u_grid
                    v_laps_fg(i,j,k) = v_grid
                enddo
                enddo
                enddo

                i4_elapsed = ishow_timer()

        endif

        if(l_grid_north_anal)then       
                write(6,*)' rotating obs arrays to grid north'

                do k = 1, nz_l
                do j = 1, ny_l
                do i = 1, nx_l
                    if(     grid_laps_u(i,j,k) .ne. r_missing_data
     1                .and. grid_laps_v(i,j,k) .ne. r_missing_data )then
                        call uvtrue_to_uvgrid(
     1                           grid_laps_u(i,j,k),grid_laps_v(i,j,k)
     1                          ,u_grid   ,v_grid
     1                          ,lon(i,j)           )
                        grid_laps_u(i,j,k) = u_grid
                        grid_laps_v(i,j,k) = v_grid
                    endif
                enddo
                enddo
                enddo

                write(6,*)' rotating obs structure to grid north'
                do iob = 1,nobs_point
                    i = obs_point(iob)%i
                    j = obs_point(iob)%j
                    call uvtrue_to_uvgrid(obs_point(iob)%valuef(1)
     1                                   ,obs_point(iob)%valuef(2)
     1                                   ,u_grid   ,v_grid
     1                                   ,lon(i,j)           )
                    obs_point(iob)%valuef(1) = u_grid
                    obs_point(iob)%valuef(2) = v_grid
                enddo ! iob

                i4_elapsed = ishow_timer()

        endif ! l_grid_north

!sms$serial end

        call laps_anl(grid_laps_u,grid_laps_v
     1       ,obs_point,max_wind_obs,nobs_point                         ! i
     1       ,n_radars,istat_radar_vel                                  ! i
     1       ,grid_ra_vel,grid_ra_nyq,v_nyquist_in,idx_radar_a          ! i
!    1       ,upass1,vpass1                                             ! o
     1       ,n_var                                                     ! i
     1       ,uanl,vanl                                                 ! o
     1       ,grid_laps_wt,weight_bkg_const_wind,rms_thresh_wind        ! i/l
     1       ,r0_barnes_max_m,brns_conv_rate_wind                       ! i
     1       ,max_radars                                                ! i
     1       ,n_vel_grids,rlat_radar,rlon_radar,rheight_radar           ! i
     1       ,thresh_2_radarobs_lvl_unfltrd                             ! i
     1       ,thresh_4_radarobs_lvl_unfltrd                             ! i
     1       ,thresh_9_radarobs_lvl_unfltrd                             ! i
     1       ,thresh_25_radarobs_lvl_unfltrd                            ! i
     1       ,u_laps_fg,v_laps_fg                                       ! i/l
     1       ,nx_l,ny_l,nz_l,lat,lon                                    ! i
     1       ,nx_r,ny_r,ioffset,joffset                                 ! i
     1       ,i4time_sys,grid_spacing_m                                 ! i
     1       ,r_missing_data                                            ! i
     1       ,heights_3d                                                ! i
     1       ,i4_loop_total                                             ! o
     1       ,l_derived_output,l_grid_north_anal,l_3pass                ! i
     1       ,l_correct_unfolding                                       ! i
!    1       ,n_iter_wind
     1       ,weight_cdw,weight_sfc,weight_pirep,weight_prof            ! i
     1       ,weight_radar                                              ! i
     1       ,istatus)                                                  ! o

!sms$serial(default=ignore)  begin              

        if(istatus .ne. 1)then
                write(6,*)' error in wind analysis (laps_anl)'
                goto 999
        endif

        if(l_grid_north .and. .not. l_grid_north_out)then

                i4_elapsed = ishow_timer()

                write(6,*)' rotating analyzed wind back to true north'
                do k = 1, nz_l
                do j = 1, ny_l
                do i = 1, nx_l
                        call uvgrid_to_uvtrue(
     1                           uanl(i,j,k),vanl(i,j,k)
     1                          ,u_true   ,v_true
     1                          ,lon(i,j)           )
                    uanl(i,j,k) = u_true
                    vanl(i,j,k) = v_true
                enddo
                enddo
                enddo

                i4_elapsed = ishow_timer()

        endif ! rotate back to true north

        write(6,*)'uanl(nx_l/2+1,ny_l/2+1,1) = '
     1            ,uanl(nx_l/2+1,ny_l/2+1,1)
        write(6,*)'vanl(nx_l/2+1,ny_l/2+1,1) = '
     1            ,vanl(nx_l/2+1,ny_l/2+1,1)

        i4_elapsed = ishow_timer()

#ifdef multi_exec
!       need to allocate actual datatype member, can not use wanl pointer for 
!       allocation
        allocate( wind%wanl(nx_l,ny_l,nk_laps), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate wind%wanl'
            stop
        endif
#endif
        wanl => wind%wanl

        call wind_post_process(i4time_sys
     1                        ,uanl,vanl                            ! i
     1                        ,wanl                                 ! o
     1                        ,nx_l,ny_l,nz_l,n_3d_fields           ! i
     1                        ,heights_3d                           ! i
     1                        ,uanl_sfcitrp,vanl_sfcitrp            ! i
     1                        ,topo,lat,lon,grid_spacing_m          ! i
     1                        ,r_missing_data,l_grid_north_out      ! i
     1                        ,istat_lw3)

        i4_elapsed = ishow_timer()

        if(i4_elapsed .gt. 0)then
            pct_loop_total = 
     1          float(i4_loop_total)/float(i4_elapsed) * 100.
        else
            pct_loop_total = 0.
        endif
        write(6,998)           i4_loop_total                 
     1                        ,i4_elapsed - i4_loop_total
     1                        ,i4_elapsed     
     1                        ,pct_loop_total
998     format(' i4time (final) loop/nonloop/total/% = ',3i5,f8.2)

        istat_wind = 1

        write(6,*)' end of laps wind analysis'

        close(15)
!sms$serial end

 999    continue
   
        deallocate (grid_ra_vel, grid_ra_nyq)

        return

        end

