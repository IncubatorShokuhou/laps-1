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
  

        subroutine get_wind_3d_obs(
     1            nx_l,ny_l,nz_l,                                   ! i
     1            r_missing_data,i2_missing_data,                   ! i
     1            i4time_lapswind,heights_3d,heights_1d,            ! i
     1            max_pr,max_pr_levels,weight_prof,l_use_raob,      ! i
     1            l_use_cdw,                                        ! i
     1            n_sfc,n_pirep,                                    ! i
     1            lat,lon,topo,                                     ! i
     1            ntmin,ntmax,                                      ! i
     1            u_mdl_bkg_4d, v_mdl_bkg_4d,                       ! i
     1            grid_laps_u,grid_laps_v,grid_laps_wt,             ! i/o
     1            max_obs,obs_point,nobs_point,                     ! i/o
     1            rlat_radar,rlon_radar,rheight_radar,              ! i
     1            istat_radar_vel,n_vel_grids,                      ! i
     1            istatus_remap_pro,                                ! o
     1            istatus                )                          ! o

!       1997 jun     ken dritz       added nx_l, ny_l, nz_l as dummy arguments,
!                                    making arrays with those dimensions
!                                    automatic.
!       1997 jun     ken dritz       added r_missing_data and i2_missing_data
!                                    as dummy arguments.
!       1997 jun     ken dritz       removed include of 'lapsparms.for'.

        include 'barnesob.inc'
        type (barnesob) :: obs_point(max_obs)                           

!       laps grid dimensions
                                                   
        real lat(nx_l,ny_l)                      
        real lon(nx_l,ny_l)                      
        real topo(nx_l,ny_l)

!       profiler stuff
        real lat_pr(max_pr)                        
        real lon_pr(max_pr)                        
        character*8 obstype(max_pr)
        character*5 c5_name_a(max_pr)
        integer i4time_ob_pr(max_pr)

!       profiler observations

        integer nlevels_obs_pr(max_pr)             
        real ob_pr_u (max_pr,nz_l) ! vertically interpolated profiler wind
        real ob_pr_v (max_pr,nz_l) ! vertically interpolated profiler wind

!       laps analysis grids
        real grid_laps_wt(nx_l,ny_l,nz_l)                               
        real grid_laps_u(nx_l,ny_l,nz_l)                                
        real grid_laps_v(nx_l,ny_l,nz_l)                                

!       data flags for laps analysis
        logical l_profiler
        parameter (l_profiler = .true.)

        real u_mdl_bkg_4d(nx_l,ny_l,nz_l,ntmin:ntmax)
        real v_mdl_bkg_4d(nx_l,ny_l,nz_l,ntmin:ntmax)

        real u_laps_fg(nx_l,ny_l,nz_l)
        real v_laps_fg(nx_l,ny_l,nz_l)

        real heights_3d(nx_l,ny_l,nz_l)
        real heights_1d(nz_l)

        character*3 ext_in

        logical l_use_raob,l_use_cdw,l_use_all_nontower_lvls

        l_use_all_nontower_lvls = .false.

!  ***  read in model first guess data  **************************************

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' error getting laps_cycle_time'
            return
        endif

!  ***  read in profiler data  ********************************************

        write(6,*)' calling read_profiles'

        call read_profiles(
     1            i4time_lapswind,heights_3d,                       ! i
     1            lat_pr,lon_pr,obstype,c5_name_a,                  ! o
     1            lat,lon,topo,i4time_ob_pr,                        ! i
     1            max_pr,max_pr_levels,                             ! i
     1            l_use_raob,l_use_all_nontower_lvls,               ! i
     1            ob_pr_u , ob_pr_v ,                               ! o
     1            max_obs,obs_point,nobs_point,weight_prof,         ! i/o
     1            nlevels_obs_pr,n_profiles,                        ! o
     1            rlat_radar,rlon_radar,rheight_radar,              ! i
     1            n_vel_grids,                                      ! i
     1            u_mdl_bkg_4d,v_mdl_bkg_4d,ntmin,ntmax,            ! i
     1            ilaps_cycle_time,r_missing_data,                  ! i
     1            nx_l,ny_l,nz_l,                                   ! i
     1            istatus                )                          ! o

        if(istatus .ne. 1)then
            write(6,*)' abort read_profiles'
            return
        endif

        i4_elapsed = ishow_timer()

! ***   remapping + barnes analysis of profiler data in u & v *****************

        call remap_profiles(
     1           ob_pr_u,ob_pr_v                                  ! i
     1          ,grid_laps_u,grid_laps_v,grid_laps_wt             ! i/o
     1          ,max_obs,obs_point,nobs_point                     ! i/o
     1          ,lat,lon                                          ! i
     1          ,nx_l,ny_l,nz_l,max_pr                            ! i
     1          ,nlevels_obs_pr,lat_pr,lon_pr,obstype,n_profiles  ! i
     1          ,c5_name_a,i4time_ob_pr                           ! i
     1          ,r_missing_data                                   ! i
     1          ,weight_prof                                      ! o
     1          ,l_profiler,l_use_raob,l_use_all_nontower_lvls    ! i
     1          ,istatus_remap_pro)                               ! o

!  ***  read in sfc wind data   *******************************************

	call get_maxstns(maxstns,istatus)
	if (istatus .ne. 1) then
	   write (6,*) 'error obtaining maxstns'
           return
	endif

        call rdsfc(i4time_lapswind,heights_3d                            ! i
     1            ,n_sfc,maxstns                                         ! i
     1            ,lat,lon,r_missing_data                                ! i
     1            ,n_sfc_obs                                             ! o
     1            ,grid_laps_wt,grid_laps_u,grid_laps_v                  ! i/o
     1            ,max_obs,obs_point,nobs_point                          ! i/o
     1            ,nx_l,ny_l,nz_l                                        ! i
     1            ,istatus)                                              ! o
        if(istatus .ne. 1)then
            write(6,*)
     1     ' aborting from get_wind_obs - error reading sfc obs'
            return
!       else
!           i4time_array(n_sag) = i4time_lapswind
!           j_status(n_sag) = ss_normal
        endif

        i4_elapsed = ishow_timer()

!  ***  read in pirep data   ********************************************

        n_pirep_obs = 0

        ext_in = 'pin'

        call rdpoint(i4time_lapswind,heights_3d
     1  ,n_pirep,n_pirep_obs,ext_in
     1  ,nx_l,ny_l,nz_l
     1  ,u_mdl_bkg_4d,v_mdl_bkg_4d,ntmin,ntmax                         ! i
     1  ,lat,lon
!    1  ,pirep_i,pirep_j,pirep_k,pirep_u,pirep_v
     1  ,grid_laps_wt,grid_laps_u,grid_laps_v                          ! i/o
     1  ,max_obs,obs_point,nobs_point                                  ! i/o
     1  ,istatus)                                                      ! o

        if(istatus .ne. 1)then
            write(6,*)
     1       ' aborting from get_wind_obs - error reading pireps'
            return
!       else
!           i4time_array(n_pig) = i4time_lapswind
!           j_status(n_pig) = ss_normal
        endif

!  ***  read in cloud drift wind data   ***********************************

        if(l_use_cdw)then

            ext_in = 'cdw'

            call rdpoint(i4time_lapswind,heights_3d
     1      ,n_pirep,n_pirep_obs,ext_in
     1      ,nx_l,ny_l,nz_l
     1      ,u_mdl_bkg_4d,v_mdl_bkg_4d,ntmin,ntmax                     ! i
     1      ,lat,lon
!    1      ,pirep_i,pirep_j,pirep_k,pirep_u,pirep_v
     1      ,grid_laps_wt,grid_laps_u,grid_laps_v                      ! i/o
     1      ,max_obs,obs_point,nobs_point                              ! i/o
     1      ,istatus)                                                  ! o

            if(istatus .ne. 1)then
                write(6,*)
     1           ' aborting from laps wind anal - error reading cdw'
                return
!           else
!               i4time_array(n_pig) = i4time_lapswind
!               j_status(n_pig) = ss_normal
            endif

        endif ! l_use_cdw


        return
        end


        subroutine remap_profiles(
     1           ob_pr_u,ob_pr_v                                     ! i
     1          ,grid_laps_u,grid_laps_v,grid_laps_wt                ! i/o
     1          ,max_obs,obs_point,nobs_point                        ! i/o
     1          ,lat,lon                                             ! i
     1          ,ni,nj,nk,max_pr                                     ! i
     1          ,nlevels_obs_pr,lat_pr,lon_pr,obstype,n_profiles     ! i
     1          ,c5_name_a,i4time_ob_pr                              ! i
     1          ,r_missing_data                                      ! i
     1          ,weight_prof                                         ! o
     1          ,l_profiler,l_use_raob,l_use_all_nontower_lvls       ! i
     1          ,istatus)                                            ! o

!       perform horizontal remapping of profile obs onto laps grid
!       they have already been vertically interpolated
!	2006	yuanfu xie	use of the fraction grid values of obs_point.

        include 'barnesob.inc'
        type (barnesob) :: obs_point(max_obs)                           

!       profile observations

        integer nlevels_obs_pr(max_pr)
        real lat_pr(max_pr)
        real lon_pr(max_pr)
        character*8 obstype(max_pr)
        character*12 c_obstype
        character*5 c5_name_a(max_pr)
        character*9 a9time
        integer i4time_ob_pr(max_pr)
        real ob_pr_u (max_pr,nk) ! vertically interpolated profile wind
        real ob_pr_v (max_pr,nk) ! vertically interpolated profile wind

!       barnes profile analysis

        real lat(ni,nj),lon(ni,nj)

        real grid_laps_u(ni,nj,nk)
        real grid_laps_v(ni,nj,nk)
        real grid_laps_wt(ni,nj,nk)

        logical l_profiler, l_use_all_nontower_lvls, l_use_raob

        write(6,*)
        write(6,*)' subroutine remap_profiles: # of profiles = '
     1           ,n_profiles
        write(6,*)' u/v are true north with time trending applied...'

        do i_pr = 1,n_profiles ! max_pr

            if(.not. l_use_all_nontower_lvls  .and.
     1        obstype(i_pr)(1:5) .ne. 'tower' .and.
     1        obstype(i_pr)(1:5) .ne. 'sodar'      )then       

                if(nlevels_obs_pr(i_pr) .gt. 0)then
                    call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr)
     1                                      ,lat,lon,ni,nj,ri,rj
     1                                      ,istatus)     
                    if(istatus .ne. 1)then
                        write(6,*)' note... profile is outside domain'  
     1                           ,i_pr,i_ob,j_ob,' ',obstype(i_pr)

                    else ! inside domain
                        i_ob = nint(ri)
                        j_ob = nint(rj)

                        call make_fnam_lp(i4time_ob_pr(i_pr),a9time
     1                                   ,istatus)
                        if(istatus .ne. 1)then
                            write(6,*)
     1                      ' error in remap_profiles - invalid obtime'       
                            return
                        endif
              
                        write(6,311)i_pr,i_ob,j_ob,nlevels_obs_pr(i_pr)
     1                           ,obstype(i_pr),c5_name_a(i_pr),a9time
 311                    format(1x,' remapping profile ',4i6,1x,a8,1x,a5
     1                        ,1x,a9,' (intrp laps lvls)')      

                        do k = 1,nk
                            if(ob_pr_u(i_pr,k) .ne. r_missing_data)then

                                ob_u = ob_pr_u (i_pr,k)
                                ob_v = ob_pr_v (i_pr,k)

!                 ***           map observation onto laps grid   ***
                                if(l_profiler)then
!                                   grid_laps_u(i_ob,j_ob,k) = ob_u
!                                   grid_laps_v(i_ob,j_ob,k) = ob_v
!                                   grid_laps_wt(i_ob,j_ob,k) = weight_prof

!                                   add to data structure (still is subsampling)
                                    nobs_point = nobs_point + 1

                                    write(6,11)k,nobs_point,ob_u,ob_v
 11                                 format(10x,i4,i6,2f8.1)

                                    obs_point(nobs_point)%i = i_ob
                                    obs_point(nobs_point)%j = j_ob
                                    obs_point(nobs_point)%k = k
                                    obs_point(nobs_point)%ri = ri    ! yuanfu
                                    obs_point(nobs_point)%rj = rj    ! yuanfu
                                    obs_point(nobs_point)%rk = k
                                    obs_point(nobs_point)%valuef(1)=ob_u       
                                    obs_point(nobs_point)%valuef(2)=ob_v
                                    obs_point(nobs_point)%weight = 
     1                                                    weight_prof       
                                    obs_point(nobs_point)%vert_rad_rat =
     1                                                    1.0   
                                    call downcase(obstype(i_pr)
     1                                           ,c_obstype)    
                                    obs_point(nobs_point)%type = 
     1                                            c_obstype
                                    obs_point(nobs_point)%file = 'pro'      
                                    obs_point(nobs_point)%i4time =
     1                                               i4time_ob_pr(i_pr)       
                                    if(obstype(i_pr)(1:4) 
     1                                                  .eq. 'raob')then
                                      if(.not. l_use_raob)then
                                        write(6,*)' withholding this ob'
                                        obs_point(nobs_point)%l_withhold       
     1                                      = .true.
                                      endif
                                    endif

                                endif ! l_profiler

                            endif ! not missing data
                        enddo ! k
                    endif ! istatus (in/out of domain)

                else
                    if(i_pr .le. 100)then
                        write(6,*)' zero levels / outside domain'
     1                           ,i_pr,i_ob,j_ob,' ',obstype(i_pr)
                    endif

                endif ! data present

            endif ! remapping with vertical interpolation

        enddo ! i_pr

        i4_elapsed = ishow_timer()

        istatus = 1

        return
        end


