cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis

        subroutine get_vis(i4time,solar_alt,l_use_vis,l_use_vis_add ! I
     1                    ,l_use_vis_partial,lat                    ! I
     1                    ,i4_sat_window,i4_sat_window_offset       ! I
     1                    ,rlaps_land_frac,topo                     ! I
     1                    ,cloud_frac_vis_a,vis_albedo,ihist_alb    ! O
     1                    ,static_albedo                            ! O
     1                    ,comment                                  ! O
     1                    ,ni,nj,nk,r_missing_data                  ! I
     1                    ,istat_vis_a                              ! O
     1                    ,istatus)                                 ! O

!       Steve Albers 1997

        integer*4 ihist_alb(-10:20)
        integer*4 ihist_alb_sfc(-10:20)
        integer*4 ihist_frac_sat(-10:20)
        integer*4 istat_vis_a(ni,nj)            ! Cloud mask based on VIS

        real*4 lat(ni,nj)
        real*4 sfc_albedo(ni,nj), sfc_albedo_lwrb(ni,nj)
        real*4 static_albedo(ni,nj)   ! Static albedo database
        real*4 vis_albedo(ni,nj)
        real*4 rlaps_land_frac(ni,nj)
        real*4 topo(ni,nj)

!       This stuff is for reading VIS data from LVD file
        real*4 solar_alt(ni,nj)
        real*4 cloud_frac_vis_a(ni,nj)
        integer*4 mxstn
        parameter (mxstn = 100)       ! max number of "stations" in data file
        character*9 filename
        character*31 ext

        character*3 lvd_ext
        data lvd_ext /'lvd'/

        character var*3,comment*125,units*10
        logical l_use_vis, l_use_vis_add, l_use_vis_partial

!       l_use_vis_add = .false.
        icount_vis_add_potl = 0

!       Initialize histograms
        do i = -10,20
            ihist_alb(i) = 0
            ihist_alb_sfc(i) = 0
            ihist_frac_sat(i) = 0
        enddo ! i

!       Initialize
        do i = 1,ni
        do j = 1,nj
            vis_albedo(i,j) = r_missing_data
            cloud_frac_vis_a(i,j) = r_missing_data
            istat_vis_a(i,j) = 0
        enddo ! j
        enddo ! i

        call get_sfc_albedo(ni,nj,lat,r_missing_data,i4time              ! I
     1                     ,rlaps_land_frac,topo                         ! I
     1                     ,sfc_albedo,sfc_albedo_lwrb                   ! O
     1                     ,static_albedo,istat_sfc_alb)                 ! O

!       Determine whether to use VIS / ALBEDO data
        if(.not. l_use_vis)then
            write(6,*)' Note: l_use_vis set to not use vis data'
            istatus = 0
            return
        endif

!       Read in albedo data
        write(6,*)' Getting the VIS data from LVD file'
        call make_fnam_lp(i4time,filename,istatus)

        ext = lvd_ext
        var = 'ALB'
        ilevel = 0
        call get_laps_2dvar(i4time+i4_sat_window_offset,i4_sat_window
     1                     ,i4time_nearest,EXT,var,units
     1                     ,comment,ni,nj,vis_albedo,ilevel,istatus)
        write(6,*)' istatus from vis_albedo data = ',istatus
        if(istatus .ne. 1 .and. istatus .ne. -1)then
            write(6,*)' No VIS / ALBEDO available'
            write(6,*)' return from get_vis'
            istatus = 0
            return
        endif

!       Initial test for missing albedo (and partial data coverage)
        do i = 1,ni
        do j = 1,nj
            if(vis_albedo(i,j) .eq. r_missing_data .and. 
     1                          (.not. l_use_vis_partial)      )then
                write(6,*)
     1              ' No VIS / ALBEDO available (missing data found)'          
                write(6,*)' return from get_vis'
                vis_albedo = r_missing_data
                istatus = 0
                return
            endif
        enddo ! j
        enddo ! i

        n_missing_albedo = 0

!       Horizontal array loop
        do i = 1,ni
        do j = 1,nj

!         We will now only use the VIS data if the solar alt exceeds 15 deg
!         7 degrees now used to allow 30 min slack in data availability
          if(solar_alt(i,j) .lt. 7.0)then
              if(vis_albedo(i,j) .ne. r_missing_data)then
                  write(6,*)' Error -  vis_albedo not missing:'
     1                     ,solar_alt(i,j)
                  stop
              endif
!             vis_albedo(i,j) = r_missing_data
          endif

          if(vis_albedo(i,j) .ne. r_missing_data)then

!           Translate the vis_albedo into cloud fraction

!           Store histogram information for satellite data
            iscr_alb  = nint(vis_albedo(i,j)*10.)
            iscr_alb  = min(max(iscr_alb,-10),20)
            ihist_alb(iscr_alb) = ihist_alb(iscr_alb) + 1

            if(sfc_albedo_lwrb(i,j) .ne. r_missing_data)then
                iscr_alb_sfc  = nint(sfc_albedo_lwrb(i,j)*10.)
                iscr_alb_sfc  = min(max(iscr_alb_sfc,-10),20)
                ihist_alb_sfc(iscr_alb_sfc) = 
     1          ihist_alb_sfc(iscr_alb_sfc) + 1
            endif

            clear_albedo = sfc_albedo_lwrb(i,j)
            cloud_frac_vis = albedo_to_cloudfrac2(clear_albedo
     1                                           ,vis_albedo(i,j))

            iscr_frac_sat = nint(cloud_frac_vis*10.)
            iscr_frac_sat = min(max(iscr_frac_sat,-10),20)
            ihist_frac_sat(iscr_frac_sat) = 
     1      ihist_frac_sat(iscr_frac_sat) + 1

!           Make sure satellite cloud fraction is between 0 and 1
            if(cloud_frac_vis .le. 0.0)cloud_frac_vis = 0.0
            if(cloud_frac_vis .ge. 1.0)cloud_frac_vis = 1.0

            cloud_frac_vis_a(i,j) = cloud_frac_vis

!           Is there enough of a signal from the VIS to say a cloud is present?
            if(       cloud_frac_vis_a(i,j) .gt. 0.3 
     1          .and. sfc_albedo(i,j) .ne. r_missing_data
     1          .and. sfc_albedo(i,j) .le. 0.3              
     1          .and. l_use_vis_add                         )then
                istat_vis_a(i,j) = 1
                icount_vis_add_potl = icount_vis_add_potl + 1
            endif

          else
            n_missing_albedo =  n_missing_albedo + 1
            cloud_frac_vis_a(i,j) = r_missing_data

          endif

        enddo ! i
        enddo ! j

        write(6,*)
        write(6,*)' N_MISSING_ALBEDO = ',n_missing_albedo
        write(6,*)

        write(6,*)' Number of potential visible clear/add = '
     1            ,ni*nj-n_missing_albedo,icount_vis_add_potl

        if(n_missing_albedo .eq. ni*nj)then ! Return with status = 0
            write(6,*)' All albedos were missing - return from get_vis'
            istatus = 0
            return
        endif

        write(6,*)
        write(6,*)'              HISTOGRAMS'
        write(6,*)' I       ',
     1  ' Sat Albedo  Sfc Albedo  Cld Frac Sat'
        do i = -5,15
            write(6,11)i,ihist_alb(i),ihist_alb_sfc(i),ihist_frac_sat(i)       
11          format(i4,i12,i12,i12)
        enddo ! i

        write(6,*)
        istatus = 1

        return
        end

        subroutine get_sfc_albedo(ni,nj,lat,r_missing_data,i4time    ! I
     1                           ,rlaps_land_frac,topo               ! I
     1                           ,sfc_albedo,sfc_albedo_lwrb         ! O
     1                           ,static_albedo                      ! O
     1                           ,istat_sfc_alb)                     ! O

!       This returns the surface albedo. This is from the static database
!       to yield a less confident "lower bound". If we are confident that
!       snow/ice is not present, then the 'sfc_albedo' array is also 
!       populated and can be used more at face value.

        character*3 var
        real*4 lat(ni,nj)
        real*4 sfc_albedo(ni,nj)      ! Populated with "reliable" values that
                                      ! may include land/sea snow/ice

        real*4 sfc_albedo_lwrb(ni,nj) ! Populated with lower bound (i.e. from
                                      ! static database)

        real*4 static_albedo(ni,nj)   ! Static albedo database

        real*4 rlaps_land_frac(ni,nj)
        real*4 topo(ni,nj)

        write(6,*)' Subroutine get_sfc_albedo...'

        istat_sfc_alb = 0
        static_albedo = r_missing_data

        call get_static_field_interp('albedo',i4time,ni,nj
     1                               ,static_albedo,istat_sfc_alb)       
        if(istat_sfc_alb .ne. 1)then ! Read sfc albedo from fixed database
            write(6,*)' Monthly Albedo Data N/A, look for fixed data'
            var = 'ALB'
            call read_static_grid(ni,nj,var,static_albedo,istat_sfc_alb)
        endif

        icount_albedo = 0
        icount_albedo_lwrb = 0

        do i = 1,ni
        do j = 1,nj
            if(istat_sfc_alb .eq. 1 .and. 
     1         static_albedo(i,j) .ne. r_missing_data)then ! static data avalbl
                sfc_albedo_lwrb(i,j) = static_albedo(i,j)

!               Consider lat/topo now - later also SST (tgd) > 5C ? 

                if(rlaps_land_frac(i,j) .le. 0.25)then         ! over water
                  if((abs(lat(i,j)) .le. 40. .and. topo(i,j) .le.  100.)
     1                                     .OR.
     1               (abs(lat(i,j)) .le. 20. .and. topo(i,j) .le. 3000.)
     1                                          )then          ! it's reliable
                      sfc_albedo(i,j) = static_albedo(i,j)
                  else                               ! sea/lake ice possible?
                      sfc_albedo(i,j) = r_missing_data
                  endif

                else                                           ! over land
                  sfc_albedo(i,j)      = r_missing_data

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

        write(6,*)' Number of sfc albedo and lower bound points = '       
     1           ,icount_albedo,icount_albedo_lwrb

        return
        end

        function albedo_to_cloudfrac2(clear_albedo,albedo)

        cloud_albedo = .4485300

        arg = albedo

        call stretch2(clear_albedo,cloud_albedo,0.,1.,arg)

        albedo_to_cloudfrac2 = arg

        return
        end
