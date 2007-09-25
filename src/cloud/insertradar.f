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
        subroutine insert_radar(i4time,cldcv,cld_hts
     1         ,temp_3d,temp_sfc_k,td_sfc_k                          ! I
     1         ,grid_spacing_m,ni,nj,nk,kcloud                       ! I
     1         ,cloud_base,ref_base                                  ! I
     1         ,topo,r_missing_data                                  ! I
     1         ,grid_ra_ref,dbz_max_2d                               ! I/O
     1         ,vis_radar_thresh_dbz                                 ! I
     1         ,l_unresolved                                         ! O
     1         ,heights_3d                                           ! I
     1         ,istatus)                                             ! O

!       Insert radar data into cloud grid and reconcile cloud/radar fields.

!       For wideband radar, echoes can be added in as cloud if the echo top 
!       is >2000m AGL and the echo is located above an existing "nominal" 
!       cloud base. The nominal cloud base is taken (if present) from the 
!       pre-existing cloud field, either at the horizontal grid location 
!       in question or from a neighboring grid point that has existing cloud.

!       High echo tops with consistent cloud bases are treated as "resolved".
!       Low echo tops are treated as "unresolved" as are echoes that lack
!       pre-existing cloud.

!       For all "unresolved" echoes no cloud is added and the echo is blanked 
!       out.

!       Narrowband radar is treated similarly, except that the echo top
!       threshold does not apply. Instead, this subroutine will be actively
!       adding clouds (between the pre-existing cloud layers) only for 
!       "trusted" narrowband radars such as NOWRAD. For other narrowband 
!       radars, the radar is blanked out if there is no pre-existing cloud 
!       in the calling routine 'laps_cloud_sub'.

!       Further reconciling of radar and satellite/cloud information is done 
!       later on in subroutine 'insert_vis'.

!       1992 - 2003        Steve Albers

        real temp_3d(ni,nj,nk)
        real heights_3d(ni,nj,nk)
        real temp_sfc_k(ni,nj)
        real td_sfc_k(ni,nj)

        real cldcv(ni,nj,kcloud)
        real cld_hts(kcloud)

        real grid_ra_ref(ni,nj,nk)
        real dbz_max_2d(ni,nj)
        real cloud_base(ni,nj)
        real cloud_base_buf(ni,nj) ! Lowest SAO/IR base within search radius
        real topo(ni,nj)

        real echo_top(ni,nj)           ! L
        real echo_top_agl(ni,nj)       ! L
        real lcl_agl(ni,nj)            ! L

!       Cloud not filled in unless radar echo is higher than base calculated
!       with THIS threshold.
        real     thresh_cvr
        parameter (thresh_cvr = 0.10)

        logical l_below_base, l_inserted
        logical l_unresolved(ni,nj)

!       Calculate Cloud Bases
        unlimited_ceiling = 200000.

!       Echo top threshold
        echo_agl_thr = 2000.

!       Search radius for cloud bases
        search_radius_m = 120000.

        do j = 1,nj
        do i = 1,ni
            cloud_base(i,j) = unlimited_ceiling
            echo_top(i,j) = r_missing_data
            echo_top_agl(i,j) = r_missing_data

            do k = kcloud-1,1,-1
                if(cldcv(i,j,k  ) .lt. thresh_cvr .and.
     1             cldcv(i,j,k+1) .ge. thresh_cvr       )then
                    cloud_base(i,j) = 0.5 * (cld_hts(k) + cld_hts(k+1))
                endif
            enddo ! k

            cloud_base_buf(i,j) = cloud_base(i,j)

            l_unresolved(i,j) = .false.

!           Approximation based on a 4 deg F dewpoint depression for each
!           1000 ft of cloud base above the ground
            lcl_agl(i,j) = (temp_sfc_k(i,j) - td_sfc_k(i,j)) * 137.16 ! meters

        enddo ! i
        enddo ! j

        icount_below = 0
        isearch_base = 0
        insert_count_tot = 0
        iwrite_inserted = 0

        isearch_radius = nint(search_radius_m / grid_spacing_m)
        intvl = max((isearch_radius / 4),1)

        write(6,*)' isearch_radius/intvl = ',isearch_radius,intvl

        do k = nk,1,-1 ! Essential that this go downward to detect radar tops
                       ! in time to search for a new cloud base

            kp1 = min(k+1,nk)

            icount_radar_lvl = 0
            insert_count_lvl = 0

            klow = 10000
            khigh = -1

           !Find range of cloud grid levels for this LAPS grid level
           !Find lowest level
            do kk = 1,kcloud
c               write(6,*)k,kk,height_to_zcoord(cld_hts(kk),istatus)
                if(nint(height_to_zcoord(cld_hts(kk),istatus)) .eq. k
     1                                                        )then
c                   write(6,*)' klow = ',kk
                    klow = kk
                    goto400
                endif
            enddo
400         do kk = kcloud,1,-1
                if(nint(height_to_zcoord(cld_hts(kk),istatus)) .eq. k
     1                                                        )then
c                   write(6,*)' khigh = ',kk
                    khigh = kk
                    goto500
                endif
            enddo

            if(istatus .ne. 1)then
                write(6,*)' Error in insert_radar'
                return
            endif

            if(klow .eq. 10000 .or. khigh .eq. -1)goto600

500         do i = 1,ni
            do j = 1,nj

                l_inserted = .false.

                if(temp_3d(i,j,k) .lt. 278.)then
                    ref_thresh = 0.0
                else
                    ref_thresh = 0.0 ! Lowered from 20. and 13. 9/94
                endif

                if(grid_ra_ref(i,j,k) .gt. ref_thresh)then
                    icount_radar_lvl = icount_radar_lvl + 1
                    l_below_base = .false.

!                   Test if we are at echo top
                    if(k .eq. nk   .or. 
     1                 grid_ra_ref(i,j,kp1) .lt. ref_thresh)then

                        echo_top(i,j) = heights_3d(i,j,k)
                        echo_top_agl(i,j) = echo_top(i,j) - topo(i,j)

!                       Test if we are below the cloud base
                        if(echo_top(i,j) .lt. cloud_base(i,j))then

!                           Radar Echo Top is below analyzed cloud base
!                           Search for Modified Cloud Base, i.e. other neighboring
!                           bases lower than the current buffer value

                            ilow =  max(i-isearch_radius,1)
                            ihigh = min(i+isearch_radius,ni)
                            jlow =  max(j-isearch_radius,1)
                            jhigh = min(j+isearch_radius,nj)

                            do jj = jlow,jhigh,intvl
                            do ii = ilow,ihigh,intvl
                                cloud_base_buf(i,j)
     1                    = min(cloud_base(ii,jj),cloud_base_buf(i,j))
                            enddo ! ii
                            enddo ! jj

                            if(cloud_base_buf(i,j) .lt. echo_top(i,j)
     1                                       .AND. 
     1                         echo_top_agl(i,j) .gt. echo_agl_thr)then       

                              isearch_base = isearch_base + 1
                              if(isearch_base .lt. 50)then ! limit log output
                                write(6,71)i,j,k
     1                                    ,nint(echo_top(i,j))
     1                                    ,nint(cloud_base(i,j))
     1                                    ,nint(cloud_base_buf(i,j))
71                              format(' Rdr Top > Bse ',2i4,i3,3i7
     1                                ,' Resolved')
                              endif

                            else ! Potentially Unresolved base
                                if(cloud_base(i,j) 
     1                                      .eq. unlimited_ceiling    ! No clds
     1                                      .OR.
     1                             echo_top_agl(i,j) .lt. echo_agl_thr! Gnd Clut
     1                                                             )then       

!                                   We will want to reconcile cloud/radar
                                    l_unresolved(i,j) = .true.

                                    write(6,72)i,j,k
     1                                    ,nint(echo_top(i,j))
     1                                    ,nint(cloud_base(i,j))
     1                                    ,nint(cloud_base_buf(i,j))
72                                  format(' Rdr Top < Bse ',2i4,i3,3i7
     1                                    ,' Unresolved      - CLD_RDR')

                                else
                                    write(6,73)i,j,k
     1                                    ,nint(echo_top(i,j))
     1                                    ,nint(cloud_base(i,j))
     1                                    ,nint(cloud_base_buf(i,j))
73                                  format(' Rdr Top < Bse ',2i4,i3,3i7
     1                                    ,' Potl Unresolved - CLD_RDR')

                                endif

                            endif ! Resolved Base tests

                        endif ! Below Cloud Base

                    endif ! At Echo Top

!                   Loop through range of cloud grid levels for this LAPS level
                    do kk = klow,khigh
!                       Insert radar if we are above cloud base
                        if(cld_hts(kk) .gt. cloud_base_buf(i,j)
     1                                .and. .not. l_unresolved(i,j)
     1                                                           )then
                            cldcv(i,j,kk) = 1.0
                            insert_count_lvl = insert_count_lvl + 1
                            insert_count_tot = insert_count_tot + 1
                            l_inserted = .true.
                        else ! Radar Echo below cloud base
                            l_below_base = .true.
                        endif
                    enddo ! kk

                    if(l_below_base)then
                        icount_below = icount_below + 1

                        if(icount_below .le. 100)then
                            write(6,81)i,j,k,nint(cld_hts(klow))
     1                                  ,nint(cloud_base_buf(i,j))
81                          format(' Rdr     < Bse ',2i4,i3,2i7)
                        endif
                    endif ! below base

                endif ! Reflectivity > thresh

                if(l_inserted .and. iwrite_inserted .le. 200)then
                    write(6,591)i,j,k,l_unresolved(i,j)
591                 format(' Inserted radar',2i4,i3,l2)
                    iwrite_inserted = iwrite_inserted + 1
                endif

            enddo ! j
            enddo ! i

            write(6,592)k,klow,khigh
     1          ,icount_radar_lvl,insert_count_lvl,insert_count_tot
592         format(' Inserted radar',3i3,3i8)

600         continue

        enddo ! k (LAPS grid level)

        write(6,*)' Total cloud grid points modified by radar = '
     1                                          ,insert_count_tot

        do i = 1,ni
        do j = 1,nj
            if(echo_top(i,j) .ne. r_missing_data)then
                if(echo_top_agl(i,j) .lt. echo_agl_thr .and. 
     1             .not. l_unresolved(i,j)       )then
!                   Should we set these to unresolved here or better yet above?
                    write(6,610)i,j       
     1                        ,nint(echo_top_agl(i,j))
     1                        ,nint(echo_top(i,j))
     1                        ,nint(cloud_base(i,j))
     1                        ,nint(cloud_base_buf(i,j))
 610                format('CLD_RDR - Low echo top yet resolved '
     1                    ,2i5,' ET',2i5,' Base',2i7)
                endif
            endif

            if(l_unresolved(i,j))then

!               Reconcile radar and satellite 
!               (Block out "unresolved" radar echoes)

                if(dbz_max_2d(i,j) .lt. vis_radar_thresh_dbz)then
!                   Blank out radar
                    write(6,601)i,j,dbz_max_2d(i,j)
     1                         ,nint(echo_top_agl(i,j))
601                 format(' CLD_RDR - insert_radar: '
     1                             ,'Blank out radar < '
     1                             ,2x,2i4,f6.1,' dbz',i6,' agl')

                    if(.true.)then ! Block out the radar
                        do k = 1,nk
                            grid_ra_ref(i,j,k) = ref_base
                        enddo ! k
                        dbz_max_2d(i,j) = ref_base
                    endif

                else

!                   This situation is undesirable because we want to
!                   avoid cases where the radar reflectivity is high
!                   but is "unresolved" (e.g. analyzed (SAO/IR) clouds at 
!                   this point in) the analysis. We should check as to why 
!                   it is unresolved (e.g. no clouds analyzed). For now we 
!                   blank out the radar, but we could also create a cloud at 
!                   the radar echo location.

                    write(6,602)i,j,dbz_max_2d(i,j)
     1                         ,nint(echo_top_agl(i,j))       
602                 format(' CLD_RDR - insert_radar: '
     1                             ,'Blank out radar > *'
     1                             ,1x,2i4,f6.1,' dbz',i6,' agl')

                    if(.true.)then ! Block out the radar
                        do k = 1,nk
                            grid_ra_ref(i,j,k) = ref_base
                        enddo ! k
                        dbz_max_2d(i,j) = ref_base
                    endif

                endif

            endif
        enddo ! j
        enddo ! i

999     return

        end
