cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis 
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps 
cdis 
cdis    this software and its documentation are in the public domain and 
cdis    are furnished "as is."  the united states government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  they assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  if significant modifications or enhancements 
cdis    are made to this software, the fsl software policy manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
        subroutine rd_sfc_t(i4time,heights_3d,temp_bkg_3d           ! i
     1                       ,pres_3d                               ! i
     1                       ,max_sfc                               ! i
     1                       ,n_good_sfc                            ! o
     1                       ,ext_in                                ! i
!    1                       ,u_maps_inc,v_maps_inc                 ! i
     1                       ,ni,nj,nk                              ! i
     1                       ,lat,lon,r_missing_data                ! i
     1                       ,temp_obs,max_obs,n_obs                ! i/o
     1                       ,istatus)                              ! o

!       2001        steve albers, fsl      can be called for sfc temps

!******************************************************************************

        use mem_namelist, only: iwrite_output

        include 'tempobs.inc'

!       laps grid dimensions
        real lat(ni,nj)
        real lon(ni,nj)

!       sfc stations

        integer sfc_i(max_sfc)     ! i sfc gridpoint
        integer sfc_j(max_sfc)     ! j sfc gridpoint
        integer sfc_ht(max_sfc)    ! ht sfc
        real    sfc_temp(max_sfc)  ! sfc temp

        real sfc_lat(max_sfc)
        real sfc_lon(max_sfc)
        real sfc_elev(max_sfc)

c
        integer i4time, jstatus
c
        character filetime*9, infile*256, btime*24
        character stations(max_sfc)*20, provider(max_sfc)*11


!******************************************************************************

        real heights_3d(ni,nj,nk)
        real temp_bkg_3d(ni,nj,nk)
        real pres_3d(ni,nj,nk)

        real u_maps_inc(ni,nj,nk)
        real v_maps_inc(ni,nj,nk)

        character*9 asc9_tim_sfc
        character ext*31, ext_in*3

        logical l_eof

        write(6,*)
        write(6,*)' subroutine rd_sfc_t...'

        n_sfc_read = 0
        n_sfc_obs = 0
        n_good_sfc = 0
        n_bad_sfc = 0

        lun_tmg = 32
        ext = 'tmg'
        close(lun_tmg)

!       open output intermediate graphics file
        if(iwrite_output .ge. 1)then
            call open_lapsprd_file_append(lun_tmg,i4time,ext,istatus)       
            if(istatus .ne. 1)go to 888
        endif

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' error getting laps_cycle_time'
            return
        endif

        call get_tempob_time_window('sfc',i4_window_sfc,istatus)
        if(istatus .eq. 1)then
            write(6,*)' i4_window_sfc = ',i4_window_sfc
        else
            write(6,*)' error getting i4_window_sfc'
            return
        endif

        write(6,*)
        write(6,*)'             reading sfc obs: ',ext_in
        write(6,*)
     1  '   n   i  j  k   t-raw  t-grid  bias  '

10      i_qc = 1

        if(n_sfc_read .le. 200)then
            iwrite = 1
            write(6,*)
        else
            iwrite = 0
        endif

!       call read_sfc_ob(lun_in,'temp',xlat,xlon,elev_in,temp_ob,arg2       
!    1                                  ,asc9_tim_sfc,iwrite,l_eof)

	call read_sfc_temp(i4time,btime,n_obs_g,n_obs_b 
     1                    ,stations,provider
     1                    ,sfc_lat,sfc_lon,sfc_elev
     1                    ,sfc_temp,max_sfc
     1                    ,jstatus)

        do i = 1,n_obs_b

          temp_ob = sfc_temp(i)

          elev_in = sfc_elev(i)

          if(elev_in .eq. 0.)i_qc = 0

          n_sfc_read = n_sfc_read + 1

          call cv_asc_i4time(asc9_tim_sfc,i4time_sfc)

          if(abs(i4time_sfc - i4time) .le. i4_window_sfc)then ! in time window

            rcycles = float(i4time - i4time_sfc) 
     1              / float(ilaps_cycle_time)

!           climo qc check
            if(temp_ob .lt. 500. .and. i_qc .eq. 1)then

                call latlon_to_rlapsgrid(sfc_lat(i),sfc_lon(i)
     1                                  ,lat,lon,ni,nj
     1                                  ,ri,rj,istatus)
                i_grid = nint(ri)
                j_grid = nint(rj)

                if(i_grid .ge.  1 .and. j_grid .ge. 1 .and.
     1             i_grid .le. ni .and. j_grid .le. nj)then

!                   sfc is in horizontal domain

                    if(ext_in .eq. 'pin')then
!                       assume sfc elev is pressure altitude msl
                        elev_std = elev_in

                        if(abs(elev_std) .lt. 90000.)then ! within flag value
                            pres_mb = ztopsa(elev_std)
                            pres_pa = pres_mb * 100.
                            call pressure_to_height(pres_pa,heights_3d       
     1                                             ,ni,nj,nk
     1                                             ,i_grid,j_grid
     1                                             ,elev_geo
     1                                             ,istatus_rk)      
                        else
                            istatus_rk = 0
                        endif

                        if(istatus_rk .eq. 1)then
                            rk = height_to_zcoord2(elev_geo,heights_3d       
     1                           ,ni,nj,nk,i_grid,j_grid,istatus_rk)       
                        endif

                        if(istatus_rk .ne. 1)then
                            write(6,*)' warning: rejecting sfc ',
     1                      'elevation questionable ',elev_std,elev_geo       
                        endif

                    else 
                        istatus = 0
                        return

                    endif

                    k_grid = nint(rk)

                    if(istatus_rk .eq. 1
     1             .and. k_grid .le. nk
     1             .and. k_grid .ge. 1    )then ! sfc is in vertical domain

                        n_sfc_obs = n_sfc_obs + 1

                        if(n_sfc_obs .gt. max_sfc)then
                           write(6,*)' warning: too many sfcs, '
     1                              ,'limit is ',max_sfc
                           istatus = 0
                           return
                        endif

                        sfc_i(n_sfc_obs) = i_grid
                        sfc_j(n_sfc_obs) = j_grid

                        sfc_ht(n_sfc_obs) = elev_geo

                        t_diff = 0.
!                       call get_time_term(t_diff)

                        pres_ob = r_missing_data

                        call interp_tobs_to_laps(
     1                             elev_geo,temp_ob,                     ! i
     1                             pres_ob,                              ! i
     1                             t_diff,temp_bkg_3d,                   ! i
     1                             t_interp,                             ! o
     1                             1,iwrite,k_grid,.true.,               ! i
     1                             1,                                    ! i
     1                             lat_pr,lon_pr,i_grid,j_grid,          ! i
     1                             ni,nj,nk,                             ! i
     1                             1,1,r_missing_data,                   ! i
     1                             pres_3d,                              ! i
     1                             heights_3d)                           ! i

                        if(t_interp .eq. r_missing_data)then
                            write(6,*)' error: t_interp = ',t_interp
                            istatus = 0
                            return
                        endif

                        if(iwrite_output .ge. 1)then
                            write(lun_tmg,*)ri,rj,k_grid
     1                                     ,t_interp,'sfc'
                        endif

!                       calculate observation bias
                        bias = t_interp - 
     1                         temp_bkg_3d(i_grid,j_grid,k_grid)

!                       qc check of bias
                        if(abs(bias) .le. 10.)then
                            n_good_sfc = n_good_sfc + 1            
                            n_obs = n_obs + 1

                            if(n_obs .gt. max_obs)then
                                write(6,*)
     1                        ' error - too many obs in data structure'       
                                write(6,*)
     1                        ' increase max_obs parameter from',max_obs     
                                istatus = 0
                                return
                            endif

!                           insert ob into data structure
                            temp_obs(n_obs,i_ri) = i_grid
                            temp_obs(n_obs,i_rj) = j_grid
                            temp_obs(n_obs,i_rk) = rk
                            temp_obs(n_obs,i_ob_raw) = temp_ob
                            temp_obs(n_obs,i_i) = i_grid
                            temp_obs(n_obs,i_j) = j_grid
                            temp_obs(n_obs,i_k) = k_grid
                            temp_obs(n_obs,i_ob_grid) = t_interp
                            temp_obs(n_obs,i_wt) = 1.0
                            temp_obs(n_obs,i_bias) = bias
                            temp_obs(n_obs,i_inst_err) = 1.0
 

                        else
                            n_bad_sfc = n_bad_sfc + 1            
          
                        endif


                        if(iwrite .eq. 1)write(6,20,err=21)
     1                                    n_sfc_read,n_sfc_obs       
     1                                   ,i_grid,j_grid,k_grid       
     1                                   ,rk ! ,rk_pspace
     1                                   ,temp_ob,t_interp,bias
20                      format(2i5,1x,3i4,2x,f8.3,2x,3f7.1)
21                      continue

                    else
                        write(6,*)' note: out of vertical bounds'
     1                             ,n_sfc_read       

                    endif ! in vertical bounds

                else
                    write(6,*)' out of horizontal bounds'
     1                        ,n_sfc_read,i_grid,j_grid        

                endif ! in horizontal bounds
            endif ! good data

          else
            write(6,*)' out of temporal bounds',n_sfc_read
     1                              ,abs(i4time_sfc - i4time)

          endif ! in temporal bounds

900     enddo ! i

        write(6,*)' end of sfc ',ext_in,' file'

        write(6,*)' # of sfc read in = ',n_sfc_read
        write(6,*)' # of sfc passing bounds checks = ',n_sfc_obs      
        write(6,*)' # of sfc passing qc check = ',n_good_sfc
        write(6,*)' # of sfc failing qc check = ',n_bad_sfc
        write(6,*)' % of sfc failing qc check = ',
     1                     pct_rejected(n_good_sfc,n_bad_sfc)

        close(lun_tmg)

        istatus = 1

        return

999     write(6,*)' no sfc data present'
        istatus = 1
        return


888     write(6,*)' open error for tmg file'
        istatus = 0
        close(lun_tmg)
        return

        end


