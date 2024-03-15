 
        subroutine read_snd_data2(lun,i4time_snd,ext                   ! i
     1                         ,max_pr,max_pr_levels                   ! i
     1                         ,lat,lon,topo,imax,jmax,kmax            ! i
     1                         ,vgrid_3d,l_fill_ht                     ! i
     1                         ,mode                                   ! i
     1                         ,n_profiles                             ! i/o
     1                         ,elev_pr                                ! o
     1                         ,nlevels_obs_pr                         ! o
     1                         ,c5_name,obstype                        ! o
     1                         ,ob_pr_ht_obs,ob_pr_pr_obs              ! o
     1                         ,ob_pr_u_obs,ob_pr_v_obs                ! o
     1                         ,ob_pr_t_obs,ob_pr_td_obs               ! o
     1                         ,lat_pr_lvl,lon_pr_lvl,i4time_ob_pr_lvl ! o
     1                         ,istatus)                               ! o

cdoc    returns sounding wind, t, td data from the snd file
cdoc    also returns the lat/lon/time info from all the levels

!       use mem_grid, only: topo
        real topo(imax,jmax)

!       profile stuff
        real lat_pr(max_pr),lat_pr_lvl(max_pr,max_pr_levels)
        real lon_pr(max_pr),lon_pr_lvl(max_pr,max_pr_levels)
        real elev_pr(max_pr)

        real ob_pr_ht_obs(max_pr,max_pr_levels)         ! o
        real ob_pr_pr_obs(max_pr,max_pr_levels)         ! o (mb)
        real ob_pr_di_obs(max_pr_levels)                ! l
        real ob_pr_sp_obs(max_pr_levels)                ! l
        real ob_pr_u_obs(max_pr,max_pr_levels)          ! o (m/s)
        real ob_pr_v_obs(max_pr,max_pr_levels)          ! o (m/s)
        real ob_pr_t_obs(max_pr,max_pr_levels)          ! o (deg c)
        real ob_pr_td_obs(max_pr,max_pr_levels)         ! o (deg c)

        integer i4time_ob_pr(max_pr)
        integer i4time_ob_pr_lvl(max_pr,max_pr_levels)
        integer nlevels_obs_pr(max_pr)

        character*5 c5_name(max_pr)
        character*8 obstype(max_pr)
        character*9 a9time_ob, a9time_in
        character ext*(*)

        real lat_in,lon_in

        real lat(imax,jmax)
        real lon(imax,jmax)
        real vgrid_3d(imax,jmax,kmax)
        real pres_3d(imax,jmax,kmax)
        real heights_3d(imax,jmax,kmax)

        logical l_good_level, l_fill_ht, ltest_vertical_grid

        ipr_times_lvls = max_pr*max_pr_levels

        write(6,*)' subroutine read_snd_data2: ',max_pr,max_pr_levels   
     1                                          ,ipr_times_lvls

        if(ipr_times_lvls .gt. 1000000)then
            write(6,*)' warning, large allocation for ipr_times_lvls' 
     1               ,ipr_times_lvls
        endif

        istatus = 0

        rcycles_pr = 0

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)then
          write(6,*) ' error getting r_missing_data'
          go to 600
      endif

      if(ext .eq. 'snd' .or. ext .eq. 'snd')then
          call open_lapsprd_file_read(lun,i4time_snd,ext,istatus)
          if(istatus .ne. 1)then
              write(6,*) ' error opening snd file'
              go to 600
          endif
      else ! ext is treated as a full filename
          call s_len(ext,len_ext)
          istatus = 0
          open(lun,file=ext(1:len_ext),status='old',err=600)
          write(6,*)' successful open of ',ext(1:len_ext)
      endif

      if(ltest_vertical_grid('sigma_ht'))then
        print*,'topo: ',topo(3,3)
        call get_ht_3d(imax,jmax,kmax,topo,heights_3d,istatus) 
        pres_3d = vgrid_3d
      else
        heights_3d = vgrid_3d
      endif
!     print*,'hello'

      do i_pr = n_profiles+1,max_pr

        if(i_pr .le. 200 .or. i_pr .eq. (i_pr/100)*100)then
            iwrite = 1
        else
            iwrite = 0
        endif

        read(lun,511,err=530,end=550)ista,nlevels_in,
     1       lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr),
     1       c5_name(i_pr),a9time_ob,obstype(i_pr)
  511   format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)

!       obstype(i_pr) = 'raob'

        call i4time_fname_lp(a9time_ob,i4time_ob,istatus)
        i4time_ob_pr(i_pr) = i4time_ob

        if(iwrite .eq. 1)write(6,512)i_pr,ista,nlevels_in,lat_pr(i_pr)       
     1             ,lon_pr(i_pr)
     1             ,elev_pr(i_pr),rcycles_pr,c5_name(i_pr),a9time_ob       
     1             ,obstype(i_pr)     
 512    format(1x,' snd #',i4,i6,i5,2f8.2,e10.3,f8.2,1x,a5,3x,a9,1x,a8)       

        n_good_levels = 0
        nlevels_obs_pr(i_pr) = 0

        ht_prev = -9999.

        do level = 1,nlevels_in

          read(lun,*,err=515)ht_in,pr_in,t_in,td_in,di_in,sp_in ! (sp = m/s)
     1                      ,a9time_in,lat_in,lon_in

!         test this by deliberately setting ht_in to missing
!         ht_in = r_missing_data

!         determine whether we need to supply our own height (only pres given)
          if(ht_in .eq. r_missing_data .and. 
     1       pr_in .ne. r_missing_data .and. l_fill_ht)then

              call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr),lat,lon       
     1                                ,imax,jmax,ri,rj,istatus)

              if(istatus .ne. 1)goto505

              i_ob = nint(ri)
              j_ob = nint(rj)
              if(i_ob .ge. 1 .and. i_ob .le. imax .and.
     1           j_ob .ge. 1 .and. j_ob .le. jmax      )then

                  pr_in_pa = pr_in * 100.

                  if(ltest_vertical_grid('sigma_ht'))then
                      call pres_to_ht(pr_in_pa,pres_3d,heights_3d
     1                     ,imax,jmax,kmax,i_ob,j_ob,ht_buff,istatus)
                      if(istatus .ne. 1)goto505
                  else
                      call pressure_to_height(pr_in_pa,heights_3d
     1                     ,imax,jmax,kmax,i_ob,j_ob,ht_buff,istatus)
                      if(istatus .ne. 1)goto505
                  endif

                  ht_in = ht_buff
                  if(iwrite .eq. 1)
     1            write(6,*)' pressure was given, height was derived:'       
     1                     ,pr_in,ht_in
              endif
          endif

 505      if(mode .eq. 1 .or. mode .eq. 2)then   ! key good levels to wind data
              if(ht_in .ne. r_missing_data .and.
     1           di_in .ne. r_missing_data .and.
     1           sp_in .ne. r_missing_data          )then           
                  l_good_level = .true.
              else
                  l_good_level = .false.
              endif

          elseif(mode .eq. 3)then
              if(ht_in .ne. r_missing_data)then           
                  l_good_level = .true.
              else
                  l_good_level = .false.
              endif

          endif ! mode

          if(l_good_level)then           

              n_good_levels = n_good_levels + 1

              if(n_good_levels .gt. max_pr_levels)then
                  write(6,*)' error: too many sounding (.snd) levels '       
     1                     ,i_pr,n_good_levels,max_pr_levels
                  goto515
              endif

              if(ht_in .le. ht_prev)then
                  write(6,*)
     1             ' warning: sounding (.snd) ht levels out of sequence'       
     1                     ,i_pr,n_good_levels,ht_prev,ht_in
                  write(6,*)' keep just the lower part of this sounding'
                  goto 516
              endif

              ht_prev = ht_in

              nlevels_obs_pr(i_pr) = n_good_levels
              ob_pr_ht_obs(i_pr,n_good_levels) = ht_in
              ob_pr_pr_obs(i_pr,n_good_levels) = pr_in
              ob_pr_di_obs(n_good_levels) = di_in
              ob_pr_sp_obs(n_good_levels) = sp_in

              ! yuanfu adds an 'if' for wind in mode = 3 case:
	      if ((ob_pr_di_obs(n_good_levels) .eq. r_missing_data) 
     1	        .or. (ob_pr_sp_obs(n_good_levels) .eq. r_missing_data)) 
     1          then
                ob_pr_u_obs(i_pr,n_good_levels) = r_missing_data
                ob_pr_v_obs(i_pr,n_good_levels) = r_missing_data
              else
                call disp_to_uv(ob_pr_di_obs(n_good_levels),
     1                          ob_pr_sp_obs(n_good_levels),
     1                          ob_pr_u_obs(i_pr,n_good_levels),
     1                          ob_pr_v_obs(i_pr,n_good_levels))
              endif

!             write(6,311,err=312)ista,i_pr
!    1        ,ob_snd_ht_obs(i_pr,level)
!    1        ,ob_snd_t_obs(i_pr,level)
!311          format(1x,i6,i4,5f8.1)

              ob_pr_t_obs(i_pr,n_good_levels)  = t_in
              ob_pr_td_obs(i_pr,n_good_levels) = td_in

              lat_pr_lvl(i_pr,n_good_levels) = lat_in
              lon_pr_lvl(i_pr,n_good_levels) = lon_in

              call i4time_fname_lp(a9time_in,i4time_ob,istatus)
              i4time_ob_pr_lvl(i_pr,n_good_levels) = i4time_ob

          endif ! good data at the level

          go to 516

  515     write(6,*)' error reading raob, raw level # =',level
          write(6,*)' while reading sounding # ',(i_pr-n_profiles)
          n_profiles=i_pr-1
          istatus = 0
          go to 600

  516   end do                   ! level
      end do ! i_pr

      write(6,*)' error: # of soundings+profilers has reached'
     1         ,' dimensions of max_pr ',max_pr

      n_profiles=max_pr
      istatus = 0
      go to 600
c
c     sounding reading error handling
c
  530 write(6,*) ' error during read of snd file (sounding header)'       
      write(6,*) ' while reading sounding number ',(i_pr-n_profiles)
      n_profiles=i_pr-1
      istatus = 0
      go to 600

  550 continue ! used for end of file
      n_profiles=i_pr-1
      istatus = 1

  600 continue 

      close(lun)

      return
      end



        subroutine read_snd_data(lun,i4time_snd,ext                    ! i
     1                         ,max_pr,max_pr_levels                   ! i
     1                         ,lat,lon,topo,imax,jmax,kmax            ! i
     1                         ,heights_3d,l_fill_ht                   ! i
     1                         ,mode                                   ! i
     1                         ,n_profiles                             ! i/o
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr   ! o
     1                         ,c5_name,i4time_ob_pr,obstype           ! o
     1                         ,ob_pr_ht_obs,ob_pr_pr_obs              ! o
     1                         ,ob_pr_u_obs,ob_pr_v_obs                ! o
     1                         ,ob_pr_t_obs,ob_pr_td_obs               ! o
     1                         ,istatus)                               ! o

cdoc    returns sounding wind, t, td data from the snd file

        real topo(imax,jmax)

!       profile stuff
        real lat_pr(max_pr),lat_pr_lvl(max_pr,max_pr_levels)
        real lon_pr(max_pr),lon_pr_lvl(max_pr,max_pr_levels)
        real elev_pr(max_pr)

        real ob_pr_ht_obs(max_pr,max_pr_levels)         ! o
        real ob_pr_pr_obs(max_pr,max_pr_levels)         ! o (mb)
        real ob_pr_di_obs(max_pr_levels)                ! l
        real ob_pr_sp_obs(max_pr_levels)                ! l
        real ob_pr_u_obs(max_pr,max_pr_levels)          ! o (m/s)
        real ob_pr_v_obs(max_pr,max_pr_levels)          ! o (m/s)
        real ob_pr_t_obs(max_pr,max_pr_levels)          ! o (deg c)
        real ob_pr_td_obs(max_pr,max_pr_levels)         ! o (deg c)

        integer i4time_ob_pr(max_pr)
        integer i4time_ob_pr_lvl(max_pr,max_pr_levels)
        integer nlevels_obs_pr(max_pr)

        character*5 c5_name(max_pr)
        character*8 obstype(max_pr)
        character*9 a9time_ob
        character ext*(*)

        real lat(imax,jmax)
        real lon(imax,jmax)
        real heights_3d(imax,jmax,kmax)

        logical l_good_level, l_fill_ht

        n_profiles_in = n_profiles

        call read_snd_data2(    lun,i4time_snd,ext                     ! i
     1                         ,max_pr,max_pr_levels                   ! i
     1                         ,lat,lon,topo,imax,jmax,kmax            ! i
     1                         ,heights_3d,l_fill_ht                   ! i
     1                         ,mode                                   ! i
     1                         ,n_profiles                             ! i/o
     1                         ,elev_pr                                ! o
     1                         ,nlevels_obs_pr                         ! o
     1                         ,c5_name,obstype                        ! o
     1                         ,ob_pr_ht_obs,ob_pr_pr_obs              ! o
     1                         ,ob_pr_u_obs,ob_pr_v_obs                ! o
     1                         ,ob_pr_t_obs,ob_pr_td_obs               ! o
     1                         ,lat_pr_lvl,lon_pr_lvl,i4time_ob_pr_lvl ! o
     1                         ,istatus)                               ! o

        if(n_profiles_in .lt. n_profiles)then
            do i_pr = n_profiles_in+1,n_profiles
                lat_pr(i_pr) = lat_pr_lvl(i_pr,1)
                lon_pr(i_pr) = lon_pr_lvl(i_pr,1)
                i4time_ob_pr(i_pr) = i4time_ob_pr_lvl(i_pr,1)
            enddo ! i_pr
        endif

        return
        end



        subroutine read_snd_metadata(lun,i4time_snd,ext                ! i
     1                         ,max_pr,max_pr_levels                   ! i
     1                         ,lat,lon,imax,jmax                      ! i
     1                         ,n_profiles                             ! o
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr   ! o
     1                         ,c5_name,i4time_ob_pr,obstype           ! o
     1                         ,cloud_base_temp,cloud_liquid_int       ! o
     1                         ,istatus)                               ! o

cdoc    returns sounding metadata from the snd file

!       profile stuff
        real lat_pr(max_pr)
        real lon_pr(max_pr)
        real elev_pr(max_pr)

        real cloud_base_temp(max_pr)
        real cloud_liquid_int(max_pr)

        integer i4time_ob_pr(max_pr)
        integer nlevels_obs_pr(max_pr)

        character*5 c5_name(max_pr)
        character*8 obstype(max_pr)
        character*9 a9time_ob
        character ext*(*)
        character*200 line

        real lat(imax,jmax)
        real lon(imax,jmax)

        istatus = 0

        rcycles_pr = 0

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)then
          write(6,*) ' error getting r_missing_data'
          go to 600
      endif

      call open_lapsprd_file_read(lun,i4time_snd,ext,istatus)
      if(istatus .ne. 1)then
          write(6,*) ' error opening snd file'
          istatus = 0
          go to 600
      endif

      do i_pr = 1,max_pr

        if(i_pr .le. 200 .or. i_pr .eq. (i_pr/10)*10)then
            iwrite = 1
        else
            iwrite = 0
        endif

        read(lun,501,err=530,end=550)line
 501    format(a)
 
        read(line,511,err=530,end=550)ista,nlevels_in,
     1           lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr),
     1           c5_name(i_pr),a9time_ob,obstype(i_pr)
 
        if(obstype(i_pr) .eq. 'radiomtr')then
            read(line,511,err=530)ista,nlevels_in,
     1           lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr),
     1           c5_name(i_pr),a9time_ob,obstype(i_pr),
     1           cloud_base_temp(i_pr),cloud_liquid_int(i_pr)
        else
            read(line,511,err=530)ista,nlevels_in,
     1           lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr),
     1           c5_name(i_pr),a9time_ob,obstype(i_pr)
            cloud_base_temp(i_pr) = r_missing_data     
            cloud_liquid_int(i_pr) = r_missing_data     
        endif

  511   format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8,f8.0,f8.0)

        call i4time_fname_lp(a9time_ob,i4time_ob,istatus)
        i4time_ob_pr(i_pr) = i4time_ob

        if(iwrite .eq. 1)write(6,512)i_pr,ista,nlevels_in,lat_pr(i_pr)       
     1             ,lon_pr(i_pr)
     1             ,elev_pr(i_pr),rcycles_pr,c5_name(i_pr),a9time_ob
     1             ,obstype(i_pr)     
 512    format(/' snd #',i4,i6,i5,2f8.2,e10.3,f8.2,1x,a5,3x,a9,1x,a8)       

        n_good_levels = 0

        do level = 1,nlevels_in

          read(lun,*,end=515)

!         test this by deliberately setting ht_in to missing
!         ht_in = r_missing_data


          go to 516

  515     write(6,*)' error reading level # =',level
          write(6,*)' while reading sounding # ',(i_pr-n_profiles)
          n_profiles=i_pr-1
          istatus = 0
          go to 600

  516   end do                   ! level
      end do ! i_pr

      write(6,*)' error: # of soundings has reached'
     1         ,' dimensions of max_pr ',max_pr

      n_profiles=max_pr
      istatus = 0
      go to 600
c
c     sounding reading error handling
c
  530 write(6,*) ' error during read of snd file (sounding header)'
      write(6,*) ' while reading sounding number ',(i_pr-n_profiles)
      n_profiles=i_pr-1
      istatus = 0
      go to 600

  550 continue ! used for end of file
      n_profiles=i_pr-1

      istatus = 1

  600 continue 

      close(lun)

      return
      end

	subroutine read_sfc_snd(i4time,atime_s,
     &       n_obs_g,n_obs_b,                                            ! i/o
     &       obstime,wmoid,stations,provider,wx_s,reptype,autostntype,
     &       lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &       alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &       pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &       td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &       sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,maxsta,       
     &       lat,lon,imax,jmax,kmax,                                     ! i
     &       max_pr,max_pr_levels,                                       ! i
     &       topo,                                                       ! i
     &       istatus)

        include 'constants.inc'

        include 'read_sfc.inc'

        real topo(imax,jmax)

!       declarations for 'read_snd_data' call
        integer max_pr,max_pr_levels

        real lat_pr(max_pr)
        real lon_pr(max_pr)
        real elev_pr(max_pr)

        real ob_pr_ht_obs(max_pr,max_pr_levels)                        ! o
        real ob_pr_pr_obs(max_pr,max_pr_levels)                        ! o
        real ob_pr_u_obs(max_pr,max_pr_levels)                         ! o
        real ob_pr_v_obs(max_pr,max_pr_levels)                         ! o
        real ob_pr_t_obs(max_pr,max_pr_levels)                         ! o
        real ob_pr_td_obs(max_pr,max_pr_levels)                        ! o

        real lat(imax,jmax)
        real lon(imax,jmax)

        real heights_3d_dum(imax,jmax,kmax)

        integer i4time_ob_pr(max_pr)
        integer nlevels_obs_pr(max_pr)

        character*5 c5_name(max_pr)
        character*8 obstype(max_pr)
        character*9 a9time_ob

!       other declarations
        logical l_good_snd
        character*3 ext

        write(6,*)' subroutine read_sfc_snd...'

        call get_r_missing_data(r_missing_data,istatus)
        call get_sfc_badflag(badflag,istatus)
        elev_msg = -999.

        n_profiles = 0
        n_good_snd = 0
        ext = 'snd'
        lun = 44
        mode = 3

        call read_snd_data(lun,i4time,ext                              ! i
     1                         ,max_pr,max_pr_levels                   ! i
     1                         ,lat,lon,topo,imax,jmax,kmax            ! i
     1                         ,heights_3d_dum,.false.                 ! i
     1                         ,mode                                   ! i
     1                         ,n_profiles                             ! i/o
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr   ! o
     1                         ,c5_name,i4time_ob_pr,obstype           ! o
     1                         ,ob_pr_ht_obs,ob_pr_pr_obs              ! o
     1                         ,ob_pr_u_obs,ob_pr_v_obs                ! o
     1                         ,ob_pr_t_obs,ob_pr_td_obs               ! o
     1                         ,istatus)                               ! o
        if(istatus .ne. 1)then
            write(6,*)' bad istatus from read_snd_data'
            go to 999
        endif

        if(n_profiles .eq. 0)then
            write(6,*)' no profiles available via read_snd_data'
            go to 999
        endif

!       append sounding obs to surface arrays
        do i_pr = 1,n_profiles
            l_good_snd = .false.

            ilevel_best_wind      = 0
            ilevel_best_temp      = 0
            ilevel_best_mslp      = 0
            ilevel_best_stnp      = 0

            if(nlevels_obs_pr(i_pr) .gt. 0    .and.
     1         obstype(i_pr)(1:4) .ne. 'goes' .and.
     1         obstype(i_pr)(1:4) .ne. 'poes'       )then 

              if(elev_pr(i_pr) .ne. elev_msg)then ! assume station elev present
                height_best_wind_diff = 99999.
                height_best_wind      = r_missing_data
                ob_best_wind          = r_missing_data

                height_best_temp_diff = 99999.
                height_best_temp      = r_missing_data
                ob_best_temp          = r_missing_data
                ob_best_dwpt          = r_missing_data

                ob_best_mslp          = r_missing_data

                ob_best_stnp          = r_missing_data

!               find closest temp and wind obs to levels to 2m and 10m above
!               station elevation, respectivtly
                do il = 1,nlevels_obs_pr(i_pr)
                    rlevel_temp_diff = abs( 
     1                   (ob_pr_ht_obs(i_pr,il)-elev_pr(i_pr)) - 2.0 )

                    rlevel_wind_diff = abs( 
     1                   (ob_pr_ht_obs(i_pr,il)-elev_pr(i_pr)) - 10.0 )

                    if(rlevel_temp_diff .lt. height_best_temp_diff)then
                      if(ob_pr_t_obs(i_pr,il) .ne. r_missing_data)then
                        height_best_temp_diff = rlevel_temp_diff
                        ilevel_best_temp      = il
                        height_best_temp      = ob_pr_ht_obs(i_pr,il)
                      endif
                    endif

                    if(rlevel_wind_diff .lt. height_best_wind_diff)then       
                      if(ob_pr_u_obs(i_pr,il) .ne. r_missing_data
     1             .and. ob_pr_v_obs(i_pr,il) .ne. r_missing_data)then
                        height_best_wind_diff = rlevel_wind_diff
                        height_best_wind      = ob_pr_ht_obs(i_pr,il)
                        ilevel_best_wind      = il
                      endif
                    endif

                enddo ! il

!               test whether best diffs represent an overall good level for 
!               temp/wind measurement
                if(height_best_wind_diff .le. 6.0 .or. 
     1             height_best_temp_diff .le. 2.0)then
                    write(6,*)' good sounding elev     ',i_pr
     1                       ,c5_name(i_pr)
     1                       ,height_best_temp,height_best_temp_diff       
     1                       ,height_best_wind,height_best_wind_diff
                    l_good_snd = .true.
                    n_good_snd = n_good_snd + 1
                elseif(i_pr .le. 100)then
                    write(6,*)' unusable sounding elev ',i_pr
     1                       ,c5_name(i_pr)
     1                       ,height_best_temp,height_best_temp_diff       
     1                       ,height_best_wind,height_best_wind_diff
                endif

              else ! station elevation is missing (e.g. for dropsondes)

!               these steps assume station elevation is -999. (missing value)
                il = 1

!               step 2, fill in mslp if lowest level is near zero
                if(abs(ob_pr_ht_obs(i_pr,il)) .le. 1.)then
                    l_good_snd = .true.
                    n_good_snd = n_good_snd + 1
                    ilevel_best_mslp = il
                    ob_best_mslp = ob_pr_pr_obs(i_pr,ilevel_best_mslp)
                    write(6,*)' good sounding for mslp ',i_pr
     1                       ,c5_name(i_pr)
     1                       ,ilevel_best_mslp
     1                       ,ob_best_mslp
                endif

!               step 3, fill in station pressure if lowest level is near topo,
!                       topo is being passed in for this

                call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr)
     1                                  ,lat,lon       
     1                                  ,imax,jmax,ri,rj,istatus)
                if(istatus .ne. 1)then
                    write(6,*)
     1                     ' warning, laps grid location not calculated'
                else
                    call bilinear_laps(ri,rj,imax,jmax,topo,topo_sta)
                    if(abs(ob_pr_ht_obs(i_pr,il) - topo_sta) .le. 100.
     1                                                            )then                    
                        l_good_snd = .true.
                        n_good_snd = n_good_snd + 1
                        ilevel_best_stnp = il
                        ob_best_stnp = 
     1                      ob_pr_pr_obs(i_pr,ilevel_best_stnp)     
                        write(6,*)' good sounding for stnp ',i_pr
     1                       ,c5_name(i_pr)
     1                       ,ilevel_best_stnp,topo_sta
     1                       ,ob_best_stnp
                    endif
                endif

              endif ! station elevation present

            endif ! nlevels > 0 as well as valid obstype

!           if(.false.)then
            if(l_good_snd)then ! this ob is good enough to append
                n_obs_g = n_obs_g + 1
                n_obs_b = n_obs_b + 1

                lat_s(n_obs_b) = lat_pr(i_pr)
                lon_s(n_obs_b) = lon_pr(i_pr)
                call sfci4_to_sfchhmm(i4time_ob_pr(i_pr)
     1                               ,obstime(n_obs_b),istatus) ! assign obstime
                elev_s(n_obs_b) = elev_pr(i_pr)
                stations(n_obs_b) = c5_name(i_pr)
                reptype(n_obs_b) = obstype(i_pr)

!               initialize values to badflag
                t_s(n_obs_b) = badflag
                td_s(n_obs_b) = badflag
                rh_s(n_obs_b) = badflag
                dd_s(n_obs_b) = badflag
                ff_s(n_obs_b) = badflag
                ddg_s(n_obs_b) = badflag
                ffg_s(n_obs_b) = badflag
                alt_s(n_obs_b) = badflag
                pstn_s(n_obs_b) = badflag
                pmsl_s(n_obs_b) = badflag
                delpch(n_obs_b) = badflag
                delp(n_obs_b) = badflag
                vis_s(n_obs_b) = badflag
                solar_s(n_obs_b) = badflag
                sfct(n_obs_b) = badflag
                sfcm(n_obs_b) = badflag
                pcp1(n_obs_b) = badflag
                pcp3(n_obs_b) = badflag
                pcp6(n_obs_b) = badflag
                pcp24(n_obs_b) = badflag
                snow(n_obs_b) = badflag
                kloud_s(n_obs_b) = badflag
                max24t(n_obs_b) = badflag
                min24t(n_obs_b) = badflag

                t_ea(n_obs_b) = badflag
                td_ea(n_obs_b) = badflag
                rh_ea(n_obs_b) = badflag
                dd_ea(n_obs_b) = badflag
                ff_ea(n_obs_b) = badflag
                alt_ea(n_obs_b) = badflag
                p_ea(n_obs_b) = badflag
                vis_ea(n_obs_b) = badflag
                solar_ea(n_obs_b) = badflag
                sfct_ea(n_obs_b) = badflag
                sfcm_ea(n_obs_b) = badflag
                pcp_ea(n_obs_b)  = badflag
                snow_ea(n_obs_b) = badflag
                store_amt(n_obs_b,:) = '     '
                store_hgt(n_obs_b,:) = badflag

!               add temperature ob into arrays 
                if(ilevel_best_temp .gt. 0)then
                    if(ob_pr_t_obs(i_pr,ilevel_best_temp) 
     1                                         .ne. r_missing_data)then
                        t_ob_c = ob_pr_t_obs(i_pr,ilevel_best_temp)     
                        t_s(n_obs_b) = c_to_f(t_ob_c)
                        write(6,*)' good temp ob',t_s(n_obs_b)
                    endif      

!                   add dewpoint ob into arrays 
                    if(ob_pr_td_obs(i_pr,ilevel_best_temp) 
     1                                         .ne. r_missing_data)then       
                        td_ob_c = ob_pr_td_obs(i_pr,ilevel_best_temp)     
                        td_s(n_obs_b) = c_to_f(td_ob_c)
                        write(6,*)' good dwpt ob',td_s(n_obs_b)
                    endif
                endif ! good temp level determined

!               add wind ob into arrays
                if(ilevel_best_wind .gt. 0)then
                  if(ob_pr_u_obs(i_pr,ilevel_best_wind)
     1                                          .ne.r_missing_data      
     1       .and.   ob_pr_v_obs(i_pr,ilevel_best_wind)
     1                                          .ne.r_missing_data
     1                                                         )then
                    call uv_to_disp(ob_pr_u_obs(i_pr,ilevel_best_wind)
     1                             ,ob_pr_v_obs(i_pr,ilevel_best_wind)
     1                             ,dd_s(n_obs_b)
     1                             ,speed_ms )
                    ff_s(n_obs_b) = speed_ms / mspkt    
                    write(6,*)' good wind ob',dd_s(n_obs_b)
     1                                       ,ff_s(n_obs_b)      
                  endif

                endif ! good wind level determined

!               add mslp ob into arrays
                if(ilevel_best_mslp .gt. 0)then
                  if(ob_pr_pr_obs(i_pr,ilevel_best_mslp)
     1                                          .ne.r_missing_data      
     1                                                         )then
                    pmsl_s(n_obs_b) = 
     1                  ob_pr_pr_obs(i_pr,ilevel_best_mslp)
                    write(6,*)' good mslp ob',pmsl_s(n_obs_b)
                  endif

                endif ! good mslp level determined

!               add stnp ob into arrays
                if(ilevel_best_stnp .gt. 0)then
                  if(ob_pr_pr_obs(i_pr,ilevel_best_stnp)
     1                                          .ne.r_missing_data      
     1                                                         )then
                    pstn_s(n_obs_b) = 
     1                  ob_pr_pr_obs(i_pr,ilevel_best_stnp)
                    elev_s(n_obs_b) = topo_sta
                    write(6,*)' good stnp ob',pstn_s(n_obs_b)
                  endif

                endif ! good stnp level determined

            endif ! good sounding

        enddo ! i_pr

 999	write(6,*)' # of good soundings = ',n_good_snd,n_obs_g,n_obs_b       

        return
        end

