
      subroutine ingest_drpsnd(path_to_raw_drpsnd,c8_drpsnd_format
     1                        ,l_fill_ht,lun_out)                  ! i       

!     steve albers fsl   2003     original version

!     input file 
      character*200 filename_in
!     character*200 dropsonde_in
      character*9 a9_time
      character*180 dir_in
      character*255 c_filespec
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)
      character*9 a8_to_a9
      character*8 a8_time,a8_time_orig(max_files)
      character*8 c8_drpsnd_format

!     output file
      character*13 filename13, cvt_i4time_wfo_fname13
      character*31    ext
      integer       len_dir

      character*40 c_vars_req
      character*(*) path_to_raw_drpsnd

      logical l_parse, l_fill_ht

!     define interval to be used (between timestamps) for creation of snd files
      integer i4_snd_interval
      parameter (i4_snd_interval = 3600)

      iopen = 0

      call getenv('laps_a9time',a9_time)
      call s_len(a9_time,ilen)

      if(ilen .eq. 9)then
          write(6,*)' systime (from env) = ',a9_time
          call i4time_fname_lp(a9_time,i4time_sys,istatus)
      else
          call get_systime(i4time_sys,a9_time,istatus)
          if(istatus .ne. 1)go to 999
          write(6,*)' systime = ',a9_time
      endif


!     i4time_sys = (i4time_sys/i4_snd_interval) * i4_snd_interval ! for testing only

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

      call get_laps_cycle_time(ilaps_cycle_time,istatus)
      if(istatus .eq. 1)then
          write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
      else
          write(6,*)' error getting laps_cycle_time'
          go to 999
      endif

      dir_in = path_to_raw_drpsnd

      call s_len(dir_in,len_dir_in)

      if(c8_drpsnd_format(1:6) .eq. 'nimbus')then
          c_filespec = dir_in(1:len_dir_in)//'/*0300o'
          call get_file_times(c_filespec,max_files,c_fnames
     1                       ,i4times,i_nbr_files_ret,istatus)

!tt -- copy and paste of the previous. i'll change the file names of the
!tt    avaps dropsondes to follow the nimbus/raob convention.

      elseif(c8_drpsnd_format(1:5) .eq. 'avaps')then
          c_filespec = dir_in(1:len_dir_in)//'/*0300o'
          call get_file_times(c_filespec,max_files,c_fnames
     1                       ,i4times,i_nbr_files_ret,istatus)

      elseif(c8_drpsnd_format(1:3) .eq. 'wfo' .or.
     1       c8_drpsnd_format(1:3) .eq. 'rsa' .or.
     1       c8_drpsnd_format(1:3) .eq. 'snd'      )then
          c_filespec = dir_in(1:len_dir_in)
          call get_file_times(c_filespec,max_files,c_fnames
     1                       ,i4times,i_nbr_files_ret,istatus)

      elseif(c8_drpsnd_format(1:3) .eq. 'cwb')then
          c_filespec = dir_in(1:len_dir_in)//'/drpsnd*'
          call get_file_names(c_filespec,i_nbr_files_ret,c_fnames
     1                      ,max_files,istatus)

!         obtain file times from file names
          do i = 1,i_nbr_files_ret
              call s_len(c_fnames(i),len_fname)
              call get_directory_length(c_fnames(i),len_dir)
              a8_time = c_fnames(i)(len_dir+7:len_fname)
              a8_time_orig(i) = a8_time

              a9_time = a8_to_a9(a8_time)
              call i4time_fname_lp(a9_time,i4times(i),istatus)
              write(6,*)c_fnames(i)(1:len_fname),i4times(i)
          enddo ! i

      else
          write(6,*)' error - invalid c8_drpsnd_format '
     1             ,c8_drpsnd_format      
          istatus = 0
          goto 999

      endif

!     get dropsonde time window
      call get_windob_time_window('raob',i4_wind_ob,istatus)
      if(istatus .ne. 1)goto 997

      call get_tempob_time_window('raob',i4_temp_ob,istatus)
      if(istatus .ne. 1)goto 997

      i4_drpsnd_window = max(i4_wind_ob,i4_temp_ob)
      isnd_staname = 0

!     loop through dropsonde files and choose ones in time window
      write(6,*)' # of files using filename format ',c8_drpsnd_format      
     1                                        ,' = ',i_nbr_files_ret
      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)

          if(c8_drpsnd_format(1:6) .eq. 'nimbus')then
              filename_in = dir_in(1:len_dir_in)//'/'//a9_time//'0300o'       
!             i4_drpsnd_window = 60000  ! temporary for testing
              i4_contains_early = 7200
              i4_contains_late  = 3600

          elseif(c8_drpsnd_format(1:5) .eq. 'avaps')then
              filename_in = dir_in(1:len_dir_in)//'/'//a9_time//'0300o'       
              i4_contains_early = 0
              i4_contains_late  = 10800

          elseif(c8_drpsnd_format(1:5) .eq. 'snd')then
              filename_in = dir_in(1:len_dir_in)//'/'//a9_time//'.snd'       
              i4_contains_early = 3600
              i4_contains_late  = 3600

          elseif(c8_drpsnd_format(1:3) .eq. 'wfo')then
              filename13 = cvt_i4time_wfo_fname13(i4times(i))
              filename_in = dir_in(1:len_dir_in)//'/'//filename13      
              i4_contains_early = 43200
              i4_contains_late  = 0

          elseif(c8_drpsnd_format(1:3) .eq. 'rsa')then
              filename13 = cvt_i4time_wfo_fname13(i4times(i))
              filename_in = dir_in(1:len_dir_in)//'/'//filename13      
              i4_contains_early = 43200
              i4_contains_late  = 0

          elseif(c8_drpsnd_format(1:3) .eq. 'cwb')then 
!             filename_in = dir_in(1:len_dir_in)//'/temp'//
!    1                      a8_time_orig(i)//'.dat'
!             this may need to be adjusted
              filename_in = dir_in(1:len_dir_in)//'/drpsnd'//
     1                       a8_time_orig(i)//'.dat'
              i4_contains_early = 19800         
              i4_contains_late  = 23400       

          else
              write(6,*)' error - invalid c8_drpsnd_format '
     1                 ,c8_drpsnd_format    
              istatus = 0
              goto 999

          endif

!         filename_in = 'test.nc                                 '

!         define limits of dropsonde data times we are interested in
          i4time_drpsnd_latest =   i4time_sys + i4_drpsnd_window 
          i4time_drpsnd_earliest = i4time_sys - i4_drpsnd_window 

!         define limits of file times we are interested in. 
          i4_filetime_latest =   i4time_drpsnd_latest+i4_contains_early
          i4_filetime_earliest = i4time_drpsnd_earliest-i4_contains_late     
          
          if(i .eq. 1)then
              write(6,*)' i4 drpsnd sys/window'
     1                   ,i4time_sys,i4_drpsnd_window
              write(6,*)' i4 drpsnd range     '
     1                   ,i4time_drpsnd_earliest,i4time_drpsnd_latest
              write(6,*)' i4 file range     '
     1                   ,i4_filetime_earliest,i4_filetime_latest
          endif

          if(i4times(i) .lt. i4_filetime_earliest)then
              write(6,*)' file is too early ',a9_time,i

          elseif(i4times(i) .gt. i4_filetime_latest)then
              write(6,*)' file is too late ',a9_time,i

          else
              write(6,*)
              write(6,*)' file is in time window ',a9_time,i
              write(6,*)' input file ',filename_in
              write(6,*)

!             read from the raw file and write to the opened snd file
              if(c8_drpsnd_format(1:6) .eq. 'nimbus' .or.
     1           c8_drpsnd_format(1:3) .eq. 'wfo'         )then

                  write(6,*)' calling get_drpsnd_data...'

                  call get_drpsnd_data(i4time_sys,ilaps_cycle_time
     1                ,nx_l,ny_l
     1                ,i4time_drpsnd_earliest,i4time_drpsnd_latest
     1                ,filename_in,isnd_staname,lun_out,l_fill_ht
     1                ,istatus)

              elseif(c8_drpsnd_format(1:5) .eq. 'avaps')then
                  call avapsread_sub(filename_in, lun_out
     1                          ,i4time_drpsnd_earliest
     1                          ,i4time_drpsnd_latest,istatus)

              elseif(c8_drpsnd_format(1:3) .eq. 'rsa')then

                  write(6,*)' dropsonde access routine not yet set up'       

!                 call get_raob_data   (i4time_sys,ilaps_cycle_time
!    1                ,nx_l,ny_l
!    1                ,i4time_drpsnd_earliest,i4time_drpsnd_latest
!    1                ,filename_in,istatus)

              elseif(c8_drpsnd_format(1:3) .eq. 'snd')then

                  write(6,*)' calling combine_snd_file...'       

                  call combine_snd_file(i4times(i),i4time_sys
     1                                 ,i4time_drpsnd_earliest
     1                                 ,i4time_drpsnd_latest
     1                                 ,filename_in
     1                                 ,nx_l,ny_l,nz_l
     1                                 ,lun_out,istatus)

              elseif(c8_drpsnd_format(1:3) .eq. 'cwb')then
                  call get_drpsnd_data_cwb(i4time_sys, ilaps_cycle_time,       
     ~                 nx_l, ny_l, 
     ~                 i4time_drpsnd_earliest,i4time_drpsnd_latest,
     ~                 a9_time, filename_in, lun_out, istatus)

              else
                  write(6,*)' error - invalid c8_drpsnd_format '
     1                     ,c8_drpsnd_format
                  istatus = 0
                  goto 999

              endif

          endif

      enddo

      go to 999

 997  write(6,*)' error in dropsonde ingest (ob time windows)'

      go to 999

 998  write(6,*)' error opening output sounding file: '

 999  continue

      write(6,*)' end of dropsonde ingest routine'

      return
      end

      subroutine combine_snd_file(i4time_file,i4time_sys
     1                           ,i4time_drpsnd_earliest
     1                           ,i4time_drpsnd_latest
     1                           ,filename_in,ni,nj,nk
     1                           ,lun_out,istatus)       

      use mem_grid, only: topo

      integer max_pr, max_pr_levels
      parameter (max_pr=100)
      parameter (max_pr_levels=1000)

      integer nlevels_obs_pr(max_pr)
      integer i4time_ob_pr(max_pr)
      integer i4time_ob_pr_lvl(max_pr,max_pr_levels)
      integer iwmostanum(max_pr)

      real stalat(max_pr),stalat_lvl(max_pr,max_pr_levels)
      real stalon(max_pr),stalon_lvl(max_pr,max_pr_levels)
      real staelev(max_pr)

      real height_m(max_pr,max_pr_levels)
      real pressure_mb(max_pr,max_pr_levels)
      real ob_pr_u_obs(max_pr,max_pr_levels)
      real ob_pr_v_obs(max_pr,max_pr_levels)
      real dir_deg(max_pr,max_pr_levels)
      real spd_mps(max_pr,max_pr_levels)
      real temp_c(max_pr,max_pr_levels)
      real dewpoint_c(max_pr,max_pr_levels)

      real heights_3d(ni,nj,nk)
      real lat(ni,nj)
      real lon(ni,nj)

      character*5 c5_name(max_pr)
      character*8 c8_obstype(max_pr)
      character*9 a9time_ob(max_pr),a9time_ob_lvl(max_pr,max_pr_levels)       

      character*(*) filename_in

      logical l_fill_ht, l_snd2

      l_snd2 = .true.

      lun_in = 59
      mode = 1
      n_profiles = 0

      l_fill_ht = .false. ! decide whether to fill in heights that are read in
                          ! from the snd file 

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)return
      heights_3d = r_missing_data
      lat = r_missing_data
      lon = r_missing_data

      staelev = -999. ! missing value for dropsondes

      write(6,*)' call read_snd_data for ',filename_in

      if(.not. l_snd2)then
!         call read routine for input snd file
          call read_snd_data(lun_in,i4time_file,filename_in            ! i
     1                         ,max_pr,max_pr_levels                   ! i
     1                         ,lat,lon,topo,imax,jmax,kmax            ! i
     1                         ,heights_3d,l_fill_ht                   ! i
     1                         ,mode                                   ! i
     1                         ,n_profiles                             ! i/o
     1                         ,nlevels_obs_pr,stalat,stalon,staelev   ! o
     1                         ,c5_name,i4time_ob_pr,c8_obstype        ! o
     1                         ,height_m,pressure_mb                   ! o
     1                         ,ob_pr_u_obs,ob_pr_v_obs                ! o
     1                         ,temp_c,dewpoint_c                      ! o
     1                         ,istatus)                               ! o

      else
          call read_snd_data2(lun_in,i4time_file,filename_in           ! i
     1                         ,max_pr,max_pr_levels                   ! i
     1                         ,lat,lon,topo,imax,jmax,kmax            ! i
     1                         ,heights_3d,l_fill_ht                   ! i
     1                         ,mode                                   ! i
     1                         ,n_profiles                             ! i/o
     1                         ,staelev                                ! o
     1                         ,nlevels_obs_pr                         ! o
     1                         ,c5_name,c8_obstype                     ! o
     1                         ,height_m,pressure_mb                   ! o
     1                         ,ob_pr_u_obs,ob_pr_v_obs                ! o
     1                         ,temp_c,dewpoint_c                      ! o
     1                         ,stalat_lvl,stalon_lvl,i4time_ob_pr_lvl ! o
     1                         ,istatus)                               ! o
      endif ! l_snd2

      if(istatus .ne. 1)then
          write(6,*)' warning: bad istatus from read_snd_data, '
     1             ,'abort execution of combine_snd'
          return
      endif

      do i = 1,n_profiles
          
          write(6,*)' number of levels = ',i,nlevels_obs_pr(i)

          do j = 1,nlevels_obs_pr(i)
              if(ob_pr_u_obs(i,j) .ne. r_missing_data .and.
     1           ob_pr_v_obs(i,j) .ne. r_missing_data       )then
                  call uv_to_disp(ob_pr_u_obs(i,j),ob_pr_v_obs(i,j)
     1                           ,dir_deg(i,j),spd_mps(i,j))
              else
                  dir_deg(i,j) = r_missing_data
                  spd_mps(i,j) = r_missing_data
              endif

              if(.not. l_snd2)then
                  call make_fnam_lp(i4time_ob_pr(i),a9time_ob_lvl(i,j)
     1                             ,istatus)       
              else
                  call make_fnam_lp(i4time_ob_pr_lvl(i,j)
     1                             ,a9time_ob_lvl(i,j),istatus)       
!                 write(6,*)i,j,a9time_ob_lvl(i,j),i4time_ob_pr_lvl(i,j)
              endif

              if(istatus .ne. 1)then
                  write(6,*)' error in i4time ',i4time_ob_pr(i)
                  return
              endif

              if(.not. l_snd2)then
                  stalat_lvl(i,j) = stalat(i)
                  stalon_lvl(i,j) = stalon(i)
              endif

          enddo ! j

!         determine effective sounding time (closest level to systime)
          call get_closest_a9time(i4time_sys,a9time_ob_lvl(i,:)
     1                           ,i4time_closest,nlevels_obs_pr(i)
     1                           ,istatus)

          i4time_diff = abs(i4time_sys - i4time_closest)

          if(i4time_closest .ge. i4time_drpsnd_earliest .and.
     1       i4time_closest .le. i4time_drpsnd_latest         )then

              write(6,*)' i/i4time_diff = ',i,i4time_diff
     1                 ,' inside time bounds'

!             call write routine to append this sounding to output snd file
              call open_ext(lun_out,i4time_sys,'snd',istatus)
              if(istatus .ne. 1)then
                  write(6,*)
     1                 ' could not open output snd file with open_ext'   
                  return
              endif

              call write_snd(lun_out                        ! i
     1                    ,1,max_pr_levels,1                ! i
     1                    ,iwmostanum(i)                    ! i
     1                    ,stalat_lvl(i,:),stalon_lvl(i,:)  ! i
     1                    ,staelev(i)                       ! i
     1                    ,c5_name(i),a9time_ob_lvl(i,:)    ! i
     1                    ,c8_obstype(i)                    ! i
     1                    ,nlevels_obs_pr(i)                ! i
     1                    ,height_m(i,:)                    ! i
     1                    ,pressure_mb(i,:)                 ! i
     1                    ,temp_c(i,:)                      ! i
     1                    ,dewpoint_c(i,:)                  ! i
     1                    ,dir_deg(i,:)                     ! i
     1                    ,spd_mps(i,:)                     ! i
     1                    ,istatus)                         ! o

          else
              write(6,*)' i/i4time_diff = ',i,i4time_diff
     1                 ,' outside time bounds'

          endif ! write out this sounding

      enddo ! i

      return
      end

      subroutine get_closest_a9time(i4time_ref,a9times,i4time_closest
     1                             ,ntimes,istatus)

      character*9 a9times(ntimes)

      i4_min_diff = 99999
      i4time_closest = 99999

      do i = 1,ntimes
          call i4time_fname_lp(a9times(i),i4time_i,istatus) 
          i4_diff = abs(i4time_i - i4time_ref)
          if(i4_diff .lt. i4_min_diff)then
              i4_min_diff = i4_diff
              i4time_closest = i4time_i
          endif
      enddo ! i      

      return
      end
