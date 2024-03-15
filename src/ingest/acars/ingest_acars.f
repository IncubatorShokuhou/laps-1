
      subroutine ingest_acars(istatus)

!     input file 
      character*170 filename_in
      character*9 a9_time
      character*180 dir_in
      character*3 ext_in,fname_fmt
      character*8 c8_file_fmt,c8_project
      character*255 c_filespec
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)
      character*13 filename13, cvt_i4time_wfo_fname13
      logical l_use_tamdar ! applies to non-wfo, non_rsa runs

!     output file
      character*31    ext

      character*40 c_vars_req
      character*180 c_values_req

      call get_systime(i4time_sys,a9_time,istatus)
      if(istatus .ne. 1)go to 999

!     i4time = (i4time/3600) * 3600

      call get_grid_dim_xy(nx_l,ny_l,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'error getting horizontal domain dimensions'
          go to 999
      endif

      call get_c8_project(c8_project,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'error getting c8_project'
          go to 999
      endif

      call get_l_use_tamdar(l_use_tamdar,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'error getting l_use_tamdar'
          go to 999
      endif

      lag_time_report = 3600

!     get path to input files (/public netcdf format)
      c_vars_req = 'path_to_qc_acars'
      call get_static_info(c_vars_req,c_values_req,1,istatus)
      if(istatus .eq. 1)then
          write(6,*)c_vars_req(1:30),' = ',c_values_req
          dir_in = c_values_req
      else
          write(6,*)' error getting ',c_vars_req
          return
      endif

      if(c8_project .eq. 'nimbus' .or. c8_project .eq. 'cwb')then
          ext_in = 'cdf'
          call s_len(dir_in,len_dir_in)
          c_filespec = dir_in(1:len_dir_in)//'/'//'*00q.'//ext_in

!         wait for latest input data (only for nimbus format)
          i4time_now = i4time_now_gg()
          i4_hour = (i4time_now/3600) * 3600        ! i4time at the top of the hour
          minutes_now = (i4time_now - i4_hour) / 60

          if(minutes_now .ge. 19 .and. minutes_now .lt. 22)then
              i4time_desired = i4_hour
              i4_check_interval = 10
              i4time_stop_waiting = i4_hour + 22*60
              i4_total_wait = min(i4time_stop_waiting - i4time_now, 120)
              i4_thresh_age = 3600

              call wait_for_data(c_filespec,i4time_desired
     1                   ,i4_check_interval,i4_total_wait
     1                   ,i4_thresh_age       ! only loop through the waiting
                                              ! if data is younger than this
                                              ! threshold
     1                   ,istatus)
          endif ! within time range to wait for data

          c8_file_fmt = 'nimbus'

      elseif(c8_project .eq. 'afwa')then
          ext_in = 'ac'
          c8_file_fmt = 'afwa'
          call s_len(dir_in,len_dir_in)
          c_filespec = dir_in(1:len_dir_in)//'/'//'*00q.'//ext_in

      else ! create filespec for wfo format
          ext_in = 'wfo'
          c8_file_fmt = 'wfo'
          call s_len(dir_in,len_dir_in)
          c_filespec = dir_in(1:len_dir_in)

      endif

      fname_fmt = ext_in

!     get list of files
      call get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_ret,istatus)

!     check to see if files are named (presumably) with the wfo convention 
      if(i_nbr_files_ret .eq. 0 .and. fname_fmt .ne. 'wfo')then
          write(6,*)' no raw files with filename convention of '
     1              ,fname_fmt
          write(6,*)
     1        ' try for wfo filename convention (e.g. madis acars)...'       
          ext_in = 'wfo'
          fname_fmt = 'wfo'
          call s_len(dir_in,len_dir_in)
          c_filespec = dir_in(1:len_dir_in)
          call get_file_times(c_filespec,max_files,c_fnames
     1                          ,i4times,i_nbr_files_ret,istatus)
      endif

!     open output pin file for appending
      if(i_nbr_files_ret .gt. 0 .and. istatus .eq. 1)then
          ext = 'pin'
      else
          write(6,*)' no raw data files identified of fname fmt '
     1              ,fname_fmt
          goto999
      endif

!     get acars time window
      call get_windob_time_window('acars',i4_wind_ob,istatus)
      if(istatus .ne. 1)goto 997

      call get_tempob_time_window('acars',i4_temp_ob,istatus)
      if(istatus .ne. 1)goto 997

      i4_acars_window = max(i4_wind_ob,i4_temp_ob)

!     loop through acars files and choose ones in time window
      write(6,*)' # of raw acars files = ',i_nbr_files_ret
      write(6,*)' file format = ',c8_file_fmt
      write(6,*)' filename format = ',fname_fmt

      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)
          filename_in = c_fnames(i)
!         test whether we want the netcdf file for this time
          i4time_file_earliest = i4time_sys - i4_acars_window
     1                                      - lag_time_report
          i4time_file_latest =   i4time_sys + i4_acars_window
          
          if(i4times(i) .lt. i4time_file_earliest)then
              write(6,*)
              write(6,*)' file is too early ',a9_time,i
          elseif(i4times(i) .gt. i4time_file_latest)then
              write(6,*)
              write(6,*)' file is too late ',a9_time,i
          else
              write(6,*)
              write(6,*)' file is in time window ',a9_time,i
              write(6,*)' input file ',filename_in

              if(c8_file_fmt .eq. 'nimbus')then ! nimbus netcdf
!                 read from the acars file 
!                 write to the opened pin file
!                 ingest_acars_sub.f
                  call get_acars_data(i4time_sys,i4_acars_window
     1                                      ,nx_l,ny_l
     1                                      ,c8_file_fmt
     1                                      ,ext
     1                                      ,l_use_tamdar
     1                                      ,filename_in,istatus)
              elseif(c8_file_fmt .eq. 'wfo')then ! wfo netcdf
!                 read from the acars file 
!                 write to the opened pin file
                  filename13 = cvt_i4time_wfo_fname13(i4times(i))
                  call get_acars_data(i4time_sys,i4_acars_window
     1                                      ,nx_l,ny_l
     1                                      ,c8_file_fmt
     1                                      ,ext
     1                                      ,l_use_tamdar
     1                                      ,filename_in,status)
              elseif(c8_file_fmt .eq. 'afwa')then ! afwa ascii
!                 read from the acars file 
!                 write to the opened pin file
                  call get_acars_afwa(i4time_sys,i4_acars_window
     1                                      ,nx_l,ny_l
     1                                      ,ext
     1                                      ,filename_in,istatus)
              else
                  write(6,*)' error, invalid c8_file_fmt: ',c8_file_fmt       
                  istatus = 0
                  return
              endif

          endif
      enddo

!990  close(11) ! output pin file

      go to 999

 997  write(6,*)' error in acars ingest (ob time window)'

 999  continue

      write(6,*)' end of acars ingest routine'
 
      return
      end
 
