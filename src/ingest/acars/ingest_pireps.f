
      subroutine ingest_pireps(istatus)

!     steve albers                      original version
!     ken dritz     15-jul-1997         added call to get_grid_dim_xy to
!                                       get the values of nx_l, ny_l
!     ken dritz     15-jul-1997         pass nx_l, ny_l to get_pirep_data

!     input file 
      character*170 filename_in
      character*9 a9_time
      character*180 dir_in
      character*255 c_filespec
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)
      character*8 c8_project

!     output file
      character*13 filename13
      character*31    ext
      integer       len_dir

      character*40 c_vars_req
      character*180 c_values_req
      character*5 cftype

      call get_systime(i4time,a9_time,istatus)
      if(istatus .ne. 1)go to 999

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

      call get_laps_cycle_time(ilaps_cycle_time,istatus)
      if(istatus .eq. 1)then
          write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
      else
          write(6,*)' error getting laps_cycle_time'
          return
      endif

      lag_time_report = 3600

!     get list of input /public netcdf files
      c_vars_req = 'path_to_raw_pirep'
      call get_static_info(c_vars_req,c_values_req,1,istatus)
      if(istatus .eq. 1)then
          write(6,*)c_vars_req(1:30),' = ',c_values_req
          dir_in = c_values_req
      else
          write(6,*)' error getting ',c_vars_req
          return
      endif

!     dir_in = path_to_raw_pirep

      call s_len(dir_in,len_dir_in)

      cftype = '0100o'
!     cftype = '0005r'

      if(c8_project .eq. 'wfo')then
        c_filespec = dir_in(1:len_dir_in)
      else
        c_filespec = dir_in(1:len_dir_in)//'*'//cftype
      endif

      call get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_ret,istatus)

!     open output pin file
      if(i_nbr_files_ret .gt. 0 .and. istatus .eq. 1)then
          ext = 'pin'
          if(istatus .ne. 1)then
              write(6,*)' error opening output file'
              go to 999
          endif
      else
          write(6,*)' no raw data files identified'
          go to 999
      endif

!     determine range of ob time within file relative to file time
      if(cftype .eq. '0100o')then ! hourly obs time
          ifile_ob_b = 0
          ifile_ob_e = 3599
      elseif(cftype .eq. '0005r')then ! every 5 min reciept time     
          ifile_ob_b = -3600
          ifile_ob_e = 299
      endif

!     range of desired obs times
      i4time_ob_earliest = i4time - (ilaps_cycle_time / 2)
      i4time_ob_latest   = i4time + (ilaps_cycle_time / 2)

      write(6,*)'i4time/ob range  ',i4time,i4time_ob_earliest
     1                                    ,i4time_ob_latest

!     range of potential file time stamps having obs within time window
      i4time_file_earliest = i4time_ob_earliest - ifile_ob_e
      i4time_file_latest   = i4time_ob_latest   - ifile_ob_b

      write(6,*)'i4time/file range',i4time,i4time_file_earliest
     1                                    ,i4time_file_latest

!     loop through /public netcdf files and choose ones in time window
      write(6,*)' type / # of files = ',cftype,i_nbr_files_ret
      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)
          if(c8_project .eq. 'wfo')then
            filename_in = c_fnames(i)
          else
            filename_in = dir_in(1:len_dir_in)//a9_time//cftype
          endif
!         filename_in = 'test.nc                                 '

!         test whether we want the netcdf file for this time
          if(i4times(i) .lt. i4time_file_earliest)then
              write(6,*)' file is too early ',a9_time,i
          elseif(i4times(i) .gt. i4time_file_latest)then
              write(6,*)' file is too late ',a9_time,i
          else
              write(6,*)' file is in time window ',a9_time,i
              write(6,*)' input file ',filename_in

!             read from the netcdf pirep file and write to the opened pin file
              if(c8_project .eq. 'wfo')then
                write(6,*)' beginning awips pirep ingest - ',filename_in
                write(6,*)
                call get_pirep_data_wfo(i4time,ilaps_cycle_time
     1                                      ,filename_in
     1                                      ,ext
     1                                      ,nx_l,ny_l
     1                                      ,istatus)
                write(6,*)
                write(6,*)' end of awips pirep ingest - ', filename_in
                write(6,*)
              else
                call get_pirep_data(i4time,ilaps_cycle_time,filename_in
     1                                      ,ext
     1                                      ,nx_l,ny_l
     1                                      ,istatus)
              endif
          endif
      enddo


!     close(11) ! output pin file

 999  continue

      write(6,*)' end of pirep ingest routine'

      return
      end
 
