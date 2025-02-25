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
       program vrc_driver
       
       character*80  c_grid_fname
       character*3   c_raddat_type_vrc
       character*200 c_dataroot

       call get_grid_dim_xy(nx_l,ny_l,istatus)
       if(istatus.eq.1)then
         write(*,*)'laps parameters obtained'
       else
          write(*,*)'istatus = ',istatus,'error - get_grid_dim_xy'
          write(*,*)'terminating laps-vrc. wsi remapping'
          stop
       endif
       call get_raddat_type(c_raddat_type_vrc,istatus)
       if(istatus.ne.1)then
          write(*,*)'istatus = ',istatus,'error - get_radar_type'
          write(*,*)'terminating laps-vrc. wsi remapping'
          stop
       endif
       call find_domain_name(c_dataroot,c_grid_fname,istatus)
       if(istatus.ne.1)then
          write(*,*)'istatus = ',istatus,'error - find_domain_name'
          write(*,*)'terminating laps-vrc. wsi remapping'
          stop
       endif

       call vrc_driver_sub(nx_l,ny_l,c_raddat_type_vrc,
     +c_grid_fname)
 
       stop
       end

       subroutine vrc_driver_sub(nx_l,ny_l,c_raddat_type,
     +c_grid_fname)
c
c program drives transformation of wsi-nowrad high density (hd) radar to laps
c domain (subroutine nowrad_to_laps). 'hd' files are assumed to reside in
c /public/data/radar/wsi/nowrad/netcdf. program also handles wsi-nowrad in
c wfo (c_raddat_type = 'wfo').
c
       integer extreme_thrsh_70
       integer extreme_thrsh_47
       integer maxradars
       integer max_files

       include 'lapsparms.for' ! max_radar_files
       parameter (max_files=max_radar_files)

       parameter (extreme_thrsh_47=0.30,
     &            extreme_thrsh_70=0.10,
     &            maxradars=200)

       character*150 dir_vrc
       character*31 ext_vrc
       character*200 dir_static
       character*(*) c_grid_fname
       character*200 aoml_path_in
       character*7 vrc_outdir

       character*125 comment_ll(2)
       character*125 comment_vrc(2)
       character*10  units_ll(2)
       character*10  units_vrc(2)
       character*3   var_ll(2)
       character*3   var_vrc(2)

       character*4 lvl_coord_2d(2)

       logical     l_archive

       real lat(nx_l,ny_l)
       real lon(nx_l,ny_l)
       real rdbz(nx_l,ny_l)
       real grid_spacing
       real data(nx_l,ny_l,2)
       real radar_dist_min(nx_l,ny_l)
       real rdum

       real percent_extreme_47
       real percent_extreme_70

       integer i4time_cur
       integer i4time_latest_diff
       integer i4time_data
       integer i4time_latest_vrc
       integer i4time_latest_wsi
       integer i4time_now_gg
       integer i4_validtime
       integer istatus
       integer laps_cycle_time
       integer lvl_2d(2)
       integer n,nn,nd, len
       integer n_vars_req
       integer irad
       integer msngrad,i4_check_interval
       integer i4_total_wait,i4_thresh_age
       integer i4times_recent(max_files)
       integer i4times_proc(max_files)

       character*100 c_values_req
       character*40  c_vars_req

       character*13 c_fname_cur
       character*13 c_fname_cur_temp
       character*13 cvt_i4time_wfo_fname13
       character*14 c_filetime
       character*9 wfo_fname13_to_fname9
       character*9 c_filename
       character*9 c_fname_sys
       character*200 wsi_dir_path
       character*255 c_filespec
       character*200 c_filenames_proc(max_files)
       character*200 c_fnames_recent(max_files)
       character*3   c_raddat_type
     
       data lvl_2d/0,0/

       l_archive = .false.
c
c get vrc runtime parameters
c
       call read_vrc_nl(wsi_dir_path,msngrad,i4_check_interval,
     +i4_total_wait,i4_thresh_age,aoml_path_in,vrc_outdir,istatus)

c
c set filename. wfo data is 13 character. however, if reading from wfo-type
c data from /public then we must fool the filename to still be 9 characters.
c we do this by checking for "public" when the radar type is wfo.
c
       iwait = 1

 50    i4time_cur = i4time_now_gg()
       call get_systime(i4time_sys,c_fname_sys,istatus)
       call make_fnam_lp(i4time_cur,c_fname_cur,istatus)

c this is designed to allow archive data runs!
       call get_laps_cycle_time(laps_cycle_time,istatus)
       if(i4time_cur-i4time_sys .gt. 4*laps_cycle_time)then
          print*,'set current time to contents of systime.dat'
          c_fname_cur=c_fname_sys
          i4time_cur=i4time_sys
          l_archive = .true.
       endif

       c_fname_cur_temp = cvt_i4time_wfo_fname13(i4time_cur)

       if(c_raddat_type.eq.'wfo'.and.wsi_dir_path(2:7).ne.
     1'public')then
          irad = 2
          c_fname_cur = c_fname_cur_temp
       elseif(c_raddat_type.eq.'wfo'.and.wsi_dir_path(2:7).eq.
     1'public')then
          c_fname_cur = wfo_fname13_to_fname9(c_fname_cur_temp)
          irad = 2
       else
          c_fname_cur = wfo_fname13_to_fname9(c_fname_cur_temp)
          irad = 1
       endif
       write(6,*)'current (nominal) time: ',c_fname_cur
c
       n=index(wsi_dir_path,' ')
       write(6,*)'wsi_dir_path = ',wsi_dir_path(1:n-1)
c
c
c convert to fname9 and determine if this time has already been processed
c
!      call get_directory('vrc',dir_vrc,nd)
       i_vrc = 1
       write(6,*)
       call get_vrc_full_path(i_vrc,vrc_outdir,dir_vrc,nd,istatus)
       write(6,*)' wsi vrc output path is ',dir_vrc(1:nd)

c       dir_vrc = '../lapsprd/vrc/'    !this also use below for output
c       nd = index(dir_vrc,' ')-1
       c_filespec = dir_vrc(1:nd)//'*'
       write(6,*)'latest time for vrc files'

       call get_latest_file_time(c_filespec,i4time_latest_vrc)

       write(6,*)' returned from get_latest_file_time for '
     1            ,dir_vrc(1:nd)
       write(6,*)
c
c find files matching filespec and that have yet to be processed      
c
       n=index(wsi_dir_path,' ')
       if(c_raddat_type.eq.'wfo')then
          c_filespec=wsi_dir_path(1:n-1)//'*'
       elseif(.not.l_archive)then
          c_filespec=wsi_dir_path(1:n-1)//'*_hd'
       else
          c_filespec=wsi_dir_path(1:n-1)//c_fname_cur(1:8)//'*_hd'
       endif

       call get_file_times(c_filespec,max_files,c_fnames_recent
     1                    ,i4times_recent,i_nbr_files_out,istatus)

       n_new_files = 0

       if(i_nbr_files_out .eq. 0)then
           write(6,*)' no hd files meeting filespec: ',trim(c_filespec)
           goto 60
       endif

       do k = 1,i_nbr_files_out
           if(l_archive .or. 
     1       ((.not. l_archive) .and. 
     1        (i4times_recent(k) .ge. i4time_cur-5400) ) 
     1                                                   )then
             if(i4times_recent(k) .gt. i4time_latest_vrc)then
               write(6,*)' recent unprocessed file is: '
     1                  ,trim(c_fnames_recent(k))
               n_new_files = n_new_files + 1
               c_filenames_proc(n_new_files) = c_fnames_recent(k)
               i4times_proc(n_new_files) = i4times_recent(k)
             else
               write(6,*)' recent processed file is:   '
     1                  ,trim(c_fnames_recent(k))
             endif                                             

           else
             write(6,*)' file is too old for processing: '
!    1                 ,i4times_recent(k),i4time_cur-5400
     1                 ,trim(c_fnames_recent(k))

           endif ! recent enough to consider processing
       enddo ! k

       write(6,*)' n_new_files is ',n_new_files

60     if(n_new_files .eq. 0 .and. iwait .le. 4)then
           write(6,*)' waiting for 60 seconds...'
           write(6,*)
           call sleep(60)
           iwait = iwait + 1
           goto 50
       endif

       nfiles = n_new_files

       do ifile = 1,nfiles

         i4time_latest_wsi = i4times_proc(ifile)         

         if(c_raddat_type.eq.'wfo'.and.wsi_dir_path(2:7).ne.'public'       
     1)then
           c_filetime=cvt_i4time_wfo_fname13(i4time_latest_wsi)
           nn=index(c_filetime,' ')-1
           c_filenames_proc(ifile)=wsi_dir_path(1:n-1)//c_filetime(1:nn)      
         else
           call make_fnam_lp (i4time_latest_wsi, c_filename, istatus)
           c_filetime=c_filename
           nn=index(c_filetime,' ')-1
           c_filenames_proc(ifile)=wsi_dir_path(1:n-1)//c_filetime(1:nn)
     1                                                //'_hd'
         endif

       enddo ! ifile
c
c this for output.  laps vrc files as indicated.
c
       ext_vrc = 'vrc'
       var_vrc(1) = 'ref'
       var_vrc(2) = 'dis'
       units_vrc(1) = 'dbz'
       units_vrc(2) = 'm'
       comment_vrc(1) = 'wsi nowrad data remapped to laps domain'
       comment_vrc(2) = 'wsi nowrad data: minimum radar distance'
c
c definitions needed for acquiring laps latitude and longitude arrays.
c
       call get_directory('static',dir_static,len)
       var_ll(1)='lat'
       var_ll(2)='lon'

       call rd_laps_static(dir_static, c_grid_fname, nx_l, ny_l, 2,
     &     var_ll, units_ll, comment_ll, data, grid_spacing,
     &     istatus)

       if(istatus.eq.1)then
          write(*,*)'laps lat/lon grid obtained'
          write(*,*)
          do j=1,ny_l
          do i=1,nx_l
             lat(i,j)=data(i,j,1)
             lon(i,j)=data(i,j,2)
          end do
          end do
       else
          write(*,*)'unable to get lat/lon data'
          write(*,*)'nowrad-vrc process terminating'
          stop
       end if
c
c *****************************************************************************
c process the nowrad high density data. remap to laps domain.  generate output.
c *****************************************************************************
c
c need routine to read wsi-nowrad header and return nlines and nelems.
c      call get_wsi_parms_vrc(irad,nlines,nelems,
c    +rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum,istatus)

       do k=1,nfiles

          write(6,*)
          write(6,*)' processing: ',trim(c_filenames_proc(k))

          call read_nowrad_dims(c_raddat_type,c_filenames_proc(k),
     &         nelems,nlines)

          write(6,*)' nowrad dims are: ',nelems,nlines

          call nowradwsi_to_laps(c_raddat_type,
     &                     c_filenames_proc(k),
     &                     nlines,nelems,
     &                     nx_l,ny_l,
     &                     lat,lon,
     &                     i4_validtime,
     &                     rdbz,
     &                     maxradars,
     &                     radar_dist_min,
     &                     istatus )

         if(istatus .eq. 1)then
            write(6,*)'wsi data properly remapped'
         else
            goto 18
         endif
c
c quick qc check
c
         n_extreme=0
         do j=1,ny_l
         do i=1,nx_l
            if(rdbz(i,j) .ge. 70.)then
               n_extreme=n_extreme+1
            end if
         end do
         end do
         percent_extreme_70 = n_extreme/(nx_l*ny_l)
         n_extreme=0
         do j=1,ny_l
         do i=1,nx_l
            if(rdbz(i,j) .ge. 47.)then
               n_extreme=n_extreme+1
            end if
         end do
         end do
         percent_extreme_47 = n_extreme/(nx_l*ny_l)

         n=index(c_filenames_proc(k),' ')
         if(percent_extreme_70 .le. extreme_thrsh_70 .and.
     &      percent_extreme_47 .le. extreme_thrsh_47)then
c
c   output ... adjust i4time to 1960.
c
            i4time_data=i4_validtime+315619200

            data(:,:,1)=rdbz
            data(:,:,2)=radar_dist_min

            call write_laps_data(i4time_data,
     &                         dir_vrc,
     &                         ext_vrc,
     &                         nx_l,ny_l,1,2,
     &                         var_vrc,
     &                         lvl_2d,
     &                         lvl_coord_2d,
     &                         units_vrc,
     &                         comment_vrc,
     &                         data,
     &                         istatus)

            if(istatus.eq.1)then
               write(*,*)'vrc file successfully written'
               write(*,*)'for: ',trim(c_filenames_proc(k)(1:n-1))
               write(*,*)'i4 time: ',i4time_data
            else
               goto 14
            end if

         else

            write(6,*)'either'
            write(6,*)'more than 10% of dbz >= 70 or'
            write(6,*)'more than 25% of dbz >= 47; thus,'
            write(6,*)'not writing this vrc'
            write(*,*)'i4time: ',i4time_data,
     &              ' ftime: ',trim(c_filenames_proc(k)(1:n-1))

         end if

         goto 25

18       write(6,*)'error with the wsi data. not writing a'
         write(6,*)'vrc file for this time'
         write(6,*)'filename: ',c_filenames_proc(k)(1:n-1)

         goto 25

14       write(*,*)'error writing vrc file - terminating'
         write(*,*)'i4time: ',i4time_data,
     &           ' ftime: ',c_filenames_proc(k)(1:n-1)

25     enddo

990    write(6,*)'normal completion of wsi ingest'
       goto 1001

898    write(6,*)'error getting raw data path'
       goto 16

995    write(6,*)'error opening wait.parms'
       goto 16

999    write(*,*)'error reading systime.dat - terminating'
       write(*,*)
       goto 16
998    write(6,*)'no wsi data was found'

1001   write(6,*)
       write(6,*)' call map_aoml_sub'
       if(.true.)then
          call map_aoml_sub(nx_l,ny_l,aoml_path_in,vrc_outdir         ! i
     1                     ,i4time_sys,laps_cycle_time                ! i
     1                     ,istatus)                                  ! o
       endif

       write(6,*)'normal completion of vrc_driver'

16     stop
       end
