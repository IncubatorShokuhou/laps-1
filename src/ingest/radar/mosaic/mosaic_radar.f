      program  mosaic_radar

      include  'radar_mosaic_dim.inc'

      integer   istatus
      integer   n_radars
      integer   i_window
      integer   imosaic_3d
      integer   i
      character c_radar_ext(max_radars_mosaic)*4
      character c_radar_mosaic*3, ext_dum*3

      call get_laps_config('nest7grid',istatus)

      istat = init_timer()

      n_radars_wideband = 0
      n_radars_narrowband = 0

      call get_grid_dim_xy(nx_l,ny_l,istatus)
      if(istatus .ne. 1)go to 1000

      call get_laps_dimensions(nz_l,istatus)
      if(istatus .ne. 1)go to 1000

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)go to 1000

      call get_laps_cycle_time(laps_cycle_time,istatus)
      if(istatus .ne. 1)go to 1000

      call get_max_radar_files(max_radar_files,istatus)
      if(istatus .ne. 1)go to 1000

      call get_grid_spacing_cen(grid_spacing_cen_m,istatus)
      if(istatus .ne. 1)go to 1000
c
c radar mosaic. there are two types:  wideband remapper (vxx series - rrv output)
c and rpg narrow band (wfo arena: rdr/./vrc/).
c
c namelist items

      call mosaic_radar_nl(c_radar_mosaic,n_radars,c_radar_ext,
     & i_window,mosaic_cycle_time,imosaic_3d,n_radars_wideband,
     & n_radars_narrowband,
     & istatus)      

      if(n_radars_wideband .eq. 0. .and. n_radars_narrowband .eq. 0)then    
          c_radar_mosaic = 'mrm'
      endif

      if(n_radars_wideband .ne. 0)then  
          write(6,*)     
          write(6,*)' wideband scenario'

          if(n_radars_wideband .eq. -1)then
              call get_max_radars(n_radars_wideband,istatus)
              if(istatus .ne. 1)goto1000
              write(6,*)' setting n_radars_wideband = max_radars = '      
     1                 ,n_radars_wideband
          endif

          n_radars = n_radars_wideband
          c_radar_mosaic='vxx'

!         set c_radar_ext
          do i = 1,n_radars
              if(i .le. 99)then
                  write(c_radar_ext(i),11)i
 11               format('v',i2.2)
              else
                  write(c_radar_ext(i),12)i
 12               format('v',i3.3)
              endif
          enddo ! i

!         duplicate code section

          if(n_radars .gt. max_radars_mosaic)then
              print*,'the namelist item n_radars exceeds',
     +               'the maximum number of radars allowed'
              print*,'aborting ','n_radars = ',n_radars,' max = ',
     +               max_radars_mosaic
              goto 1000
          endif

          print*,'process parameters'
          print*,'mosaic type: ',c_radar_mosaic
          print*,'n radars to mosaic: ',n_radars
          write(6,50)'radar extensions: '
     1              ,(trim(c_radar_ext(i)),i=1,n_radars)
50        format(1x,a,10(:,a,1x))
 
          call mosaic_radar_sub(nx_l,ny_l,nz_l,
     1         max_radars_mosaic,grid_spacing_cen_m,
     1         n_radars,c_radar_ext,c_radar_mosaic,i_window,
     1         r_missing_data,laps_cycle_time,mosaic_cycle_time,
     1         imosaic_3d,max_radar_files,istatus)
 
          if(istatus.ne.1)then
             print*,'error in mosaic_radar_sub'
          endif

      endif

      if(n_radars_narrowband .ne. 0)then       
          write(6,*)     
          write(6,*)' narrowband scenario'
          n_radars = n_radars_narrowband
          c_radar_mosaic='rdr'

          if(imosaic_3d .eq. 2)then
              write(6,*)' warning: inconsistent parameter settings'
              write(6,*)' n_radars_narrowband = ',n_radars_narrowband
              write(6,*)' imosaic_3d = ',imosaic_3d
              write(6,*)
          endif

!         set c_radar_ext
          do i = 1,n_radars
              write(c_radar_ext(i),21)i
 21           format(i3.3)
          enddo ! i

!         duplicate code section

          if(n_radars .gt. max_radars_mosaic)then
              print*,'the namelist item n_radars exceeds',
     +               'the maximum number of radars allowed'
              print*,'aborting ','n_radars = ',n_radars,' max = ',
     +               max_radars_mosaic
              goto 1000
          endif

          print*,'process parameters'
          print*,'mosaic type: ',c_radar_mosaic
          print*,'n radars to mosaic: ',n_radars
          write(6,50)'radar extensions: '
     1                    ,(trim(c_radar_ext(i)),i=1,n_radars)
!50       format(1x,a,10(:,a3,1x))
 
          call mosaic_radar_sub(nx_l,ny_l,nz_l,
     1         max_radars_mosaic,grid_spacing_cen_m,
     1         n_radars,c_radar_ext,c_radar_mosaic,i_window,
     1         r_missing_data,laps_cycle_time,mosaic_cycle_time,
     1         imosaic_3d,max_radar_files,istatus)
 
          if(istatus.ne.1)then
             print*,'error in mosaic_radar_sub'
          endif

      endif

      if(c_radar_mosaic .eq. 'mrm')then
          call process_mrms(nx_l,ny_l,nz_l)
      endif

1000  stop
      end

      subroutine mosaic_radar_sub(nx_l,ny_l,nz_l,mx_radars,
     1     grid_spacing_cen_m,
     1     n_radars,c_radar_ext,c_mosaic_type,i_window_size,
     1     r_missing_data,laps_cycle_time,mosaic_cycle_time,
     1     imosaic_3d,maxfiles,istatus)
c
c a laps vrc file is generated by mosaicing v'xx' type files
c (where xx = 01, 02, ..., 20; corresponding to radar ingest files v01, v02,
c  ..., v20). the v'xx' files are remapped 3d laps grid single radar files
c corresponding to a wsr88d within the laps domain.
c
c mosaics can also be made from 2d radar files using switch 'rdr' which is
c wfo type configuration.
c
c there is no path to data. the files processed are in the lapsprd subdirectory
c and get_directory satisfies the pathway requirements.
c

      integer       maxfiles

      integer       x,y,z,record

!     real        grid_ra_ref(nx_l,ny_l,nz_l,n_radars)
      real, allocatable, dimension(:,:,:,:) :: grid_ra_ref
      real, allocatable, dimension(:,:,:,:) :: grid_ra_ref_offset

      real        grid_mosaic_3dref(nx_l,ny_l,nz_l)
      real        grid_mosaic_2dref(nx_l,ny_l)
      real        grid_mosaic_2dvrc(nx_l,ny_l,2)
      real        lat(nx_l,ny_l)
      real        lon(nx_l,ny_l)
      real        topo(nx_l,ny_l)
      real        rheight_laps(nx_l,ny_l,nz_l)
      real        rlat_radar(n_radars)
      real        rlon_radar(n_radars)
      real        rheight_radar(n_radars)
      real        dist_multiradar_2d(nx_l,ny_l,n_radars)
      real        closest_radar_m(nx_l,ny_l)
      real        radius_r, grid_spacing_cen_m

      real        zcoord_of_level
      integer       lvl_3d(nz_l)
c
      character     c_filename_vxx(maxfiles,n_radars)*200
      character     path_rdr*200
      character     path*200
      character     atime*24
      character     a9time*9
      character     c_radar_ext(mx_radars)*4
      character     c_ra_ext(n_radars)*4
      character     c_directory*256
      character     c_mosaic_type*(*)
      character     c_rad_types(n_radars)

      character     cradars*3

      integer       nfiles_vxx(mx_radars)
      integer       ioffset(n_radars),joffset(n_radars)
      integer       i_ra_count
      integer       i_window_size
      integer       len_dir
      integer       nx_r,ny_r,igrid_r

      logical       first_time
      logical       found_data
      logical       l_low_level
      logical       l_offset_vxx, l_offset
!     parameter (l_offset_vxx = .true.)

      integer       i4time_mos
      integer       i4time_pre
      integer       i4time_now_gg
      integer       i4time_diff
      integer       i4time_window_beg
      integer       i4time_window_end
      integer       i4timefile_vxx(maxfiles,n_radars)
      integer       i4timefile_proc(maxfiles,n_radars)
      integer       i4_file_closest(n_radars)
      integer       i4time_nearest
c
c vrc output definitions
c
!     character     dir_vrc*50
      character     ext_vrc*31
      character*125 comment_vrc(2)
      character*10  units_vrc(2)
      character*3   var_vrc(2)
      character     lvl_coord_2d*4
c
c vrz output definitions
c
      character     dir_vrz*50
      character     ext_vrz*31
      character     comment_vrz*240
      character     units_vrz*10
      character     var_vrz*3
c
c for getting laps heights & vrc
c
      character     ext*31
     +             ,var_3d(nz_l)*3
     +             ,lvl_coord_3d(nz_l)*4
     +             ,units_3d(nz_l)*10
     +             ,comment_3d(nz_l)*200
     +             ,comment_tmp*40

      character     units_2d*10
      character     var_2d*3, var_dis*3
      character     comment_2d*125

      character     c_radar_id(n_radars)*4
      character     c_ra_ftime(n_radars)*9
      character     c_ftime_data*9

      data          lvl_2d/0/
c
c---------------------------------------------------------
c start
c
      call get_l_offset_radar(nx_l,ny_l,grid_spacing_cen_m,     ! i
     1                        nx_r,ny_r,igrid_r,l_offset_vxx)   ! o

      if(l_offset_vxx .and. c_mosaic_type.eq.'vxx')then
         l_offset = .true.
      else
         l_offset = .false.
      endif

      if(.not. l_offset)then

         n_radars_grid = n_radars
         n_radars_offset = 0

      else 
         radius_r = 500000. ! 500km max radar radius
!        igrid_r = int(radius_r / grid_spacing_cen_m) + 1
!        nx_r = min(((2 * igrid_r) + 1),nx_l)
!        ny_r = min(((2 * igrid_r) + 1),ny_l)


         n_radars_grid = 0
         n_radars_offset = n_radars

      endif

      allocate(grid_ra_ref(nx_l,ny_l,nz_l,n_radars_grid)
     1        ,stat=istat_alloc)
      if(istat_alloc .ne. 0)then
         write(6,*)' error: could not allocate grid_ra_ref'      
         stop
      endif

      allocate(grid_ra_ref_offset(nx_r,ny_r,nz_l,n_radars_offset)
     1        ,stat=istat_alloc)             
      if(istat_alloc .ne. 0)then
         write(6,*)' error: could not allocate grid_ra_ref'      
                     stop
      else
         write(6,*)' allocated grid_ra_ref_offset ',nx_r,ny_r
      endif

      istatus = 0
      l_low_level=.false.

      closest_radar_m = r_missing_data ! initialize as a placeholder
c
c get current time. make the time window.
c --------------------------------------------------------
      call get_systime(i4time_sys,a9time,istatus)
c
c get lat/lon/topo data
c
      call get_laps_domain(nx_l,ny_l,'nest7grid',
     &                     lat,lon,topo,istatus)
      if(istatus .ne. 1)then
          write(6,*)'error reading static file'
          goto 1000
      end if
      write(6,*)
c
c determine if data from any radar is current. count the number of
c radars with current radar data. terminate if no current data available.
c
      if(c_mosaic_type.eq.'rdr')then
         call get_directory('rdr',path_rdr,lprdr)
      endif

      do i=1,n_radars

         if(c_mosaic_type.eq.'vxx')then
            call get_directory(trim(c_radar_ext(i)),path,lenp)
         elseif(c_mosaic_type.eq.'rdr')then
            path=path_rdr(1:lprdr)//trim(c_radar_ext(i))//'/vrc/'
            call s_len(path,lenp)
         endif

!        should this be simplified with a call to 'get_file_times'?
         call get_file_names(path,
     &                    numoffiles,
     &                    c_filename_vxx(1,i),
     &                    maxfiles,
     &                    istatus)
         if(istatus.eq.1)then
            print*,'success in get_file_names. numoffiles = ',numoffiles
            if(numoffiles .le. 0)then
               write(6,*)'no data available in: ',trim(c_radar_ext(i))
               nfiles_vxx(i)=0
               goto 333
            end if
c
c laps internal filename convention always applies (yyjjjhhmm.ext).
c
            do l=1,numoffiles
               nn=index(c_filename_vxx(l,i),' ')
c              write(6,*)c_filename_vxx(l,i)(1:nn)
               if(i .ge. 100)then
                  call cv_asc_i4time(c_filename_vxx(l,i)(nn-14:nn-6),
     &                               i4timefile_vxx(l,i))
               else
                  call cv_asc_i4time(c_filename_vxx(l,i)(nn-13:nn-5),
     &                               i4timefile_vxx(l,i))
               endif
            end do
            nfiles_vxx(i)=numoffiles
         else
            write(6,*)'istatus ne 1 in getfilenames - abort'
            stop
         end if

333   enddo ! i

      i4_elapsed = ishow_timer()

      if(mosaic_cycle_time .ge. 300 .and. 
     1   mosaic_cycle_time .le. laps_cycle_time)then

          write(6,*)' we have a valid input mosaic_cycle_time of'
     1             ,mosaic_cycle_time
      else
          mosaic_cycle_time = laps_cycle_time
      endif

      num_mosaics = laps_cycle_time / mosaic_cycle_time

      if(laps_cycle_time .eq. mosaic_cycle_time * num_mosaics)then
          write(6,*)' mosaics divide evenly ',num_mosaics
      else
          mosaic_cycle_time = laps_cycle_time
          num_mosaics = 1
      endif

      i4time_start = i4time_sys - mosaic_cycle_time * (num_mosaics-1)       
      i4time_end   = i4time_sys

      do i4time_mos = i4time_start,i4time_end,mosaic_cycle_time

          i_ra_count = 0
          write(6,*)' mosaic time loop = ',i4time_start
     1                                    ,i4time_mos,i4time_end

          i4time_window_beg = i4time_mos-i_window_size
          i4time_window_end = i4time_mos+i_window_size

          found_data=.false.
          do i=1,n_radars
              first_time=.true.
              do l=1,nfiles_vxx(i)
                  if(i4timefile_vxx(l,i).ge.i4time_window_beg.and.
     &               i4timefile_vxx(l,i).le.i4time_window_end)then

                    if(first_time)then
                      found_data=.true.
                      first_time=.false.
                      i_ra_count=i_ra_count+1
                      c_ra_ext(i_ra_count) = c_radar_ext(i)

                      i4_file_closest(i_ra_count) = i4timefile_vxx(l,i)
                      i4_diff_min = abs(i4timefile_vxx(l,i)-i4time_mos)       
                      call make_fnam_lp(i4timefile_vxx(l,i),
     &                                  c_ra_ftime(i_ra_count),istatus)      

                    else   ! this switch= if more than one file for same radar 
                           ! within window. determine closest
                      i4_diff = abs(i4timefile_vxx(l,i)-i4time_mos)
                      if(i4_diff .lt. i4_diff_min)then
                        i4_file_closest(i_ra_count)=i4timefile_vxx(l,i)      
                        i4_diff_min = i4_diff
                        call make_fnam_lp(i4timefile_vxx(l,i),
     &                                  c_ra_ftime(i_ra_count),istatus)      
                      endif

                    endif ! first_time

                  endif ! in time window

              enddo ! file

              write(6,*)trim(c_radar_ext(i)), ' found radar = ', 
     1                  .not. first_time       
          enddo ! radar

          if(.not.found_data)then
             write(6,*)
     1           'no current files in any radar directories of type '    
     1           ,c_mosaic_type
             write(6,*)'no data will be processed for this time'
             goto 895
          endif
c
c determine appropriate i4time for these data
c
          if(i_ra_count.gt.1)then
             i4time_data = i4time_mos
             call make_fnam_lp(i4time_data,c_ftime_data,istatus)
             do i=1,i_ra_count
                print*,'radar info: ',i,' ',trim(c_ra_ext(i)),' '
     1                ,i4_file_closest(i_ra_count),c_ra_ftime(i)
             enddo
             print*,'data filetime: ',c_ftime_data
          elseif(i_ra_count.eq.1)then
             i4time_data = i4time_mos
             call make_fnam_lp(i4time_data,c_ftime_data,istatus)
             print*,'radar filetime: ',c_ra_ftime(1)
             print*,'data filetime: ',c_ftime_data
          else
             write(6,*)'ooops, i_ra_count = 0!'
             goto 1000
          endif

          write(6,*)'get the data'
c
c read analysis heights
c
          ext = 'lt1'
          call get_directory(ext,c_directory,len_dir)
          do k = 1,nz_l
             lvl_3d(k) = nint(zcoord_of_level(k))/100
             lvl_coord_3d(k) = 'hpa'
             var_3d(k) = 'ht'
             units_3d(k) = 'm'
          enddo

          call get_file_time(c_directory,i4time_data,i4time_nearest)
          if(i4time_nearest-i4time_data .gt. 3600)then
             write(6,*)'no current hgts available'
             l_low_level=.false.
          else
             write(6,*)'reading analysis heights'

             var_2d = 'ht'
             i4_tol = 7200
             call get_laps_3dgrid(i4time_data,i4_tol,i4time_nearest,
     1                            nx_l,ny_l,nz_l,ext,var_2d,units_2d,
     1                            comment_2d,rheight_laps,istatus)
c
c            call read_laps_data(i4time_nearest,c_directory,ext,
c    &                           nx_l, ny_l, nz_l,2d,rheight_laps,istatus) 
c    &                           var_3d, lvl_3d, lvl_coord_3d,
c    &                           units_3d, comment_3d, rheight_laps, istatus)

             if(istatus.ne.1)then
                write(6,*)'error reading heights '
                write(6,*)'setting l_low_level = false'
                l_low_level=.false.
             endif
          endif

          i4_elapsed = ishow_timer()
c
c these subroutines could be in loop (do i=1,n_mosaics).
c ----------------------------------------------------------
          if(c_mosaic_type(1:3).eq.'vxx')then

!            read 'vxx' file within time window around 'i4time_mos'

             write(6,*)' calling getlapsvxx'

             call getlapsvxx(nx_l,ny_l,nz_l,n_radars,c_radar_id,         ! i
     &          i_ra_count,c_ra_ext,i4time_mos,i_window_size,            ! i
     &          rheight_laps,lat,lon,topo,i4_file_closest,               ! i
     &          nx_r,ny_r,igrid_r,                                       ! i
     &          rlat_radar,rlon_radar,rheight_radar,n_valid_radars,      ! o
     &          grid_ra_ref,n_radars_grid,                               ! o
     &          grid_ra_ref_offset,ioffset,joffset,n_radars_offset,      ! o
     &          l_offset,                                                ! i
     &          istatus)                                                 ! o

             i4_elapsed = ishow_timer()

          elseif(c_mosaic_type(1:3).eq.'rdr')then ! rd 'vrc' files on laps grid
             ext = 'vrc'
             var_2d = 'ref'
             var_dis = 'dis'
             ilevel = 0

             n_valid_radars = 0 ! being reset for each radar time
                                ! needed (as input) only in this if block

             do i = 1,i_ra_count
                path=path_rdr(1:lprdr)//trim(c_ra_ext(i))//'/vrc/'
                grid_ra_ref(:,:,:,i) = r_missing_data  ! initialize this 3d ref

                write(6,*)

                call get_2dgrid_dname(path
     1           ,i4_file_closest(i),0,i4time_nearest
     1           ,ext,var_2d,units_2d
     1           ,comment_2d,nx_l,ny_l,grid_ra_ref(1,1,1,i)
     1           ,ilevel,istatus)

                if(istatus.ne.1 .and. istatus.ne.-1)then
                  print*,'error reading radar ',i

                else
                  n_valid_radars = n_valid_radars + 1

                  if(comment_2d(1:3) .eq. 'wsi')then
                      c_radar_id(i)(1:4)='wsi '
                      rlat_radar(i) = r_missing_data
                      rlon_radar(i) = r_missing_data
                      rheight_radar(i) = r_missing_data

                      write(6,*)
     1                   ' wsi radar data detected, read distance array'      

                      call get_2dgrid_dname(path
     1                   ,i4_file_closest(i),0,i4time_nearest
     1                   ,ext,var_dis,units_2d
     1                   ,comment_2d,nx_l,ny_l,dist_multiradar_2d(1,1,i)     
     1                   ,ilevel,istatus)

                      if(istatus .ne. 1 .and. istatus.ne.-1)then
                          write(6,*)' error reading dis array ',i
                          istatus = 0
                          goto 1000
                      endif

                  else
                      if(comment_2d(1:9)   .eq. '        ' .or. 
     1                   comment_2d(10:18) .eq. '        '      )then
                          write(6,*)' lat/lon header missing for radar '
     1                             ,i     
                          write(6,*)' reading in distance array instead'
                          rlat_radar(i) = r_missing_data
                          rlon_radar(i) = r_missing_data
                          rheight_radar(i) = r_missing_data

                          call get_2dgrid_dname(path
     1                   ,i4_file_closest(i),0,i4time_nearest
     1                   ,ext,var_dis,units_2d
     1                   ,comment_2d,nx_l,ny_l,dist_multiradar_2d(1,1,i)     
     1                   ,ilevel,istatus)

                          if(istatus .ne. 1 .and. istatus.ne.-1)then
                              write(6,*)' error reading dis array ',i
                              istatus = 0
                              goto 1000
                          endif

                      else
                          read(comment_2d(1:9),'(f9.3)')
     1                         rlat_radar(i)
                          read(comment_2d(10:18),'(f9.3)')
     1                         rlon_radar(i)
                          read(comment_2d(19:26),'(f8.0)')
     1                         rheight_radar(i)       

                      endif

                      c_radar_id(i)(1:4)=comment_2d(34:37)
                  endif

                endif

            enddo ! i

          endif

          i4_elapsed = ishow_timer()

c
c determine max reflectivity 2-d field as composite of all radar files for
c the given time. test n_valid_radars > 1. if not then no need to mosaic!
c -------------------------------------------------------------------------
          if(n_valid_radars .ge. 1)then

c this subroutine does not yet use the imosaic_3d parameter.

             call mosaic_ref_multi(n_valid_radars,n_radars,l_low_level,   ! i
     &                         c_radar_id,lat,lon,nx_l,ny_l,nz_l,         ! i
     &                         rlat_radar,rlon_radar,rheight_radar,       ! i
     &                         topo,rheight_laps,grid_ra_ref,             ! i
     &                         grid_ra_ref_offset,ioffset,joffset,        ! i
     &                         nx_r,ny_r,                                 ! i
     &                         n_radars_grid,n_radars_offset,             ! i
     &                         imosaic_3d,                                ! i
     &                         dist_multiradar_2d,                        ! i  
     &                         l_offset,                                  ! i
     &                         grid_mosaic_2dref,grid_mosaic_3dref,       ! i/o
     &                         closest_radar_m,istatus)                   ! o
             if(istatus .ne. 1)goto 1000

             i4_elapsed = ishow_timer()

          elseif(.false.)then

             print*,'only 1 radar - no mosaic'

             if(imosaic_3d.eq.0.or.imosaic_3d.eq.2)then
                call move(grid_ra_ref(1,1,1,1),grid_mosaic_3dref(1,1,1),
     &                    nx_l,ny_l)
             elseif(imosaic_3d.eq.1.or.imosaic_3d.eq.3)then
                do k=1,nz_l
                   call move(grid_ra_ref(1,1,k,1),
     &                       grid_mosaic_3dref(1,1,k),nx_l,ny_l)
                enddo
             endif

          else
              print*,'no radars?'
          endif

c check it out
c
          print*,'------------------'
          print*,'the mosaiced field'
          print*,'------------------'

          do j=1,ny_l,10
          do i=1,nx_l,10
             write(6,31)i,j,grid_mosaic_3dref(i,j,1)
          enddo
          enddo
31        format(2(2x,i4),1x,f8.1)

          write(6,*)'mosaiced column',grid_mosaic_3dref(nx_l/2,ny_l/2,:)
 
          i4_elapsed = ishow_timer()

          write(cradars,100)n_valid_radars
100       format(i3)

          do i=1,3
             if(cradars(i:i).eq.' ')cradars(i:i)='0'
          enddo
c
c vrc output... 
c
          if(imosaic_3d.eq.0.or.imosaic_3d.eq.2)then
             ext_vrc = 'vrc'

             var_vrc(1) = 'ref'
             var_vrc(2) = 'dis'

             units_vrc(1) = 'dbz'
             units_vrc(2) = 'm'

!comment next line...shouldn't reset n_radars lw 4-2-03
!            read(cradars,*)n_radars

             do ic = 1,2
                 comment_vrc(ic)='radar mosaic. type = '//c_mosaic_type       
     1                           //' '//cradars
             enddo ! ic

             call move(grid_mosaic_2dref,grid_mosaic_2dvrc(1,1,1)
     1                ,nx_l,ny_l)       
             call move(closest_radar_m,grid_mosaic_2dvrc(1,1,2)
     1                ,nx_l,ny_l)       

             call put_laps_multi_2d(i4time_data,ext_vrc,var_vrc
     1                             ,units_vrc
     1                             ,comment_vrc,grid_mosaic_2dvrc(1,1,1)       
     1                             ,nx_l,ny_l,2,istatus)

             if(istatus.eq.1)then
                write(*,*)'vrc file successfully written'
                call cv_i4tim_asc_lp(i4time_data,atime,istatus)
                write(6,*)'for: ',atime
                write(*,*)'i4 time: ',i4time_data
             else
                write(6,*)'vrc not written!'
             end if

             write(6,*)'column after vrc section'
     1                  ,grid_mosaic_3dref(nx_l/2,ny_l/2,:)

          endif

          i4_elapsed = ishow_timer()
c
c vrz output. 
c
          if(imosaic_3d.eq.1.or.imosaic_3d.eq.2)then
             call clean_radar(grid_mosaic_3dref,nx_l,ny_l,nz_l)

             write(6,*)'cleaned column'
     1                 ,grid_mosaic_3dref(nx_l/2,ny_l/2,:)

             write(6,*)' output vrz file, n_valid_radars = '
     1                ,n_valid_radars       
             ext_vrz = 'vrz'
             var_vrz = 'ref'
             units_vrz = 'dbz'
             comment_vrz='radar mosaic. type = '//c_mosaic_type//' '
     1                   //cradars

             if(.true.)then ! write radar info into comments
                 call get_directory(ext_vrz,path,len_dir)
                 write(6,11)path,ext_vrz,var_vrz
11               format(' writing 3d ',a50,1x,a5,1x,a3)

                 do k = 1,nz_l
                     units_3d(k) = units_vrz
                     lvl_3d(k) = nint(zcoord_of_level(k))/100
                     lvl_coord_3d(k) = 'hpa'
                     var_3d(k) = var_vrz
                     comment_3d(k) = " " ! initialize
                 enddo ! k

                 comment_3d(1) = comment_vrz

                 n_ref = 0

!                write comments in 6 columns each 40 characters wide
                 nch=40
                 do i_radar = 1,n_valid_radars
                     if(i_radar .le. (nz_l-1) )then ! write in 1st column
                         ii = i_radar + 1
                         write(comment_3d(ii),1)rlat_radar(i_radar)
     1                                         ,rlon_radar(i_radar)
     1                                         ,rheight_radar(i_radar)
     1                                         ,n_ref
     1                                         ,c_radar_id(i_radar)
1                        format(2f9.3,f8.0,i7,a4)

                         do ip = 1,40
                            write(6,*)' comment is ',ip,' ',
     1                             comment_3d(ii)(ip:ip)      
                         enddo ! ip

                     elseif(i_radar .le. 2*(nz_l-1) )then ! write in 2nd column
                         write(comment_tmp,1)rlat_radar(i_radar)
     1                                      ,rlon_radar(i_radar)
     1                                      ,rheight_radar(i_radar)
     1                                      ,n_ref
     1                                      ,c_radar_id(i_radar)
                         ii = i_radar - (nz_l-1) + 1
                         comment_3d(ii)(1*nch+1:2*nch-3) = comment_tmp    

                     elseif(i_radar .le. 3*(nz_l-1) )then ! write in 3rd column
                         write(comment_tmp,1)rlat_radar(i_radar)
     1                                      ,rlon_radar(i_radar)
     1                                      ,rheight_radar(i_radar)
     1                                      ,n_ref
     1                                      ,c_radar_id(i_radar)
                         ii = i_radar - (2*(nz_l-1)) + 1
                         comment_3d(ii)(2*nch+1:3*nch-3) = comment_tmp

                     elseif(i_radar .le. 4*(nz_l-1) )then ! write in 4th column
                         write(comment_tmp,1)rlat_radar(i_radar)
     1                                      ,rlon_radar(i_radar)
     1                                      ,rheight_radar(i_radar)
     1                                      ,n_ref
     1                                      ,c_radar_id(i_radar)
                         ii = i_radar - (3*(nz_l-1)) + 1
                         comment_3d(ii)(3*nch+1:4*nch-3) = comment_tmp

                     elseif(i_radar .le. 5*(nz_l-1) )then ! write in 5th column
                         write(comment_tmp,1)rlat_radar(i_radar)
     1                                      ,rlon_radar(i_radar)
     1                                      ,rheight_radar(i_radar)
     1                                      ,n_ref
     1                                      ,c_radar_id(i_radar)
                         ii = i_radar - (4*(nz_l-1)) + 1
                         comment_3d(ii)(4*nch+1:5*nch-3) = comment_tmp

                     elseif(i_radar .le. 6*(nz_l-1) )then ! write in 6th column
                         write(comment_tmp,1)rlat_radar(i_radar)
     1                                      ,rlon_radar(i_radar)
     1                                      ,rheight_radar(i_radar)
     1                                      ,n_ref
     1                                      ,c_radar_id(i_radar)
                         ii = i_radar - (5*(nz_l-1)) + 1
                         comment_3d(ii)(5*nch+1:6*nch-3) = comment_tmp

                     else
                         write(6,*)
     1                   ' error: too many radars for comment output'
                         write(6,*)' limit is ',i_radar-1
                         istatus = 0
                         goto 1000

                     endif

                 enddo ! i

                 call write_laps_data(i4time_data,path,ext_vrz,
     1                                nx_l,ny_l,nz_l,nz_l,
     1                                var_3d,lvl_3d,lvl_coord_3d,
     1                                units_3d,comment_3d,
     1                                grid_mosaic_3dref,istatus)       

            endif

          endif

          i4_elapsed = ishow_timer()

          go to 900

895       write(6,*)'no data for this time'

900   enddo ! i4time_mos (looping through times to mosaic)

      goto 1000

998   write(6,*)'error using systime.dat'

1000  if(allocated(grid_ra_ref))       deallocate(grid_ra_ref)
      if(allocated(grid_ra_ref_offset))deallocate(grid_ra_ref_offset)

      return
      end
