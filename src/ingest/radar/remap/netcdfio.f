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

 
       subroutine radar_init(i_radar,path_to_radar,path_to_vrc,itimes  ! i
     1                      ,b_missing_data                            ! i
     1                      ,i_tilt_proc                               ! i/o
     1                      ,l_realtime                                ! o
     1                      ,i_last_scan,istatus)                      ! o

       use mem_vol
 
!      open/read polar netcdf file for the proper time
       integer max_files
       parameter(max_files=20000)  ! max_radar_files

       character*150 path_to_radar,c_filespec,filename,directory
     1              ,c_fnames(max_files),c_fnames_in(max_files)
       character*15 path_to_vrc
       character*9 a9_time
       integer i4times_raw(max_files),i4times_lapsprd(max_files)
       character*2 c2_tilt
       character*4 laps_radar_ext, c3_radar_subdir
       character*8 radar_subdir                             
       logical l_multi_tilt,l_exist,l_output,l_realtime
       character*13 a13_time
       integer cvt_wfo_fname13_i4time
       integer z_bin, v_bin, radial

       character*31 station

!      integer gater, gater_hi, gatev, gatev_hi, radialr, radialr_hi,
!    +     radialv, radialv_hi, scanr, scanr_hi, scanv,
!    +     scanv_hi,nf_fid, nf_vid, nf_status

       save a9_time, i_nbr_files_raw, i_nbr_files_2nd, i_nbr_files_vol
       save i4times_raw, i4times_lapsprd, i_nbr_lapsprd_files

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
       include 'remap_constants.dat' ! for debugging only? (also for structure)
!      include 'remap.cmn' ! for debugging only

!      make this allocatable?
!      integer   max_tilts   
!      parameter (max_tilts=10) 
!      integer   z_vol(max_ref_gates*max_ray_tilt*max_tilts)
!      integer   v_vol(max_vel_gates*max_ray_tilt*max_tilts)

!      this call is still needed for return of 'laps_radar_ext/c3_radar_subdir'
!      we could change this to pass these in through the 'radar_init' call
       call get_remap_parms(i_radar,n_radars_remap                        ! i/o
     1                    ,max_times,path_to_radar                        ! o
     1                    ,laps_radar_ext,c3_radar_subdir                 ! o
     1                    ,path_to_vrc                                    ! o
     1                    ,ref_min,min_ref_samples,min_vel_samples,dgr    ! o
     1                    ,namelist_parms,istatus)                        ! o
       if(istatus .ne. 1)then
           write(6,*)'warning: bad status return from get_remap_parms'       
           return
       endif

!      this is nominally set to zero if input filenames are in utc
       i4_input_offset = 0 ! 9*3600

       call get_systime_i4(i4time_sys,istatus)
       if(istatus .ne. 1)then ! use wall clock time if systime is not available
           i4time_sys = i4time_now_gg()
           l_realtime = .true.

       else ! success reading 'systime.dat'
           i4_age = i4time_now_gg() - i4time_sys
           if(i4_age .gt. 43200)then ! archive case assumed
               l_realtime = .false.
           else
               l_realtime = .true.
           endif
       endif

       call get_laps_cycle_time(laps_cycle_time,istatus)   
       if(istatus .ne. 1)return   

       call get_r_missing_data(r_missing_data,istatus)
       if(istatus .ne. 1)then
           write(6,*)' error calling get_r_missing_data'
           return
       endif

!      c8_fname_format = 'unknown'
c
c      determine output filename extension
       write(6,*)' radar_init: laps_ext = ',laps_radar_ext
       if(laps_radar_ext(1:1) .ne. 'v')then ! sanity check
           laps_radar_ext = 'v01'
       endif

       if(laps_radar_ext .eq. 'vrc')then 
           radar_subdir = c3_radar_subdir
           write(6,*)' radar_init: radar_subdir = ',radar_subdir
       endif

       i_last_scan = 0

       if(i_tilt_proc .lt. 10)then
           write(c2_tilt,101)i_tilt_proc
 101       format('0',i1)
       else
           write(c2_tilt,102)i_tilt_proc
 102       format(i2)
       endif

       if(i_tilt_proc .eq. 1)then
           call s_len(path_to_radar,len_path)
           if(len_path .le. 0)then
               write(6,*)' len_path <= 0, return from radar_init'
               istatus = 0
               return
           endif
 
!          get filecount of 02 elevation raw files
           c2_tilt = '02'
           c_filespec = path_to_radar(1:len_path)//'/*elev'//c2_tilt       

           write(6,*)' itimes = ',itimes

           if(itimes .eq. 1)then
               i_nbr_files_vol = 0

               call get_file_names(c_filespec,i_nbr_files_2nd,c_fnames
     1                            ,max_files,istatus)
               if(istatus .ne. 1)then
                   write(6,*)' istatus returned from get_file_names ='        
     1                      ,istatus
                   write(6,*)' assume no directory present (tilt)'
                   return
               endif

               if(i_nbr_files_2nd .eq. 0)then ! try for volume files
                   c_filespec = path_to_radar(1:len_path)//'/*.nc'
                   write(6,*)' volume filespec = ',trim(c_filespec)
                   call get_file_names(c_filespec,i_nbr_files_vol
     1                                ,c_fnames_in,max_files,istatus)
                   if(istatus .ne. 1)then
                       write(6,*)
     1                       ' istatus returned from get_file_names ='  
     1                        ,istatus
                       write(6,*)' assume no directory present (vol)'
                       return
                   else
                       write(6,*)' i_nbr_files_vol = ',i_nbr_files_vol
                   endif
               endif

           endif

           call s_len(c_filespec,lenspec)
           write(6,*)' input filespec = ',c_filespec(1:lenspec)

           write(6,*)' # of 2nd tilt raw files = ',i_nbr_files_2nd

           i4_elapsed = ishow_timer()

           if(i_nbr_files_vol .gt. 0)then
               write(6,*)' we have volume data, # files ='
     1                                        ,i_nbr_files_vol
               c8_fname_format = 'volume'
               c_filespec = path_to_radar(1:len_path)//'/*.nc'

               if(itimes .eq. 1)then ! determine file times
                   do i = 1,i_nbr_files_vol
                       if(i .eq. 1)then
                           call get_directory_length(c_fnames_in(i)
     1                                              ,lend)
                           call s_len(c_fnames_in(i),len_fname)
                           lenf = len_fname - lend
                       endif
                       a13_time = c_fnames_in(i)(lend+5:lend+17)
                       write(6,*)' vol a13_time = ',a13_time
                       i4times_raw(i) = cvt_wfo_fname13_i4time(a13_time)
                   enddo ! i
               endif

               i_nbr_files_raw = i_nbr_files_vol

               write(6,*)' # of (volume) raw files = ',i_nbr_files_raw

           else
               if(i_nbr_files_2nd .gt. 0)then
                   l_multi_tilt = .true.
                   write(6,*)' we have multiple tilt data'
               else
                   l_multi_tilt = .false.
                   write(6,*)' we have single tilt data'
               endif
 
!              get i4times of 01 elevation raw files
               c2_tilt = '01'
               c_filespec = path_to_radar(1:len_path)//'/*elev'//c2_tilt       

               if(itimes .eq. 1)then
                   call get_file_times(c_filespec,max_files,c_fnames_in
     1                                ,i4times_raw,i_nbr_files_raw
     1                                ,istatus)
               endif

               i4times_raw(1:i_nbr_files_raw) =
     1         i4times_raw(1:i_nbr_files_raw) - i4_input_offset

               write(6,*)' # of (1st tilt) raw files = ',i_nbr_files_raw

           endif

           i4_elapsed = ishow_timer()

           if(istatus .ne. 1 .or. i_nbr_files_raw .eq. 0)then
               istatus = 0
               return
           endif

!          get output filespec
           if(laps_radar_ext .ne. 'vrc')then
               call get_filespec(laps_radar_ext,1,c_filespec,istatus)

           else ! laps_radar_ext = 'vrc', now check path_to_vrc
               if(path_to_vrc .eq. 'rdr')then
                   call get_directory('rdr',directory,len_dir)
                   c_filespec = directory(1:len_dir)//radar_subdir(1:3)
     1                          //'/vrc'     
               else ! path_to_vrc = 'lapsprd'
                   call get_filespec(laps_radar_ext,1,c_filespec
     1                              ,istatus)      
               endif

           endif

           call s_len(c_filespec,lenspec)
           write(6,*)' output filespec = ',c_filespec(1:lenspec)

!          get i4times of output files 
           if(itimes .eq. 1)then
               call get_file_times(c_filespec,max_files,c_fnames
     1                            ,i4times_lapsprd,i_nbr_lapsprd_files
     1                            ,istatus)

           else ! more efficiency for subsequent radars
               i_nbr_lapsprd_files = i_nbr_lapsprd_files + 1
               i4times_lapsprd(i_nbr_lapsprd_files) = i4time_process
               if(i4time_process .eq. 0)then
                   write(6,*)' warning: i4time_process = 0'
               else
                   write(6,*)' adding to i4times_lapsprd:'
     1                      ,i4times_lapsprd(i_nbr_lapsprd_files)
               endif

           endif
               
           write(6,*)' # of output files = ',i_nbr_lapsprd_files

           if(i_nbr_lapsprd_files .ge. 1)then ! write latest output filetime
               call make_fnam_lp(i4times_lapsprd(i_nbr_lapsprd_files)
     1                          ,a9_time,istatus)
               write(6,*)' latest output filetime = ',a9_time
           endif
           
           i4_elapsed = ishow_timer()

           i4time_process = 0

           if(l_realtime .and. i_nbr_files_vol .eq. 0)then
               needed_raw_files = 2
               latest_raw_file = i_nbr_files_raw-1
           else
               needed_raw_files = 1
               latest_raw_file = i_nbr_files_raw
           endif

           write(6,*)' l_realtime/needed/latest'
     1                ,l_realtime,needed_raw_files,latest_raw_file

!          get input filetime to process
           if(i_nbr_files_raw .ge. needed_raw_files)then ! use earliest time
               i4_earliest_window = i4time_sys - laps_cycle_time - 1800       
               call make_fnam_lp(i4_earliest_window,a9_time,istatus)

               write(6,*)
     1           ' looking for earliest unprocessed input file back to '       
     1           ,a9_time

               if(l_realtime)then
                   i4_latest_window = i4time_now_gg()
               else
                   i4_latest_window = i4time_sys + laps_cycle_time + 1800
               endif
               call make_fnam_lp(i4_latest_window,a9_time,istatus)

               write(6,*)' files accepted up to ',a9_time

               i_process = 0
               do i = latest_raw_file,1,-1
                   if(i .eq. latest_raw_file)then ! write latest raw filetime
                       call make_fnam_lp(i4times_raw(i),a9_time,istatus)
                       write(6,*)' latest raw filetime = ',a9_time
                   endif

                   l_output = .false.
                   do j = 1,i_nbr_lapsprd_files
                       if(i4times_raw(i) .eq. i4times_lapsprd(j))then
                           l_output = .true.
                       endif
                   enddo ! j

                   if( (.not. l_output)                      .and. 
     1                i4times_raw(i) .ge. i4_earliest_window       
     1                              .and.
     1                i4times_raw(i) .le. i4_latest_window
     1                                                           )then
                       i4time_process = i4times_raw(i)
                       i_process = i
                   endif

               enddo

           endif

           if(i4time_process .ne. 0)then
               call make_fnam_lp(i4time_process+i4_input_offset,
     1                           a9_time,istatus)
               write(6,*)' processing file/a9time ',i_process,a9_time
           else
               write(6,*)' no new filetimes to process'
               istatus = 0
               return
           endif

       endif ! i_tilt_proc = 1

!      pull in housekeeping data from 1st tilt
       write(6,*)' radar_init: looking for file for tilt... '
     1          ,i_tilt_proc

       i_skip = 0

!      test existence of raw 'yyjjjhhmm_elevtt / yyyymmdd_hhmm.elevtt' input
 200   if(i_nbr_files_vol .eq. 0)then ! tilt data
           c8_fname_format = 'unknown'
           call check_input_file(path_to_radar,a9_time,i_tilt_proc      ! i
     1                          ,c8_fname_format                        ! i/o
     1                          ,filename,l_exist)                      ! o     
       else                           ! volume data
         l_exist = .true.
         if(i_tilt_proc .eq. 1)then   ! first tilt
!          filename = path_to_radar(1:len_path)//'/'//a9_time//'.nc'
           write(6,*)' processing volume directory file # ',i_process
           if(i_process .lt. 0 .or. i_process .gt. max_files)then
               write(6,*)' error: i_process out of bounds in radar_init'
               stop
           endif
           filename = c_fnames_in(i_process)
           write(6,*)' c_fnames array',(c_fnames_in(i),i=1,i_process)
         endif ! first tilt
       endif

       if(l_exist)then ! these calls will fill the variables in 
                       ! 'netcdfio_radar_common.inc'

         write(6,*)' c8_fname_format = ',trim(c8_fname_format)

         if(c8_fname_format .ne. 'volume')then

           write(6,*)filename
           call get_tilt_netcdf_hdr  (filename,nf_fid
     1                               ,radarname  ! is this returned?
     1                               ,sitelat                        
     1                               ,sitelon                        
     1                               ,sitealt                        
     1                               ,elevationangle
     1                               ,numradials 
     1                               ,elevationnumber
     1                               ,vcp
     1                               ,radialazim
     1                               ,resolutionv
     1                               ,gatesizev,gatesizez
     1                               ,firstgaterangev,firstgaterangez
     1                               ,max_vel_gates, max_ref_gates ! i
     1                               ,max_ray_tilt                 ! i
     1                               ,v_bin,     z_bin,     radial ! o
     1                               ,istatus)

!          equivalent names, particularly in 'netcdfio_radar_common.inc'
           numradials = radial
           ngates_ref_cdf = z_bin
           ngates_vel_cdf = v_bin

!          initialize in case they aren't present in the netcdf file
           z_scale  = r_missing_data
           z_offset = r_missing_data
           v_scale  = r_missing_data
           v_offset = r_missing_data

!          note that file remains open from call to 'get_tilt_netcdf_hdr'
           call get_tilt_netcdf_data(filename,nf_fid
     1                               ,radarname
     1                               ,sitelat                        
     1                               ,sitelon                        
     1                               ,sitealt                        
     1                               ,elevationangle
     1                               ,numradials
     1                               ,elevationnumber
     1                               ,vcp
     1                               ,r_nyquist
     1                               ,radialazim
     1                               ,z
     1                               ,v
     1                               ,resolutionv
     1                               ,gatesizev,gatesizez
     1                               ,firstgaterangev,firstgaterangez
     1                               ,z_scale, z_offset
     1                               ,v_scale, v_offset
     1                               ,v_bin,     z_bin,     radial    ! i
     1                               ,istatus)

!          use default values if needed
           if(abs(z_scale ) .gt. 1e10)z_scale  = 2.0
           if(abs(z_offset) .gt. 1e10)z_offset = 66.
           if(abs(v_scale ) .gt. 1e10)v_scale  = 2.0
           if(abs(v_offset) .gt. 1e10)v_offset = 129.

           write(6,*)' z/v scale/offset = ',z_scale,z_offset
     1                                     ,v_scale,v_offset

         else ! c8_fname_format = 'volume'

           if(i_tilt_proc .eq. 1)then
               write(6,*)' call get_vol_netcdf_hdr'
               write(6,*)trim(filename)

               call get_vol_netcdf_hdr(filename,
     +            gater, gater_hi, gatev, gatev_hi, radialr, radialr_hi,
     +            radialv, radialv_hi, scanr, scanr_hi, scanv,
     +            scanv_hi,nf_fid, nf_vid, nf_status)

               write(6,*)' call get_attribute_vol'  

               call get_attribute_vol(nf_fid,sitelat,sitelon,sitealt
     +            ,station
     +            ,istatus)

               radarname = station

               write(6,*)' lat/lon/alt/name ',sitelat,sitelon,sitealt
     1                                       ,trim(radarname)

!              deallocate old volume scan (if needed)
               if(allocated(reflectivity))deallocate(reflectivity)
               if(allocated(reflectivity_hi))deallocate(reflectivity_hi)
               if(allocated(radialvelocity))deallocate(radialvelocity)
               if(allocated(radialvelocity_hi)) 
     1                                     deallocate(radialvelocity_hi)

               if(allocated(elevationr))deallocate(elevationr)
               if(allocated(elevationr_hi))deallocate(elevationr_hi)
               if(allocated(elevationv))deallocate(elevationv)
               if(allocated(elevationv_hi))deallocate(elevationv_hi)

               if(allocated(azimuthr))deallocate(azimuthr)
               if(allocated(azimuthr_hi))deallocate(azimuthr_hi)
               if(allocated(azimuthv))deallocate(azimuthv)
               if(allocated(azimuthv_hi))deallocate(azimuthv_hi)

               if(allocated(distancer))deallocate(distancer)
               if(allocated(distancer_hi))deallocate(distancer_hi)
               if(allocated(distancev))deallocate(distancev)
               if(allocated(distancev_hi))deallocate(distancev_hi)
               if(allocated(nyquistvelocityv))
     1           deallocate(nyquistvelocityv)
               if(allocated(nyquistvelocityv_hi))
     1           deallocate(nyquistvelocityv_hi)

!              allocate new volume scan
               write(6,*)'gater,radialr,scanr = ',gater,radialr,scanr
               allocate(reflectivity(gater,radialr,scanr))
               allocate(elevationr(radialr,scanr))
               allocate(azimuthr(radialr,scanr))
               allocate(distancer(gater))

               write(6,*)'gater_hi,radialr_hi,scanr_hi = ',
     1                    gater_hi,radialr_hi,scanr_hi
               allocate(reflectivity_hi(gater_hi,radialr_hi,scanr_hi))
               allocate(elevationr_hi(radialr_hi,scanr_hi))
               allocate(azimuthr_hi(radialr_hi,scanr_hi))
               allocate(distancer_hi(gater_hi))

               write(6,*)'gatev,radialv,scanv = ',gatev,radialv,scanv
               allocate(radialvelocity(gatev,radialv,scanv))
               allocate(elevationv(radialv,scanv))
               allocate(azimuthv(radialv,scanv))
               allocate(distancev(gatev))              
               allocate(nyquistvelocityv(gatev))              

               write(6,*)'gatev_hi,radialv_hi,scanv_hi = ',
     1                    gatev_hi,radialv_hi,scanv_hi
               allocate(radialvelocity_hi(gatev_hi,radialv_hi,scanv_hi))
               allocate(elevationv_hi(radialv_hi,scanv_hi))
               allocate(azimuthv_hi(radialv_hi,scanv_hi))
               allocate(distancev_hi(gatev_hi))                   
               allocate(nyquistvelocityv_hi(gatev_hi))                   

               write(6,*)' call get_vol_netcdf_data'
               call get_vol_netcdf_data(nf_fid, gater, gater_hi, 
     +              gatev, gatev_hi,
     +              radialr, radialr_hi, radialv, radialv_hi, 
     +              scanr, scanr_hi,
     +              scanv, scanv_hi, 
     +              reflectivity, reflectivity_hi,
     +              radialvelocity, radialvelocity_hi,
     +              elevationr, elevationr_hi,
     +              elevationv, elevationv_hi,
     +              azimuthr, azimuthr_hi,
     +              azimuthv, azimuthv_hi, 
     +              distancer, distancer_hi,
     +              distancev, distancev_hi,
     +              nyquistvelocityv, nyquistvelocityv_hi)

!              use default values
               z_scale  = 2.0
               z_offset = 66.
               v_scale  = 2.0
               v_offset = 129.
               resolutionv = 1. / v_scale

               write(6,*)' z/v scale/offset = ',z_scale,z_offset
     1                                         ,v_scale,v_offset

               write(6,*)' elevationr = '   ,elevationr(1,:)
               write(6,*)' elevationr_hi = ',elevationr_hi(1,:)
               write(6,*)' elevationv = '   ,elevationv(1,:)
               write(6,*)' elevationv_hi = ',elevationv_hi(1,:)

           endif ! i_tilt_proc = 1

!          transfer tilt to tilt arrays
           write(6,*)' transfer tilt to tilt arrays'

!          radarname = 'kftg'
!          c4_radarname = 'kftg'

           elevationnumber = i_tilt_proc

!          note the _hi scans are lowest elevations

!          high elevation (r)
           if(i_tilt_proc .gt. scanr_hi .and. 
     1        i_tilt_proc .le. scanr_hi + scanr)then
               i_array = i_tilt_proc - scanr_hi
               iscr = 0   
!              irmax = 0                                ! debug
               do j = 1,radialr
                 do i = 1,gater
                   iscr = iscr + 1
                   z(iscr) = reflectivity(i,j,i_array)
!                  if(z(iscr) .gt. 0)then               ! debug
!                    irmax = i                          ! debug
!                  endif                                ! debug
                 enddo ! i
!                if(irmax .gt. 0)then                   ! debug
!                 write(6,*)j,azimuthr(j,i_array),irmax ! debug
!    1                     ,distancer(irmax)/1000.      ! debug
!                endif                                  ! debug
                 radialazim(j) = azimuthr(j,i_array)
               enddo ! j
               elevationangle = elevationr(1,i_array)
               numradials = radialr
               ngates_ref_cdf = gater
               firstgaterangez = distancer(1) / 1000.
               gatesizez = (distancer(2) - distancer(1)) / 1000.
               write(6,*)' r i_tilt_proc/i_array/elev/gsp = ',
     1                     i_tilt_proc,i_array,elevationangle,gatesizez
           endif

!          low elevation (r_hi)
           if(i_tilt_proc .le. scanr_hi)then
               i_array = i_tilt_proc                   
               iscr = 0   
!              irmax = 0                                   ! debug
               if(radialr_hi .gt. max_ray_tilt)then
                 write(6,*)' error: radialr_hi > max_ray_tilt',
     1                     radialr_hi,max_ray_tilt
                 istatus = 0
                 return
               endif
               do j = 1,radialr_hi
                 do i = 1,gater_hi
                   iscr = iscr + 1
                   z(iscr) = reflectivity_hi(i,j,i_array)
!                  if(z(iscr) .gt. 0)then                  ! debug
!                    irmax = i                             ! debug
!                  endif                                   ! debug
                 enddo ! i
!                if(irmax .gt. 0)then                      ! debug
!                 write(6,*)j,azimuthr_hi(j,i_array),irmax ! debug
!    1                     ,distancer_hi(irmax)/1000.      ! debug
!                endif                                     ! debug
                 radialazim(j) = azimuthr_hi(j,i_array)
               enddo ! j
               elevationangle = elevationr_hi(1,i_array)
               numradials = radialr_hi
               ngates_ref_cdf = gater_hi
               firstgaterangez = distancer_hi(1) / 1000.
               gatesizez = (distancer_hi(2) -  distancer_hi(1)) / 1000.
               write(6,*)' r_hi i_tilt_proc/i_array/elev/gsp = ',
     1                     i_tilt_proc,i_array,elevationangle,gatesizez
           endif
 
!          high elevation (v)
           i_v_match = 0
           do i_v = 1,scanv    
               if(elevationv(1,i_v) .eq. elevationangle)then
                   i_v_match = i_v
                   write(6,*)' v elevation match ',i_v
               endif
           enddo

           if(i_v_match .gt. 0)then                              
               i_array = i_v_match                     
               iscr = 0   
               do j = 1,radialv
                 do i = 1,gatev
                   iscr = iscr + 1
                   v(iscr) = radialvelocity(i,j,i_array)
                 enddo ! i
                 radialazim(j) = azimuthv(j,i_array)         
               enddo ! j
               elevationangle = elevationv(1,i_array)
               numradials = radialv
               ngates_vel_cdf = gatev
               firstgaterangev = distancev(1) / 1000.
               gatesizev = (distancev(2) - distancev(1)) / 1000.
               r_nyquist = nyquistvelocityv(i_array)
               write(6,*)' v i_tilt_proc/i_array/elev/nyq = ',
     1                     i_tilt_proc,i_array,elevationangle,r_nyquist       
           else
               write(6,*)' no v match at elev ',elevationangle
           endif

!          low elevation (v_hi)
           i_v_match = 0
           do i_v = 1,scanv_hi    
               if(elevationv_hi(1,i_v) .eq. elevationangle)then
                   i_v_match = i_v
                   write(6,*)' v_hi elevation match ',i_v
               endif
           enddo

           if(i_v_match .gt. 0)then                              
               i_array = i_v_match                     
               iscr = 0   
               do j = 1,radialv_hi
                 do i = 1,gatev_hi
                   iscr = iscr + 1
                   v(iscr) = radialvelocity_hi(i,j,i_array)
                 enddo ! i
                 radialazim(j) = azimuthv_hi(j,i_array)         
               enddo ! j
               elevationangle = elevationv_hi(1,i_array)
               numradials = radialv_hi
               ngates_vel_cdf = gatev_hi
               firstgaterangev = distancev_hi(1) / 1000.
               gatesizev = (distancev_hi(2) -  distancev_hi(1)) / 1000. 
               r_nyquist = nyquistvelocityv_hi(i_array)
               write(6,*)' v_hi i_tilt_proc/i_array/elev/nyq = ',
     1                     i_tilt_proc,i_array,elevationangle,r_nyquist       
           else
               write(6,*)' no v_hi match at elev ',elevationangle
           endif

           if(i_tilt_proc .eq. scanr + scanr_hi)then
               i_last_scan = 1
           endif

!          if(i_tilt_proc .eq. scanv + scanv_hi)then
!              i_last_scan = 1
!          endif

           istatus = 1

           write(6,*)' end of volume processing for tilt ',i_tilt_proc       
     1                                                    ,i_last_scan

         endif ! end of volume processing for this tilt

       elseif(i_tilt_proc .le. 20)then ! l_exist is .false.
           i_tilt_proc = i_tilt_proc + 1
           i_skip = 1
           goto 200 ! check for existence of next tilt

       else
           istatus = 0

       endif ! file exists

       if(i_tilt_proc .gt. 50)then
           write(6,*)' error in radar_init, i_tilt_proc = ',i_tilt_proc
           istatus = 0
           return
       endif

       if(istatus .eq. 1
     1             .and.
     1     .not. (laps_radar_ext .eq. 'vrc' .and. i_tilt_proc .gt. 1)
     1                                                             )then       
           if(i_tilt_proc .eq. 1)then
               write(6,201)elevationnumber, i_tilt_proc
 201           format(' elevationnumber, i_tilt_proc',2i4)

           else
               write(6,202)elevationnumber, i_tilt_proc
 202           format(' elevationnumber, i_tilt_proc',2i4
     1               ,' (upcoming tilt)')

               if(i_skip .eq. 1)then ! tilt is found and prior one didn't exist
                   if(.false.)then
                       write(6,*)
     1               ' warning: we had to skip past some missing tilts'       
                   else
                       call s_len(radarname,lenr)
                       write(6,*)' error: missing tilts present for '
     1                          ,radarname(1:lenr),' at ',a9_time
     1                          ,' ',trim(laps_radar_ext)
                       istatus = 0
                       return
                   endif
               endif

           endif

       else
           write(6,*)' could not read tilt # ',i_tilt_proc
           i_last_scan = 1

       endif

       write(6,*)

       istatus = 1
       return
       end
 
       function get_altitude()
       integer get_altitude          ! site altitude (meters)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_altitude = nint(sitealt)

       return
       end
 
 
       function get_latitude()
       integer get_latitude          ! site latitude (degrees * 100000)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_latitude = nint(sitelat*100000)
       return
       end
 
 
       function get_longitude()
       integer get_longitude         ! site longitude (degrees * 100000)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'

       character*8 c8_project

!      call get_c8_project(c8_project,istatus)
!      if(istatus .ne. 1)then
!          write(6,*)' error, no c8_project'
!          stop
!      endif

!      if(c8_project(1:3) .ne. 'cwb')then

       if(.false.)then
           get_longitude = -abs(nint(sitelon*100000))
       else
           get_longitude =      nint(sitelon*100000)
       endif

       return
       end

       subroutine get_radarname(c4_radarname,istatus)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
       character*4 c4_radarname

       c4_radarname = radarname
       call upcase(c4_radarname,c4_radarname)
       write(6,*)' c4_radarname = ',c4_radarname

       istatus = 1
       return
       end
       
 
       function get_field_num(c3_field)
       integer get_field_num
       character*3 c3_field
 
       if(c3_field .eq. 'dbz')get_field_num = 1
       if(c3_field .eq. 'vel')get_field_num = 2
 
       return
       end
 
 
       function read_radial()
 
       read_radial = 0
       return
       end
 
 
       function get_status()
       integer get_status
 
       get_status = 0
       return
       end
 
 
       function get_fixed_angle()
       integer get_fixed_angle     ! beam tilt angle (degrees * 100)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'

       if( abs(elevationangle) .le. 1e10 )then
           get_fixed_angle = nint(elevationangle * 100.)
       else
           write(6,*)' warning in get_fixed_angle, invalid value'
           get_fixed_angle = -999 ! i_missing_data
       endif

 
       return
       end
 
 
       function get_scan()
       integer get_scan            ! scan #

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_scan = elevationnumber
       return
       end
 
 
       function get_tilt()
       integer get_tilt            ! tilt #

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_tilt = elevationnumber
       return
       end
 
 
       function get_num_rays()
       integer get_num_rays

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_num_rays = numradials
       return
       end
 
 
       subroutine get_volume_time(i4time_process_ret)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       i4time_process_ret = i4time_process
       return
       end
 
 
       function get_vcp()
       integer get_vcp

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_vcp = vcp
       return
       end
 
 
       function get_azi(iray) ! azimuth * 100.
       integer get_azi

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       if( abs(radialazim(iray)) .le. 1e10 )then
           get_azi = nint(radialazim(iray)*100.)
       else
           write(6,*)' warning in get_azi, azimuth = ',iray
     1                                   , radialazim(iray)
           get_azi = -999 ! i_missing_data
       endif

       return
       end
 
 
       function get_nyquist()
       real get_nyquist                 ! nyquist velocity of the radial (m/s)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_nyquist = r_nyquist
       return
       end
 
 
       function get_number_of_gates(index)
       integer get_number_of_gates

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_number_of_gates = 0

       if(index .eq. 1)then
           get_number_of_gates = ngates_ref_cdf
       endif

       if(index .eq. 2)then
           get_number_of_gates = ngates_vel_cdf
       endif

       return
       end
 
 
       subroutine get_first_gate(index,first_gate_m,gate_spacing_m)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'

       call get_r_missing_data(r_missing_data,istatus)
       if(istatus .ne. 1)then
           write(6,*)' error calling get_r_missing_data'
           stop
       endif
 
       if(index .eq. 1)then
           if(firstgaterangez .ge. -10. .and. 
     1        firstgaterangez .lt. 1000.)then     
               first_gate_m = firstgaterangez * 1000.
           else
               write(6,*)' warning: firstgaterangez is outside range'
     1                  ,firstgaterangez
               first_gate_m = r_missing_data
           endif

           if(gatesizez .ge. 0. .and. gatesizez .lt. 1000.)then     
               gate_spacing_m = gatesizez * 1000.
           else
               write(6,*)' warning: gatesizez is outside range'
     1                  ,gatesizez
               gate_spacing_m = r_missing_data
           endif

       elseif(index .eq. 2)then
           if(firstgaterangev .ge. -10. .and. 
     1        firstgaterangev .lt. 1000.)then     
               first_gate_m = firstgaterangev * 1000.
           else
!              write(6,*)' warning: firstgaterangev is outside range'
!    1                  ,firstgaterangev
               first_gate_m = r_missing_data
           endif

           if(gatesizev .ge. 0. .and. gatesizev .lt. 1000.)then     
               gate_spacing_m = gatesizev * 1000.
           else
!              write(6,*)' warning: gatesizev is outside range'
!    1                  ,gatesizev      
               gate_spacing_m = r_missing_data
           endif

       endif

       return
       end
 
 
       function get_data_field(index, data, n_ptr, n_gates
     1                                           , b_missing_data)
       integer get_data_field

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'

!      note that z and v are integer arrays that have been read in from the
!      netcdf files as byte values. these are here converted to dbz and 
!      meters per second, respectively. this is done for one radial for
!      each call to this routine.
 
       real data(n_gates)

       if(index .eq. 1)then ! reflectivity
           do i = 1,n_gates
               call counts_to_dbz(z(n_ptr + (i-1))                       ! i
     1                           ,z_scale,z_offset,b_missing_data        ! i
     1                           ,data(i),istatus)                       ! o
           enddo

           if(istatus .ne. 1)then
               get_data_field = 0
               return
           endif

       elseif(index .eq. 2)then ! velocity
           do i = 1,n_gates
               call counts_to_vel(v(n_ptr + (i-1))                       ! i
     1                           ,b_missing_data,v_scale,v_offset        ! i
     1                           ,resolutionv,data(i))                   ! o
           enddo

       endif

       get_data_field = 1
       return
       end
 
       subroutine counts_to_dbz(zcounts,z_scale,z_offset,b_missing_data  ! i
     1                         ,dbz,istatus)                             ! o

!      convert integer z count value to dbz

!      from the netcdf header
!      z:valid_range = 2b, -2b ;     (2 through 254)
!      z:below_threshold = 0b ;      (0)
!      z:range_ambiguous = 1b ;      (1)
!      z:_fillvalue = -1b ;          (255)

       integer zcounts
       real dbz,dbz_hold,b_missing_data                         

       dbz_hold = zcounts

!      convert from signed to unsigned
       if(dbz_hold .gt. 127.) then
           print *, 'error in reflectivity: ',dbz_hold,zcounts
           istatus = 0
       endif

       if(dbz_hold .lt. 0.) then
           dbz_hold = 256. + dbz_hold
       endif

!      if(dbz_hold .eq. 1.)then 
!          dbz_hold = b_missing_data  ! range ambiguous
!      endif

       if(dbz_hold .ne. b_missing_data)then ! scale
           dbz_hold = (dbz_hold - z_offset) / z_scale
       endif

       dbz = dbz_hold

       istatus = 1
       return
       end


       subroutine counts_to_vel(vcounts,b_missing_data,v_scale          ! i
     1                         ,v_offset,resolutionv,vel_ms)            ! o

!      convert integer v count value to radial velocity (meters/sec)

       integer vcounts
       real vel_ms,vel_hold,b_missing_data                         

       vel_hold = vcounts

!      convert from signed to unsigned
       if(vel_hold .gt. 127.) then
           print *, 'error in velocity: ',vel_hold
           stop
       endif

       if(vel_hold .lt. 0.) then
           vel_hold = 256. + vel_hold
       endif

       if(vel_hold .eq. 1. .or. vel_hold .eq. 0.)then 
           vel_hold = b_missing_data  ! invalid measurement
       endif

       if(resolutionv .eq. 0.)then ! qc check
           vel_hold = b_missing_data
       endif

       if(vel_hold .ne. b_missing_data)then ! scale valid v
           vel_hold = (vel_hold - v_offset) / v_scale
       endif

       vel_ms = vel_hold
   
       return
       end

       function cvt_fname_data()
 
       cvt_fname_data = 0
       return
       end
 
 
       subroutine check_input_file(path_to_radar,a9_time,i_tilt           ! i
     1                            ,c8_fname_format                        ! i/o
     1                            ,filename,l_exist)                      ! o

       logical l_exist

       character*(*) path_to_radar
       character*(*) filename
       character*2 c2_tilt
       character*9 a9_time
       character*8 c8_fname_format
       character*13 a13_time, fname9_to_wfo_fname13

!      determine type of filename (if not yet known) and test for existence
!      of radar file for given time/tilt

       if(i_tilt .lt. 10)then
           write(c2_tilt,101)i_tilt
 101       format('0',i1)
       else
           write(c2_tilt,102)i_tilt
 102       format(i2)
       endif

       call s_len(path_to_radar,len_path)

       if(c8_fname_format .eq. 'unknown')then
           write(6,*)'check_input_file: resolving filename format'
       endif

       if(c8_fname_format .eq. 'nimbus' .or. 
     1    c8_fname_format .eq. 'unknown')then
           filename = path_to_radar(1:len_path)//'/'//a9_time//'_elev'
     1                //c2_tilt

!          test existence of radar file
           inquire(file=filename,exist=l_exist)

           if(l_exist)c8_fname_format = 'nimbus'
       endif

       if(c8_fname_format .eq. 'wfo' .or. 
     1    c8_fname_format .eq. 'unknown')then
           a13_time = fname9_to_wfo_fname13(a9_time)
           filename = path_to_radar(1:len_path)//'/'//a13_time//'.elev'
     1                //c2_tilt

!          test existence of radar file
           inquire(file=filename,exist=l_exist)

           if(l_exist)c8_fname_format = 'wfo'
       endif

!      if(c8_fname_format .eq. 'volume')then
!          a15_time = fname9_to_wfo_fname15(a9_time)
!          filename = path_to_radar(1:len_path)//'/kftg'//a15_time//'_v03.nc'
!    1                //c2_tilt
!      endif

       call s_len(filename,len_file)
       write(6,*)' check_input_file: ',filename(1:len_file),' ',l_exist       

       return
       end
