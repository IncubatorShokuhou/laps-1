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
      program lga

c     use laps_static
      implicit none
      include 'bgdata.inc'

c
c
c------------------> background model designation <-----------------------------
c
c *** the following variable designates which model to use as the background.
c *** read from standard input.
c
c        bgmodel = 1 ---> ruc (60 km native grid)
c        bgmodel = 2 ---> eta (48 km conus-c grid)
c        bgmodel = 3 ---> cwb fa (5km) model
c        bgmodel = 4 ---> sbn conus211 (eta or ruc)
c        bgmodel = 5 ---> ruc (40 km native grid)
c        bgmodel = 6 ---> avn (360 x 181 lat-lon grid)
c        bgmodel = 7 ---> eta (48 km from grib file)
c        bgmodel = 8 ---> nogaps (1.0 deg)
c        bgmodel = 9 ---> nws conus (ruc, eta, ngm, avn)
c        bgmodel =10 ---> unidata default netcdf format from gribtonc  !wni
c        bgmodel =11 ---> wrf-arw raw netcdf files (single time per file) !wni
c        bgmodel =12 ---> ecmwf with two netcdf options: esrl and fmi
c        bgmodel =13 ---> grib1 and grib2-formatted: gfs, nam, ruc, ecmwf, etc
c
c------------------> grid dimension specification <-----------------------------
c
c *** the following variables specify the various model grid dimensions.
c *** read from standard input.
c *** the first line specifies the laps grid dimensions and p level specs.
c *** the second line specifies the laps cycle time.
c *** the third line specifies the laps root path directory.
c *** the fourth line specifies the laps domain designation.
c *** each following line specifies the grid dimensions for each corresponding
c        bgmodel defined above.
c
      integer nx_laps,ny_laps,nz_laps,       !laps grid dimensions
     .          laps_cycle_time,             !laps cycle time
     .          nx_bg,ny_bg,nz_bg,
     .          lga_status                   !status returned from lga_driver
c                                            ! 1=good 0=bad
      integer   np,ntmin,ntmax
c
c
c------------------> background grid data path <--------------------------------
c
c *** the following variables specify the various model data directory paths.
c *** read from standard input.
c *** each line contains the background data directory path for each
c        corresponding bgmodel defined above.
c
      character*256 cfilespec

      character*2   gproj
      real          dlat,dlon
      real          lat0,lat1,lon0
      real          sw(2),ne(2)

c
c  this is the max number of paths allowed in nest7grid.parms
c  and should match the value in lib/lapsgrid.f 
c 
      character*256 bgpaths(maxbgmodels)
      character*256 bgpath
      integer bgmodels(maxbgmodels)
      integer bgmodel
      integer istatus, init_timer, mode, mlen
 
      character*13 cvt_i4time_wfo_fname13
c
      character*9 a9
      character*10 c_mode
      integer i4time_now,i4time_latest
      integer i4time_now_lga,i4time_now_gg
      integer i4time_lga
      integer max_files,bg_files
      integer itime_inc
      integer itime
      parameter (max_files=20000)             
      character*256 names(max_files)            ! should be basenames
      character*256 reject_names(max_files)
      integer reject_cnt
      integer accepted_files
      integer lbgp
      integer nbgmodels
      integer idum,jdum,kdum
      integer n_written
      character*14 c_ftimes_written(100)
c
c-------------------------------------------------------------------------------
c
c *** comments used in writing netcdf files only.
c
      integer istat, i,l, no_infinite_loops
c
c cmodel is really only 12 chars but the sbn netcdf carries 132
c
      character*132 cmodels(maxbgmodels)
      character*132 cmodel
      integer forecast_length
      logical use_analysis, use_forecast, use_systime
      logical ltime(-1:1)
      logical luse_sfc_bkgd
      logical smooth_fields
      logical lgb_only
c_______________________________________________________________________________
c
c read information from static/nest7grid.parms

      print *,'istatus ',istatus
      istatus=init_timer()
      call get_grid_dim_xy(nx_laps,ny_laps,istatus)
      call get_laps_dimensions(nz_laps,istatus)
      call get_laps_cycle_time(laps_cycle_time,istatus)
c
c read information from static/background.nl
c
      call get_background_info(bgpaths,bgmodels
     +,forecast_length,use_analysis,use_forecast,cmodels
     +,itime,smooth_fields,luse_sfc_bkgd,ntmin,ntmax,lgb_only)

      nbgmodels=0
      do i=1,maxbgmodels
         if(index(bgpaths(i),' ').ne.0)then
            nbgmodels=nbgmodels+1
         endif
      enddo

c
c *** initialize esat table.
c
      call es_ini
c
c *** get current time from systime.dat
c
      use_systime=.true.
      if(use_systime)then
         call get_systime(i4time_now,a9,lga_status)
         print*,'systime: ',a9,' ',i4time_now
      else
         i4time_now = i4time_now_gg()
         print*,'using i4time now'
      endif

!     set mode
!     1 - spatial interp only
!     2 - temporal interp only (not yet supported)
!     3 - both spatial and temporal interp

      call getenv('lga_mode',c_mode)
      call s_len(c_mode,mlen)
      if(mlen .gt. 0)then
          write(6,*)' obtaining lga mode from environment variable'
          read(c_mode,*)mode
      else
          mode = 3 
      endif

      write(6,*)' lga mode = ',mode

      if(mode .eq. 1)then ! spatial interp only, process one cycle ahead
          i4time_now = i4time_now + laps_cycle_time
      endif

      if(lgb_only)then
          ntmin = 0
          ntmax = 0
!     else
!         ntmin = -1
!         ntmax = +1
!         ntmin = 0
!         ntmax = +6
      endif

      ltime=.true.

      n_written = 0

      do itime_inc=ntmin,ntmax

       print*
       print*,'----------------------------------------------'
       if(itime_inc.lt.0)then
        print*,'start time-minus-one cycle time interpolation'
       elseif(itime_inc.eq.0)then
        print*,'start cycle time interpolation'
       else
        print*,'start time-plus',itime_inc,'cycle time interpolation'
       endif
       print*,'----------------------------------------------'
       print*

       i4time_now_lga=i4time_now+(itime_inc*laps_cycle_time)
       call make_fnam_lp(i4time_now_lga,a9,istatus)
       print*,'processing background data for ',a9
       print*

       bg_files=0
       i=1

       bgmodel = bgmodels(i)
       bgpath =  bgpaths(i)
       cmodel =  cmodels(i)
       lga_status = -99 

       reject_cnt=0
       bg_files=0
       no_infinite_loops=0

       do while((lga_status.le.0 .and. i.le.nbgmodels)
     +          .and. (no_infinite_loops.le.nbgmodels))


         print *,'looping through models (lga-1):',lga_status, i, 
     1                                           bg_files,reject_cnt

c        if(i.eq.nbgmodels)i=0
         if(reject_cnt.gt.bg_files )then    !.and. bg_files.gt.0) then
            i=i+1
            bgmodel = bgmodels(i)
            bgpath =  bgpaths(i)
            cmodel =  cmodels(i)
            reject_cnt = 0
         endif
         if (bgmodel .lt. 0 .or. bgmodel .gt. maxbgmodels) then
            print*
            print*,' cannot proceed with model specification in lga'
            print*,' check bgpaths in static/background.nl'
            print*,'bgmodel=',bgmodel,' maxbgmodels=,',maxbgmodels
            print*,' lga process aborted...'
            stop
         endif


         print *,'looping through models (lga-2):',lga_status, i, 
     1                                           bg_files,reject_cnt

c
c this commented subroutine is under development to replace
c get_acceptable_files.  at the moment get_acceptable_files
c is still doing the job even though it is difficult software
c to work with. get_acceptable_files is also used in dprep.f 
c
c        call get_bkgd_files(i4time_now_lga,bgpath,bgmodel
         if (bgmodel .ne. 11) then !wni
           print*,'calling get_acceptable_files '
           print*,'bgpath:   ',trim(bgpath)
           print*,'cmodel:   ',trim(cmodel)
           print*,'bgmodel:  ',bgmodel
           print*,'lgb_only: ',lgb_only
           print*

           call get_acceptable_files(i4time_now_lga,bgpath,bgmodel
     +        ,names,max_files,use_analysis,use_forecast,bg_files
     +        ,accepted_files
     +        ,forecast_length,cmodel
     +        ,nx_bg,ny_bg,nz_bg,reject_names,reject_cnt)

           if(accepted_files.eq.0.and.bg_files.eq.0) then

             print*,'no acceptable files found for background model:'
             print*,'bgpath =  ',trim(bgpath)
             print*,'bgmodel = ',bgmodel 

             no_infinite_loops=no_infinite_loops+1
             reject_cnt=reject_cnt+1
             if(i.eq.nbgmodels)lga_status = 0 

           elseif(accepted_files.eq.0 .and. cmodel .eq. 'laps')then

             print*,'time interp previous laps analyses',
     +'and write result for current and cycle_time+1 backgrounds'
             print*,'code not yet available'
 
c             i4time_lga=i4time_now_lga+laps_cycle_time
c             call advance_analyses(i4time_now_lga,i4time_lga,nx_laps
c    +,ny_laps,nz_laps)
c             print*,'finished in advance_analyses'

           elseif(accepted_files.eq.0.and.reject_cnt.eq.bg_files)then
             write(6,*)' accepted files is 0, incrementing reject_cnt'
             reject_cnt=reject_cnt+1

           elseif(accepted_files.eq.0)then
             write(6,*)'note: accepted files is 0, skip lga_driver call'

           else

             print *, ' '
             print *, 'input parameters'
             print *, '----------------'
             print *
             print *, ' analysis setup: '
             print *, nx_laps,ny_laps,nz_laps
             print *, laps_cycle_time
             print *

             if(luse_sfc_bkgd .and. cmodel.ne.'eta48_conus'
!    1                        .and. cmodel.ne.'fmi_netcdf_ll'
     1                        .and. cmodel.ne.'laps_fua'
     1                        .and. cmodel.ne.'laps'
     1                                                       )then
                print*,'warning: ',cmodel,
     1             ' is experimental when turning on luse_sfc_bkgd'
             endif
c
c *** call lga driver and, if necessary, interpolate acceptable files.
c
             write(6,*)' names(1) = ',trim(names(1))

             call lga_driver(nx_laps,ny_laps,nz_laps,luse_sfc_bkgd,
     .          mode,
     .          laps_cycle_time,bgmodel,bgpath,cmodel,reject_cnt,
     .          reject_names,names,max_files,accepted_files,
     .          n_written,c_ftimes_written,
     .          i4time_now_lga, smooth_fields,lgb_only,lga_status)

           endif ! test for acceptable files

         else  
cc
c wrf-arw2.1 netcdf data processing
c ---------------------------------
           call lga_driver_wrfarw(nx_laps,ny_laps,nz_laps, ! add nx_laps,ny_laps,nz_laps by wei-ting(130326)
     .             bgpath,cmodel(1:12),use_analysis,forecast_length,
     .             luse_sfc_bkgd,
     .             i4time_now_lga,smooth_fields,lga_status)

         endif     

c these constructs force t-1 and t+1 cycle time background generation.
c these should be removed when lga runs outside sched.pl
c
         if(lga_status.eq.1.and.ltime(itime_inc))then
c             lga_status = -99
            ltime(itime_inc) = .false.
            if(itime_inc.lt.0)then
              print*,'completed time-minus-one cycle time interp'
            elseif(itime_inc.eq.0)then
              print*,'completed cycle time interpolation'
            else
               print*,'completed time-plus-one cycle time interp'
            endif
         endif
         print*
                                                                                                                                              
         if(lga_status.ne.1) then 
           reject_cnt = 1
           bg_files = 0
         endif

       enddo !do while

      enddo !itime_inc = -1,+1

      if(no_infinite_loops.gt.nbgmodels) then
         print*,'error: lga infinite loop condition found'
      endif
c
      print*,'number of hinterp lga/lgb files written = ',n_written
      do i = 1,n_written
          print*,c_ftimes_written(i)
      enddo ! i

      if(lga_status.eq.0) goto 965
      print *,'lga...normal completion.'
      stop
c
c *** error trapping.
c
900   continue
      print *,'error reading background model number.'
      stop
c
910   continue
      print *,'error reading laps grid information.'
      stop
c
920   continue
      print *,'error reading laps cycle time.'
      stop
c
930   continue
      print *,'error reading laps root directory path.'
      stop
c
940   continue
      print *,'error reading laps domain specification.'
      stop
c
950   continue
      print *,'error reading background model grid dimensions.'
      stop
c
960   continue
      print *,'error reading background data directory paths.'
      stop
965   continue
      print *,'no acceptable background model found.'
      
c
      end
c
