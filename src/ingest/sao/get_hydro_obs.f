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
c
        subroutine get_hydro_obs(maxobs,maxsta,i4time_sys,
     &                 path_to_local_data,local_format,
     &                 itime_before,itime_after,
     &                 itest_madis_qc,l_multiple_reports,
     &                 lat,lon,ni,nj,grid_spacing,
     &                 nn,n_obs_g,n_obs_b,stations,
     &                 reptype,atype,weather,wmoid,
     &                 store_1,!store_2,store_2ea,
!    &                 store_3,store_3ea,store_4,store_4ea,
     &                 store_5,store_5ea,store_6,store_6ea,
!    &                 store_7,store_cldht,store_cldamt,
     &                 provider, laps_cycle_time, 
     &                 local_obs_thresh, i4wait_local_obs_max, jstatus)       

c*****************************************************************************
c
c.....  input variables/arrays
c
        integer maxsta ! processed stations for lso file
        character*(*) path_to_local_data, local_format
c
c.....  local variables/arrays
c
	integer    maxobs
	integer    rtime
        integer  i4time_ob_a(maxobs), before, after
        real    lat(ni,nj), lon(ni,nj)
        real  k_to_f
        character*9 a9time_before, a9time_after, a9time_a(maxobs)
        logical l_reject(maxobs), ltest_madis_qc, ltest_madis_qcb
        logical l_multiple_reports, l_same_stn
        logical l_good_global,l_first_solar
c
	integer  wmoid(maxsta)
c
	character  save_stn(maxobs)*6
	character  timech*9, time*4
	character  stations(maxsta)*20
	character  provider(maxsta)*11
	character  presweather(maxobs)*25, weather(maxsta)*25
	character  reptype(maxsta)*6, atype(maxsta)*6
	character  store_cldamt(maxsta,5)*4 
        character*13 filename13, cvt_i4time_wfo_fname13
        character*150 data_file 
        character*40 string
c
c.....  declarations for call to netcdf reading routine (from gennet)

      include 'netcdf.inc'
      integer icchecknum, qcchecknum, maxstaticids, ninventorybins,
     +     recnum,nf_fid, nf_vid, nf_status
      parameter (icchecknum=100)       ! manually added
      parameter (qcchecknum=100)       ! manually added
      parameter (maxstaticids=30000)   ! manually added
      parameter (ninventorybins=24)    ! manually added
      integer filtersetnum, firstinbin(ninventorybins), firstoverflow,
     +     globalinventory, invtime(maxobs), inventory(maxstaticids),
     +     isoverflow(maxobs), lastinbin(ninventorybins),
     +     lastrecord(maxstaticids), nstaticids,
     +     numericwmoid(maxobs), precip12hrica(maxobs),
     +     precip12hricr(maxobs), precip12hrqca(maxobs),
     +     precip12hrqcr(maxobs), precip1hrica(maxobs),
     +     precip1hricr(maxobs), precip1hrqca(maxobs),
     +     precip1hrqcr(maxobs), precip24hrica(maxobs),
     +     precip24hricr(maxobs), precip24hrqca(maxobs),
     +     precip24hrqcr(maxobs), precip3hrica(maxobs),
     +     precip3hricr(maxobs), precip3hrqca(maxobs),
     +     precip3hrqcr(maxobs), precip5minica(maxobs),
     +     precip5minicr(maxobs), precip5minqca(maxobs),
     +     precip5minqcr(maxobs), precip6hrica(maxobs),
     +     precip6hricr(maxobs), precip6hrqca(maxobs),
     +     precip6hrqcr(maxobs), precipaccumica(maxobs),
     +     precipaccumicr(maxobs), precipaccumqca(maxobs),
     +     precipaccumqcr(maxobs), prevrecord(maxobs),
     +     secondsstage1_2(maxobs), secondsstage3(maxobs)
      real elevation(maxobs), latitude(maxobs), longitude(maxobs),
     +     precip12hr(maxobs), precip12hrqcd( qcchecknum, maxobs),
     +     precip1hr(maxobs), precip1hrqcd( qcchecknum, maxobs),
     +     precip24hr(maxobs), precip24hrqcd( qcchecknum, maxobs),
     +     precip3hr(maxobs), precip3hrqcd( qcchecknum, maxobs),
     +     precip5min(maxobs), precip5minqcd( qcchecknum, maxobs),
     +     precip6hr(maxobs), precip6hrqcd( qcchecknum, maxobs),
     +     precipaccum(maxobs), precipaccumqcd( qcchecknum, maxobs),
     +     riverflow(maxobs), riverstage(maxobs)
      double precision observationtime(maxobs), receivedtime(maxobs),
     +     reporttime(maxobs), riverreportchangetime(maxobs)
      character precip5mindd(maxobs)
      character precip12hrdd(maxobs)
      character*72 ict(icchecknum)
      character precip1hrdd(maxobs)
      character*11 dataprovider(maxobs)
      character*60 qct(qcchecknum)
      character*11 handbook5id(maxobs)
      character precip24hrdd(maxobs)
      character precipaccumdd(maxobs)
      character*51 stationname(maxobs)
      character precip3hrdd(maxobs)
      character*11 stationtype(maxobs)
      character precip6hrdd(maxobs)
      character*256 rawmessage(maxobs)
      character*24 staticids(maxstaticids)
      character*12 providerid(maxobs)
      character*4 homewfo(maxobs)
      character*11 stationid(maxobs)

      real seasurfacetemp(maxobs) ! manually added
c
c.....  output arrays.
c
	real  store_1(maxsta,4), 
     &          store_2(maxsta,3), store_2ea(maxsta,3),
     &          store_3(maxsta,4), store_3ea(maxsta,2),
     &          store_4(maxsta,5), store_4ea(maxsta,2),
     &          store_5(maxsta,4), store_5ea(maxsta,4),
     &          store_6(maxsta,5), store_6ea(maxsta,2),
     &          store_7(maxsta,3),
     &          store_cldht(maxsta,5)

        integer ibmask(8)
c
c.....  start.
c
 
        l_first_solar = .true.

        if(itest_madis_qc .gt. 0)then
            if(itest_madis_qc .eq. 15)then  ! call dd & qcr checking routines
                ltest_madis_qc  = .true.    ! for subjective qc reject list
                ltest_madis_qcb = .true.
                ibmask(1) = 0
                ibmask(2) = 1               ! validity check applied
                ibmask(3) = 0
                ibmask(4) = 0
                ibmask(5) = 0
                ibmask(6) = 1               ! statistical spatial consistency check
                ibmask(7) = 0
                ibmask(8) = 0
                level_qc = 0                ! subjective qc (reject list) only
            else                            ! values of 1-2 (dd flag check)
                ltest_madis_qc  = .true.
                ltest_madis_qcb = .false.
                level_qc = itest_madis_qc
            endif
        else                                ! value of 0 (neither check routine)
            ltest_madis_qc  = .false.
            ltest_madis_qcb = .false.
        endif

             
        write(6,*)' subroutine get_hydro_obs:' 
        write(6,*)
     1      ' itest_madis_qc/ltest_madis_qc/ltest_madis_qcb/level_qc = '   
     1       ,itest_madis_qc,ltest_madis_qc,ltest_madis_qcb,level_qc       

c
c.....	set jstatus flag for the local data to bad until we find otherwise.
c
	jstatus = -1

        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return

        call get_box_size(box_size,istatus)
        if(istatus .ne. 1)return
c
c.....  figure out the size of the "box" in gridpoints.  user defines
c.....  the 'box_size' variable in degrees, then we convert that to an
c.....  average number of gridpoints based on the grid spacing.
c
        box_length = box_size * 111.137 !km/deg lat (close enough for lon)
        ibox_points = box_length / (grid_spacing / 1000.) !in km
        box_low = 1. - float(ibox_points)    !buffer on west/south side
        box_idir = float( ni + ibox_points)  !buffer on east
        box_jdir = float( nj + ibox_points)  !buffer on north

        nn_in = nn
c
c.....	zero out the counters.
c
 10     nn = nn_in
        n_obs_g = 0	        ! # of local obs in the laps grid
        n_obs_ng = 0	        ! # of local obs not in the laps grid
        n_obs_b = 0	        ! # of local obs in the box
c
c.....  get the data from the netcdf file.  first, open the file.
c.....  if not there, return to obs_driver.
c
        ix = 1
c
c.....  set up the time window.
c
	before = i4time_sys - itime_before
	after  = i4time_sys + itime_after

!       ob times contained in each file
        i4_contains_early = 0 
        i4_contains_late = 3599

        call get_filetime_range(before,after                
     1                         ,i4_contains_early,i4_contains_late       
     1                         ,3600                                     
     1                         ,i4time_file_b,i4time_file_a)              

        i4_elapsed = ishow_timer()

        do i4time_file = i4time_file_a, i4time_file_b, -3600

            call s_len(path_to_local_data,len_path)
            filename13= cvt_i4time_wfo_fname13(i4time_file)
 	    data_file = path_to_local_data(1:len_path)//filename13

            write(6,*)' mesonet file = ',data_file(1:len_path+13)

!           goto 590 ! debugging test
	    nf_status = nf_open(data_file,nf_nowrite,nf_fid)

	    if(nf_status.ne.nf_noerr) then
	       print *, nf_strerror(nf_status)
	       print *, data_file
	       go to 590
	    endif
c
c.....  get the dimension of some of the variables.
c
c.....  "recnum"
c
	    nf_status = nf_inq_dimid(nf_fid,'recnum',nf_vid)
	    if(nf_status.ne.nf_noerr) then
	       print *, nf_strerror(nf_status)
	       print *,'dim recnum'
	    endif
	    nf_status = nf_inq_dimlen(nf_fid,nf_vid,recnum)
	    if(nf_status.ne.nf_noerr) then
	       print *, nf_strerror(nf_status)
	       print *,'dim recnum'
	    endif

            if(recnum .gt. maxobs-ix+1)then
                write(6,*)
     1              ' error: exceeded maxobs limits in get_local_obs'
     1              ,ix-1,recnum,(ix-1)+recnum,maxobs
                write(6,*)' try increasing "maxobs" in obs_driver.nl'
                go to 590
            endif

c
c.....  call the read routine.
c
        call read_madis_hydro_netcdf(nf_fid, icchecknum, qcchecknum, 
     +     maxstaticids, ninventorybins, recnum, elevation(ix), 
     +     latitude(ix), longitude(ix), precip12hr(ix), 
     +     precip12hrqcd(1,ix), precip1hr(ix), precip1hrqcd(1,ix), 
     +     precip24hr(ix), precip24hrqcd(1,ix), precip3hr(ix), 
     +     precip3hrqcd(1,ix), precip5min(ix), precip5minqcd(1,ix), 
     +     precip6hr(ix), precip6hrqcd(1,ix), precipaccum(ix), 
     +     precipaccumqcd(1,ix), riverflow(ix), riverstage(ix), ict, 
     +     qct, dataprovider(ix), handbook5id(ix), homewfo(ix), 
     +     precip12hrdd(ix), precip1hrdd(ix), precip24hrdd(ix), 
     +     precip3hrdd(ix), precip5mindd(ix), precip6hrdd(ix), 
     +     precipaccumdd(ix), providerid(ix), rawmessage(ix), 
     +     staticids, stationid(ix), stationname(ix), 
     +     stationtype(ix), observationtime(ix), receivedtime(ix), 
     +     reporttime(ix), riverreportchangetime(ix), filtersetnum, 
     +     firstinbin, firstoverflow, globalinventory, invtime(ix), 
     +     inventory, isoverflow(ix), lastinbin, lastrecord, 
     +     nstaticids, numericwmoid(ix), precip12hrica(ix), 
     +     precip12hricr(ix), precip12hrqca(ix), precip12hrqcr(ix), 
     +     precip1hrica(ix), precip1hricr(ix), precip1hrqca(ix), 
     +     precip1hrqcr(ix), precip24hrica(ix), precip24hricr(ix), 
     +     precip24hrqca(ix), precip24hrqcr(ix), precip3hrica(ix), 
     +     precip3hricr(ix), precip3hrqca(ix), precip3hrqcr(ix), 
     +     precip5minica(ix), precip5minicr(ix), precip5minqca(ix), 
     +     precip5minqcr(ix), precip6hrica(ix), precip6hricr(ix), 
     +     precip6hrqca(ix), precip6hrqcr(ix), precipaccumica(ix), 
     +     precipaccumicr(ix), precipaccumqca(ix), 
     +     precipaccumqcr(ix), prevrecord(ix), secondsstage1_2(ix), 
     +     secondsstage3(ix),badflag)

            n_local_file = recnum
            write(6,*)'     n_local_file = ',n_local_file

            ix = ix + n_local_file

            i4_elapsed = ishow_timer()

590     enddo                  ! i4time_file

        n_local_all = ix - 1
        write(6,*)' n_local_all = ',n_local_all
c
c.....  first check the data coming from the netcdf files.  there can be
c.....  "floatinf" (used as fill value) in some of the variables.  these
c.....  are not handled the same by different operating systems.  for 
c.....  example, ibm systems make "floatinf" into "nan" and store them that
c.....  way in the file, which messes up other laps routines.  this code
c.....  checks for "floatinf" and sets the variable to 'badflag'.  if the
c.....  "floatinf" is in the lat, lon, elevation, or time of observation,
c.....  we toss the whole ob since we can't be sure where it is.
c
        max_write = 100

        if(n_local_all .gt. maxobs)then
           write(6,*)' error in get_hydro_obs: n_local_all is ',
     1                                         n_local_all
           write(6,*)' try increasing obs_driver.nl/maxobs from ',maxobs
           stop
        endif
c
c..................................
c.....	first qc loop over all the obs.
c..................................
c
	do i=1,n_local_all
           l_reject(i) = .false.
c
c........  toss the ob if lat/lon/elev or observation time are bad by setting 
c........  lat to badflag (-99.9), which causes the bounds check to think that
c........  the ob is outside the laps domain.
	   if( nanf( latitude(i) ) .eq. 1 ) l_reject(i) = .true.
	   if( nanf( longitude(i) ) .eq. 1 ) l_reject(i) = .true.
	   if( nanf( elevation(i) ) .eq. 1 ) l_reject(i) = .true.
	   if( nanf( observationtime(i) ) .eq. 1 ) l_reject(i) = .true.

c
c.....  bounds check: is station in the box?  find the ob i,j location
c.....  on the laps grid, then check if outside past box boundary.
c
!          test for invalid latitude
           if(latitude(i) .lt. -90 .or. latitude(i) .gt. +90.)then
               if(.true.)then
                   write(6,81,err=82)i,n_local_all
     1                               ,stationid(i)
     1                               ,latitude(i)
 81                format(2i7,1x,a8,' invalid latitude ',e12.5)
               endif
 82            l_reject(i) = .true.
               go to 105
           endif

!          test for badflag or (close to s pole but not quite at it)
!          check can also be generalized in 'latlon_to_rlapsgrid' for 'lambert'
           if(latitude(i) .lt. -89.999 .and. latitude(i) .ne. -90.)then
               l_reject(i) = .true.
               go to 105
           endif

           call latlon_to_rlapsgrid(latitude(i),longitude(i),lat,lon,       
     &                              ni,nj,ri_loc,rj_loc,istatus)
           if(ri_loc.lt.box_low .or. ri_loc.gt.box_idir
     1   .or. rj_loc.lt.box_low .or. rj_loc.gt.box_jdir) then
               if(i .le. max_write)then
                   write(6,91,err=92)i,stationid(i)
     1                               ,nint(ri_loc),nint(rj_loc)
 91                format(i6,1x,a8,' out of box ',2i12)
               endif
 92            l_reject(i) = .true.
               go to 105
           endif
c
c.....  elevation ok?
c
	   if(elevation(i).gt.5200. .or. elevation(i).lt.-400.)then
               l_reject(i) = .true.
               go to 105
           endif

!          end of geographic location check

	   i4time_ob_a(i) = nint(observationtime(i)) + 315619200
	   call make_fnam_lp(i4time_ob_a(i),a9time_a(i),istatus)

           call filter_string(stationid(i))
c
c........  check to see if its in the desired time window.
c
	   if(i4time_ob_a(i) .lt. before 
     1   .or. i4time_ob_a(i) .gt. after) then
               if(i .le. max_write)then
                   write(6,71,err=72)i,stationid(i)
     1                               ,a9time_a(i),i4time_ob_a(i)
     1                               ,before,after
 71		   format(i6,1x,a8,' out of time ',a11,3i12)
               endif
 72            l_reject(i) = .true.
               go to 105
           endif

!          end of time check

!          pick closest station if multiple stations are in time window
           if(.not. l_multiple_reports)then
             do k = 1,i-1
               if(stationid(i) .eq. stationid(k))then ! possibly the same stn
                 l_same_stn = .true.
                 if(latitude(i)  .ne. latitude(k) .or.
     1              longitude(i) .ne. longitude(k)      )then
                     l_same_stn = .false.
                 endif

                 if(l_same_stn)then ! added to if block for efficiency
                   if( (.not. l_reject(i)) .and. (.not. l_reject(k)) 
     1                                                             )then       
                     i_diff = abs(i4time_ob_a(i) - i4time_sys)
                     k_diff = abs(i4time_ob_a(k) - i4time_sys)

                     if(i_diff .ge. k_diff)then
                         i_reject = i
                     else
                         i_reject = k
                     endif

                     write(6,51)i,k,stationid(i),a9time_a(i),a9time_a(k)       
     1                         ,i_reject
 51		     format(' duplicate detected ',2i6,1x,a6,1x,a9,1x,a9
     1                     ,1x,i6)

                     l_reject(i_reject) = .true.
                   endif ! both weren't rejected
                 endif ! same station
               endif ! possibly the same stn
             enddo ! k
           endif ! l_multiple_reports
c
c
!          if( nanf( soilmoisturepercent(i)) .eq. 1 ) 
!    1               soilmoisturepercent(i) = badflag

 105       continue
c
	enddo !i

        write(6,*)' completed 1st qc loop'

        i4_elapsed = ishow_timer()
c
c..................................
c.....	second qc loop over all the obs.
c..................................
c
	jfirst = 1
c
	do i=1,n_local_all

           if(l_reject(i))go to 125
c
c.....  right time, right location...

           timech = a9time_a(i)
	   time = timech(6:9)
	   read(time,*) rtime
c
c.....  check if station is reported more than once this
c.....  time period.
c
           if(.false.)then ! we may not need this second dupe check
	     if(jfirst .eq. 1) then
	       icount = 1
	       save_stn(1) = stationid(i)
	       jfirst = 0
	       go to 150
	     endif
c
	     do k=1,icount
               if(stationid(i) .eq. save_stn(k)) then
                 write(6,*)' rejecting duplicate ',i,k,stationid(i)
     1                    ,' ',a9time_a(i),' ',a9time_a(k)
                 go to 125
               endif
	     enddo !k
c
	     icount = icount + 1
	     save_stn(icount) = stationid(i)  ! only one...save for checking
c 
           endif ! second dupe check (set to .false.)

 150	   nn = nn + 1

           if(nn .gt. maxsta)then
              write(6,*)' error in get_hydro_obs: increase maxsta '
     1                 ,nn,maxsta
              stop
           endif
 
           n_obs_b = n_obs_b + 1     !station is in the box

           call s_len2(dataprovider(i),lenp)  
c
c.....  check if its in the laps grid.
c
           call latlon_to_rlapsgrid(latitude(i),longitude(i),lat,lon,       
     &                              ni,nj,ri_loc,rj_loc,istatus)
           if(  (ri_loc.lt.1. .or. ri_loc.gt.float(ni)) .or.
     1          (rj_loc.lt.1. .or. rj_loc.gt.float(nj))     )then
               n_obs_ng = n_obs_ng + 1 ! outside the grid (inside the box)
           else
               n_obs_g = n_obs_g + 1   ! on grid...count it
           endif
c
c
c.....  convert units for storage.  for those variables with a "change
c.....  time", check to make sure the variable was observed within the
c.....  last cycle (and that they're not just carrying an old ob for the 
c.....  current time).
c
c
c
c..... precip 1hr
         pcp1 = badflag
         if(.true.                            .and.
     1      precip1hr(i) .ge. 0.              .and.    
     1      precip1hr(i) .lt. 10000.          .and.    
     1      precip1hr(i) .ne. badflag               )then
             pcp1 = precip1hr(i) / 25.4 ! convert mm to inches
             write(6,*)' found a hydro 1hr precip ob: '
     1                 ,pcp1,' ',dataprovider(i)(1:lenp),' '
     1                 ,precip1hrdd(i)
         endif
c
c..... precip 3hr
         pcp3 = badflag
         if(.true.                            .and.
     1      precip3hr(i) .ge. 0.              .and.    
     1      precip3hr(i) .lt. 10000.          .and.    
     1      precip3hr(i) .ne. badflag               )then
             pcp3 = precip3hr(i) / 25.4 ! convert mm to inches
             write(6,*)' found a hydro 3hr precip ob: '
     1                 ,pcp3,' ',dataprovider(i)(1:lenp),' '
     1                 ,precip3hrdd(i)
         endif
c
c..... precip 6hr
         pcp6 = badflag
         if(.true.                            .and.
     1      precip6hr(i) .ge. 0.              .and.    
     1      precip6hr(i) .lt. 10000.          .and.    
     1      precip6hr(i) .ne. badflag               )then
             pcp6 = precip6hr(i) / 25.4 ! convert mm to inches
             write(6,*)' found a hydro 6hr precip ob: '
     1                 ,pcp6,' ',dataprovider(i)(1:lenp),' '
     1                 ,precip6hrdd(i)
         endif
c
c..... precip 24hr
         pcp24 = badflag
         if(.true.                            .and.
     1      precip24hr(i) .ge. 0.             .and.    
     1      precip24hr(i) .lt. 10000.         .and.    
     1      precip24hr(i) .ne. badflag               )then
             pcp24 = precip24hr(i) / 25.4 ! convert mm to inches
             write(6,*)' found a hydro 24hr precip ob: '
     1                 ,pcp24,' ',dataprovider(i)(1:lenp),' '
     1                 ,precip24hrdd(i)
         endif
c
c..... other stuff.  
c
c
	 store_6ea(nn,1) = 0.0             ! precipitation (in)
	 store_6ea(nn,2) = 0.0             ! snow cover (in) 
c
c
c..... output the data to the storage arrays
c
	 call s_len(stationid(i), len)
         if(len .ne. 0)then
             stations(nn)(1:len) = stationid(i)(1:len) ! station name
         else
             write(6,*)' warning in get_local_obs: blank station name.'
     1                ,' assigning name ',i
             write(stations(nn),101)i
 101	     format(i5,15x)
         endif
c
         if(lenp .ne. 0) then
	     provider(nn)(1:lenp) = dataprovider(i)(1:lenp)    ! data provider
         endif
         call filter_string(provider(nn))
c
         call s_len(stationtype(i), len)
         if(len .ne. 0) then
            ilen = min(len, 6)
            atype(nn)(1:ilen) = stationtype(i)(1:ilen) ! auto stn type
         endif
         call filter_string(atype(nn))
c
	 reptype(nn)(1:6) = 'ldad  '            ! report type
	 wmoid(nn) = ibadflag                   ! wmo id
c
	 store_1(nn,1) = latitude(i)            ! station latitude
	 store_1(nn,2) = longitude(i)           ! station longitude
	 store_1(nn,3) = elevation(i)           ! station elevation
	 store_1(nn,4) = rtime                  ! observation time
c
         store_6(nn,1) = pcp1                   ! 1-h precipitation
         store_6(nn,2) = pcp3                   ! 3-h precipitation
         store_6(nn,3) = pcp6                   ! 6-h precipitation
         store_6(nn,4) = pcp24                  ! 24-h precipitation
         store_6(nn,5) = badflag                ! snow cover
c
 125     continue
       enddo !i

       i4_elapsed = ishow_timer()

       i4t_since_sys = i4time_now_gg() - i4time_sys

       write(6,*)' i4t_since_sys / nobs = ',i4t_since_sys,n_obs_b

       if(n_obs_b       .lt. local_obs_thresh .and. 
     1    i4t_since_sys .lt. i4wait_local_obs_max       )then
           write(6,*)' waiting 60 sec for more obs'
           call snooze_gg(60.,istatus)
           go to 10
       endif
c
c..... that's it...lets go home.
c
       print *,' found ',n_obs_b,' hydro obs in the laps box'
       print *,' found ',n_obs_g,' hydro obs in the laps grid'
       print *,'       ',n_obs_ng,' hydro obs outside the laps grid'
       print *,' '
       jstatus = 1            ! everything's ok...
       return
c
 990   continue               ! no data available
       jstatus = 0
       print *,' warning: no data available from get_hydro_obs'
       return
c
       end

