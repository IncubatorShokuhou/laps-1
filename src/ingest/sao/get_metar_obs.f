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
c
	subroutine get_metar_obs(maxobs,maxsta,i4time_sys,
     &                      path_to_metar,metar_format,
     &                      minutes_to_wait_for_metars,
     &                      ick_metar_time,itime_before,itime_after,
     &                      eastg,westg,anorthg,southg,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_sao_g,n_sao_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, jstatus)

c
c*****************************************************************************
c
cdoc    routine to gather metar (and synop) data for laps.   
c
c	changes:
c		p. stamus  10-30-96  original version (from get_sao_obs).
c                          11-13-96  pass in full path to data.
c                          12-06-96  put atype (01 or 02) in obstype var.
c	                   06-09-97  changes for new metar cdl.
c                          01-29-98  upgrades for ls2.
c                          05-01-98  add soil moisture variables.
c                          06-21-99  change ob location check to gridpt space.
c                                      figure box size in gridpoint space from
c                                      user-defined size (deg) and grid_spacing.
c                          07-22-99  add elev ck for badflag.  fix wmoid var.
c
c
c*****************************************************************************
c
	include 'netcdf.inc'
c
c.....  read arrays.
c
        integer maxobs ! raw data file
        integer maxsta ! output lso file
	real*8  timeobs(maxobs)
	real  lats(maxobs), lons(maxobs), elev(maxobs)
	real  t(maxobs), td(maxobs), tt(maxobs), ttd(maxobs)
	real  dd(maxobs), ddg(maxobs), ff(maxobs), ffg(maxobs)
	real  stnp(maxobs), mslp(maxobs), alt(maxobs)
	real  ht(6,maxobs), vis(maxobs), dp(maxobs)
        real  rh(maxobs), sr(maxobs), st(maxobs)
	real  pcp1(maxobs), pcp3(maxobs), pcp6(maxobs), pcp24(maxobs)
	real  max24t(maxobs), min24t(maxobs), snowcvr(maxobs)
        real    lat(ni,nj), lon(ni,nj), k_to_f
c
	real  store_1(maxsta,4), 
     &          store_2(maxsta,3), store_2ea(maxsta,3),
     &          store_3(maxsta,4), store_3ea(maxsta,2),
     &          store_4(maxsta,5), store_4ea(maxsta,2),
     &          store_5(maxsta,4), store_5ea(maxsta,4),
     &          store_6(maxsta,5), store_6ea(maxsta,2),
     &          store_7(maxsta,3),
     &          store_cldht(maxsta,5)
c
	integer  i4time_ob, wmoid(maxobs), wmoid_in(maxobs)
	integer  i4time_before, i4time_after
	integer    rtime, dpchar(maxobs)
	integer    maxskycover, recnum, nf_fid, nf_vid, nf_status
c
	character  stname(maxobs)*5, save_stn(maxobs)*5
	character  data_file*150, timech*9, time*4, metar_format*(*)
	character  filename13*13, fname9_to_wfo_fname13*13
        character  path_to_metar*(*), a9time*9, a8time*8, a9_to_a8*8
        character  path_to_local_cwb*150
	character  cvr(6,maxobs)*8
	character  stations(maxsta)*20, provider(maxsta)*11
     1                                , c11_provider*11
	character  weather(maxobs)*25, wx(maxobs)*25
	character  reptype(maxobs)*6, atype(maxobs)*6
	character  reptype_in(maxobs)*6, atype_in(maxobs)*6
	character  store_cldamt(maxsta,5)*4
c
        integer cnt
	logical exists, l_parse
        logical l_dupe_time(maxobs)
        data exists/.false./
        data cnt/0/
c
c
c.....	set jstatus flag for the sao data to bad until we find otherwise.
c
	jstatus = -1

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return

        call get_box_size(box_size,istatus)
        if(istatus .ne. 1)return

        l_dupe_time = .false.

!       initialize variables mainly used for cwb combined synop/local data
        ddg  = badflag
        stnp = badflag
        rh   = badflag
        sr   = badflag
        st   = badflag
c
c.....  figure out the size of the "box" in gridpoints.  user defines
c.....  the 'box_size' variable in degrees, then we convert that to an
c.....  average number of gridpoints based on the grid spacing.
c
        box_length = box_size * 111.137 !km/deg lat (close enough for lon)
        ibox_points = box_length / (grid_spacing / 1000.) !in km
c
c.....	zero out the counters.
c
	n_sao_g = 0		! # of saos in the laps grid
	n_sao_b = 0		! # of saos in the box
c
c.....  set up the time window.
c
	i4time_before = i4time_sys - itime_before
	i4time_after  = i4time_sys + itime_after

        call s_len(metar_format,len_metar_format)

        if(     l_parse(metar_format,'fsl') 
     1     .or. l_parse(metar_format,'nimbus')
     1     .or. l_parse(metar_format,'wfo')     
     1     .or. l_parse(metar_format,'madis')     
     1                                             )then ! netcdf format

            ix = 1

            i4_contains_early = 900
            i4_contains_late = 2699
            i4_file_interval = 3600

            call get_filetime_range(i4time_before,i4time_after                
     1                             ,i4_contains_early,i4_contains_late       
     1                             ,i4_file_interval                         
     1                             ,i4time_file_b,i4time_file_a)              

            call s_len(path_to_metar,len_path)

            n_sao_all = 0

            do i4time_file = i4time_file_b, i4time_file_a
     1                     , i4_file_interval

                call make_fnam_lp(i4time_file,a9time,istatus)

                if(metar_format(1:len_metar_format) .eq. 'nimbus')then
	            data_file = path_to_metar(1:len_path)
     1                            //a9time// '0100o'

                elseif(metar_format(1:len_metar_format) .eq. 'wfo'
     1            .or. metar_format(1:len_metar_format) .eq. 'madis'
     1                                                           )then      
                    filename13=fname9_to_wfo_fname13(a9time)       
                    data_file = path_to_metar(1:len_path) // filename13       

                else
                    write(6,*)' error: unknown metar format '
     1                       ,metar_format          
                    istatus = 0
                    return

                endif

                do while(.not. exists 
     &                           .and. 
     &                    cnt .le. minutes_to_wait_for_metars
     &                           .and.
     &                    i4time_file .eq. i4time_file_a        
     &                           .and.
     &                    i4time_file .le. i4time_now_gg()        )
c
	            inquire(file=data_file,exist=exists)
                    if(.not. exists) then
                        if(cnt .lt. minutes_to_wait_for_metars)then
                            print*,'waiting for file ', data_file
                            call waiting_c(60)
                        endif
                        cnt = cnt+1               
	            endif

	        enddo ! while in waiting loop

c
c.....          get the data from the netcdf file.  first, open the file.
c 
	        nf_status = nf_open(data_file,nf_nowrite,nf_fid)

	        if(nf_status.ne.nf_noerr) then ! no file found to open
	           print *, nf_strerror(nf_status)
	           print *, data_file
                   write(6,*)' warning: not found in get_metar_obs - '
     1                      ,data_file       
                   n_metar_file = 0
                   goto580
                else
                   write(6,*)' file found - ',data_file
	        endif
c
c.....          get the dimension of some of the variables.
c.....          "maxskycover"
c
	        nf_status = nf_inq_dimid(nf_fid,'maxskycover',nf_vid)
	        if(nf_status.ne.nf_noerr) then
	           print *, nf_strerror(nf_status)
	           print *,'dim maxskycover'
	        endif
	        nf_status = nf_inq_dimlen(nf_fid,nf_vid,maxskycover)
	        if(nf_status.ne.nf_noerr) then
	           print *, nf_strerror(nf_status)
	           print *,'dim maxskycover'
	        endif
c
c.....          "recnum"
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
     1              ' error: exceeded maxobs limits in get_metar_obs'
                    go to 590
                endif
c
c.....          call the read routine.
c
	        call read_metar(nf_fid , maxskycover, recnum, alt(ix),       
     &             atype_in(ix), td(ix), ttd(ix), elev(ix),
     &             lats(ix), lons(ix), max24t(ix), min24t(ix),
     &             pcp1(ix), pcp24(ix), pcp3(ix), pcp6(ix),
     &             wx(ix), dp(ix), dpchar(ix),
     &             reptype_in(ix), mslp(ix), cvr(1,ix), ht(1,ix),
     &             snowcvr(ix), stname(ix), tt(ix), t(ix),
     &             timeobs(ix), vis(ix), dd(ix), ffg(ix), ff(ix),
     &             wmoid_in(ix), badflag, istatus)

                if(istatus .ne. 1)then
                    write(6,*)
     1              '     warning: bad status return from read_metar'       
                    n_metar_file = 0
                else
                    n_metar_file = recnum
                    write(6,*)'     n_metar_file = ',n_metar_file
                endif

 580            ix = ix + n_metar_file

 590	    enddo             ! i4time_file (in time window)

            n_sao_all = ix - 1
            i4time_offset = 315619200
            c11_provider = 'nws        '

        else ! read cwb metar and synop obs
            i4time_file = (i4time_sys/3600) * 3600
            call make_fnam_lp(i4time_file,a9time,istatus)
            a8time = a9_to_a8(a9time(1:9))

            call s_len(path_to_metar,len_path)

!           read metar obs
            maxskycover=6
            recnum = maxobs

	    data_file = 
     1            path_to_metar(1:len_path)//'metar'//a8time//'.dat'

            call s_len(data_file,len_file)
            write(6,*)' cwb metar data: ',data_file(1:len_file)

            call read_metar_cwb(data_file , maxskycover, recnum, alt,    
     &         atype_in, td, ttd, elev,
     &         lats, lons, max24t, min24t,
     &         pcp1, pcp24, pcp3, pcp6,
     &         wx, dp, dpchar,
     &         reptype_in, mslp, cvr, ht,
     &         snowcvr, stname, tt, t,
     &         timeobs, vis, dd, ffg, ff,
     &         wmoid_in, badflag, n_metar_cwb, istatus)

            write(6,*)' n_metar_cwb = ',n_metar_cwb

            ix = n_metar_cwb + 1

!           read synop obs
            i4time_file = (i4time_sys/3600) * 3600
            call make_fnam_lp(i4time_file,a9time,istatus)
            a8time = a9_to_a8(a9time(1:9))

            maxskycover=6
            recnum = maxobs-ix+1

	    data_file = 
     1          path_to_metar(1:len_path)//'synop'//a8time//'.dat'

            path_to_local_cwb = path_to_metar(1:len_path)//'../loc/'

            call s_len(data_file,len_file)
            write(6,*)' cwb synop/local (via read_synop_cwb): '
     1                ,data_file(1:len_file)       

            n_synop_cwb = 0

            call read_synop_cwb(data_file , maxskycover, recnum, 
     &         i4time_sys, path_to_local_cwb,
     &         alt(ix), atype_in(ix), td(ix), ttd(ix), elev(ix),
     &         lats(ix), lons(ix), max24t(ix), min24t(ix),
     &         pcp1(ix), pcp24(ix), pcp3(ix), pcp6(ix),
     &         wx(ix), stnp(ix), dp(ix), dpchar(ix),
     &         reptype_in(ix), rh(ix),
     &         mslp(ix), cvr(1,ix), ht(1,ix),
     &         snowcvr(ix), sr(ix), st(ix), 
     &         stname(ix), tt(ix), t(ix), timeobs(ix), vis(ix), 
     &         dd(ix), ddg(ix), ff(ix), ffg(ix), wmoid_in(ix), badflag,       
     &         n_synop_cwb, istatus)
            if(istatus .ne. 1)then
                n_synop_cwb = 0
            endif

            write(6,*)' n_synop_cwb = ',n_synop_cwb

            n_sao_all = n_metar_cwb + n_synop_cwb
            if(n_sao_all .le. 0) go to 990

            i4time_offset = 0

            c11_provider = 'cwb        '

        endif
c
        write(6,*)' n_sao_all = ',n_sao_all

	if(n_sao_all .gt. maxobs)then
            write(6,*)' error, n_sao_all > maxobs ',n_sao_all,maxobs
            return
        endif
c
c.....  first check the data coming from the netcdf file.  there can be
c.....  "floatinf" (used as fill value) in some of the variables.  these
c.....  are not handled the same by different operating systems.  for 
c.....  example, ibm systems make "floatinf" into "nan" and store them that
c.....  way in the lso, which messes up other laps routines.  this code
c.....  checks for "floatinf" and sets the variable to 'badflag'.  if the
c.....  "floatinf" is in the lat, lon, elevation, or time of observation,
c.....  we toss the whole ob since we can't be sure where it is.
c
	do i=1,n_sao_all
c
c.....  toss the ob if lat/lon/elev or observation time are bad by setting 
c.....  lat to badflag (-99.9), which causes the bounds check to think that
c.....  the ob is outside the laps domain.
c
	   if( nanf( lats(i) ) .eq. 1 ) lats(i)  = badflag
	   if( nanf( lons(i) ) .eq. 1 ) lats(i)  = badflag
	   if( nanf( elev(i) ) .eq. 1 ) lats(i)  = badflag
c
	   if( nanf( timeobs(i) ) .eq. 1 ) lats(i) = badflag
c
	   if( nanf( vis(i)  ) .eq. 1 ) vis(i)   = badflag
	   if( nanf( mslp(i) ) .eq. 1 ) mslp(i)  = badflag
	   if( nanf( t(i)    ) .eq. 1 ) t(i)     = badflag
	   if( nanf( td(i)   ) .eq. 1 ) td(i)    = badflag
	   if( nanf( tt(i)   ) .eq. 1 ) tt(i)    = badflag
	   if( nanf( ttd(i)  ) .eq. 1 ) ttd(i)   = badflag
	   if( nanf( dd(i)   ) .eq. 1 ) dd(i)    = badflag
	   if( nanf( ff(i)   ) .eq. 1 ) ff(i)    = badflag
	   if( nanf( ffg(i)  ) .eq. 1 ) ffg(i)   = badflag
	   if( nanf( alt(i)  ) .eq. 1 ) alt(i)   = badflag
	   if( nanf( pcp1(i) ) .eq. 1 ) pcp1(i)  = badflag
	   if( nanf( pcp3(i) ) .eq. 1 ) pcp3(i)  = badflag
	   if( nanf( pcp6(i) ) .eq. 1 ) pcp6(i)  = badflag
	   if( nanf( pcp24(i)) .eq. 1 ) pcp24(i) = badflag
	   if( nanf( dp(i)   ) .eq. 1 ) dp(i)    = badflag
c
	   if( nanf( snowcvr(i) ) .eq. 1 ) snowcvr(i) = badflag
	   if( nanf( max24t(i)  ) .eq. 1 ) max24t(i)  = badflag
	   if( nanf( min24t(i)  ) .eq. 1 ) min24t(i)  = badflag
c
	   do j=1,5
	      if( nanf( ht(j,i) ) .eq. 1 ) ht(j,i) = badflag
	   enddo !j
c
	enddo !i
c
c..................................
c.....	now loop over all the obs.
c..................................
c
	jfirst = 1
        box_low = 1. - float(ibox_points)    !buffer on west/south side
        box_idir = float( ni + ibox_points)  !buffer on east
        box_jdir = float( nj + ibox_points)  !buffer on north
c 
	do 125 i=1,n_sao_all
c
c.....  bounds check: is station in the box?  find the ob i,j location
c.....  on the laps grid, then check if outside past box boundary.
c
           if(lats(i) .lt. -90.) go to 125   !badflag (-99.9)...from nan ck
           call latlon_to_rlapsgrid(lats(i),lons(i),lat,lon,ni,nj,
     &                              ri_loc,rj_loc,istatus)
           if(ri_loc.lt.box_low .or. ri_loc.gt.box_idir) go to 125
           if(rj_loc.lt.box_low .or. rj_loc.gt.box_jdir) go to 125
c
c.....  elevation ok?
c
          if(elev(i) .eq. badflag) go to 125 !from read_metar (missing elev)
          if(elev(i).gt.5200. .or. elev(i).lt.-400.) go to 125
c
c.....  check to see if its in the desired time window (if the flag
c.....  says to check this).
c
	  i4time_ob = nint(timeobs(i)) + i4time_offset
c
	  if(ick_metar_time .eq. 1) then
	    if(i4time_ob.lt.i4time_before 
     1    .or. i4time_ob.gt.i4time_after) go to 125
	  endif
c
c.....  right time, right location...

 	  call make_fnam_lp(i4time_ob,timech,istatus)
	  time = timech(6:9)
	  read(time,*) rtime

!         write(6,*)'debug ',i,stname(i),timech
c
c.....  check if station is reported more than once this
c.....  time period.
c
          istn_reject = 0
          istn_keep = i

	  if(jfirst .eq. 1) then ! initialize with first station in dataset
	     isaved_sta = 1
	     save_stn(1) = stname(i)
	     jfirst = 0

             nn = nn + 1
             np = nn 

          else ! all subsequent stations
             if(isaved_sta .eq. i)then ! this violates assumed logic of code
                 write(6,*)'error isaved_sta = i',i
                 stop
             endif

!            do duplication check unless station names are 'unk'
	     do k=1,isaved_sta
	        if(stname(i) .eq. save_stn(k))then
                   if(stname(i)(1:3) .ne. 'unk')then

!                      i4time_ob_k = nint(timeobs(k)) + i4time_offset
                       int_obtime = store_1(k,4)
                       call get_sfc_obtime(int_obtime,i4time_sys
     1                                    ,i4time_ob_k,istatus)
                       if(istatus .ne. 1)then
                           write(6,*)
     1                            ' error returned from get_sfc_obtime'       
                           return
                       endif

!                      alternatively set l_dupe_time based on abs(obstime-systime)
                       i_diff = abs(i4time_ob   - i4time_sys)
                       k_diff = abs(i4time_ob_k - i4time_sys)

                       if(i_diff .ge. k_diff)then
                          istn_reject = i ! current metar isn't closer to systime
                          istn_keep   = k 
                       else
                          istn_reject = k ! current metar is closer to systime
                          istn_keep   = i 
                       endif

                       call filter_string(stname(i))

                       write(6,*)' dupe metar at ',stname(i),timech,i,k
     1                          ,i_diff,k_diff,istn_keep,istn_reject

                       l_dupe_time(i) = .true.

                   endif
                endif
	     enddo ! k
c
             if(l_dupe_time(i))then
                 if(istn_keep .eq. i)then 
                    np = istn_reject
                    write(6,*)' backfill metar at closer time '
     1                        ,stname(i),i,istn_keep,np

                 else ! skip processing of duplicate station
                    go to 125

                 endif

             else ! not flagged as dupe
                 isaved_sta = isaved_sta + 1
	         save_stn(isaved_sta) = stname(i) ! only one...save for checking

                 nn = nn + 1
                 np = nn 
!                write(6,*)' advancing isaved_sta/np to ',isaved_sta,np       

             endif

          endif ! first station

          if(np .gt. maxsta)then
              write(6,*)' error in get_metar_obs: increase maxsta '
     1                 ,np,maxsta
              stop
          endif

          call init_station(np
     1                      ,stations,provider,weather,reptype,atype      
     1                      ,store_1,store_2,store_3,store_4,store_5
     1                      ,store_6,store_7
     1                      ,store_2ea,store_3ea,store_4ea,store_5ea
     1                      ,store_6ea,dpchar,wmoid
     1                      ,store_cldht,store_cldamt,maxsta,badflag)

	  n_sao_b = n_sao_b + 1	!station is in the box
c
c.....  check if its in the laps grid.
c
          if(ri_loc.lt.1. .or. ri_loc.gt.float(ni)) go to 151  !off grid
          if(rj_loc.lt.1. .or. rj_loc.gt.float(nj)) go to 151  !off grid
	  n_sao_g = n_sao_g + 1  !on grid...count it
 151	  continue
c
c.....	figure out the cloud data.
c
	  k_layers = 0               ! number of cloud layers
c
!         this requires that the first array element has a cloud layer in it
	  if(cvr(1,i)(1:1) .eq. ' ' .and. .false.) then
	     k_layers = 0

!            test section to help determine whether this if block is needed
             k_layers1 = 0
	     do k=1,5
		if(cvr(k,i)(1:1) .ne. ' ') k_layers1 = k_layers1 + 1       
	     enddo !k

             if(k_layers1 .ne. 0)then
                 write(6,*)' warning: k_layers1 is > 0',k_layers1
             endif

	  else
	     do k=1,5
                call s_len(cvr(k,i),lenc)   
!               if(cvr(k,i)(1:1) .ne. ' ' .and.        ! valid layer coverage
		if(lenc .gt. 0            .and.        ! valid layer coverage
     1             cvr(k,i)(1:2) .ne. '-9'       )then ! cwb error check
                    k_layers = k_layers + 1       
                endif
	     enddo !k
	  endif
c
	  if(k_layers .eq. 0) then   ! no cloud data...probably amos station
	    go to 126		     ! skip rest of cloud stuff
	  endif

          iihigh = k_layers

	  do ii=1,iihigh
	    if(ht(ii,i) .gt. 25000.0) ht(ii,i) = badflag
c
	    if(cvr(ii,i)(1:3) .eq. 'skc') then
	       ht(ii,i) = 22500.0    ! manual ob

	    elseif(cvr(ii,i)(1:3) .eq. 'clr') then
	       ht(ii,i) = 3657.4     ! automatic ob

            elseif(ht(ii,i) .gt. 17000.0    .or.
     1             ht(ii,i) .eq. badflag        ) then ! check for bad height
               write(6,*)' warning in get_metar_obs: '      
     1                   ,' reject cloud ob, height = '
     1                   ,ht(ii,i),stname(i),k_layers
               ht(ii,i) = badflag
               k_layers = 0

            elseif(ii .ge. 2)then              ! check for out of order heights
               if(ht(ii,i) .lt. ht(ii-1,i))then 
                   write(6,*)' warning in get_metar_obs: '      
     1                      ,' reject cloud ob, out of order heights = '       
     1                      ,ht(ii-1,i),ht(ii,i),stname(i)
                   ht(ii,i) = badflag
                   k_layers = 0
               endif

	    endif

	  enddo !ii

 126	  continue
c
c.....	check cloud info for very high heights...set to max if greater.
c.....	also convert agl cloud heights to msl by adding elevation.
c
	if(k_layers .gt. 0) then
	  do ii=1,k_layers
	    if(ht(ii,i) .ge. 0.) then
	      ht(ii,i) = ht(ii,i) + elev(i)		! conv agl to msl
	      if(ht(ii,i) .gt. 22500.) ht(ii,i) = 22500.
	    endif
	  enddo !ii
	endif
c
c
c.....  convert units for storage. some qc checking.
c
c.....  temperature and dewpoint
c
	temp_k = tt(i)                        !set to temp_from_tenths
	if(temp_k .eq. badflag) temp_k = t(i) !no temp_from_tenths, set to t
	if(temp_k .eq. badflag) then          !t bad too?
	   temp_f = badflag                   !          bag
	else
           temp_f = k_to_f(temp_k)
	endif
        call sfc_climo_qc_r('t_f',temp_f)
c
	dewp_k = ttd(i)                        !set to dp_from_tenths
	if(dewp_k .eq. badflag) dewp_k = td(i) !no dp_from_tenths, set to dp
	if(dewp_k .eq. badflag) then           !dp bad too?
	   dewp_f = badflag                    !         bag
	else
	   dewp_f = k_to_f(dewp_k)
	endif
c
c..... wind speed and direction
c
        call sfc_climo_qc_r('dir_deg',dd(i))
        call sfc_climo_qc_r('spd_ms',ff(i))
	if(ff(i)  .ne. badflag) ff(i)  = 1.94254 * ff(i)   !m/s to kt

        call sfc_climo_qc_r('dir_deg',ddg(i))
        call sfc_climo_qc_r('spd_ms',ffg(i))
	if(ffg(i) .ne. badflag) then
	   ffg(i) = 1.94254 * ffg(i) !m/s to kt
	   if(ddg(i) .eq. badflag)ddg(i) = dd(i)
	endif
c
c..... pressure...msl and altimeter, 3-h pressure change
c
	if(alt(i)  .ne. badflag)  alt(i) =  alt(i) * 0.01   !pa to mb
        call sfc_climo_qc_r('alt_mb',alt(i))

	if(mslp(i) .ne. badflag) mslp(i) = mslp(i) * 0.01   !pa to mb
        call sfc_climo_qc_r('mslp_mb',mslp(i))

	if(stnp(i) .ne. badflag) stnp(i) = stnp(i) * 0.01   !pa to mb
        call sfc_climo_qc_r('stnp_mb',stnp(i))

	if(dp(i)   .ne. badflag)   dp(i) =   dp(i) * 0.01   !pa to mb
c
c..... visibility
c
        if(vis(i) .lt. 0.0 .or. vis(i) .ge. 1e6) then ! qc check
           vis(i) = badflag
        endif
	if(vis(i) .ne. badflag) then
	   vis(i) = vis(i) * .001      !m to km
	   vis(i) = 0.621371 * vis(i)  !km to miles
	endif
c
c..... climo-type stuff...precip and snow cover, max/min temps
c
	if(pcp1(i)  .ne. badflag)  pcp1(i) =  pcp1(i) * 39.370079 ! m to in
	if(pcp3(i)  .ne. badflag)  pcp3(i) =  pcp3(i) * 39.370079 ! m to in
	if(pcp6(i)  .ne. badflag)  pcp6(i) =  pcp6(i) * 39.370079 ! m to in
	if(pcp24(i) .ne. badflag) pcp24(i) = pcp24(i) * 39.370079 ! m to in

        if(pcp1(i)  .gt. 50.)pcp1(i) = badflag
        if(pcp3(i)  .gt. 50.)pcp3(i) = badflag
        if(pcp6(i)  .gt. 50.)pcp6(i) = badflag
        if(pcp24(i) .gt. 50.)pcp24(i) = badflag

	if(snowcvr(i) .ne. badflag) snowcvr(i) = snowcvr(i) * 39.370079 ! m to in
	if(max24t(i) .ne. badflag) then
	   max24t(i) = k_to_f(max24t(i))
	endif
	if(min24t(i) .ne. badflag) then
	   min24t(i) = k_to_f(min24t(i))
	endif
c
c
c..... fill the expected accuracy arrays.  values are based on information
c..... in the 'federal meteorological handbook no. 1' for the metars, 
c..... appendix c (http://www.nws.noaa.gov/oso/oso1/oso12/fmh1/fmh1appc.htm)
c..... note that we convert the units in appendix c to match what we're 
c..... using here.
c
c..... temperature (deg f)
c
	fon = 9. / 5.  !ratio when converting c to f
	store_2ea(np,1) = 5.0 * fon        ! start...we don't know what we have
	if(temp_f .ne. badflag) then
	   if(temp_f.ge.c2f(-62.) .and. temp_f.le.c2f(-50.)) then
	      store_2ea(np,1) = 1.1 * fon  ! conv to deg f
	   elseif(temp_f.gt.c2f(-50.) .and. temp_f.lt.c2f(50.)) then
	      store_2ea(np,1) = 0.6 * fon  ! conv to deg f
	   elseif(temp_f.ge.c2f(50.) .and. temp_f.le.c2f(54.)) then
	      store_2ea(np,1) = 1.1 * fon  ! conv to deg f
	   endif
	endif
c
c..... dew point (deg f).  also estimate a rh accuracy based on the dew point.
c..... estimates for the rh expected accuracy are from playing around with the
c..... psychrometric tables for various t/td combinations (including their
c..... accuracies from the fmh-1 appendix c).
c
	 store_2ea(np,2) = 5.0 * fon       ! start...don't know what we have 
	 store_2ea(np,3) = 50.0            ! relative humidity %
	 if(dewp_f .ne. badflag) then
	    if(dewp_f.ge.c2f(-34.) .and. dewp_f.lt.c2f(-24.)) then
	       store_2ea(np,2) = 2.2 * fon ! conv to deg f
	       store_2ea(np,3) = 20.0      ! rh (%) 
	    elseif(dewp_f.ge.c2f(-24.) .and. dewp_f.lt.c2f(-1.)) then
	       store_2ea(np,2) = 1.7 * fon ! conv to deg f
	       store_2ea(np,3) = 12.0      ! rh (%) 
	    elseif(dewp_f.ge.c2f(-1.) .and. dewp_f.le.c2f(30.)) then
	       store_2ea(np,2) = 1.1 * fon ! conv to deg f
	       store_2ea(np,3) = 8.0       ! rh (%) 
	    endif
	 endif
c
c..... wind direction (deg) and speed (kts)
c
	 store_3ea(np,1) = 10.0    ! deg 
	 store_3ea(np,2) =  1.0    ! kt
	 if(ff(i) .ne. badflag) then
	    if(ff(i).ge.1.0 .and. ff(i).le.10.0) then
	       store_3ea(np,2) = 1.0          ! kt
	    elseif((ff(i) .gt. 10.0).and.(ff(i) .lt. 200)) then
	       store_3ea(np,2) = ff(i) * 0.1  ! 10% of speed (kts)
	    endif
c
	    if(ff(i) .ge. 5.0) then    ! dir check
	       store_3ea(np,1) = 5.0   ! deg
	    endif
	 endif
c
c..... pressure and altimeter (mb)
c
	 store_4ea(np,1) = 0.68            ! pressure (mb)
	 store_4ea(np,2) = 0.68            ! altimeter (mb)
c
c..... visibility (miles).  for automated stations use a guess based 
c..... on table c-2 in appendix c of fmh-1.  for manual stations, use
c..... a guess based on the range between reportable values (e.g., for
c..... reported visibility between 0 and 3/8th mile, set accuracy to 
c..... 1/16th mile).  this isn't ideal, but its a start.
c
	 store_5ea(np,1) = 10.00         ! start with this (miles)
	 if(vis(i) .ne. badflag) then
	    if(atype_in(i)(1:2) .eq. 'a0') then   ! have an auto station
	       if(vis(i) .lt. 2.0) then
		  store_5ea(np,1) = 0.25         ! miles
	       elseif(vis(i).ge.2.0 .and. vis(i).lt.3.0) then
		  store_5ea(np,1) = 0.50         ! miles
	       elseif(vis(i) .gt. 3.0) then
		  store_5ea(np,1) = 1.00         ! miles
	       endif
	    else		! have a manual station
	       if(vis(i) .le. 0.375) then
		  store_5ea(np,1) = 0.0625       ! miles
	       elseif(vis(i).gt.0.375 .and. vis(i).le.2.0) then
		  store_5ea(np,1) = 0.125        ! miles
	       elseif(vis(i).gt.2.0 .and. vis(i).le.3.0) then
		  store_5ea(np,1) = 0.25         ! miles
	       elseif(vis(i).gt.3.0 .and. vis(i).le.15.0) then
		  store_5ea(np,1) = 1.00         ! miles
	       elseif(vis(i) .gt. 15.0) then
		  store_5ea(np,1) = 5.00         ! miles
	       endif
	    endif
	 endif
c
c..... other stuff.  don't really know about the precip, but probably
c..... worse that this guess.
c
	 store_5ea(np,2) = 0.0             ! solar radiation 
	 store_5ea(np,3) = 0.0             ! soil/water temperature
	 store_5ea(np,4) = 0.0             ! soil moisture 
c
	 store_6ea(np,1) = 0.01            ! precipitation (in)
	 store_6ea(np,2) = 1.0             ! snow cover (in) 
c
c
c..... output the data to the storage arrays
c
!        write(6,*)' saving station ',i,np,stname(i),rtime

	 call s_len(stname(i), len)
	 stations(np)(1:len) = stname(i)(1:len) ! station name
c
	 call s_len(atype_in(i), len)
	 if(len .ne. 0) then
	    atype(np)(1:len) = atype_in(i)(1:len) ! auto stn type
	 endif
c
	 call s_len(reptype_in(i), len)
	 if(len .ne. 0) then
	    reptype(np)(1:len) = reptype_in(i)(1:len) ! report type
	 endif
c
 	 weather(np)(1:25) = wx(i)(1:25)        ! present weather
         call filter_string(weather(np))

	 provider(np)(1:11) = c11_provider      ! data provider 
	 wmoid(np) = wmoid_in(i)
c
	 store_1(np,1) = lats(i)                ! station latitude
	 store_1(np,2) = lons(i)                ! station longitude
	 store_1(np,3) = elev(i)                ! station elevation
	 store_1(np,4) = rtime                  ! observation time
c
	 store_2(np,1) = temp_f                 ! temperature (deg f)
	 store_2(np,2) = dewp_f                 ! dew point (deg f)
	 store_2(np,3) = rh(i)                  ! relative humidity
c
	 store_3(np,1) = dd(i)                  ! wind dir (deg)
	 store_3(np,2) = ff(i)                  ! wind speed (kt)
	 store_3(np,3) = ddg(i)                 ! wind gust dir (deg)
	 store_3(np,4) = ffg(i)                 ! wind gust speed (kt)
c
	 store_4(np,1) = alt(i)                 ! altimeter setting (mb)
	 store_4(np,2) = stnp(i)                ! station pressure (mb)
	 store_4(np,3) = mslp(i)                ! msl pressure (mb)
	 store_4(np,4) = float(dpchar(i))       ! 3-h press change character
         store_4(np,5) = dp(i)                  ! 3-h press change (mb)
c
	 store_5(np,1) = vis(i)                 ! visibility (miles)
	 store_5(np,2) = sr(i)                  ! solar radiation 
	 store_5(np,3) = st(i)                  ! soil/water temperature
	 store_5(np,4) = badflag                ! soil moisture
c
	 store_6(np,1) = pcp1(i)                ! 1-h precipitation
	 store_6(np,2) = pcp3(i)                ! 3-h precipitation
	 store_6(np,3) = pcp6(i)                ! 6-h precipitation
	 store_6(np,4) = pcp24(i)               ! 24-h precipitation
	 store_6(np,5) = snowcvr(i)             ! snow cover
c
	 store_7(np,1) = float(k_layers)        ! number of cloud layers
	 store_7(np,2) = max24t(i)              ! 24-h max temperature
	 store_7(np,3) = min24t(i)              ! 24-h min temperature
c
c.....	store cloud info if we have any. 
c
	 if(k_layers .gt. 0) then
	   do ii=1,k_layers
	     store_cldht(np,ii) = ht(ii,i)
	     store_cldamt(np,ii)(1:1) = ' '
	     store_cldamt(np,ii)(2:4) = cvr(ii,i)(1:3)
	   enddo !ii
	 endif
c
c
  125	 continue

         write(6,*)' n_sao_all/np = ',n_sao_all,np
c
c
c.....  that's it...lets go home.
c
	 print *,' found ',n_sao_b,' metar/synops in the laps box'
	 print *,' found ',n_sao_g,' metar/synops in the laps grid'
	 print *,' '
	 jstatus = 1		! everything's ok...
	 return
c
 990	 continue		! no data available
	 jstatus = 0
	 print *,' warning: no data available from read_metar.'
	 return
c
	 end
c
c
	function c2f(temp_c)
c
c       takes a single value in deg c and returns deg f.
c
	c2f = (temp_c * 9./5.) + 32.
c
	return
	end


        subroutine sfc_climo_qc_r(c_var,arg)

        character(*) c_var

        call get_sfc_badflag(badflag,istatus)
        
        if(c_var .eq. 'alt_mb' .or. c_var .eq. 'mslp_mb')then
            if(arg .gt. 1100.)arg = badflag
            if(arg .lt.  850.)arg = badflag

        elseif(c_var .eq. 'stnp_mb')then
            if(arg .gt. 1100.)arg = badflag
            if(arg .lt.  400.)arg = badflag

        elseif(c_var .eq. 'alt_pa' .or. c_var .eq. 'mslp_pa')then
            if(arg .gt. 110000.)arg = badflag
            if(arg .lt.  85000.)arg = badflag

        elseif(c_var .eq. 'stnp_pa')then
            if(arg .gt. 110000.)arg = badflag
            if(arg .lt.  40000.)arg = badflag

        elseif(c_var .eq. 't_f')then
            if(arg .gt. +145.)arg = badflag
            if(arg .lt. -145.)arg = badflag

        elseif(c_var .eq. 'tgd_k')then
            if(arg .gt. 340.)arg = badflag
            if(arg .lt. 200.)arg = badflag

        elseif(c_var .eq. 'sst_k')then
            if(arg .gt. 313.)arg = badflag
            if(arg .lt. 200.)arg = badflag

        elseif(c_var .eq. 'td_k')then
            if(arg .gt. 320.)arg = badflag
            if(arg .lt. 210.)arg = badflag

        elseif(c_var .eq. 'dir_deg')then
            if(arg .gt. 360.)arg = badflag
            if(arg .lt.   0.)arg = badflag

        elseif(c_var .eq. 'spd_ms')then
            if(arg .gt. 500.)arg = badflag
            if(arg .lt.   0.)arg = badflag

        else
            write(6,*)' warning: unknown variable in sfc_climo_qc_r'
     1	             ,c_var

        endif

        return
        end
