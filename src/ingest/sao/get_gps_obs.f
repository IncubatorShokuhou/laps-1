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
	subroutine get_gps_obs(maxobs,maxsta,i4time_sys,
     &                      path_to_gps_data,gps_format,
     &                      itime_before,itime_after,
     &                      eastg,westg,anorthg,southg,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_obs_g,n_obs_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, jstatus)

c
c*****************************************************************************
c
cdoc	routine to gather data from the gps files for laps.   
c
c*****************************************************************************
c
	include 'netcdf.inc'
c
c.....  input variables/arrays.
c
        integer maxobs ! raw data file
        integer maxsta ! processed stations for lso file
        character  path_to_gps_data*(*)
c
c.....  local variables/arrays
c
        real    lat(ni,nj), lon(ni,nj)
	real*8  timeobs(maxobs)
	real  lats(maxobs), lons(maxobs), elev(maxobs)
	real  t(maxobs), rh(maxobs), stnp(maxobs)

	integer  i4time_ob, before, after, wmoid(maxobs)
	integer    rtime, dpchar(maxobs), iplat_type(maxobs)

        character*9 a9time_file
	character  stname(maxobs)*5, save_stn(maxobs)*8
	character  data_file*255, timech*9, time*4, gps_format*(*)
	character  weather(maxobs)*25, wx(maxsta)*25
	character  reptype(maxobs)*6, atype(maxobs)*6
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
c
	integer    recnum, nf_fid, nf_vid, nf_status
	character  stations(maxsta)*20, provider(maxsta)*11
	character  store_cldamt(maxsta,5)*4
c
c
c.....  start.
c
	ibadflag = int( badflag )
c
c.....	set jstatus flag for the gps data to bad until we find otherwise.
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
        ibox_points = box_length / (grid_spacing / 1000.)  !in km
c
c.....	zero out the counters.
c
	n_obs_g = 0		! # of gps obs in the laps grid
	n_obs_b = 0		! # of gps obs in the box

        call s_len(gps_format, len_gps_format)
        if(gps_format(1:len_gps_format) .eq. 'nimbus' .or. ! fsl netcdf format
     1     gps_format(1:len_gps_format) .eq. 'cwb'   )then ! cwb netcdf format
            call s_len(path_to_gps_data,len_path)

            call get_file_time(path_to_gps_data(1:len_path)
     1                        ,i4time_sys,i4time_nearest)
            if(abs(i4time_sys - i4time_nearest) .le. 900)then
                i4time_file = i4time_nearest
            else
                write(6,*)'no data available within 900s of systime'
                goto990
            endif
c
!           i4time_file = i4time_sys - 900
            call make_fnam_lp(i4time_file,a9time_file,istatus)

            if(gps_format(1:len_gps_format) .eq. 'nimbus') then
               data_file = path_to_gps_data(1:len_path)//a9time_file
     1                                                 //'0030o.nc'
            elseif(gps_format(1:len_gps_format) .eq. 'cwb') then
               data_file = path_to_gps_data(1:len_path)//a9time_file
     1                                                 //'00.cwb.nc'
            endif
c
c.....      get the data from the netcdf file.  first, open the file.
c.....      if not there, return.
c
	    nf_status = nf_open(data_file,nf_nowrite,nf_fid)

	    if(nf_status.ne.nf_noerr) then
	       print *, nf_strerror(nf_status)
	       print *, data_file
               go to 990
            else
               write(6,*)' file opened successfully'
!              write(6,*)' returning since gps code is not tested yet'
!              go to 990
	    endif
c
c.....      get the dimension of some of the variables.
c
c.....      "recnum"
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

	    call read_gps(nf_fid , recnum, 
     &         stnp, elev, lats, lons, 
     &         t, rh, timeobs, stname)

            i4time_offset=315619200
c
        else
            write(6,*)' unknown gps format ',gps_format
            go to 990

        endif

	if(istatus .ne. 1) go to 990
	n_gps_raw = recnum
c
c.....  first check the data coming from the netcdf file.  there can be
c.....  "floatinf" (used as fill value) in some of the variables.  these
c.....  are not handled the same by different operating systems.  for 
c.....  example, ibm systems make "floatinf" into "nan" and store them that
c.....  way in the file, which messes up other laps routines.  this code
c.....  checks for "floatinf" and sets the variable to 'badflag'.  if the
c.....  "floatinf" is in the lat, lon, elevation, or time of observation,
c.....  we toss the whole ob since we can't be sure where it is.
c
	do i=1,n_gps_raw
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
	   if( nanf( stnp(i) ) .eq. 1 ) stnp(i)  = badflag
	   if( nanf( t(i)    ) .eq. 1 ) t(i)     = badflag
	   if( nanf( rh(i)   ) .eq. 1 ) rh(i)    = badflag
c
	enddo !i
c
c.....  set up the time window.
c
	before = i4time_sys - itime_before
	after  = i4time_sys + itime_after
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
	do 125 i=1,n_gps_raw
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
          if(elev(i).gt.5200. .or. elev(i).lt.-400.) go to 125
c
c.....  check to see if its in the desired time window.
c
	  i4time_ob = nint(timeobs(i)) + i4time_offset
	  if(i4time_ob.lt.before .or. i4time_ob.gt.after) go to 125
c
c.....  right time, right location...

 	  call make_fnam_lp(i4time_ob,timech,istatus)
	  time = timech(6:9)
	  read(time,*) rtime
c
c.....  check if station is reported more than once this
c.....  time period.
c
	  if(jfirst .eq. 1) then
	     icount = 1
	     save_stn(1) = stname(i)
	     jfirst = 0
	     go to 150
	  endif
c
	  do k=1,icount
	     if(stname(i) .eq. save_stn(k)) go to 125
	  enddo !k
c
	  icount = icount + 1
	  save_stn(icount) = stname(i)  ! only one...save for checking
c
 150	  nn = nn + 1
	  n_obs_b = n_obs_b + 1     !station is in the box
c
c.....  check if its in the laps grid.
c
          if(ri_loc.lt.1. .or. ri_loc.gt.float(ni)) go to 151  !off grid
          if(rj_loc.lt.1. .or. rj_loc.gt.float(nj)) go to 151  !off grid
	  n_obs_g = n_obs_g + 1  !on grid...count it
 151	  continue
c
c.....	figure out the cloud data.
c.....  note: not reading cloud data from ship/gps file.  the data
c.....        is too ambiguous for laps use at this time.
c
	  kkk = 0               ! number of cloud layers
c
c
c.....  convert units for storage.
c
c.....  temperature and dewpoint
c
	temp_c = t(i)                         
	if(temp_c .lt. -85. .or. temp_c .gt. +70.)temp_c = badflag
	if(temp_c .eq. badflag) then                 ! t bad?
	   temp_f = badflag                          !          bag
	else
	   temp_f = c_to_f(temp_c)                   ! c to f
	endif
c
	dewp_f = badflag
c
c..... pressure...station pressure
c
	if(stnp(i).lt.400. .or. stnp(i).gt.1200.) then
	   stnp(i) = badflag
	endif
c
c..... fill the expected accuracy arrays.  values are based on information
c..... in the 'coastal-marine automated network (c-man) users guide', and 
c..... we assume that they are about the same for gpss, ships, and c-man
c..... sites.  note that we convert the units to match what we're using here.
c
c..... temperature (deg f)
c
	fon = 9. / 5.  !ratio when converting c to f
	store_2ea(nn,1) = 5.0 * fon        ! start...we don't know what we have
	if(temp_f .ne. badflag) then
	   if(temp_f.ge.c2f(-40.) .and. temp_f.le.c2f(50.)) then
	      store_2ea(nn,1) = 1.0 * fon  ! conv to deg f
	   endif
	endif
c
c..... dew point (deg f).  also estimate a rh accuracy based on the dew point.
c..... estimates for the rh expected accuracy are from playing around with the
c..... psychrometric tables for various t/td combinations.
c
         store_2ea(nn,2) = 0.0             ! dew point not reported
         store_2ea(nn,3) = 10.0            ! relative humidity %
c
c..... wind (not reported)
c
         store_3ea(nn,1) = 0.00            ! deg
         store_3ea(nn,2) = 0.00            ! kt
c
c..... pressure and altimeter (mb)
c
         store_4ea(nn,1) = 1.00            ! pressure (mb)
         store_4ea(nn,2) = 0.00            ! altimeter (mb)
c
c..... other stuff (don't report these). 
c 
         store_5ea(nn,1) = 0.0             ! visibility 
         store_5ea(nn,2) = 0.0             ! solar radiation       
         store_5ea(nn,3) = 0.0             ! soil/water temperature
         store_5ea(nn,4) = 0.0             ! soil moisture
c
         store_6ea(nn,1) = 0.0             ! precipitation (in)
         store_6ea(nn,2) = 0.0             ! snow cover (in) 
c
c
c..... output the data to the storage arrays
c
	 call s_len(stname(i), len)
	 stations(nn)(1:len) = stname(i)(1:len) ! station name
	 provider(nn)(1:11) = 'nws        '     ! data provider (all from nws)
	 reptype(nn)(1:6) = 'gps   '            ! station type
	 atype(nn)(1:6) ='unk   '               ! used here for moving/fixed stns
	 wmoid(nn) = ibadflag                   ! wmo id...not applicable here
	 weather(nn)(1:25) = 
     1              'unk                      ' ! present weather
c       
	 store_1(nn,1) = lats(i)                ! station latitude
	 store_1(nn,2) = lons(i)                ! station longitude
	 store_1(nn,3) = elev(i)                ! station elevation
	 store_1(nn,4) = rtime                  ! observation time
c
	 store_2(nn,1) = temp_f                 ! temperature (deg f)
	 store_2(nn,2) = badflag                ! dew point (deg f)
	 store_2(nn,3) = rh(i)                  ! relative humidity
c
	 store_3(nn,1) = badflag                ! wind dir (deg)
	 store_3(nn,2) = badflag                ! wind speed (kt)
	 store_3(nn,3) = badflag                ! wind gust dir (deg)
	 store_3(nn,4) = badflag                ! wind gust speed (kt)
c
	 store_4(nn,1) = badflag                ! altimeter setting (mb)
	 store_4(nn,2) = stnp(i)                ! station pressure (mb)
	 store_4(nn,3) = badflag                ! msl pressure (mb)
         store_4(nn,4) = badflag
         store_4(nn,5) = badflag                ! 3-h press change (mb)
c
	 store_5(nn,1) = badflag                ! visibility (miles)
	 store_5(nn,2) = badflag                ! solar radiation 
	 store_5(nn,3) = badflag                ! soil/water temperature (f)
	 store_5(nn,4) = badflag                ! soil moisture 
c
	 store_6(nn,1) = badflag                ! 1-h precipitation (in)
	 store_6(nn,2) = badflag                ! 3-h precipitation
	 store_6(nn,3) = badflag                ! 6-h precipitation
	 store_6(nn,4) = badflag                ! 24-h precipitation
	 store_6(nn,5) = badflag                ! snow cover
c
	 store_7(nn,1) = 0                      ! number of cloud layers
	 store_7(nn,2) = badflag                ! 24-h max temperature
	 store_7(nn,3) = badflag                ! 24-h min temperature
c
  125	 continue
c
c
c.....  that's it...lets go home.
c
	 print *,' found ',n_obs_b,' gps obs in the laps box'
	 print *,' found ',n_obs_g,' gps obs in the laps grid'
	 print *,' '
	 jstatus = 1		! everything's ok...
	 return
c
 990	 continue		! no data available
	 jstatus = 0
	 print *,' warning.  no data available from read_gps'
	 return
c
	 end
