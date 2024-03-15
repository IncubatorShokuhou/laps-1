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
	subroutine get_sao_obs_af(filename,path_to_obs,maxsta,
     &                      eastg,westg,anorthg,southg,nn,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      n_sao_g,n_sao_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider,jstatus)
c
c*****************************************************************************
c
c	routine to gather sao data for laps.  air force version.
c
c	changes:
c		p. stamus  10-28-94  original version (from get_surface_obs).
c                          12-08-94  adjust cld amts in character string.
c                          03-14-95  filter out multiples of same obs.
c                          11-10-95  fixed bug in first station id (dups).
c                          06-05-96  code to handle floatinf in data for ibms.
c                          06-24-96  air force version.
c                          03-06-97  pass c*5 instead of c*3 stn names.
c                          10-28-97  changes for dynamic laps.
c                          11-18-98  changes for new lso format.
c                          03-26-99  better drib/ship check.
c	                   04-07-99  change stname to left justify.
c                          06-11-99  change ob location check to gridpt space
c                                      to better handle +/- longitudes.
c                          06-17-99  figure box size in gridpoint space from
c                                      user-defined size (deg) and grid_spacing.
c
c*****************************************************************************
c
	real  timeobs(maxsta)
	real  lats(maxsta), lons(maxsta), elev(maxsta)
	real  t(maxsta), td(maxsta), dd(maxsta), ff(maxsta)
	real  mslp(maxsta), alt(maxsta), ht(5,maxsta), ffg(maxsta)
	real  vis(maxsta)
	real    lat(ni,nj), lon(ni,nj)
c
	real	store_1(maxsta,4),
     &          store_2(maxsta,3), store_2ea(maxsta,3), 
     &          store_3(maxsta,4), store_3ea(maxsta,2), 
     &          store_4(maxsta,5), store_4ea(maxsta,2), 
     &          store_5(maxsta,4), store_5ea(maxsta,4), 
     &          store_6(maxsta,5), store_6ea(maxsta,2), 
     &          store_7(maxsta,3),
     &          store_cldht(maxsta,5)
c
	integer  i4time_ob, wmoid(maxsta)
	integer  rtime, dpchar(maxsta), dp(maxsta)
c
	character  path_to_obs*(*) 
	character  drib*2, ship*2
	character  stname(maxsta)*5, save_stn(maxsta)*5, timech*9, time*4
	character  cvr(5,maxsta)*8 
	character  stations(maxsta)*20, provider(maxsta)*11
	character  weather(maxsta)*25, wx(maxsta)*8
	character  reptype(maxsta)*6, atype(maxsta)*6, filename*9
	character  reptype_in(maxsta)*6, atype_in(maxsta)*6
	character  store_cldamt(maxsta,5)*4 
c
c.....	set jstatus flag for the sao data to bad until we find otherwise.
c
	jstatus = -1

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return
c
c.....  figure out the size of the "box" in gridpoints.  user defines
c.....  the 'box_size' variable in degrees, then we convert that to an
c.....  average number of gridpoints based on the grid spacing.
c
        call get_box_size(box_size,istatus)
        if(istatus .ne. 1)return

	box_length = box_size * 111.137 !km/deg lat (close enough for lon)
	ibox_points = box_length / (grid_spacing / 1000.) !in km
c
c.....	zero out the counters.
c
	n_sao_g = 0		! # of saos in the laps grid
	n_sao_b = 0		! # of saos in the box
c
c.....  call the routine that reads the netcdf sao data files.
c
        call read_sao_af(filename,path_to_obs,n_sao_all,stname,
     &                   lats,lons,elev,
     &                   timeobs,reptype_in,atype_in,
     &                   cvr,ht,vis,weather,
     &                   mslp,t,td,dd,ff,ffg,alt,
     &                   dpchar,dp,maxsta,badflag,istatus)
c
	if(istatus .ne. 1) go to 990
c
c.....  first check the data coming from the netcdf file.  there can be
c.....  "floatinf" (used as fill value) in some of the variables.  these
c.....  are not handled the same by different operating systems.  in our
c.....  case, ibm systems make "floatinf" into "nan" and store them that
c.....  way in the lso, which messes up other laps routines.  this code
c.....  checks for "floatinf" and sets the variable to 'badflag'.  if the
c.....  "floatinf" is in the lat, lon, or elevation, we toss the whole ob
c.....  since we don't know where it is.
c
	do i=1,n_sao_all
c
c.....  toss the ob if lat/lon/elev bad by setting lat to badflag (-99.9),
c.....  which causes the bounds check to think its outside the laps domain.
c
	   if( nanf( lats(i) ) .eq. 1 ) lats(i) = badflag
	   if( nanf( lons(i) ) .eq. 1 ) lats(i) = badflag
	   if( nanf( elev(i) ) .eq. 1 ) lats(i) = badflag
c
	   if( nanf( timeobs(i) ) .eq. 1 ) lats(i) = badflag
c
	   if( nanf( vis(i)   ) .eq. 1 ) vis(i)  = badflag
	   if( nanf( mslp(i)  ) .eq. 1 ) mslp(i) = badflag
	   if( nanf( t(i)     ) .eq. 1 ) t(i)    = badflag
	   if( nanf( td(i)    ) .eq. 1 ) td(i)   = badflag
	   if( nanf( dd(i)    ) .eq. 1 ) dd(i)   = badflag
	   if( nanf( ff(i)    ) .eq. 1 ) ff(i)   = badflag
	   if( nanf( ffg(i)   ) .eq. 1 ) ffg(i)  = badflag
	   if( nanf( alt(i)   ) .eq. 1 ) alt(i)  = badflag
c
	   do j=1,5
	      if( nanf(ht(j,i) ) .eq. 1 ) ht(j,i) = badflag
	   enddo !j
c
	enddo !i
c
c..................................
c.....	now loop over all the obs.
c..................................
c
	idrib_cnt = 0
	iship_cnt = 0
	jfirst = 1
	box_low = 1. - float(ibox_points)    !buffer on west/south side
	box_idir = float( ni + ibox_points)  !buffer on east
	box_jdir = float( nj + ibox_points)  !buffer on north
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
	   if(elev(i).gt.5200. .or. elev(i).lt.-400.) go to 125
c
c.....  check for drifting buoys ('drib') and ship reports ('ship').
c.....  change the names so we can track more than just one of each.
c.....  for ships, change reptype so we can filter it later.
c
	  if(stname(i)(1:4) .eq. 'drib') then
	     idrib_cnt = idrib_cnt + 1
	     write(drib,911) idrib_cnt
	     stname(i)(1:4) = 'by' // drib(1:2)
	  endif
c
	  if(stname(i)(1:4) .eq. 'ship') then
	     iship_cnt = iship_cnt + 1
	     write(ship,911) iship_cnt
	     stname(i)(1:4) = 'sh' // ship(1:2)
	     reptype_in(i)(1:6) = 'ship  '
	  endif
 911	  format(i2.2)
c
c.....  in the box, so lets check if its reported more than once this
c.....  time period.
c
	  if(jfirst .eq. 1) then
	     icount = 1
	     save_stn(1) = stname(i)(1:5)
	     jfirst = 0
	     go to 150
	  endif
c
	  do k=1,icount
	     if(stname(i)(1:5) .eq. save_stn(k)) go to 125
	  enddo !k
c
	  icount = icount + 1
	  save_stn(icount) = stname(i)(1:5)   ! only one...save for checking
c
 150	  nn = nn + 1
	  n_sao_b = n_sao_b + 1		!station is in the box
c
c.....  check if the ob is on the laps grid.
c
cc	 if(lons(i).gt.eastg .or. lons(i).lt.westg) go to 151
cc	 if(lats(i).gt.anorthg .or. lats(i).lt.southg) go to 151
c
	  if(ri_loc.lt.1. .or. ri_loc.gt.float(ni)) go to 151  !off grid
	  if(rj_loc.lt.1. .or. rj_loc.gt.float(nj)) go to 151  !off grid
	  n_sao_g = n_sao_g + 1    !on grid...count it
c
c.....  figure out the report time.
c
 151	 continue
cc 	 i4time_ob = nint(timeobs(i)) + 315619200
cc 	 call make_fnam_lp(i4time_ob,timech,istatus)
cc	 time = timech(6:9)
cc	 read(time,*) rtime
	 rtime = timeobs(i)
c
c.....	figure out the cloud data.
c
	  kkk = 0               ! number of cloud layers
c
	  if(cvr(1,i)(1:1) .eq. ' ') then
	     kkk = 0
	  else
	     do k=1,5
		if(cvr(k,i)(1:1) .ne. ' ') kkk = kkk + 1
	     enddo !k
	  endif
c
	  if(kkk .eq. 0) then	! no cloud data...probably amos station
	    go to 126		! skip rest of cloud stuff
	  endif
	  do ii=1,kkk
	    if(ht(ii,i) .gt. 25000.0) ht(ii,i) = badflag
c
	    if(cvr(ii,i)(1:3) .eq. 'clr') then
	       ht(ii,i) = 3657.4
	    endif
c	    
	    if(cvr(ii,i)(1:3) .eq. 'skc') then
	       ht(ii,i) = 22500.
	    endif
	  enddo !ii
 126	  continue
c
c.....	check cloud info for very high heights...set to max if greater
c.....	also convert agl cloud heights to msl by adding elevation.
c
	if(kkk .gt. 0) then
	  do ii=1,kkk
	     if(ht(ii,i) .ge. 0.) then
		ht(ii,i) = ht(ii,i) + elev(i) ! conv agl to msl
		if(ht(ii,i) .gt. 22500.) ht(ii,i) = 22500.
		if(ht(ii,i) .lt.     0.) ht(ii,i) = badflag
	     endif
	  enddo !ii
	endif
c
c
c.....  convert units for storage, while doing a quick and dirty qc on the data.
c
	if(alt(i).lt.900. .or. alt(i).gt.1150.) alt(i) = badflag !alt (mb)
	if(mslp(i).lt.900. .or. mslp(i).gt.1150.) mslp(i) = badflag !msl p (mb)
c
	temp_k = t(i)
	if(temp_k .eq. badflag) then
	   temp_f = badflag
	else
	   temp_f = ((temp_k - 273.16) * 9./5.) + 32.
	endif
c
	dewp_k = td(i)
	if(dewp_k .eq. badflag) then
	   dewp_f = badflag
	else
	   dewp_f = ((dewp_k - 273.16) * 9./5.) + 32.
	endif
c
	if(dd(i).lt.0. .or. dd(i).gt.360.)dd(i) = badflag
	ff(i) = 1.94254 * ff(i)                             !m/s to kt
	ffg(i) = 1.94254 * ffg(i)                           !m/s to kt
	if(ff(i).lt.0. .or. ff(i).gt.40.)ff(i) = badflag
	if(ffg(i).lt.0. .or. ffg(i).gt.55.)ffg(i) = badflag
	ddg = dd(i)
	if(ffg(i) .eq. badflag) ddg = badflag
c
	vis(i) = 0.621371 * vis(i)                      !km to miles
	if(vis(i) .gt. 200.) vis(i) = badflag			! visibility
	if(vis(i) .lt.   0.) vis(i) = badflag
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
	store_2ea(nn,1) = 5.0 * fon        ! start...we don't know what we have
	if(temp_f .ne. badflag) then
	   if(temp_f.ge.c2f(-62.) .and. temp_f.le.c2f(-50.)) then
	      store_2ea(nn,1) = 1.1 * fon  ! conv to deg f
	   elseif(temp_f.gt.c2f(-50.) .and. temp_f.lt.c2f(50.)) then
	      store_2ea(nn,1) = 0.6 * fon  ! conv to deg f
	   elseif(temp_f.ge.c2f(50.) .and. temp_f.le.c2f(54.)) then
	      store_2ea(nn,1) = 1.1 * fon  ! conv to deg f
	   endif
	endif
c
c..... dew point (deg f).  also estimate a rh accuracy based on the dew point.
c..... estimates for the rh expected accuracy are from playing around with the
c..... psychrometric tables for various t/td combinations (including their
c..... accuracies from the fmh-1 appendix c).
c
	 store_2ea(nn,2) = 5.0 * fon       ! start...don't know what we have 
	 store_2ea(nn,3) = 50.0            ! relative humidity %
	 if(dewp_f .ne. badflag) then
	    if(dewp_f.ge.c2f(-34.) .and. dewp_f.lt.c2f(-24.)) then
	       store_2ea(nn,2) = 2.2 * fon ! conv to deg f
	       store_2ea(nn,3) = 20.0      ! rh (%) 
	    elseif(dewp_f.ge.c2f(-24.) .and. dewp_f.lt.c2f(-1.)) then
	       store_2ea(nn,2) = 1.7 * fon ! conv to deg f
	       store_2ea(nn,3) = 12.0      ! rh (%) 
	    elseif(dewp_f.ge.c2f(-1.) .and. dewp_f.le.c2f(30.)) then
	       store_2ea(nn,2) = 1.1 * fon ! conv to deg f
	       store_2ea(nn,3) = 8.0       ! rh (%) 
	    endif
	 endif
c
c..... wind direction (deg) and speed (kts)
c
	 store_3ea(nn,1) = 10.0    ! deg 
	 store_3ea(nn,2) =  1.0    ! kt
	 if(ff(i) .ne. badflag) then
	    if(ff(i).ge.1.0 .and. ff(i).le.10.0) then
	       store_3ea(nn,2) = 1.0          ! kt
	    elseif(ff(i) .gt. 10.0) then
	       store_3ea(nn,2) = ff(i) * 0.1  ! 10% of speed (kts)
	    endif
c
	    if(ff(i) .ge. 5.0) then    ! dir check
	       store_3ea(nn,1) = 5.0   ! deg
	    endif
	 endif
c
c..... pressure and altimeter (mb)
c
	 store_4ea(nn,1) = 0.68            ! pressure (mb)
	 store_4ea(nn,2) = 0.68            ! altimeter (mb)
c
c..... visibility (miles).  for automated stations use a guess based 
c..... on table c-2 in appendix c of fmh-1.  for manual stations, use
c..... a guess based on the range between reportable values (e.g., for
c..... reported visibility between 0 and 3/8th mile, set accuracy to 
c..... 1/16th mile).  this isn't ideal, but its a start.
c
	 store_5ea(nn,1) = 10.00         ! start with this (miles)
	 if(vis(i) .ne. badflag) then
	    if(atype_in(i)(1:1) .eq. 'a') then   ! have an auto station
	       if(vis(i) .lt. 2.0) then
		  store_5ea(nn,1) = 0.25         ! miles
	       elseif(vis(i).ge.2.0 .and. vis(i).lt.3.0) then
		  store_5ea(nn,1) = 0.50         ! miles
	       elseif(vis(i) .gt. 3.0) then
		  store_5ea(nn,1) = 1.00         ! miles
	       endif
	    else		! have a manual station
	       if(vis(i) .le. 0.375) then
		  store_5ea(nn,1) = 0.0625       ! miles
	       elseif(vis(i).gt.0.375 .and. vis(i).le.2.0) then
		  store_5ea(nn,1) = 0.125        ! miles
	       elseif(vis(i).gt.2.0 .and. vis(i).le.3.0) then
		  store_5ea(nn,1) = 0.25         ! miles
	       elseif(vis(i).gt.3.0 .and. vis(i).le.15.0) then
		  store_5ea(nn,1) = 1.00         ! miles
	       elseif(vis(i) .gt. 15.0) then
		  store_5ea(nn,1) = 5.00         ! miles
	       endif
	    endif
	 endif
c
c..... other stuff.  don't really know about the precip, but probably
c..... worse that this guess.
c
	 store_5ea(nn,2) = 0.0             ! solar radiation 
	 store_5ea(nn,3) = 0.0             ! soil/water temperature
	 store_5ea(nn,4) = 0.0             ! soil moisture 
c
	 store_6ea(nn,1) = 0.01            ! precipitation (in)
	 store_6ea(nn,2) = 1.0             ! snow cover (in) 
c
c
c..... output the data to the storage arrays
c
cc	 call s_len(stname(i), len)
	 len = 5
	 stations(nn)(1:len) = stname(i)(1:len) ! station name
c
	 call s_len(atype_in(i), len)
	 if(len .ne. 0) then
	    atype(nn)(1:len) = atype_in(i)(1:len) ! auto stn type
	 endif
c
	 call s_len(reptype_in(i), len)
	 if(len .ne. 0) then
	    reptype(nn)(1:len) = reptype_in(i)(1:len) ! report type
	 endif
c
	 weather(nn)(1:8) = wx(i)(1:8)          ! present weather
	 provider(nn)(1:11) = 'afwa       '     ! data provider (all from afwa)
	 wmoid(nn) = 0                          ! wmo id number (don't have)
c
	 store_1(nn,1) = lats(i)                ! station latitude
	 store_1(nn,2) = lons(i)                ! station longitude
	 store_1(nn,3) = elev(i)                ! station elevation
	 store_1(nn,4) = rtime                  ! observation time
c
	 store_2(nn,1) = temp_f                 ! temperature (deg f)
	 store_2(nn,2) = dewp_f                 ! dew point (deg f)
	 store_2(nn,3) = badflag                ! relative humidity
c
	 store_3(nn,1) = dd(i)                  ! wind dir (deg)
	 store_3(nn,2) = ff(i)                  ! wind speed (kt)
	 store_3(nn,3) = ddg                    ! wind gust dir (deg)
	 store_3(nn,4) = ffg(i)                 ! wind gust speed (kt)
c
	 store_4(nn,1) = alt(i)                 ! altimeter setting (mb)
	 store_4(nn,2) = badflag                ! station pressure (mb)
	 store_4(nn,3) = mslp(i)                ! msl pressure (mb)
	 store_4(nn,4) = float(dpchar(i))       ! 3-h press change character
         store_4(nn,5) = float(dp(i))           ! 3-h press change (mb)
c
	 store_5(nn,1) = vis(i)                 ! visibility (miles)
	 store_5(nn,2) = badflag                ! solar radiation 
	 store_5(nn,3) = badflag                ! soil/water temperature
	 store_5(nn,4) = badflag                ! soil moisture
c
	 store_6(nn,1) = badflag                ! 1-h precipitation
	 store_6(nn,2) = badflag                ! 3-h precipitation
	 store_6(nn,3) = badflag                ! 6-h precipitation
	 store_6(nn,4) = badflag                ! 24-h precipitation
	 store_6(nn,5) = badflag                ! snow cover
c
	 store_7(nn,1) = float(kkk)             ! number of cloud layers
	 store_7(nn,2) = badflag                ! 24-h max temperature
	 store_7(nn,3) = badflag                ! 24-h min temperature
c
c.....	store cloud info if we have any. 
c
	 if(kkk .gt. 0) then
	   do ii=1,kkk
	     store_cldht(nn,ii) = ht(ii,i)
	     store_cldamt(nn,ii)(1:1) = ' '
	     store_cldamt(nn,ii)(2:4) = cvr(ii,i)(1:3)
	   enddo !ii
	 endif
c
c
  125	 continue
c
c
c.....  that's it...lets go home.
c
	 print *,' found ',idrib_cnt,' drifting buoys in the laps box'
	 print *,' found ',iship_cnt,' "no-name" ships in the laps box'
	 print *,' '
	 print *,' found ',n_sao_b,' saos in the laps box'
	 print *,' found ',n_sao_g,' saos in the laps grid'
	 print *,' '
	 jstatus = 1		! everything's ok...
	 return
c
 990	 continue		! no data available
	 jstatus = 0
	 print *,' error.  no data available from read_sao_af.'
	 return
c
	 end
c
c
