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
	subroutine read_master_file_af(master_file,maxsta,n_master,
     &               stn_name,stn_lat,stn_lon,stn_elev,stn_ii,
     &               stn_jj,n_updates,stn_t,trend_t,stn_td,trend_td,
     &               stn_u,trend_u,stn_v,trend_v,stn_alt,trend_alt,
     &               jstatus)
c
c*****************************************************************************
c
c	routine to read the laps "master" data file.   the file contains
c       the station id, obs, and trend, for the past several hours (usually
c       three).
c
c	changes:
c		p. stamus  02-24-97  original version.
c
c*****************************************************************************
c
	real stn_lat(maxsta), stn_lon(maxsta), stn_elev(maxsta)
        real stn_t(maxsta), trend_t(maxsta)
        real stn_td(maxsta), trend_td(maxsta)
        real stn_u(maxsta), trend_u(maxsta)
        real stn_v(maxsta), trend_v(maxsta)
        real stn_alt(maxsta), trend_alt(maxsta)
c
        integer stn_ii(maxsta), stn_jj(maxsta)
        integer n_updates(maxsta)
	integer jstatus
c
        character master_file*70
	character stn_name(maxsta)*5
c
c
c.....	open the file.
c
	jstatus = 0
	open(11,file=master_file,status='unknown',err=500)
c
c.....	read the header.
c
	read(11,900,err=500,end=500) n_master      ! # of stations in this file.
 900    format(i7)
c
c.....	read the data.
c
	do k=1,n_master
	  read(11,910) stn_name(k)(1:5), stn_lat(k), stn_lon(k),
     &             stn_elev(k), stn_ii(k), stn_jj(k), n_updates(k)
c
	  read(11,920) stn_t(k),   trend_t(k)
	  read(11,920) stn_td(k),  trend_td(k)
	  read(11,920) stn_u(k),   trend_u(k)
	  read(11,920) stn_v(k),   trend_v(k)
	  read(11,920) stn_alt(k), trend_alt(k)
c
	enddo !k
c
 910      format(1x,a5,1x,2(f8.2,1x),f6.0,1x,i4,1x,i4,1x,i3)
 920      format(4x,f10.1,1x,f10.2)
c
c.....  end of data reading.  let's go home...
c
	close(11)	
	jstatus = 1
	return
c
c.....  error opening the file.
c
 500	continue
	close(11)
	print *,' error. unable to open master file.'
	jstatus = -1
	return
c
	end
c
c
	subroutine write_master_file_af(master_file,maxsta,n_master,
     &               stn_name,stn_lat,stn_lon,stn_elev,stn_ii,
     &               stn_jj,n_updates,stn_t,trend_t,stn_td,trend_td,
     &               stn_u,trend_u,stn_v,trend_v,stn_alt,trend_alt)
c
c*****************************************************************************
c
c	routine to write the laps "master" data file.   the file contains
c       the station id, obs, and trend, for the past several hours (usually
c       three).
c
c	changes:
c		p. stamus  02-24-97  original version.
c
c*****************************************************************************
c
	real stn_lat(maxsta), stn_lon(maxsta), stn_elev(maxsta)
        real stn_t(maxsta), trend_t(maxsta)
        real stn_td(maxsta), trend_td(maxsta)
        real stn_u(maxsta), trend_u(maxsta)
        real stn_v(maxsta), trend_v(maxsta)
        real stn_alt(maxsta), trend_alt(maxsta)
c
        integer stn_ii(maxsta), stn_jj(maxsta)
        integer n_updates(maxsta)
c
        character master_file*70
	character stn_name(maxsta)*5
c
c
c.....	open the file.
c
	open(12,file=master_file,status='unknown')
c
c.....	write the header.
c
	write(12,900) n_master           ! # of stations in this file.
 900    format(i7)
c
c.....	write the data.
c
	do k=1,n_master
	  write(12,910) stn_name(k)(1:5), stn_lat(k), stn_lon(k),
     &             stn_elev(k), stn_ii(k), stn_jj(k), n_updates(k)
c
	  write(12,920) stn_t(k),   trend_t(k)
	  write(12,920) stn_td(k),  trend_td(k)
	  write(12,920) stn_u(k),   trend_u(k)
	  write(12,920) stn_v(k),   trend_v(k)
	  write(12,920) stn_alt(k), trend_alt(k)
c
	enddo !k
c
 910      format(1x,a5,1x,2(f8.2,1x),f6.0,1x,i4,1x,i4,1x,i3)
 920      format(4x,f10.1,1x,f10.2)
c
	endfile(12)
	close(12)	
c
c..... end of data writing.  let's go home...
c
	return
	end
c
c
	subroutine windconv_af(uwind,vwind,direction,speed)
c
c-----  given wind components, calculate the corresponding speed and direction.

c
c argument	i/o	type			description
c --------	---	----	-----------------------------------------------
c uwind		 i	r*4	u-component of wind
c vwind		 i	r*4	v-component of wind
c direction	 o	r*4	wind direction (meteorological degrees)
c speed		 o	r*4	wind speed (same units as input arguments)
c
c
	parameter(flag = -99.9)
c
	real		uwind,vwind,direction,speed
c
	if(uwind.eq.flag .or. vwind.eq.flag) then
	   speed = flag
	   direction = flag
	elseif(abs(uwind).gt.200. .or. abs(vwind).gt.200.) then
	   speed = flag
	   direction = flag
	elseif(uwind.eq.0.0 .and. vwind.eq.0.0) then
	   speed = 0.0
	   direction = 0.0	!undefined
	else
	   speed = sqrt(uwind*uwind + vwind*vwind) !speed
	   direction = 57.2957795 * (atan2(uwind,vwind)) + 180.	!dir
	endif
c
	return
	end
c
c
	subroutine decomp_wind_af(dd,ff,ucomp,vcomp,status)
c***decompose vector wind into u and v

c	j. wakefield	16 sep 83	original version
c	p. stamus       29 apr 93	unix version

c argument	i/o	type			description
c --------	---	----	-----------------------------------------------
c dd		 i	r*4	wind direction (meteorological degrees)
c ff		 i	r*4	wind speed
c ucomp		 o	r*4	u-component of wind
c vcomp		 o	r*4	v-component of wind
c status	 o	i*4	standard system status

	parameter	(flag = -99.9)

	real		dd,ff,ucomp,vcomp
	integer	istatus

	istatus = 1

	if(dd.eq.flag .or. ff.eq.flag) then
	   ucomp = flag
	   vcomp = flag
	elseif(ff .eq. 0.) then
	   ucomp = 0.0
	   vcomp = 0.0
	elseif(dd.ge.0. .and. dd.le.360.) then
	   angle = .01745239 * dd            !conv to radians
	   ucomp = -ff * sin(angle)
	   vcomp = -ff * cos(angle)
	else
	   ucomp = flag
	   vcomp = flag
	   istatus = 0
	endif
c
	return
	end
c
c
	subroutine find_ij_af(lat_s,lon_s,lat,lon,numsta,mxsta,
     &                     ni,nj,ii,jj,rii,rjj)
c
c.....	routine to find the i,j locations for each station.  do not "round"
c.....  the ii,jj's "up"...straight truncation puts the ob at the proper
c.....  grid point on the major grid.
c
	real lat_s(mxsta), lon_s(mxsta)
        real lat(ni,nj), lon(ni,nj)
	integer ii(mxsta), jj(mxsta)
        real rii(mxsta), rjj(mxsta)
c
	do ista=1,numsta
          call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista),lat,lon,
     &       ni,nj,rii(ista),rjj(ista),istatus)
	  ii(ista) = rii(ista)
	  jj(ista) = rjj(ista)
	enddo !ista
c
	return
	end
