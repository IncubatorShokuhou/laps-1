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
      subroutine read_metar(nf_fid , maxskycover, recnum, altimeter,
     &     autostationtype, dewpoint, dpfromtenths, elevation,
     &     latitude, longitude, maxtemp24hour, mintemp24hour,
     &     precip1hour, precip24hour, precip3hour, precip6hour,
     &     presweather, presschange3hour, presschangechar,
     &     reporttype, sealevelpress, skycover, skylayerbase,
     &     snowcover, stationname, tempfromtenths, temperature,
     &     timeobs, visibility, winddir, windgust, windspeed, wmoid,
     &     badflag, istatus)
c
c======================================================================
c
c     routine to read the metar netcdf files at fsl.
c     code created with 'xgennet.pl' by j. edwards, noaa/fsl.
c     
c     original:  p. stamus, noaa/fsl  12 mar 1998
c     changes:
c        p. stamus, noaa/fsl  01 sep 1998  fix for differences in
c           /public vs wfo cdls (no snowcover in wfo at this time.)
c
c======================================================================
c
      include 'netcdf.inc'
c
      integer maxskycover, recnum, nf_fid, nf_vid, nf_status
      integer ifilval
c
      character*6 autostationtype(recnum)
      character*25 presweather(recnum)
      character*6 reporttype(recnum)
      character*8 skycover( maxskycover, recnum)
      character*5 stationname(recnum)
      integer   presschangechar(recnum), wmoid(recnum)

      double precision timeobs(recnum), dfilval

      real altimeter(recnum), dewpoint(recnum), dpfromtenths(recnum)
      real elevation(recnum), latitude(recnum), longitude(recnum)
      real maxtemp24hour(recnum), mintemp24hour(recnum)
      real precip1hour(recnum), precip24hour(recnum)
      real precip3hour(recnum), precip6hour(recnum)
      real presschange3hour(recnum), sealevelpress(recnum)
      real skylayerbase( maxskycover, recnum), snowcover(recnum)
      real tempfromtenths(recnum), temperature(recnum)
      real visibility(recnum), winddir(recnum), windgust(recnum)
      real windspeed(recnum)

      real filval
c
c
c..... start.
c
      istatus = 0
c
c     variable        netcdf long name
c      autostationtype"automated station type" 
c
        nf_status = nf_inq_varid(nf_fid,'autostationtype',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var autostationtype'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,autostationtype)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ autostationtype '
      endif
c
c     variable        netcdf long name
c      presweather  "present weather" 
c
        nf_status = nf_inq_varid(nf_fid,'presweather',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var presweather'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,presweather)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ presweather '
      endif
c
c     variable        netcdf long name
c      reporttype   "report type" 
c
        nf_status = nf_inq_varid(nf_fid,'reporttype',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reporttype'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,reporttype)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ reporttype '
      endif
c
c     variable        netcdf long name
c      skycover     "sky cover" 
c
        nf_status = nf_inq_varid(nf_fid,'skycover',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var skycover'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,skycover)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ skycover '
      endif
c
c     variable        netcdf long name
c      stationname  "alphanumeric station identification" 
c
        nf_status = nf_inq_varid(nf_fid,'stationname',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationname'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,stationname)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stationname '
      endif
c
c     variable        netcdf long name
c      presschangechar"character of pressure change" 
c
        nf_status = nf_inq_varid(nf_fid,'presschangechar',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var presschangechar'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,presschangechar)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ presschangechar '
      endif
        nf_status = nf_get_att_int(nf_fid,nf_vid,'_fillvalue',ifilval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var presschangechar'
      endif
      do i=1,recnum
         if(presschangechar(i) .eq. ifilval) 
     &                         presschangechar(i) = int(badflag)
      enddo !i
c
c     variable        netcdf long name
c      wmoid        "numeric wmo identification" 
c
        nf_status = nf_inq_varid(nf_fid,'wmoid',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wmoid'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,wmoid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ wmoid '
      endif
        nf_status = nf_get_att_int(nf_fid,nf_vid,'_fillvalue',ifilval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var wmoid'
      endif
      do i=1,recnum
         if(wmoid(i) .eq. ifilval) wmoid(i) = int(badflag)
      enddo !i
c
c     variable        netcdf long name
c      altimeter    "altimeter setting" 
c
        nf_status = nf_inq_varid(nf_fid,'altimeter',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var altimeter'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,altimeter)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ altimeter '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var altimeter'
      endif
      call ck_array_real(altimeter, recnum, filval, badflag)
c
c     variable        netcdf long name
c      dewpoint     "dewpoint" 
c
        nf_status = nf_inq_varid(nf_fid,'dewpoint',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dewpoint'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,dewpoint)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ dewpoint '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var dewpoint'
      endif
      call ck_array_real(dewpoint, recnum, filval, badflag)
c
c     variable        netcdf long name
c      dpfromtenths "dewpoint from tenths of a degree celsius" 
c
        nf_status = nf_inq_varid(nf_fid,'dpfromtenths',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dpfromtenths'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,dpfromtenths)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ dpfromtenths '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
        print *,'in var dpfromtenths'
      endif
      call ck_array_real(dpfromtenths, recnum, filval, badflag)
c
c     variable        netcdf long name
c      elevation    "elevation" 
c
        nf_status = nf_inq_varid(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevation'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,elevation)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ elevation '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in var elevation'
      endif
      call ck_array_real(elevation, recnum, filval, badflag)
c
c     variable        netcdf long name
c      latitude     "latitude" 
c
        nf_status = nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,latitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ latitude '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var latitude' 
      endif
      call ck_array_real(latitude, recnum, filval, badflag)
c
c     variable        netcdf long name
c      longitude    "longitude" 
c
        nf_status = nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,longitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ longitude '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var longitude' 
      endif
      call ck_array_real(longitude, recnum, filval, badflag)
c
c     variable        netcdf long name
c      maxtemp24hour"24 hour max temperature" 
c
        nf_status = nf_inq_varid(nf_fid,'maxtemp24hour',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var maxtemp24hour'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,maxtemp24hour)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ maxtemp24hour '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in var maxtemp24hour'
      endif
      call ck_array_real(maxtemp24hour, recnum, filval, badflag)
c
c     variable        netcdf long name
c      mintemp24hour"24 hour min temperature" 
c
        nf_status = nf_inq_varid(nf_fid,'mintemp24hour',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var mintemp24hour'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,mintemp24hour)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ mintemp24hour '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in var mintemp24hour'
      endif
      call ck_array_real(mintemp24hour, recnum, filval, badflag)
c
c     variable        netcdf long name
c      precip1hour  "1 hour precipitation" 
c
        nf_status = nf_inq_varid(nf_fid,'precip1hour',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var precip1hour'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,precip1hour)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ precip1hour '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var precip1hour'
      endif
      call ck_array_real(precip1hour, recnum, filval, badflag)
c
c     variable        netcdf long name
c      precip24hour "24 hour precipitation" 
c
        nf_status = nf_inq_varid(nf_fid,'precip24hour',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var precip24hour'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,precip24hour)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ precip24hour '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var precip24hour'
      endif
      call ck_array_real(precip24hour, recnum, filval, badflag)
c
c     variable        netcdf long name
c      precip3hour  "3 hour precipitation" 
c
        nf_status = nf_inq_varid(nf_fid,'precip3hour',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var precip3hour'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,precip3hour)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ precip3hour '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var precip3hour'
      endif
      call ck_array_real(precip3hour, recnum, filval, badflag)
c
c     variable        netcdf long name
c      precip6hour  "6 hour precipitation" 
c
        nf_status = nf_inq_varid(nf_fid,'precip6hour',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var precip6hour'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,precip6hour)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ precip6hour '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var precip6hour'
      endif
      call ck_array_real(precip6hour, recnum, filval, badflag)
c
c     variable        netcdf long name
c      presschange3hour"3 hour pressure change" 
c
        nf_status = nf_inq_varid(nf_fid,'presschange3hour',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var presschange3hour'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,presschange3hour)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ presschange3hour '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var presschange3hour'
      endif
      call ck_array_real(presschange3hour, recnum, filval, badflag)
c
c     variable        netcdf long name
c      sealevelpress"sea level pressure" 
c
        nf_status = nf_inq_varid(nf_fid,'sealevelpress',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sealevelpress'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,sealevelpress)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ sealevelpress '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var sealevelpress'
      endif
      call ck_array_real(sealevelpress, recnum, filval, badflag)
c
c     variable        netcdf long name
c      skylayerbase "sky cover layer base" 
c
        nf_status = nf_inq_varid(nf_fid,'skylayerbase',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var skylayerbase'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,skylayerbase)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ skylayerbase '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var skylayerbase'
      endif
      call ck_array_real(skylayerbase, recnum, filval, badflag)
c
c     variable        netcdf long name
c      snowcover    "snow cover" 
c
        nf_status = nf_inq_varid(nf_fid,'snowcover',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var snowcover'
        do i=1,recnum
           snowcover(i) = badflag
        enddo !i
        go to 500
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,snowcover)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ snowcover '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var snowcover'
      endif
      call ck_array_real(snowcover, recnum, filval, badflag)
c
c     variable        netcdf long name
c      tempfromtenths"temperature from tenths of a degree celsius" 
c
 500  nf_status = nf_inq_varid(nf_fid,'tempfromtenths',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tempfromtenths'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,tempfromtenths)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ tempfromtenths '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var tempfromtenths'
      endif
      call ck_array_real(tempfromtenths, recnum, filval, badflag)
c
c     variable        netcdf long name
c      temperature  "temperature" 
c
        nf_status = nf_inq_varid(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,temperature)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ temperature '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var temperature'
      endif
      call ck_array_real(temperature, recnum, filval, badflag)
c
c     variable        netcdf long name
c      visibility   "visibility" 
c
        nf_status = nf_inq_varid(nf_fid,'visibility',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var visibility'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,visibility)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ visibility '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var visibility'
      endif
      call ck_array_real(visibility, recnum, filval, badflag)
c
c     variable        netcdf long name
c      winddir      "wind direction" 
c
        nf_status = nf_inq_varid(nf_fid,'winddir',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,winddir)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ winddir '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var winddir'
      endif
      call ck_array_real(winddir, recnum, filval, badflag)
c
c     variable        netcdf long name
c      windgust     "wind gust" 
c
        nf_status = nf_inq_varid(nf_fid,'windgust',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windgust'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,windgust)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ windgust '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var windgust'
      endif
      call ck_array_real(windgust, recnum, filval, badflag)
c
c     variable        netcdf long name
c      windspeed    "wind speed" 
c
        nf_status = nf_inq_varid(nf_fid,'windspeed',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeed'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,windspeed)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ windspeed '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var windspeed'
      endif
      call ck_array_real(windspeed, recnum, filval, badflag)
c
c     variable        netcdf long name
c      timeobs      "time of observation" 
c
        nf_status = nf_inq_varid(nf_fid,'timeobs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timeobs'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,timeobs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ timeobs '
      endif
       nf_status = nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dfilval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var timeobs'
      endif
      do i=1,recnum
         if(timeobs(i) .eq. dfilval) timeobs(i) = badflag
      enddo !i
c
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif
c
c..... that's it.
c
      istatus = 1
c
      return
      end
