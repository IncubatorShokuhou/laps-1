c
c
      subroutine read_local_obs(nf_fid , recnum, altimeter, 
     +     dataprovider, solarradiation, seasurfacetemp,    
     +     soiltemperature,
     +     dewpoint, elevation, latitude, longitude, observationtime,
     +     presweather, relhumidity, rhchangetime, sealevelpressure,
     +     stationid, stationpresschangetime, stationpressure,
     +     stationtype, tempchangetime, temperature, visibility,
     +     winddir, winddirchangetime, winddirmax, windgust,
     +     windgustchangetime, windspeed, windspeedchangetime,
     &     badflag, istatus)
c
c**********************************************************************
c
c     routine to read the ldad netcdf mesonet observation files at fsl.
c     code created with 'xgennet.pl' by j. edwards, noaa/fsl.
c     
c     original:  p. stamus, noaa/fsl  28 aug 1998
c	         09-30-98  p. stamus
c                     housekeeping changes.
c                10-18-99  p. stamus
c                     add check for missing values.
c
c**********************************************************************
c
      include 'netcdf.inc'
      integer recnum, nf_fid, nf_vid, nf_status !, ifilval

      character*11 dataprovider(recnum)
      character*25 presweather(recnum)
      character*6 stationid(recnum)
      character*11 stationtype(recnum)
      real altimeter(recnum), dewpoint(recnum), elevation(recnum),
     +     latitude(recnum), longitude(recnum), relhumidity(recnum),
     +     sealevelpressure(recnum), stationpressure(recnum),
     +     temperature(recnum), visibility(recnum), winddir(recnum),
     +     winddirmax(recnum), windgust(recnum), windspeed(recnum),
     +     solarradiation(recnum), seasurfacetemp(recnum),
     +     soiltemperature(recnum)
      real filval, misval

      double precision observationtime(recnum), rhchangetime(recnum),
     +     stationpresschangetime(recnum), tempchangetime(recnum),
     +     winddirchangetime(recnum), windgustchangetime(recnum),
     +     windspeedchangetime(recnum), dfilval, dmisval
c
c..... start here.
c
      istatus = 0
c
c     variable        netcdf long name
c      dataprovider "ldad data provider" 
c
        nf_status = nf_inq_varid(nf_fid,'dataprovider',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dataprovider'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,dataprovider)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ dataprovider '
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
c      stationid    "alphanumeric station id" 
c
        nf_status = nf_inq_varid(nf_fid,'stationid',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationid'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,stationid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stationid '
      endif
c
c     variable        netcdf long name
c      stationtype  "ldad station type" 
c
        nf_status = nf_inq_varid(nf_fid,'stationtype',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationtype'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,stationtype)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stationtype '
      endif
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
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var altimeter'
      endif
      call ck_array_real(altimeter, recnum, misval, badflag)
c
c     variable        netcdf long name
c      dewpoint     "dew point temperature" 
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
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var dewpoint'
      endif
      call ck_array_real(dewpoint, recnum, misval, badflag)
c
c     variable        netcdf long name
c      dewpoint     "solar radiation" 
c
        nf_status = nf_inq_varid(nf_fid,'solarradiation',nf_vid)       
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var solarradiation'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,solarradiation)       
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ solarradiation '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var solarradiation'
      endif
      call ck_array_real(solarradiation, recnum, filval, badflag)       
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var solarradiation'
      endif
      call ck_array_real(solarradiation, recnum, misval, badflag)       
c
c     variable        netcdf long name
c      dewpoint     "sea surface temperature" 
c
        nf_status = nf_inq_varid(nf_fid,'seasurfacetemp',nf_vid)       
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var seasurfacetemp'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,seasurfacetemp)       
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ seasurfacetemp '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var seasurfacetemp'
      endif
      call ck_array_real(seasurfacetemp, recnum, filval, badflag)       
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var seasurfacetemp'
      endif
      call ck_array_real(seasurfacetemp, recnum, misval, badflag)       
c
c     variable        netcdf long name
c      "soil temperature" 
c
        nf_status = nf_inq_varid(nf_fid,'soiltemperature',nf_vid)       
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var soiltemperature'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,soiltemperature)       
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ soiltemperature '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var soiltemperature'
      endif
      call ck_array_real(soiltemperature, recnum, filval, badflag)       
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var soiltemperature'
      endif
      call ck_array_real(soiltemperature, recnum, misval, badflag)       
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
         print *, ' in var elevation'
      endif
      call ck_array_real(elevation, recnum, filval, badflag)
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var elevation'
      endif
      call ck_array_real(elevation, recnum, misval, badflag)
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
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var latitude'
      endif
      call ck_array_real(latitude, recnum, misval, badflag)
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
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var longitude'
      endif
      call ck_array_real(longitude, recnum, misval, badflag)
c
c     variable        netcdf long name
c      relhumidity  "relative humidity" 
c
        nf_status = nf_inq_varid(nf_fid,'relhumidity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidity'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,relhumidity)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ relhumidity '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var relhumidity'
      endif
      call ck_array_real(relhumidity, recnum, filval, badflag)
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var relhumidity'
      endif
      call ck_array_real(relhumidity, recnum, misval, badflag)
c
c     variable        netcdf long name
c      sealevelpressure"sea level pressure" 
c
        nf_status = nf_inq_varid(nf_fid,'sealevelpressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sealevelpressure'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,sealevelpressure)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ sealevelpressure '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var sealevelpressure'
      endif
      call ck_array_real(sealevelpressure, recnum, filval, badflag)
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var sealevelpressure'
      endif
      call ck_array_real(sealevelpressure, recnum, misval, badflag)
c
c     variable        netcdf long name
c      stationpressure"station pressure" 
c
        nf_status = nf_inq_varid(nf_fid,'stationpressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationpressure'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,stationpressure)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stationpressure '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var stationpressure'
      endif
      call ck_array_real(stationpressure, recnum, filval, badflag)
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var stationpressure'
      endif
      call ck_array_real(stationpressure, recnum, misval, badflag)
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
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var temperature'
      endif
      call ck_array_real(temperature, recnum, misval, badflag)
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
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var visibility'
      endif
      call ck_array_real(visibility, recnum, misval, badflag)
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
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var winddir'
      endif
      call ck_array_real(winddir, recnum, misval, badflag)
c
c     variable        netcdf long name
c      winddirmax   "wind direction at gust" 
c
        nf_status = nf_inq_varid(nf_fid,'winddirmax',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddirmax'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,winddirmax)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ winddirmax '
      endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',filval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var winddirmax'
      endif
      call ck_array_real(winddirmax, recnum, filval, badflag)
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var winddirmax'
      endif
      call ck_array_real(winddirmax, recnum, misval, badflag)
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
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var windgust'
      endif
      call ck_array_real(windgust, recnum, misval, badflag)
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
       nf_status = nf_get_att_real(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var windspeed'
      endif
      call ck_array_real(windspeed, recnum, misval, badflag)
c
c     variable        netcdf long name
c      observationtime"time of observation" 
c
        nf_status = nf_inq_varid(nf_fid,'observationtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var observationtime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,observationtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ observationtime '
      endif
       nf_status = nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dfilval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var observationtime'
      endif
      do i=1,recnum
        if(observationtime(i) .eq. dfilval) observationtime(i) = badflag
      enddo !i
      nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var observationtime'
      endif
      do i=1,recnum
        if(observationtime(i) .eq. dmisval) observationtime(i) = badflag
      enddo !i
c
c     variable        netcdf long name
c      rhchangetime "relative humidity time of last change" 
c
        nf_status = nf_inq_varid(nf_fid,'rhchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rhchangetime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,rhchangetime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ rhchangetime '
      endif
       nf_status = nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dfilval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var rhchangetime'
      endif
      do i=1,recnum
         if(rhchangetime(i) .eq. dfilval) rhchangetime(i) = badflag
      enddo !i
      nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var rhchangetime'
      endif
      do i=1,recnum
         if(rhchangetime(i) .eq. dmisval) rhchangetime(i) = badflag
      enddo !i
c
c     variable        netcdf long name
c      stationpresschangetime"station press time of last change" 
c
        nf_status = nf_inq_varid(nf_fid,'stationpresschangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationpresschangetime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,
     &                                        stationpresschangetime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stationpresschangetime '
      endif
       nf_status = nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dfilval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var stationpresschangetime'
      endif
      do i=1,recnum
         if(stationpresschangetime(i) .eq. dfilval) 
     &                             stationpresschangetime(i) = badflag
      enddo !i
      nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var stationpresschangetime'
      endif
      do i=1,recnum
         if(stationpresschangetime(i) .eq. dmisval) 
     &                             stationpresschangetime(i) = badflag
      enddo !i
c
c     variable        netcdf long name
c      tempchangetime"temperature time of last change" 
c
        nf_status = nf_inq_varid(nf_fid,'tempchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tempchangetime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,tempchangetime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ tempchangetime '
      endif
       nf_status = nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dfilval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var tempchangetime'
      endif
      do i=1,recnum
         if(tempchangetime(i) .eq. dfilval) tempchangetime(i) = badflag
      enddo !i
      nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var tempchangetime'
      endif
      do i=1,recnum
         if(tempchangetime(i) .eq. dmisval) tempchangetime(i) = badflag
      enddo !i
c
c     variable        netcdf long name
c      winddirchangetime"wind direction time of last change" 
c
        nf_status = nf_inq_varid(nf_fid,'winddirchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddirchangetime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,winddirchangetime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ winddirchangetime '
      endif
       nf_status = nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dfilval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var winddirchangetime'
      endif
      do i=1,recnum
         if(winddirchangetime(i) .eq. dfilval) 
     &                                 winddirchangetime(i) = badflag
      enddo !i
      nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var winddirchangetime'
      endif
      do i=1,recnum
         if(winddirchangetime(i) .eq. dmisval) 
     &                                 winddirchangetime(i) = badflag
      enddo !i
c
c     variable        netcdf long name
c      windgustchangetime"wind gust time of last change" 
c
        nf_status = nf_inq_varid(nf_fid,'windgustchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windgustchangetime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,windgustchangetime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ windgustchangetime '
      endif
       nf_status = nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dfilval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var windgustchangetime'
      endif
      do i=1,recnum
         if(windgustchangetime(i) .eq. dfilval) 
     &                               windgustchangetime(i) = badflag
      enddo !i
      nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var windgustchangetime'
      endif
      do i=1,recnum
         if(windgustchangetime(i) .eq. dmisval) 
     &                               windgustchangetime(i) = badflag
      enddo !i
c
c     variable        netcdf long name
c      windspeedchangetime"wind speed time of last change" 
c
        nf_status = nf_inq_varid(nf_fid,'windspeedchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeedchangetime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,windspeedchangetime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ windspeedchangetime '
      endif
       nf_status = nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dfilval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var windspeedchangetime'
      endif
      do i=1,recnum
         if(windspeedchangetime(i) .eq. dfilval) 
     &                                windspeedchangetime(i) = badflag
      enddo !i
      nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. nf_noerr) then
         print *, nf_strerror(nf_status)
         print *, ' in var windspeedchangetime'
      endif
      do i=1,recnum
         if(windspeedchangetime(i) .eq. dmisval) 
     &                                windspeedchangetime(i) = badflag
      enddo !i
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif
c
      istatus = 1
c
      return
      end
