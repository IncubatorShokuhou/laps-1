c
c  subroutine to read the file "ldad automated mesonet data " 
c
      subroutine read_ldad_madis_netcdf(nf_fid, maxsensor, recnum, 
     +     maxpstentries, code1pst, code2pst, namepst,
     +     altimeterqcr, dewpointqcr, firstoverflow, globalinventory, 
     +     nstaticids, numericwmoid, precipaccumqcr, precipintensity, 
     +     preciprateqcr, preciptype, presschange3hourqcr, 
     +     presschangechar, relhumidityqcr, sealevelpressureqcr, 
     +     stationpressureqcr, temperatureqcr, visibilityqcr, 
     +     winddirqcr, windspeedqcr, altimeter, dewpoint, elevation, 
     +     latitude, longitude, meanweightedtemperature, precipaccum, 
     +     preciprate, presschange3hour, relhumidity, 
     +     sealevelpressure, seasurfacetemp, soilmoisturepercent, 
     +     soiltemperature, solarradiation, stationpressure, 
     +     temperature, visibility, winddir, winddirmax, windgust, 
     +     windspeed, altimeterdd, dataprovider, dewpointdd, 
     +     precipaccumdd, precipratedd, presweather, 
     +     presschange3hourdd, providerid, relhumiditydd, 
     +     sealevelpressuredd, stationid, stationname, 
     +     stationpressuredd, stationtype, temperaturedd, 
     +     visibilitydd, winddirdd, windspeeddd, observationtime, 
     +     receivedtime, reporttime, rhchangetime, 
     +     stationpresschangetime, tempchangetime, winddirchangetime, 
     +     windgustchangetime, windspeedchangetime,badflag)
c
      include 'netcdf.inc'
      integer maxsensor, recnum,nf_fid, nf_vid, nf_status
      integer code1pst(maxpstentries), code2pst(maxpstentries), 
     +     code4pst(maxpstentries),
     +     altimeterqcr(recnum), dewpointqcr(recnum),
     +     firstoverflow, globalinventory, nstaticids,
     +     numericwmoid(recnum), precipaccumqcr(recnum),
     +     precipintensity( maxsensor, recnum),
     +     preciprateqcr(recnum), preciptype( maxsensor, recnum),
     +     presschange3hourqcr(recnum), presschangechar(recnum),
     +     relhumidityqcr(recnum), sealevelpressureqcr(recnum),
     +     stationpressureqcr(recnum), temperatureqcr(recnum),
     +     visibilityqcr(recnum), winddirqcr(recnum),
     +     windspeedqcr(recnum)
      real altimeter(recnum), dewpoint(recnum), elevation(recnum),
     +     latitude(recnum), longitude(recnum),
     +     meanweightedtemperature(recnum), precipaccum(recnum),
     +     preciprate(recnum), presschange3hour(recnum),
     +     relhumidity(recnum), sealevelpressure(recnum),
     +     seasurfacetemp(recnum), soilmoisturepercent(recnum),
     +     soiltemperature(recnum), solarradiation(recnum),
     +     stationpressure(recnum), temperature(recnum),
     +     visibility(recnum), winddir(recnum), winddirmax(recnum),
     +     windgust(recnum), windspeed(recnum)
      double precision observationtime(recnum), receivedtime(recnum),
     +     reporttime(recnum), rhchangetime(recnum),
     +     stationpresschangetime(recnum), tempchangetime(recnum),
     +     winddirchangetime(recnum), windgustchangetime(recnum),
     +     windspeedchangetime(recnum)
      character temperaturedd(recnum)
      character precipaccumdd(recnum)
      character dewpointdd(recnum)
      character*6 stationid(recnum)
      character presschange3hourdd(recnum)
      character relhumiditydd(recnum)
      character*25 presweather(recnum)
      character*12 providerid(recnum)
      character visibilitydd(recnum)
      character winddirdd(recnum)
      character*11 dataprovider(recnum)
      character*11 namepst(maxpstentries)
      character altimeterdd(recnum)
      character*51 stationname(recnum)
      character precipratedd(recnum)
      character*11 stationtype(recnum)
      character stationpressuredd(recnum)
      character sealevelpressuredd(recnum)
      character windspeeddd(recnum)


c   variables of type real
c
c     variable        netcdf long name
c     altimeter     "altimeter setting"
c
      nf_status=nf_inq_varid(nf_fid,'altimeter',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for altimeter'
       print *,'set altimeter to badflag'
       altimeter = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,altimeter)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for altimeter'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - altimeter'
       else
        call ck_array_real(altimeter,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - altimeter'
       else
        call ck_array_real(altimeter,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     dewpoint      "dew point temperature"
c
      nf_status=nf_inq_varid(nf_fid,'dewpoint',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dewpoint'
       print *,'set dewpoint to badflag'
       dewpoint = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,dewpoint)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dewpoint'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - dewpoint'
       else
        call ck_array_real(dewpoint,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - dewpoint'
       else
        call ck_array_real(dewpoint,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     elevation     "elevation"
c
      nf_status=nf_inq_varid(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for elevation'
       print *,'set elevation to badflag'
       elevation = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,elevation)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for elevation'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - elevation'
       else
        call ck_array_real(elevation,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - elevation'
       else
        call ck_array_real(elevation,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     latitude      "latitude"
c
      nf_status=nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for latitude'
       print *,'set latitude to badflag'
       latitude = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,latitude)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for latitude'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - latitude'
       else
        call ck_array_real(latitude,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - latitude'
       else
        call ck_array_real(latitude,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     longitude     "longitude"
c
      nf_status=nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for longitude'
       print *,'set longitude to badflag'
       longitude = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,longitude)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for longitude'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - longitude'
       else
        call ck_array_real(longitude,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - longitude'
       else
        call ck_array_real(longitude,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     meanweightedtemperature"mean weighted temperature"
c
      nf_status=nf_inq_varid(nf_fid,'meanweightedtemperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for meanweightedtemperature'
       print *,'set meanweightedtemperature to badflag'
       meanweightedtemperature = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,meanweightedtemperature)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for meanweightedtemperature'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - meanweightedtemperature'
       else
        call ck_array_real(meanweightedtemperature,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - meanweightedtemperature'
       else
        call ck_array_real(meanweightedtemperature,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precipaccum   "precip accumulation"
c
      nf_status=nf_inq_varid(nf_fid,'precipaccum',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipaccum'
       print *,'set precipaccum to badflag'
       precipaccum = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precipaccum)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipaccum'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precipaccum'
       else
        call ck_array_real(precipaccum,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precipaccum'
       else
        call ck_array_real(precipaccum,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     preciprate    "precipitation rate"
c
      nf_status=nf_inq_varid(nf_fid,'preciprate',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for preciprate'
       print *,'set preciprate to badflag'
       preciprate = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,preciprate)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for preciprate'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - preciprate'
       else
        call ck_array_real(preciprate,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - preciprate'
       else
        call ck_array_real(preciprate,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     presschange3hour"3 hour pressure change"
c
      nf_status=nf_inq_varid(nf_fid,'presschange3hour',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for presschange3hour'
       print *,'set presschange3hour to badflag'
       presschange3hour = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,presschange3hour)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for presschange3hour'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - presschange3hour'
       else
        call ck_array_real(presschange3hour,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - presschange3hour'
       else
        call ck_array_real(presschange3hour,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     relhumidity   "relative humidity"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidity',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for relhumidity'
       print *,'set relhumidity to badflag'
       relhumidity = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,relhumidity)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for relhumidity'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - relhumidity'
       else
        call ck_array_real(relhumidity,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - relhumidity'
       else
        call ck_array_real(relhumidity,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     sealevelpressure"sea level pressure"
c
      nf_status=nf_inq_varid(nf_fid,'sealevelpressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for sealevelpressure'
       print *,'set sealevelpressure to badflag'
       sealevelpressure = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,sealevelpressure)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for sealevelpressure'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - sealevelpressure'
       else
        call ck_array_real(sealevelpressure,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - sealevelpressure'
       else
        call ck_array_real(sealevelpressure,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     seasurfacetemp"sea surface temperature"
c
      nf_status=nf_inq_varid(nf_fid,'seasurfacetemp',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for seasurfacetemp'
       print *,'set seasurfacetemp to badflag'
       seasurfacetemp = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,seasurfacetemp)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for seasurfacetemp'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - seasurfacetemp'
       else
        call ck_array_real(seasurfacetemp,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - seasurfacetemp'
       else
        call ck_array_real(seasurfacetemp,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     soilmoisturepercent  "soil moisture"
c
      nf_status=nf_inq_varid(nf_fid,'soilmoisturepercent',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for soilmoisturepercent'
       print *,'set soilmoisturepercent to badflag'
       soilmoisturepercent = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,soilmoisturepercent)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for soilmoisturepercent'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - soilmoisturepercent'
       else
        call ck_array_real(soilmoisturepercent,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - soilmoisturepercent'
       else
        call ck_array_real(soilmoisturepercent,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     soiltemperature"soil temperature"
c
      nf_status=nf_inq_varid(nf_fid,'soiltemperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for soiltemperature'
       print *,'set soiltemperature to badflag'
       soiltemperature = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,soiltemperature)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for soiltemperature'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - soiltemperature'
       else
        call ck_array_real(soiltemperature,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - soiltemperature'
       else
        call ck_array_real(soiltemperature,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     solarradiation"solar radiation"
c
      nf_status=nf_inq_varid(nf_fid,'solarradiation',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for solarradiation'
       print *,'set solarradiation to badflag'
       solarradiation = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,solarradiation)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for solarradiation'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - solarradiation'
       else
        call ck_array_real(solarradiation,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - solarradiation'
       else
        call ck_array_real(solarradiation,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     stationpressure"station pressure"
c
      nf_status=nf_inq_varid(nf_fid,'stationpressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationpressure'
       print *,'set stationpressure to badflag'
       stationpressure = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,stationpressure)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationpressure'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - stationpressure'
       else
        call ck_array_real(stationpressure,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - stationpressure'
       else
        call ck_array_real(stationpressure,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     temperature   "temperature"
c
      nf_status=nf_inq_varid(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for temperature'
       print *,'set temperature to badflag'
       temperature = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,temperature)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for temperature'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - temperature'
       else
        call ck_array_real(temperature,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - temperature'
       else
        call ck_array_real(temperature,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     visibility    "visibility"
c
      nf_status=nf_inq_varid(nf_fid,'visibility',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for visibility'
       print *,'set visibility to badflag'
       visibility = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,visibility)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for visibility'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - visibility'
       else
        call ck_array_real(visibility,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - visibility'
       else
        call ck_array_real(visibility,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     winddir       "wind direction"
c
      nf_status=nf_inq_varid(nf_fid,'winddir',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for winddir'
       print *,'set winddir to badflag'
       winddir = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,winddir)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for winddir'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - winddir'
       else
        call ck_array_real(winddir,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - winddir'
       else
        call ck_array_real(winddir,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     winddirmax    "wind direction at gust"
c
      nf_status=nf_inq_varid(nf_fid,'winddirmax',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for winddirmax'
       print *,'set winddirmax to badflag'
       winddirmax = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,winddirmax)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for winddirmax'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - winddirmax'
       else
        call ck_array_real(winddirmax,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - winddirmax'
       else
        call ck_array_real(winddirmax,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     windgust      "wind gust"
c
      nf_status=nf_inq_varid(nf_fid,'windgust',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for windgust'
       print *,'set windgust to badflag'
       windgust = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,windgust)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for windgust'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - windgust'
       else
        call ck_array_real(windgust,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - windgust'
       else
        call ck_array_real(windgust,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     windspeed     "wind speed"
c
      nf_status=nf_inq_varid(nf_fid,'windspeed',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for windspeed'
       print *,'set windspeed to badflag'
       windspeed = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,windspeed)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for windspeed'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - windspeed'
       else
        call ck_array_real(windspeed,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - windspeed'
       else
        call ck_array_real(windspeed,recnum,valmis
     1                    ,badflag)
       endif
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c     altimeterqcr  "altimeter setting qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'altimeterqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for altimeterqcr'
       print *,'set altimeterqcr to -99'
       altimeterqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,altimeterqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for altimeterqcr'
       endif
      endif
c
c     variable        netcdf long name
c     dewpointqcr   "dew point qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'dewpointqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dewpointqcr'
       print *,'set dewpointqcr to -99'
       dewpointqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,dewpointqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dewpointqcr'
       endif
      endif
c
c     variable        netcdf long name
c     firstoverflow 
c
      nf_status=nf_inq_varid(nf_fid,'firstoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for firstoverflow'
       print *,'set firstoverflow to -99'
       firstoverflow = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,firstoverflow)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for firstoverflow'
       endif
      endif
c
c     variable        netcdf long name
c     globalinventory
c
      nf_status=nf_inq_varid(nf_fid,'globalinventory',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for globalinventory'
       print *,'set globalinventory to -99'
       globalinventory = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,globalinventory)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for globalinventory'
       endif
      endif
c
c     variable        netcdf long name
c     nstaticids    
c
      nf_status=nf_inq_varid(nf_fid,'nstaticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for nstaticids'
       print *,'set nstaticids to -99'
       nstaticids = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,nstaticids)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for nstaticids'
       endif
      endif
c
c     variable        netcdf long name
c     numericwmoid  "numeric wmo identification"
c
      nf_status=nf_inq_varid(nf_fid,'numericwmoid',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for numericwmoid'
       print *,'set numericwmoid to -99'
       numericwmoid = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,numericwmoid)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for numericwmoid'
       endif
      endif
c
c     variable        netcdf long name
c     precipaccumqcr"precip amount qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precipaccumqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipaccumqcr'
       print *,'set precipaccumqcr to -99'
       precipaccumqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precipaccumqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipaccumqcr'
       endif
      endif
c
c     variable        netcdf long name
c     precipintensity"precipitation intensity"
c
      nf_status=nf_inq_varid(nf_fid,'precipintensity',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipintensity'
       print *,'set precipintensity to -99'
       precipintensity = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precipintensity)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipintensity'
       endif
      endif
c
c     variable        netcdf long name
c     preciprateqcr "precip rate qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'preciprateqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for preciprateqcr'
       print *,'set preciprateqcr to -99'
       preciprateqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,preciprateqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for preciprateqcr'
       endif
      endif
c
c     variable        netcdf long name
c     preciptype    "precipitation type"
c
      nf_status=nf_inq_varid(nf_fid,'preciptype',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for preciptype'
       print *,'set preciptype to -99'
       preciptype = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,preciptype)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for preciptype'
       endif
      endif
c
c     variable        netcdf long name
c     presschange3hourqcr"3h pressure change qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'presschange3hourqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for presschange3hourqcr'
       print *,'set presschange3hourqcr to -99'
       presschange3hourqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,presschange3hourqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for presschange3hourqcr'
       endif
      endif
c
c     variable        netcdf long name
c     presschangechar"character of pressure change"
c
      nf_status=nf_inq_varid(nf_fid,'presschangechar',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for presschangechar'
       print *,'set presschangechar to -99'
       presschangechar = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,presschangechar)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for presschangechar'
       endif
      endif
c
c     variable        netcdf long name
c     relhumidityqcr"relative humidity qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidityqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for relhumidityqcr'
       print *,'set relhumidityqcr to -99'
       relhumidityqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,relhumidityqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for relhumidityqcr'
       endif
      endif
c
c     variable        netcdf long name
c     sealevelpressureqcr"sea level pressure qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'sealevelpressureqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for sealevelpressureqcr'
       print *,'set sealevelpressureqcr to -99'
       sealevelpressureqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,sealevelpressureqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for sealevelpressureqcr'
       endif
      endif
c
c     variable        netcdf long name
c     stationpressureqcr"station pressure qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'stationpressureqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationpressureqcr'
       print *,'set stationpressureqcr to -99'
       stationpressureqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,stationpressureqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationpressureqcr'
       endif
      endif
c
c     variable        netcdf long name
c     temperatureqcr"temperature qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'temperatureqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for temperatureqcr'
       print *,'set temperatureqcr to -99'
       temperatureqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,temperatureqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for temperatureqcr'
       endif
      endif
c
c     variable        netcdf long name
c     visibilityqcr "visibility qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'visibilityqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for visibilityqcr'
       print *,'set visibilityqcr to -99'
       visibilityqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,visibilityqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for visibilityqcr'
       endif
      endif
c
c     variable        netcdf long name
c     winddirqcr    "wind direction qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'winddirqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for winddirqcr'
       print *,'set winddirqcr to -99'
       winddirqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,winddirqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for winddirqcr'
       endif
      endif
c
c     variable        netcdf long name
c     windspeedqcr  "wind speed qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'windspeedqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for windspeedqcr'
       print *,'set windspeedqcr to -99'
       windspeedqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,windspeedqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for windspeedqcr'
       endif
      endif
c
c     variable        netcdf long name
c     code1pst      "precip variable definition"
c
      nf_status=nf_inq_varid(nf_fid,'code1pst',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for code1pst'
       print *,'set code1pst to -99'
       code1pst = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,code1pst)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for code1pst'
       endif
      endif

c
c     variable        netcdf long name
c     code2pst      "solarradiation variable definition"
c
      nf_status=nf_inq_varid(nf_fid,'code2pst',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for code2pst'
       print *,'set code2pst to -99'
       code2pst = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,code2pst)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for code2pst'
       endif
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c     observationtime"time of observation"
c
      nf_status=nf_inq_varid(nf_fid,'observationtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for observationtime'
       print *,'set observationtime to badflag'
       observationtime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,observationtime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for observationtime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - observationtime'
       else
        call ck_array_dble(observationtime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - observationtime'
       else
        call ck_array_dble(observationtime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     receivedtime  "time data was received"
c
      nf_status=nf_inq_varid(nf_fid,'receivedtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for receivedtime'
       print *,'set receivedtime to badflag'
       receivedtime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,receivedtime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for receivedtime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - receivedtime'
       else
        call ck_array_dble(receivedtime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - receivedtime'
       else
        call ck_array_dble(receivedtime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     reporttime    "time data was processed by the provider"
c
      nf_status=nf_inq_varid(nf_fid,'reporttime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for reporttime'
       print *,'set reporttime to badflag'
       reporttime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,reporttime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for reporttime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - reporttime'
       else
        call ck_array_dble(reporttime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - reporttime'
       else
        call ck_array_dble(reporttime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     rhchangetime  "relative humidity time of last change"
c
      nf_status=nf_inq_varid(nf_fid,'rhchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for rhchangetime'
       print *,'set rhchangetime to badflag'
       rhchangetime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,rhchangetime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for rhchangetime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - rhchangetime'
       else
        call ck_array_dble(rhchangetime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - rhchangetime'
       else
        call ck_array_dble(rhchangetime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     stationpresschangetime"station press time of last change"
c
      nf_status=nf_inq_varid(nf_fid,'stationpresschangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationpresschangetime'
       print *,'set stationpresschangetime to badflag'
       stationpresschangetime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,stationpresschangetime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationpresschangetime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - stationpresschangetime'
       else
        call ck_array_dble(stationpresschangetime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - stationpresschangetime'
       else
        call ck_array_dble(stationpresschangetime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     tempchangetime"temperature time of last change"
c
      nf_status=nf_inq_varid(nf_fid,'tempchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for tempchangetime'
       print *,'set tempchangetime to badflag'
       tempchangetime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,tempchangetime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for tempchangetime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - tempchangetime'
       else
        call ck_array_dble(tempchangetime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - tempchangetime'
       else
        call ck_array_dble(tempchangetime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     winddirchangetime"wind direction time of last change"
c
      nf_status=nf_inq_varid(nf_fid,'winddirchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for winddirchangetime'
       print *,'set winddirchangetime to badflag'
       winddirchangetime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,winddirchangetime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for winddirchangetime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - winddirchangetime'
       else
        call ck_array_dble(winddirchangetime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - winddirchangetime'
       else
        call ck_array_dble(winddirchangetime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     windgustchangetime"wind gust time of last change"
c
      nf_status=nf_inq_varid(nf_fid,'windgustchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for windgustchangetime'
       print *,'set windgustchangetime to badflag'
       windgustchangetime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,windgustchangetime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for windgustchangetime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - windgustchangetime'
       else
        call ck_array_dble(windgustchangetime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - windgustchangetime'
       else
        call ck_array_dble(windgustchangetime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     windspeedchangetime"wind speed time of last change"
c
      nf_status=nf_inq_varid(nf_fid,'windspeedchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for windspeedchangetime'
       print *,'set windspeedchangetime to badflag'
       windspeedchangetime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,windspeedchangetime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for windspeedchangetime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - windspeedchangetime'
       else
        call ck_array_dble(windspeedchangetime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - windspeedchangetime'
       else
        call ck_array_dble(windspeedchangetime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c     altimeterdd   "altimeter setting qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'altimeterdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for altimeterdd'
       print *,'set altimeterdd to " "'
       altimeterdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,altimeterdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for altimeterdd'
       endif
      endif
c
c     variable        netcdf long name
c     dataprovider  "local data provider"
c
      nf_status=nf_inq_varid(nf_fid,'dataprovider',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dataprovider'
       print *,'set dataprovider to " "'
       dataprovider = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,dataprovider)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dataprovider'
       endif
      endif
c
c     variable        netcdf long name
c     dewpointdd    "dew point qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'dewpointdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dewpointdd'
       print *,'set dewpointdd to " "'
       dewpointdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,dewpointdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dewpointdd'
       endif
      endif
c
c     variable        netcdf long name
c     namepst       "pst provider or subprovider name"
c
      nf_status=nf_inq_varid(nf_fid,'namepst',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for namepst'
       print *,'set namepst to " "'
       namepst = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,namepst)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for namepst'
       endif
      endif
c
c     variable        netcdf long name
c     precipaccumdd "precip amount qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'precipaccumdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipaccumdd'
       print *,'set precipaccumdd to " "'
       precipaccumdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,precipaccumdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipaccumdd'
       endif
      endif
c
c     variable        netcdf long name
c     precipratedd  "precip rate qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'precipratedd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipratedd'
       print *,'set precipratedd to " "'
       precipratedd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,precipratedd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipratedd'
       endif
      endif
c
c     variable        netcdf long name
c     presweather   "present weather"
c
      nf_status=nf_inq_varid(nf_fid,'presweather',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for presweather'
       print *,'set presweather to " "'
       presweather = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,presweather)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for presweather'
       endif
      endif
c
c     variable        netcdf long name
c     presschange3hourdd"3h pressure change qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'presschange3hourdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for presschange3hourdd'
       print *,'set presschange3hourdd to " "'
       presschange3hourdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,presschange3hourdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for presschange3hourdd'
       endif
      endif
c
c     variable        netcdf long name
c     providerid    "data provider station id"
c
      nf_status=nf_inq_varid(nf_fid,'providerid',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for providerid'
       print *,'set providerid to " "'
       providerid = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,providerid)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for providerid'
       endif
      endif
c
c     variable        netcdf long name
c     relhumiditydd "relative humidity qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'relhumiditydd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for relhumiditydd'
       print *,'set relhumiditydd to " "'
       relhumiditydd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,relhumiditydd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for relhumiditydd'
       endif
      endif
c
c     variable        netcdf long name
c     sealevelpressuredd"sea level pressure qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'sealevelpressuredd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for sealevelpressuredd'
       print *,'set sealevelpressuredd to " "'
       sealevelpressuredd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,sealevelpressuredd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for sealevelpressuredd'
       endif
      endif
c
c     variable        netcdf long name
c     stationid     "alphanumeric station id"
c
      nf_status=nf_inq_varid(nf_fid,'stationid',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationid'
       print *,'set stationid to " "'
       stationid = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,stationid)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationid'
       endif
      endif
c
c     variable        netcdf long name
c     stationname   "alphanumeric station name"
c
      nf_status=nf_inq_varid(nf_fid,'stationname',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationname'
       print *,'set stationname to " "'
       stationname = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,stationname)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationname'
       endif
      endif
c
c     variable        netcdf long name
c     stationpressuredd"station pressure qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'stationpressuredd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationpressuredd'
       print *,'set stationpressuredd to " "'
       stationpressuredd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,stationpressuredd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationpressuredd'
       endif
      endif
c
c     variable        netcdf long name
c     stationtype   "ldad station type"
c
      nf_status=nf_inq_varid(nf_fid,'stationtype',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationtype'
       print *,'set stationtype to " "'
       stationtype = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,stationtype)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationtype'
       endif
      endif
c
c     variable        netcdf long name
c     temperaturedd "temperature qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'temperaturedd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for temperaturedd'
       print *,'set temperaturedd to " "'
       temperaturedd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,temperaturedd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for temperaturedd'
       endif
      endif
c
c     variable        netcdf long name
c     visibilitydd  "visibility qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'visibilitydd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for visibilitydd'
       print *,'set visibilitydd to " "'
       visibilitydd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,visibilitydd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for visibilitydd'
       endif
      endif
c
c     variable        netcdf long name
c     winddirdd     "wind direction qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'winddirdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for winddirdd'
       print *,'set winddirdd to " "'
       winddirdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,winddirdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for winddirdd'
       endif
      endif
c
c     variable        netcdf long name
c     windspeeddd   "wind speed qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'windspeeddd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for windspeeddd'
       print *,'set windspeeddd to " "'
       windspeeddd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,windspeeddd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for windspeeddd'
       endif
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
