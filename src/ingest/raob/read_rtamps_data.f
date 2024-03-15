c
c  subroutine to read the file 
c
      subroutine read_rtamps_netcdf_dum(nfid, manlevel, maxstaticids,        
     +     ninventorybins, rawlevel, recnum, stdlevel, termlevel, 
     +     troplevel, editflag, firstinbin, firstoverflow, 
     +     globalinventory, indxrefr, invtime, inventory, irman, 
     +     isoverflow, lastinbin, lastrecord, nstaticids, oiman, 
     +     optindxrefr, prevrecord, storedobs, abshumidity, 
     +     airdensity, barompressure, bpman, bpstd, bpterm, bptrop, 
     +     dewpt, direction, dpman, dpstd, dptrop, drman, drstd, 
     +     elevation, geomheight, geopheight, ghman, ghterm, ghtrop, 
     +     gpstd, gpterm, latitude, longitude, precipwater, 
     +     relhumidity, rhman, rhstd, riserate, shear, sheardir, 
     +     shearmagx, shearmagy, spman, spstd, speed, temperature, 
     +     tpman, tpstd, tptrop, vaporpressure, velerror, velsound, 
     +     observationtime, receivedtime, reporttime, dataprovider, 
     +     providerid, staticids, stationname)
c
      include 'netcdf.inc'
      integer manlevel, maxstaticids, ninventorybins, rawlevel, 
     +     recnum, stdlevel, termlevel, troplevel,nfid, nf_vid, 
     +     nf_status
      integer editflag( rawlevel, recnum), firstinbin(ninventorybins),
     +     firstoverflow, globalinventory, indxrefr( rawlevel,
     +     recnum), invtime(recnum), inventory(maxstaticids), irman(
     +     manlevel, recnum), isoverflow(recnum),
     +     lastinbin(ninventorybins), lastrecord(maxstaticids),
     +     nstaticids, oiman( manlevel, recnum), optindxrefr(
     +     rawlevel, recnum), prevrecord(recnum), storedobs(recnum)
      real abshumidity( rawlevel, recnum), airdensity( rawlevel,
     +     recnum), barompressure( rawlevel, recnum), bpman(
     +     manlevel, recnum), bpstd( stdlevel, recnum), bpterm(
     +     termlevel, recnum), bptrop( troplevel, recnum), dewpt(
     +     rawlevel, recnum), direction( rawlevel, recnum), dpman(
     +     manlevel, recnum), dpstd( stdlevel, recnum), dptrop(
     +     troplevel, recnum), drman( manlevel, recnum), drstd(
     +     stdlevel, recnum), elevation(recnum), geomheight(
     +     rawlevel, recnum), geopheight( rawlevel, recnum), ghman(
     +     manlevel, recnum), ghterm( termlevel, recnum), ghtrop(
     +     troplevel, recnum), gpstd( stdlevel, recnum), gpterm(
     +     termlevel, recnum), latitude(recnum), longitude(recnum),
     +     precipwater( rawlevel, recnum), relhumidity( rawlevel,
     +     recnum), rhman( manlevel, recnum), rhstd( stdlevel,
     +     recnum), riserate( rawlevel, recnum), shear( rawlevel,
     +     recnum), sheardir( rawlevel, recnum), shearmagx( rawlevel,
     +     recnum), shearmagy( rawlevel, recnum), spman( manlevel,
     +     recnum), spstd( stdlevel, recnum), speed( rawlevel,
     +     recnum), temperature( rawlevel, recnum), tpman( manlevel,
     +     recnum), tpstd( stdlevel, recnum), tptrop( troplevel,
     +     recnum), vaporpressure( rawlevel, recnum), velerror(
     +     rawlevel, recnum), velsound( rawlevel, recnum)
      double precision observationtime(recnum), receivedtime(recnum),
     +     reporttime(recnum)
      character*12 providerid(recnum)
      character*24 staticids(maxstaticids)
      character*51 stationname(recnum)
      character*11 dataprovider(recnum)


c   variables of type real
c
c     variable        netcdf long name
c      abshumidity  "absolute humidity"
c
        nf_status = nf_inq_varid(nfid,'abshumidity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var abshumidity'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,abshumidity)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var abshumidity'
      endif
c
c     variable        netcdf long name
c      airdensity   "density of air"
c
        nf_status = nf_inq_varid(nfid,'airdensity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var airdensity'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,airdensity)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var airdensity'
      endif
c
c     variable        netcdf long name
c      barompressure"pressure"
c
        nf_status = nf_inq_varid(nfid,'barompressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var barompressure'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,barompressure)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var barompressure'
      endif
c
c     variable        netcdf long name
c      bpman        "pressure - mandatory level"
c
        nf_status = nf_inq_varid(nfid,'bpman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpman'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,bpman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpman'
      endif
c
c     variable        netcdf long name
c      bpstd        "pressure - standard level"
c
        nf_status = nf_inq_varid(nfid,'bpstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpstd'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,bpstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpstd'
      endif
c
c     variable        netcdf long name
c      bpterm       "pressure - termination"
c
        nf_status = nf_inq_varid(nfid,'bpterm',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpterm'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,bpterm)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpterm'
      endif
c
c     variable        netcdf long name
c      bptrop       "pressure - tropopause level"
c
        nf_status = nf_inq_varid(nfid,'bptrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bptrop'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,bptrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bptrop'
      endif
c
c     variable        netcdf long name
c      dewpt        "dew point temperature"
c
        nf_status = nf_inq_varid(nfid,'dewpt',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dewpt'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,dewpt)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dewpt'
      endif
c
c     variable        netcdf long name
c      direction    "wind direction"
c
        nf_status = nf_inq_varid(nfid,'direction',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var direction'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,direction)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var direction'
      endif
c
c     variable        netcdf long name
c      dpman        "dew point temperature - mandatory level"
c
        nf_status = nf_inq_varid(nfid,'dpman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dpman'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,dpman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dpman'
      endif
c
c     variable        netcdf long name
c      dpstd        "dew point temperature - standard level"
c
        nf_status = nf_inq_varid(nfid,'dpstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dpstd'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,dpstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dpstd'
      endif
c
c     variable        netcdf long name
c      dptrop       "dew point temperature - tropopause level"
c
        nf_status = nf_inq_varid(nfid,'dptrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dptrop'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,dptrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dptrop'
      endif
c
c     variable        netcdf long name
c      drman        "wind direction - mandatory level"
c
        nf_status = nf_inq_varid(nfid,'drman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var drman'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,drman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var drman'
      endif
c
c     variable        netcdf long name
c      drstd        "wind direction - standard level"
c
        nf_status = nf_inq_varid(nfid,'drstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var drstd'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,drstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var drstd'
      endif
c
c     variable        netcdf long name
c      elevation    "station elevation"
c
        nf_status = nf_inq_varid(nfid,'elevation',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevation'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,elevation)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevation'
      endif
c
c     variable        netcdf long name
c      geomheight   "geometric height"
c
        nf_status = nf_inq_varid(nfid,'geomheight',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var geomheight'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,geomheight)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var geomheight'
      endif
c
c     variable        netcdf long name
c      geopheight   "geopotential height"
c
        nf_status = nf_inq_varid(nfid,'geopheight',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var geopheight'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,geopheight)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var geopheight'
      endif
c
c     variable        netcdf long name
c      ghman        "geometric - mandatory level"
c
        nf_status = nf_inq_varid(nfid,'ghman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghman'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,ghman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghman'
      endif
c
c     variable        netcdf long name
c      ghterm       "geometric - termination"
c
        nf_status = nf_inq_varid(nfid,'ghterm',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghterm'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,ghterm)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghterm'
      endif
c
c     variable        netcdf long name
c      ghtrop       "geometric - tropopause level"
c
        nf_status = nf_inq_varid(nfid,'ghtrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghtrop'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,ghtrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghtrop'
      endif
c
c     variable        netcdf long name
c      gpstd        "geopotential - standard level"
c
        nf_status = nf_inq_varid(nfid,'gpstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gpstd'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,gpstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gpstd'
      endif
c
c     variable        netcdf long name
c      gpterm       "geopotential - termination"
c
        nf_status = nf_inq_varid(nfid,'gpterm',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gpterm'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,gpterm)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gpterm'
      endif
c
c     variable        netcdf long name
c      latitude     "station latitude"
c
        nf_status = nf_inq_varid(nfid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,latitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
c
c     variable        netcdf long name
c      longitude    "station longitude"
c
        nf_status = nf_inq_varid(nfid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,longitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
c
c     variable        netcdf long name
c      precipwater  "precipitable water"
c
        nf_status = nf_inq_varid(nfid,'precipwater',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var precipwater'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,precipwater)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var precipwater'
      endif
c
c     variable        netcdf long name
c      relhumidity  "relative humidity"
c
        nf_status = nf_inq_varid(nfid,'relhumidity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidity'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,relhumidity)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidity'
      endif
c
c     variable        netcdf long name
c      rhman        "relative humidity - mandatory level"
c
        nf_status = nf_inq_varid(nfid,'rhman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rhman'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,rhman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rhman'
      endif
c
c     variable        netcdf long name
c      rhstd        "relative humidity - standard level"
c
        nf_status = nf_inq_varid(nfid,'rhstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rhstd'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,rhstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rhstd'
      endif
c
c     variable        netcdf long name
c      riserate     "rise rate"
c
        nf_status = nf_inq_varid(nfid,'riserate',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var riserate'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,riserate)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var riserate'
      endif
c
c     variable        netcdf long name
c      shear        "shear"
c
        nf_status = nf_inq_varid(nfid,'shear',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shear'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,shear)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shear'
      endif
c
c     variable        netcdf long name
c      sheardir     "shear direction"
c
        nf_status = nf_inq_varid(nfid,'sheardir',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sheardir'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,sheardir)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sheardir'
      endif
c
c     variable        netcdf long name
c      shearmagx    "shear magnitude x-direction"
c
        nf_status = nf_inq_varid(nfid,'shearmagx',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shearmagx'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,shearmagx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shearmagx'
      endif
c
c     variable        netcdf long name
c      shearmagy    "shear magnitude y-direction"
c
        nf_status = nf_inq_varid(nfid,'shearmagy',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shearmagy'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,shearmagy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shearmagy'
      endif
c
c     variable        netcdf long name
c      spman        "wind speed - mandatory level"
c
        nf_status = nf_inq_varid(nfid,'spman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var spman'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,spman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var spman'
      endif
c
c     variable        netcdf long name
c      spstd        "wind speed - standard level"
c
        nf_status = nf_inq_varid(nfid,'spstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var spstd'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,spstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var spstd'
      endif
c
c     variable        netcdf long name
c      speed        "wind speed"
c
        nf_status = nf_inq_varid(nfid,'speed',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var speed'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,speed)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var speed'
      endif
c
c     variable        netcdf long name
c      temperature  "temperature"
c
        nf_status = nf_inq_varid(nfid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,temperature)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
c
c     variable        netcdf long name
c      tpman        "temperature - mandatory level"
c
        nf_status = nf_inq_varid(nfid,'tpman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpman'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,tpman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpman'
      endif
c
c     variable        netcdf long name
c      tpstd        "temperature - standard level"
c
        nf_status = nf_inq_varid(nfid,'tpstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpstd'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,tpstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpstd'
      endif
c
c     variable        netcdf long name
c      tptrop       "temperature - tropopause level"
c
        nf_status = nf_inq_varid(nfid,'tptrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tptrop'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,tptrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tptrop'
      endif
c
c     variable        netcdf long name
c      vaporpressure"pressure"
c
        nf_status = nf_inq_varid(nfid,'vaporpressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vaporpressure'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,vaporpressure)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vaporpressure'
      endif
c
c     variable        netcdf long name
c      velerror     "velocity error"
c
        nf_status = nf_inq_varid(nfid,'velerror',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var velerror'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,velerror)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var velerror'
      endif
c
c     variable        netcdf long name
c      velsound     "velocity of sound"
c
        nf_status = nf_inq_varid(nfid,'velsound',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var velsound'
      endif
        nf_status = nf_get_var_real(nfid,nf_vid,velsound)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var velsound'
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      editflag     "edit flag"
c
        nf_status = nf_inq_varid(nfid,'editflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var editflag'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,editflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var editflag'
      endif
c
c     variable        netcdf long name
c      firstinbin   
c
        nf_status = nf_inq_varid(nfid,'firstinbin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstinbin'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,firstinbin)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstinbin'
      endif
c
c     variable        netcdf long name
c      firstoverflow
c
        nf_status = nf_inq_varid(nfid,'firstoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstoverflow'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,firstoverflow)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstoverflow'
      endif
c
c     variable        netcdf long name
c      globalinventory
c
        nf_status = nf_inq_varid(nfid,'globalinventory',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var globalinventory'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,globalinventory)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var globalinventory'
      endif
c
c     variable        netcdf long name
c      indxrefr     "microwave index of refraction"
c
        nf_status = nf_inq_varid(nfid,'indxrefr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var indxrefr'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,indxrefr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var indxrefr'
      endif
c
c     variable        netcdf long name
c      invtime      
c
        nf_status = nf_inq_varid(nfid,'invtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var invtime'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,invtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var invtime'
      endif
c
c     variable        netcdf long name
c      inventory    
c
        nf_status = nf_inq_varid(nfid,'inventory',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var inventory'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,inventory)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var inventory'
      endif
c
c     variable        netcdf long name
c      irman        "microwave index of refraction - mandatory level"
c
        nf_status = nf_inq_varid(nfid,'irman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var irman'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,irman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var irman'
      endif
c
c     variable        netcdf long name
c      isoverflow   
c
        nf_status = nf_inq_varid(nfid,'isoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var isoverflow'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,isoverflow)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var isoverflow'
      endif
c
c     variable        netcdf long name
c      lastinbin    
c
        nf_status = nf_inq_varid(nfid,'lastinbin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lastinbin'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,lastinbin)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lastinbin'
      endif
c
c     variable        netcdf long name
c      lastrecord   
c
        nf_status = nf_inq_varid(nfid,'lastrecord',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lastrecord'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,lastrecord)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lastrecord'
      endif
c
c     variable        netcdf long name
c      nstaticids   
c
        nf_status = nf_inq_varid(nfid,'nstaticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nstaticids'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,nstaticids)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nstaticids'
      endif
c
c     variable        netcdf long name
c      oiman        "optical index of refraction - mandatory level"
c
        nf_status = nf_inq_varid(nfid,'oiman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var oiman'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,oiman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var oiman'
      endif
c
c     variable        netcdf long name
c      optindxrefr  "optical index of refraction"
c
        nf_status = nf_inq_varid(nfid,'optindxrefr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var optindxrefr'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,optindxrefr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var optindxrefr'
      endif
c
c     variable        netcdf long name
c      prevrecord   
c
        nf_status = nf_inq_varid(nfid,'prevrecord',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prevrecord'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,prevrecord)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prevrecord'
      endif
c
c     variable        netcdf long name
c      storedobs    "stored \'raw\' profile observations"
c
        nf_status = nf_inq_varid(nfid,'storedobs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var storedobs'
      endif
        nf_status = nf_get_var_int(nfid,nf_vid,storedobs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var storedobs'
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c      observationtime"observation time"
c
        nf_status = nf_inq_varid(nfid,'observationtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var observationtime'
      endif
        nf_status = nf_get_var_double(nfid,nf_vid,observationtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var observationtime'
      endif
c
c     variable        netcdf long name
c      receivedtime "received time"
c
        nf_status = nf_inq_varid(nfid,'receivedtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var receivedtime'
      endif
        nf_status = nf_get_var_double(nfid,nf_vid,receivedtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var receivedtime'
      endif
c
c     variable        netcdf long name
c      reporttime   "report time"
c
        nf_status = nf_inq_varid(nfid,'reporttime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reporttime'
      endif
        nf_status = nf_get_var_double(nfid,nf_vid,reporttime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reporttime'
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c      dataprovider "local data provider"
c
        nf_status = nf_inq_varid(nfid,'dataprovider',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dataprovider'
      endif
        nf_status = nf_get_var_text(nfid,nf_vid,dataprovider)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dataprovider'
      endif
c
c     variable        netcdf long name
c      providerid   "data provider station id"
c
        nf_status = nf_inq_varid(nfid,'providerid',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var providerid'
      endif
        nf_status = nf_get_var_text(nfid,nf_vid,providerid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var providerid'
      endif
c
c     variable        netcdf long name
c      staticids    
c
        nf_status = nf_inq_varid(nfid,'staticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staticids'
      endif
        nf_status = nf_get_var_text(nfid,nf_vid,staticids)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staticids'
      endif
c
c     variable        netcdf long name
c      stationname  "alphanumeric station name"
c
        nf_status = nf_inq_varid(nfid,'stationname',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationname'
      endif
        nf_status = nf_get_var_text(nfid,nf_vid,stationname)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationname'
      endif

      nf_status = nf_close(nfid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
