      subroutine get_wisdom_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l
     1                   ,i4time_earliest          
     1                   ,i4time_latest           
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer qcchecknum, maxstaticids, ninventorybins, recnum,nf_fid,
     +     nf_vid, nf_status
c
c  open netcdf file for reading
c
      nf_status=nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),filename
        istatus=0
        return
      endif
c
c  fill all dimension values
c
c
c get size of qcchecknum
c
      nf_status=nf_inq_dimid(nf_fid,'qcchecknum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim qcchecknum'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,qcchecknum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim qcchecknum'
      endif
c
c get size of maxstaticids
c
      nf_status=nf_inq_dimid(nf_fid,'maxstaticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxstaticids'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,maxstaticids)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxstaticids'
      endif
c
c get size of ninventorybins
c
      nf_status=nf_inq_dimid(nf_fid,'ninventorybins',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim ninventorybins'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,ninventorybins)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim ninventorybins'
      endif
c
c get size of recnum
c
      nf_status=nf_inq_dimid(nf_fid,'recnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,recnum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
      endif
      call read_wisdom_data(nf_fid, qcchecknum, maxstaticids,
     +     ninventorybins, recnum, i4time_sys, ilaps_cycle_time,
     +     nx_l, ny_l, i4time_earliest, i4time_latest, lun_out,
     +     istatus)

      return
      end
c
c
      subroutine read_wisdom_data(nf_fid, qcchecknum, maxstaticids,
     +     ninventorybins, recnum, i4time_sys, ilaps_cycle_time,
     +     nx_l, ny_l, i4time_earliest, i4time_latest, lun_out,
     +     istatus)


      include 'netcdf.inc'
      integer qcchecknum, maxstaticids, ninventorybins, recnum,nf_fid,
     +     nf_vid, nf_status
      integer cutdownflag(recnum), firstinbin(ninventorybins),
     +     firstoverflow, globalinventory, invtime(recnum),
     +     inventory(maxstaticids), isoverflow(recnum),
     +     lastinbin(ninventorybins), lastrecord(maxstaticids),
     +     latitudeqca(recnum), latitudeqcr(recnum),
     +     longitudeqca(recnum), longitudeqcr(recnum), nstaticids,
     +     pressureqca(recnum), pressureqcr(recnum),
     +     prevrecord(recnum), relhumidityqca(recnum),
     +     relhumidityqcr(recnum), secondsstage1_2(recnum),
     +     temperatureqca(recnum), temperatureqcr(recnum),
     +     winddirqca(recnum), winddirqcr(recnum),
     +     windspeedqca(recnum), windspeedqcr(recnum)
      real circularerror(recnum), elevation(recnum), latitude(recnum),
     +     latitudeqcd( qcchecknum, recnum), longitude(recnum),
     +     longitudeqcd( qcchecknum, recnum), mobileelev(recnum),
     +     mobilelat(recnum), mobilelon(recnum), pressure(recnum),
     +     pressureqcd( qcchecknum, recnum), relhumidity(recnum),
     +     relhumidityqcd( qcchecknum, recnum), temperature(recnum),
     +     temperatureqcd( qcchecknum, recnum), winddir(recnum),
     +     winddirqcd( qcchecknum, recnum), windspeed(recnum),
     +     windspeedqcd( qcchecknum, recnum)
      double precision observationtime(recnum), receivedtime(recnum),
     +     reporttime(recnum)
      character pressuredd(recnum)
      character*11 dataprovider(recnum)
      character longitudedd(recnum)
      character winddirdd(recnum)
      character*60 qct(qcchecknum)
      character*6 handbook5id(recnum)
      character*51 stationname(recnum)
      character*11 stationtype(recnum)
      character latitudedd(recnum)
      character relhumiditydd(recnum)
      character*512 rawmessage(recnum)
      character temperaturedd(recnum)
      character*24 staticids(maxstaticids)
      character windspeeddd(recnum)
      character*12 providerid(recnum)
      character*6 stationid(recnum)
      character*4 homewfo(recnum)

!     declarations for 'write_pin' call
      integer iwmostanum(recnum)
      character a9time_ob_r(recnum)*9
      logical l_closest_time, l_closest_time_i, l_in_domain, l_geoalt
      logical l_debug
      real*4 lat_a(nx_l,ny_l)
      real*4 lon_a(nx_l,ny_l)
      real*4 topo_a(nx_l,ny_l)

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'error getting r_missing_data'
          return
      endif
      call get_domain_perimeter(nx_l,ny_l,'nest7grid',lat_a,lon_a,
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error in get_domain_perimeter'
          return
      endif

      call read_wisdom_netcdf(nf_fid, qcchecknum, maxstaticids, 
     +     ninventorybins, recnum, cutdownflag, firstinbin, 
     +     firstoverflow, globalinventory, invtime, inventory, 
     +     isoverflow, lastinbin, lastrecord, latitudeqca, 
     +     latitudeqcr, longitudeqca, longitudeqcr, nstaticids, 
     +     pressureqca, pressureqcr, prevrecord, relhumidityqca, 
     +     relhumidityqcr, secondsstage1_2, temperatureqca, 
     +     temperatureqcr, winddirqca, winddirqcr, windspeedqca, 
     +     windspeedqcr, circularerror, elevation, latitude, 
     +     latitudeqcd, longitude, longitudeqcd, mobileelev, 
     +     mobilelat, mobilelon, pressure, pressureqcd, relhumidity, 
     +     relhumidityqcd, temperature, temperatureqcd, winddir, 
     +     winddirqcd, windspeed, windspeedqcd, qct, dataprovider, 
     +     handbook5id, homewfo, latitudedd, longitudedd, pressuredd, 
     +     providerid, rawmessage, relhumiditydd, staticids, 
     +     stationid, stationname, stationtype, temperaturedd, 
     +     winddirdd, windspeeddd, observationtime, receivedtime, 
     +     reporttime)
c
c the netcdf variables are filled - your pin write call may go here
c
      write(6,*)' # of wisdom reports read in = ',recnum

!     initial loop through obs to get times and stanums
      do iob = 1,recnum
          iwmostanum(iob) = 0
          if(abs(observationtime(iob)) .le. 1e10)then
              i4time_ob = idint(observationtime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

      icount_ob_written = 0

      do iob = 1,recnum
          l_geoalt = .true.

!         madis qc flag checks can be added here if desired

          if(iob .le. 10)then
              l_debug = .true.
          else
              l_debug = .false.
          endif

          call write_aircraft_sub(lun_out,'pin'
     1                           ,a9time_ob_r(iob),a9time_ob_r(iob)
     1                           ,i4time_sys
     1                           ,i4time_earliest          
     1                           ,i4time_latest            
     1                           ,latitude(iob),longitude(iob)
     1                           ,elevation(iob)
     1                           ,winddir(iob),windspeed(iob)
     1                           ,temperature(iob),relhumidity(iob)
     1                           ,l_geoalt                          ! i
     1                           ,l_debug                           ! i
     1                           ,istat_ob)                         ! o

          icount_ob_written = icount_ob_written + istat_ob

      enddo ! iob

      write(6,*)' # of wisdom reports written to pin file = '
     1          ,icount_ob_written        

      return
      end
c
c  subroutine to read the file "madis - wisdom" 
c
      subroutine read_wisdom_netcdf(nf_fid, qcchecknum, maxstaticids, 
     +     ninventorybins, recnum, cutdownflag, firstinbin, 
     +     firstoverflow, globalinventory, invtime, inventory, 
     +     isoverflow, lastinbin, lastrecord, latitudeqca, 
     +     latitudeqcr, longitudeqca, longitudeqcr, nstaticids, 
     +     pressureqca, pressureqcr, prevrecord, relhumidityqca, 
     +     relhumidityqcr, secondsstage1_2, temperatureqca, 
     +     temperatureqcr, winddirqca, winddirqcr, windspeedqca, 
     +     windspeedqcr, circularerror, elevation, latitude, 
     +     latitudeqcd, longitude, longitudeqcd, mobileelev, 
     +     mobilelat, mobilelon, pressure, pressureqcd, relhumidity, 
     +     relhumidityqcd, temperature, temperatureqcd, winddir, 
     +     winddirqcd, windspeed, windspeedqcd, qct, dataprovider, 
     +     handbook5id, homewfo, latitudedd, longitudedd, pressuredd, 
     +     providerid, rawmessage, relhumiditydd, staticids, 
     +     stationid, stationname, stationtype, temperaturedd, 
     +     winddirdd, windspeeddd, observationtime, receivedtime, 
     +     reporttime)
c
      include 'netcdf.inc'
      integer qcchecknum, maxstaticids, ninventorybins, recnum,nf_fid, 
     +     nf_vid, nf_status
      integer cutdownflag(recnum), firstinbin(ninventorybins),
     +     firstoverflow, globalinventory, invtime(recnum),
     +     inventory(maxstaticids), isoverflow(recnum),
     +     lastinbin(ninventorybins), lastrecord(maxstaticids),
     +     latitudeqca(recnum), latitudeqcr(recnum),
     +     longitudeqca(recnum), longitudeqcr(recnum), nstaticids,
     +     pressureqca(recnum), pressureqcr(recnum),
     +     prevrecord(recnum), relhumidityqca(recnum),
     +     relhumidityqcr(recnum), secondsstage1_2(recnum),
     +     temperatureqca(recnum), temperatureqcr(recnum),
     +     winddirqca(recnum), winddirqcr(recnum),
     +     windspeedqca(recnum), windspeedqcr(recnum)
      real circularerror(recnum), elevation(recnum), latitude(recnum),
     +     latitudeqcd( qcchecknum, recnum), longitude(recnum),
     +     longitudeqcd( qcchecknum, recnum), mobileelev(recnum),
     +     mobilelat(recnum), mobilelon(recnum), pressure(recnum),
     +     pressureqcd( qcchecknum, recnum), relhumidity(recnum),
     +     relhumidityqcd( qcchecknum, recnum), temperature(recnum),
     +     temperatureqcd( qcchecknum, recnum), winddir(recnum),
     +     winddirqcd( qcchecknum, recnum), windspeed(recnum),
     +     windspeedqcd( qcchecknum, recnum)
      double precision observationtime(recnum), receivedtime(recnum),
     +     reporttime(recnum)
      character pressuredd(recnum)
      character*11 dataprovider(recnum)
      character longitudedd(recnum)
      character winddirdd(recnum)
      character*60 qct(qcchecknum)
      character*6 handbook5id(recnum)
      character*51 stationname(recnum)
      character*11 stationtype(recnum)
      character latitudedd(recnum)
      character relhumiditydd(recnum)
      character*512 rawmessage(recnum)
      character*24 staticids(maxstaticids)
      character temperaturedd(recnum)
      character*12 providerid(recnum)
      character windspeeddd(recnum)
      character*4 homewfo(recnum)
      character*6 stationid(recnum)


c   variables of type real
c
c     variable        netcdf long name
c     circularerror "position fix error"
c
      nf_status=nf_inq_varid(nf_fid,'circularerror',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for circularerror'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,circularerror)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for circularerror'
       endif
      endif
c
c     variable        netcdf long name
c     elevation     "elevation"
c
      nf_status=nf_inq_varid(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for elevation'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,elevation)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for elevation'
       endif
      endif
c
c     variable        netcdf long name
c     latitude      "latitude"
c
      nf_status=nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for latitude'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,latitude)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for latitude'
       endif
      endif
c
c     variable        netcdf long name
c     latitudeqcd   "latitude qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'latitudeqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for latitudeqcd'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,latitudeqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for latitudeqcd'
       endif
      endif
c
c     variable        netcdf long name
c     longitude     "longitude"
c
      nf_status=nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for longitude'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,longitude)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for longitude'
       endif
      endif
c
c     variable        netcdf long name
c     longitudeqcd  "longitude qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'longitudeqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for longitudeqcd'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,longitudeqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for longitudeqcd'
       endif
      endif
c
c     variable        netcdf long name
c     mobileelev    "mobile elevation"
c
      nf_status=nf_inq_varid(nf_fid,'mobileelev',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for mobileelev'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,mobileelev)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for mobileelev'
       endif
      endif
c
c     variable        netcdf long name
c     mobilelat     "mobile latitude"
c
      nf_status=nf_inq_varid(nf_fid,'mobilelat',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for mobilelat'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,mobilelat)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for mobilelat'
       endif
      endif
c
c     variable        netcdf long name
c     mobilelon     "mobile longitude"
c
      nf_status=nf_inq_varid(nf_fid,'mobilelon',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for mobilelon'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,mobilelon)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for mobilelon'
       endif
      endif
c
c     variable        netcdf long name
c     pressure      "pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pressure'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pressure)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pressure'
       endif
      endif
c
c     variable        netcdf long name
c     pressureqcd   "pressure qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'pressureqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pressureqcd'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pressureqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pressureqcd'
       endif
      endif
c
c     variable        netcdf long name
c     relhumidity   "relhumidity"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidity',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for relhumidity'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,relhumidity)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for relhumidity'
       endif
      endif
c
c     variable        netcdf long name
c     relhumidityqcd"relhumidity qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidityqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for relhumidityqcd'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,relhumidityqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for relhumidityqcd'
       endif
      endif
c
c     variable        netcdf long name
c     temperature   "temperature"
c
      nf_status=nf_inq_varid(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for temperature'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,temperature)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for temperature'
       endif
      endif
c
c     variable        netcdf long name
c     temperatureqcd"temperature qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'temperatureqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for temperatureqcd'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,temperatureqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for temperatureqcd'
       endif
      endif
c
c     variable        netcdf long name
c     winddir       "wind direction"
c
      nf_status=nf_inq_varid(nf_fid,'winddir',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for winddir'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,winddir)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for winddir'
       endif
      endif
c
c     variable        netcdf long name
c     winddirqcd    "wind direction qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'winddirqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for winddirqcd'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,winddirqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for winddirqcd'
       endif
      endif
c
c     variable        netcdf long name
c     windspeed     "wind speed"
c
      nf_status=nf_inq_varid(nf_fid,'windspeed',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for windspeed'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,windspeed)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for windspeed'
       endif
      endif
c
c     variable        netcdf long name
c     windspeedqcd  "wind speed qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'windspeedqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for windspeedqcd'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,windspeedqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for windspeedqcd'
       endif
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c     cutdownflag   "cutdown flag"
c
      nf_status=nf_inq_varid(nf_fid,'cutdownflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for cutdownflag'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,cutdownflag)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for cutdownflag'
       endif
      endif
c
c     variable        netcdf long name
c     firstinbin    
c
      nf_status=nf_inq_varid(nf_fid,'firstinbin',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for firstinbin'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,firstinbin)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for firstinbin'
       endif
      endif
c
c     variable        netcdf long name
c     firstoverflow 
c
      nf_status=nf_inq_varid(nf_fid,'firstoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for firstoverflow'
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
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,globalinventory)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for globalinventory'
       endif
      endif
c
c     variable        netcdf long name
c     invtime       
c
      nf_status=nf_inq_varid(nf_fid,'invtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for invtime'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,invtime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for invtime'
       endif
      endif
c
c     variable        netcdf long name
c     inventory     
c
      nf_status=nf_inq_varid(nf_fid,'inventory',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for inventory'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,inventory)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for inventory'
       endif
      endif
c
c     variable        netcdf long name
c     isoverflow    
c
      nf_status=nf_inq_varid(nf_fid,'isoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for isoverflow'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,isoverflow)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for isoverflow'
       endif
      endif
c
c     variable        netcdf long name
c     lastinbin     
c
      nf_status=nf_inq_varid(nf_fid,'lastinbin',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for lastinbin'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,lastinbin)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for lastinbin'
       endif
      endif
c
c     variable        netcdf long name
c     lastrecord    
c
      nf_status=nf_inq_varid(nf_fid,'lastrecord',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for lastrecord'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,lastrecord)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for lastrecord'
       endif
      endif
c
c     variable        netcdf long name
c     latitudeqca   "latitude qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'latitudeqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for latitudeqca'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,latitudeqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for latitudeqca'
       endif
      endif
c
c     variable        netcdf long name
c     latitudeqcr   "latitude qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'latitudeqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for latitudeqcr'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,latitudeqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for latitudeqcr'
       endif
      endif
c
c     variable        netcdf long name
c     longitudeqca  "longitude qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'longitudeqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for longitudeqca'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,longitudeqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for longitudeqca'
       endif
      endif
c
c     variable        netcdf long name
c     longitudeqcr  "longitude qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'longitudeqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for longitudeqcr'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,longitudeqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for longitudeqcr'
       endif
      endif
c
c     variable        netcdf long name
c     nstaticids    
c
      nf_status=nf_inq_varid(nf_fid,'nstaticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for nstaticids'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,nstaticids)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for nstaticids'
       endif
      endif
c
c     variable        netcdf long name
c     pressureqca   "pressure qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'pressureqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pressureqca'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,pressureqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pressureqca'
       endif
      endif
c
c     variable        netcdf long name
c     pressureqcr   "pressure qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'pressureqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pressureqcr'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,pressureqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pressureqcr'
       endif
      endif
c
c     variable        netcdf long name
c     prevrecord    
c
      nf_status=nf_inq_varid(nf_fid,'prevrecord',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for prevrecord'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,prevrecord)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for prevrecord'
       endif
      endif
c
c     variable        netcdf long name
c     relhumidityqca"relhumidity qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidityqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for relhumidityqca'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,relhumidityqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for relhumidityqca'
       endif
      endif
c
c     variable        netcdf long name
c     relhumidityqcr"relhumidity qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidityqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for relhumidityqcr'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,relhumidityqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for relhumidityqcr'
       endif
      endif
c
c     variable        netcdf long name
c     secondsstage1_2
c
      nf_status=nf_inq_varid(nf_fid,'secondsstage1_2',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for secondsstage1_2'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,secondsstage1_2)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for secondsstage1_2'
       endif
      endif
c
c     variable        netcdf long name
c     temperatureqca"temperature qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'temperatureqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for temperatureqca'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,temperatureqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for temperatureqca'
       endif
      endif
c
c     variable        netcdf long name
c     temperatureqcr"temperature qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'temperatureqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for temperatureqcr'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,temperatureqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for temperatureqcr'
       endif
      endif
c
c     variable        netcdf long name
c     winddirqca    "wind direction qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'winddirqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for winddirqca'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,winddirqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for winddirqca'
       endif
      endif
c
c     variable        netcdf long name
c     winddirqcr    "wind direction qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'winddirqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for winddirqcr'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,winddirqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for winddirqcr'
       endif
      endif
c
c     variable        netcdf long name
c     windspeedqca  "wind speed qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'windspeedqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for windspeedqca'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,windspeedqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for windspeedqca'
       endif
      endif
c
c     variable        netcdf long name
c     windspeedqcr  "wind speed qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'windspeedqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for windspeedqcr'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,windspeedqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for windspeedqcr'
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
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,observationtime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for observationtime'
       endif
      endif
c
c     variable        netcdf long name
c     receivedtime  "time data was received"
c
      nf_status=nf_inq_varid(nf_fid,'receivedtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for receivedtime'
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,receivedtime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for receivedtime'
       endif
      endif
c
c     variable        netcdf long name
c     reporttime    "time data was processed by the provider"
c
      nf_status=nf_inq_varid(nf_fid,'reporttime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for reporttime'
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,reporttime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for reporttime'
       endif
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c     qct           "list of possible qc checks"
c
      nf_status=nf_inq_varid(nf_fid,'qct',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for qct'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,qct)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for qct'
       endif
      endif
c
c     variable        netcdf long name
c     dataprovider  "local data provider"
c
      nf_status=nf_inq_varid(nf_fid,'dataprovider',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dataprovider'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,dataprovider)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dataprovider'
       endif
      endif
c
c     variable        netcdf long name
c     handbook5id   "handbook5 id (afos or shef id)"
c
      nf_status=nf_inq_varid(nf_fid,'handbook5id',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for handbook5id'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,handbook5id)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for handbook5id'
       endif
      endif
c
c     variable        netcdf long name
c     homewfo       "home wfo id"
c
      nf_status=nf_inq_varid(nf_fid,'homewfo',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for homewfo'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,homewfo)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for homewfo'
       endif
      endif
c
c     variable        netcdf long name
c     latitudedd    "latitude qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'latitudedd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for latitudedd'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,latitudedd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for latitudedd'
       endif
      endif
c
c     variable        netcdf long name
c     longitudedd   "longitude qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'longitudedd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for longitudedd'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,longitudedd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for longitudedd'
       endif
      endif
c
c     variable        netcdf long name
c     pressuredd    "pressure qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'pressuredd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pressuredd'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,pressuredd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pressuredd'
       endif
      endif
c
c     variable        netcdf long name
c     providerid    "data provider station id"
c
      nf_status=nf_inq_varid(nf_fid,'providerid',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for providerid'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,providerid)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for providerid'
       endif
      endif
c
c     variable        netcdf long name
c     rawmessage    "raw text ldad mesonet message"
c
      nf_status=nf_inq_varid(nf_fid,'rawmessage',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for rawmessage'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,rawmessage)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for rawmessage'
       endif
      endif
c
c     variable        netcdf long name
c     relhumiditydd "relhumidity qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'relhumiditydd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for relhumiditydd'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,relhumiditydd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for relhumiditydd'
       endif
      endif
c
c     variable        netcdf long name
c     staticids     
c
      nf_status=nf_inq_varid(nf_fid,'staticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for staticids'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,staticids)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for staticids'
       endif
      endif
c
c     variable        netcdf long name
c     stationid     "alphanumeric station id"
c
      nf_status=nf_inq_varid(nf_fid,'stationid',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationid'
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
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,stationname)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationname'
       endif
      endif
c
c     variable        netcdf long name
c     stationtype   "ldad station type"
c
      nf_status=nf_inq_varid(nf_fid,'stationtype',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationtype'
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
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,temperaturedd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for temperaturedd'
       endif
      endif
c
c     variable        netcdf long name
c     winddirdd     "wind direction qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'winddirdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for winddirdd'
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
