      subroutine get_ldad_madis_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*170 filename

      integer maxsensor, recnum,nf_fid, nf_vid, nf_status
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
c get size of maxsensor
c
      nf_status=nf_inq_dimid(nf_fid,'maxsensor',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxsensor'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,maxsensor)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxsensor'
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
      call read_ldad_madis_data(nf_fid, maxsensor, recnum, i4time_sys,
     +     ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, istatus)

      return
      end
c
c
      subroutine read_ldad_madis_data(nf_fid, maxsensor, recnum,
     +     i4time_sys, ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, istatus)


      include 'netcdf.inc'
      integer maxsensor, recnum,nf_fid, nf_vid, nf_status
      integer altimeterqcr(recnum), dewpointqcr(recnum),
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
     +     seasurfacetemp(recnum), soilmoisture(recnum),
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
      character altimeterdd(recnum)
      character*51 stationname(recnum)
      character precipratedd(recnum)
      character*11 stationtype(recnum)
      character stationpressuredd(recnum)
      character sealevelpressuredd(recnum)
      character windspeeddd(recnum)

!     declarations for 'write_lso' call
      integer iwmostanum(recnum)
      logical l_closest_time, l_closest_time_i
      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)

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

      call read_ldad_madis_netcdf(nf_fid, maxsensor, recnum, 
     +     altimeterqcr(ix), dewpointqcr(ix), firstoverflow, 
     +     globalinventory, nstaticids, numericwmoid, 
     +     precipaccumqcr(ix), precipintensity, 
     +     preciprateqcr(ix), preciptype, presschange3hourqcr(ix), 
     +     presschangechar, relhumidityqcr(ix), sealevelpressureqcr(ix),       
     +     stationpressureqcr(ix), temperatureqcr(ix), 
     +     visibilityqcr(ix),       
     +     winddirqcr(ix), windspeedqcr(ix), altimeter(ix), 
     +     dewpoint(ix), 
     +     elevation(ix), latitude(ix), longitude(ix), 
     +     meanweightedtemperature(ix), precipaccum(ix), 
     +     preciprate(ix), presschange3hour(ix), relhumidity(ix), 
     +     sealevelpressure(ix), seasurfacetemp(ix), 
     +     soilmoisture(ix), soiltemperature(ix), solarradiation(ix), 
     +     stationpressure(ix), temperature(ix), visibility(ix), 
     +     winddir(ix), winddirmax(ix), windgust(ix), windspeed(ix), 
     +     altimeterdd(ix), dataprovider(ix), dewpointdd(ix), 
     +     precipaccumdd(ix), precipratedd(ix), presweather(ix), 
     +     presschange3hourdd(ix), providerid(ix), relhumiditydd(ix), 
     +     sealevelpressuredd(ix), stationid(ix), stationname(ix), 
     +     stationpressuredd(ix), stationtype(ix), temperaturedd(ix), 
     +     visibilitydd(ix), winddirdd(ix), windspeeddd(ix), 
     +     observationtime(ix), receivedtime(ix), reporttime(ix), 
     +     rhchangetime(ix), stationpresschangetime(ix), 
     +     tempchangetime(ix), winddirchangetime(ix), 
     +     windgustchangetime(ix), windspeedchangetime(ix),badflag)
c
c the netcdf variables are filled - your lso write call may go here
c
!     initial loop through obs to get times and stanums
      do iob = 1,recnum
          iwmostanum(iob) = 0
          if(abs(observationtime(iob)) .le. 1e10)then
              i4time_ob = idint(observationtime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

!     c8_obstype = 

      do iob = 1,recnum
          l_closest_time = .true.

      enddo ! iob
      return
      end
