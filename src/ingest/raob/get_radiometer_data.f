      subroutine get_radiometer_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer icchecknum, qcchecknum, level, maxstaticids,
     +     ninventorybins, recnum,nf_fid, nf_vid, nf_status
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
c get size of icchecknum
c
      nf_status=nf_inq_dimid(nf_fid,'icchecknum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim icchecknum'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,icchecknum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim icchecknum'
      endif
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
c get size of level
c
      nf_status=nf_inq_dimid(nf_fid,'level',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim level'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,level)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim level'
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
      call read_radiometer_data(nf_fid, icchecknum, qcchecknum, level,
     +     maxstaticids, ninventorybins, recnum, i4time_sys,
     +     ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, istatus)

      return
      end
c
c
      subroutine read_radiometer_data(nf_fid, icchecknum, qcchecknum,
     +     level, maxstaticids, ninventorybins, recnum, i4time_sys,
     +     ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, istatus)


      include 'netcdf.inc'
      integer icchecknum, qcchecknum, level, maxstaticids,
     +     ninventorybins, recnum,nf_fid, nf_vid, nf_status
      integer cloudbasetempica(recnum), cloudbasetempicr(recnum),
     +     cloudbasetempqca(recnum), cloudbasetempqcr(recnum),
     +     firstinbin(ninventorybins), firstoverflow,
     +     globalinventory, integratedliquidica(recnum),
     +     integratedliquidicr(recnum), integratedliquidqca(recnum),
     +     integratedliquidqcr(recnum), integratedvaporica(recnum),
     +     integratedvaporicr(recnum), integratedvaporqca(recnum),
     +     integratedvaporqcr(recnum), invtime(recnum),
     +     inventory(maxstaticids), isoverflow(recnum),
     +     lastinbin(ninventorybins), lastrecord(maxstaticids),
     +     liquiddensityica( level, recnum), liquiddensityicr( level,
     +     recnum), liquiddensityqca( level, recnum),
     +     liquiddensityqcr( level, recnum), nstaticids,
     +     prevrecord(recnum), rainflag(recnum), relhumidityica(
     +     level, recnum), relhumidityicr( level, recnum),
     +     relhumidityqca( level, recnum), relhumidityqcr( level,
     +     recnum), stationtype(recnum), temperatureica( level,
     +     recnum), temperatureicr( level, recnum), temperatureqca(
     +     level, recnum), temperatureqcr( level, recnum),
     +     vapordensityica( level, recnum), vapordensityicr( level,
     +     recnum), vapordensityqca( level, recnum), vapordensityqcr(
     +     level, recnum), wmostanum(recnum)
      real cloudbasetemp(recnum), elevation(recnum),
     +     integratedliquid(recnum), integratedvapor(recnum),
     +     latitude(recnum), levels( level, recnum), liquiddensity(
     +     level, recnum), longitude(recnum), pressure(recnum),
     +     relhumidity( level, recnum), relhumiditysfc(recnum),
     +     temperature( level, recnum), temperaturesfc(recnum),
     +     vapordensity( level, recnum)
      double precision observationtime(recnum)
      character*51 stationname(recnum)
      character*6 providerid(recnum)
      character vapordensitydd( level, recnum)
      character*30 staticids(maxstaticids)
      character temperaturedd( level, recnum)
      character relhumiditydd( level, recnum)
      character*11 dataprovider(recnum)
      character*72 ict(icchecknum)
      character cloudbasetempdd(recnum)
      character*60 qct(qcchecknum)
      character liquiddensitydd( level, recnum)
      character integratedliquiddd(recnum)
      character integratedvapordd(recnum)

!     declarations for 'write_snd' call
      integer iwmostanum(recnum)
      real stalat(level),stalon(level)
      character a9time_ob_r(recnum)*9,a9time_ob_l(level)*9
      character c8_obstype*8
      real height_m(level)
      real pressure_mb(level)
      real temp_c(level)
      real dewpoint_c(level)
      real dir_deg(level)
      real spd_mps(level)

      logical l_closest_time, l_closest_time_i
      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)

      integer max_lvls
      parameter (max_lvls=200)
      real liquid_a(max_lvls)

      common /write_snd_data/ cloud_base_temp,cloud_integrated_liquid
     1                       ,liquid_a

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

      call read_radiometer_netcdf(nf_fid, icchecknum, qcchecknum, 
     +     level, maxstaticids, ninventorybins, recnum, 
     +     cloudbasetempica, cloudbasetempicr, cloudbasetempqca, 
     +     cloudbasetempqcr, firstinbin, firstoverflow, 
     +     globalinventory, integratedliquidica, integratedliquidicr, 
     +     integratedliquidqca, integratedliquidqcr, 
     +     integratedvaporica, integratedvaporicr, 
     +     integratedvaporqca, integratedvaporqcr, invtime, 
     +     inventory, isoverflow, lastinbin, lastrecord, 
     +     liquiddensityica, liquiddensityicr, liquiddensityqca, 
     +     liquiddensityqcr, nstaticids, prevrecord, rainflag, 
     +     relhumidityica, relhumidityicr, relhumidityqca, 
     +     relhumidityqcr, stationtype, temperatureica, 
     +     temperatureicr, temperatureqca, temperatureqcr, 
     +     vapordensityica, vapordensityicr, vapordensityqca, 
     +     vapordensityqcr, wmostanum, cloudbasetemp, elevation, 
     +     integratedliquid, integratedvapor, latitude, levels, 
     +     liquiddensity, longitude, pressure, relhumidity, 
     +     relhumiditysfc, temperature, temperaturesfc, vapordensity, 
     +     observationtime, ict, qct, cloudbasetempdd, dataprovider, 
     +     integratedliquiddd, integratedvapordd, liquiddensitydd, 
     +     providerid, relhumiditydd, staticids, stationname, 
     +     temperaturedd, vapordensitydd)
c
c the netcdf variables are filled - your snd write call may go here
c
!     initial loop through obs to get times and stanums
      do iob = 1,recnum
          read(providerid(iob),'(3x,i2)',err=101)iwmostanum(iob)
          goto 102
101       iwmostanum(iob) = iob
          write(6,*)' warning: unreadable providerid '
     1             ,trim(providerid(iob))
102       continue

          if(abs(observationtime(iob)) .le. 1e10)then
              i4time_ob = idint(observationtime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

      c8_obstype = 'radiomtr'

      do iob = 1,recnum
          call convert_array(levels(:,iob),height_m,level
     1                      ,'none',r_missing_data,istatus)

          call addcon_miss(height_m,elevation(iob),height_m,level,1)

          stalat = latitude(iob)
          stalon = longitude(iob)

!         convert arrays for a single sounding
          a9time_ob_l = a9time_ob_r(iob)

          pressure_mb = r_missing_data

          call convert_array(pressure(iob),pressure_mb(1),1
     1                      ,'pa_to_mb',r_missing_data,istatus)

          call convert_array(temperature(:,iob),temp_c,level
     1                      ,'k_to_c',r_missing_data,istatus)

          do ilvl = 1,level
              dewpoint_c(ilvl) =
     1            dwpt_laps(temp_c(ilvl),relhumidity(ilvl,iob))
          enddo

          dir_deg = r_missing_data

          spd_mps = r_missing_data


          call get_nlevels_snd(pressure_mb,height_m,r_missing_data
     +                        ,level,nlevels_snd)

          if(integratedvapor(iob) .gt. .06)nlevels_snd=0

          cloud_base_temp = cloudbasetemp(iob)
          cloud_integrated_liquid = integratedliquid(iob)

          l_closest_time = l_closest_time_i(iwmostanum,a9time_ob_r
     1                                     ,recnum,iob,i4time_sys
     1                                     ,istatus)

          if(nlevels_snd .gt. 0 .and. l_closest_time)then
              call filter_string(providerid(iob))
              write(6,*)' valid radiometer near analysis time'
     1                 ,providerid(iob)

              write(6,*)' cloud_base_temp/liq = ',cloud_base_temp
     1                  ,cloud_integrated_liquid

              write(6,*)' liquid density:'
              do ilvl = 1,nlevels_snd
                  write(6,*)ilvl,height_m(ilvl),liquiddensity(ilvl,iob) 
                  liquid_a(ilvl) = liquiddensity(ilvl,iob)       
              enddo ! ilvl

!             call 'write_snd' for a single profile
              call open_ext(lun_out,i4time_sys,'snd',istatus)

              call write_snd(lun_out
     +                      ,1,nlevels_snd,1
     +                      ,iwmostanum
     +                      ,stalat,stalon,elevation(iob)
     +                      ,providerid(iob)
     +                      ,a9time_ob_l,c8_obstype
     +                      ,nlevels_snd
     +                      ,height_m
     +                      ,pressure_mb
     +                      ,temp_c
     +                      ,dewpoint_c
     +                      ,dir_deg
     +                      ,spd_mps
     +                      ,istatus)
          endif ! valid profile

      enddo ! iob
      return
      end
c
c  subroutine to read the file 
c
      subroutine read_radiometer_netcdf(nf_fid, icchecknum, 
     +     qcchecknum, level, maxstaticids, ninventorybins, recnum, 
     +     cloudbasetempica, cloudbasetempicr, cloudbasetempqca, 
     +     cloudbasetempqcr, firstinbin, firstoverflow, 
     +     globalinventory, integratedliquidica, integratedliquidicr, 
     +     integratedliquidqca, integratedliquidqcr, 
     +     integratedvaporica, integratedvaporicr, 
     +     integratedvaporqca, integratedvaporqcr, invtime, 
     +     inventory, isoverflow, lastinbin, lastrecord, 
     +     liquiddensityica, liquiddensityicr, liquiddensityqca, 
     +     liquiddensityqcr, nstaticids, prevrecord, rainflag, 
     +     relhumidityica, relhumidityicr, relhumidityqca, 
     +     relhumidityqcr, stationtype, temperatureica, 
     +     temperatureicr, temperatureqca, temperatureqcr, 
     +     vapordensityica, vapordensityicr, vapordensityqca, 
     +     vapordensityqcr, wmostanum, cloudbasetemp, elevation, 
     +     integratedliquid, integratedvapor, latitude, levels, 
     +     liquiddensity, longitude, pressure, relhumidity, 
     +     relhumiditysfc, temperature, temperaturesfc, vapordensity, 
     +     observationtime, ict, qct, cloudbasetempdd, dataprovider, 
     +     integratedliquiddd, integratedvapordd, liquiddensitydd, 
     +     providerid, relhumiditydd, staticids, stationname, 
     +     temperaturedd, vapordensitydd)
c
      include 'netcdf.inc'
      integer icchecknum, qcchecknum, level, maxstaticids, 
     +     ninventorybins, recnum,nf_fid, nf_vid, nf_status
      integer cloudbasetempica(recnum), cloudbasetempicr(recnum),
     +     cloudbasetempqca(recnum), cloudbasetempqcr(recnum),
     +     firstinbin(ninventorybins), firstoverflow,
     +     globalinventory, integratedliquidica(recnum),
     +     integratedliquidicr(recnum), integratedliquidqca(recnum),
     +     integratedliquidqcr(recnum), integratedvaporica(recnum),
     +     integratedvaporicr(recnum), integratedvaporqca(recnum),
     +     integratedvaporqcr(recnum), invtime(recnum),
     +     inventory(maxstaticids), isoverflow(recnum),
     +     lastinbin(ninventorybins), lastrecord(maxstaticids),
     +     liquiddensityica( level, recnum), liquiddensityicr( level,
     +     recnum), liquiddensityqca( level, recnum),
     +     liquiddensityqcr( level, recnum), nstaticids,
     +     prevrecord(recnum), rainflag(recnum), relhumidityica(
     +     level, recnum), relhumidityicr( level, recnum),
     +     relhumidityqca( level, recnum), relhumidityqcr( level,
     +     recnum), stationtype(recnum), temperatureica( level,
     +     recnum), temperatureicr( level, recnum), temperatureqca(
     +     level, recnum), temperatureqcr( level, recnum),
     +     vapordensityica( level, recnum), vapordensityicr( level,
     +     recnum), vapordensityqca( level, recnum), vapordensityqcr(
     +     level, recnum), wmostanum(recnum)
      real cloudbasetemp(recnum), elevation(recnum),
     +     integratedliquid(recnum), integratedvapor(recnum),
     +     latitude(recnum), levels( level, recnum), liquiddensity(
     +     level, recnum), longitude(recnum), pressure(recnum),
     +     relhumidity( level, recnum), relhumiditysfc(recnum),
     +     temperature( level, recnum), temperaturesfc(recnum),
     +     vapordensity( level, recnum)
      double precision observationtime(recnum)
      character*51 stationname(recnum)
      character*6 providerid(recnum)
      character vapordensitydd( level, recnum)
      character*30 staticids(maxstaticids)
      character temperaturedd( level, recnum)
      character relhumiditydd( level, recnum)
      character*11 dataprovider(recnum)
      character*72 ict(icchecknum)
      character cloudbasetempdd(recnum)
      character*60 qct(qcchecknum)
      character liquiddensitydd( level, recnum)
      character integratedliquiddd(recnum)
      character integratedvapordd(recnum)


c   variables of type real
c
c     variable        netcdf long name
c      cloudbasetemp"infrared cloud base temperature"
c
      nf_status=nf_inq_varid(nf_fid,'cloudbasetemp',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetemp'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,cloudbasetemp)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetemp'
      endif
c
c     variable        netcdf long name
c      elevation    "elevation above msl"
c
      nf_status=nf_inq_varid(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevation'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,elevation)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevation'
      endif
c
c     variable        netcdf long name
c      integratedliquid"integrated liquid"
c
      nf_status=nf_inq_varid(nf_fid,'integratedliquid',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquid'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,integratedliquid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquid'
      endif
c
c     variable        netcdf long name
c      integratedvapor"integrated vapor"
c
      nf_status=nf_inq_varid(nf_fid,'integratedvapor',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvapor'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,integratedvapor)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvapor'
      endif
c
c     variable        netcdf long name
c      latitude     "station latitude"
c
      nf_status=nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,latitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
c
c     variable        netcdf long name
c      levels       "instrument level, height above station"
c
      nf_status=nf_inq_varid(nf_fid,'levels',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var levels'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,levels)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var levels'
      endif
c
c     variable        netcdf long name
c      liquiddensity"liquid density"
c
      nf_status=nf_inq_varid(nf_fid,'liquiddensity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensity'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,liquiddensity)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensity'
      endif
c
c     variable        netcdf long name
c      longitude    "station longitude"
c
      nf_status=nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,longitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
c
c     variable        netcdf long name
c      pressure     "station pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var pressure'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,pressure)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var pressure'
      endif
c
c     variable        netcdf long name
c      relhumidity  "relative humidity"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidity'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,relhumidity)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidity'
      endif
c
c     variable        netcdf long name
c      relhumiditysfc"surface relative humidity"
c
      nf_status=nf_inq_varid(nf_fid,'relhumiditysfc',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumiditysfc'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,relhumiditysfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumiditysfc'
      endif
c
c     variable        netcdf long name
c      temperature  "temperature"
c
      nf_status=nf_inq_varid(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,temperature)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
c
c     variable        netcdf long name
c      temperaturesfc"surface temperature"
c
      nf_status=nf_inq_varid(nf_fid,'temperaturesfc',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperaturesfc'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,temperaturesfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperaturesfc'
      endif
c
c     variable        netcdf long name
c      vapordensity "vapor density"
c
      nf_status=nf_inq_varid(nf_fid,'vapordensity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensity'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,vapordensity)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensity'
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      cloudbasetempica"cloud base temp ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'cloudbasetempica',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetempica'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,cloudbasetempica)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetempica'
      endif
c
c     variable        netcdf long name
c      cloudbasetempicr"cloud base temp ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'cloudbasetempicr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetempicr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,cloudbasetempicr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetempicr'
      endif
c
c     variable        netcdf long name
c      cloudbasetempqca"cloud base temp qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'cloudbasetempqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetempqca'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,cloudbasetempqca)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetempqca'
      endif
c
c     variable        netcdf long name
c      cloudbasetempqcr"cloud base temp qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'cloudbasetempqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetempqcr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,cloudbasetempqcr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetempqcr'
      endif
c
c     variable        netcdf long name
c      firstinbin   
c
      nf_status=nf_inq_varid(nf_fid,'firstinbin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstinbin'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,firstinbin)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstinbin'
      endif
c
c     variable        netcdf long name
c      firstoverflow
c
      nf_status=nf_inq_varid(nf_fid,'firstoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstoverflow'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,firstoverflow)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstoverflow'
      endif
c
c     variable        netcdf long name
c      globalinventory
c
      nf_status=nf_inq_varid(nf_fid,'globalinventory',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var globalinventory'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,globalinventory)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var globalinventory'
      endif
c
c     variable        netcdf long name
c      integratedliquidica"integrated liquid ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'integratedliquidica',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquidica'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,integratedliquidica)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquidica'
      endif
c
c     variable        netcdf long name
c      integratedliquidicr"integrated liquid ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'integratedliquidicr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquidicr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,integratedliquidicr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquidicr'
      endif
c
c     variable        netcdf long name
c      integratedliquidqca"integrated liquid qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'integratedliquidqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquidqca'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,integratedliquidqca)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquidqca'
      endif
c
c     variable        netcdf long name
c      integratedliquidqcr"integrated liquid qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'integratedliquidqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquidqcr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,integratedliquidqcr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquidqcr'
      endif
c
c     variable        netcdf long name
c      integratedvaporica"integrated vapor ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'integratedvaporica',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvaporica'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,integratedvaporica)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvaporica'
      endif
c
c     variable        netcdf long name
c      integratedvaporicr"integrated vapor ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'integratedvaporicr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvaporicr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,integratedvaporicr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvaporicr'
      endif
c
c     variable        netcdf long name
c      integratedvaporqca"integrated vapor qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'integratedvaporqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvaporqca'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,integratedvaporqca)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvaporqca'
      endif
c
c     variable        netcdf long name
c      integratedvaporqcr"integrated vapor qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'integratedvaporqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvaporqcr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,integratedvaporqcr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvaporqcr'
      endif
c
c     variable        netcdf long name
c      invtime      
c
      nf_status=nf_inq_varid(nf_fid,'invtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var invtime'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,invtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var invtime'
      endif
c
c     variable        netcdf long name
c      inventory    
c
      nf_status=nf_inq_varid(nf_fid,'inventory',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var inventory'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,inventory)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var inventory'
      endif
c
c     variable        netcdf long name
c      isoverflow   
c
      nf_status=nf_inq_varid(nf_fid,'isoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var isoverflow'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,isoverflow)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var isoverflow'
      endif
c
c     variable        netcdf long name
c      lastinbin    
c
      nf_status=nf_inq_varid(nf_fid,'lastinbin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lastinbin'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,lastinbin)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lastinbin'
      endif
c
c     variable        netcdf long name
c      lastrecord   
c
      nf_status=nf_inq_varid(nf_fid,'lastrecord',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lastrecord'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,lastrecord)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lastrecord'
      endif
c
c     variable        netcdf long name
c      liquiddensityica"liquid density ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'liquiddensityica',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensityica'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,liquiddensityica)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensityica'
      endif
c
c     variable        netcdf long name
c      liquiddensityicr"liquid density ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'liquiddensityicr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensityicr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,liquiddensityicr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensityicr'
      endif
c
c     variable        netcdf long name
c      liquiddensityqca"liquid density qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'liquiddensityqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensityqca'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,liquiddensityqca)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensityqca'
      endif
c
c     variable        netcdf long name
c      liquiddensityqcr"liquid density qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'liquiddensityqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensityqcr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,liquiddensityqcr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensityqcr'
      endif
c
c     variable        netcdf long name
c      nstaticids   
c
      nf_status=nf_inq_varid(nf_fid,'nstaticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nstaticids'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,nstaticids)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nstaticids'
      endif
c
c     variable        netcdf long name
c      prevrecord   
c
      nf_status=nf_inq_varid(nf_fid,'prevrecord',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prevrecord'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,prevrecord)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prevrecord'
      endif
c
c     variable        netcdf long name
c      rainflag     "rain flag"
c
      nf_status=nf_inq_varid(nf_fid,'rainflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rainflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,rainflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rainflag'
      endif
c
c     variable        netcdf long name
c      relhumidityica"relative humidity ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidityica',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidityica'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,relhumidityica)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidityica'
      endif
c
c     variable        netcdf long name
c      relhumidityicr"relative humidity ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidityicr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidityicr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,relhumidityicr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidityicr'
      endif
c
c     variable        netcdf long name
c      relhumidityqca"relative humidity qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidityqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidityqca'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,relhumidityqca)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidityqca'
      endif
c
c     variable        netcdf long name
c      relhumidityqcr"relative humidity qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'relhumidityqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidityqcr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,relhumidityqcr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidityqcr'
      endif
c
c     variable        netcdf long name
c      stationtype  "station type"
c
      nf_status=nf_inq_varid(nf_fid,'stationtype',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationtype'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,stationtype)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationtype'
      endif
c
c     variable        netcdf long name
c      temperatureica"temperature ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'temperatureica',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperatureica'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,temperatureica)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperatureica'
      endif
c
c     variable        netcdf long name
c      temperatureicr"temperature ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'temperatureicr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperatureicr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,temperatureicr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperatureicr'
      endif
c
c     variable        netcdf long name
c      temperatureqca"temperature qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'temperatureqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperatureqca'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,temperatureqca)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperatureqca'
      endif
c
c     variable        netcdf long name
c      temperatureqcr"temperature qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'temperatureqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperatureqcr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,temperatureqcr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperatureqcr'
      endif
c
c     variable        netcdf long name
c      vapordensityica"vapor density ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'vapordensityica',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensityica'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,vapordensityica)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensityica'
      endif
c
c     variable        netcdf long name
c      vapordensityicr"vapor density ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'vapordensityicr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensityicr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,vapordensityicr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensityicr'
      endif
c
c     variable        netcdf long name
c      vapordensityqca"vapor density qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'vapordensityqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensityqca'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,vapordensityqca)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensityqca'
      endif
c
c     variable        netcdf long name
c      vapordensityqcr"vapor density qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'vapordensityqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensityqcr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,vapordensityqcr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensityqcr'
      endif
c
c     variable        netcdf long name
c      wmostanum    "wmo numeric station id"
c
      nf_status=nf_inq_varid(nf_fid,'wmostanum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wmostanum'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,wmostanum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wmostanum'
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c      observationtime"time of observation"
c
      nf_status=nf_inq_varid(nf_fid,'observationtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var observationtime'
      endif
      nf_status=nf_get_var_double(nf_fid,nf_vid,observationtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var observationtime'
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c      ict          "list of possible ic checks"
c
      nf_status=nf_inq_varid(nf_fid,'ict',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ict'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,ict)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ict'
      endif
c
c     variable        netcdf long name
c      qct          "list of possible qc checks"
c
      nf_status=nf_inq_varid(nf_fid,'qct',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var qct'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,qct)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var qct'
      endif
c
c     variable        netcdf long name
c      cloudbasetempdd"cloud base temp qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'cloudbasetempdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetempdd'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,cloudbasetempdd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbasetempdd'
      endif
c
c     variable        netcdf long name
c      dataprovider "name of organization responsible for delivering the data"
c
      nf_status=nf_inq_varid(nf_fid,'dataprovider',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dataprovider'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,dataprovider)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dataprovider'
      endif
c
c     variable        netcdf long name
c      integratedliquiddd"integrated liquid qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'integratedliquiddd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquiddd'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,integratedliquiddd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedliquiddd'
      endif
c
c     variable        netcdf long name
c      integratedvapordd"integrated vapor qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'integratedvapordd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvapordd'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,integratedvapordd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var integratedvapordd'
      endif
c
c     variable        netcdf long name
c      liquiddensitydd"liquid density qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'liquiddensitydd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensitydd'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,liquiddensitydd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var liquiddensitydd'
      endif
c
c     variable        netcdf long name
c      providerid   "alphanumeric station name"
c
      nf_status=nf_inq_varid(nf_fid,'providerid',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var providerid'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,providerid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var providerid'
      endif
c
c     variable        netcdf long name
c      relhumiditydd"relative humidity qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'relhumiditydd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumiditydd'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,relhumiditydd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumiditydd'
      endif
c
c     variable        netcdf long name
c      staticids    
c
      nf_status=nf_inq_varid(nf_fid,'staticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staticids'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,staticids)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staticids'
      endif
c
c     variable        netcdf long name
c      stationname  "alphanumeric station name"
c
      nf_status=nf_inq_varid(nf_fid,'stationname',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationname'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,stationname)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationname'
      endif
c
c     variable        netcdf long name
c      temperaturedd"temperature qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'temperaturedd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperaturedd'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,temperaturedd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperaturedd'
      endif
c
c     variable        netcdf long name
c      vapordensitydd"vapor density qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'vapordensitydd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensitydd'
      endif
      nf_status=nf_get_var_text(nf_fid,nf_vid,vapordensitydd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vapordensitydd'
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
