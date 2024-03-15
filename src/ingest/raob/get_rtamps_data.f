      subroutine get_rtamps_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*170 filename

      integer manlevel, maxstaticids, ninventorybins, rawlevel,
     +     recnum, stdlevel, termlevel, troplevel,nf_fid, nf_vid,
     +     nf_status
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
c get size of manlevel
c
      nf_status=nf_inq_dimid(nf_fid,'manlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim manlevel'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,manlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim manlevel'
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
c get size of rawlevel
c
      nf_status=nf_inq_dimid(nf_fid,'rawlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim rawlevel'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,rawlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim rawlevel'
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
c
c get size of stdlevel
c
      nf_status=nf_inq_dimid(nf_fid,'stdlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim stdlevel'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,stdlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim stdlevel'
      endif
c
c get size of termlevel
c
      nf_status=nf_inq_dimid(nf_fid,'termlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim termlevel'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,termlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim termlevel'
      endif
c
c get size of troplevel
c
      nf_status=nf_inq_dimid(nf_fid,'troplevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim troplevel'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,troplevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim troplevel'
      endif
      call read_rtamps_data(nf_fid, manlevel, maxstaticids,
     +     ninventorybins, rawlevel, recnum, stdlevel, termlevel,
     +     troplevel, i4time_sys, ilaps_cycle_time, nx_l, ny_l,
     +     i4time_earliest, i4time_latest, lun_out, istatus)

      return
      end
c
c
      subroutine read_rtamps_data(nf_fid, manlevel, maxstaticids,
     +     ninventorybins, rawlevel, recnum, stdlevel, termlevel,
     +     troplevel, i4time_sys, ilaps_cycle_time, nx_l, ny_l,
     +     i4time_earliest, i4time_latest, lun_out, istatus)


      include 'netcdf.inc'
      integer manlevel, maxstaticids, ninventorybins, rawlevel,
     +     recnum, stdlevel, termlevel, troplevel,nf_fid, nf_vid,
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
      character*11 dataprovider(recnum)
      character*51 stationname(recnum)
      character*24 staticids(maxstaticids)

!     declarations for 'write_snd' call
      integer iwmostanum(recnum)
      real stalat(rawlevel),stalon(rawlevel)
      character a9time_ob_r(recnum)*9,a9time_ob_l(rawlevel)*9
      character c8_obstype*8
      real height_m(rawlevel)
      real pressure_mb(rawlevel)
      real temp_c(rawlevel)
      real dewpoint_c(rawlevel)
      real dir_deg(rawlevel)
      real spd_mps(rawlevel)

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

      call read_rtamps_netcdf(nf_fid, manlevel, maxstaticids, 
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
c the netcdf variables are filled - your snd write call may go here
c
!     initial loop through obs to get times and stanums
      do iob = 1,recnum
          read(providerid(iob),*)iwmostanum(iob)
          if(abs(observationtime(iob)) .le. 1e10)then
              i4time_ob = idint(observationtime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

      c8_obstype = 'raob    '

      do iob = 1,recnum
          call convert_array(geopheight(:,iob),height_m,rawlevel
     1                      ,'none',r_missing_data,istatus)

          call addcon_miss(height_m,elevation(iob),height_m,rawlevel,1)

          stalat = latitude(iob)
          stalon = longitude(iob)

!         convert arrays for a single sounding
          a9time_ob_l = a9time_ob_r(iob)

          call convert_array(barompressure(:,iob),pressure_mb,rawlevel
     1                      ,'none',r_missing_data,istatus)

          call convert_array(temperature(:,iob),temp_c,rawlevel
     1                      ,'k_to_c',r_missing_data,istatus)

          call convert_array(dewpt(:,iob),dewpoint_c,rawlevel
     1                      ,'k_to_c',r_missing_data,istatus)

          call convert_array(direction(:,iob),dir_deg,rawlevel
     1                      ,'none',r_missing_data,istatus)

          call convert_array(speed(:,iob),spd_mps,rawlevel
     1                      ,'none',r_missing_data,istatus)


          call get_nlevels_snd(pressure_mb,height_m,r_missing_data
     +                        ,rawlevel,nlevels_snd)

!         apply qc editflag
          do i = 1,nlevels_snd
              if(editflag(i,iob) .eq. 3)then ! set wind to missing
                  direction(iob,3) = r_missing_data
                  speed(iob,3) = r_missing_data
              endif
          enddo ! i

          l_closest_time = .true.

          if(nlevels_snd .gt. 0 .and. l_closest_time)then
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
      subroutine read_rtamps_netcdf(nf_fid, manlevel, maxstaticids, 
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
     +     recnum, stdlevel, termlevel, troplevel,nf_fid, nf_vid, 
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
      character*11 dataprovider(recnum)
      character*51 stationname(recnum)
      character*24 staticids(maxstaticids)


c   variables of type real
c
c     variable        netcdf long name
c      abshumidity  "absolute humidity"
c
      nf_status=nf_inq_varid(nf_fid,'abshumidity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var abshumidity'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,abshumidity)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var abshumidity'
      endif
c
c     variable        netcdf long name
c      airdensity   "density of air"
c
      nf_status=nf_inq_varid(nf_fid,'airdensity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var airdensity'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,airdensity)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var airdensity'
      endif
c
c     variable        netcdf long name
c      barompressure"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'barompressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var barompressure'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,barompressure)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var barompressure'
      endif
c
c     variable        netcdf long name
c      bpman        "pressure - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'bpman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpman'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,bpman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpman'
      endif
c
c     variable        netcdf long name
c      bpstd        "pressure - standard level"
c
      nf_status=nf_inq_varid(nf_fid,'bpstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpstd'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,bpstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpstd'
      endif
c
c     variable        netcdf long name
c      bpterm       "pressure - termination"
c
      nf_status=nf_inq_varid(nf_fid,'bpterm',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpterm'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,bpterm)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bpterm'
      endif
c
c     variable        netcdf long name
c      bptrop       "pressure - tropopause level"
c
      nf_status=nf_inq_varid(nf_fid,'bptrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bptrop'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,bptrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bptrop'
      endif
c
c     variable        netcdf long name
c      dewpt        "dew point temperature"
c
      nf_status=nf_inq_varid(nf_fid,'dewpt',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dewpt'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,dewpt)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dewpt'
      endif
c
c     variable        netcdf long name
c      direction    "wind direction"
c
      nf_status=nf_inq_varid(nf_fid,'direction',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var direction'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,direction)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var direction'
      endif
c
c     variable        netcdf long name
c      dpman        "dew point temperature - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'dpman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dpman'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,dpman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dpman'
      endif
c
c     variable        netcdf long name
c      dpstd        "dew point temperature - standard level"
c
      nf_status=nf_inq_varid(nf_fid,'dpstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dpstd'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,dpstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dpstd'
      endif
c
c     variable        netcdf long name
c      dptrop       "dew point temperature - tropopause level"
c
      nf_status=nf_inq_varid(nf_fid,'dptrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dptrop'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,dptrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dptrop'
      endif
c
c     variable        netcdf long name
c      drman        "wind direction - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'drman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var drman'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,drman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var drman'
      endif
c
c     variable        netcdf long name
c      drstd        "wind direction - standard level"
c
      nf_status=nf_inq_varid(nf_fid,'drstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var drstd'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,drstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var drstd'
      endif
c
c     variable        netcdf long name
c      elevation    "station elevation"
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
c      geomheight   "geometric height"
c
      nf_status=nf_inq_varid(nf_fid,'geomheight',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var geomheight'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,geomheight)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var geomheight'
      endif
c
c     variable        netcdf long name
c      geopheight   "geopotential height"
c
      nf_status=nf_inq_varid(nf_fid,'geopheight',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var geopheight'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,geopheight)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var geopheight'
      endif
c
c     variable        netcdf long name
c      ghman        "geometric - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'ghman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghman'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,ghman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghman'
      endif
c
c     variable        netcdf long name
c      ghterm       "geometric - termination"
c
      nf_status=nf_inq_varid(nf_fid,'ghterm',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghterm'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,ghterm)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghterm'
      endif
c
c     variable        netcdf long name
c      ghtrop       "geometric - tropopause level"
c
      nf_status=nf_inq_varid(nf_fid,'ghtrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghtrop'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,ghtrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ghtrop'
      endif
c
c     variable        netcdf long name
c      gpstd        "geopotential - standard level"
c
      nf_status=nf_inq_varid(nf_fid,'gpstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gpstd'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,gpstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gpstd'
      endif
c
c     variable        netcdf long name
c      gpterm       "geopotential - termination"
c
      nf_status=nf_inq_varid(nf_fid,'gpterm',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gpterm'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,gpterm)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gpterm'
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
c      precipwater  "precipitable water"
c
      nf_status=nf_inq_varid(nf_fid,'precipwater',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var precipwater'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,precipwater)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var precipwater'
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
c      rhman        "relative humidity - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'rhman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rhman'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,rhman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rhman'
      endif
c
c     variable        netcdf long name
c      rhstd        "relative humidity - standard level"
c
      nf_status=nf_inq_varid(nf_fid,'rhstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rhstd'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,rhstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rhstd'
      endif
c
c     variable        netcdf long name
c      riserate     "rise rate"
c
      nf_status=nf_inq_varid(nf_fid,'riserate',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var riserate'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,riserate)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var riserate'
      endif
c
c     variable        netcdf long name
c      shear        "shear"
c
      nf_status=nf_inq_varid(nf_fid,'shear',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shear'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,shear)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shear'
      endif
c
c     variable        netcdf long name
c      sheardir     "shear direction"
c
      nf_status=nf_inq_varid(nf_fid,'sheardir',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sheardir'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,sheardir)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sheardir'
      endif
c
c     variable        netcdf long name
c      shearmagx    "shear magnitude x-direction"
c
      nf_status=nf_inq_varid(nf_fid,'shearmagx',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shearmagx'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,shearmagx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shearmagx'
      endif
c
c     variable        netcdf long name
c      shearmagy    "shear magnitude y-direction"
c
      nf_status=nf_inq_varid(nf_fid,'shearmagy',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shearmagy'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,shearmagy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var shearmagy'
      endif
c
c     variable        netcdf long name
c      spman        "wind speed - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'spman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var spman'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,spman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var spman'
      endif
c
c     variable        netcdf long name
c      spstd        "wind speed - standard level"
c
      nf_status=nf_inq_varid(nf_fid,'spstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var spstd'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,spstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var spstd'
      endif
c
c     variable        netcdf long name
c      speed        "wind speed"
c
      nf_status=nf_inq_varid(nf_fid,'speed',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var speed'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,speed)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var speed'
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
c      tpman        "temperature - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'tpman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpman'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,tpman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpman'
      endif
c
c     variable        netcdf long name
c      tpstd        "temperature - standard level"
c
      nf_status=nf_inq_varid(nf_fid,'tpstd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpstd'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,tpstd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpstd'
      endif
c
c     variable        netcdf long name
c      tptrop       "temperature - tropopause level"
c
      nf_status=nf_inq_varid(nf_fid,'tptrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tptrop'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,tptrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tptrop'
      endif
c
c     variable        netcdf long name
c      vaporpressure"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'vaporpressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vaporpressure'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,vaporpressure)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vaporpressure'
      endif
c
c     variable        netcdf long name
c      velerror     "velocity error"
c
      nf_status=nf_inq_varid(nf_fid,'velerror',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var velerror'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,velerror)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var velerror'
      endif
c
c     variable        netcdf long name
c      velsound     "velocity of sound"
c
      nf_status=nf_inq_varid(nf_fid,'velsound',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var velsound'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,velsound)
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
      nf_status=nf_inq_varid(nf_fid,'editflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var editflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,editflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var editflag'
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
c      indxrefr     "microwave index of refraction"
c
      nf_status=nf_inq_varid(nf_fid,'indxrefr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var indxrefr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,indxrefr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var indxrefr'
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
c      irman        "microwave index of refraction - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'irman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var irman'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,irman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var irman'
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
c      oiman        "optical index of refraction - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'oiman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var oiman'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,oiman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var oiman'
      endif
c
c     variable        netcdf long name
c      optindxrefr  "optical index of refraction"
c
      nf_status=nf_inq_varid(nf_fid,'optindxrefr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var optindxrefr'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,optindxrefr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var optindxrefr'
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
c      storedobs    "stored \'raw\' profile observations"
c
      nf_status=nf_inq_varid(nf_fid,'storedobs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var storedobs'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,storedobs)
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
c
c     variable        netcdf long name
c      receivedtime "received time"
c
      nf_status=nf_inq_varid(nf_fid,'receivedtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var receivedtime'
      endif
      nf_status=nf_get_var_double(nf_fid,nf_vid,receivedtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var receivedtime'
      endif
c
c     variable        netcdf long name
c      reporttime   "report time"
c
      nf_status=nf_inq_varid(nf_fid,'reporttime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reporttime'
      endif
      nf_status=nf_get_var_double(nf_fid,nf_vid,reporttime)
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
c      providerid   "data provider station id"
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

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
