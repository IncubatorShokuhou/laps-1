      subroutine get_sodar_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*170 filename

      integer level, maxstaticids, ninventorybins, recnum,nf_fid,
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
      call read_sodar_data(nf_fid, level, maxstaticids,
     +     ninventorybins, recnum, i4time_sys, ilaps_cycle_time,
     +     nx_l, ny_l, i4time_earliest, i4time_latest, lun_out,
     +     istatus)

      return
      end
c
c
      subroutine read_sodar_data(nf_fid, level, maxstaticids,
     +     ninventorybins, recnum, i4time_sys, ilaps_cycle_time,
     +     nx_l, ny_l, i4time_earliest, i4time_latest, lun_out,
     +     istatus)


      include 'netcdf.inc'
      integer level, maxstaticids, ninventorybins, recnum,nf_fid,
     +     nf_vid, nf_status
      integer assetid(recnum), firstinbin(ninventorybins),
     +     firstoverflow, globalinventory, invtime(recnum),
     +     inventory(maxstaticids), isoverflow(recnum),
     +     lastinbin(ninventorybins), lastrecord(maxstaticids),
     +     nstaticids, prevrecord(recnum), uqcflag( level, recnum),
     +     ustdqcflag( level, recnum), ugqcflag( level, recnum),
     +     vqcflag( level, recnum), vstdqcflag( level, recnum),
     +     vgqcflag( level, recnum), wqcflag( level, recnum),
     +     wstdqcflag( level, recnum), wdqcflag( level, recnum),
     +     winddir( level, recnum), wsqcflag( level, recnum)
      real elevation(recnum), latitude(recnum), levels( level,
     +     recnum), longitude(recnum), ucomponent( level, recnum),
     +     ugustcomponent( level, recnum), ustddevcomponent( level,
     +     recnum), vcomponent( level, recnum), vgustcomponent(
     +     level, recnum), vstddevcomponent( level, recnum),
     +     wcomponent( level, recnum), wstddevcomponent( level,
     +     recnum), windspeed( level, recnum)
      double precision observationtime(recnum), receipttime(recnum),
     +     reporttime(recnum)
      character*6 providerid(recnum)
      character*24 dataprovider(recnum)
      character*30 staticids(maxstaticids)

!     declarations for 'write_pro' call
      integer iwmostanum(recnum)
      character a9time_ob_r(recnum)*9
      character c8_obstype*8
      real height_m(level)
      real dir_deg(level)
      real spd_mps(level)

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

      call read_sodar_netcdf(nf_fid, level, maxstaticids, 
     +     ninventorybins, recnum, assetid, firstinbin, 
     +     firstoverflow, globalinventory, invtime, inventory, 
     +     isoverflow, lastinbin, lastrecord, nstaticids, prevrecord, 
     +     uqcflag, ustdqcflag, ugqcflag, vqcflag, vstdqcflag, 
     +     vgqcflag, wqcflag, wstdqcflag, wdqcflag, winddir, 
     +     wsqcflag, elevation, latitude, levels, longitude, 
     +     ucomponent, ugustcomponent, ustddevcomponent, vcomponent, 
     +     vgustcomponent, vstddevcomponent, wcomponent, 
     +     wstddevcomponent, windspeed, observationtime, receipttime, 
     +     reporttime, dataprovider, providerid, staticids)
c
c the netcdf variables are filled - your pro write call may go here
c
!     initial loop through obs to get times and stanums
      do iob = 1,recnum
          iwmostanum(iob) = assetid(iob)
          if(abs(observationtime(iob)) .le. 1e10)then
              i4time_ob = idint(observationtime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

      c8_obstype = 'sodar   '

      do iob = 1,recnum
          call convert_array(levels(:,iob),height_m,level
     1                      ,'none',r_missing_data,istatus)

          call addcon_miss(height_m,elevation(iob),height_m,level,1)

          call convert_array_i2r(winddir(:,iob),dir_deg,level
     1                      ,'none',r_missing_data,istatus)
          call apply_qc_rsa(wdqcflag(:,iob),dir_deg,level)
          call convert_array(dir_deg,dir_deg,level
     1                      ,'none',r_missing_data,istatus)

          call convert_array(windspeed(:,iob),spd_mps,level
     1                      ,'none',r_missing_data,istatus)
          call apply_qc_rsa(wsqcflag(:,iob),spd_mps,level)
          call convert_array(spd_mps,spd_mps,level
     1                      ,'none',r_missing_data,istatus)

          l_closest_time = l_closest_time_i(iwmostanum,a9time_ob_r
     1                                     ,recnum,iob,i4time_sys
     1                                     ,istatus)

          if(l_closest_time)then
!             call 'write_pro' for a single profile
              call open_ext(lun_out,i4time_sys,'pro',istatus)

              call write_pro(lun_out
     +                      ,1,level,1
     +                      ,assetid(iob)
     +                      ,latitude(iob),longitude(iob),elevation(iob)
     +                      ,providerid(iob)
     +                      ,a9time_ob_r(iob),c8_obstype
     +                      ,level
     +                      ,height_m
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
      subroutine read_sodar_netcdf(nf_fid, level, maxstaticids, 
     +     ninventorybins, recnum, assetid, firstinbin, 
     +     firstoverflow, globalinventory, invtime, inventory, 
     +     isoverflow, lastinbin, lastrecord, nstaticids, prevrecord, 
     +     uqcflag, ustdqcflag, ugqcflag, vqcflag, vstdqcflag, 
     +     vgqcflag, wqcflag, wstdqcflag, wdqcflag, winddir, 
     +     wsqcflag, elevation, latitude, levels, longitude, 
     +     ucomponent, ugustcomponent, ustddevcomponent, vcomponent, 
     +     vgustcomponent, vstddevcomponent, wcomponent, 
     +     wstddevcomponent, windspeed, observationtime, receipttime, 
     +     reporttime, dataprovider, providerid, staticids)
c
      include 'netcdf.inc'
      integer level, maxstaticids, ninventorybins, recnum,nf_fid, 
     +     nf_vid, nf_status
      integer assetid(recnum), firstinbin(ninventorybins),
     +     firstoverflow, globalinventory, invtime(recnum),
     +     inventory(maxstaticids), isoverflow(recnum),
     +     lastinbin(ninventorybins), lastrecord(maxstaticids),
     +     nstaticids, prevrecord(recnum), uqcflag( level, recnum),
     +     ustdqcflag( level, recnum), ugqcflag( level, recnum),
     +     vqcflag( level, recnum), vstdqcflag( level, recnum),
     +     vgqcflag( level, recnum), wqcflag( level, recnum),
     +     wstdqcflag( level, recnum), wdqcflag( level, recnum),
     +     winddir( level, recnum), wsqcflag( level, recnum)
      real elevation(recnum), latitude(recnum), levels( level,
     +     recnum), longitude(recnum), ucomponent( level, recnum),
     +     ugustcomponent( level, recnum), ustddevcomponent( level,
     +     recnum), vcomponent( level, recnum), vgustcomponent(
     +     level, recnum), vstddevcomponent( level, recnum),
     +     wcomponent( level, recnum), wstddevcomponent( level,
     +     recnum), windspeed( level, recnum)
      double precision observationtime(recnum), receipttime(recnum),
     +     reporttime(recnum)
      character*6 providerid(recnum)
      character*24 dataprovider(recnum)
      character*30 staticids(maxstaticids)


c   variables of type real
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
c      ucomponent   "u (eastward) component"
c
      nf_status=nf_inq_varid(nf_fid,'ucomponent',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ucomponent'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,ucomponent)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ucomponent'
      endif
c
c     variable        netcdf long name
c      ugustcomponent"gust u (eastward) component"
c
      nf_status=nf_inq_varid(nf_fid,'ugustcomponent',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ugustcomponent'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,ugustcomponent)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ugustcomponent'
      endif
c
c     variable        netcdf long name
c      ustddevcomponent"std dev u (eastward) component"
c
      nf_status=nf_inq_varid(nf_fid,'ustddevcomponent',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ustddevcomponent'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,ustddevcomponent)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ustddevcomponent'
      endif
c
c     variable        netcdf long name
c      vcomponent   "v (northward) component"
c
      nf_status=nf_inq_varid(nf_fid,'vcomponent',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vcomponent'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,vcomponent)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vcomponent'
      endif
c
c     variable        netcdf long name
c      vgustcomponent"gust v (northward) component"
c
      nf_status=nf_inq_varid(nf_fid,'vgustcomponent',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vgustcomponent'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,vgustcomponent)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vgustcomponent'
      endif
c
c     variable        netcdf long name
c      vstddevcomponent"std dev v (northward) component"
c
      nf_status=nf_inq_varid(nf_fid,'vstddevcomponent',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vstddevcomponent'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,vstddevcomponent)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vstddevcomponent'
      endif
c
c     variable        netcdf long name
c      wcomponent   "w (upward) component"
c
      nf_status=nf_inq_varid(nf_fid,'wcomponent',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wcomponent'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,wcomponent)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wcomponent'
      endif
c
c     variable        netcdf long name
c      wstddevcomponent"std dev w (upward) component"
c
      nf_status=nf_inq_varid(nf_fid,'wstddevcomponent',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wstddevcomponent'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,wstddevcomponent)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wstddevcomponent'
      endif
c
c     variable        netcdf long name
c      windspeed    "wind speed (scalar)"
c
      nf_status=nf_inq_varid(nf_fid,'windspeed',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeed'
      endif
      nf_status=nf_get_var_real(nf_fid,nf_vid,windspeed)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeed'
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      assetid      "rsa asset identifier"
c
      nf_status=nf_inq_varid(nf_fid,'assetid',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var assetid'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,assetid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var assetid'
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
c      uqcflag      "rsa mini-sodar wind u-component quality control flag"
c
      nf_status=nf_inq_varid(nf_fid,'uqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var uqcflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,uqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var uqcflag'
      endif
c
c     variable        netcdf long name
c      ustdqcflag   "rsa mini-sodar wind standard deviation u-component quality control flag"
c
      nf_status=nf_inq_varid(nf_fid,'ustdqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ustdqcflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,ustdqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ustdqcflag'
      endif
c
c     variable        netcdf long name
c      ugqcflag     "rsa mini-sodar wind gust u-component quality control flag"
c
      nf_status=nf_inq_varid(nf_fid,'ugqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ugqcflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,ugqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ugqcflag'
      endif
c
c     variable        netcdf long name
c      vqcflag      "rsa mini-sodar wind v-component quality control flag"
c
      nf_status=nf_inq_varid(nf_fid,'vqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vqcflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,vqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vqcflag'
      endif
c
c     variable        netcdf long name
c      vstdqcflag   "rsa mini-sodar wind standard deviation v-component quality control flag"
c
      nf_status=nf_inq_varid(nf_fid,'vstdqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vstdqcflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,vstdqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vstdqcflag'
      endif
c
c     variable        netcdf long name
c      vgqcflag     "rsa mini-sodar wind gust v-component quality control flag"
c
      nf_status=nf_inq_varid(nf_fid,'vgqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vgqcflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,vgqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vgqcflag'
      endif
c
c     variable        netcdf long name
c      wqcflag      "rsa mini-sodar wind w-component quality control flag"
c
      nf_status=nf_inq_varid(nf_fid,'wqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wqcflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,wqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wqcflag'
      endif
c
c     variable        netcdf long name
c      wstdqcflag   "rsa mini-sodar wind standard deviation w-component quality control flag"
c
      nf_status=nf_inq_varid(nf_fid,'wstdqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wstdqcflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,wstdqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wstdqcflag'
      endif
c
c     variable        netcdf long name
c      wdqcflag     "rsa wind direction quality control flag"
c
      nf_status=nf_inq_varid(nf_fid,'wdqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdqcflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,wdqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdqcflag'
      endif
c
c     variable        netcdf long name
c      winddir      "wind direction (scalar)"
c
      nf_status=nf_inq_varid(nf_fid,'winddir',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,winddir)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      endif
c
c     variable        netcdf long name
c      wsqcflag     "rsa wind speed quality control flag"
c
      nf_status=nf_inq_varid(nf_fid,'wsqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wsqcflag'
      endif
      nf_status=nf_get_var_int(nf_fid,nf_vid,wsqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wsqcflag'
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
c
c     variable        netcdf long name
c      receipttime  "file time stamp (time file was received)"
c
      nf_status=nf_inq_varid(nf_fid,'receipttime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var receipttime'
      endif
      nf_status=nf_get_var_double(nf_fid,nf_vid,receipttime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var receipttime'
      endif
c
c     variable        netcdf long name
c      reporttime   "time of observation"
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

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
