      subroutine get_poes_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer maxlevels, maxstaticids, ninventorybins, recnum,nf_fid,
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
c get size of maxlevels
c
      nf_status=nf_inq_dimid(nf_fid,'maxlevels',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxlevels'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,maxlevels)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxlevels'
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

      write(6,*)' get_poes_data: number of records is ',recnum
      call read_poes_data(nf_fid, maxlevels, maxstaticids,
     +     ninventorybins, recnum, i4time_sys, ilaps_cycle_time,
     +     nx_l, ny_l, i4time_earliest, i4time_latest, lun_out,
     +     istatus)

      return
      end
c
c
      subroutine read_poes_data(nf_fid, maxlevels, maxstaticids,
     +     ninventorybins, recnum, i4time_sys, ilaps_cycle_time,
     +     nx_l, ny_l, i4time_earliest, i4time_latest, lun_out,
     +     istatus)


      include 'netcdf.inc'
      integer maxlevels, maxstaticids, ninventorybins, recnum,nf_fid,
     +     nf_vid, nf_status
      integer cloudamount(recnum), daynight(recnum),
     +     fieldofviewnum(recnum), firstinbin(ninventorybins),
     +     firstoverflow, globalinventory, invtime(recnum),
     +     inventory(maxstaticids), isoverflow(recnum),
     +     landsea(recnum), lastinbin(ninventorybins),
     +     lastrecord(maxstaticids), mixingratioqca( maxlevels,
     +     recnum), mixingratioqcr( maxlevels, recnum), nstaticids,
     +     numlevels(recnum), prevrecord(recnum), satproc(recnum),
     +     satelliteid(recnum), superadiabatic(recnum),
     +     temperatureqca( maxlevels, recnum), temperatureqcr(
     +     maxlevels, recnum), terrain(recnum)
      real cloudtoppressure(recnum), cloudtoptemperature(recnum),
     +     mixingratio( maxlevels, recnum), precipwater(recnum),
     +     pressure( maxlevels, recnum), skintemp(recnum),
     +     solarelev(recnum), staelev(recnum), stalat(recnum),
     +     stalon(recnum), stapress(recnum), temperature( maxlevels,
     +     recnum), zenithangle(recnum)
      double precision validtime(recnum)
      character temperaturedd( maxlevels, recnum)
      character*8 staticids(maxstaticids)
      character*7 staname(recnum)
      character mixingratiodd( maxlevels, recnum)

!     declarations for 'write_snd' call
      integer iwmostanum(recnum)
      real stalatl(maxlevels),stalonl(maxlevels)
      character a9time_ob_r(recnum)*9,a9time_ob_l(maxlevels)*9
      character staname_o(recnum)*5
      character c8_obstype*8
      real height_m(maxlevels)
      real pressure_mb(maxlevels)
      real temp_c(maxlevels)
      real dewpoint_c(maxlevels)
      real dir_deg(maxlevels)
      real spd_mps(maxlevels)

      integer iob_tot
      save iob_tot
      data iob_tot /0/

      logical l_closest_time, l_closest_time_i, l_in_domain
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

      call read_poes_netcdf(nf_fid, maxlevels, maxstaticids, 
     +     ninventorybins, recnum, cloudamount, daynight, 
     +     fieldofviewnum, firstinbin, firstoverflow, 
     +     globalinventory, invtime, inventory, isoverflow, landsea, 
     +     lastinbin, lastrecord, mixingratioqca, mixingratioqcr, 
     +     nstaticids, numlevels, prevrecord, satproc, satelliteid, 
     +     superadiabatic, temperatureqca, temperatureqcr, terrain, 
     +     cloudtoppressure, cloudtoptemperature, mixingratio, 
     +     precipwater, pressure, skintemp, solarelev, staelev, 
     +     stalat, stalon, stapress, temperature, zenithangle, 
     +     mixingratiodd, staname, staticids, temperaturedd, 
     +     validtime)
c
c the netcdf variables are filled - your snd write call may go here
c
!     initial loop through obs to get times and stanums
      do iob = 1,recnum
          iwmostanum(iob) = 0
          if(abs(validtime(iob)) .le. 1e10)then
              i4time_ob = idint(validtime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

!         create station name from ob count (try hex later if more needed)
          iob_tot = iob_tot + 1
          write(staname_o(iob),44) iob_tot
 44       format(i5.5)

      enddo ! iob

      c8_obstype = 'poessnd '

      do iob = 1,recnum
          height_m = r_missing_data
          stalatl = stalat(iob)
          stalonl = stalon(iob)

!         convert arrays for a single sounding
          a9time_ob_l = a9time_ob_r(iob)

          call convert_array(pressure(:,iob),pressure_mb,maxlevels
     1                      ,'pa_to_mb',r_missing_data,istatus)

          call convert_array(temperature(:,iob),temp_c,maxlevels
     1                      ,'k_to_c',r_missing_data,istatus)

          dewpoint_c = r_missing_data
          do ilvl = 1,maxlevels
              if(        mixingratio(ilvl,iob) .gt. 0.
     1             .and. mixingratio(ilvl,iob) .le. 1.)then
                  dewpoint_c(ilvl) =
     1              tmr(mixingratio(ilvl,iob)*1000.,pressure_mb(ilvl))
              endif
          enddo

          dir_deg = r_missing_data

          spd_mps = r_missing_data


          call get_nlevels_snd(pressure_mb,height_m,r_missing_data
     +                        ,maxlevels,nlevels_snd)

          l_closest_time = .true.

          if(   stalatl(1) .le. rnorth .and. stalatl(1) .ge. south
     1    .and. stalonl(1) .le. east   .and. stalonl(1) .ge. west
     1                                                       )then
              l_in_domain = .true.
          else
              l_in_domain = .false.
          endif

          if(nlevels_snd .gt. 0 .and. l_closest_time
     1                           .and. l_in_domain)then
!             call 'write_snd' for a single profile
              call open_ext(lun_out,i4time_sys,'snd',istatus)

              call write_snd(lun_out
     +                      ,1,nlevels_snd,1
     +                      ,iwmostanum
     +                      ,stalatl,stalonl,staelev(iob)
     +                      ,staname_o(iob)
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
      subroutine read_poes_netcdf(nf_fid, maxlevels, maxstaticids, 
     +     ninventorybins, recnum, cloudamount, daynight, 
     +     fieldofviewnum, firstinbin, firstoverflow, 
     +     globalinventory, invtime, inventory, isoverflow, landsea, 
     +     lastinbin, lastrecord, mixingratioqca, mixingratioqcr, 
     +     nstaticids, numlevels, prevrecord, satproc, satelliteid, 
     +     superadiabatic, temperatureqca, temperatureqcr, terrain, 
     +     cloudtoppressure, cloudtoptemperature, mixingratio, 
     +     precipwater, pressure, skintemp, solarelev, staelev, 
     +     stalat, stalon, stapress, temperature, zenithangle, 
     +     mixingratiodd, staname, staticids, temperaturedd, 
     +     validtime)
c
      include 'netcdf.inc'
      integer maxlevels, maxstaticids, ninventorybins, recnum,nf_fid, 
     +     nf_vid, nf_status
      integer cloudamount(recnum), daynight(recnum),
     +     fieldofviewnum(recnum), firstinbin(ninventorybins),
     +     firstoverflow, globalinventory, invtime(recnum),
     +     inventory(maxstaticids), isoverflow(recnum),
     +     landsea(recnum), lastinbin(ninventorybins),
     +     lastrecord(maxstaticids), mixingratioqca( maxlevels,
     +     recnum), mixingratioqcr( maxlevels, recnum), nstaticids,
     +     numlevels(recnum), prevrecord(recnum), satproc(recnum),
     +     satelliteid(recnum), superadiabatic(recnum),
     +     temperatureqca( maxlevels, recnum), temperatureqcr(
     +     maxlevels, recnum), terrain(recnum)
      real cloudtoppressure(recnum), cloudtoptemperature(recnum),
     +     mixingratio( maxlevels, recnum), precipwater(recnum),
     +     pressure( maxlevels, recnum), skintemp(recnum),
     +     solarelev(recnum), staelev(recnum), stalat(recnum),
     +     stalon(recnum), stapress(recnum), temperature( maxlevels,
     +     recnum), zenithangle(recnum)
      double precision validtime(recnum)
      character temperaturedd( maxlevels, recnum)
      character*8 staticids(maxstaticids)
      character*7 staname(recnum)
      character mixingratiodd( maxlevels, recnum)


c   variables of type real
c
c     variable        netcdf long name
c     cloudtoppressure"cloud top pressure"
c
      nf_status=nf_inq_varid(nf_fid,'cloudtoppressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for cloudtoppressure'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,cloudtoppressure)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for cloudtoppressure'
       endif
      endif
c
c     variable        netcdf long name
c     cloudtoptemperature"cloud top temperature"
c
      nf_status=nf_inq_varid(nf_fid,'cloudtoptemperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for cloudtoptemperature'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,cloudtoptemperature)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for cloudtoptemperature'
       endif
      endif
c
c     variable        netcdf long name
c     mixingratio   "mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'mixingratio',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for mixingratio'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,mixingratio)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for mixingratio'
       endif
      endif
c
c     variable        netcdf long name
c     precipwater   "total precipitable water"
c
      nf_status=nf_inq_varid(nf_fid,'precipwater',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipwater'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precipwater)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipwater'
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
c     skintemp      "skin temperature"
c
      nf_status=nf_inq_varid(nf_fid,'skintemp',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for skintemp'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,skintemp)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for skintemp'
       endif
      endif
c
c     variable        netcdf long name
c     solarelev     "grid point solar elevation"
c
      nf_status=nf_inq_varid(nf_fid,'solarelev',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for solarelev'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,solarelev)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for solarelev'
       endif
      endif
c
c     variable        netcdf long name
c     staelev       "grid point elevation"
c
      nf_status=nf_inq_varid(nf_fid,'staelev',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for staelev'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,staelev)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for staelev'
       endif
      endif
c
c     variable        netcdf long name
c     stalat        "grid point latitude"
c
      nf_status=nf_inq_varid(nf_fid,'stalat',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stalat'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,stalat)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stalat'
       endif
      endif
c
c     variable        netcdf long name
c     stalon        "grid point longitude"
c
      nf_status=nf_inq_varid(nf_fid,'stalon',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stalon'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,stalon)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stalon'
       endif
      endif
c
c     variable        netcdf long name
c     stapress      "grid point pressure"
c
      nf_status=nf_inq_varid(nf_fid,'stapress',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stapress'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,stapress)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stapress'
       endif
      endif
c
c     variable        netcdf long name
c     temperature   "retrieved temperature"
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
c     zenithangle   "satellite zenith angle"
c
      nf_status=nf_inq_varid(nf_fid,'zenithangle',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for zenithangle'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,zenithangle)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for zenithangle'
       endif
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c     cloudamount   "cloud amount"
c
      nf_status=nf_inq_varid(nf_fid,'cloudamount',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for cloudamount'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,cloudamount)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for cloudamount'
       endif
      endif
c
c     variable        netcdf long name
c     daynight      "day/night qualifier"
c
      nf_status=nf_inq_varid(nf_fid,'daynight',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for daynight'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,daynight)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for daynight'
       endif
      endif
c
c     variable        netcdf long name
c     fieldofviewnum"field of view number"
c
      nf_status=nf_inq_varid(nf_fid,'fieldofviewnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for fieldofviewnum'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,fieldofviewnum)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for fieldofviewnum'
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
c     landsea       "land/sea mask"
c
      nf_status=nf_inq_varid(nf_fid,'landsea',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for landsea'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,landsea)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for landsea'
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
c     mixingratioqca"mixingratio qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'mixingratioqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for mixingratioqca'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,mixingratioqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for mixingratioqca'
       endif
      endif
c
c     variable        netcdf long name
c     mixingratioqcr"mixingratio qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'mixingratioqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for mixingratioqcr'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,mixingratioqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for mixingratioqcr'
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
c     numlevels     "number of sounding levels"
c
      nf_status=nf_inq_varid(nf_fid,'numlevels',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for numlevels'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,numlevels)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for numlevels'
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
c     satproc       "satellite processing technique used"
c
      nf_status=nf_inq_varid(nf_fid,'satproc',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for satproc'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,satproc)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for satproc'
       endif
      endif
c
c     variable        netcdf long name
c     satelliteid   "satellite identifier"
c
      nf_status=nf_inq_varid(nf_fid,'satelliteid',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for satelliteid'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,satelliteid)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for satelliteid'
       endif
      endif
c
c     variable        netcdf long name
c     superadiabatic"superadiabatic indicator"
c
      nf_status=nf_inq_varid(nf_fid,'superadiabatic',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for superadiabatic'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,superadiabatic)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for superadiabatic'
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
c     terrain       "terrain type"
c
      nf_status=nf_inq_varid(nf_fid,'terrain',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for terrain'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,terrain)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for terrain'
       endif
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c     validtime     "sounding valid time"
c
      nf_status=nf_inq_varid(nf_fid,'validtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for validtime'
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,validtime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for validtime'
       endif
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c     mixingratiodd "mixingratio qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'mixingratiodd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for mixingratiodd'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,mixingratiodd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for mixingratiodd'
       endif
      endif
c
c     variable        netcdf long name
c     staname       "grid point station identifier"
c
      nf_status=nf_inq_varid(nf_fid,'staname',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for staname'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,staname)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for staname'
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

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
