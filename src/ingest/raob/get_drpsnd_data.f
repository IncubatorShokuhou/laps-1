      subroutine get_drpsnd_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l
     +                   ,i4time_drpsnd_earliest,i4time_drpsnd_latest       
     +                   ,filename,isnd_staname
     +                   ,lun_out,l_fill_ht
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      logical l_fill_ht

      integer manlevel, recnum, sigtlevel, sigwlevel,
     +     troplevel,nf_fid, nf_vid, nf_status, nlvl_out
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
c get size of sigtlevel
c
      nf_status=nf_inq_dimid(nf_fid,'sigtlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigtlevel'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,sigtlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigtlevel'
      endif
c
c get size of sigwlevel
c
      nf_status=nf_inq_dimid(nf_fid,'sigwlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigwlevel'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,sigwlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigwlevel'
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
      nlvl_out = manlevel + sigtlevel + sigwlevel
      call read_drpsnd_data(nf_fid, manlevel, recnum, sigtlevel,
     +     sigwlevel, troplevel, i4time_sys, ilaps_cycle_time, nx_l,
     +     ny_l, i4time_drpsnd_earliest, i4time_drpsnd_latest, 
     +     nlvl_out, lun_out, l_fill_ht, isnd_staname, istatus)

      return
      end
c
c
      subroutine read_drpsnd_data(nf_fid, manlevel, recnum, sigtlevel,
     +     sigwlevel, troplevel, i4time_sys, ilaps_cycle_time, nx_l,
     +     ny_l, i4time_drpsnd_earliest, i4time_drpsnd_latest, nlvl_out,
     +     lun_out, l_fill_ht, isnd_staname, istatus)


      include 'netcdf.inc'
      integer manlevel, recnum, sigtlevel, sigwlevel,
     +     troplevel,nf_fid, nf_vid, nf_status
      integer marsdensquare(recnum), nummand(recnum), numsigt(recnum),
     +     numsigw(recnum), numtrop(recnum), prman( manlevel,
     +     recnum), prtrop( troplevel, recnum), wdman( manlevel,
     +     recnum), wdsigw( sigwlevel, recnum), wdtrop( troplevel,
     +     recnum)
      real htman( manlevel, recnum), latitude(recnum),
     +     longitude(recnum), prsigt( sigtlevel, recnum), prsigw(
     +     sigwlevel, recnum), tdman( manlevel, recnum), tdsigt(
     +     sigtlevel, recnum), tdtrop( troplevel, recnum), tpman(
     +     manlevel, recnum), tpsigt( sigtlevel, recnum), tptrop(
     +     troplevel, recnum), wsman( manlevel, recnum), wssigw(
     +     sigwlevel, recnum), wstrop( troplevel, recnum)
      double precision timenominal(recnum), timeobs(recnum)
      character*12 dropsondelocation(recnum)
      character*2048 rawdropsonde(recnum)

      integer wmostanum(recnum)
      real staelev(recnum)
      character*1 staname(6,recnum)
      character*6 c6_staname
      real wdman_r(manlevel,recnum)
      real prman_r(manlevel,recnum)
      real wdsigw_r(sigwlevel,recnum)
      real htsigw(sigwlevel, recnum)

      logical l_closest_time, l_closest_time_i, l_in_domain, l_fill_ht       
      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)
      character*8 c8_obstype
      character*9 a9time_sys,a9time_release,a9time_syn,a9time_drpsnd       

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

      call read_drpsnd_netcdf(nf_fid, manlevel, recnum, sigtlevel, 
     +     sigwlevel, troplevel, marsdensquare, nummand, numsigt, 
     +     numsigw, numtrop, prman, prtrop, wdman, wdsigw, wdtrop, 
     +     htman, latitude, longitude, prsigt, prsigw, tdman, tdsigt, 
     +     tdtrop, tpman, tpsigt, tptrop, wsman, wssigw, wstrop, 
     +     dropsondelocation, rawdropsonde, timenominal, timeobs)
c
c the netcdf variables are filled - your snd write call may go here
c
!     write all dropsondes to laps snd file

      r_nc_missing_data = 1e20

      htsigw = r_missing_data

      n_snd = recnum

      do isnd = 1,n_snd

!         qc and write out the sounding
          i4time_drpsnd = 0

          if(abs(timeobs(isnd)) .lt. 1e10 .and.
     1       abs(timeobs(isnd)) .gt.    0.      )then
              i4time_release = idint(timeobs(isnd))+315619200

              i4time_diff = i4time_release - i4time_sys

              if(abs(i4time_diff) .gt. 20000)then
                  write(6,*)' warning: i4time_release is not '
     1                     ,'consistent with i4time_diff'
     1                     ,i4time_release,i4time_sys
              endif

!             correction for balloon fall time to mid-sounding
              i4time_drpsnd = i4time_release + 100

          else
              i4time_release = 0
              i4time_diff = 0

          endif

          write(6,*)
          write(6,*)' drpsnd #',isnd,i4time_sys,i4time_release
     1                         ,i4time_diff,i4time_syn

          call make_fnam_lp(i4time_sys    , a9time_sys    , istatus)
          call make_fnam_lp(i4time_release, a9time_release, istatus)
          call make_fnam_lp(i4time_syn    , a9time_syn    , istatus)
          call make_fnam_lp(i4time_drpsnd , a9time_drpsnd , istatus)

          write(6,*)' times - sys/release/syn/drpsnd: '
     1             ,a9time_sys,' ',a9time_release,' '
     1             ,a9time_syn,' ',a9time_drpsnd

          if(latitude(isnd) .ge. r_nc_missing_data)then
              write(6,*)' missing first latitude',isnd
              goto 999
          endif

          if(longitude(isnd) .ge. r_nc_missing_data)then
              write(6,*)' missing first longitude',isnd
              goto 999
          endif

          if(latitude(isnd) .le. rnorth .and. latitude(isnd) .ge. south 
     1                                .and.      
     1       longitude(isnd) .ge. west   .and. longitude(isnd) .le. east      
     1                                                            )then       

!         if(.true.)then      ! for testing

              write(6,*)' drpsnd is inside domain lat/lon perimeter'
          else
              write(6,*)
     1            ' outside domain lat/lon perimeter - reject'
              goto 999
          endif

          if(i4time_drpsnd .ne. 0)then ! test window
              if(i4time_drpsnd .ge. i4time_drpsnd_earliest .and.
     1           i4time_drpsnd .le. i4time_drpsnd_latest)then
                  write(6,*)' inside time window'
              else
                  write(6,*)' outside time window - reject'
                  goto 999
              endif
          endif

          staelev(isnd) = -999.
          c8_obstype = 'dropsnd'

          isnd_staname = isnd_staname + 1
          write(c6_staname,1)isnd_staname
 1        format('d',i4.4,1x)
          do i = 1,6
              staname(i,isnd) = c6_staname(i:i)
          enddo 

!         convert arrays from integer to real
          do ilev = 1,nummand(isnd) 
              wdman_r(ilev,isnd) = wdman(ilev,isnd)
              prman_r(ilev,isnd) = prman(ilev,isnd)
          enddo ! l

          do ilev = 1,numsigw(isnd) 
              wdsigw_r(ilev,isnd) = wdsigw(ilev,isnd)
          enddo ! l

          call sort_and_write(i4time_sys,lun_out,l_fill_ht
     1                       ,recnum,isnd,r_missing_data,a9time_drpsnd
     1                       ,c8_obstype
     1                       ,wmostanum,staname,latitude,longitude
     1                       ,staelev
     1                       ,nummand,htman,prman_r,tpman,tdman      
     1                       ,wdman_r,wsman
     1                       ,numsigt,prsigt,tpsigt,tdsigt
     1                       ,numsigw,prsigw,htsigw,wdsigw_r,wssigw
     1                       ,nlvl_out 
     1                       ,manlevel,sigtlevel,sigwlevel,istatus)

          go to 999

 998      write(6,*)' error writing out drpsnd'

 999      continue

      enddo ! i

      return
      end
c
c  subroutine to read the file "drop sonde data : selected by ob time : time range from 1201190400 to 1201201200" 
c
      subroutine read_drpsnd_netcdf(nf_fid, manlevel, recnum, 
     +     sigtlevel, sigwlevel, troplevel, marsdensquare, nummand, 
     +     numsigt, numsigw, numtrop, prman, prtrop, wdman, wdsigw, 
     +     wdtrop, htman, latitude, longitude, prsigt, prsigw, tdman, 
     +     tdsigt, tdtrop, tpman, tpsigt, tptrop, wsman, wssigw, 
     +     wstrop, dropsondelocation, rawdropsonde, timenominal, 
     +     timeobs)
c
      include 'netcdf.inc'
      integer manlevel, recnum, sigtlevel, sigwlevel, 
     +     troplevel,nf_fid, nf_vid, nf_status
      integer marsdensquare(recnum), nummand(recnum), numsigt(recnum),
     +     numsigw(recnum), numtrop(recnum), prman( manlevel,
     +     recnum), prtrop( troplevel, recnum), wdman( manlevel,
     +     recnum), wdsigw( sigwlevel, recnum), wdtrop( troplevel,
     +     recnum)
      real htman( manlevel, recnum), latitude(recnum),
     +     longitude(recnum), prsigt( sigtlevel, recnum), prsigw(
     +     sigwlevel, recnum), tdman( manlevel, recnum), tdsigt(
     +     sigtlevel, recnum), tdtrop( troplevel, recnum), tpman(
     +     manlevel, recnum), tpsigt( sigtlevel, recnum), tptrop(
     +     troplevel, recnum), wsman( manlevel, recnum), wssigw(
     +     sigwlevel, recnum), wstrop( troplevel, recnum)
      double precision timenominal(recnum), timeobs(recnum)
      character*12 dropsondelocation(recnum)
      character*2048 rawdropsonde(recnum)


c   variables of type real
c
c     variable        netcdf long name
c     htman         "geopotential - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'htman',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for htman'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,htman)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for htman'
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
c     prsigt        "pressure - significant level wrt t"
c
      nf_status=nf_inq_varid(nf_fid,'prsigt',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for prsigt'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,prsigt)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for prsigt'
       endif
      endif
c
c     variable        netcdf long name
c     prsigw        "pressure - significant level wrt w"
c
      nf_status=nf_inq_varid(nf_fid,'prsigw',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for prsigw'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,prsigw)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for prsigw'
       endif
      endif
c
c     variable        netcdf long name
c     tdman         "dew point depression - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'tdman',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for tdman'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,tdman)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for tdman'
       endif
      endif
c
c     variable        netcdf long name
c     tdsigt        "dew point depression - significant level wrt t"
c
      nf_status=nf_inq_varid(nf_fid,'tdsigt',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for tdsigt'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,tdsigt)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for tdsigt'
       endif
      endif
c
c     variable        netcdf long name
c     tdtrop        "dew point depression - tropopause levels"
c
      nf_status=nf_inq_varid(nf_fid,'tdtrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for tdtrop'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,tdtrop)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for tdtrop'
       endif
      endif
c
c     variable        netcdf long name
c     tpman         "temperature - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'tpman',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for tpman'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,tpman)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for tpman'
       endif
      endif
c
c     variable        netcdf long name
c     tpsigt        "temperature - significant level wrt t"
c
      nf_status=nf_inq_varid(nf_fid,'tpsigt',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for tpsigt'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,tpsigt)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for tpsigt'
       endif
      endif
c
c     variable        netcdf long name
c     tptrop        "temperature - tropopause level"
c
      nf_status=nf_inq_varid(nf_fid,'tptrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for tptrop'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,tptrop)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for tptrop'
       endif
      endif
c
c     variable        netcdf long name
c     wsman         "wind speed"
c
      nf_status=nf_inq_varid(nf_fid,'wsman',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for wsman'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,wsman)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for wsman'
       endif
      endif
c
c     variable        netcdf long name
c     wssigw        "wind speed - significant level wrt w"
c
      nf_status=nf_inq_varid(nf_fid,'wssigw',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for wssigw'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,wssigw)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for wssigw'
       endif
      endif
c
c     variable        netcdf long name
c     wstrop        "wind speed - tropopause levels"
c
      nf_status=nf_inq_varid(nf_fid,'wstrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for wstrop'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,wstrop)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for wstrop'
       endif
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c     marsdensquare "marsden square"
c
      nf_status=nf_inq_varid(nf_fid,'marsdensquare',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for marsdensquare'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,marsdensquare)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for marsdensquare'
       endif
      endif
c
c     variable        netcdf long name
c     nummand       "number of mandatory levels"
c
      nf_status=nf_inq_varid(nf_fid,'nummand',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for nummand'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,nummand)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for nummand'
       endif
      endif
c
c     variable        netcdf long name
c     numsigt       "number of significant levels wrt t"
c
      nf_status=nf_inq_varid(nf_fid,'numsigt',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for numsigt'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,numsigt)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for numsigt'
       endif
      endif
c
c     variable        netcdf long name
c     numsigw       "number of significant levels wrt w"
c
      nf_status=nf_inq_varid(nf_fid,'numsigw',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for numsigw'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,numsigw)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for numsigw'
       endif
      endif
c
c     variable        netcdf long name
c     numtrop       "number of tropopause levels"
c
      nf_status=nf_inq_varid(nf_fid,'numtrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for numtrop'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,numtrop)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for numtrop'
       endif
      endif
c
c     variable        netcdf long name
c     prman         "pressure - mandatory level"
c
      nf_status=nf_inq_varid(nf_fid,'prman',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for prman'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,prman)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for prman'
       endif
      endif
c
c     variable        netcdf long name
c     prtrop        "pressure - tropopause level"
c
      nf_status=nf_inq_varid(nf_fid,'prtrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for prtrop'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,prtrop)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for prtrop'
       endif
      endif
c
c     variable        netcdf long name
c     wdman         "wind direction"
c
      nf_status=nf_inq_varid(nf_fid,'wdman',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for wdman'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,wdman)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for wdman'
       endif
      endif
c
c     variable        netcdf long name
c     wdsigw        "wind direction - significant level wrt w"
c
      nf_status=nf_inq_varid(nf_fid,'wdsigw',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for wdsigw'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,wdsigw)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for wdsigw'
       endif
      endif
c
c     variable        netcdf long name
c     wdtrop        "wind direction - tropopause levels"
c
      nf_status=nf_inq_varid(nf_fid,'wdtrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for wdtrop'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,wdtrop)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for wdtrop'
       endif
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c     timenominal   "drop sonde data hour"
c
      nf_status=nf_inq_varid(nf_fid,'timenominal',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for timenominal'
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,timenominal)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for timenominal'
       endif
      endif
c
c     variable        netcdf long name
c     timeobs       "observation time"
c
      nf_status=nf_inq_varid(nf_fid,'timeobs',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for timeobs'
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,timeobs)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for timeobs'
       endif
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c     dropsondelocation"dropsonde location or observation type (from 62626)"
c
      nf_status=nf_inq_varid(nf_fid,'dropsondelocation',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dropsondelocation'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,dropsondelocation)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dropsondelocation'
       endif
      endif
c
c     variable        netcdf long name
c     rawdropsonde  "raw drop sonde ascii message"
c
      nf_status=nf_inq_varid(nf_fid,'rawdropsonde',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for rawdropsonde'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,rawdropsonde)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for rawdropsonde'
       endif
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
