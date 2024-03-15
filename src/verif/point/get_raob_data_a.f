      subroutine get_raob_data_a(nx,ny,i4time_sys,i4time_raob_earliest,
     1                           i4time_raob_latest,filename,maxraob,
     1                           maxm,maxt,maxw,lat,lon,timesyn,timerel,
     1                           numsigt, numsigw, wmostanum,staname,
     1                           typew,typet,prsigt, tsigt,tdsigt,
     1                           htsigt, htsigw, wdsigw,wssigw, 
     1                           stalat, stalon, staelev,
     1                           max_ht_m_proc,min_pres_mb_proc,numraob, 
     1                           n_raobs_avail,verif_missing_data, 
     1                           raob_missing_data,istatus)

      implicit none
      include 'netcdf.inc'

!     input variables
      integer     nx,ny,i4time_sys,i4time_raob_earliest,
     1              i4time_raob_latest,maxraob,maxm,maxt,maxw
      character*(*) filename
      real        lat(nx,ny),lon(nx,ny)

!     variables filled for output
      integer     timesyn(maxraob),timerel(maxraob),
     1              numsigt(maxraob),numsigw(maxraob),
     1              istatus
      integer       numraob, n_raobs_avail
      real	    verif_missing_data,raob_missing_data
      character*6   staname(maxraob)
      integer	    wmostanum(maxraob),wmostanum_f
      character*1   typew(maxw,maxraob),typet(maxt,maxraob)
      real        prsigt(maxt,maxraob),tsigt(maxt,maxraob),
     1              tdsigt(maxt,maxraob),
     1              htsigt(maxt,maxraob),
     1              htsigw(maxw,maxraob),wdsigw(maxw,maxraob),
     1              wssigw(maxw,maxraob),staelev(maxraob), 
     1              max_ht_m_proc, min_pres_mb_proc,
     1              stalat(maxraob), stalon(maxraob) 

!     other variables
      real        htman(maxm,maxraob), prman(maxm,maxraob),
     1              tman(maxm,maxraob), tdman(maxm,maxraob),
     1              wsman(maxm,maxraob), wdman(maxm,maxraob),
     1              prsigti(maxt,maxraob),tsigti(maxt,maxraob),
     1              tdsigti(maxt,maxraob),
     1              htsigwi(maxw,maxraob),wdsigwi(maxw,maxraob),
     1              wssigwi(maxw,maxraob)

      double precision d_timerel, d_timesyn
      character*9   a9time_syn, a9time_release, a9time_raob, 
     1              a9time_sys
      character*8   c8_project
      character*6   staname_f
      integer     i, j, recnum, sigtlevel, sigwlevel, stanamelen, 
     1              stanamelen_f, nf_fid, nf_vid, nf_status,
     1              nobs, i4time_syn, i4time_release, i4time_diff,
     1              i4time_raob, manlevel, numman(maxraob),
     1              itime_delay,status
      real        stalat_f, stalon_f, staelev_f, ri, rj,
     1              r_nc_missing_data
      integer       index_1(1), start(2), count(2)

!.............................................................................
! begin

      n_raobs_avail = 0

      istatus = 1    ! assume a good return
      stanamelen = len(staname(1))

c
c  open netcdf file for reading
c
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        write(6,*) nf_strerror(nf_status)
        write(6,*)'nf_open ',filename
        istatus = 0
        return
      else
        write(6,*) 'opening file ',filename
      endif
c
c  read dimension values recnum, manlevel, sigtlevel, sigwlevel, stanamelen
c
c
c get size of recnum
c
      nf_status = nf_inq_dimid(nf_fid,'recnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
        istatus = 0
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,recnum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
        istatus = 0
        return
      endif
c
c get size of manlevel
c
      nf_status = nf_inq_dimid(nf_fid,'manlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim manlevel'
        istatus = 0
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,manlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim manlevel'
        istatus = 0
        return
      endif
c
c get size of sigtlevel
c
      nf_status = nf_inq_dimid(nf_fid,'sigtlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigtlevel'
        istatus = 0
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,sigtlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigtlevel'
        istatus = 0
        return
      endif
c
c get size of sigwlevel
c
      nf_status = nf_inq_dimid(nf_fid,'sigwlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigwlevel'
        istatus = 0
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,sigwlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigwlevel'
        istatus = 0
        return
      endif
c
c get size of stanamelen
c
      nf_status = nf_inq_dimid(nf_fid,'stanamelen',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim stanamelen'
        istatus = 0
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,stanamelen_f)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim stanamelen'
        istatus = 0
        return
      endif

      if (stanamelen_f .gt. stanamelen) then
        print *, 'staname truncated to ',stanamelen,' characters.'
        stanamelen_f = stanamelen
      endif
c
c     read missing data value for stalat and stalon
c
        nf_status = nf_inq_varid(nf_fid,'stalat',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var stalat'
          istatus = 0
          return
        endif
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',
     1                              r_nc_missing_data)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in att for stalat'
          istatus = 0
          return
        endif

! read data from netcdf file
      nobs = 0   !number of raobs stored to returning arrays
      start(1) = 1
      count(2) = 1

      do i = 1, recnum
        index_1(1) = i
        start(2) = i
c
c       read staname      "station identifier"
c
        nf_status = nf_inq_varid(nf_fid,'staname',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var staname'
        endif

        count(1) = stanamelen_f
        nf_status = nf_get_vara_text(nf_fid,nf_vid,start,
     1                               count,staname_f)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var staname'
        endif
c
c       read wmostanum    "wmo station number"
c
        nf_status = nf_inq_varid(nf_fid,'wmostanum',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var wmostanum'
        endif
        nf_status = nf_get_var1_int(nf_fid,nf_vid,index_1,wmostanum_f)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var wmostanum'
        endif
c
c       read stalat       "station latitude"
c
        nf_status = nf_inq_varid(nf_fid,'stalat',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var stalat'
        endif
        nf_status = nf_get_var1_real(nf_fid,nf_vid,index_1,stalat_f)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var stalat'
        endif
c
c       read stalon       "station longitude"
c
        nf_status = nf_inq_varid(nf_fid,'stalon',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var stalon'
        endif
        nf_status = nf_get_var1_real(nf_fid,nf_vid,index_1,stalon_f)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var stalon'
        endif

        if(stalat_f .ge. r_nc_missing_data)then
          goto 888
        endif

        if(stalon_f .ge. r_nc_missing_data)then
          goto 888
        endif

        call latlon_to_rlapsgrid(stalat_f,stalon_f,lat,lon,nx,ny,
     1                           ri,rj,status)
        if (status .ne. 1) then   !raob is not in laps domain
          goto 888
        else
          write(6,*)
          write(6,*) 'raob ',wmostanum_f,' is in laps domain',
     1               stalat_f,' ',stalon_f
          n_raobs_avail = n_raobs_avail + 1
        endif

!       read syntime and reltime and see if raob is in time window
c
c       read reltime      "sounding release time"

        nf_status = nf_inq_varid(nf_fid,'reltime',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var reltime'
        endif
        nf_status = nf_get_var1_double(nf_fid,nf_vid,index_1,d_timerel)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var reltime'
        endif
c
c       read syntime      "synoptic time"
c
        nf_status = nf_inq_varid(nf_fid,'syntime',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var syntime'
        endif
        nf_status = nf_get_var1_double(nf_fid,nf_vid,index_1,d_timesyn)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var syntime'
        endif

        i4time_raob = 0

        if(abs(d_timesyn) .lt. 1e10)then
          i4time_syn  = idint(d_timesyn)+315619200
          i4time_raob = i4time_syn
        else
          i4time_syn = 0
        endif
c
c somewhat bogus but allows synch of delays with time of data
c
        call get_c8_project(c8_project,istatus)
        call upcase(c8_project,c8_project)
        if(c8_project.eq.'airdrop')itime_delay=0  !2*3600

        if(abs(d_timerel) .lt. 1e10)then
          i4time_release = idint(d_timerel)+315619200

          i4time_diff = i4time_release - i4time_sys

          if(abs(i4time_diff) .gt. 20000)then
            write(6,*)' warning: i4time_release is not '
     1               ,'consistent with i4time_diff'
     1               ,i4time_release,i4time_sys
            call make_fnam_lp(i4time_sys,a9time_sys,istatus)
!            if (a9time_sys(6:9) .eq. '0000') then
!              i4time_release = i4time_release - 210600
!            else
!              i4time_release = i4time_release - 166200
!            endif
          endif

!         correction for balloon rise time to mid-troposphere
!         and time delay due to cron and sched inputs.
          i4time_raob = i4time_release + 1800 + itime_delay
        else
          i4time_release = 0
          i4time_diff = 0
        endif

        if(i4time_raob .ne. 0)then ! test window
          if(i4time_raob .ge. i4time_raob_earliest .and.
     1       i4time_raob .le. i4time_raob_latest)then
c.or.(abs(i4time_raob-i4time_sys).gt.10800))then
          else   !outside time window - reject
            write(6,*) 'raob ',wmostanum_f,' is outside of time window '
            write(6,*) 'raob:',i4time_raob,' window:',
     1                 i4time_raob_earliest, '-',i4time_raob_latest
            n_raobs_avail = n_raobs_avail - 1
            goto 888
          endif
        elseif (i4time_raob .eq. 0) then  !not a valid time
          goto 888
        endif
        i4time_diff = abs(i4time_raob-i4time_sys)
        if(i4time_diff.gt.10800)then
           print*,'raob not valid verification for this cycle'
           print*,'i4time_raob/i4time_sys: ',i4time_raob,i4time_sys
           n_raobs_avail = n_raobs_avail - 1
           goto 888
        endif
! if you get to here, raob is in domain and within time window, save it
! write what have to return arrays: staname, stalat, stalon, timesyn, timerel

        nobs = nobs + 1
        staname(nobs) = staname_f
        wmostanum(nobs) = wmostanum_f
        stalat(nobs) = stalat_f
        stalon(nobs) = stalon_f
        timesyn(nobs) = i4time_syn
        timerel(nobs) = i4time_release

! read staelev, numman, numsigt, numsigw, prsigt, tpsigt, tdsigt, htsigw, 
!               wdsigw, wssigw, prman, htman, tman, tdman, wsman, wdman 
!   into return arrays
c
c       read staelev      "station elevation"
c
        nf_status = nf_inq_varid(nf_fid,'staelev',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var staelev'
        endif
        nf_status = nf_get_var1_real(nf_fid,nf_vid,index_1,
     1                               staelev(nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var staelev'
        endif
c
c       read numman      "number of mandatory levels"
c
        nf_status = nf_inq_varid(nf_fid,'nummand',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var nummand'
        endif
        nf_status = nf_get_var1_int(nf_fid,nf_vid,index_1,numman(nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var numman'
        endif

!       read in prman and tpman
        count(1) = numman(nobs)
c
c       read prman       "pressure - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'prman',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var prman' 
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               prman(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var prman'
        endif
c
c       read missing for prman...use as missing for rest of data too
c
        nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue',
     1                              raob_missing_data)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in _fillvalue for prman'
          istatus = 0
          return
        endif
c
c       read tpman       "temperature - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'tpman',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var tpman'
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               tman(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var tpman'
        endif
c
c       read tdman       "dewpoint - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'tdman',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var tcman'
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               tdman(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var tdman'
        endif
c
c       read htman       "geopotential - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'htman',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var htman'
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               htman(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var htman'
        endif
c
c       read wdman       "wind direction - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'wdman',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var wdman'
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               wdman(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var wdman'
        endif
c
c       read wsman       "wind speed - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'wsman',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var wsman'
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               wsman(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var wsman'
        endif

c
c       read numsigt      "number of significant levels wrt t"
c
        nf_status = nf_inq_varid(nf_fid,'numsigt',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var numsigt'
        endif
        nf_status = nf_get_var1_int(nf_fid,nf_vid,index_1,numsigt(nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var numsigt'
        endif

!       read in prsigt and tpsigt
        count(1) = numsigt(nobs)
c
c       read prsigt       "pressure - significant level wrt t"
c
        nf_status = nf_inq_varid(nf_fid,'prsigt',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var prsigt' 
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               prsigti(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var prsigt'
        endif
c
c       read tpsigt       "temperature - significant level wrt t"
c
        nf_status = nf_inq_varid(nf_fid,'tpsigt',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var tpsigt'
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               tsigti(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var tpsigt'
        endif

c
c       read tdsigt       "dewpoint - significant level wrt t"
c
        nf_status = nf_inq_varid(nf_fid,'tdsigt',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var tdsigt'
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               tdsigti(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var tdsigt'
        endif

c
c       read numsigw      "number of significant levels wrt w"
c
        nf_status = nf_inq_varid(nf_fid,'numsigw',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var numsigw'
        endif
        nf_status = nf_get_var1_int(nf_fid,nf_vid,index_1,numsigw(nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var numsigw'
        endif

!       read in htsigw, wdsigw and wssigw
        count(1) = numsigw(nobs)
c
c       read htsigw       "geopotential - significant level wrt w"
c
        nf_status = nf_inq_varid(nf_fid,'htsigw',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var htsigw'
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               htsigwi(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var htsigw'
        endif
c
c       read wdsigw       "wind direction - significant level wrt w"
c
        nf_status = nf_inq_varid(nf_fid,'wdsigw',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var wdsigw'
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               wdsigwi(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var wdsigw'
        endif
c
c       read wssigw       "wind speed - significant level wrt w"
c
        nf_status = nf_inq_varid(nf_fid,'wssigw',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var wssigw'
        endif
        nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,
     1                               wssigwi(1,nobs))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status),'in var wssigw'
        endif

888     continue
      enddo
      numraob = nobs
c
c     close netcdf file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close ',filename
      endif

c     convert dewpoint depression to dewpoint
      do i = 1, numraob
        do j = 1, numman(i)
          if ((tdman(j,i) .ne. raob_missing_data) .and.
     1        (tman(j,i) .ne. raob_missing_data))
     1      tdman(j,i) = tman(j,i) - tdman(j,i)
        enddo
 
        do j = 1, numsigt(i)
          if ((tdsigti(j,i) .ne. raob_missing_data) .and.
     1        (tsigti(j,i) .ne. raob_missing_data))
     1      tdsigti(j,i) = tsigti(j,i) - tdsigti(j,i)
        enddo
      enddo

          

c
c     interleave man obs into sigt and sigw obs
c
      call minglemansig(maxm,maxt,maxw,maxraob,numraob,
     1                  max_ht_m_proc, min_pres_mb_proc,
     1                  numsigt,numsigw,numman,typew,typet,
     1                  htman,prman,tman,tdman,wsman,wdman,
     1                  staelev,raob_missing_data, 
     1                  prsigti,tsigti,tdsigti,htsigwi,
     1                  wssigwi,wdsigwi,prsigt,tsigt,tdsigt,
     1                  htsigt,htsigw,wssigw,wdsigw,
     1                  verif_missing_data,istatus)

      if (istatus .ne. 1) then
        write(6,*) 'error mingling man and sig obs '
        istatus = 0
      endif

      return
      end
!1............................................................................
      subroutine minglemansig(maxm,maxt,maxw,maxraob,numraob,
     1                  max_ht_m_proc, min_pres_mb_proc,
     1                  numsigt,numsigw,numman,typew,typet,
     1                  htman,prman,tman,tdman,wsman,wdman,
     1                  staelev,raob_missing_data,
     1                  prsigti,tsigti,tdsigti,htsigwi,
     1                  wssigwi,wdsigwi,prsigt,tsigt,tdsigt,
     1                  htsigt,htsigw,wssigw,wdsigw,
     1                  verif_missing_data,istatus)

      implicit none

      integer     maxm,maxt,maxw,maxraob,numraob,
     1              numsigt(maxraob),numsigw(maxraob),
     1              numman(maxraob)
      real	    max_ht_m_proc, min_pres_mb_proc
      character*1   typew(maxw,maxraob),typet(maxt,maxraob)
      real        htman(maxm,maxraob), prman(maxm,maxraob),
     1              tman(maxm,maxraob), tdman(maxm,maxraob),
     1              wsman(maxm,maxraob), wdman(maxm,maxraob),
     1              staelev(maxraob), raob_missing_data,
     1              prsigti(maxt,maxraob),tsigti(maxt,maxraob),
     1              tdsigti(maxt,maxraob),
     1              htsigwi(maxw,maxraob),wdsigwi(maxw,maxraob),
     1              wssigwi(maxw,maxraob)
      real        prsigt(maxt,maxraob),tsigt(maxt,maxraob),
     1              tdsigt(maxt,maxraob),
     1              htsigt(maxt,maxraob),
     1              htsigw(maxw,maxraob),wdsigw(maxw,maxraob),
     1              wssigw(maxw,maxraob), htsfc, prsfc,
     1              verif_missing_data
      integer	    istatus

      integer     mptr, sptr, jptr, i, j

c     begin

c     loop through raobs to interleave man with sigt and man with sigw
      do i = 1, numraob

c       debug  print out data
        if(.true.)then

        write(6,*) 'raob ',i
        write(6,*) 'htman, wsman, wdman'
        do j = 1, numman(i)
          write(6,*) j,htman(j,i),wsman(j,i),wdman(j,i)
        enddo

        write(6,*) 'htsigwi, wssigwi, wdsigwi'
        do j = 1, numsigw(i)
          write(6,*) j,htsigwi(j,i),wssigwi(j,i),wdsigwi(j,i)
        enddo

        write(6,*) 'prman, tman, tdman'
        do j = 1, numman(i)
          write(6,*) j,prman(j,i),tman(j,i),tdman(j,i)
        enddo

        write(6,*) 'prsigti, tsigti, tdsigti'
        do j = 1, numsigt(i)
          write(6,*) j,prsigti(j,i),tsigti(j,i),tdsigti(j,i)
        enddo
        write(6,*)

        endif

c       get surface ht and pr setup
        jptr = 1
        mptr = 1
        sptr = 1

c       see if first level of man is missing
        if (htman(1,i) .eq. raob_missing_data) then  !see if it's in htsigw
          if ((htsigwi(1,i) .eq. raob_missing_data) .or.
     1        ((htsigwi(1,i) .eq. 0.0) .and.
     1         (wssigwi(1,i) .eq. raob_missing_data) .and.
     1         (wdsigwi(1,i) .eq. raob_missing_data))) then  !skip surface level
            htsfc = raob_missing_data 
            mptr = 2
            sptr = 2
          else
            htsfc = htsigwi(1,i)
            htsigw(jptr,i) = htsigwi(sptr,i)
            wssigw(jptr,i) = wssigwi(sptr,i)
            wdsigw(jptr,i) = wdsigwi(sptr,i)
            typew(jptr,i) = 'w'
            mptr = mptr + 1  !skip sfc in man
            sptr = sptr + 1
            jptr = jptr + 1
          endif
        else
          htsfc = htman(1,i)
          htsigw(jptr,i) = htman(mptr,i)
          wssigw(jptr,i) = wsman(mptr,i)
          wdsigw(jptr,i) = wdman(mptr,i)
          typew(jptr,i) = 'm'
          mptr = mptr + 1
          sptr = sptr + 1  !skip sfc in sigw
          jptr = jptr + 1
        endif

c       set mptr and sptr so ht .gt. htsfc 
        do while ((mptr .le. numman(i)) .and.
     1            (htman(mptr,i) .le. htsfc))
          mptr = mptr + 1
        enddo
        do while ((sptr .le. numsigw(i)) .and.
     1            (htsigwi(sptr,i) .le. htsfc))
          sptr = sptr + 1
        enddo

c       mingle man with sigw - go until sigw levels exhausted
        do while (sptr .le. numsigw(i))  

c         make sure mptr is past any missing htman
          do while ((mptr .le. numman(i)) .and.
     1          ((htman(mptr,i) .gt. max_ht_m_proc)
     1      .or. (htman(mptr,i) .eq. raob_missing_data)
     1      .or. (wsman(mptr,i) .eq. raob_missing_data) 
     1      .or. (wdman(mptr,i) .eq. raob_missing_data)))
            mptr = mptr + 1
          enddo
       
c         make sure sptr is past any missing htsigw
          do while ((sptr .le. numsigw(i)) .and.
     1          ((htsigwi(sptr,i) .gt. max_ht_m_proc)
     1      .or. (htsigwi(sptr,i) .eq. raob_missing_data)
     1      .or. (wssigwi(mptr,i) .eq. raob_missing_data) 
     1      .or. (wdsigwi(mptr,i) .eq. raob_missing_data)))
            sptr = sptr + 1
          enddo

          if (sptr .gt. numsigw(i)) goto 777  !no sigw obs to mingle
       
          if (mptr .le. numman(i)) then
            if (htsigwi(sptr,i) .lt. htman(mptr,i)) then
              htsigw(jptr,i) = htsigwi(sptr,i)
              wssigw(jptr,i) = wssigwi(sptr,i)
              wdsigw(jptr,i) = wdsigwi(sptr,i)
              typew(jptr,i) = 'w'
              sptr = sptr + 1
              jptr = jptr + 1
            else 
              if (htman(mptr,i) .lt. htsigwi(sptr,i)) then
                htsigw(jptr,i) = htman(mptr,i)
                wssigw(jptr,i) = wsman(mptr,i)
                wdsigw(jptr,i) = wdman(mptr,i)
                typew(jptr,i) = 'm'
                mptr = mptr + 1
                jptr = jptr + 1
              else
                if (htman(mptr,i) .eq. htsigwi(sptr,i)) then
                  htsigw(jptr,i) = htman(mptr,i)
                  wssigw(jptr,i) = wsman(mptr,i)
                  wdsigw(jptr,i) = wdman(mptr,i)
                  typew(jptr,i) = 'm'
                  mptr = mptr + 1
                  sptr = sptr + 1
                  jptr = jptr + 1
                endif
              endif
            endif
          else
            if (htsigw(jptr,i) .le. max_ht_m_proc) then
              htsigw(jptr,i) = htsigwi(sptr,i)
              wssigw(jptr,i) = wssigwi(sptr,i)
              wdsigw(jptr,i) = wdsigwi(sptr,i)
              typew(jptr,i) = 'w'
              sptr = sptr + 1
              jptr = jptr + 1
            endif
          endif

777       continue 

        enddo

c       finish rest of man if any are left
        do while (mptr .le. numman(i)) 
          if ((htman(mptr,i) .ne. raob_missing_data) 
     1      .and. (htman(mptr,i) .le. max_ht_m_proc) 
     1      .and. (wsman(mptr,i) .ne. raob_missing_data) 
     1      .and. (wdman(mptr,i) .ne. raob_missing_data))
     1      then
            htsigw(jptr,i) = htman(mptr,i)
            wssigw(jptr,i) = wsman(mptr,i)
            wdsigw(jptr,i) = wdman(mptr,i)
            typew(jptr,i) = 'm'
            jptr = jptr + 1
          endif
          mptr = mptr + 1
        enddo

c       set numsigw to mingled value
        numsigw(i) = jptr - 1

        write(6,*) numsigw(i)
        write(6,*) 'htsigw, typew, wssigw, wdsigw '        
        do j = 1, numsigw(i)
          write(6,*) j, htsigw(j,i), typew(j,i), wssigw(j,i), 
     1               wdsigw(j,i)         
        enddo

c       mingle man with sigt
c       get surface pr setup
        jptr = 1
        mptr = 1
        sptr = 1

c       see if first level of man is missing
        if (prman(1,i) .eq. raob_missing_data) then  !see if it's in prsigt
          if ((prsigti(1,i) .eq. raob_missing_data) .or.
     1        (tsigti(1,i) .eq. raob_missing_data))then  !skip surface level
            prsfc = raob_missing_data 
            mptr = 2
            sptr = 2
          else
            prsfc = prsigti(1,i)
            prsigt(jptr,i) = prsigti(sptr,i)
            tsigt(jptr,i) = tsigti(sptr,i)
            tdsigt(jptr,i) = tdsigti(sptr,i)
            htsigt(jptr,i) = verif_missing_data
            typet(jptr,i) = 't'
            mptr = mptr + 1  !skip sfc in man
            sptr = sptr + 1
            jptr = jptr + 1
          endif
        else
          prsfc = prman(1,i)
          prsigt(jptr,i) = prman(mptr,i)
          tsigt(jptr,i) = tman(mptr,i)
          tdsigt(jptr,i) = tdman(mptr,i)
          htsigt(jptr,i) = htman(mptr,i)
          typet(jptr,i) = 'm'
          mptr = mptr + 1
          sptr = sptr + 1  !skip sfc in sigt
          jptr = jptr + 1
        endif

c       set mptr and sptr so pr .lt. prsfc
        do while ((mptr .le. numman(i)) .and.
     1            (prman(mptr,i) .ge. prsfc))
          mptr = mptr + 1
        enddo
        do while ((sptr .le. numsigt(i)) .and.
     1            (prsigti(sptr,i) .ge. prsfc))
          sptr = sptr + 1
        enddo

c       mingle man with sigt - go until sigt levels exhausted
        do while (sptr .le. numsigt(i))  

c         make sure mptr is past any missing prman
          do while ((mptr .le. numman(i)) .and.
     1          ((prman(mptr,i) .lt. min_pres_mb_proc)
     1      .or. (prman(mptr,i) .eq. raob_missing_data)
     1      .or. (tman(mptr,i) .eq. raob_missing_data)))
            mptr = mptr + 1
          enddo
       
c         make sure sptr is past any missing prsigt
          do while ((sptr .le. numsigt(i)) .and.
     1         ((prsigti(sptr,i) .lt. min_pres_mb_proc)
     1      .or.(prsigti(sptr,i) .eq. raob_missing_data)
     1      .or.(tsigti(mptr,i) .eq. raob_missing_data)))
            sptr = sptr + 1
          enddo

          if (sptr .gt. numsigt(i)) goto 778

          if (mptr .le. numman(i)) then
            if (prsigti(sptr,i) .gt. prman(mptr,i)) then
              prsigt(jptr,i) = prsigti(sptr,i)
              tsigt(jptr,i) = tsigti(sptr,i)
              tdsigt(jptr,i) = tdsigti(sptr,i)
              htsigt(jptr,i) = verif_missing_data
              typet(jptr,i) = 't'
              sptr = sptr + 1
              jptr = jptr + 1
            else 
              if (prman(mptr,i) .gt. prsigti(sptr,i)) then
                prsigt(jptr,i) = prman(mptr,i)
                tsigt(jptr,i) = tman(mptr,i)
                tdsigt(jptr,i) = tdman(mptr,i)
                htsigt(jptr,i) = htman(mptr,i)
                typet(jptr,i) = 'm'
                mptr = mptr + 1
                jptr = jptr + 1
              else
                if (prman(mptr,i) .eq. prsigti(sptr,i)) then
                  prsigt(jptr,i) = prman(mptr,i)
                  tsigt(jptr,i) = tman(mptr,i)
                  tdsigt(jptr,i) = tdman(mptr,i)
                  htsigt(jptr,i) = htman(mptr,i)
                  typet(jptr,i) = 'm'
                  mptr = mptr + 1
                  sptr = sptr + 1
                  jptr = jptr + 1
                endif
              endif
            endif
          else
            prsigt(jptr,i) = prsigti(sptr,i)
            tsigt(jptr,i) = tsigti(sptr,i)
            tdsigt(jptr,i) = tdsigti(sptr,i)
            htsigt(jptr,i) = verif_missing_data
            typet(jptr,i) = 't'
            sptr = sptr + 1
            jptr = jptr + 1
          endif

778       continue

        enddo

c       finish rest of man if any are left
        if (mptr .le. numman(i)) then
          do while (mptr .le. numman(i)) 
            if ((prman(mptr,i) .ne. raob_missing_data).and. 
     1          (tman(mptr,i) .ne. raob_missing_data)) then
              prsigt(jptr,i) = prman(mptr,i)
              tsigt(jptr,i) = tman(mptr,i)
              tdsigt(jptr,i) = tdman(mptr,i)
              htsigt(jptr,i) = htman(mptr,i)
              typet(jptr,i) = 'm'
              jptr = jptr + 1
            endif
            mptr = mptr + 1
          enddo
        endif

c       set numsigt to mingled value
        numsigt(i) = jptr - 1

        write(6,*) 'raob ',i,' after mingle'
        write(6,*) numsigt(i)
        write(6,*) 'prsigt, typet tsigt, tdsigt htsigt'        
        do j = 1, numsigt(i)
          write(6,*) j,prsigt(j,i),typet(j,i),tsigt(j,i), 
     1               tdsigt(j,i), htsigt(j,i) 
        enddo

      enddo

      return
      end
!2............................................................................
      subroutine calc_domain_perim(nx, ny, lat, lon, north, south,
     1                             east,west,r_buffer)

      integer     nx, ny, i, j
      real        north, south, east, west, r_buffer,
     1              lat(nx,ny), lon(nx,ny)

!     calculate domain perimeter
        north = -90.
        south  = +90.
        west = +1000.
        east = -1000.

        do i = 1,nx
          north = max(north,lat(i,1),lat(i,ny))
          south = min(south,lat(i,1),lat(i,ny))
          east  = max(east ,lon(i,1),lon(i,ny))
          west  = min(west ,lon(i,1),lon(i,ny))
        enddo ! i

        do j = 1,ny
          north = max(rnorth,lat(1,j),lat(nx,j))
          south = min(south ,lat(1,j),lat(nx,j))
          east  = max(east  ,lon(1,j),lon(nx,j))
          west  = min(west  ,lon(1,j),lon(nx,j))
        enddo ! j

        north = north + r_buffer
        south = south - r_buffer
        east  = east  + r_buffer
        west  = west  - r_buffer

       write(6,101)north,south,east,west
101     format(1x,' box around laps grid - nsew ',4f9.2)

       return
       end
!.............................................................................
