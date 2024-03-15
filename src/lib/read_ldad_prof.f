
      subroutine read_ldad_prof(i4time_sys,i4_prof_window
     1                                    ,nx_l,ny_l
     1                                    ,ext,lun
     1                                    ,filename,n_good_obs,istatus)

      character*(*) filename,ext

!.............................................................................

      include 'netcdf.inc'
      integer level, maxstaticids, ninventorybins, recnum,nf_fid,
     +     nf_vid, nf_status

      write(6,*)
      write(6,*)' subroutine read_ldad_prof...'

c
c  open netcdf file for reading
c
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open filename'
      endif
c
c  fill all dimension values
c
c
c get size of level
c
      nf_status = nf_inq_dimid(nf_fid,'level',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim level'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,level)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim level'
      endif
c
c get size of maxstaticids
c
      nf_status = nf_inq_dimid(nf_fid,'maxstaticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxstaticids'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,maxstaticids)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxstaticids'
      endif
c
c get size of ninventorybins
c
      nf_status = nf_inq_dimid(nf_fid,'ninventorybins',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim ninventorybins'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,ninventorybins)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim ninventorybins'
      endif
c
c get size of recnum
c
      nf_status = nf_inq_dimid(nf_fid,'recnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,recnum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
      endif
      call read_prof(nf_fid, level, maxstaticids, ninventorybins,
     +     recnum,ext,lun,
!.............................................................................
     1     i4time_sys,i4_prof_window,nx_l,ny_l,n_good_obs,istatus)

      return
!.............................................................................

      end
c
c
      subroutine read_prof(nf_fid, level, maxstaticids, ninventorybins,       
     +     recnum,ext,lun,
!.............................................................................
     1     i4time_sys,i4_prof_window,nx_l,ny_l,n_good_obs,istatus)
!.............................................................................

      include 'netcdf.inc'
      integer level, maxstaticids, ninventorybins, recnum,nf_fid,
     +     nf_vid, nf_status
      integer assetid(recnum), averageminutes(recnum), 
     +     firstinbin(ninventorybins),
     +     firstoverflow, globalinventory, invtime(recnum),
     +     inventory(maxstaticids), isoverflow(recnum),
     +     lastinbin(ninventorybins), lastrecord(maxstaticids),
     +     nstaticids, prevrecord(recnum), tempqcflag( level,
     +     recnum), wdqcflag( level, recnum), winddir( level,
     +     recnum), wsqcflag( level, recnum)
      real elevation(recnum), latitude(recnum), levels( level,
     +     recnum), longitude(recnum), temperature( level, recnum),
     +     windspeed( level, recnum)
      double precision observationtime(recnum), receipttime(recnum),
     +     reporttime(recnum)
      character*6 stationid(recnum)
      character*4 homewfo(recnum)
      character*24 dataprovider(recnum)
      character*6 providerid(recnum)
      integer stationtype(recnum)
      character*30 staticids(maxstaticids)
      character*51 stationname(recnum)
!.............................................................................

      character*9 a9_timeobs,a9_recpttime,a9_closest,a9time_ob
      character*8 c8_project,c8_format
      character*6 provider_ref
      character*(*)ext
      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)
      real ht_out(200),di_out(200),sp_out(200),temp_out(200)
      integer iqc1_out(200),iqc2_out(200)
      integer assetid_ref
      real mspkt
      data mspkt/.518/

!............................................................................

!     initialize arrays
      averageminutes = -999
      stationid = 'unk   '

      call get_c8_project(c8_project,istatus)

      if(c8_project .eq. 'rsa')then
          c8_format = 'ldad'
      else
          c8_format = 'madis'
      endif

      write(6,*)' reading profiler/rass data, format = ',c8_format

      call read_ldad_prof_netcdf(nf_fid, level, maxstaticids, 
     +     ninventorybins,       
     +     recnum, assetid, averageminutes, firstinbin, firstoverflow, 
     +     globalinventory, invtime, inventory, isoverflow, 
     +     lastinbin, lastrecord, nstaticids, prevrecord, tempqcflag, 
     +     wdqcflag, winddir, wsqcflag, elevation, latitude, levels, 
     +     longitude, temperature, windspeed, observationtime, 
     +     receipttime, reporttime, dataprovider, homewfo, 
     +     providerid, staticids, stationid, stationname, stationtype,
     +     c8_format)
c
c the netcdf variables are filled - your code goes here
c
!............................................................................

      call get_latlon_perimeter(nx_l,ny_l,1.0
     1                         ,lat_a,lon_a,topo_a
     1                         ,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error in get_latlon_perimeter'
          return
      endif

      write(6,*)' # of profilers = ',nstaticids
      write(6,*)' # of records = ',recnum

      do i_sta = 1,nstaticids
        ilast_rec = lastrecord(i_sta) + 1 ! offset going from c to fortran
        provider_ref = providerid(ilast_rec)

        call s_len(staticids(i_sta),len_id)

        write(6,*)
        write(6,*)' looping for profiler ',i_sta
     1                                    ,staticids(i_sta)(1:len_id)

        write(6,*)' station / last record / provider '
     1           ,i_sta,ilast_rec,provider_ref         

        i4_resid_closest = 999999
        i4_avg_window = -999
        a9_closest = '---------'

        do irec = 1,recnum

          call s_len(providerid(irec),len_prov)

          if(providerid(irec) .eq. provider_ref)then

              write(6,*)
              write(6,*)' providerid match = ',irec
     1                 ,providerid(irec)(1:len_prov)

              rlat = latitude(irec)
              rlon = longitude(irec)

              if(rlat .le. rnorth .and. rlat .ge. south .and.
     1           rlon .ge. west   .and. rlon .le. east            )then        

                  write(6,*)irec,providerid(irec)(1:len_prov)
     1                     ,' is in box'       

                  elev = elevation(irec)

!                 convert u_std, v_std to rms

!                 test observation time
                  if(abs(observationtime(irec))      .lt. 3d9)then
                      ictime_ob = nint(observationtime(irec))

                      if(averageminutes(irec) .ne. -999)then
                          i4_avg_window = averageminutes(irec) * 60
                          i4_avg_window_half = i4_avg_window / 2
                          ictime_ob = ictime_ob - i4_avg_window_half
                      endif

                      call c_time2fname(ictime_ob,a9_timeobs)

                  else
                      write(6,*)' bad observation time - reject record'         
     1                           ,observationtime(irec)
                      goto 300

                  endif

!                 test number of good levels
                  n_good_levels = 0

                  if(ext(1:3) .eq. 'pro')then ! test wind profile
                      do i = 1,level
                          if(winddir(i,irec) .ge. 0      .and.
     1                       winddir(i,irec) .le. 360    .and.
     1                 iqc_rsa(wdqcflag(i,irec)) .ne. -1 .and.       
     1                 iqc_rsa(wsqcflag(i,irec)) .ne. -1     
     1                                                          )then ! good qc
                              n_good_levels = n_good_levels + 1
                          endif
                      enddo ! i

                  elseif(ext(1:3) .eq. 'lrs')then ! test temperature profile
                      do i = 1,level
                          if(iqc_rsa(tempqcflag(i,irec)) .ne. -1
     1                            .and.
     1                        temperature(i,irec) .gt. 200.
     1                            .and.
     1                        temperature(i,irec) .lt. 400.
     1                                                   )then ! good qc
                              n_good_levels = n_good_levels + 1
                          endif
                      enddo ! i
                  endif

                  if(n_good_levels .le. 0)then
                      write(6,*)' no good levels - reject record'         
     1                           ,observationtime(irec)
                      goto 300

                  else ! good levels detected
                      write(6,*)' good levels = ',n_good_levels
     1                           ,observationtime(irec)
                  endif

!                 determine if this report is closest to the analysis time
                  call cv_asc_i4time(a9_timeobs,i4time_ob)
                  i4_resid = abs(i4time_ob - i4time_sys)
                  if(i4_resid .lt. i4_resid_closest)then
                      i4_resid_closest = i4_resid
                      i_pr_cl = irec
                      a9_closest = a9_timeobs
                  endif
                  write(6,*)'i4_resid/closest/avg_window = '
     1                      ,i4_resid,i4_resid_closest,i4_avg_window       

              else !
                  write(6,*)irec,providerid(irec)(1:len_prov)
     1                     ,' is outside of domain'

                  go to 900 ! loop back to next station

              endif ! in box

!         else
!             write(6,*)' providerid = ',irec
!    1                 ,providerid(irec)(1:len_prov)
 
          endif ! correct provider id

 300      continue

        enddo ! irec 

        write(6,*)
        write(6,*)' evaluating profile closest in time'

        if(i4_resid_closest .gt. i4_prof_window)then ! outside time window
            write(6,*)' outside time window - reject '
     1               ,a9_closest,i4_resid,i4_prof_window
        
        else
!           lun=1 ! for both 'pro' and 'lrs'
c
c           open intermediate output file.
c
            call open_ext(lun,i4time_sys,ext(1:3),istatus)
            if(istatus .ne. 1)then
                write(6,*)' error opening product file',ext
                goto980
            endif

            n_good_obs = n_good_obs + 1

            a9time_ob = a9_closest

            n_good_levels = 0

            call filter_string(stationid(i_pr_cl))

            if(ext(1:3) .eq. 'pro')then

                rms = 1.0

                write(6,*)'i/winddir/wdqcflag/wsqcflag, i_pr_cl='
     1                                                 ,i_pr_cl     

                do i = 1,level
                    write(6,*)i,winddir(i,i_pr_cl)
     1                         ,wdqcflag(i,i_pr_cl),wsqcflag(i,i_pr_cl)       
                    if(winddir(i,i_pr_cl) .ge. 0      .and.
     1                 winddir(i,i_pr_cl) .le. 360    .and.
     1                 iqc_rsa(wdqcflag(i,i_pr_cl)) .ne. -1 .and.       
     1                 iqc_rsa(wsqcflag(i,i_pr_cl)) .ne. -1     
     1                                                          )then ! good qc
                        n_good_levels = n_good_levels + 1
                        ht_out(n_good_levels) = levels(i,i_pr_cl)
                        di_out(n_good_levels) = winddir(i,i_pr_cl)
                        sp_out(n_good_levels) = windspeed(i,i_pr_cl)
     1                                        * mspkt
                        iqc1_out(n_good_levels) = wdqcflag(i,i_pr_cl)     
                        iqc2_out(n_good_levels) = wsqcflag(i,i_pr_cl)     
                    endif
                enddo ! i

                call filter_string(provider_ref)

                if(c8_format .eq. 'ldad')then
                    write(6,401)provider_ref
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,stationid(i_pr_cl)(1:6)
     1                         ,a9time_ob,'profiler'
                    write(lun,401)provider_ref
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,stationid(i_pr_cl)(1:6)
     1                         ,a9time_ob,'profiler'
                else ! madis
                    write(6,402)i_sta
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,provider_ref(1:6)
     1                         ,a9time_ob,'profiler'
                    write(lun,402)i_sta
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,provider_ref(1:6)
     1                         ,a9time_ob,'profiler'
                endif

401             format(a12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)
402             format(i12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)

                do i = 1,n_good_levels
                    write(lun,411,err=421)ht_out(i)
     1                                   ,di_out(i),sp_out(i)
     1                                   ,rms
                    write(6  ,411,err=421)ht_out(i)
     1                                   ,di_out(i),sp_out(i)
     1                                   ,rms
     1                                   ,iqc1_out(i)
     1                                   ,iqc2_out(i)
411                 format(1x,f6.0,f6.0,2f6.1,2i7)
421                 continue
                enddo ! i

            elseif(ext(1:3) .eq. 'lrs')then

                rms = 1.0
                iqc = 1

                do i = 1,level
                    if(iqc_rsa(tempqcflag(i,i_pr_cl)) .ne. -1
     1                            .and.
     1                 temperature(i,i_pr_cl) .gt. 200.
     1                            .and.
     1                 temperature(i,i_pr_cl) .lt. 400.
     1                                                   )then ! good qc
                        n_good_levels = n_good_levels + 1
                        ht_out(n_good_levels) = levels(i,i_pr_cl)
                        temp_out(n_good_levels) = temperature(i,i_pr_cl)
                        iqc1_out(n_good_levels) = tempqcflag(i,i_pr_cl)
                    endif
                enddo ! i

                call filter_string(provider_ref)

                if(c8_format .eq. 'ldad')then
                    write(6,501)provider_ref
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,stationid(i_pr_cl)(1:5)
     1                         ,a9time_ob,'rass    '
                    write(lun,501)provider_ref
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,stationid(i_pr_cl)(1:5)
     1                         ,a9time_ob,'rass    '
                else ! madis
                    write(6,502)i_sta
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,provider_ref(1:5)
     1                         ,a9time_ob,'rass    '
                    write(lun,502)i_sta
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,provider_ref(1:5)
     1                         ,a9time_ob,'rass    '
                endif

501             format(a12,i12,f11.3,f15.3,f15.0,5x,a5,3x,a9,1x,a8)
502             format(i12,i12,f11.3,f15.3,f15.0,5x,a5,3x,a9,1x,a8)

                do i = 1,n_good_levels
                    write(lun,511,err=521)ht_out(i)
     1                                   ,temp_out(i)
     1                                   ,iqc
     1                                   ,rms
                    write(6  ,511,err=521)ht_out(i)
     1                                   ,temp_out(i)
     1                                   ,iqc1_out(i)
     1                                   ,rms
511                 format(1x,f6.0,f6.1,i6,f6.1)
521                 continue
                enddo ! i

            endif ! ext

        endif ! in time window

980     continue

900   enddo ! ista

!............................................................................
      return
      end
c
c  subroutine to read the file 
c
      subroutine read_ldad_prof_netcdf(nf_fid, level, maxstaticids, 
     +     ninventorybins, recnum, assetid, averageminutes, firstinbin,        
     +     firstoverflow, globalinventory, invtime, inventory, 
     +     isoverflow, lastinbin, lastrecord, nstaticids, prevrecord, 
     +     tempqcflag, wdqcflag, winddir, wsqcflag, elevation, 
     +     latitude, levels, longitude, temperature, windspeed, 
     +     observationtime, receipttime, reporttime, dataprovider, 
     +     homewfo, providerid, staticids, stationid, stationname, 
     +     stationtype,c8_format)
c
      include 'netcdf.inc'
      integer level, maxstaticids, ninventorybins, recnum,nf_fid, 
     +     nf_vid, nf_status
      integer assetid(recnum), averageminutes(recnum), 
     +     firstinbin(ninventorybins),
     +     firstoverflow, globalinventory, invtime(recnum),
     +     inventory(maxstaticids), isoverflow(recnum),
     +     lastinbin(ninventorybins), lastrecord(maxstaticids),
     +     nstaticids, prevrecord(recnum), tempqcflag( level,
     +     recnum), wdqcflag( level, recnum), winddir( level,
     +     recnum), wsqcflag( level, recnum)
      real elevation(recnum), latitude(recnum), levels( level,
     +     recnum), longitude(recnum), temperature( level, recnum),
     +     windspeed( level, recnum)
      double precision observationtime(recnum), receipttime(recnum),
     +     reporttime(recnum)
      character*6 stationid(recnum)
      character*4 homewfo(recnum)
      character*24 dataprovider(recnum)
      character*6 providerid(recnum)
      integer stationtype(recnum)
      character*30 staticids(maxstaticids)
      character*51 stationname(recnum)

      character*8 c8_format


c   variables of type real
c
c     variable        netcdf long name
c      elevation    "elevation above msl"
c
      nf_status = nf_inq_varid(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevation'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,elevation)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var elevation'
        endif
      endif
c
c     variable        netcdf long name
c      latitude     "station latitude"
c
      nf_status = nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,latitude)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var latitude'
        endif
      endif
c
c     variable        netcdf long name
c      levels       "instrument level, height above station"
c
      nf_status = nf_inq_varid(nf_fid,'levels',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var levels'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,levels)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var levels'
        endif
      endif
c
c     variable        netcdf long name
c      longitude    "station longitude"
c
      nf_status = nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,longitude)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var longitude'
        endif
      endif
c
c     variable        netcdf long name
c      temperature  "temperature"
c
      nf_status = nf_inq_varid(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,temperature)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var temperature'
        endif
      endif
c
c     variable        netcdf long name
c      windspeed    "wind speed (scalar)"
c
      nf_status = nf_inq_varid(nf_fid,'windspeed',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeed'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,windspeed)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var windspeed'
        endif
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      assetid      "rsa asset identifier"
c
!       nf_status = nf_inq_varid(nf_fid,'assetid',nf_vid)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'in var assetid'
!     endif
!       nf_status = nf_get_var_int(nf_fid,nf_vid,assetid)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'in var assetid'
!     endif
c
c     variable        netcdf long name
c      firstinbin   
c
      nf_status = nf_inq_varid(nf_fid,'firstinbin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstinbin'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,firstinbin)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var firstinbin'
        endif
      endif
c
c     variable        netcdf long name
c      averageminutes
c
      nf_status = nf_inq_varid(nf_fid,'averageminutes',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var averageminutes'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,averageminutes)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var averageminutes'
        endif
      endif
c
c     variable        netcdf long name
c      firstoverflow
c
      nf_status = nf_inq_varid(nf_fid,'firstoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstoverflow'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,firstoverflow)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var firstoverflow'
        endif
      endif
c
c     variable        netcdf long name
c      globalinventory
c
      nf_status = nf_inq_varid(nf_fid,'globalinventory',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var globalinventory'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,globalinventory)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var globalinventory'
        endif
      endif
c
c     variable        netcdf long name
c      invtime      
c
      nf_status = nf_inq_varid(nf_fid,'invtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var invtime'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,invtime)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var invtime'
        endif
      endif
c
c     variable        netcdf long name
c      inventory    
c
      nf_status = nf_inq_varid(nf_fid,'inventory',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var inventory'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,inventory)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var inventory'
        endif
      endif
c
c     variable        netcdf long name
c      isoverflow   
c
      nf_status = nf_inq_varid(nf_fid,'isoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var isoverflow'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,isoverflow)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var isoverflow'
        endif
      endif
c
c     variable        netcdf long name
c      lastinbin    
c
      nf_status = nf_inq_varid(nf_fid,'lastinbin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lastinbin'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,lastinbin)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var lastinbin'
        endif
      endif
c
c     variable        netcdf long name
c      lastrecord   
c
      nf_status = nf_inq_varid(nf_fid,'lastrecord',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lastrecord'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,lastrecord)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var lastrecord'
        endif
      endif
c
c     variable        netcdf long name
c      nstaticids   
c
      nf_status = nf_inq_varid(nf_fid,'nstaticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nstaticids'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,nstaticids)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var nstaticids'
        endif
      endif
c
c     variable        netcdf long name
c      prevrecord   
c
      nf_status = nf_inq_varid(nf_fid,'prevrecord',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prevrecord'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,prevrecord)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var prevrecord'
        endif
      endif
c
c     variable        netcdf long name
c      tempqcflag   "rsa temperature quality control flag"
c
      nf_status = nf_inq_varid(nf_fid,'tempqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tempqcflag'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,tempqcflag)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var tempqcflag'
        endif
      endif
c
c     variable        netcdf long name
c      wdqcflag     "rsa wind direction quality control flag"
c
      nf_status = nf_inq_varid(nf_fid,'wdqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdqcflag'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,wdqcflag)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var wdqcflag'
        endif
      endif
c
c     variable        netcdf long name
c      winddir      "wind direction (scalar)"
c
      nf_status = nf_inq_varid(nf_fid,'winddir',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,winddir)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var winddir'
        endif
      endif
c
c     variable        netcdf long name
c      wsqcflag     "rsa wind speed quality control flag"
c
      nf_status = nf_inq_varid(nf_fid,'wsqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wsqcflag'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,wsqcflag)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var wsqcflag'
        endif
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c      observationtime"time of observation"
c
      nf_status = nf_inq_varid(nf_fid,'observationtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var observationtime'
      else
        nf_status = nf_get_var_double(nf_fid,nf_vid,observationtime)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var observationtime'
        endif
      endif
c
c     variable        netcdf long name
c      receipttime  "file time stamp (time file was received)"
c
!       nf_status = nf_inq_varid(nf_fid,'receipttime',nf_vid)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'in var receipttime'
!     endif
!       nf_status = nf_get_var_double(nf_fid,nf_vid,receipttime)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'in var receipttime'
!     endif
c
c     variable        netcdf long name
c      reporttime   "time of observation"
c
!       nf_status = nf_inq_varid(nf_fid,'reporttime',nf_vid)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'in var reporttime'
!     endif
!       nf_status = nf_get_var_double(nf_fid,nf_vid,reporttime)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'in var reporttime'
!      endif


c   variables of type char
c
c
c     variable        netcdf long name
c      dataprovider "name of organization responsible for delivering the data"
c
      nf_status = nf_inq_varid(nf_fid,'dataprovider',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dataprovider'
      else
        nf_status = nf_get_var_text(nf_fid,nf_vid,dataprovider)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var dataprovider'
        endif
      endif
c
c     variable        netcdf long name
c      homewfo      "home wfo id"
c
!       nf_status = nf_inq_varid(nf_fid,'homewfo',nf_vid)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'in var homewfo'
!     endif
!       nf_status = nf_get_var_text(nf_fid,nf_vid,homewfo)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'in var homewfo'
!     endif
c
c     variable        netcdf long name
c      providerid   "alphanumeric station name"
c
      nf_status = nf_inq_varid(nf_fid,'providerid',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var providerid'
      else
        nf_status = nf_get_var_text(nf_fid,nf_vid,providerid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var providerid'
        endif
      endif
c
c     variable        netcdf long name
c      staticids    
c
      nf_status = nf_inq_varid(nf_fid,'staticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staticids'
      else
        nf_status = nf_get_var_text(nf_fid,nf_vid,staticids)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var staticids'
        endif
      endif
c
c     variable        netcdf long name
c      stationid    "alphanumeric station id"
c
      nf_status = nf_inq_varid(nf_fid,'stationid',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stationid'
      else
        nf_status = nf_get_var_text(nf_fid,nf_vid,stationid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var stationid'
        endif
      endif
c
c     variable        netcdf long name
c      stationname  "alphanumeric station name"
c
!       nf_status = nf_inq_varid(nf_fid,'stationname',nf_vid)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'in var stationname'
!     endif
!       nf_status = nf_get_var_text(nf_fid,nf_vid,stationname)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'in var stationname'
!     endif


!     for rsa, this variable doesn't contain any useful information.
!     this could be useful if 'c8_format' is madis, if we read it in 
!     given its "short" declaration as an int variable.

c
c     variable        netcdf long name
c      stationtype  "ldad station type"
c
      if(c8_format .eq. 'madis')then
        write(6,*)' read stationtype as short/integer'
        nf_status = nf_inq_varid(nf_fid,'stationtype',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var stationtype'
        else
          nf_status = nf_get_var_int(nf_fid,nf_vid,stationtype)
          if(nf_status.ne.nf_noerr) then
            print *, nf_strerror(nf_status)
            print *,'in var stationtype'
          endif
        endif

      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end


      function iqc_rsa(iflag)

!     possible outputs
!     -1 bad data
!      0 unknown quality or unknown input flag
!     +1 good data

      iqc_rsa = 0

      if(iflag .eq.     0)iqc_rsa = +1 ! ok
      if(iflag .eq.     1)iqc_rsa = -1 ! out of range
      if(iflag .eq.     2)iqc_rsa = -1 ! questionable
      if(iflag .eq.     3)iqc_rsa =  0 ! not tested
      if(iflag .eq.     4)iqc_rsa = -1 ! missing data
      if(iflag .eq.     5)iqc_rsa = +1 ! mffg auto mode algorithm good data
      if(iflag .eq.     6)iqc_rsa = +1 ! mffg manual mode algorithm good data 
                                       ! operator no action 
      if(iflag .eq.     7)iqc_rsa = -1 ! mffg auto mode algorithm bad data
      if(iflag .eq.     8)iqc_rsa = -1 ! mffg auto mode algorithm bad data
                                       ! operator no action 
      if(iflag .eq.     9)iqc_rsa =  0 ! reserved for mffg but not used
      if(iflag .eq.    10)iqc_rsa = -1 ! mffg manual mode algorithm good data
                                       ! operator flagged as bad data
      if(iflag .eq.    11)iqc_rsa =  0 ! reserved for mffg but not used
      if(iflag .eq.    12)iqc_rsa = -1 ! mffg manual mode algorithm bad data 
                                       ! operator flagged as bad data
      if(iflag .eq.    13)iqc_rsa =  0 ! reserved for mffg but not used
      if(iflag .eq.    14)iqc_rsa = +1 ! mffg manual mode algorithm good data 
                                       ! operator no action
      if(iflag .eq.    15)iqc_rsa =  0 ! reserved for mffg but not used
      if(iflag .eq.    16)iqc_rsa = -1 ! mffg manual mode algorithm bad data 
                                       ! operator no action
      if(iflag .eq.    17)iqc_rsa =  0 ! reserved for mffg but not used
      if(iflag .eq.    18)iqc_rsa = -1 ! mffg manual mode algorithm good data 
                                       ! operator flagged as bad data
      if(iflag .eq.    19)iqc_rsa =  0 ! reserved for mffg but not used
      if(iflag .eq.    20)iqc_rsa = -1 ! mffg manual mode algorithm bad data 
                                       ! operator flagged as bad data
      if(iflag .eq. -9999)iqc_rsa = -1 ! missing data

      return
      end
