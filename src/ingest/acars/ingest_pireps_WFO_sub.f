
      subroutine get_pirep_data_wfo(i4time_sys,ilaps_cycle_time,
     1                              filename,ext,nx_l,ny_l,istatus)

      character*(*) filename
      character*(*) ext

!.............................................................................

      include 'netcdf.inc'
      integer maxnumcloudlayers, recnum,nf_fid, nf_vid, nf_status
c
c  open netcdf file for reading
c
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',filename
      endif
c
c  fill all dimension values
c
c
c get size of maxnumcloudlayers
c
      nf_status = nf_inq_dimid(nf_fid,'maxnumcloudlayers',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxnumcloudlayers'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,maxnumcloudlayers)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxnumcloudlayers'
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

      call pireps_wfo_sub(nf_fid, maxnumcloudlayers, recnum,
!.............................................................................
     1     ext,i4time_sys,ilaps_cycle_time,nx_l,ny_l,istatus)
      return
!.............................................................................
      end
c
c
      subroutine pireps_wfo_sub(nf_fid, maxnumcloudlayers, recnum,
!.............................................................................
     1     ext,i4time_sys,ilaps_cycle_time,nx_l,ny_l,istatus)
!.............................................................................

      include 'netcdf.inc'
      integer maxnumcloudlayers, recnum,nf_fid, nf_vid, 
     +        nf_status, cloudamt(maxnumcloudlayers,recnum),
     +        written
      real lat(recnum), lon(recnum), 
     +     cloudbaseht(maxnumcloudlayers,recnum),
     +     cloudtopht(maxnumcloudlayers,recnum)
      double precision timeobs(recnum)

!.............................................................................

      integer i4time_sys, ilaps_cycle_time, nx_l, ny_l, istatus
      character*(*) ext
      character*9 a9_timeobs 
      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)

!.............................................................................
c
c      read latitude          "latitude of report"
c
        nf_status = nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,lat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
c
c      read longitude          "longitude of report"
c
        nf_status = nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,lon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
c
c      read timeobs      "time of report"
c
        nf_status = nf_inq_varid(nf_fid,'timeobs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timeobs'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,timeobs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timeobs'
      endif
c
c      read cloudbaseheight
c
        nf_status = nf_inq_varid(nf_fid,'cloudbaseheight',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbaseheight'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,cloudbaseht)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudbaseheight'
      endif
c
c      read cloudtopheight
c
        nf_status = nf_inq_varid(nf_fid,'cloudtopheight',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudtopheight'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,cloudtopht)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudtopheight'
      endif
c
c      read cloudamount
c
        nf_status = nf_inq_varid(nf_fid,'cloudamount',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudamount'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,cloudamt)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var cloudamount'
      endif
c
c the netcdf variables are filled 
c
!.............................................................................

      call get_domain_perimeter(nx_l,ny_l,'nest7grid',lat_a,lon_a,
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error in get_domain_perimeter'
          return
      endif

!     write all pireps to laps pin file
      altitude = 0.

      r_nc_missing_data = 1e20

      write(6,*)' recnum = ',recnum

      num_pireps = recnum

      do i = 1,num_pireps

          write(6,*)
          write(6,*)' pirep #',i

          if(lat(i) .ge. r_nc_missing_data)then
              write(6,*)' missing first latitude',i
              goto 999
          endif
          if(lon(i) .ge. r_nc_missing_data)then
              write(6,*)' missing first longitude',i
              goto 999
          endif
          if(timeobs(i) .ge. r_nc_missing_data)then
              write(6,*)' missing ob time',i
              goto 999
          endif

!         test to see how many lat/lons are present
          n_latitude_present = 0
          if(lat(i) .lt. r_nc_missing_data)then
            write(6,*)' lat/lon = ',lat(i),lon(i)
            n_latitude_present = n_latitude_present + 1
          endif

          write(6,*)' num locations = ',n_latitude_present       

          if(n_latitude_present .gt. 1)then
              write(6,*)' multiple locations, reject ob',i
              goto 999
          endif

          if(lat(i) .le. rnorth .and. lat(i) .ge. south .and.
     1       lon(i) .ge. west   .and. lon(i) .le. east      )then        
              write(6,*)' pirep is inside domain lat/lon perimeter'
          else
              write(6,*)
     1            ' outside domain lat/lon perimeter - reject'
              goto 999
          endif

!         write individual pirep to laps pin file
          call c_time2fname(nint(timeobs(i)),a9_timeobs)

          call cv_asc_i4time(a9_timeobs,i4time_ob)
          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. (ilaps_cycle_time / 2) )then
              write(6,*)' outside time window - reject '
     1                              ,a9_timeobs,i4_resid
              goto 999        
          endif


!         write out cloud base/top in feet and cloud amount in eighths
!         written = 0 if timeobs/lat/lon not written out, set to 1 once written to file
          written = 0
          do ilyr = 1,maxnumcloudlayers

!             use -1000. for missing value of cloud base/top
              if ((cloudbaseht(ilyr,i) .ge. 1e10) .or.
     +         (cloudbaseht(ilyr,i) .eq. -9999.)) then
                  rbase = -1000.
              else
                  rbase = cloudbaseht(ilyr,i)
              endif

              if((cloudtopht(ilyr,i) .ge. 1e10) .or.
     +         (cloudtopht(ilyr,i) .eq. -9999.)) then
                  rtop = -1000.
              else
                  rtop = cloudtopht(ilyr,i)
              endif

!             test for missing or unusable cloud cover
!             wfo bases cloudamount on wfm bufr table 020011
              if((cloudamt(ilyr,i) .gt. 8) .or.
     +           (cloudamt(ilyr,i) .lt. 0)) then 
                  ieighths = -999
              else 
                  ieighths = cloudamt(ilyr,i)
              endif                  

              write(6,3)rbase,rtop,ieighths
 3            format(' cloud layer',2f8.0,i5)

              if ((rbase .ne. -1000.) .and. (rtop .ne. -1000.) .and.
     +            (ieighths .ne. -999)) then
                if (written .eq. 0) then   ! write timeobs, lat, lon, etc
                  call open_ext(31,i4time_sys,ext(1:3),istatus)

                  write(6,1)a9_timeobs 
                  write(31,1)a9_timeobs, a9_timeobs 
 1                format(' time - prp/rcvd:'/1x,a9,2x,a9) 

                  write(6,2)lat(i),lon(i),altitude
                  write(31,2)lat(i),lon(i),altitude
 2                format(' lat, lon, altitude'/f8.3,f10.3,f8.0)  

                  write(6,33)
                  write(31,33)
 33               format(' cloud layer')
                  written = 1
                endif
!               write(6,*)' above layer written to pin file'
                write(31,3)rbase,rtop,ieighths
              endif

          enddo ! ilyr

          if (written .eq. 0) then ! none of layers in pirep contained info
            write(6,*)' above pirep not written to pin file - ',
     +'no usable data.'
          endif

 999      continue

      enddo ! i

!.............................................................................

      return
      end
