
      subroutine get_pirep_data(i4time_sys,ilaps_cycle_time,filename
     1                                                     ,ext
     1                                                     ,nx_l,ny_l
     1                                                     ,istatus)

      character*170 filename
      character*(*) ext

!.............................................................................

      include 'netcdf.inc'
      integer maxicinglvls, maxlocs, maxskylvls, maxturbelements,
     +     maxturblvls, recnum,nf_fid, nf_vid, nf_status
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
c get size of maxicinglvls
c
      nf_status = nf_inq_dimid(nf_fid,'maxicinglvls',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxicinglvls'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,maxicinglvls)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxicinglvls'
      endif
c
c get size of maxlocs
c
      nf_status = nf_inq_dimid(nf_fid,'maxlocs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxlocs'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,maxlocs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxlocs'
      endif
c
c get size of maxskylvls
c
      nf_status = nf_inq_dimid(nf_fid,'maxskylvls',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxskylvls'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,maxskylvls)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxskylvls'
      endif
c
c get size of maxturbelements
c
      nf_status = nf_inq_dimid(nf_fid,'maxturbelements',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxturbelements'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,maxturbelements)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxturbelements'
      endif
c
c get size of maxturblvls
c
      nf_status = nf_inq_dimid(nf_fid,'maxturblvls',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxturblvls'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,maxturblvls)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxturblvls'
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
      call pireps_sub(nf_fid, maxicinglvls, maxlocs, maxskylvls,
     +     maxturbelements, maxturblvls, recnum,
!.............................................................................
     1     ext,i4time_sys,ilaps_cycle_time,nx_l,ny_l,istatus)
      return
!.............................................................................
      end
c
c
      subroutine pireps_sub(nf_fid, maxicinglvls, maxlocs, maxskylvls,
     +     maxturbelements, maxturblvls, recnum,
!.............................................................................
     1     ext,i4time_sys,ilaps_cycle_time,nx_l,ny_l,istatus)
!.............................................................................

      include 'netcdf.inc'
      integer maxicinglvls, maxlocs, maxskylvls, maxturbelements,
     +     maxturblvls, recnum,nf_fid, nf_vid, nf_status
      integer icingintens( maxicinglvls, recnum),
     +     lowlvlwndshr(recnum), temperature(recnum), turbintens(
     +     maxturbelements,  maxturblvls, recnum), vis(recnum),
     +     winddir(recnum), windspd(recnum)
      real fltlvlbottom(recnum), fltlvltop(recnum), icingbottom(
     +     maxicinglvls, recnum), icingtop( maxicinglvls, recnum),
     +     lat( maxlocs, recnum), lon( maxlocs, recnum),
     +     skycvrbottom( maxskylvls, recnum), skycvrtop( maxskylvls,
     +     recnum), turbbottom( maxturblvls, recnum), turbtop(
     +     maxturblvls, recnum)
      double precision recpttime(recnum), timeobs(recnum)
      character*4 icinghtind( maxicinglvls, recnum)
      character*8 skycvramt( maxskylvls, recnum)
      character*4 turbhtind( maxturblvls, recnum)
      character*16 locstr( maxlocs, recnum)
      character*4 skycvrhtind( maxskylvls, recnum)
      character*4 turbfreq( maxturbelements,  maxturblvls, recnum)
      character*33 fltwthr(recnum)
      character*5 icingtype( maxicinglvls, recnum)
      character*201 origreport(recnum)
      character*5 aircrafttype(recnum)
      character*4 reporttype(recnum)
      character*4 fltlvlind(recnum)
      character*4 collsite(recnum)
      character*5 turbtype( maxturblvls, recnum)

!.............................................................................

      character*170 filename
      character*(*) ext
      character*9 a9_timeobs,a9_recpttime 
      character*7 c7_skycover
      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)

!.............................................................................

      call read_pireps_netcdf(nf_fid, maxicinglvls, maxlocs, maxskylvls,     
     +     maxturbelements, maxturblvls, recnum, icingintens, 
     +     lowlvlwndshr, temperature, turbintens, vis, winddir, 
     +     windspd, fltlvlbottom, fltlvltop, icingbottom, icingtop, 
     +     lat, lon, skycvrbottom, skycvrtop, turbbottom, turbtop, 
     +     recpttime, timeobs, aircrafttype, collsite, fltlvlind, 
     +     fltwthr, icinghtind, icingtype, locstr, origreport, 
     +     reporttype, skycvramt, skycvrhtind, turbfreq, turbhtind, 
     +     turbtype)
c
c the netcdf variables are filled - your code goes here
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

          if(lat(1,i) .ge. r_nc_missing_data)then
              write(6,*)' missing first latitude',i
              goto 999
          endif
          if(lon(1,i) .ge. r_nc_missing_data)then
              write(6,*)' missing first longitude',i
              goto 999
          endif
          if(timeobs(i) .ge. r_nc_missing_data)then
              write(6,*)' missing ob time',i
              goto 999
          endif
          if(recpttime(i) .ge. r_nc_missing_data)then
              write(6,*)' missing received time',i
              goto 999
          endif

!         test to see how many lat/lons are present
          n_latitude_present = 0
          do j = 1,4
            if(lat(j,i) .lt. r_nc_missing_data)then
              write(6,*)' lat/lon = ',lat(j,i),lon(j,i)
              n_latitude_present = n_latitude_present + 1
            endif
          enddo ! j

          write(6,*)' num locations = ',n_latitude_present       

          if(n_latitude_present .gt. 1)then
              write(6,*)' multiple locations, reject ob',i
              goto 999
          endif

          if(lat(1,i) .le. rnorth .and. lat(1,i) .ge. south .and.
     1       lon(1,i) .ge. west   .and. lon(1,i) .le. east      )then        
              write(6,*)' pirep is inside domain lat/lon perimeter'
          else
              write(6,*)
     1            ' outside domain lat/lon perimeter - reject'
              goto 999
          endif

!         write individual pirep to laps pin file
          call c_time2fname(nint(timeobs(i)),a9_timeobs)
          call c_time2fname(nint(recpttime(i)),a9_recpttime)

          call cv_asc_i4time(a9_timeobs,i4time_ob)
          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. (ilaps_cycle_time / 2) )then
              write(6,*)' outside time window - reject '
     1                              ,a9_timeobs,i4_resid
              goto 999        
          endif

          call open_ext(31,i4time_sys,ext(1:3),istatus)

          write(6,1)a9_timeobs,a9_recpttime 
          write(31,1)a9_timeobs,a9_recpttime 
 1        format(' time - prp/rcvd:'/1x,a9,2x,a9) 

          write(6,2)lat(1,i),lon(1,i),altitude
          write(31,2)lat(1,i),lon(1,i),altitude
 2        format(' lat, lon, altitude'/f8.3,f10.3,f8.0)  

          write(6,33)
          write(31,33)
 33       format(' cloud layer')

!         write out cloud base/top in feet and cloud amount in eighths
          do ilyr = 1,3

!             use -1000. for missing value of cloud base/top
              if(skycvrbottom(ilyr,i) .ge. 1e10)then
                  rbase = -1000.
              else
                  rbase = skycvrbottom(ilyr,i) * 3.281
              endif
              if(skycvrtop(ilyr,i) .ge. 1e10)then
                  rtop = -1000.
              else
                  rtop = skycvrtop(ilyr,i) * 3.281
              endif

              c7_skycover(1:7) = skycvramt(ilyr,i)(1:7)

              call skycover_to_frac(c7_skycover,fraction,istatus)

!             test for missing or unusable cloud cover
              if(istatus .ne. 1)then 
                  ieighths = -999
              else 
                  ieighths = nint(fraction * 8.0)
              endif                  

              write(6,3)rbase,rtop,ieighths,skycvramt(ilyr,i)
     1                                     ,fraction
 3            format(' cloud layer',2f8.0,i5,1x,a8,f5.2)

!             if(istatus .eq. 1)then
!                 write(6,*)' above layer written to pin file'
                  write(31,3)rbase,rtop,ieighths
!             endif

          enddo ! ilyr

 999      continue

      enddo ! i

!.............................................................................

      return
      end
c
c  subroutine to read the file "pirep data : selected by receipt time : time range from 885926400 to 885926700" 
c
      subroutine read_pireps_netcdf(nf_fid, maxicinglvls, maxlocs, 
     +     maxskylvls, maxturbelements, maxturblvls, recnum, 
     +     icingintens, lowlvlwndshr, temperature, turbintens, vis, 
     +     winddir, windspd, fltlvlbottom, fltlvltop, icingbottom, 
     +     icingtop, lat, lon, skycvrbottom, skycvrtop, turbbottom, 
     +     turbtop, recpttime, timeobs, aircrafttype, collsite, 
     +     fltlvlind, fltwthr, icinghtind, icingtype, locstr, 
     +     origreport, reporttype, skycvramt, skycvrhtind, turbfreq, 
     +     turbhtind, turbtype)
c
      include 'netcdf.inc'
      integer maxicinglvls, maxlocs, maxskylvls, maxturbelements, 
     +     maxturblvls, recnum,nf_fid, nf_vid, nf_status
      integer icingintens( maxicinglvls, recnum),
     +     lowlvlwndshr(recnum), temperature(recnum), turbintens(
     +     maxturbelements,  maxturblvls, recnum), vis(recnum),
     +     winddir(recnum), windspd(recnum)
      real fltlvlbottom(recnum), fltlvltop(recnum), icingbottom(
     +     maxicinglvls, recnum), icingtop( maxicinglvls, recnum),
     +     lat( maxlocs, recnum), lon( maxlocs, recnum),
     +     skycvrbottom( maxskylvls, recnum), skycvrtop( maxskylvls,
     +     recnum), turbbottom( maxturblvls, recnum), turbtop(
     +     maxturblvls, recnum)
      double precision recpttime(recnum), timeobs(recnum)
      character*4 icinghtind( maxicinglvls, recnum)
      character*8 skycvramt( maxskylvls, recnum)
      character*4 turbhtind( maxturblvls, recnum)
      character*16 locstr( maxlocs, recnum)
      character*4 skycvrhtind( maxskylvls, recnum)
      character*33 fltwthr(recnum)
      character*4 turbfreq( maxturbelements,  maxturblvls, recnum)
      character*5 icingtype( maxicinglvls, recnum)
      character*201 origreport(recnum)
      character*5 aircrafttype(recnum)
      character*4 reporttype(recnum)
      character*4 fltlvlind(recnum)
      character*4 collsite(recnum)
      character*5 turbtype( maxturblvls, recnum)


c   variables of type real
c
c     variable        netcdf long name
c      fltlvlbottom "bottom of pireps flight level"
c
        nf_status = nf_inq_varid(nf_fid,'fltlvlbottom',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var fltlvlbottom'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,fltlvlbottom)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var fltlvlbottom'
      endif
c
c     variable        netcdf long name
c      fltlvltop    "top of pireps flight level"
c
        nf_status = nf_inq_varid(nf_fid,'fltlvltop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var fltlvltop'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,fltlvltop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var fltlvltop'
      endif
c
c     variable        netcdf long name
c      icingbottom  "bottom of icing layer"
c
        nf_status = nf_inq_varid(nf_fid,'icingbottom',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var icingbottom'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,icingbottom)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var icingbottom'
      endif
c
c     variable        netcdf long name
c      icingtop     "top of icing layer"
c
        nf_status = nf_inq_varid(nf_fid,'icingtop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var icingtop'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,icingtop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var icingtop'
      endif
c
c     variable        netcdf long name
c      lat          "latitude of report"
c
        nf_status = nf_inq_varid(nf_fid,'lat',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lat'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,lat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lat'
      endif
c
c     variable        netcdf long name
c      lon          "longitude of report"
c
        nf_status = nf_inq_varid(nf_fid,'lon',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lon'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,lon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lon'
      endif
c
c     variable        netcdf long name
c      skycvrbottom "bottom of sky cover layer"
c
        nf_status = nf_inq_varid(nf_fid,'skycvrbottom',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var skycvrbottom'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,skycvrbottom)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var skycvrbottom'
      endif
c
c     variable        netcdf long name
c      skycvrtop    "bottom of sky cover layer"
c
        nf_status = nf_inq_varid(nf_fid,'skycvrtop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var skycvrtop'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,skycvrtop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var skycvrtop'
      endif
c
c     variable        netcdf long name
c      turbbottom   "bottom of turbulence layer"
c
        nf_status = nf_inq_varid(nf_fid,'turbbottom',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbbottom'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,turbbottom)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbbottom'
      endif
c
c     variable        netcdf long name
c      turbtop      "top of turbulence layer"
c
        nf_status = nf_inq_varid(nf_fid,'turbtop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbtop'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,turbtop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbtop'
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      icingintens  "icing intensity"
c
        nf_status = nf_inq_varid(nf_fid,'icingintens',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var icingintens'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,icingintens)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var icingintens'
      endif
c
c     variable        netcdf long name
c      lowlvlwndshr "low level wind shear indicator"
c
        nf_status = nf_inq_varid(nf_fid,'lowlvlwndshr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lowlvlwndshr'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,lowlvlwndshr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lowlvlwndshr'
      endif
c
c     variable        netcdf long name
c      temperature  "temperature"
c
        nf_status = nf_inq_varid(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,temperature)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
c
c     variable        netcdf long name
c      turbintens   "turbulence intensity"
c
        nf_status = nf_inq_varid(nf_fid,'turbintens',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbintens'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,turbintens)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbintens'
      endif
c
c     variable        netcdf long name
c      vis          "visibility"
c
        nf_status = nf_inq_varid(nf_fid,'vis',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vis'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,vis)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vis'
      endif
c
c     variable        netcdf long name
c      winddir      "wind direction"
c
        nf_status = nf_inq_varid(nf_fid,'winddir',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,winddir)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      endif
c
c     variable        netcdf long name
c      windspd      "wind speed"
c
        nf_status = nf_inq_varid(nf_fid,'windspd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspd'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,windspd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspd'
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c      recpttime    "time report was received"
c
        nf_status = nf_inq_varid(nf_fid,'recpttime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var recpttime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,recpttime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var recpttime'
      endif
c
c     variable        netcdf long name
c      timeobs      "time of report"
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


c   variables of type char
c
c
c     variable        netcdf long name
c      aircrafttype "type of aircraft issuing report"
c
        nf_status = nf_inq_varid(nf_fid,'aircrafttype',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var aircrafttype'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,aircrafttype)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var aircrafttype'
      endif
c
c     variable        netcdf long name
c      collsite     "data collection site"
c
        nf_status = nf_inq_varid(nf_fid,'collsite',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var collsite'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,collsite)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var collsite'
      endif
c
c     variable        netcdf long name
c      fltlvlind    "flight level indicator"
c
        nf_status = nf_inq_varid(nf_fid,'fltlvlind',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var fltlvlind'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,fltlvlind)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var fltlvlind'
      endif
c
c     variable        netcdf long name
c      fltwthr      "flight weather"
c
        nf_status = nf_inq_varid(nf_fid,'fltwthr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var fltwthr'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,fltwthr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var fltwthr'
      endif
c
c     variable        netcdf long name
c      icinghtind   "icing height indicator"
c
        nf_status = nf_inq_varid(nf_fid,'icinghtind',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var icinghtind'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,icinghtind)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var icinghtind'
      endif
c
c     variable        netcdf long name
c      icingtype    "type of icing observed"
c
        nf_status = nf_inq_varid(nf_fid,'icingtype',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var icingtype'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,icingtype)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var icingtype'
      endif
c
c     variable        netcdf long name
c      locstr       "location string for each derived lat/lon"
c
        nf_status = nf_inq_varid(nf_fid,'locstr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var locstr'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,locstr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var locstr'
      endif
c
c     variable        netcdf long name
c      origreport   "original pireps report"
c
        nf_status = nf_inq_varid(nf_fid,'origreport',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var origreport'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,origreport)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var origreport'
      endif
c
c     variable        netcdf long name
c      reporttype   "report type"
c
        nf_status = nf_inq_varid(nf_fid,'reporttype',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reporttype'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,reporttype)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reporttype'
      endif
c
c     variable        netcdf long name
c      skycvramt    "amount of sky cover"
c
        nf_status = nf_inq_varid(nf_fid,'skycvramt',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var skycvramt'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,skycvramt)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var skycvramt'
      endif
c
c     variable        netcdf long name
c      skycvrhtind  "sky cover height indicator"
c
        nf_status = nf_inq_varid(nf_fid,'skycvrhtind',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var skycvrhtind'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,skycvrhtind)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var skycvrhtind'
      endif
c
c     variable        netcdf long name
c      turbfreq     "weather qualifier"
c
        nf_status = nf_inq_varid(nf_fid,'turbfreq',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbfreq'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,turbfreq)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbfreq'
      endif
c
c     variable        netcdf long name
c      turbhtind    "turbulence height indicator"
c
        nf_status = nf_inq_varid(nf_fid,'turbhtind',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbhtind'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,turbhtind)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbhtind'
      endif
c
c     variable        netcdf long name
c      turbtype     "type of turbulence observed"
c
        nf_status = nf_inq_varid(nf_fid,'turbtype',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbtype'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,turbtype)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var turbtype'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
