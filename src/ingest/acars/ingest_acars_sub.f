
      subroutine get_acars_data(i4time_sys,i4_acars_window
     1                                    ,nx_l,ny_l
     1                                    ,c8_project
     1                                    ,ext
     1                                    ,l_use_tamdar
     1                                    ,filename,istatus)

      character*(*) filename,ext
      logical l_use_tamdar

!.............................................................................

      include 'netcdf.inc'
      integer recnum,nf_fid, nf_vid, nf_status
      character*8 c8_project

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
      call acars_sub(nf_fid, recnum, c8_project, ext, l_use_tamdar,
!.............................................................................
     1              i4time_sys,i4_acars_window,nx_l,ny_l,istatus)
      return
!.............................................................................
      end
c
c
      subroutine acars_sub(nf_fid, recnum, c8_project, ext, 
     1                     l_use_tamdar,
!.............................................................................
     1              i4time_sys,i4_acars_window,nx_l,ny_l,istatus)
!.............................................................................

      include 'netcdf.inc'
      character*8 c8_project
      character*(*) ext
      integer recnum,nf_fid, nf_vid, nf_status
      integer airline(recnum), bounceerror(recnum),
     +     correctedflag(recnum), datadescriptor(recnum),
     +     datasource(recnum),
     +     errortype(recnum), interpolatedll(recnum),
     +     interpolatedtime(recnum), missinginputminutes,
     +     rollflag(recnum), speederror(recnum), temperror(recnum),
     +     watervaporqc(recnum), winddirerror(recnum),
     +     windspeederror(recnum)
      real altitude(recnum), heading(recnum), latitude(recnum),
     +     longitude(recnum), maxturbulence(recnum),
     +     medturbulence(recnum), temperature(recnum),
     +     vertaccel(recnum), downlinkedrh(recnum), winddir(recnum),
     +     windspeed(recnum)
      double precision maxsecs, minsecs, timeobs(recnum),
     +     timereceived(recnum)
      character*4 rptstation(recnum)
      character*6 destairport(recnum)
      character*30 mindate
      character*6 origairport(recnum)
      character*9 tailnumber(recnum)
      character*30 maxdate
      character*13 flight(recnum)

!.............................................................................

      character*9 a9_timeobs,a9_recpttime 
      character*7 c7_skycover
      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)
      logical l_use_tamdar, l_debug

!............................................................................

      if (c8_project(1:3) .eq. 'wfo' .or. 
     1    c8_project(1:3) .eq. 'rsa'      ) then     
        call read_acars_netcdf_wfo(nf_fid, recnum, airline, bounceerror, 
     +     correctedflag, datadescriptor, errortype, interpolatedll, 
     +     interpolatedtime, missinginputminutes, rollflag, 
     +     speederror, temperror, watervaporqc, winddirerror, 
     +     windspeederror, altitude, heading, latitude, longitude, 
     +     maxturbulence, medturbulence, temperature, vertaccel, 
     +     downlinkedrh, winddir, windspeed, maxsecs, minsecs, 
     +     timeobs, timereceived, destairport, flight, maxdate, 
     +     mindate, origairport, rptstation, tailnumber)

         datasource = 0

      else
        call read_acars_netcdf(nf_fid, recnum, airline, bounceerror, 
     +     correctedflag, datadescriptor, datasource,
     +     errortype, interpolatedll, 
     +     interpolatedtime, missinginputminutes, rollflag, 
     +     speederror, temperror, watervaporqc, winddirerror, 
     +     windspeederror, altitude, heading, latitude, longitude, 
     +     maxturbulence, medturbulence, temperature, vertaccel, 
     +     downlinkedrh, winddir, windspeed, maxsecs, minsecs, 
     +     timeobs, timereceived, destairport, flight, maxdate, 
     +     mindate, origairport, rptstation, tailnumber)

      endif
c
c the netcdf variables are filled - your code goes here
c
!............................................................................

      call get_domain_perimeter(nx_l,ny_l,'nest7grid',lat_a,lon_a, 
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error in get_laps_perimeter'
          return
      endif

      num_acars = recnum
      write(6,*)' # of acars = ',recnum

      do i = 1,num_acars

          if(i .le. 1000 .or. i .eq. ((i/10)*10) )then
              l_debug = .true.
          else
              l_debug = .false.
          endif

          if(l_debug)write(6,*)
          if(l_debug)write(6,*)' acars #',i,'  ',char(datadescriptor(i))       
     1                                          ,char(errortype(i))
     1                                          ,datasource(i)
!         write(6,*)' location = '
!    1             ,latitude(i),longitude(i),altitude(i)

          if((.not. l_use_tamdar) .and. datasource(i) .eq. 4)then
              if(l_debug)write(6,*)' tamdar observation - reject'
              goto 900
          endif

          if(latitude(i) .le. rnorth .and. latitude(i)  .ge. south .and.
     1       longitude(i) .ge. west  .and. longitude(i) .le. east      
     1                                                             )then       
              continue
          else ! outside lat/lon perimeter - reject
              if(l_debug)write(6,*)' lat/lon - reject'       
!    1                 ,latitude(i),longitude(i)
              goto 900
          endif

          if(nanf(altitude(i)) .eq. 1)then
              if(l_debug)write(6,*)' altitude failed nan test - reject'       
     1                            ,altitude(i)
              goto 900
          endif

          if(altitude(i) .gt. 20000. .or. altitude(i) .lt. -1000.
     1                               .or. altitude(i) .eq.     0.)then
              if(l_debug)write(6,*)' altitude is suspect - reject'
     1                            ,altitude(i)
              goto 900
          endif

          if(abs(timeobs(i))      .lt. 3d9       .and.
     1       abs(timereceived(i)) .lt. 3d9              )then
              call c_time2fname(nint(timeobs(i)),a9_timeobs)
              call c_time2fname(nint(timereceived(i)),a9_recpttime)
          else
              if(l_debug)write(6,*)' bad observation time - reject'       
     1                   ,timeobs(i),timereceived(i)
              goto 900
          endif


          call cv_asc_i4time(a9_timeobs,i4time_ob)
          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. i4_acars_window)then ! outside time window
              if(l_debug)write(6,*)' time - reject '
     1           ,a9_timeobs,i4_resid,i4_acars_window
              goto 900        
          endif

          lun = 31

          call open_ext(lun,i4time_sys,ext(1:3),istatus)       

          if(l_debug)write(6,1)a9_timeobs,a9_recpttime 
                     write(lun,1)a9_timeobs,a9_recpttime 
 1        format(' time - prp/rcvd:'/1x,a9,2x,a9) 

!         l_geoalt is implicitly false with pressure altitude data
          if(l_debug)write(6,2)latitude(i),longitude(i),altitude(i)
          write(lun,2)          latitude(i),longitude(i),altitude(i)
 2        format(' lat, lon, altitude'/f8.3,f10.3,f8.0)  

!         test for bad winds
          if(char(datadescriptor(i)) .eq. 'x')then
            if(char(errortype(i)) .eq. 'w' .or. 
     1         char(errortype(i)) .eq. 'b'                         )then
              if(l_debug)write(6,*)' qc flag is bad - reject wind'
     1                 ,char(datadescriptor(i)),char(errortype(i))
              goto 850
            endif
          endif

          if(abs(windspeed(i)) .gt. 250.)then
              if(l_debug)write(6,*)' wind speed is suspect - reject'
     1                              ,windspeed(i)

          elseif(int(winddir(i)).lt.0 .or. int(winddir(i)).gt.360)then     
              if(l_debug)write(6,*)' wind direction is suspect - reject'       
     1                              ,winddir(i)

          else ! write out valid wind
              if(l_debug)write(6 ,3)int(winddir(i)),windspeed(i)
                         write(lun,3)int(winddir(i)),windspeed(i)
 3            format(' wind:'/' ', i3, ' deg @ ', f6.1, ' m/s')

          endif

 850      continue

!         test for bad temps
          if(char(datadescriptor(i)) .eq. 'x')then
            if(char(errortype(i)) .eq. 't' .or. 
     1         char(errortype(i)) .eq. 'b'                         )then
              if(l_debug)write(6,*)' qc flag is bad - reject temp'
     1                 ,char(datadescriptor(i)),char(errortype(i))
              goto 860
            endif
          endif

          if(abs(temperature(i)) .lt. 400.)then
              if(l_debug)write(6,13)temperature(i)
                         write(lun,13)temperature(i)
 13           format(' temp:'/1x,f10.1)
       
          else
              if(l_debug)write(6,*)' temperature is suspect - reject'
     1                             , temperature(i)

          endif

 860      continue

          if(downlinkedrh(i) .ge. 0.   .and. 
     1       downlinkedrh(i) .le. 1.00 .and.
     1       watervaporqc(i) .le. 2    .and.
     1       watervaporqc(i) .ge. 0             )then
              if(l_debug)write(6,23)downlinkedrh(i)
              if(l_debug)write(6,*)' rh qc value = ',watervaporqc(i)
              write(lun,23)downlinkedrh(i)
 23           format(' rh:'/1x,f10.3)

          else
              if(l_debug)write(6,*)' rh rejected: '
     1                             ,downlinkedrh(i),watervaporqc(i)

          endif

 900  enddo ! i

!............................................................................

      return
      end
c
c  subroutine to read the file "acars data" 
c
      subroutine read_acars_netcdf(nf_fid, recnum, airline, bounceerror, 
     +     correctedflag, datadescriptor, datasource,
     +     errortype, interpolatedll, 
     +     interpolatedtime, missinginputminutes, rollflag, 
     +     speederror, temperror, watervaporqc, winddirerror, 
     +     windspeederror, altitude, heading, latitude, longitude, 
     +     maxturbulence, medturbulence, temperature, vertaccel, 
     +     downlinkedrh, winddir, windspeed, maxsecs, minsecs, 
     +     timeobs, timereceived, destairport, flight, maxdate, 
     +     mindate, origairport, rptstation, tailnumber)
c
      include 'netcdf.inc'
      integer recnum,nf_fid, nf_vid, nf_status
      integer airline(recnum), bounceerror(recnum),
     +     correctedflag(recnum), datadescriptor(recnum),
     +     datasource(recnum),
     +     errortype(recnum), interpolatedll(recnum),
     +     interpolatedtime(recnum), missinginputminutes,
     +     rollflag(recnum), speederror(recnum), temperror(recnum),
     +     watervaporqc(recnum), winddirerror(recnum),
     +     windspeederror(recnum)
      real altitude(recnum), heading(recnum), latitude(recnum),
     +     longitude(recnum), maxturbulence(recnum),
     +     medturbulence(recnum), temperature(recnum),
     +     vertaccel(recnum), downlinkedrh(recnum), winddir(recnum),
     +     windspeed(recnum)
      double precision maxsecs, minsecs, timeobs(recnum),
     +     timereceived(recnum)
      character*4 rptstation(recnum)
      character*6 destairport(recnum)
      character*30 mindate
      character*6 origairport(recnum)
      character*9 tailnumber(recnum)
      character*30 maxdate
      character*13 flight(recnum)


c   variables of type real
c
c     variable        netcdf long name
c      altitude     
c
        nf_status = nf_inq_varid(nf_fid,'altitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var altitude'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,altitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var altitude'
      endif
c
c     variable        netcdf long name
c      heading      "heading of flight path over ground"
c
        nf_status = nf_inq_varid(nf_fid,'heading',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var heading'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,heading)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var heading'
      endif
c
c     variable        netcdf long name
c      latitude     
c
        nf_status = nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,latitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
c
c     variable        netcdf long name
c      longitude    
c
        nf_status = nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,longitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
c
c     variable        netcdf long name
c      maxturbulence"maximum eddy dissipation rate"
c
        nf_status = nf_inq_varid(nf_fid,'maxturbulence',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var maxturbulence'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,maxturbulence)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var maxturbulence'
      endif
c
c     variable        netcdf long name
c      medturbulence"median eddy dissipation rate"
c
        nf_status = nf_inq_varid(nf_fid,'medturbulence',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var medturbulence'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,medturbulence)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var medturbulence'
      endif
c
c     variable        netcdf long name
c      temperature  
c
        nf_status = nf_inq_varid(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,temperature)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
c
c     variable        netcdf long name
c      vertaccel    "peak vertical acceleration"
c
        nf_status = nf_inq_varid(nf_fid,'vertaccel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vertaccel'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,vertaccel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vertaccel'
      endif
c
c     variable        netcdf long name
c      downlinkedrh "downlinked relative humidity"
c
        nf_status = nf_inq_varid(nf_fid,'downlinkedrh',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var downlinkedrh'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,downlinkedrh)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var downlinkedrh'
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
        nf_status = nf_get_var_real(nf_fid,nf_vid,winddir)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      endif
c
c     variable        netcdf long name
c      windspeed    "wind speed"
c
        nf_status = nf_inq_varid(nf_fid,'windspeed',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeed'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,windspeed)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeed'
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      airline      "airline"
c
        nf_status = nf_inq_varid(nf_fid,'airline',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var airline'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,airline)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var airline'
      endif
c
c     variable        netcdf long name
c      bounceerror  "aircraft altitude variance error"
c
        nf_status = nf_inq_varid(nf_fid,'bounceerror',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bounceerror'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,bounceerror)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var bounceerror'
      endif
c
c     variable        netcdf long name
c      correctedflag"corrected data indicator"
c
        nf_status = nf_inq_varid(nf_fid,'correctedflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var correctedflag'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,correctedflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var correctedflag'
      endif
c
c     variable        netcdf long name
c      datadescriptor"awips-type data descriptor"
c
        nf_status = nf_inq_varid(nf_fid,'datadescriptor',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var datadescriptor'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,datadescriptor)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var datadescriptor'
      endif
c
c     variable        netcdf long name
c      datasource"awips-type data source"
c
        nf_status = nf_inq_varid(nf_fid,'datasource',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var datasource'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,datasource)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var datasource'
      endif
c
c     variable        netcdf long name
c      errortype    
c
        nf_status = nf_inq_varid(nf_fid,'errortype',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var errortype'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,errortype)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var errortype'
      endif
c
c     variable        netcdf long name
c      interpolatedll"ups ascent/descent lat&lon interpolation indicator"
c
        nf_status = nf_inq_varid(nf_fid,'interpolatedll',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var interpolatedll'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,interpolatedll)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var interpolatedll'
      endif
c
c     variable        netcdf long name
c      interpolatedtime"ups ascent/descent time interpolation indicator"
c
        nf_status = nf_inq_varid(nf_fid,'interpolatedtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var interpolatedtime'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,interpolatedtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var interpolatedtime'
      endif
c
c     variable        netcdf long name
c      missinginputminutes"missing minutes of input data"
c
        nf_status = nf_inq_varid(nf_fid,'missinginputminutes',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var missinginputminutes'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,missinginputminutes)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var missinginputminutes'
      endif
c
c     variable        netcdf long name
c      rollflag     "aircraft roll angle flag "
c
        nf_status = nf_inq_varid(nf_fid,'rollflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rollflag'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,rollflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rollflag'
      endif
c
c     variable        netcdf long name
c      speederror   "aircraft ground speed error"
c
        nf_status = nf_inq_varid(nf_fid,'speederror',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var speederror'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,speederror)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var speederror'
      endif
c
c     variable        netcdf long name
c      temperror    
c
        nf_status = nf_inq_varid(nf_fid,'temperror',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperror'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,temperror)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperror'
      endif
c
c     variable        netcdf long name
c      watervaporqc "water vapor mixing ratio qc character"
c
        nf_status = nf_inq_varid(nf_fid,'watervaporqc',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var watervaporqc'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,watervaporqc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var watervaporqc'
      endif
c
c     variable        netcdf long name
c      winddirerror 
c
        nf_status = nf_inq_varid(nf_fid,'winddirerror',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddirerror'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,winddirerror)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddirerror'
      endif
c
c     variable        netcdf long name
c      windspeederror
c
        nf_status = nf_inq_varid(nf_fid,'windspeederror',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeederror'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,windspeederror)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeederror'
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c      maxsecs      "maximum observation time"
c
        nf_status = nf_inq_varid(nf_fid,'maxsecs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var maxsecs'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,maxsecs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var maxsecs'
      endif
c
c     variable        netcdf long name
c      minsecs      "minimum observation time"
c
        nf_status = nf_inq_varid(nf_fid,'minsecs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var minsecs'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,minsecs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var minsecs'
      endif
c
c     variable        netcdf long name
c      timeobs      "time of observation"
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
c     variable        netcdf long name
c      timereceived "time data was received at ground station"
c
        nf_status = nf_inq_varid(nf_fid,'timereceived',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timereceived'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,timereceived)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timereceived'
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c      destairport  "destination airport"
c
        nf_status = nf_inq_varid(nf_fid,'destairport',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var destairport'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,destairport)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var destairport'
      endif
c
c     variable        netcdf long name
c      flight       "flight number"
c
        nf_status = nf_inq_varid(nf_fid,'flight',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var flight'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,flight)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var flight'
      endif
c
c     variable        netcdf long name
c      maxdate      "maximum observation date"
c
        nf_status = nf_inq_varid(nf_fid,'maxdate',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var maxdate'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,maxdate)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var maxdate'
      endif
c
c     variable        netcdf long name
c      mindate      "minimum observation date"
c
        nf_status = nf_inq_varid(nf_fid,'mindate',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var mindate'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,mindate)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var mindate'
      endif
c
c     variable        netcdf long name
c      origairport  "originating airport"
c
        nf_status = nf_inq_varid(nf_fid,'origairport',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var origairport'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,origairport)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var origairport'
      endif
c
c     variable        netcdf long name
c      rptstation   "station reporting through"
c
        nf_status = nf_inq_varid(nf_fid,'rptstation',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rptstation'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,rptstation)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rptstation'
      endif
c
c     variable        netcdf long name
c      tailnumber   "tail number"
c
        nf_status = nf_inq_varid(nf_fid,'tailnumber',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tailnumber'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,tailnumber)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tailnumber'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
c
c  subroutine to read wfo format file "acars data" 
c
      subroutine read_acars_netcdf_wfo(nf_fid, recnum, airline, 
     +     bounceerror, 
     +     correctedflag, datadescriptor, errortype, interpolatedll, 
     +     interpolatedtime, missinginputminutes, rollflag, 
     +     speederror, temperror, watervaporqc, winddirerror, 
     +     windspeederror, altitude, heading, latitude, longitude, 
     +     maxturbulence, medturbulence, temperature, vertaccel, 
     +     downlinkedrh, winddir, windspeed, maxsecs, minsecs, 
     +     timeobs, timereceived, destairport, flight, maxdate, 
     +     mindate, origairport, rptstation, tailnumber)
c
      include 'netcdf.inc'
      integer recnum,nf_fid, nf_vid, nf_status
      integer airline(recnum), bounceerror(recnum),
     +     correctedflag(recnum), datadescriptor(recnum),
     +     errortype(recnum), interpolatedll(recnum),
     +     interpolatedtime(recnum), missinginputminutes,
     +     rollflag(recnum), speederror(recnum), temperror(recnum),
     +     watervaporqc(recnum), winddirerror(recnum),
     +     windspeederror(recnum)
      real altitude(recnum), heading(recnum), latitude(recnum),
     +     longitude(recnum), maxturbulence(recnum),
     +     medturbulence(recnum), temperature(recnum),
     +     vertaccel(recnum), downlinkedrh(recnum), winddir(recnum),
     +     windspeed(recnum)
      double precision maxsecs, minsecs, timeobs(recnum),
     +     timereceived(recnum)
      character*4 rptstation(recnum)
      character*6 destairport(recnum)
      character*30 mindate
      character*6 origairport(recnum)
      character*9 tailnumber(recnum)
      character*30 maxdate
      character*13 flight(recnum)
      integer i, wind_err, temp_err
      real ws, temp, alt, max_temp, min_temp, wmax

c   variables of type real
c
c     variable        netcdf long name
c      altitude     
c
        nf_status = nf_inq_varid(nf_fid,'indaltitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var indaltitude'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,altitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var altitude'
      endif
c
c     variable        netcdf long name
c      heading      "heading of flight path over ground"
c      heading not available in wfo file....fill with 99999.0
c
      do i = 1, recnum
        heading(i) = 99999.0
      enddo
c
c     variable        netcdf long name
c      latitude     
c
        nf_status = nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,latitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latitude'
      endif
c
c     variable        netcdf long name
c      longitude    
c
        nf_status = nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,longitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var longitude'
      endif
c
c     variable        netcdf long name
c      maxturbulence"maximum eddy dissipation rate"
c      maxturbulence not available in wfo file....fill with -9.99 
c
      do i = 1, recnum
        maxturbulence(i) = -9.99
      enddo
c
c     variable        netcdf long name
c      medturbulence"median eddy dissipation rate"
c      medturbulence not available in wfo file....fill with -9.99 
c
      do i = 1, recnum
        medturbulence(i) = -9.99
      enddo
c
c     variable        netcdf long name
c      temperature  
c
        nf_status = nf_inq_varid(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,temperature)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
c
c     variable        netcdf long name
c      vertaccel    "peak vertical acceleration"
c      vertaccel not available in wfo file....fill with 0.0  
c
      do i = 1, recnum
        vertaccel(i) = 0.0 
      enddo
c
c     variable        netcdf long name
c      downlinkedrh "downlinked relative humidity"
c
        nf_status = nf_inq_varid(nf_fid,'relhumidity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidity'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,downlinkedrh)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relhumidity'
      endif
c
c reset values in downlinkedrh to 0-1 range from percent
c
      do i = 1, recnum
        downlinkedrh(i) = downlinkedrh(i)/100.0
      enddo
c
c     variable        netcdf long name
c      winddir      "wind direction"
c
        nf_status = nf_inq_varid(nf_fid,'winddir',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,winddir)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      endif
c
c     variable        netcdf long name
c      windspeed    "wind speed"
c
        nf_status = nf_inq_varid(nf_fid,'windspeed',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeed'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,windspeed)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeed'
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      airline      "airline"
c      airline not available in wfo file....fill with 0
c
      do i = 1, recnum
        airline(i) = 0
      enddo
c
c     variable        netcdf long name
c      bounceerror  "aircraft altitude variance error"
c      bounceerror not available in wfo file....fill with 45
c
      do i = 1, recnum
        bounceerror(i) = 45
      enddo
c
c     variable        netcdf long name
c      correctedflag"corrected data indicator"
c      correctedflag not available in wfo file....fill with 114
c
      do i = 1, recnum
        correctedflag(i) = 114
      enddo
c
c     variable        netcdf long name
c      interpolatedll "ups ascent/descent lat&lon interpolation indicator"
c      interpolatedll not available in wfo file....fill with 114
c
      do i = 1, recnum
        interpolatedll(i) = 114
      enddo
c
c     variable        netcdf long name
c      interpolatedtime "ups ascent/descent time interpolation indicator"
c      interpolatedtime not available in wfo file....fill with 114
c
      do i = 1, recnum
        interpolatedtime(i) = 114
      enddo
c
c     variable        netcdf long name
c      missinginputminutes "missing minutes of input data"
c      missinginputminutes not available in wfo file....fill with 0
c
      missinginputminutes = 0
c
c     variable        netcdf long name
c      rollflag     "aircraft roll angle flag "
c
        nf_status = nf_inq_varid(nf_fid,'rollquality',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rollquality'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,rollflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rollflag'
      endif

      do i = 1, recnum
        if (rollflag(i) .eq. 0) then
          rollflag(i) = 71
        elseif (rollflag(i) .eq. 1) then
          rollflag(i) = 66
        else
          rollflag(i) = 78
        endif
      enddo
 
c
c     variable        netcdf long name
c      speederror   "aircraft ground speed error"
c      speederror not available in wfo file....fill with 45
c
      do i = 1, recnum
        speederror(i) = 45
      enddo
c
c     variable        netcdf long name
c      temperror    
c      temperror not available in wfo file....fill with 45
c
      do i = 1, recnum
        temperror(i) = 45
      enddo
c
c     variable        netcdf long name
c      winddirerror 
c      winddirerror not available in wfo file....fill with 45
c
      do i = 1, recnum
        winddirerror(i) = 45
      enddo
c
c     variable        netcdf long name
c      windspeederror
c      windspeederror not available in wfo file....fill with 45
c
      do i = 1, recnum
        windspeederror(i) = 45
      enddo
c
c   variables of type double
c
c
c     variable        netcdf long name
c      maxsecs      "maximum observation time"
c      maxsecs not available in wfo file....fill with 0
c
      maxsecs = 0
c
c     variable        netcdf long name
c      minsecs      "minimum observation time"
c      minsecs not available in wfo file....fill with 0
c
      minsecs = 0
c
c     variable        netcdf long name
c      timeobs      "time of observation"
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
c     variable        netcdf long name
c      timereceived "time data was received at ground station"
c
        nf_status = nf_inq_varid(nf_fid,'timereceived',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timereceived'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,timereceived)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timereceived'
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c      destairport  "destination airport"
c      destairport not available in wfo file....fill with "   " (3 spaces)
c
      do i = 1, recnum
        destairport(i) = "   "
      enddo
c
c     variable        netcdf long name
c      flight       "flight number"
c
        nf_status = nf_inq_varid(nf_fid,'flight',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var flight'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,flight)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var flight'
      endif
c
c     variable        netcdf long name
c      maxdate      "maximum observation date"
c      maxdate not available in wfo file....fill with 30 blank spaces
c
c     variable        netcdf long name
c      mindate      "minimum observation date"
c      mindate not available in wfo file....fill with 30 blank spaces
c
      maxdate = "                              "
      mindate = "                              "
c
c     variable        netcdf long name
c      origairport  "originating airport"
c      origairport not available in wfo file....fill with "   " (3 spaces)
c
      do i = 1, recnum
        origairport(i) = "   "
      enddo
c
c     variable        netcdf long name
c      rptstation   "station reporting through"
c      rptstation not available in wfo file....fill with "    " (4 spaces)
c
      do i = 1, recnum
        rptstation(i) = "    "
      enddo
c
c     variable        netcdf long name
c      tailnumber   "tail number"
c
        nf_status = nf_inq_varid(nf_fid,'tailnumber',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tailnumber'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,tailnumber)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tailnumber'
      endif

c     used by acars program...not available on wfo
c       use fsl acars quality controls to fill
c
c     max temp:  if altitude > 35000ft, t < -20c
c                else   t < 60 - 80 * (altitude / 35000)
c     min temp:  if altitude < 18000ft, t > -60c
c                if altitude > 35000ft, t > -100c
c                else t > 60 - 40 * (altitude -18000) / 17000
c     wind dir:  0 <= winddir <= 360
c     wind spd:  windspeed(knots) >= 0
c     max wind spd:  bad = 0
c                    if (altitude < 30000.) {
c                      wmax = 70. + 230.*altitude / 30000.;
c                    } else if (altitude < 40000.) {
c                      wmax = 300.;
c                    } else if (altitude < 45000.) {
c                      wmax = 300. - 100 * (altitude - 40000.) / 5000. ;
c                    } else {
c                      wmax = 200.;
c                    }
c                    if (windspeed > wmax) {
c                      bad = 1;
c                    }
c
c     datadescriptor: 'r' if temp and wind within ranges above
c                     'x' if temp and wind failed any test    
c     if datadescriptor = 'x':
c     errortype:      'w' for wind
c                     't' for temperature
c                     'b' for both
c
      do i = 1, recnum
        wind_err = 0  ! no error
        temp_err = 0  ! no error
        temp = temperature(i) - 273.15  ! convert k to c

        if(abs(windspeed(i)) .lt. 1000.)then
            ws = windspeed(i) * 1.9438 ! convert m/s to knots
        endif

        if(abs(altitude(i)) .lt. 1e10)then
            alt = altitude(i) * 3.280839895 ! convert m to ft
        endif

        if (alt .gt. 35000.0) then
          max_temp = -20.0
        else
          max_temp = 80 * (alt / 35000.)
        endif
        
        if (alt .lt. 18000.0) then
          min_temp = -60.0
        elseif (alt .gt. 35000.0) then
          min_temp = -100.0
        else
          min_temp = 40 * (alt - 18000.) / 17000.
        endif
        
        if ((temp .le. min_temp) .or. (temp .ge. max_temp)) 
     1    temp_err = 1    ! bad temp

        if ((winddir(i) .lt. 0.0) .or. (winddir(i) .gt. 360.0))
     1    wind_err = 1    ! bad wind dir

        if (alt .lt. 30000.) then
          wmax = 70. + 230.*alt / 30000.
        elseif (alt .lt. 40000.) then
          wmax = 300.
        elseif (alt .lt. 45000.) then
          wmax = 300. - 100 * (alt - 40000.) / 5000. 
        else 
          wmax = 200.
        endif 

        if ((ws .lt. 0.0) .or. (ws .gt. wmax))
     1    wind_err = 1  ! bad wind speed

        if ((temp_err .eq. 0) .and. (wind_err .eq. 0)) then
          datadescriptor(i) = 82
          errortype(i) = 46
        elseif ((temp_err .eq. 1) .and. (wind_err .eq. 1)) then
          datadescriptor(i) = 88
          errortype(i) = 66
        elseif (temp_err .eq. 1) then
          datadescriptor(i) = 88
          errortype(i) = 84
        else  ! (wind_err .eq. 1)
          datadescriptor(i) = 88
          errortype(i) = 87
        endif
 
      enddo
c
c     variable        netcdf long name
c      watervaporqc "water vapor mixing ratio qc character"
c      watervaporqc not available in wfo file....fill with 0
c
      do i = 1, recnum
        watervaporqc(i) = 0
      enddo
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
