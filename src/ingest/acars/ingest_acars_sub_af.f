
      subroutine get_acars_afwa(i4time_sys,i4_acars_window
     1                                    ,nx_l,ny_l
     1                                    ,ext
     1                                    ,filename,istatus)

      character*(*) filename,ext

!.............................................................................

      character*6 c6_a1acid
      character*9 a9_timeobs,a9_recpttime 
!     character*7 c7_skycover
      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)
      real latitude,longitude

!............................................................................

      lun_in = 21

      open(lun_in,file=filename,status='old')
  
      i = 0

      do while (.true.)

          read(lun_in,101,err=890,end=999) !    name          units & factor
     1         i_a1cycc,                       
     1         i_a1type,
     1         i_a1dpd,
     1         i_a1gwc,
     1         i_a1jul,                ! julian hour          hr since 673650000
     1         i_a1min,                ! time-report-minutes
     1         i_a1lat,                ! latitude             deg * 100
     1         i_a1lon,                ! longitude            deg * -100 
     1         i_a1kind,
     1         i_a1alt,                ! flight-alt (pressure) meters msl 
     1         i_a1pla,
     1         i_a1dval,
     1         i_a1holm,
     1         i_a1fltp,               ! temperature          kelvins * 10
     1         i_a1wd,                 ! wind-direction       deg
     1         i_a1wfls,               ! wind-speed           m/s * 10
     1         c6_a1acid 
 101      format(16(i9,2x),a6)

          i = i + 1

          write(6,*)
          write(6,*)' acars #',i

          latitude  =  float(i_a1lat)/100.
          longitude = +float(i_a1lon)/100.
          altitude  =  i_a1alt

!         write(6,*)' location = '
!    1             ,latitude,longitude,altitude

!         if(latitude  .le. rnorth .and. latitude  .ge. south .and.
!    1       longitude .ge. west   .and. longitude .le. east      
!    1                                                             )then       
!             continue
!         else ! outside lat/lon perimeter - reject
!             write(6,*)' lat/lon - reject'       
!!   1                 ,latitude,longitude
!             goto 900
!         endif

          if(altitude .gt. 20000.)then
              write(6,*)' altitude is suspect - reject',altitude
              goto 900
          endif

!         if(abs(timeobs)      .lt. 3d9       .and.
!    1       abs(timereceived) .lt. 3d9              )then
!             call c_time2fname(nint(timeobs),a9_timeobs)
!             call c_time2fname(nint(timereceived),a9_recpttime)
!         else
!             write(6,*)' bad observation time - reject'       
!    1                   ,timeobs,timereceived
!             goto 900
!         endif

          call afwa_julhr_i4time(i_a1jul,i_a1min,i4time_ob)

          call make_fnam_lp(i4time_ob,a9_timeobs,istatus)
          if(istatus .ne. 1)goto900

          a9_recpttime = '         '

!         call cv_asc_i4time(a9_timeobs,i4time_ob)

          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. i4_acars_window)then ! outside time window
              write(6,*)' time - reject '
     1           ,a9_timeobs,i4_resid,i4_acars_window
              goto 900        
          endif

          call open_ext(31,i4time_sys,ext(1:3),istatus)       

          write(6,1)a9_timeobs,a9_recpttime 
          write(31,1)a9_timeobs,a9_recpttime 
 1        format(' time - prp/rcvd:'/1x,a9,2x,a9) 

          write(6,2)latitude,longitude,altitude
          write(31,2)latitude,longitude,altitude
 2        format(' lat, lon, altitude'/f8.3,f10.3,f8.0)  

!         test for bad winds
!         if(char(datadescriptor) .eq. 'x')then
!           if(char(errortype) .eq. 'w' .or. 
!    1         char(errortype) .eq. 'b'                         )then
!             write(6,*)' qc flag is bad - reject '
!    1                 ,char(datadescriptor),char(errortype)
!             goto 850
!           endif
!         endif

          windspeed = float(i_a1wfls) / 10.
          winddir = i_a1wd

          if(abs(windspeed) .gt. 250. .or. winddir .gt. 360.)then
              write(6,*)' wind is suspect - reject',winddir,windspeed

          else ! write out valid wind
              write(6,3)int(winddir),windspeed
              write(31,3)int(winddir),windspeed
 3            format(' wind:'/' ', i3, ' deg @ ', f6.1, ' m/s')

          endif

          temperature = float(i_a1fltp)/10.

 850      if(abs(temperature) .lt. 400.)then
              write(6,13)temperature
              write(31,13)temperature
 13           format(' temp:'/1x,f10.1)
       
          else
              write(6,*)' temperature is suspect - reject'
     1                , temperature

          endif

!         if(watervapormr .ge. 0. .and. 
!    1       watervapormr .le. 100.)then
!             write(6,23)watervapormr
!             write(31,23)watervapormr
!23           format(' mixr:'/1x,f10.3)

!         else
!             write(6,*)' water vapor rejected: ',watervapormr
!
!         endif

          go to 900

 890      write(6,*)' warning: read error'

 900  enddo ! read line of afwa file

!............................................................................

 999  write(6,*)' end of afwa file detected'

!     close(lun_in)
      istatus = 1
      return
      end
