          subroutine write_aircraft_sub(lun,ext
     1                          ,a9_timeobs,a9_recpttime
     1                          ,i4time_sys
     1                          ,i4time_earliest          
     1                          ,i4time_latest            
     1                          ,latitude,longitude,altitude
     1                          ,winddir,windspeed
     1                          ,temperature,relhumidity
     1                          ,l_geoalt
     1                          ,l_debug                           ! i
     1                          ,istat_ob)                         ! o

          character*(*) ext 

          character*9 a9_timeobs,a9_recpttime

          logical l_debug,l_geoalt

          real latitude, longitude

          istat_ob = 0

          call open_ext(lun,i4time_sys,ext(1:3),istatus)       

!         test the altitude
          if(nanf(altitude) .eq. 1)then
              if(l_debug)write(6,*)' altitude failed nan test - reject'       
     1                            ,altitude
              goto 900
          endif

          if(altitude .gt. 20000. .or. altitude .lt. -1000.
     1                            .or. altitude .eq.     0.)then
              if(l_debug)write(6,*)' altitude is suspect - reject'
     1                            ,altitude
              goto 900
          endif

!         test the time
          call cv_asc_i4time(a9_timeobs,i4time_ob)
          if(i4time_ob .lt. i4time_earliest .or.
     1       i4time_ob .gt. i4time_latest        )then ! outside time window
              if(l_debug)write(6,*)' time - reject '
     1           ,a9_timeobs,i4time_ob,i4time_earliest,i4time_latest
              goto 900        
          endif

          if(l_debug)write(6,1)a9_timeobs,a9_recpttime 
                     write(lun,1)a9_timeobs,a9_recpttime 
 1        format(' time - prp/rcvd:'/1x,a9,2x,a9) 
          istat_ob = 1

          if(l_geoalt)then
              if(l_debug)write(6,2)latitude,longitude,altitude
              write(lun,2)          latitude,longitude,altitude
 2            format(' lat, lon, geoalt  '/f8.3,f10.3,f8.0)  
          else
              if(l_debug)write(6,3)latitude,longitude,altitude
              write(lun,3)          latitude,longitude,altitude
 3            format(' lat, lon, altitude'/f8.3,f10.3,f8.0)  
          endif

!         test for bad winds
!         if(char(datadescriptor) .eq. 'x')then
!           if(char(errortype) .eq. 'w' .or. 
!    1         char(errortype) .eq. 'b'                         )then
!             if(l_debug)write(6,*)' qc flag is bad - reject wind'
!    1                 ,char(datadescriptor),char(errortype)
!             goto 850
!           endif
!         endif

          if(abs(windspeed) .gt. 250.)then
              if(l_debug)write(6,*)' wind speed is suspect - reject'
     1                              ,windspeed

          elseif(int(winddir).lt.0 .or. int(winddir).gt.360)then     
              if(l_debug)write(6,*)' wind direction is suspect - reject'       
     1                              ,winddir

          else ! write out valid wind
              if(l_debug)then
                  write(6,4)int(winddir),windspeed
              endif

              write(lun,4)int(winddir),windspeed
 4            format(' wind:'/' ', i3, ' deg @ ', f6.1, ' m/s')     

          endif

!         test/write temperature
          if(abs(temperature) .lt. 400.)then
              if(l_debug)write(6,13)temperature
                         write(lun,13)temperature
 13           format(' temp:'/1x,f10.1)
       
          else
              if(l_debug)write(6,*)' temperature is suspect - reject'
     1                             , temperature

          endif

!         test/write relative humidity
          if(relhumidity     .ge. 0.   .and. 
     1       relhumidity     .le. 1.00 
!    1       watervaporqc(i) .le. 2    .and.
!    1       watervaporqc(i) .ge. 0             
     1                                       )then
              if(l_debug)write(6,23)relhumidity
!             if(l_debug)write(6,*)' rh qc value = ',watervaporqc(i)
              write(lun,23)relhumidity
 23           format(' rh:'/1x,f10.3)

          else
              if(l_debug)write(6,*)' rh rejected: '
     1                             ,relhumidity ! ,watervaporqc(i)

          endif

 900      continue

!............................................................................

      return
      end
