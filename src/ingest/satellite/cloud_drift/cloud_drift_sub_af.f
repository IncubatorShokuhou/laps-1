cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis

      subroutine get_cloud_drift_afwa(i4time_sys,i4_ob_window
     1                                ,nx_l,ny_l
     1                                ,filename,istatus)

      character*(*) filename

!.............................................................................

      character*9 a9_timeobs
!     character*7 c7_skycover
      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)
      real latitude,longitude

!............................................................................

      lun_in = 21
      lun_out = 11

      open(lun_in,file=filename,status='old')

      r_mspkt = .518  

      i = 0

      do while (.true.)

          read(lun_in,101,err=890,end=999) !    name          units & factor
     1         i_g9cyc,                ! db-cycle       
     1         i_g9gwcg,               ! gwc-region
     1         i_g9julg,               ! julian hour          hr since 673650000
     1         i_g9latg,               ! latitude             deg * 100
     1         i_g9long,               ! longitude            deg * 100 
     1         i_g9optc,               ! options-code-goes
     1         i_g9otg,                ! hhmm
     1         i_g9prgo,               ! pressure             pa * 10 (hpa*1000)
     1         i_g9rtgo,               ! report-type-id
     1         i_g9satg,               ! satellite operator
     1         i_g9sng,                ! satellite number
     1         i_g9tmgo,               ! temperature
     1         i_g9wdgo,               ! wind-direction       deg
     1         i_g9wsgo,               ! wind-speed           m/s * 10
     1         i_g9wuig                ! wind-speed units     1=kt, 2=m/s
 101      format(14(i9,2x),i9)

          i_g9min = 0

          i = i + 1

          write(6,*)
          write(6,*)' cdw #',i

          latitude  =  float(i_g9latg)/100.
          longitude = +float(i_g9long)/100.
          pres_pa   =  i_g9prgo / 10.

          write(6,2)latitude,longitude,pres_pa
 2        format(' lat, lon, pres_pa'/f8.3,f10.3,f8.0)  

          if(pres_pa .gt. 1e10)then
              write(6,*)' pres_pa is suspect - reject',pres_pa
              goto 900
          endif

          i_hr_ob  = i_g9otg / 100
          i_min_ob = i_g9otg - i_hr_ob*100
          i_hr_jul = i_g9julg - (i_g9julg/24)*24

          if(i_hr_ob .ne. i_hr_jul)then
              write(6,*)' error, i_hr_ob .ne. i_hr_jul',i_hr_ob,i_hr_jul
              goto900
          endif

          call afwa_julhr_i4time(i_g9julg,i_min_ob,i4time_ob)

          call make_fnam_lp(i4time_ob,a9_timeobs,istatus)
          if(istatus .ne. 1)goto900

          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. i4_ob_window)then ! outside time window
              write(6,*)' time - reject '
     1           ,a9_timeobs,i4_resid,i4_ob_window
              goto 900        
          endif

          winddir = i_g9wdgo

!         calculate wind speed in m/s
          if(i_g9wuig .eq. 1)then     ! kt
!             write(6,*)' convert wind speed from kt to m/s ',i_g9wuig       
              windspeed = (float(i_g9wsgo) / 10.) * r_mspkt

          elseif(i_g9wuig .eq. 2)then ! m/s
              windspeed = float(i_g9wsgo) / 10.

          else ! invalid units identifier
              write(6,*)' warning: wind speed units are undefined '
     1                 ,i_g9wuig,i_g9wsgo     
              goto900

          endif

          if(abs(windspeed) .gt. 250. .or. winddir .gt. 360.)then
              write(6,*)' wind is suspect - reject',winddir,windspeed
              goto900

          else ! write out valid wind
              write(6      ,21)latitude,longitude,pres_pa
     1                        ,winddir,windspeed,a9_timeobs       
              write(lun_out,21)latitude,longitude,pres_pa
     1                        ,winddir,windspeed,a9_timeobs       
 21           format(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)

          endif

          go to 900

 890      write(6,*)' warning: read error'

 900  enddo ! read line of afwa file

!............................................................................

 999  write(6,*)' end of afwa file detected'

      close(lun_in)
      istatus = 1
      return
      end
