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
      subroutine read_ship_cwb ( filename, maxskycover, recnum, 
     +     dataplatformtype, dewpoint, elevation, equivwindspeed10m,
     +     latitude, longitude, precip1hour, precip24hour, precip3hour,        
     +     precip6hour, presweather, presschange3hour, presschangechar,
     +     sealevelpress, seasurfacetemp, skycover, skylayerbase, 
     +     stationname, temperature, timeobs, visibility, 
     +     wetbulbtemperature, winddir, windgust, windspeed, wmoid,
     +     badflag, n, istatus )
 
      integer  recnum

      character*(*)  filename
      character*25   presweather(recnum)
      character*8    skycover(maxskycover,recnum), stationname(recnum)

      integer  dataplatformtype(recnum), presschangechar(recnum)
      integer  wmoid(recnum)

      real  dewpoint(recnum), elevation(recnum)
      real  equivwindspeed10m(recnum), latitude(recnum)
      real  longitude(recnum), precip1hour(recnum)
      real  precip24hour(recnum), precip3hour(recnum)
      real  precip6hour(recnum), presschange3hour(recnum)
      real  sealevelpress(recnum), seasurfacetemp(recnum)
      real  skylayerbase(maxskycover,recnum), temperature(recnum)
      real  visibility(recnum), wetbulbtemperature(recnum)
      real  winddir(recnum), windgust(recnum), windspeed(recnum)

      double precision  timeobs(recnum)

      character*3   reportflag(recnum)
      character*2   yy(recnum), mo(recnum), dd(recnum)
      character*2   hh(recnum), mn(recnum)
      character*10  time(recnum)
      character*9   a10_to_a9

      integer  windqua(recnum), temperaturequa(recnum)
      integer  sealevelpressqua(recnum), seasurfacetempqua(recnum)
      integer  dewpointqua(recnum), presschange3hourqua(recnum)
      integer  wetbulbtemperaturequa(recnum)

      real  relahumility(recnum), tempdewdiff(recnum)

      istatus = 0
      n= 0
 
      open ( 1, file=filename, status='old', err=1000 )

      istatus= 1

      do j= 1,recnum
         read ( 1, 10, end=99, err=999 ) reportflag(j),
     ~               stationname(j), latitude(j), longitude(j),
     ~               yy(j), mo(j), dd(j), hh(j), mn(j)
         read (1,20) winddir(j), windspeed(j), windqua(j),
     ~               visibility(j), presweather(j),
     ~               sealevelpress(j), sealevelpressqua(j),
     ~               temperature(j), temperaturequa(j),
     ~               skycover(1,j), skylayerbase(1,j)
         read (1,30) tempdewdiff(j), dewpointqua(j), presschangechar(j),       
     ~               presschange3hour(j), presschange3hourqua(j),
     ~               precip3hour(j), wetbulbtemperature(j),
     ~               wetbulbtemperaturequa(j), windgust(j)
         read (1,40) skycover(2,j), skylayerbase(2,j), precip24hour(j)
         read (1,50) seasurfacetemp(j), seasurfacetempqua(j)

         if ( reportflag(j) .ne. '*34' )  then
            write (6,*) 'read ship data heading error'
            go to 1000
         endif

         n= n+1
      enddo

!hj: w>=d+3. 10/14/2013
10    format ( a3, a5, 4x, 2f5.2, 2x, 5a2 )
20    format ( 2x, 2f3.0, i1, f23.0, a2, 3x, f5.1, i1, f4.1, i1, a2, 2x,
     ~         f3.0 )
30    format ( f4.1, i1, 1x, i2, f4.1, i1, 1x, f4.1, 10x, f4.1, i1, 3x,
     ~         f3.0 )
40    format ( a2, 2x, f3.0, 14x, f4.1 )
50    format ( 16x, f4.1, i1 )

c      ----------       examing data quality and changing units       ---------
99    do j= 1,n
         if ( windqua(j) .eq. 9 )  then
            winddir(j)= badflag
            windspeed(j)= badflag
         endif

         if ( windgust(j) .eq. -99. )  windgust(j)= badflag
         if ( elevation(j) .eq. -999. )  elevation(j)= badflag
         if ( precip3hour(j) .eq. -99.9 )  precip3hour(j)= badflag
         if ( presschangechar(j).eq.-9 ) presschangechar(j)=int(badflag)
         if ( presweather(j) .eq. '-9' )  presweather(j)= 'unk'

         if ( sealevelpressqua(j) .eq. 1 )  then
               sealevelpress(j)= sealevelpress(j) *100.   ! millibar -> pascal
            else
               sealevelpress(j)= badflag
         endif
  
         if ( presschange3hourqua(j) .eq. 1 )  then       ! millibar -> pascal
               presschange3hour(j)= presschange3hour(j) *100.
            else
               presschange3hour(j)= badflag
         endif

         if ( seasurfacetempqua(j) .eq. 1 )  then
               seasurfacetemp(j)= seasurfacetemp(j) +273.15     ! degc -> degk
            else
               seasurfacetemp(j)= badflag
         endif

         if ( temperaturequa(j) .eq. 1 )  then
               temperature(j)= temperature(j) +273.15           ! degc -> degk
            else
               temperature(j)= badflag
         endif

         if ( dewpointqua(j) .eq. 1 )  then
               dewpoint(j)= temperature(j) -tempdewdiff(j)        ! unit: degk
            else
               dewpoint(j)= badflag
         endif

         if ( wetbulbtemperaturequa(j) .eq. 1 )  then           ! degc -> degk
               wetbulbtemperature(j)= wetbulbtemperature(j) +273.15
            else
               wetbulbtemperature(j)= badflag
         endif

         if ( precip24hour(j) .eq. -99.9 )  then
               precip24hour(j)= badflag
            else
               precip24hour(j)= precip24hour(j) *0.001   ! millimeter -> meter
         endif

         if ( yy(j)(1:1) .eq. ' ' )  yy(j)= '0'//yy(j)(2:2)
         if ( mo(j)(1:1) .eq. ' ' )  mo(j)= '0'//mo(j)(2:2)
         if ( dd(j)(1:1) .eq. ' ' )  dd(j)= '0'//dd(j)(2:2)
         if ( hh(j)(1:1) .eq. ' ' )  hh(j)= '0'//hh(j)(2:2)
         if ( mn(j)(1:1) .eq. ' ' )  mn(j)= '0'//mn(j)(2:2)
         time(j)= yy(j)//mo(j)//dd(j)//hh(j)//mn(j)
         call cv_asc_i4time( a10_to_a9(time(j),istatus), i4time )
         timeobs(j)= dble( i4time )                       ! seconds since 1960
      enddo

c    -------     code figure transformed into visibility ( unit: m )  -------
      do j= 1,n
         if     ( visibility(j) .eq.  0. )  then
               visibility(j)= 50.
         elseif ( visibility(j) .gt.  0.  .and.
     ~            visibility(j) .lt. 51. )  then
               visibility(j)= visibility(j) *100.
         elseif ( visibility(j) .gt. 55.  .and.
     ~            visibility(j) .lt. 81. )  then
               visibility(j)= ( visibility(j) -50. ) *1000.
         elseif ( visibility(j) .gt. 80.  .and.
     ~            visibility(j) .lt. 90. )  then
               visibility(j)= ( visibility(j) -74. ) *5000.
         elseif ( visibility(j) .eq. 90. )  then
               visibility(j)= 25.
         elseif ( visibility(j) .eq. 91. )  then
               visibility(j)= 50.
         elseif ( visibility(j) .eq. 92. )  then
               visibility(j)= 200.
         elseif ( visibility(j) .eq. 93. )  then
               visibility(j)= 500.
         elseif ( visibility(j) .eq. 94. )  then
               visibility(j)= 1000.
         elseif ( visibility(j) .eq. 95. )  then
               visibility(j)= 2000.
         elseif ( visibility(j) .eq. 96. )  then
               visibility(j)= 4000.
         elseif ( visibility(j) .eq. 97. )  then
               visibility(j)= 10000.
         elseif ( visibility(j) .eq. 98. )  then
               visibility(j)= 20000.
         elseif ( visibility(j) .eq. 99. )  then
               visibility(j)= 50000.
         else
               visibility(j)= badflag
         endif

c  -----  code figure transformed into base of lowest cloud ( unit: m )  -----
         if     ( skylayerbase(1,j) .eq. 0. )  then
               skylayerbase(1,j)= 25.
         elseif ( skylayerbase(1,j) .eq. 1. )  then
               skylayerbase(1,j)= 75.
         elseif ( skylayerbase(1,j) .eq. 2. )  then
               skylayerbase(1,j)= 150.
         elseif ( skylayerbase(1,j) .eq. 3. )  then
               skylayerbase(1,j)= 250.
         elseif ( skylayerbase(1,j) .eq. 4. )  then
               skylayerbase(1,j)= 450.
         elseif ( skylayerbase(1,j) .eq. 5. )  then
               skylayerbase(1,j)= 800.
         elseif ( skylayerbase(1,j) .eq. 6. )  then
               skylayerbase(1,j)= 1250.
         elseif ( skylayerbase(1,j) .eq. 7. )  then
               skylayerbase(1,j)= 1750.
         elseif ( skylayerbase(1,j) .eq. 8. )  then
               skylayerbase(1,j)= 2250.
         elseif ( skylayerbase(1,j) .eq. 9. )  then
               skylayerbase(1,j)= 3000.
         else
               skylayerbase(1,j)= badflag
         endif

c ---- code figure transformed into base of cloud layer indicated (unit: m) ----
         if     ( skylayerbase(2,j) .eq.  0. )  then
               skylayerbase(2,j)= 15.
         elseif ( skylayerbase(2,j) .gt.  0.  .and.
     ~            skylayerbase(2,j) .lt. 51. )  then
               skylayerbase(2,j)= skylayerbase(2,j) *30.
         elseif ( skylayerbase(2,j) .gt. 55.  .and.
     ~            skylayerbase(2,j) .lt. 81. )  then
               skylayerbase(2,j)= ( skylayerbase(2,j) -50. ) *300.
         elseif ( skylayerbase(2,j) .gt. 80.  .and.
     ~            skylayerbase(2,j) .lt. 90. )  then
               skylayerbase(2,j)= ( skylayerbase(2,j) -74. ) *1500.
         elseif ( skylayerbase(2,j) .eq. 90. )  then
               skylayerbase(2,j)= 25.
         elseif ( skylayerbase(2,j) .eq. 91. )  then
               skylayerbase(2,j)= 75.
         elseif ( skylayerbase(2,j) .eq. 92. )  then
               skylayerbase(2,j)= 150.
         elseif ( skylayerbase(2,j) .eq. 93. )  then
               skylayerbase(2,j)= 250.
         elseif ( skylayerbase(2,j) .eq. 94. )  then
               skylayerbase(2,j)= 450.
         elseif ( skylayerbase(2,j) .eq. 95. )  then
               skylayerbase(2,j)= 800.
         elseif ( skylayerbase(2,j) .eq. 96. )  then
               skylayerbase(2,j)= 1250.
         elseif ( skylayerbase(2,j) .eq. 97. )  then
               skylayerbase(2,j)= 1750.
         elseif ( skylayerbase(2,j) .eq. 98. )  then
               skylayerbase(2,j)= 2250.
         elseif ( skylayerbase(2,j) .eq. 99. )  then
               skylayerbase(2,j)= 3000.
         else
               skylayerbase(2,j)= badflag
         endif
      enddo

c --- code figure of total cloud cover or      transformed into metar format ---
c                    special cloud layer cover
      do 100 j= 1,n
      do 100 i= 1,maxskycover
         if     ( skycover(i,j) .eq. ' 0' )  then
            skycover(i,j)= 'clr'
         elseif ( skycover(i,j) .eq. ' 1'  .or.
     ~            skycover(i,j) .eq. ' 2' )  then
            skycover(i,j)= 'few'
         elseif ( skycover(i,j) .eq. ' 3'  .or.
     ~            skycover(i,j) .eq. ' 4' )  then
            skycover(i,j)= 'sct'
         elseif ( skycover(i,j) .eq. ' 5'  .or.
     ~            skycover(i,j) .eq. ' 6'  .or.
     ~            skycover(i,j) .eq. ' 7' )  then
            skycover(i,j)= 'bkn'
         elseif ( skycover(i,j) .eq. ' 8' )  then
            skycover(i,j)= 'ovc'
         else
            skycover(i,j)= '   '
         endif
100   continue

c               -------      dealing with lacking of data      -------
      do j= 1,n
         presweather(j)= "unk"

         dataplatformtype(j)= 1

         elevation(j)= 0.
         equivwindspeed10m(j)= badflag
         precip1hour(j)= badflag
         precip6hour(j)= badflag 
         wmoid(j)= badflag 
      enddo

      go to 1000 

999   write (6,*) ' error reading ship file'
      do j= 1,n
         write(6,*)  reportflag(j),
     ~               stationname(j), latitude(j), longitude(j),
     ~               yy(j), mo(j), dd(j), hh(j), mn(j)
         write(6,*)  winddir(j), windspeed(j), windqua(j),
     ~               visibility(j), presweather(j),
     ~               sealevelpress(j), sealevelpressqua(j), 
     ~               temperature(j), temperaturequa(j),
     ~               skycover(1,j), skylayerbase(1,j)
         write(6,*)  dewpoint(j), dewpointqua(j), presschangechar(j),       
     ~               presschange3hour(j), presschange3hourqua(j),
     ~               precip3hour(j), wetbulbtemperature(j),
     ~               wetbulbtemperaturequa(j), windgust(j)
         write(6,*)  skycover(2,j), skylayerbase(2,j), precip24hour(j)
      enddo

1000  return
      end
