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
      subroutine read_buoy_cwb ( filename ,maxskycover, recnum,
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

      integer  logicrecnum(recnum)
      integer  windqua(recnum), temperaturequa(recnum)
      integer  sealevelpressqua(recnum), seasurfacetempqua(recnum)
      integer  dewpointqua(recnum)

      real  relahumility(recnum)

      istatus = 0
      n= 0
 
      open ( 1, file=filename, status='old', err=1000 )

      istatus= 1

      do j= 1,recnum
         read ( 1, 10, end=99, err=999 )
     ~               reportflag(j), wmoid(j), latitude(j), longitude(j),       
     ~               yy(j), mo(j), dd(j), hh(j), mn(j), logicrecnum(j)
         read (1,20) winddir(j), windspeed(j), windqua(j),
     ~               temperature(j), temperaturequa(j),
     ~               sealevelpress(j), sealevelpressqua(j),
     ~               presschangechar(j), presschange3hour(j)
         read (1,30) seasurfacetemp(j), seasurfacetempqua(j),
     ~               dewpoint(j), dewpointqua(j)
         read (1,40) relahumility(j)

         do 5 i= 1,logicrecnum(j)-4
5           read (1,*)

         if ( reportflag(j) .ne. '*81' )  then
            write (6,*) 'read buoy data heading error'
            go to 1000
         endif

         n= n+1
      enddo

10    format ( a3, i5, 4x, 2f5.2, 2x, 5a2, i3 )
!hj: w>=d+3. 10/14/2013
20    format ( 2f3.0, i1, f4.1, i1, f5.1, 2i1, f4.1 ) 
30    format ( f4.1, i1, 22x, f5.1, i1 )
40    format ( 30x, f3.0 )

c      ----------       examing data quality and changing units       ---------
99    do j= 1,n
         if ( windqua(j) .eq. 9 )  then
            winddir(j)= badflag
            windspeed(j)= badflag
         endif

         if ( sealevelpressqua(j) .eq. 1 )  then
               sealevelpress(j)= sealevelpress(j) *100.         ! mb -> pascal
            else
               sealevelpress(j)= badflag
         endif
  
         if ( presschange3hour(j) .eq. 1 )  then
               presschange3hour(j)= presschange3hour(j) *100.   ! mb -> pascal
            else
               presschange3hour(j)= badflag
         endif

         if ( seasurfacetempqua(j) .eq. 1 )  then
               seasurfacetemp(j)= seasurfacetemp(j) +273.15     ! degc -> degk
            else
               seasurfacetemp(j)= badflag
         endif

         if ( dewpointqua(j) .eq. 1 )  then
               dewpoint(j)= dewpoint(j) +273.15                 ! degc -> degk
            else
               if ( temperature(j) .ne. -999.9   .and. 
     ~              relahumility(j) .ne. -99. )  then
                  dewpoint(j)= dwpt( temperature(j), relahumility(j) )
                  dewpoint(j)= dewpoint(j) +273.15              ! degc -> degk
                else
                  dewpoint(j)= badflag
                  relahumility(j)= badflag
               endif
         endif

         if ( temperaturequa(j) .eq. 1 )  then
               temperature(j)= temperature(j) +273.15           ! degc -> degk
            else
               temperature(j)= badflag
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

c               -------      dealing with lacking of data      -------
      do j= 1,n
         presweather(j)= "unk"
         skycover(1,j)=  "unk"
         skycover(2,j)=  "unk"
         stationname(j)= "unk"

         presschangechar(j)= int(badflag)
         dataplatformtype(j)= 0

         elevation(j)= 0.
         equivwindspeed10m(j)= badflag
         precip1hour(j)= badflag
         precip24hour(j)= badflag
         precip3hour(j)= badflag 
         precip6hour(j)= badflag 
         skylayerbase(1,j)= badflag 
         skylayerbase(2,j)= badflag 
         visibility(j)= badflag
         wetbulbtemperature(j)= badflag
         windgust(j)= badflag
      enddo

      go to 1000 

999   write (6,*) ' error reading buoy file'
      do j= 1,n
         write(6,*) reportflag(j), wmoid(j),latitude(j), longitude(j),
     ~              yy(j), mo(j), dd(j), hh(j), mn(j), logicrecnum(j)
         write(6,*) winddir(j), windspeed(j), windqua(j),
     ~              temperature(j), temperaturequa(j),
     ~              sealevelpress(j), sealevelpressqua(j),
     ~              presschangechar(j), presschange3hour(j)
         write(6,*) seasurfacetemp(j), seasurfacetempqua(j),
     ~              dewpoint(j), dewpointqua(j), timeobs(j)
         write(6,*) relahumility(j)
      enddo

1000  return
      end
