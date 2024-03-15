      subroutine read_metar_cwb(
     &     filename, maxskycover, recnum, altimeter,
     &     autostationtype, dewpoint, dpfromtenths, elevation,
     &     latitude, longitude, maxtemp24hour, mintemp24hour,
     &     precip1hour, precip24hour, precip3hour, precip6hour,
     &     presweather, presschange3hour, presschangechar,
     &     reporttype, sealevelpress, skycover, skylayerbase,
     &     snowcover, stationname, tempfromtenths, temperature,
     &     timeobs, visibility, winddir, windgust, windspeed, wmoid,
     &     badflag, n, istatus )

      implicit none

      integer  maxskycover, recnum

      character*6   autostationtype(recnum), reporttype(recnum)
      character*8   dmycover(10,recnum), skycover( maxskycover, recnum)
      character*25  presweather(recnum)
      character*5   stationname(recnum)

      integer  presschangechar(recnum), wmoid(recnum)
      integer  istatus

      double precision  timeobs(recnum)

      real  altimeter(recnum), dewpoint(recnum), dpfromtenths(recnum)
      real  elevation(recnum), latitude(recnum), longitude(recnum)
      real  maxtemp24hour(recnum), mintemp24hour(recnum)
      real  precip1hour(recnum), precip24hour(recnum)
      real  precip3hour(recnum), precip6hour(recnum)
      real  presschange3hour(recnum), sealevelpress(recnum)
      real  skylayerbase( maxskycover, recnum), snowcover(recnum)
      real  tempfromtenths(recnum), temperature(recnum)
      real  visibility(recnum), winddir(recnum), windgust(recnum)
      real  windspeed(recnum)
      real  badflag

      character*(*) filename
      character*2   yy(recnum), m1(recnum), dd(recnum)
      character*2   hh(recnum), m2(recnum)
      character*10  time(recnum)
      character*9   a10_to_a9

      integer windqua(recnum), windgustqua(recnum)
      integer temperaturequa(recnum), dewpointqua(recnum)
      integer altimeterqua(recnum)
      integer i, j, n, i4time

      real  dmylayerbase(10,recnum)

      integer    stnnum, len_dir
      parameter  ( stnnum=4503 )

      character*100  stn_directory, stn_filename
      character*4    stn(stnnum)

      real         stnelevation(stnnum)

      istatus= 0
      n= 0

      open ( 1, file=filename, status='old', err=1000 )
      call get_directory ( 'static', stn_directory, len_dir )
      stn_filename= stn_directory(1:len_dir) // 'metarstn.dat'
      open ( 2, file=stn_filename, status='old', err=1000 )

      istatus= 1

      do j= 1,recnum
         read ( 1, 10, end=99, err=999 )
     *   hh(j), m2(j), stationname(j), latitude(j),
     *   longitude(j), winddir(j), windspeed(j), windqua(j),
     *   windgust(j), windgustqua(j), visibility(j), presweather(j),
     *   ( dmycover(i,j), dmylayerbase(i,j), i=1,7 ), sealevelpress(j),
     *   temperature(j), temperaturequa(j), dewpoint(j), dewpointqua(j),
     *   altimeter(j), altimeterqua(j), precip1hour(j), yy(j), m1(j),
     *   dd(j)
         n= n+1
      enddo
!hj: w>=d=3. 10/14/2013
10    format( 2a2, a4, 2f5.2, 2f5.0, i1, f5.0, i1, 2x, f5.0, 5x, a2,
     *        7(a3,f5.0), 3x, f5.1, 16x, 3(f5.0,i1), 10x, f6.3, 12x, 
     *        3a2 )

c         -------       read the elevations of all stations      -------
99    read (2,*)
      read (2,*)
      read (2,'(a4,11x,f4.0)') ( stn(j), stnelevation(j), j=1,stnnum )

      do 20 j= 1,n
      do 20 i= 1,stnnum
20       if ( stationname(j) .eq. stn(i) )  elevation(j)=stnelevation(i)

c      ----------       examing data quality and changing units       ---------
      do j= 1,n

         if ( windqua(j) .ne. 1 )  then
            winddir(j)= badflag
            windspeed(j)= badflag
         endif
         if ( windgustqua(j) .ne. 1 )  windgust(j)= badflag
         if ( visibility(j) .eq. -9999. )  visibility(j)= badflag
         if ( precip1hour(j) .eq. -9.999 )  precip1hour(j)= badflag
         if ( presweather(j) .eq. '-9' )  presweather(j)= '  '

         if ( temperaturequa(j) .eq. 1 )  then
               temperature(j)= temperature(j) +273.15  ! degc -> degk
            else
               temperature(j)= badflag
         endif

         if ( dewpointqua(j) .eq. 1 )  then
               dewpoint(j)= dewpoint(j) +273.15        ! degc -> degk
            else
               dewpoint(j)= badflag
         endif

         if ( altimeterqua(j) .eq. 1 )  then
               altimeter(j)= altimeter(j) *100         ! mb -> pascal
            else
               altimeter(j)= badflag
         endif

         if ( yy(j)(1:1) .eq. ' ' )  yy(j)= '0'//yy(j)(2:2)
         if ( m1(j)(1:1) .eq. ' ' )  m1(j)= '0'//m1(j)(2:2)
         if ( dd(j)(1:1) .eq. ' ' )  dd(j)= '0'//dd(j)(2:2)
         time(j)= yy(j)//m1(j)//dd(j)//hh(j)//m2(j)
         call cv_asc_i4time( a10_to_a9(time(j),istatus), i4time )
         timeobs(j)= dble( i4time )                    ! seconds since 1960

      enddo

      do 30 j= 1,n
      do 30 i= 1,maxskycover
         skycover(i,j)= dmycover(i,j)
         skylayerbase(i,j)= dmylayerbase(i,j)

         if ( skycover(i,j) .eq. '-99' )  skycover(i,j)= '   '

         if ( skylayerbase(i,j) .eq. -9999. )  then
               skylayerbase(i,j)= badflag
	       skycover(i,j)= '   '
         elseif ( skylayerbase(i,j) .eq. 0. )  then
               skylayerbase(i,j)= 15.
         elseif ( skylayerbase(i,j) .eq. 999. )  then
               skylayerbase(i,j)= 30000.
         else
               skylayerbase(i,j)= skylayerbase(i,j) *30.        ! unit: m
         endif
30    continue

c               -------      dealing with lacking of data      -------
      do j= 1,n
         autostationtype(j)= "unk"
         reporttype(j)= "metar"

         presschangechar(j)= int(badflag)
         wmoid(j)= int(badflag)

         dpfromtenths(j)= badflag
         maxtemp24hour(j)= badflag
         mintemp24hour(j)= badflag
         precip24hour(j)= badflag
         precip3hour(j)= badflag
         precip6hour(j)= badflag
         presschange3hour(j)= badflag
         sealevelpress(j)= badflag
         snowcover(j)= badflag
         tempfromtenths(j)= badflag
      enddo

      go to 1000

999   do j= 1,n
         write(6,*)
     *   hh(j), m2(j), stationname(j), latitude(j), longitude(j),
     *   winddir(j), windspeed(j), windqua(j),
     *   windgust(j), windgustqua(j), visibility(j), presweather(j),
     *   ( dmycover(i,j), dmylayerbase(i,j), i=1,10), temperature(j),
     *   temperaturequa(j), dewpoint(j), dewpointqua(j), altimeter(j),
     *   altimeterqua(j), precip1hour(j), yy(j), m1(j), dd(j), hh(j),
     *   m2(j), time(j), timeobs(j), elevation(j)
      enddo
      write(6,*)
     *      autostationtype, reporttype,
     *      presschangechar, wmoid,
     *      dpfromtenths, maxtemp24hour, mintemp24hour,
     *      precip24hour, precip3hour, precip6hour, presschange3hour,
     *      sealevelpress, snowcover, tempfromtenths, wmoid

1000  return
      end

