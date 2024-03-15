      subroutine read_synop_cwb ( filename, maxskycvr, maxobs,
     &                          i4time_sys, path_to_local,            
     &                          altm, stntp, td, tdtths, elev,
     &                          lats, lons, t24max, t24min,
     &                          pcp1hr, pcp24hr, pcp3hr, pcp6hr, 
     &                          prswth, p, pc, pcc, rpttp, rh,
     &                          mslp, skycvr, skylyrbs, snowcvr, sr, st,       
     &                          stname, ttths, t, timeobs, vis,
     &                          dd, wgdd, ff, wgff, wmoid, badflag,
     &                          num, istatussynop )

      integer, parameter :: maxsynop = 150
      integer, parameter :: maxmso =    60

      character*(*)  filename
      character(*)   path_to_local
      character(25)  prswth(maxobs)
      character(13)  cvt_i4time_wfo_fname13,a13time_eat
      character(9)   a9time
      character(8)   skycvr(maxskycvr,maxobs)
      character(6)   rpttp(maxobs), stntp(maxobs)
      character(5)   stname(maxobs), stnno(maxobs)

      integer  pcc(maxobs), wmoid(maxobs), flag

      real  altm(maxobs), td(maxobs), tdtths(maxobs)
      real  elev(maxobs), lats(maxobs), lons(maxobs)
      real  t24max(maxobs), t24min(maxobs)
      real  rh(maxobs), pcp1hr(maxobs), pcp24hr(maxobs)
      real  pcp3hr(maxobs), pcp6hr(maxobs)
      real  p(maxobs), pc(maxobs), mslp(maxobs)
      real  skylyrbs(maxskycvr,maxobs), snowcvr(maxobs)
      real  ttths(maxobs), t(maxobs)
      real  vis(maxobs), dd(maxobs), wgdd(maxobs), wgff(maxobs)
      real  ff(maxobs), sr(maxobs), st(maxobs)

      double precision  timeobs(maxobs), timeobsmso(maxmso)

      character(6)   rpttpmso(maxobs), stntpmso(maxobs)
      character(5)   stnnomso(maxmso)
      integer  pccmso(maxmso)

      real  latsmso(maxmso), lonsmso(maxmso)
      real  elevmso(maxmso), tmso(maxmso), t24maxmso(maxmso)
      real  t24minmso(maxmso), tdmso(maxmso), rhmso(maxmso)
      real  pcp1hrmso(maxmso), pcp3hrmso(maxmso), pcp6hrmso(maxmso)
      real  pcp24hrmso(maxmso), ddmso(maxmso), ffmso(maxmso)
      real  wgddmso(maxmso), wgffmso(maxmso), pmso(maxmso)
      real  mslpmso(maxmso), pcmso(maxmso), srmso(maxmso), stmso(maxmso)           



      rpttp   = 'synop'
      stntp   = 'unk'
      skycvr  = '        '
      stname  = 'unk'
      stnno   = '     '
      wmoid   = ibadflag
      pcc     = ibadflag
      timeobs = badflag
      lats    = badflag
      lons    = badflag
      elev    = badflag
      t       = badflag
      ttths   = badflag
      t24max  = badflag
      t24min  = badflag
      td      = badflag
      tdtths  = badflag
      rh      = badflag
      pcp1hr  = badflag
      pcp3hr  = badflag
      pcp6hr  = badflag
      pcp24hr = badflag
      ff      = badflag
      dd      = badflag
      wgff    = badflag
      wgdd    = badflag
      p       = badflag
      mslp    = badflag
      pc      = badflag
      altm    = badflag
      skylyrbs= badflag
      vis     = badflag
      snowcvr = badflag
      sr      = badflag
      st      = badflag
      
      istatussynop= 0
      istatusmso=   0
      num=          0
      numsynop=     0
      nummso=       0

      nq= maxsynop

      call read_synop_cwb_sub ( filename, maxskycvr, maxsynop,
     &     altm(1:nq), stntp(1:nq), td(1:nq), tdtths(1:nq), elev(1:nq),
     &     lats(1:nq), lons(1:nq), t24max(1:nq), t24min(1:nq),
     &     pcp1hr(1:nq), pcp24hr(1:nq), pcp3hr(1:nq), pcp6hr(1:nq), 
     &     prswth(1:nq), pc(1:nq), pcc(1:nq), rpttp(1:nq), mslp(1:nq),
     &     skycvr(1:maxskycvr,1:nq), skylyrbs(1:maxskycvr,1:nq),
     &     snowcvr(1:nq), stname(1:nq), ttths(1:nq), t(1:nq),
     &     timeobs(1:nq), vis(1:nq), dd(1:nq), wgff(1:nq), ff(1:nq),
     &     wmoid(1:nq), badflag, numsynop, istatussynop )

      do i= 1,numsynop
         write(stnno(i),'(i5)') wmoid(i)      
         stname(i)(1:5)=stnno(i)(1:5)
      enddo

      np= numsynop +1
      nq= numsynop +maxmso

      call s_len ( path_to_local, len_inpath )
      path_to_local= path_to_local(1:len_inpath)//'mso/'

      call read_meso_cwb ( path_to_local, maxmso, badflag, ibadflag, 
     &                     i4time_sys, timeobsmso, rpttpmso, stntpmso,
     &                     stnnomso, latsmso, lonsmso, elevmso,
     &                     tmso, t24maxmso, t24minmso, tdmso, rhmso, 
     &                     pcp1hrmso, pcp3hrmso, pcp6hrmso, pcp24hrmso,
     &                     ddmso, ffmso, wgddmso, wgffmso, pmso, 
     &                     mslpmso, pccmso, pcmso, srmso, stmso, 
     &                     nummso, istatusmso )

c                    combine synop data and mesonet data 
      do i= 1,nummso ! added by shuyuan20100707 for checking meso rh values
         if(rhmso(i)>101) then
             rhmso(i)=  badflag       
         endif
      enddo
c

      do i= 1,numsynop ! added by shuyuan20100707 for checking synop rh data 
         if(rh(i)>101) then
             rh(i)=  badflag       
         endif
      enddo

      
      k= numsynop
      do i= 1,nummso
         flag= 0

         do j= 1,numsynop
            if ( stnnomso(i) == stnno(j) ) then
               rh(j)=     rhmso(i)
               dd(j)=     ddmso(i)
               ff(j)=     ffmso(i)
               wgdd(j)=   wgddmso(i)
               wgff(j)=   wgffmso(i)
               p(j)=      pmso(i)
               sr(j)=     srmso(i)
               st(j)=     stmso(i)
c
c modified by min-ken.hsieh
c must assign timeobsmso here
c or timeobs will be synop obs time
c also, assign mso t/td to synop t/td and ttths/tdtths
c or we will get t/td carried with synop which only update hourly
c
	       t(j)=	  tmso(i)
	       ttths(j)=  tmso(i)
	       td(j)=	  tdmso(i)
	       tdtths(j)= tdmso(i)
               timeobs(j)=timeobsmso(i)
	       pcp1hr(j)= pcp1hrmso(i)
	       pcp3hr(j)= pcp3hrmso(i)
	       pcp6hr(j)= pcp6hrmso(i)
	       pcp24hr(j)= pcp24hrmso(i)
               flag= 1
               exit
            endif
         enddo 


         if ( flag /= 1 ) then
            k= k +1

            timeobs(k)= timeobsmso(i)
            rpttp(k)=   rpttpmso(i)
            stntp(k)=   stntpmso(i)
            stnno(k)=   stnnomso(i)
            stname(k)=  stnnomso(i)
            wmoid(k)=   ibadflag
            td(k)=      tdmso(i)
            tdtths(k)=  tdmso(i)
            elev(k)=    elevmso(i)
            lats(k)=    latsmso(i)
            lons(k)=    lonsmso(i)
            t24max(k)=  t24maxmso(i)
            t24min(k)=  t24minmso(i)
            rh(k)=      rhmso(i)

            pcp1hr(k)=  pcp1hrmso(i)
            pcp24hr(k)= pcp24hrmso(i)
            pcp3hr(k)=  pcp3hrmso(i)
            pcp6hr(k)=  pcp6hrmso(i)
            p(k)=       pmso(i)
c         add if shuyuan  20100707 
            if(pcmso(i)>0.5  .and. pcmso(i)<1000)  then
              pc(k)=  pcmso(i)
            else
            pc(k)=badflag
            endif 
            if(pccmso(i)>-10000  .and. pccmso(i)<10000)  then   
              pcc(k)=     pccmso(i)
            else
              pcc(k)=badflag
            endif
            mslp(k)=    mslpmso(i)
            t(k)=       tmso(i)
            dd(k)=      ddmso(i)
            wgdd(k)=    wgddmso(i)
            ff(k)=      ffmso(i)
            wgff(k)=    wgffmso(i)
            sr(k)=      srmso(i)
            st(k)=      stmso(i)
         endif
      enddo

      num = k

      end



      subroutine read_synop_cwb_sub (
     &     filename, maxskycover, recnum, altimeter,
     &     autostationtype, dewpoint, dpfromtenths, elevation,
     &     latitude, longitude, maxtemp24hour, mintemp24hour,
     &     precip1hour, precip24hour, precip3hour, precip6hour,
     &     presweather, presschange3hour, presschangechar,
     &     reporttype, sealevelpress, skycover, skylayerbase,
     &     snowcover, stationname, tempfromtenths, temperature,
     &     timeobs, visibility, winddir, windgust, windspeed, wmoid,
     &     badflag, stanum, istatus )

      integer  maxskycover, recnum

      character*(*)  filename
      character(6)   autostationtype(recnum)
      character(25)  presweather(recnum)
      character(6)   reporttype(recnum)
      character(8)   skycover(maxskycover,recnum)
      character(5)   stationname(recnum)

      integer  presschangechar(recnum), wmoid(recnum)

      real  altimeter(recnum), dewpoint(recnum), dpfromtenths(recnum)
      real  elevation(recnum), latitude(recnum), longitude(recnum)
      real  maxtemp24hour(recnum), mintemp24hour(recnum)
      real  precip1hour(recnum), precip24hour(recnum)
      real  precip3hour(recnum), precip6hour(recnum)
      real  presschange3hour(recnum), sealevelpress(recnum)
      real  skylayerbase(maxskycover,recnum), snowcover(recnum)
      real  tempfromtenths(recnum), temperature(recnum)
      real  visibility(recnum), winddir(recnum), windgust(recnum)
      real  windspeed(recnum)
      real  badflag

      double precision  timeobs(recnum)

      character(3)   reportflag(recnum)
      character(2)   yy(recnum), mo(recnum), dd(recnum)
      character(2)   hh(recnum), mn(recnum)
      character(10)  time(recnum)
      character(9)   a10_to_a9
      character(8)   skycoverdmy(recnum)

      integer  windqua(recnum), sealevelpressqua(recnum)
      integer  temperaturequa(recnum), dewpointqua(recnum)
      integer  presschange3hourqua(recnum)
      integer  duplistation(9), stanum, dummy

c     real  tempdewdiff(recnum), lowestcloudheight(0:10)
c     real  skylayerbasedmy(recnum)
      real  skylayerbasedmy(recnum), tempdewdiff(recnum)
      real  lowestcloudheight(0:10)

      istatus= 0
      stanum= 0

      open ( 1, file=filename, status='old', err=1000 )

      istatus= 1

      data  duplistation / 58968, 58974, 59158, 59358, 59559,
     ~           	   59562, 59567, 59792, 59997 /
      data  lowestcloudheight / 0., 50., 100., 200., 300., 600., 1000.,
     ~                          1500., 2000., 2500., 0. /

c      ------   give initial values to avoid data stack problem  ------
      do 5 j= 1,recnum
      do 5 i= 1,maxskycover
         skycover(i,j)= "   "
5        skylayerbase(i,j)= badflag

      do j= 1,recnum
         read (1,20,end=99,err=9) reportflag(j), wmoid(j),
     ~                   elevation(j), latitude(j), longitude(j),
     ~                   yy(j), mo(j), dd(j), hh(j), mn(j)
         read (1,30,end=9,err=9) winddir(j), windspeed(j), windqua(j),
     ~                   visibility(j), presweather(j),
     ~                   sealevelpress(j), sealevelpressqua(j),
     ~                   temperature(j), temperaturequa(j),
     ~                   skycoverdmy(j), skylayerbasedmy(j)
         read (1,40,end=9,err=9) tempdewdiff(j), dewpointqua(j),       
     ~                   presschangechar(j), presschange3hour(j),
     ~                   presschange3hourqua(j), precip3hour(j),
     ~                   maxtemp24hour(j), mintemp24hour(j), windgust(j) 
         read (1,50,end=9,err=9) skycover(1,j), skylayerbase(1,j),
     ~                   skycover(2,j), skylayerbase(2,j), 
     ~                   skycover(3,j), skylayerbase(3,j), 
     ~                   precip24hour(j)
         read (1,60,end=9,err=9) skycover(4,j), skylayerbase(4,j)

         if ( reportflag(j) /= '*31' )  then
            write (6,*) 'read synop data heading error'
            go to 1000
         endif
	 go to 10

9        write (6,*) ' reading error of synop data'
         write (6,*) reportflag(j), wmoid(j), elevation(j),
     ~               latitude(j), longitude(j), ' ',
     ~               yy(j), mo(j), dd(j), hh(j), mn(j), ' ', timeobs(j)
         write (6,*) winddir(j), windspeed(j), windqua(j),
     ~               visibility(j), presweather(j),
     ~               sealevelpress(j), sealevelpressqua(j),
     ~               temperature(j), temperaturequa(j),
     ~               skycoverdmy(j), skylayerbasedmy(j)
         write (6,*) dewpoint(j), dewpointqua(j),
     ~               presschangechar(j), presschange3hour(j),
     ~               presschange3hourqua(j), precip3hour(j),
     ~               maxtemp24hour(j), mintemp24hour(j), windgust(j) 
         write (6,*) skycover(1,j), skylayerbase(1,j),
     ~               skycover(2,j), skylayerbase(2,j), 
     ~               skycover(3,j), skylayerbase(3,j), precip24hour(j)
         write (6,*) skycover(4,j), skylayerbase(4,j)

10       stanum= stanum +1
      enddo

20    format ( a3, i5, f4.0, 2f5.2, 2x, 5a2 )
30    format ( 2x, 2f3.0, i1, f3.0, a2, 3x, f5.1, i1, f4.1, i1, a2, 2x,
     ~         f3.0 )
40    format ( f4.1, i1, 1x, i2, f4.1, i1, 3(1x, f4.1), 8x, f3.0 )
50    format ( 3(a2, 2x, f3.0), 8x, f4.1 )
!hj: w>=d+3. 10/14/2013
60    format ( a2, 2x, f3.0 )

c     --- eliminate duplicate data coming from international broadcast ---
99    do 100 k= 1,9
      do 100 j= 1,stanum
         if ( wmoid(j) == duplistation(k) )  then
            stanum= stanum -1
            do i= j,stanum
               reportflag(i)= reportflag(i+1)
	       wmoid(i)     = wmoid(i+1)
	       elevation(i) = elevation(i+1)
	       latitude(i)  = latitude(i+1)
	       longitude(i) = longitude(i+1)
	       yy(i)= yy(i+1)
	       mo(i)= mo(i+1) 
	       dd(i)= dd(i+1)
	       hh(i)= hh(i+1)
	       mn(i)= mn(i+1)
	       winddir(i)            = winddir(i+1)
	       windspeed(i)          = windspeed(i+1)
	       windqua(i)            = windqua(i+1)
	       visibility(i)         = visibility(i+1)
	       presweather(i)        = presweather(i+1)
	       sealevelpress(i)      = sealevelpress(i+1)
	       sealevelpressqua(i)   = sealevelpressqua(i+1)
	       temperature(i)        = temperature(i+1)
	       temperaturequa(i)     = temperaturequa(i+1)
               skycoverdmy(i)        = skycoverdmy(i+1)
               skylayerbasedmy(i)    = skylayerbasedmy(i+1)
	       tempdewdiff(i)        = tempdewdiff(i+1)
	       dewpointqua(i)        = dewpointqua(i+1)
	       presschangechar(i)    = presschangechar(i+1)
	       presschange3hour(i)   = presschange3hour(i+1)
	       presschange3hourqua(i)= presschange3hourqua(i+1)
	       precip3hour(i)        = precip3hour(i+1)
	       maxtemp24hour(i)      = maxtemp24hour(i+1)
	       mintemp24hour(i)      = mintemp24hour(i+1)
	       windgust(i)           = windgust(i+1)
	       precip24hour(i)       = precip24hour(i+1)
               do l= 1,maxskycover
	          skycover(l,i)=     skycover(l,i+1)
		  skylayerbase(l,i)= skylayerbase(l,i+1)
	       enddo	
            enddo
         endif
100   continue
      write(*,*) 'synop stanum=', stanum

c      ----------       examine data quality and change units       ---------
      do j= 1,stanum
         if ( windqua(j) /= 1 )  then
            winddir(j)= badflag
            windspeed(j)= badflag
         endif

         if ( windgust(j) == -99. )  windgust(j)= badflag
         if ( elevation(j) == -999. )  elevation(j)= badflag
         if ( presschangechar(j)==-9 ) presschangechar(j)=int(badflag)
         if ( presweather(j) == '-9' )  presweather(j)= 'unk'

         if ( sealevelpressqua(j) == 1 )  then
               sealevelpress(j)= sealevelpress(j) *100.   ! millibar -> pascal
            else
               sealevelpress(j)= badflag
         endif

         if ( presschange3hourqua(j) == 1 )  then
            presschange3hour(j)= presschange3hour(j) *100. ! millibar -> pascal
          else
            presschange3hour(j)= badflag
            presschangechar(j)= int(badflag)
         endif

         if ( temperaturequa(j) == 1 )  then
               temperature(j)= temperature(j) +273.15           ! degc -> degk
            else
               temperature(j)= badflag
         endif

         tempfromtenths(j)= temperature(j)

         if ( dewpointqua(j) == 1 )  then
               dewpoint(j)= temperature(j) -tempdewdiff(j)        ! unit: degk
            else
               dewpoint(j)= badflag
         endif

         dpfromtenths(j)= dewpoint(j)

         if ( maxtemp24hour(j) == -99.9 )  then
               maxtemp24hour(j)= badflag
            else
               maxtemp24hour(j)= maxtemp24hour(j) +273.15       ! degc -> degk
         endif

         if ( mintemp24hour(j) == -99.9 )  then
               mintemp24hour(j)= badflag
            else
               mintemp24hour(j)= mintemp24hour(j) +273.15       ! degc -> degk
         endif

         if ( precip24hour(j) == -99.9 )  then
               precip24hour(j)= badflag
            else
               precip24hour(j)= precip24hour(j) *0.001   ! millimeter -> meter
         endif

         if ( precip3hour(j) == -99.9 )  then
               precip3hour(j)= badflag
            else
               precip3hour(j)= precip3hour(j) *0.001     ! millimeter -> meter
         endif

         if ( yy(j)(1:1) == ' ' )  yy(j)= '0'//yy(j)(2:2)
         if ( mo(j)(1:1) == ' ' )  mo(j)= '0'//mo(j)(2:2)
         if ( dd(j)(1:1) == ' ' )  dd(j)= '0'//dd(j)(2:2)
         if ( hh(j)(1:1) == ' ' )  hh(j)= '0'//hh(j)(2:2)
         if ( mn(j)(1:1) == ' ' )  mn(j)= '0'//mn(j)(2:2)
         time(j)= yy(j)//mo(j)//dd(j)//hh(j)//mn(j)
         call cv_asc_i4time( a10_to_a9(time(j),istatus), i4time )
         timeobs(j)= dble( i4time )                       ! seconds since 1960
      enddo

c    -------    transform code figure into visibility ( unit: m )  -------
      do j= 1,stanum
         if     ( visibility(j) ==  0. )  then
               visibility(j)= 50.
         elseif ( visibility(j) >  0.  .and. 
     ~            visibility(j) < 51. )  then
               visibility(j)= visibility(j) *100.
         elseif ( visibility(j) > 55.  .and. 
     ~            visibility(j) < 81. )  then
               visibility(j)= ( visibility(j) -50. ) *1000.
         elseif ( visibility(j) > 80.  .and. 
     ~            visibility(j) < 90. )  then
               visibility(j)= ( visibility(j) -74. ) *5000.
         elseif ( visibility(j) == 90. )  then
               visibility(j)= 25.
         elseif ( visibility(j) == 91. )  then
               visibility(j)= 50.
         elseif ( visibility(j) == 92. )  then
               visibility(j)= 200.
         elseif ( visibility(j) == 93. )  then
               visibility(j)= 500.
         elseif ( visibility(j) == 94. )  then
               visibility(j)= 1000.
         elseif ( visibility(j) == 95. )  then
               visibility(j)= 2000.
         elseif ( visibility(j) == 96. )  then
               visibility(j)= 4000.
         elseif ( visibility(j) == 97. )  then
               visibility(j)= 10000.
         elseif ( visibility(j) == 98. )  then
               visibility(j)= 20000.
         elseif ( visibility(j) == 99. )  then
               visibility(j)= 50000.
         else
               visibility(j)= badflag
         endif
      enddo

c ---- transform code figure into base of cloud layer indicated (unit: m) ----
      do 200 i= 1,maxskycover
      do 200 j= 1,stanum
         if     ( skylayerbase(i,j) ==  0. )  then
               skylayerbase(i,j)= 15.
         elseif ( skylayerbase(i,j) >  0.  .and.  
     ~            skylayerbase(i,j) < 51. )  then
               skylayerbase(i,j)= skylayerbase(i,j) *30.
         elseif ( skylayerbase(i,j) > 55.  .and.  
     ~            skylayerbase(i,j) < 81. )  then
               skylayerbase(i,j)= ( skylayerbase(i,j) -50. ) *300.
         elseif ( skylayerbase(i,j) > 80.  .and. 
     ~            skylayerbase(i,j) < 90. )  then
               skylayerbase(i,j)= ( skylayerbase(i,j) -74. ) *1500.
         elseif ( skylayerbase(i,j) == 90. )  then
               skylayerbase(i,j)= 25.
         elseif ( skylayerbase(i,j) == 91. )  then
               skylayerbase(i,j)= 75.
         elseif ( skylayerbase(i,j) == 92. )  then
               skylayerbase(i,j)= 150.
         elseif ( skylayerbase(i,j) == 93. )  then
               skylayerbase(i,j)= 250.
         elseif ( skylayerbase(i,j) == 94. )  then
               skylayerbase(i,j)= 450.
         elseif ( skylayerbase(i,j) == 95. )  then
               skylayerbase(i,j)= 800.
         elseif ( skylayerbase(i,j) == 96. )  then
               skylayerbase(i,j)= 1250.
         elseif ( skylayerbase(i,j) == 97. )  then
               skylayerbase(i,j)= 1750.
         elseif ( skylayerbase(i,j) == 98. )  then
               skylayerbase(i,j)= 2250.
         elseif ( skylayerbase(i,j) == 99. )  then
               skylayerbase(i,j)= 3000.
         else
               skylayerbase(i,j)= badflag
         endif
200   continue

      do 300 j= 1,stanum
         if (skycover(1,j) == '-9' .or. skylayerbase(1,j) == -9.)  cycle       

c          ----- eliminate duplicate skycovers and skylayerbases -----
         dummy= int( skylayerbasedmy(j) )
         skylayerbasedf= abs(skylayerbase(1,j)-lowestcloudheight(dummy))     
         if ( dummy >= 0 )  then
            if ( skycover(1,j) <= skycoverdmy(j)  .and.
     ~           (skylayerbase(1,j) >= lowestcloudheight(dummy) .or.
     ~            skylayerbasedf < 30.) )  then
               if ( dummy /= 9  .and.
     ~              skylayerbase(1,j) > lowestcloudheight(dummy+1) )
     ~            go to 250
               skycoverdmy(j)= '   '
	       skylayerbasedmy(j)= badflag
            endif
         endif

c        ----- eliminate unreasonable skycovers and skylayerbases -----
250      if ( skycover(1,j)==' 8' .and. skycoverdmy(j)==' 8' )  then
            skycoverdmy(j)= '   '
            skylayerbasedmy(j)= badflag
         endif
300   continue

c  -----  transform code figure into base of lowest cloud ( unit: m )  -----
      do j= 1,stanum
	 if     ( skylayerbasedmy(j) == 0. )  then
       	       skylayerbasedmy(j)= 25.
         elseif ( skylayerbasedmy(j) == 1. )  then
               skylayerbasedmy(j)= 75.
	 elseif ( skylayerbasedmy(j) == 2. )  then
	       skylayerbasedmy(j)= 150.
         elseif ( skylayerbasedmy(j) == 3. )  then
	       skylayerbasedmy(j)= 250.
	 elseif ( skylayerbasedmy(j) == 4. )  then
	       skylayerbasedmy(j)= 450.
	 elseif ( skylayerbasedmy(j) == 5. )  then
	       skylayerbasedmy(j)= 800.
	 elseif ( skylayerbasedmy(j) == 6. )  then
	       skylayerbasedmy(j)= 1250.
	 elseif ( skylayerbasedmy(j) == 7. )  then
	       skylayerbasedmy(j)= 1750.
	 elseif ( skylayerbasedmy(j) == 8. )  then
	       skylayerbasedmy(j)= 2250.
	 elseif ( skylayerbasedmy(j) == 9. .and. skycoverdmy(j) /= ' 0'
     ~             .and. skycoverdmy(j) /= '-9' )  then
	       skylayerbasedmy(j)= 3000.
	 else
	       skylayerbasedmy(j)= badflag
  	 endif
      enddo
 
c   assign the lowest cloud data to the first array when the latter is missing 
      do j= 1,stanum
	 if ( skylayerbasedmy(j) /= badflag .and. 
     ~        skylayerbase(1,j) == badflag )  then
            skylayerbase(1,j)= skylayerbasedmy(j)
            skycover(1,j)=     skycoverdmy(j)
         endif
      enddo
	    
c        --- transform code figure of cloud cover into metar format ---
      do 400 i= 1,maxskycover
      do 400 j= 1,stanum
         if     ( skycover(i,j) == ' 0' )  then
            skycover(i,j)= 'skc'
	    skylayerbase(i,j)= 22500.
         elseif ( skycover(i,j) == ' 1'  .or.
     ~            skycover(i,j) == ' 2' )  then
            skycover(i,j)= 'few'
         elseif ( skycover(i,j) == ' 3'  .or.
     ~            skycover(i,j) == ' 4' )  then
            skycover(i,j)= 'sct'
         elseif ( skycover(i,j) == ' 5'  .or.
     ~            skycover(i,j) == ' 6'  .or.
     ~            skycover(i,j) == ' 7' )  then
            skycover(i,j)= 'bkn'
         elseif ( skycover(i,j) == ' 8' )  then
            skycover(i,j)= 'ovc'
         else
            skycover(i,j)= '   '
         endif
400   continue

c               -------      deal with lacking of data      -------
      do j= 1,stanum
         autostationtype(j)= "unk"
c        presweather(j)= "unk"
         reporttype(j)= "synop"
         stationname(j)= "unk"

         altimeter(j)= badflag
         precip1hour(j)= badflag
         precip6hour(j)= badflag
         snowcover(j)= badflag
      enddo

1000  return
      end




      subroutine read_meso_cwb (inpath, maxobs, badflag, ibadflag,
     ~                          i4time_sys, timeobs, rpttp, stntp, 
     ~                          stname, lats, lons, elev,
     ~                          t, t24max, t24min, td, rh, 
     ~                          pcp1hr, pcp3hr, pcp6hr, pcp24hr, 
     ~                          dd, ff, wgdd, wgff,
     ~                          stnp, mslp, pcc, pc, sr, st,
     ~                          num, istatus)                    
 
c======================================================================
c
c     routine to read the cwb ascii mesonet files.
c     
c======================================================================
 
      real :: lats(maxobs), lons(maxobs), elev(maxobs)
      real :: t(maxobs), t24max(maxobs), t24min(maxobs), td(maxobs)
      real :: rh(maxobs), pcp1hr(maxobs), pcp3hr(maxobs), pcp6hr(maxobs)
      real :: pcp24hr(maxobs), dd(maxobs), ff(maxobs), wgdd(maxobs)
      real :: wgff(maxobs), stnp(maxobs), mslp(maxobs), pc(maxobs)
      real :: sr(maxobs), st(maxobs)
      integer :: pcc(maxobs), wmoid(maxobs)

      double precision  timeobs(maxobs)

      logical :: l_parse

c    larger arrays for istart and iend to read data to make processes smooth
      integer, parameter :: num40 = 40,  num70 = 70
      integer   :: istart(num70), iend(num70), hhmm, hh, flag
 
      character(*)  :: inpath
      character(13) :: cvt_i4time_wfo_fname13, a13time_eat
      character(6)  :: rpttp(maxobs), stntp(maxobs)
      character(5)  :: stname(maxobs), c5_blank
      character(3)  :: cstn_id, stn_id(maxobs)
      character     :: filename*80, line*320
 
c                      stuff for the mesonet metadata.
      real  lat_master(maxobs), lon_master(maxobs), elev_master(maxobs)
 
      character :: stn_id_master(maxobs)*3, stn_name_master(maxobs)*5
 
c               get the mesonet metadata (station information).
      call read_tmeso_stntbl (inpath, maxobs, badflag,  
     ~                        stn_id_master, stn_name_master,
     ~                        lat_master, lon_master, elev_master,
     ~                        num_master, istatus)
      if ( istatus /= 1 ) then
         write(6,*) ' error reading mesonet station table'
         return
      endif

c    fill the output arrays with something, then open the file to read.
 
      istatus=  0
      c5_blank= '     '
      rpttp=    'ldad'
      stntp=    'meso'
      stname=   c5_blank 
      pcc =     ibadflag
      t   =     badflag
      td  =     badflag
      rh  =     badflag
      stnp=     badflag
      dd  =     badflag
      ff  =     badflag
      wgdd=     badflag
      wgff=     badflag
      pc  =     badflag
      sr  =     badflag
      st  =     badflag
      timeobs=  dble( i4time_sys )

      i4time_file_eat= i4time_sys +8*3600             ! convert gmt to eat
      a13time_eat= cvt_i4time_wfo_fname13(i4time_file_eat)
c
c modified by min-ken,hsieh
c filename include hh information
c to read *_m.pri (to get each 15 min data) instead of *_h.pri
c because of the data format is different
c modified each field parsing below
c

      filename= 'data.cwb.mso.'
     ~           //a13time_eat(1:4)//'-'//a13time_eat(5:6)            ! yyyy_mm
     ~           //'-'//a13time_eat(7:11) //'00' //'_m.pri'    ! dd
      write(6,*) ' mesonet file ', filename

      call s_len ( inpath, len_inpath )
      call s_len ( filename, len_fname )
 
      num= 0
      num_keep= 0

      open (11,file=inpath(1:len_inpath)//filename(1:len_fname), 
     ~         status='old',err=980)

c.....  this starts the read loop.  since we don't know how many 
c.....  stations we have, read until we hit the end of file.
      
 100  flag= 0
 
      read (11,'(a)',end=600,err=990) line

c    find first dash in time portion (two spaces before last dash in string)
      do i= 1,300
         if (line(i:i) == ':')  exit
      enddo 
 
c                  parse the string into contiguous characters
      idash= i -9
      istart= 0
      iend=   0

      ivar= 1
      istart(1)= 1

      do i= 1,idash
         if ( i == 1 )  go to 200
         if ( line(i:i) == ' ' .and. line(i-1:i-1) /= ' ' ) then
            iend(ivar)= i-1
         endif

 200     if ( line(i:i) == ' ' .and. line(i+1:i+1) /= ' ' ) then
            ivar= ivar +1
            istart(ivar)= i+1
         endif
      enddo

      if ( istart(num40) /= 0 )  go to 100

c
c modified by min-ken,hsieh
c because m.pri file data format is different from h.pri
c each column means different data
c
      ivar= 1
      read (line(istart(ivar):iend(ivar)),'(i4)',err=399) ihrmin

      read (a13time_eat(10:13),'(i4)') hhmm
      read (a13time_eat(10:11),'(i2)') hh
      if ( ihrmin /= hhmm )  go to 100

      ivar= 2
      read (line(istart(ivar):iend(ivar)),*,err=399) cstn_id

      ivar= 9
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rstnp= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rstnp
      endif

      ivar= 13
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         slp= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) slp 
      endif

      ivar= 10
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rt= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rt
      endif

      ivar= 11
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rtd= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rtd
      endif

      ivar= 4
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         idir= ibadflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) idir
      endif

      ivar= 3
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rspd= badflag
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rspd
      endif

      ivar= 12
      if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
         rpcp= badflag
      elseif ( l_parse(line(istart(ivar):iend(ivar)),'000t') ) then
         rpcp= 0
      else
         read (line(istart(ivar):iend(ivar)),*,err=399) rpcp
      endif


      if ( num == 0 )  go to 400

      do i= 1,num
         if ( cstn_id//'  ' == stn_id(i) ) then 
            num= i
            flag= 1
            go to 500
         else
            cycle
         endif
      enddo
      go to 400

 399  write(6,*) ' read error in station/variable ', num+1, ivar
      write(6,*) ivar, line(istart(ivar):iend(ivar))
      go to 990

c    have good date/time...store ob.  adjust/scale variables while storing.
 400  num= num_keep +1                    ! add to count

      if ( num > maxobs ) then
         write(6,*) ' read_local_cwb error for too many obs: ',
     ~              num, maxobs
         istatus= 0
         return
      endif
 
c match data with metadata for this station, then store the metadata in arrays.
      imatch= 0
      do j= 1,num_master
         if ( cstn_id == stn_id_master(j) ) then
            lats(num)= lat_master(j)
            lons(num)= lon_master(j)
            elev(num)= elev_master(j)
            stn_id(num)= stn_id_master(j) 
            stname(num)= stn_name_master(j) 
            imatch=1
         endif
      enddo 

      if ( imatch == 0 ) then
         write(6,*) ' no station match ', cstn_id
      endif
 
c     stname(num)= stn_id//'  '
 
c                                quality control
 500  if ( rstnp <= 0 ) then
         stnp(num)= badflag
      else
         stnp(num)= rstnp *100.                        ! millibar -> pascal
      endif
 
      if ( slp <= 800. .or. slp > 1100. ) then
         mslp(num)= badflag
      else
         mslp(num)= slp *100.                          ! millibar -> pascal
      endif
 
      ! shuyuan & yuanfu: we found rpc and ipcc are un-initialized
      ! and are not read from the meso file. we remove the pcc and pc assignment

c      pcc(num)= ipcc 

c      if ( rpc <= 0 ) then
c         pc(num)= badflag
c      else
c         pc(num)= rpc *100.                            ! millibar -> pascal 
c      endif
 
      if ( rt <= -90 ) then
         t(num)= badflag
      else
         if ( rt > 50. )  rt= - (rt - 50.)
         t(num)= rt +273.15                            ! degc -> degk
      endif
 
      if ( rtd <= -90 ) then
         td(num)= badflag
      else
         if ( rtd > 50. ) rtd= - (rtd - 50.)
         td(num)= rtd +273.15                          ! degc -> degk
      endif
 
      if ( idir > 36 .or. idir < 0 ) then
         dd(num)= badflag
      else
         dd(num)= float(idir * 10)                     ! unit : deg
      endif
 
      if ( rspd < 0 ) then
         ff(num)= badflag
      else
         ff(num)= rspd                                 ! unit : m/s 
      endif
 
      ! shuyuan & yuanfu: rwgff/wgdd are never initialized (gust wind)
c      if ( rwgff < 0 ) then
c         wgff(num)= badflag
c      else
c         wgff(num)= rwgff                              ! unit : m/s
c      endif
 
c      if ( iwgdd > 36 .or. iwgdd < 0 ) then
c         wgdd(num)= badflag
c      else
c         wgdd(num)= float(iwgdd * 10)                  ! unit : deg
c      endif
 
      if ( rpcp < 0 ) then
         pcp1hr(num)= badflag
      else
         pcp1hr(num)= rpcp *0.001                      ! millimeter -> meter
      endif
 
      ! shuyuan & yuanfu: rsr is never initialized (gust wind)
      if ( rsr < 0 ) then
         sr(num)= badflag
      else
         sr(num)= rsr /1000. /3600.                    ! conv mj/m/m to watt/m/m
      endif
 
      ! shuyuan & yuanfu: we found another un-initialized variable - irh
      ! we have to remove the following for now:
c      if ( irh < 0 ) then
c         rh(num)= badflag
c      else
c         rh(num)= float(irh)                           ! unit : %
c      endif
 
      ! shuyuan & yuanfu: rst is never initialized variable
c      if ( rst < 0 ) then
c         st(num)= badflag
c      else
c         st(num)= rst                                  ! degc -> degk
c      endif
 
c                          go back for the next ob.
      if ( flag == 0 )  num_keep= num 
      num= num_keep
      go to 100
 
 600  call mso_t24_pcp (inpath, filename, stname, stn_id, maxobs, 
     ~     badflag, hh, num, t24max, t24min, pcp3hr, pcp6hr, pcp24hr,      
     ~     istatus)

      if ( istatus == 1 ) then
c
c        modified by min-ken hsieh
c        some stn may not have enough info to calculate t24 and pcp accum.
c        and mso_t24_pcp will return badsfc here
c        so let's check return values before unit conversion
c
         do i= 1,maxobs
	    if (t24max(i).ne.badflag) then
               t24max(i)= t24max(i) +273.15                  ! degc -> degk
	    endif
	    if (t24min(i).ne.badflag) then
               t24min(i)= t24min(i) +273.15                  ! degc -> degk
	    endif
	    if (pcp3hr(i).ne.badflag) then
               pcp3hr(i)= pcp3hr(i) *0.001                   ! millimeter -> meter
	    endif
	    if (pcp6hr(i).ne.badflag) then
               pcp6hr(i)= pcp6hr(i) *0.001                   ! millimeter -> meter
	    endif
	    if (pcp24hr(i).ne.badflag) then
               pcp24hr(i)= pcp24hr(i) *0.001                 ! millimeter -> meter
	    endif
         enddo
      else
         write(6,*) ' error estimating mso_t24_pcp '
      endif

c                        hit end of file...that's it.
      write(6,*) ' found ', num, ' mesonet stations.'
      istatus= 1
      return
      
 980  write(6,*) ' warning: could not open mesonet data file ',filename
      num= 0
      istatus= -1
      return

 990  write(6,*) ' ** error reading mesonet data.'
      num= 0
      istatus= -1
      return
      
      end
 
 
 
      subroutine read_tmeso_stntbl (inpath, maxobs, badflag, stn_id,
     ~           stn_name, lat, lon, elev, num, istatus)       
 
c======================================================================
c
c     routine to read station information for the cwb ascii mesonet 
c	data.
c     
c======================================================================
 
      real         :: lat(maxobs), lon(maxobs), elev(maxobs)
      character(3) :: stn_id(maxobs), stn_id_in
      character(5) :: stn_name(maxobs), stn_name_in
      character(*) :: inpath
 
      lat=  badflag
      lon=  badflag
      elev= badflag
      stn_id=   '   '
      stn_name= '     '
 
      call s_len ( inpath, len_inpath )
      open (13,file=inpath(1:len_inpath)//'stn-table',status='old',
     ~                                                err=990)
      num= 0

c                 skip header comments at the top of the file
      do iread= 1,2
         read (13,*,end=550,err=990)
      enddo
 
c.....  this starts the station read loop.  since we don't know how many 
c.....  stations we have, read until we hit the end of file.

 500  read (13,900,end=550,err=990) stn_id_in,stn_name_in,
     ~                              lat_deg,lat_min,lat_sec,alat_sec,       
     ~                              lon_deg,lon_min,lon_sec,alon_sec,
     ~                              elev_m
 900  format (2x,a3,1x,a5,14x,                     ! name
     ~        i2,2x,i2,1x,i2,1x,f3.0,4x,           ! lat
     ~        i3,2x,i2,1x,i2,1x,f3.0,              ! lon
     ~        f12.0)                               ! elevation
 
c         move station info to arrays for sending to calling routine.
 
      alat= float(lat_deg) +float(lat_min)/60. 
     ~                     +(float(lat_sec) +alat_sec) /3600.
      alon= float(lon_deg) +float(lon_min)/60. 
     ~                     +(float(lon_sec) +alon_sec) /3600.

      num= num +1
      stn_id(num)= stn_id_in
      stn_name(num)= stn_name_in
      lat(num)= alat
      lon(num)= alon
      elev(num)= elev_m
 
c                         go back for the next ob.
      go to 500
 
c                        hit end of file...that's it.
 550  write(6,*) ' found ', num,
     ~           ' mesonet stations in the station table.' 
      istatus= 1
      return
      
 980  write(6,*) ' warning: could not open mesonet station file ',
     ~           inpath
      istatus= -1
      return

 990  write(6,*) stn_id_in, stn_name_in,
     ~           lat_deg, lat_min, lat_sec, alat_sec,       
     ~           lon_deg, lon_min, lon_sec, alon_sec, elev_m
      write(6,*) ' ** error reading mesonet station table'
      istatus= 0
      return

      end


 
      subroutine mso_t24_pcp (inpath, filename, stname, stn_id, maxobs,
     ~           badflag, ih, num, t24max, t24min, pcp3hr, pcp6hr, 
     ~           pcp24hr, istatus)

      integer, parameter :: num24 = 24,  num40 = 40,  num70 = 70
      character(*) :: stname(maxobs), stn_id(maxobs), inpath
      character(2) :: yy, mm, dd
      character    :: filename*35, filedummy*35, line*320, stn*3
      logical :: l_parse
      integer :: istart(num70), iend(num70), d(12) 
      integer :: hr, flag
      real :: t24max(maxobs), t24min(maxobs)
      real :: tmax(maxobs,num24), tmin(maxobs,num24)
      real :: pcp3hr(maxobs), pcp6hr(maxobs), pcp24hr(maxobs)
      real :: p1hr(maxobs,-4:num24)
      data  d / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

      stn=  '   '
      tmax= badflag
      tmin= badflag
      p1hr=    0
      pcp3hr=  0
      pcp6hr=  0
      pcp24hr= 0
      maxrcd= 3000
      istatus= 0

      read (filename(14:17),'(i4)') iy
      read (filename(19:20),'(i2)') im
      read (filename(22:23),'(i2)') id

c                 open two files to read data of 24 hours
      do l= 0,1
         id= id -l

         if ( id < 1 ) then
            im= im -1
  
            if ( im < 1 ) then
               im= 12
               iy= iy -1
            endif

            id= d(im)
         endif

         iy= iy -2000
         call i2a ( iy, yy )
         call i2a ( im, mm )
         call i2a ( id, dd )
         filedummy= 'data.cwb.mso.' //'20' //yy //'-' //mm //'-' //dd
     ~                                          //'_' //'0000_h.pri'     

         iy= iy +2000

         call s_len ( inpath, len_inpath )
         call s_len ( filedummy, len_fname )

         select case ( l ) 
         case ( 0 )
c
c modified by min-ken hsieh
c we read _m file to get other data in subroutine read_mso_cwb
c but here we need to get pcp data from _h file
c

c           rewind (11)
            open (11,file=inpath(1:len_inpath)//filedummy(1:len_fname),
     ~               status='old',err=980)
         case ( 1 )
            open (11,file=inpath(1:len_inpath)//filedummy(1:len_fname),
     ~               status='old',err=980)
         end select

         do 200 k= 1,maxrcd
            read (11,'(a)',iostat=istat) line

            if ( istat == -1 .or. istat == -2 )  exit
            if ( istat > 0 )  go to 990
               
c    find first dash in time portion (two spaces before last dash in string)
            do i= 1,300
               if (line(i:i) == ':')  exit
            enddo

            idash= i -9

c                 parse the string into contiguous characters
            istart= 0
            iend=   0

            ivar= 1
            istart(1)= 1

            do i= 1,idash
               if ( i == 1 )  go to 100

               if ( line(i:i) == ' ' .and. line(i-1:i-1) /= ' ' ) then      
                  iend(ivar)= i-1
               endif

 100           if ( line(i:i) == ' ' .and. line(i+1:i+1) /= ' ' ) then
                  ivar= ivar +1
                  istart(ivar)= i+1
               endif
            enddo

c                              avoid daily data
            if ( istart(num40) /= 0 )  cycle

            ivar= 1
            read (line(istart(ivar):iend(ivar)),'(2i2)',err=199) hr, mn
            if ( l == 0 .and. hr >  ih )  cycle  
            if ( l == 1 .and. hr <= ih )  cycle  
            if ( hr == 0 )  hr= 24
           
            ivar= 2
            read (line(istart(ivar):iend(ivar)),*,err=200) stn

            do i= 1,num
            if ( stn_id(i) == stn//'  ' ) then
               ivar= 8
               if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
                 tmax(i,hr)= badflag
               else
                 read(line(istart(ivar):iend(ivar)),*,err=199)tmax(i,hr)
               endif

               ivar= 9
               if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
                 tmin(i,hr)= badflag
               else
                 read(line(istart(ivar):iend(ivar)),*,err=199)tmin(i,hr)    
               endif

               ivar= 15
               if ( l_parse(line(istart(ivar):iend(ivar)),'/') ) then
                 p1hr(i,hr)= badflag
               elseif(l_parse(line(istart(ivar):iend(ivar)),'000t'))then
                 p1hr(i,hr)= 0
               else
                 read(line(istart(ivar):iend(ivar)),*,err=199)p1hr(i,hr)
               endif
            endif
            enddo

            cycle

 199        write(6,*) ' read error in station/variable ', j, ivar
            write(6,*) ivar,line(istart(ivar):iend(ivar))
 200     enddo

         if ( ih == 23 )  exit 
 500  enddo 

      do i= 1,num
         t24max(i)= tmax(i,1)
         t24min(i)= tmin(i,1)

         do j= 2,num24
            flag= 0

            if ( tmax(i,j) > 50. .or. tmax(i,j) < -90. ) then
               tmax(i,j)= badflag
               exit
            elseif ( tmin(i,j) > 50. .or. tmin(i,j) < -90. ) then
               tmin(i,j)= badflag
               exit
            endif

            if ( t24max(i) < tmax(i,j) )  t24max(i)= tmax(i,j)
            if ( t24min(i) > tmin(i,j) )  t24min(i)= tmin(i,j)

            flag= 1
         enddo

c  once there is any data missing in tmax or tmin, flag= 0 from prior do loop
         if ( flag /= 1 ) then
            t24max(i)= badflag
            t24min(i)= badflag
            write(6,*) 'too few data to obtain tmax/tmin for ',
     ~                 stname(i), ' mesonet station ', j
            cycle
         endif
 900  enddo

c                         calculate accumulated rain gauge
      do i= 1,num
         p1hr(i, 0)= p1hr(i,24)
         p1hr(i,-1)= p1hr(i,23)
         p1hr(i,-2)= p1hr(i,22)
         p1hr(i,-3)= p1hr(i,21)
         p1hr(i,-4)= p1hr(i,20)
      enddo

      do i= 1,num
      do j= ih,ih-2,-1
         if ( p1hr(i,j) == badflag ) then
            pcp3hr(i)= badflag
            write(6,*) 'too few data to estimate pcp3hr for ',
     ~                 stname(i), ' mesonet station ', j
            exit
         else
            pcp3hr(i)= pcp3hr(i) +p1hr(i,j)
         endif
      enddo
      enddo

      do i= 1,num
      do j= ih,ih-5,-1
         if ( p1hr(i,j) == badflag ) then
            pcp6hr(i)= badflag
            write(6,*) 'too few data to estimate pcp6hr for ',
     ~                 stname(i), ' mesonet station ', j
            exit
         else
            pcp6hr(i)= pcp6hr(i) +p1hr(i,j)
         endif
      enddo
      enddo

      do i= 1,num
      do j= 1,num24
         if ( p1hr(i,j) == badflag ) then
            pcp24hr(i)= badflag
            write(6,*) 'too few data to estimate pcp24hr for ',
     ~                 stname(i), ' mesonet station ', j
            exit
         else
         endif
      enddo
      enddo

      istatus= 1
      return

 980  write(6,*) ' warning: could not open mesonet data file ',
     ~           inpath(1:len_inpath)//filedummy(1:len_fname)
      istatus= -1
      return
        
 990  write(6,*) ' ** error reading mesonet data.'
      num= 0
      istatus= -1
      return

      end
        


      subroutine  i2a (ii,aa)

      character(2) :: aa
      integer      :: ii

      if ( ii < 10 ) then
         write (aa,'(a1,i1)') '0', ii
      else
         write (aa,'(i2)') ii
      endif

      return
      end
