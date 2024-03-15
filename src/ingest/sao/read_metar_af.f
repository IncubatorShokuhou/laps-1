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
c
c
      subroutine read_sao_af(filetime,path_to_obs,n,sname,
     +                   lats,lons,elevs,
     +                   time,type,atype,
     +                   cover,ht,
     +                   vis,weather,
     +                   mslp,t,td,dd,ff,ffg,alt,
     +                   dpch,dp,nmax,badflag,istatus)
c
c-------------------------------------------------------------------------
c
c     routine to read the air force surface data files.
c
c     changes:  p. stamus  06-24-96  original from fsl version.
c                          03-19-97  changed 1st format statement.
c                          10-28-97  changes for dynamic laps.
c
c-------------------------------------------------------------------------
c
c
      real lats(nmax), lons(nmax), elevs(nmax)
      real mslp(nmax), alt(nmax)
      real t(nmax), tdd, td(nmax)
      real dd(nmax), ff(nmax), ffg(nmax)
      real vis(nmax), time(nmax)
      real ht(5,nmax)
c
      integer timeo, timen
      integer dpch(nmax), dp(nmax)
      integer cvr, base, meth
      integer rtype, astype
      integer wx1, wx2, wx3, wx4, wx5, wx6, wx7
      integer cld1a, cld1h, cld2a, cld2h, cld3a, cld3h
      integer cld4a, cld4h, cld5a, cld5h
c
      character filetime*9, af_time*8, a9_to_a8*8, af_file*256
      character id*5, dum*132
      character name_in*5, sname(nmax)*5, type(nmax)*6, atype(nmax)*6
      character weather(nmax)*40
      character cover(5,nmax)*8, amt_out*3 
      character*(*) path_to_obs
c
c..... set number of header lines to skip at top of file.
c
      num_hdr = 16
	istatus = -1
c
c..... zero out the cloud arrays and other stuff.
c
      do j=1,nmax
         do i=1,5
            cover(i,j) = '        '
            ht(i,j) = badflag
         enddo !i
         weather(j) = '                                        '
         type(j) = '      '
         sname(j) = '     '
         atype(j)(1:6) = 'unk   '
      enddo !j
c
c..... open the data file.
c
      af_time = a9_to_a8(filetime)
      len = index(path_to_obs,' ') - 1
      af_file = path_to_obs(1:len) // 'sao.' // af_time
cc      af_file = '../lapsdat/sao_af_raw2/sao.' // af_time
cc      af_file = 'sao.' // af_time
cc    af_file = '../../../lapsdat/sao_af_raw/sao.'//af_time
c
      open(11,file=af_file,status='old',err=590)
c
c..... read the file header lines.
c
      do i=1,num_hdr
         read(11,900,err=600,end=600) dum
      enddo !i
 900  format(a132)
c
c..... now read and count up each station, adjusting units/codes/etc
c..... to something laps expects.
c
      n = 0
 500  n = n + 1
c
      read(11,910,end=400) id,name_in,lats(n),lons(n),
     &             elevs(n),timeo,timen,rtype,astype,
     &             mslp(n)
 910  format(a5,2x,a5,2x,f10.2,2x,f10.2,2x,f9.0,2x,i4,5x,i4,2x,
     &       i9,2x,i9,1x,f11.1)
c
cc     write(6,910) id,name_in,lats(n),lons(n),
cc    &             elevs(n),timeo,timen,rtype,astype,
cc    &             mslp(n)
c
      read(11,920) alt(n),dpch(n),dp(n),t(n),tdd,dd(n),ff(n),
     &             ffg(n),vis(n),wx1,wx2,wx3,wx4,wx5,wx6,wx7
 920  format(2x,f11.2,2x,i9,2x,i9,2x,f10.1,2x,f11.1,3x,f8.0,3x,
     &       f9.1,3x,f9.1,3x,f9.3,2x,7i2)

cc     write(6,920) alt(n),dpch(n),dp(n),t(n),tdd,dd(n),ff(n),
cc    &             ffg(n),vis(n),wx1,wx2,wx3,wx4,wx5,wx6,wx7


c
      read(11,930) cvr,base,meth,cld1a,cld1h,cld2a,cld2h,
     &             cld3a,cld3h,cld4a
 930  format(3x,i4,3x,i9,3x,i8,7(2x,i9))
c
cc     write(6,930) cvr,base,meth,cld1a,cld1h,cld2a,cld2h,
cc    &             cld3a,cld3h,cld4a


      read(11,940) cld4h,cld5a,cld5h
 940  format(3x,i9,2x,i9,2x,i9)

cc     write(6,940) cld4h,cld5a,cld5h
c
c
c.....	check for bad station info
c
	if(name_in(3:5) .eq. '???') then
	   n = n - 1
	   go to 500
	endif
c
	if(lats(n).lt.-90.  .or. lats(n).gt.90. .or. 
     &	   lons(n).lt.-360. .or. lons(n).gt.360. .or.
     &     elevs(n).lt.-500. .or. elevs(n).gt.5000.) then
	   n = n - 1
	   go to 500
	endif
c
c..... set up final arrays and convert units as needed.
c
      if(name_in .eq. '    ') then
         sname(n)(1:5) = id(1:5)
      else
         sname(n)(1:5) = name_in(1:5)
      endif
c
c..... figure out time
c
      if(timen .gt. 24) then
         n = n - 1
         go to 500
      endif
c
      if(timeo .gt. 60) timeo = 0
c
      if(timeo .ge. 45) then
         if(timen .eq. 00) then
            timen = 23
         else
            timen = timen - 1
         endif
      endif
c
      time(n) = float(timen*100 + timeo)
c
c..... figure out report type
c
      if(rtype .eq. 2.and.astype.eq.2) type(n) = 'metar'

      if(rtype .eq. 3.and.astype.eq.2) type(n) = 'speci'

      if(rtype .eq. 2.and.astype.eq.4) type(n) = 'metar'

      if(rtype .eq. 2.and.astype.eq.9) type(n) = 'buoy'

      if(rtype .eq. 2.and.astype.eq.0) type(n) = 'synop'
c
c..... auto station?  
c
      if(astype .eq. 9) then
         atype(n)(1:4) = 'auto'
      endif
c
c..... calculate dewpoint, do some first checks, change some units.
c..... winds, vis, cld hts, taken care of in calling routine.
c
      td(n) = badflag
      if(t(n).gt.0. .and. t(n).lt.400.) then
         if(tdd.gt.0. .and. tdd.lt.100.) then
            td(n) = t(n) - tdd
         endif
      else
         t(n) = badflag
      endif
c         
      if(alt(n).gt.0. .and. alt(n).lt.100.) then
         alt(n) = alt(n) * 33.86
      else
         alt(n) = badflag
      endif
c
      if(dpch(n).lt.0 .or. dpch(n).gt.9) then !gives sign
         dpch(n) = -99
      endif
      if(dp(n).lt.0 .or. dpch(n).gt.200) then !in tenths mb
         dp(n) = -99
      endif
c
c.....  put cloud data into arrays ht(5,x), cover(5,x) 
c
      call get_cld_amt_af(cld1a, amt_out)
      call get_cld_hts_af(cld1h, ht_out)
      cover(1,n)(1:3) = amt_out
      ht(1,n) = ht_out
c
      call get_cld_amt_af(cld2a, amt_out)
      call get_cld_hts_af(cld2h, ht_out)
      cover(2,n)(1:3) = amt_out
      ht(2,n) = ht_out
c
      call get_cld_amt_af(cld3a, amt_out)
      call get_cld_hts_af(cld3h, ht_out)
      cover(3,n)(1:3) = amt_out
      ht(3,n) = ht_out
c
      call get_cld_amt_af(cld4a, amt_out)
      call get_cld_hts_af(cld4h, ht_out)
      cover(4,n)(1:3) = amt_out
      ht(4,n) = ht_out
c
      call get_cld_amt_af(cld5a, amt_out)
      call get_cld_hts_af(cld5h, ht_out)
      cover(5,n)(1:3) = amt_out
      ht(5,n) = ht_out
c
c..... back to top of loop.
c
      go to 500
c
c..... when we're done...
c
 400  continue
      n = n - 1
      write(6,955) n
 955  format(' found ',i4,' stations in surface file.')
      istatus = 1
      return
c
 590  	continue
	print *,' warning: could not open sao file.'
	istatus = 0
	return
c
 600  	continue
	print *,' error reading sao file.'
	istatus = 0
	return
c
      end
c
c
      subroutine get_cld_amt_af(amt_in, amt_out)
c
c..... routine to figure out the af cloud amounts.
c
      integer amt_in
      character amt_out*3
c
      amt_in = abs( amt_in )
      amt_out = '   '
      if(amt_in .gt. 9) return
      if(amt_in .eq. 9) amt_out = 'x  '
      if(amt_in.ge.8 .and. amt_in.lt.9) amt_out = 'ovc'
      if(amt_in.ge.4 .and. amt_in.lt.8) amt_out = 'bkn'
      if(amt_in.ge.1 .and. amt_in.lt.4) amt_out = 'sct'
      if(amt_in.ge.0 .and. amt_in.lt.1) amt_out = 'clr'
      if(amt_in .lt. 0) return
c
      return
      end
c
c    
      subroutine get_cld_hts_af(ht_in, ht_out)
c
c..... routine to figure out the af cloud heights.
c
      integer ht_in
      real ht_out
c
      ht_out = -99.9
      if(ht_in .lt.  0) return
      if(ht_in .le. 50) then
         ht_out = float(ht_in) * 30.
         return
      elseif(ht_in .le. 80) then
         ht_out = float(ht_in - 50) * 300.
      elseif(ht_in .lt. 90) then
         ht_out = (float(ht_in - 80) * 1500.) + 9000.
      endif
c
      return
      end
