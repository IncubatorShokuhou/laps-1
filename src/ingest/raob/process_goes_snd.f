
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


      subroutine process_goes_snd (path, path_len, filename, 
     1     file_len,iii,jjj,
     1     i4time_begin,i4time_end, mdf, lun_out, istatus)

c
c	modified for filename change 4/22/2013 db
c
      implicit none

c     parameter list
      integer istatus
      character*(*) path, filename
      real mdf                  !missing data flag
      integer path_len, file_len, i4time_begin,i4time_end,lun_out
      integer iii,jjj

c     variables for output writing
      integer maxsnd,maxlvl
      parameter (maxsnd = 1)
      parameter (maxlvl = 50)

      integer iwmostanum(maxsnd), nlvl (maxsnd)
      real stalat(maxsnd,maxlvl),stalon(maxsnd,maxlvl),staelev(maxsnd)
      character c5_staid(maxsnd)*5,a9time_ob(maxsnd,maxlvl)*9,
     1     c8_obstype(maxsnd)*8
      real height_m(maxsnd,maxlvl)
      real pressure_mb(maxsnd,maxlvl)
      real temp_c(maxsnd,maxlvl)
      real dewpoint_c(maxsnd,maxlvl)
      real dir_deg(maxsnd,maxlvl)
      real spd_mps(maxsnd,maxlvl)
      integer nsnd ! number of soundings per write call (1)
      integer count_level,ob_counter,id_sat



c     internal variables
      real lat,lon
      integer i4time_record,nbuf_1,nbuf_2
      integer avg_i4time_record ! to create modified time rec
      character*9 c_time_record
      character*9 m_time_record ! modified time record
      character*11 c_time_record_long
      real lat_a(iii,jjj), lon_a(iii,jjj), topo_a (iii,jjj)
      real rnorth,south,east,west

c ?
c ?
c ?  *****  data record format  *****
c ?
c ?    -----  header record (record 1)  -----
c ?
c ?   word       description                         scale
c ?  -----       -----------                         -----
c ?    1         day (yyddd)                         1
c ?    2         time (hhmmss)                       1
c ?    3         number of records in file           1
c ?    4         satid                               1
c ?    5         lat/lon of nw corner of data        1
c ?    6         lat/lon of se corner of data        1
c ?
c ?   ------   data (records 2-n)  ------
c ?   
c ?   word       description                         scale
c ?  -----       -----------                         -----
c ?    1         day (yyddd)                         1
c ?    2         time (hhmmss)                       1
c ?    3         satid                               1
c ?    4         mod flag                            1
c ?    5         latitude (deg)                      100 
c ?    6         longitude (deg)                     100
c ?    7         sam (#fovs)                         1
c ?    8         cape index                          1
c ?    9-26      observed bt (k)                     100
c ?    27        solar zenith angle (deg)            100
c ?    28        local zenith angle (deg)            100
c ?    29        retrieval type                      1
c ?    30        skin temp (k)                       100
c ?    31        total precip water (mm)             100
c ?    32        lifted index (k)                    100
c ?    33        guess total precip water (k)        100
c ?    34        guess lifted index (k)              100
c ?    35        layer precip. water (1-.9 sigma)    100
c ?    36        layer precip. water (.9-.7 sigma)   100
c ?    37        layer precip. water (.7-.3 sigma)   100
c ?    38        sfc temp (k) @ sfc pressure         100
c ?    39-78     temp profile (k) (1000-.1mb)        100
c ?    79        sfc dewpoint temp (k) @ sfc press   100
c ?    80-119    dewpoint temp (k) (1000-0.1mb)      100
c ?    120       sfc geopotential height (m)         1
c ?    121-160   geopotential height (m)             1
c ?    161       sfc presuure (mb)                   10
c ?    162-201   pressure (1000-0.1mb)               10
c ?    202-210   spares                              1
c ?
c ?   ###########################################################

      character*1 csat
      character*5 cyyjd
      character*6 cvar(37)
      character*4 chhmm
      character*80 ctext
      integer*4 nbuf(210),iscale(210)
      real*4 rbuf(210)

      data iscale /4*1,2*100,2*1,20*100,1,90*100,41*1,41*10,9*1/
      data cvar /'day   ','time  ','satid ','mod   ','lat   ',
     1'lon   ','sam   ','ca    ','bt01  ','bt02  ','bt03  ','bt04  ',
     2'bt05  ','bt06  ','bt07  ','bt08  ','bt09  ','bt10  ','bt11  ',
     3'bt12  ','bt13  ','bt14  ','bt15  ','bt16  ','bt17  ','bt18  ',
     4'sza   ','lza   ','rt    ','tskin ','tpw   ','li    ','gtpw  ',
     5'gli   ','wv1   ','wv2   ','wv3   '/

      equivalence(nbuf(1),rbuf(1))

      integer ilev,jj,js,k,lt,lp,lz,ltd,irec,iend,ibegin,nrec,mm,lll
      integer ilin


c-------------------------------------------------------------------------------------
c

      istatus = 0 ! set initially as bad

c     compute average record time for ingest
      avg_i4time_record = nint((float(i4time_begin)+
     1     float(i4time_end))/2.)
      call make_fnam_lp(avg_i4time_record,m_time_record,istatus)

c     get perimeter values for the domain
      call get_domain_perimeter(iii,jjj,'nest7grid',lat_a,
     1     lon_a, topo_a, 1.0, rnorth,south,east,west,istatus)

c     open file 


       open(70,file=path(1:path_len)//filename(1:file_len)
     1     ,access='direct',recl=840,
     1 form='unformatted',status='old',err=9000)
c
      irec = 0
c     read header information

      read(70,rec=1,err=9999) nbuf

c     endian swap if needed
      if (filename(5:6) .lt. '09') then
      do k  = 1, 210
         call endian4 (nbuf(k))
      enddo
      endif


      nrec = nbuf(1)
      write(6,*) 'header record contents:'
      write(6,*) '-----------------------'
      write(6,*) 'day  of data (yyddd)     = ',nbuf(1)
      write(6,*) 'time of data (hhmm)      = ',nbuf(2),'z'
      write(6,*) 'number of records        = ',nbuf(3)
      write(6,*) 'satellite id             = ',nbuf(4)
      write(6,*) 'lat/lon of nw corner     = ',nbuf(5)
      write(6,*) 'lat/lon of se corner     = ',nbuf(6)
      write(6,*)

      id_sat = nbuf(4)

      ibegin = 2
      iend = nbuf(3)  ! number of records

      do 200 ilin=ibegin,iend      
        do 205 lll=1,210
           nbuf(lll) = 0
  205   continue
        read(70,rec=ilin,err=9999) nbuf 

c     endian swap if needed
      if (filename(5:6) .lt. '09') then
        do k = 1, 210
           call endian4 (nbuf(k))
        enddo
      endif


        do 210 mm = 1,210
           rbuf(mm) = (nbuf(mm)/float(iscale(mm)))
  210   continue
        irec = irec + 1
c        write(6,1000) ilin 
 1000   format(20x,'***  record ',i6,' contents  ***' /)
c        write(6,1010) (cvar(k),rbuf(k),k=1,4) 
 1010   format(1x,4(a6,' = ',f7.0,1x,'|',2x),/)
c        write(6,1015) (cvar(k),rbuf(k),k=5,7) 
 1015   format(1x,a6,' = ',f6.2,2x,'|',2x,a6,' = ',
     1  f6.2,2x,'|',2x,a6,' = ',f6.0,2x,'|',/)
c        write(6,1020) (cvar(k),rbuf(k),k=8,25)
 1020   format(4(1x,a6,' = ',f8.2,2x,'|',1x),/)
c        write(6,1030) (cvar(k),rbuf(k),k=26,37)
 1030   format(/,4(1x,a6,' = ',f6.2,2x,'|',1x),/)
c        write (6,*) rbuf(37)+rbuf(36)+rbuf(35)
c        write(6,1040) 
 1040   format(/,1x,'pressure',4x,'   t    ',4x,'   td   ',4x,
     1  '   z    ')
c        write(6,1050)
 1050   format(1x,'--------',4x,'--------',4x,'--------',4x,'--------')


c     set up time window compare
        nbuf_1=int(rbuf(1))
        nbuf_2=int(rbuf(2))
        write (c_time_record_long,34) nbuf_1,nbuf_2
        c_time_record = c_time_record_long(1:9)
 34     format (i5.5,i6.6)
        call i4time_fname_lp(c_time_record,i4time_record,istatus)
        if (i4time_record.ge.i4time_begin 
     1       .and. 
     1       i4time_record.le.i4time_end) then
           continue
        else
           write (6,*) 'skipping record time bounds ',c_time_record
           go to 200 ! bail on this record
        endif


c     extract lat and lon

        lat = rbuf(5)
        lon = -1.*rbuf(6) ! our convention is east longitude

c     test whether point is in the domain

        if (lat.le.rnorth .and. lat .ge. south
     1       .and. lon .ge. west .and. lon .le. east )then
           write (6,*) 'accepting sounding ', lat,lon
           write (6,*) 'routine process_goes_snd'
           continue
        else
           ! outside the perimeter - reject
           write (6,*) 'lat lon reject', lat, lon
           go to 200            !bail on this record

        endif                   !domain test finished

c     good ob
        ob_counter=ob_counter+1

        count_level = 0
        js = 1
        ilev = 0
        write (6,*) lat,lon
        do 170 jj=38,78
          ilev = ilev + 1
          lt   = jj
          ltd  = lt + 41
          lz   = ltd + 41
          lp   = lz + 41
          if(rbuf(lt) .gt. -9.9e+3) then
c            write(6,1060) rbuf(lp),rbuf(lt)-273.15,
c     1            rbuf(ltd)-273.15  ,rbuf(lz)
            count_level=count_level+1

c     fill variables for calling the routine to write snd
            nsnd = 1
            iwmostanum(1) = 0
            stalat(1,count_level) = lat
            stalon(1,count_level) = lon
            write(c5_staid(1)(2:5),44) ob_counter
 44         format(i4.4)
            c5_staid(1)(1:1) = 's'
c            a9time_ob (1,count_level) = c_time_record
            a9time_ob (1,count_level) = m_time_record
            c8_obstype(1) = 'goes    '
            write (c8_obstype(1)(5:6),45)id_sat
 45         format(i2.2)
            height_m (1,count_level) = rbuf(lz)
            pressure_mb (1,count_level) = rbuf(lp)
            temp_c (1,count_level) = rbuf(lt)-273.15
            dewpoint_c (1,count_level) = rbuf(ltd)-273.15
            dir_deg (1,count_level) = mdf
            spd_mps (1,count_level) = mdf
c            staelev(1) = height_m(1,1)
            staelev(1)  = 0.0
          endif
 1060     format(1x,4(f8.2,4x))
  170   continue

c     finish off with stuff for writing 
        nlvl(1) = count_level ! number of levels

c     call write routine
        write(6,*) 'including satsnd in output file'
        call write_snd (lun_out,
     1       maxsnd,maxlvl,nsnd,
     1       iwmostanum,
     1       stalat,stalon,staelev,
     1       c5_staid,a9time_ob,c8_obstype,
     1       nlvl,
     1       height_m,
     1       pressure_mb,
     1       temp_c,
     1       dewpoint_c,
     1       dir_deg,
     1       spd_mps,
     1       istatus)
c
  200 continue
      go to 9001
 9000 continue
      write (6,*) 'routine process_snd failed on open, 9000'

c     call sdest('open error',0)
c
 9999 continue
      write (6,*) 'routine process_snd failed on read, 9999'
      istatus   = 0
      close(70)
      return

 9001 continue
      write (6,*) 'routine process_snd success on read'
      istatus = 1
      close (70)
      return
      end


      subroutine endian4(byte)
c
      integer*1 byte(4),tmp(4)
c
      do j = 1, 4
        tmp(j) = byte(j)
      enddo  !  j
c
      do j = 1, 4
        j1 = 5-j
        byte(j1) = tmp(j)
      enddo  !  j                      
c
      return
      end

