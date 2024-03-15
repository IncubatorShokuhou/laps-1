cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
c
        subroutine make_fnam_lp (i4time, file_name, istatus)
c
cdoc    make_fnam_lp constructs the file name string 'yyjjjhhmm' for
cdoc    the time corresponding to i4time (seconds since 1-1-1960).
c
c================================================================
c
        integer i4time, istatus
        character*9 file_name
c
        integer nyear, nmonth, nday, nhour, nmin, nsec, njulian, i
c
        integer njul_days(12)
        data njul_days/0,31,59,90,120,151,181,212,243,273,304,334/
c
c================================================================
c
c convert i4 time to string yyjjjhhmm.
c
        call cv_i4tim_int_lp (i4time, nyear, nmonth, nday, nhour, nmin,
     1                      nsec, istatus)
        if (istatus .ne. 1) go to 100

        njulian = njul_days(nmonth) + nday

        if (nmonth .gt. 2 .and. mod (nyear,4) .eq. 0) njulian=njulian +
     11

        nyear = mod(nyear,100) ! steve albers 1997

        write(file_name,1001,err=90) nyear, njulian, nhour, nmin
1001    format (i2.2,i3.3,i2.2,i2.2)

        do i = 1, 9
           if (file_name(i:i) .eq. ' ') file_name(i:i) = '0'
        end do

        istatus = 1
        return
c
c error in encode.
c
90      istatus = 0
        write( 6,* ) 'error in make_fnam_lp: error in encode.'
        return
c
c error in subroutine...
c
100     istatus = 0
        write( 6,* ) 'error in make_fnam_lp: error in subroutine.'
        return
        end
c
        subroutine cv_i4tim_int_lp (i4time,nyear,nmonth,nday,nhour,
     1                     nmin,nsec,istatus)
c
cdoc    cv_i4tim_int_lp converts i4time (seconds since 1-1-1960) to six integers
cdoc    note that nyear is the number of years since 1900
c
c================================================================
c
c       integer i4time, nyear, nmonth, nday, nhour, nmin, nsec, istatus
c
        integer nsecmo(12), nsecyr, nsecda, nsechr, nsecmn
        integer ibase, i, nn, lftovr, ndays

        parameter (ibase=60)

        data nsecyr/31536000/nsecda/86400/nsechr/3600/nsecmn/60/
        data nsecmo/2678400,0,2678400,2592000,2678400,
     1  2592000,2678400,2678400,2592000,2678400,2592000,2678400/
c
c================================================================
c
        istatus = 1
c
c verify input.
c
        if (i4time .lt. 0) then
           istatus = 0
           write(6,*) 'error in cv_i4tim_int_lp: negative time ', i4time
           return
        end if

c
c subtract out number of years.
c
        lftovr = i4time
        nn = lftovr
        do i = 0, 70
           if (mod(i,4) .eq. 0) then
              nn = nn - (nsecyr + nsecda)
           else
              nn = nn - nsecyr
           end if
           if (nn .lt. 0) go to 8
           lftovr = nn
        end do
8       nyear = ibase + i

c
c subtract out number of months.
c
        nsecmo(2) = 2419200
        if (mod(nyear,4) .ne. 0) go to 10
           nsecmo(2) = 2505600
10      nn=lftovr
        do i=1,12
           nn = nn - nsecmo(i)
           if (nn .lt. 0) go to 30
           lftovr = nn
        end do
30      nmonth = i

c
c subtract out number of days.
c
        ndays = lftovr / nsecda
        lftovr = lftovr - (ndays * nsecda)
        nday = ndays + 1

c
c subtract out number of hours.
c
        nhour = lftovr / nsechr
        lftovr = lftovr - (nhour * nsechr)

c
c subtract out number of minutes.
c
        nmin = lftovr / nsecmn
        lftovr = lftovr - (nmin * nsecmn)

c
c what's left over is number of seconds.
c
        nsec = lftovr
        return
        end

      subroutine make_fnam13_lp(initial_i4time,forecast_time,filename,
     +     status)

cdoc  converts initial time and forecast time to a 13 character filename

      integer initial_i4time, forecast_time, status
      character*13 filename

      call make_fnam_lp(initial_i4time,filename,status)
      write(filename(10:13),'(i2.2,i2.2)') forecast_time/3600,
     +     mod(forecast_time,60)
      return
      end

      subroutine c_time2fname(utime,a9time)

cdoc  convert utime to a9time. jacket routine that calls 'make_fnam_lp'.

      integer utime, i4time, istatus
      character*(*) a9time

      i4time = utime + 315619200
      call make_fnam_lp (i4time,a9time,istatus)

      return
      end


      subroutine afwa_julhr_i4time(i_a1jul,i_a1min,i4time)

cdoc  i_a1jul is number of hours since dec 31, 1967 at 00z
cdoc  this is converted to i4time, number of sec since jan 1, 1960 at 00z
      i4time_hr  = i_a1jul * 3600 + (8*365 - 1 + 2) * 86400
      i4time_min = i_a1min*60
      i4time     = i4time_hr + i4time_min

      return
      end


        subroutine cv_asc_i4time(ascii_time,i4time)

cdoc    converts 9 character ascii time into i4time (seconds since 1-1-1960)

        character*9 ascii_time
        integer i4time ! seconds since 1-1-1960

        read(ascii_time(1:2),2,err=900)iyear
        read(ascii_time(3:5),3,err=900)idoy
        read(ascii_time(6:7),2,err=900)ihour
        read(ascii_time(8:9),2,err=900)imin

2       format(i2)
3       format(i3)

!       valid for years 1960-2060
        if(iyear .lt. 60)iyear = iyear + 100

        lp = (iyear + 3 - 60) / 4

        i4time =  (iyear-60) * 31536000
     1  + (idoy-1)   * 86400
     1  + ihour      * 3600
     1  + imin       * 60

        i4time = i4time + 86400 * lp

        return

!       error return
900     write(6,*)' error in cv_asc_i4time: ascii_time = ',ascii_time

        return

        end

c
        subroutine      cv_i4tim_asc_lp(i4time,atime,istatus)
c
cdoc  takes in an i4time and returns the time as an ascii string
cdoc  (e.g. 27-mar-1990 12:30:00.00 ).  the i4time is assumed to
cdoc  be a 1960-relative time, although the starting year is easily
cdoc  changed in the code.
c
c     imports - i4time ! seconds since 1-1-1960
c
c     exports - atime, istatus
c
c================================================================
c

        implicit        none

        integer       i4time,
     1          istatus,
     1          rmndr,
     1          nsec,
     1          monthsec(12),
     1          year,
     1          month,
     1          day,
     1          hour,
     1          min,
     1          sec

        character*24    atime
        character*4     ayear
        character*3     amonth(12)
        character*2     aday
        character*2     ahour
        character*2     amin
        character*2     asec

        data            monthsec/2678400,2419200,2678400,2592000,
     1                   2678400,2592000,2678400,2678400,
     1                   2592000,2678400,2592000,2678400/

        data            amonth/'jan','feb','mar','apr','may','jun',
     1                 'jul','aug','sep','oct','nov','dec'/

c
c================================================================
c

        if (i4time .lt. 0) then
           istatus=0
           write (6,*) 'error in input to cv_i4tim_asc_lp: negative time
     1'
           return
        endif

        rmndr=i4time
        do year=1960,2100
                if (mod(year,4) .eq. 0) then
                        nsec=31622400
                else
                        nsec=31536000
                endif
                if (rmndr .lt. nsec) goto 10
                rmndr=rmndr-nsec
        enddo

10      do month=1,12
                nsec=monthsec(month)
                if (mod(year,4) .eq. 0 .and. month .eq. 2) nsec=nsec+864
     100
                if (rmndr .lt. nsec) goto 20
                rmndr=rmndr-nsec
        enddo

20      do day=1,31
                if (rmndr .lt. 86400) goto 30
                rmndr=rmndr-86400
        enddo

30      do hour=0,23
                if (rmndr .lt. 3600) goto 40
                rmndr=rmndr-3600
        enddo

40      do min=0,59
                if (rmndr .lt. 60) goto 50
                rmndr=rmndr-60
        enddo

50      sec=rmndr

!       encode(4,900,ayear) year
!       encode(2,901,aday) day
!       encode(2,901,ahour) hour
!       encode(2,901,amin) min
!       encode(2,901,asec) sec

        write(ayear,900) year
        write(aday,901) day
        write(ahour,901) hour
        write(amin,901) min
        write(asec,901) sec
900     format(i4)
901     format(i2)

!       if (day .lt. 10) aday(1:1)='0'
        if (hour .lt. 10) ahour(1:1)='0'
        if (min .lt. 10) amin(1:1)='0'
        if (sec .lt. 10) asec(1:1)='0'

        atime=aday//'-'//amonth(month)//'-'//ayear//' '//ahour//':'//
     1      amin//':'//asec//'.00 '

        istatus=1
        return
        end

      subroutine jd_to_i4time(jd,i4time,istatus)

      double precision jd, jd_1960
      integer i4time ! seconds since 1-1-1960

cdoc  converts julian day (number of days since jan 1 4713bce to i4time)

!     author: steve albers 2005

      jd_1960 = 2436934.5

      days_since_1960 = jd - jd_1960

      i4time = int(days_since_1960 * 86400.)

      istatus = 1

      return
      end

      subroutine i4time_to_jd(i4time,jd,istatus)

      double precision jd, jd_1960, days_since_1960
      integer i4time,istatus ! seconds since 1-1-1960

cdoc  converts i4time to julian day (number of days since jan 1 4713bce)

!     author: steve albers 2013

      jd_1960 = 2436934.5

      days_since_1960 = dble(i4time) / 86400.d0

!     write(6,*)' days since 1960 is ',days_since_1960

      jd = jd_1960 + days_since_1960              

      istatus = 1

      return
      end
