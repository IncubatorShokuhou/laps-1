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
        subroutine i4time_fname_lp (fname_in, i4time, istatus)

        include 'lapsparms.for'
c
cdoc    this routine converts several different file name types (i.e. yydddhhmm)
cdoc    into the corresponding i4 time.
c
c       on input
c          file_name - the file name including the directory path if wanted.
c
c       on output
c          i4time - the corresponding i4 time of the file name.
c          istatus -  the return status.
c
c================================================================
c
        character*9 wfo_fname13_to_fname9, rsa13_to_a9, a8_to_a9
        character*9 a7_to_a9_time,yr_a10_to_a9
        character*9 file_name
        character*(*) fname_in
        character*256 fname_buf
        character*20 c20_type
        integer i4time, istatus
c
        integer int_file(9), i, nyear, jday, nhour, min, month, nday
        integer i4time_int_lp
c
c================================================================
c
c       check the length of file_name to be sure that it is valid.
c 
        fname_buf = fname_in
        call get_filetime_type(fname_in,c20_type,leni,lent)

        if(c20_type .eq. 'yyjjjhhmm')then                   ! nimbus/laps type
           file_name = fname_in(leni+1:leni+lent)

        elseif(c20_type .eq. 'yyyymmdd_hhmm')then           ! wfo type
           file_name = wfo_fname13_to_fname9(fname_in(leni+1:leni+lent))

        elseif(c20_type .eq. 'yyyyjjjhhmmss')then           ! rsa type
           file_name = rsa13_to_a9(fname_in(leni+1:leni+lent))

        elseif(c20_type .eq. 'yymmddhh')then                ! afwa raob type
           file_name = a8_to_a9(fname_in(leni+1:leni+lent))

        elseif(c20_type .eq. 'ymmddhh')then               ! taiwan fa model
           file_name = a7_to_a9_time(fname_in(leni+1:leni+lent))

        elseif(c20_type .eq. 'yyyymmddhh')then
           file_name = yr_a10_to_a9(fname_in(leni+1:leni+lent))

        else                                                ! unrecognized type
c          write(6,*)'i4time_fname_lp: unable to convert to i4time',
c    1'    type = ',c20_type
           istatus = 0
           fname_in = fname_buf
           return
        endif
c
c       split the file name into individual integers, checking for invalid
c       characters in the file name.
c
        do i = 1,9
           int_file(i) = ichar(file_name(i:i))-48
           if (int_file(i) .lt. 0 .or. int_file(i) .gt.9) then
              go to 1500
           end if
        end do
c
c       convert these numbers in year, day of year (julian day), hours, 
c       and minutes.
        nyear = 10*int_file(1) + int_file(2)
        jday = 100*int_file(3) + 10*int_file(4) + int_file(5)
        nhour = 10*int_file(6) + int_file(7)
        min = 10*int_file(8) + int_file(9)

c       get the actual year, here is where we assume what century it is...

        iyear_earliest_century = (iyear_earliest/100) * 100
        iyy_cutoff = iyear_earliest - iyear_earliest_century

        if(nyear .lt. iyy_cutoff)then
            nyear = nyear + iyear_earliest_century + 100
        else
            nyear = nyear + iyear_earliest_century
        endif
c
c       convert the day of year (julian day) into month and day.
c
        call cv_jul_mmdd_lp (jday, nyear, month, nday, istatus)
        if (0 .eq. istatus)
     1    go to 2000
c
c       convert these integers into i4 time.
c
        i4time = i4time_int_lp (nyear, month, nday, nhour, min, 0, istat
     1us)
        if (0 .eq. istatus)
     1    go to 2000
c
c       successful return
c
        istatus = 1
        return
c
c       error detected in this routine.
c
 1500   continue
        istatus = 0
        write( 6,* ) 'error in i4time_fname_lp: bad digit in file name.'
        write( 6,* ) file_name
        return
c
c       error detected in called routine
c
 2000   continue
        istatus = 0
        write( 6,* ) 'error in i4time_fname_lp: error in subroutine.'
        write( 6,* ) file_name
        return

 1000   continue
        istatus = 0
        write( 6,* ) 'error in i4time_fname_lp: wrong length for file na
     1me.'
        write( 6,* ) file_name

        return
c
        end
c
        subroutine cv_jul_mmdd_lp (julian_day, year, month, day, istatus
     1)
c
cdoc    this routine converts from day of year (julian days) to month and 
cdoc    day in integer format.
c
c       on input
c          julian_day - the day of year (julian date) to be converted.
c          year - the year of the day of year (julian date), if .lt. 100, 
c                 it is assumed to be the last two digits of 19xx.
c
c       on output
c          month - the integer representation for the month (1-12)
c          day - the integer value of the day of the month.
c          istatus - the return status.
c
c================================================================
c
        integer julian_day, year, day, istatus
c
        integer mnth(12), temp_year, month, max_day
        logical*1 leap
        data mnth/31,0,31,30,31,30,31,31,30,31,30,31/
c
c================================================================
c
        if (year.lt.100) then
           temp_year = year + 1900
        else
           temp_year = year
        endif
c
c       check to see if the year is a leap year
c
        leap = ( mod(temp_year,4)  .eq. 0 )
c
        if ( leap ) then
           max_day = 366
           mnth(2) = 29
        else
           max_day = 365
           mnth(2) = 28
        end if
c
c       check day of year (julian_day) for being too large or too small
c
        if (julian_day .le. 0 .or. julian_day .gt. max_day) then
           istatus = 0
           write( 6,* ) 'error in cv_jul_mmdd_lp.f: bad julian day', jul
     1ian_day
           return
        end if
c
c       convert to month and day
c
        month = 1
        day = julian_day
        do while (day .gt. mnth(month))
           day = day - mnth(month)
           month = month + 1
        end do
        istatus = 1
c
        return
        end

c
        function i4time_int_lp (nyear,nmonth,nday,nhour,nmin,nsec
     1                         ,istatus)
c
cdoc    i4time_int_lp returns i4 time (# of seconds since 00:00 01-jan-60)
cdoc    given a 6 integers containing year, month, day, hour,
cdoc    minute, second
c
c================================================================
c
        integer i4time_int_lp
        integer nyear, nmonth, nday, nhour, nmin, nsec, istatus
c
        integer nsecyr, nsecda, nsechr, nsecmn
        integer nyr, nyrs, nleap
        integer ndays(12), nsecmo(12)
        integer ibase, isum, i
        parameter (ibase=60)

        data nsecyr/31536000/nsecda/86400/nsechr/3600/nsecmn/60/
        data nsecmo/2678400,2419200,2678400,2592000,2678400,
     1  2592000,2678400,2678400,2592000,2678400,2592000,2678400/

        data ndays/31,29,31,30,31,30,31,31,30,31,30,31/
c
c================================================================
c
        istatus = 1
c
c sum the number of years.
c
        nyr=nyear
        if (nyr.gt.1900) nyr=nyr-1900
        nyrs=nyr-ibase
        if (nyrs.lt.0.or.nyrs.gt.67) go to 1000
        isum=nyrs*nsecyr
        nleap=(nyrs/4)+1                        ! account for leap years.
        isum=isum+(nleap*nsecda)
c
c sum in the number of months.
c
        if (nmonth.lt.1.or.nmonth.gt.12) go to 1000
c
        if (nmonth.ne.1) then
                do i=1,nmonth-1
                        isum=isum + nsecmo(i)
                end do
        end if
c
c correct for jan or feb of a leap yer.
c
!       if (mod(nyrs,4).eq.0) then
        if (mod(nyr,4).eq.0) then       ! steve albers, linda wharton 1993
                if (nmonth.le.2) isum=isum-nsecda
        end if
c
c sum in the number of days.
c
        if (nday.lt.1.or.nday.gt.ndays(nmonth))
     1  go to 1000
        isum=isum+((nday-1)*nsecda)
c
c sum in the number of hours.
c
        if (nhour.lt.0.or.nhour.gt.23) go to 1000
        isum=isum+(nhour*nsechr)
c
c sum in the number of minutes.
c
        if (nmin.lt.0.or.nmin.gt.59) go to 1000
        isum=isum+(nmin*nsecmn)
c
c sum in the number of seconds.
c
        if (nsec.lt.0.or.nsec.gt.59) go to 1000
        isum=isum+nsec
c
        i4time_int_lp=isum
        return
c
c return error code.
c
1000    istatus = 0
        i4time_int_lp=0
        write(6,*) 'error in i4time_int_lp'
        return
c
        end

