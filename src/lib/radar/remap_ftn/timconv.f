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


        function int_to_i4time(iyr_in,imon,idy,ih,im,isec)

!       first generate a9time
        character*9 a9

        integer imon_a(12)

        data imon_a/0,31,59,90,120,151,181,212,243,273,304,334/

        id = imon_a(imon)

        idays = id + idy

!       determine tens and units digits for the year
        iyr = iyr_in - (iyr_in / 100) * 100

!       decide whether to add a day for leap year.
        if(iyr .eq. (iyr / 4) * 4 )then
            if(imon .ge. 3)then
                idays = idays + 1
            endif
        endif


        write(a9,2)iyr,idays,ih,im
2       format(i2,i3,i2,i2)

        if(a9(3:3) .eq. ' ')a9(3:3) = '0'
        if(a9(4:4) .eq. ' ')a9(4:4) = '0'
        if(a9(6:6) .eq. ' ')a9(6:6) = '0'
        if(a9(8:8) .eq. ' ')a9(8:8) = '0'

        call i4time_fname_lp(a9,i4time,istatus)

        int_to_i4time = i4time


        return
        end



        function a11_to_a9_dummy(a11_time)

        character*9 a11_to_a9_dummy
        character*9 a9
        character*11 a11_time

        if(a11_time(4:6) .eq. 'jan')id = 0
        if(a11_time(4:6) .eq. 'feb')id = 31
        if(a11_time(4:6) .eq. 'mar')id = 59
        if(a11_time(4:6) .eq. 'apr')id = 90
        if(a11_time(4:6) .eq. 'may')id = 120
        if(a11_time(4:6) .eq. 'jun')id = 151
        if(a11_time(4:6) .eq. 'jul')id = 181
        if(a11_time(4:6) .eq. 'aug')id = 212
        if(a11_time(4:6) .eq. 'sep')id = 243
        if(a11_time(4:6) .eq. 'oct')id = 273
        if(a11_time(4:6) .eq. 'nov')id = 304
        if(a11_time(4:6) .eq. 'dec')id = 334

        iyr = 94

        read(a11_time,1)idy,ih,im
1       format(i2,5x,i2,i2)

        idays = id + idy

        write(a9,2)iyr,idays,ih,im
2       format(i2,i3,i2,i2)

        if(a9(3:3) .eq. ' ')a9(3:3) = '0'
        if(a9(4:4) .eq. ' ')a9(4:4) = '0'
        if(a9(6:6) .eq. ' ')a9(6:6) = '0'
        if(a9(8:8) .eq. ' ')a9(8:8) = '0'

        a11_to_a9_dummy = a9

        return
        end



        function i4time_to_a11_dummy(i4time)

!       dd_mon_hhmm

        character*11 i4time_to_a11_dummy
        character*11 a11
!       character*9 a9_time
        character*24 asc_tim_vol
        character*3 c3_mon

!       call i4time_fname_lp(a9_time,i4time,istatus)

!       i4time to ascii 24 time
        call cv_i4tim_asc_lp(i4time,asc_tim_vol,istatus)

!       ascii 24 time to a11 time
        a11(1:2) = asc_tim_vol(1:2)   ! dy
        a11(3:3) = '_'

        c3_mon = asc_tim_vol(4:6)
        call downcase(c3_mon,c3_mon)

        a11(4:6) = c3_mon

        a11(7:7) = '_'
        a11(8:9) = asc_tim_vol(13:14) ! hr
        a11(10:11) = asc_tim_vol(16:17) ! mn

        i4time_to_a11_dummy = a11

        return
        end
