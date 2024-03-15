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
c given a ascii time such as 09-mar-1992 12:00:00.00, return the corresponding
c i4time (seconds since 1-1-1960). the use of the seconds field is a bit of 
c a historical relic. this function is mostly used to parse the types of dates 
c a user might type in.
c
        function        i4time_asc_gg(atime,istatus)

        implicit        none

        integer       i4time_asc_gg,  ! this function
     1          i4time_int_lp,  ! function call
     1          istatus,        ! argument (exported)
     1          year,
     1          month,
     1          day,
     1          hour,
     1          min,
     1          sec,
     1          i

        character*24    atime           ! argument (imported)
        character*4     ayear
        character*3     amonth(12)
        character*2     aday
        character*2     ahour
        character*2     amin
        character*2     asec

        data            amonth/'jan','feb','mar','apr','may','jun',
     1                 'jul','aug','sep','oct','nov','dec'/

        ayear=atime(8:11)
        aday=atime(1:2)
        ahour=atime(13:14)
        amin=atime(16:17)
        asec=atime(19:20)

        read(ayear,900)  year
        read(aday,901)  day
        read(ahour,901)  hour
        read(amin,901)  min
        read(asec,901)  sec
900     format(i4)
901     format(i2)

        do i=1,12
                if (atime(4:6) .eq. amonth(i)) then
                        month=i
                        goto 10
                endif
        enddo

 10     continue

        i4time_asc_gg=i4time_int_lp(year,month,day,hour,min,sec,istatus)

        return

        end
