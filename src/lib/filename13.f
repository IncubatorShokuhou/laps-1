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
        function filename13(i4time,ext)

!       steve albers            1990

!       this routine constructs a filename from the logical date_time and
!       a passed in extension

        character ext*3,filename13*13,cdum13*13

        character*9 asc9_time

        common /laps_diag/ no_laps_diag

        call make_fnam_lp(i4time,asc9_time,istatus)

        cdum13 = asc9_time//'.'//ext

        call downcase(cdum13,cdum13)

        if(no_laps_diag .eq. 0)then
            write(6,*)' filename13 = ',cdum13
        endif

        filename13 = cdum13

        return

        end

        function filename14(i4time,ext)

!       steve albers            1990

!       this routine constructs a filename from the logical date_time and
!       a passed in extension

        character ext*4,filename14*14,cdum14*14

        character*9 asc9_time

        common /laps_diag/ no_laps_diag

        call make_fnam_lp(i4time,asc9_time,istatus)

        cdum14 = asc9_time//'.'//ext

        call downcase(cdum14,cdum14)

        if(no_laps_diag .eq. 0)then
            write(6,*)' filename14 = ',cdum14
        endif

        filename14 = cdum14

        return

        end
