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

        subroutine wait_for_data(c_filespec,i4time_desired   !  inputs
     1               ,i4_check_interval,i4_total_wait        !  inputs
     1               ,i4_thresh_age       ! only loop through the waiting
                                          ! if most recent data is younger
                                          ! than this threshold age (input)
     1               ,istatus)            ! output

!       steve albers fsl 1995, 2001

!       note all times are in seconds

        character*255 c_filespec      ! wild card specifying the full path
!       and filenames. note: if there are very many files in the directory,
!       you will need to use a directory name only, the ls command can't
!       handle it otherwise. the wildcard is especially useful if there
!       is more than one type of file in the directory.

!       the c_filespec string should not have a particular file time within it.
!       valid examples of c_filespec are 'directory/*ext' and 'directory'.

        logical l_waited

        character*8 c8_project

        call get_c8_project(c8_project,istatus)
        if(istatus .eq. 0)return

        l_waited = .false.

        i4time_start_wait = i4time_now_gg()

        write(6,*)' wait for data:'
     1           ,' i4_check_interval,i4_total_wait,i4_thresh_age'
        write(6,*)  i4_check_interval,i4_total_wait,i4_thresh_age

 10     call get_latest_file_time(c_filespec,i4time_nearest)

        i4_data_age = i4time_desired-i4time_nearest

        write(6,*)' wait for data:'
     1           ,' i4time_desired,i4time_nearest,i4_data_age'
        write(6,*)  i4time_desired,i4time_nearest,i4_data_age

        if(i4time_nearest .eq. 0)then
            istatus = 0
            return
        endif

        if(i4_data_age .gt. 0 .and. i4_data_age .le. i4_thresh_age)then
            i4time_current = i4time_now_gg()
            i4_wait_sofar = i4time_current - i4time_start_wait

            if(i4_wait_sofar .lt. i4_total_wait)then
                write(6,*)' searching for more data, sleep '
     1                   ,i4_check_interval
                r4_check_interval = float(i4_check_interval)
                call snooze_gg(r4_check_interval,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' error in wait_for_data from snooze'
                    return
                endif
   
                l_waited = .true.

                goto 10

            endif

        endif

!       assess the latest data age
        if(i4_data_age .eq. 0)then
            write(6,*)' wait for data: found the data'
            istatus = 1
        elseif(i4_data_age .gt. i4_thresh_age)then
            write(6,*)' latest data too old, not waiting for new data'
            istatus = 0
        elseif(i4_data_age .lt. 0)then
            write(6,*)' wait for data: '
     1               ,'found the data (later than requested time)'
            istatus = 1
        else
            write(6,*)' wait for data: never did find the data'
            istatus = 0
        endif


        if(istatus .eq. 1 .and. c8_project(1:3) .eq. 'wfo')then       
            if(l_waited)then
                r_additional_wait = 60.
            else
                r_additional_wait = 50.
            endif

            write(6,*)' doing wait to allow full file update: '
     1               ,r_additional_wait
            call snooze_gg(r_additional_wait,istatus)
            if(istatus .ne. 1)then
                write(6,*)' error in wait_for_data returned from snooze'
                return
            endif
        endif

        return
        end
