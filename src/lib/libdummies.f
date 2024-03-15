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
c==========================================
c
        function rshow_timer()

        common /timer/ i4time_start,sec_start,icount_start

        real sec_elapsed
        character*24 atime

        integer init
        data init/0/
        save init

        call system_clock(icount,icount_rate,icount_max)
        t1 = float(icount) / float(icount_rate)

        if(init .eq. 0)then
            sec_start = t1
            icount_start = icount
            write(6,*)'initializing rshow_timer, count = ',icount             
            init = 1
            rshow_timer = 0.
        else
            i4time_now = i4time_now_gg()
            call cv_i4tim_asc_lp(i4time_now,atime,istatus)

            sec_elapsed = float(icount - icount_start)/
     1                    float(icount_rate)                  
!           write(6,*)' count = ',icount

!           write(6,1)sec_elapsed,1./float(icount_rate)
!1          format(' elapsed seconds = ',f10.3,5x
!    1             ,'resolution = ',f10.6)       

            min_elapsed = int(sec_elapsed / 60.)
            sec_remainder = sec_elapsed - float(min_elapsed*60)
            write(6,2,err=99)min_elapsed,sec_remainder,atime(1:20)
     1                      ,1./float(icount_rate)
! hongli jiang: w>=d+3 from f8.6 to f9.6 11/27/2013
 2          format(1x,' elapsed time -',i6,':',f6.3,5x,a20
     1            ,'     resolution = ',f9.6)

 99         rshow_timer = sec_elapsed    
        endif

        return

        end
c
c==========================================
c
        function ishow_timer()

        integer sec_elapsed

        character*24 atime

        common /timer/ i4time_start

!       save i4time_start

        i4time_now = i4time_now_gg()

        i4time_elapsed = i4time_now - i4time_start

        min_elapsed = i4time_elapsed / 60

        sec_elapsed = i4time_elapsed - min_elapsed * 60

        call cv_i4tim_asc_lp(i4time_now,atime,istatus)

        if(sec_elapsed .ge. 10)then
            write(6,1,err=99)min_elapsed,sec_elapsed,atime(1:20)
 1          format(1x,' elapsed time -',i6,':',i2,5x,a20)
        elseif(sec_elapsed .ge. 0)then
            write(6,2,err=99)min_elapsed,sec_elapsed,atime(1:20)
 2          format(1x,' elapsed time -',i6,':0',i1,5x,a20)
        else
            write(6,*)' elapsed time (sec) = ',sec_elapsed
        endif

        ishow_timer = i4time_elapsed

        return

 99     write(6,*)' error in ishow_timer: i4time_now, i4time_start = '
     1                  ,i4time_now,i4time_start
        write(6,*)' this should not happen, stopping program'
        stop

        return
        end
c
c==========================================
c
        function init_timer()

        character*24 atime

        common /timer/ i4time_start
!       save i4time_start

        i4time_start = i4time_now_gg()

        call cv_i4tim_asc_lp(i4time_start,atime,istatus)

        write(6,*)'initializing elapsed timer at ',atime(1:20)  

        init_timer = 1

        return
        end
c
c==========================================
c
        function ltest_log_gg()
c
c       write(*,20)
20      format('called lib$init_timer')
c
        ltest_log_gg = 1

        return
        end
c
c==========================================
c
        subroutine notify_exec_gs(izero,             ! input
     1                  prod_array,          ! input
     1                  i4time_array,        ! input
     1                  istatflag,           ! output?
     1                  j_status)            ! input

        integer prod_array(10)
c
cd      write(*,60)
c60     format('called lib$find_file_end')
c
        return
        end
