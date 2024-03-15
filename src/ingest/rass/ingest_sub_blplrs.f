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



        subroutine ingest_blplrs(i4time_sys,nx_l,ny_l,istatus)

c       steve albers               apr-1996    blp laps ingest
!       ken dritz                1-jul-1997  added nx_l, ny_l as dummy
!                                            arguments.
!       ken dritz                1-jul-1997  changed include of lapsparms.for
!                                            to laps_static_parameters.inc.
!       ken dritz                1-jul-1997  added call to get_r_missing_data.
!       ken dritz                8-jul-1997  replaced laps_domain_file by
!                                            'nest7grid' and removed include
!                                            of laps_static_parameters.inc.
c
c       this file shows examples of how the use prof_cdf subroutines to read
c       wpdn 60-minute rass data in netcdf files.
c
        integer cdfid,status,i,j,max_levels,max_stations
        parameter (max_levels = 100)
        parameter (max_stations = 1000)

        real temp(max_levels),prs
c       character*1 qc_flag(max_levels)
        integer i_qc_flag(max_levels)
        real rlevels(max_levels)
        integer good,bad,missing, start(2), count(2), stanamlen
        integer start_time(1), count_time(1)
        parameter (good = 0)
        parameter (bad = 8)
        parameter (missing = -1)
        character*1 qc_char(3)
        data qc_char/'g','b','m'/
        character*6 staname
c       character*1 submode
        character*6 pltc2_name
        data pltc2_name/'pltc2 '/
        character*200 fnam_in 
        character*180 dir_in
        character*255 c_filespec
        character*5 c5_data_interval

        integer wsmr_wmo_id
        integer error_code
        data wsmr_wmo_id/74533/
        data error_code/1/
        integer byte_to_i4
c
c       set error handling mode.  note that you don't have to do this, if this
c       call isn't made, default error processing will occur:
c
c       error_code              meaning
c
c       0                       return status codes -- this is the default.
c       1                       return status codes and write an error message
c                               to standard output (sys$output on vms).
c       2                       write an error message to standard output and
c                               exit the program.
c

        integer varid
        include 'netcdf.inc'
        character*(maxncnam) dimname

        character*13 filename13,c13_dum
        character*9 asc9_tim,a9time_ob

        character*40 c_vars_req
        character*180 c_values_req

        character*31    ext
        character*6 prof_name(max_stations)

        character*9 a9_timeobs
        integer timeobs

        real lat(nx_l,ny_l),lon(nx_l,ny_l)
        real topo(nx_l,ny_l)

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting r_missing_data'
           return
        endif

        call prof_cdf_set_error(error_code,status)
        if(status.ne.0)then
                write(*,*)'bad set_error ',status
                return
        endif

        c13_dum = filename13(i4time_sys,'lrs')
        asc9_tim = c13_dum(1:9)
c
c       open a 60-minute rass netcdf file for 20:00:00.00 on julian date 217,
c       1993.  both 6-minute and 60-minute files have the same filename
c       convention (you're supposed to use different directories to hold the
c       different resolutions).
c
c       yyjjjhhmmhhmmo
c
c       yy   = last 2 digits of year
c       jjj  = julian date
c       hhmm = hour, minute (utc) of the data
c       hhmm = hour, minute (utc) of the data
c              (one of these is supposed to be observation time, one receipt
c               time, but both are identical for rass files created from
c               demonstration division tapes.)
c       o    = ascii "oh" character, not "zero" -- stands for observation.
c
c       for vms systems, you must explicitly put in a period after the oh if
c       you don't want a file extension, thus the string used in the open call.
c

        c_vars_req = 'path_to_raw_blprass'
        call get_static_info(c_vars_req,c_values_req,1,istatus)
        if(istatus .eq. 1)then
            write(6,*)c_vars_req(1:30),' = ',c_values_req
            dir_in = c_values_req
        else
            write(6,*)' error getting ',c_vars_req
            return
        endif

        call s_len(dir_in,len_dir_in)

        call get_domain_perimeter(nx_l,ny_l,'nest7grid',lat,lon,
     1            topo,1.0,rnorth,south,east,west,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading domain perimeter'
            return
        endif

        call get_laps_cycle_time(laps_cycle_time,istatus)
        if (istatus .ne. 1) then
           write(6,*)'error getting laps_cycle_time'
           return
        else
           write(6,*)'laps_cycle_time = ',laps_cycle_time
        endif

!       do we want hourly or 6 minute profiler data?
        if(laps_cycle_time .le. 900)then
            c5_data_interval = '0015o'
            write(6,*)' using 15 minute data'
            i4time_desired = (i4time_sys / 900) * 900
            i4_avg_wdw_sec = 900 ! default value
        elseif(laps_cycle_time .lt. 3600)then
            c5_data_interval = '0030o'
            write(6,*)' using 30 minute data'
            i4time_desired = (i4time_sys / 1800) * 1800
            i4_avg_wdw_sec = 1800 ! default value
        else
            c5_data_interval = '0100o'
            write(6,*)' using hourly data'
            i4time_desired = i4time_sys
            i4_avg_wdw_sec = 3600 ! default value
        endif

c       wait for the data
!       c_filespec = dir_in(1:len_dir_in)//'*0100o'
        if(len_dir_in .gt. 0)then
            c_filespec = dir_in(1:len_dir_in)//'*'//c5_data_interval
        else
            write(6,*)' path_to_blp_rass has zero length'
            istatus = 0
            return
        endif

        write(6,*)c_filespec(1:80)

        i4time_desired = i4time_sys

        i4time_stop_waiting = i4time_sys + 25 * 60
        i4time_now = i4time_now_gg()
        i4_wait_period = i4time_stop_waiting - i4time_now

        i4_check_interval = 10
        i4_total_wait = min(0,i4_wait_period) ! turn off wait
        i4_thresh_age = 3600

        open(31,file='zzzz', status = 'old', err=10)
        read(31,*,err=10)i4_check_interval
        read(31,*,err=10)i4_total_wait
        read(31,*,err=10)i4_thresh_age
 10     continue
        close(31)

        call wait_for_data(c_filespec,i4time_desired
     1               ,i4_check_interval,i4_total_wait
     1               ,i4_thresh_age       ! only loop through the waiting
                                          ! if data is younger than this
                                          ! threshold
     1               ,istatus)

        if(istatus .ne. 1)then
            write(6,*)' no recent data'
            return
        endif


c       read in the raw rass data

        fnam_in = dir_in(1:len_dir_in)//asc9_tim//c5_data_interval
        call s_len(fnam_in,len_fnam_in)
        write(6,*)fnam_in(1:len_fnam_in)
        call prof_cdf_open(fnam_in(1:len_fnam_in),cdfid,status)

c
c       status of 0 means success.  positive status codes and -1 are returned
c       from netcdf routines.  any netcdf errors encountered will cause an
c       error message to be written to standard output (sys$output in vms).
c       errors -2 through -6 are from prof_cdf routines and are explained in
c       the documentation for each prof_cdf routine.
c
        if(status.ne.0)then
            write(*,*)'bad open ',status
            return
        endif

! added to read from file and set lag_time 
c 	read global attribute avgtimeperiod from input file and set lag_time
        call prof_i4_avg_wdw(i4_avg_wdw_sec,cdfid,istatus)
        if(istatus .eq. 1)then
            write(6,*)' i4_avg_wdw_sec from file = ',i4_avg_wdw_sec
        else
            i4_avg_wdw_sec = 3600
            write(6,*)' ingest_sub_blplrs: '
     1               ,'warning: could not obtain i4_avg_wdw_sec'
            write(6,*)' assuming i4_avg_wdw_sec = ',i4_avg_wdw_sec
        endif

        lag_time = i4_avg_wdw_sec/2
c
c       open an output file.
c
        ext = 'lrs'
        call open_lapsprd_file_append(1,i4time_sys,ext(1:3),istatus)
        if(istatus .ne. 1)then
            write(6,*)' error opening output file'
            istatus = 0
            return
        endif

!       get the number of levels from the netcdf file
        status = nf_inq_dimid(cdfid,'level',varid)
        status = nf_inq_dim(cdfid,varid,dimname,n_levels)

        write(6,*)' # of levels = ',n_levels

        if(n_levels .gt. max_levels)then
            write(6,*)' too many levels in the data'
            istatus = 0
            return
        endif

!       get the number of profilers from the netcdf file
        status = nf_inq_dimid(cdfid,'recnum',varid)
        status = nf_inq_dim(cdfid,varid,dimname,n_profilers)
        write(6,*)' # of blp rasss = ',n_profilers
        if (n_profilers .gt. max_stations) then
          write(6,*)' too many profilers to process'
          istatus = 0
          return
        endif

!	get the number of characters in the station name
        status = nf_inq_dimid(cdfid,'stanamlen',varid)
        status = nf_inq_dim(cdfid,varid,dimname,stanamlen)

        start(1) = 1
        count(2) = 1
        count(1) = stanamlen

!       do ista = 1,max_stations

        do ista = 1,n_profilers

!         fill profiler name into prof_name(ista)
          status = nf_inq_varid(cdfid,'staname',varid)
          start(2) = ista
          status = nf_get_vara_text(cdfid, varid, start, count,
     1                prof_name(ista))
          prof_name(ista)(6:6) = ' '

          write(6,*)
c
c         get the surface pressure.  this time we'll use the 5-character 
c         site name (plus a terminating blank) to select the station.
c

          call prof_cdf_read(cdfid,prof_name(ista),0,'stalat',2,rlat
     1                      ,status)
          if(status.ne.0)then
                write(6,*)'bad lat read ',prof_name(ista),status
                go to 900
          endif
          write(6,*)
          write(6,*)prof_name(ista),' lat ',rlat

          call prof_cdf_read(cdfid,prof_name(ista),0,'stalon',2,rlon
     1                      ,status)
          if(status.ne.0)then
                write(6,*)'bad lon read ',status
                go to 900
          endif
          write(6,*)
          write(6,*)prof_name(ista),' lon ',rlon

!         get the observation time
          status = nf_inq_varid(cdfid,'timeobs',varid)
          start_time(1) = ista
          count_time(1) = 1
          status = nf_get_vara_int(cdfid, varid, start_time, 
     1                             count_time,timeobs)
          if(status.ne.nf_noerr)then
              write(6,*)' warning: bad timeobs read ',status,timeobs

          elseif(abs(timeobs) .gt. 3d9)then
              write(6,*)' warning: bad observation time',timeobs

          else
              call c_time2fname(timeobs,a9_timeobs)

              write(6,*)
              write(6,*)' timeobs ',a9_timeobs

              call cv_asc_i4time(a9_timeobs,i4_timeobs)
              i4_resid = abs(i4_timeobs - i4time_sys)
!             if(i4_resid .gt. (ilaps_cycle_time / 2) )then ! outside time window
              if(i4_resid .gt. 0)then
                  write(6,*)' warning, time is suspect '
     1                     ,a9_timeobs,i4_resid       
              endif

          endif

          if(rlat .le. rnorth .and. rlat .ge. south .and.
     1       rlon .ge. west   .and. rlon .le. east            )then

            write(6,*)prof_name(ista),' is in box'


            call prof_cdf_read(cdfid,prof_name(ista),0,'pressure',2,prs       
     1                        ,status)
            if(status.ne.0)then
                write(*,*)'bad pressure read ',status
                prs = r_missing_data
            endif
            write(6,*)
            write(6,*)prof_name(ista),' pressure ',prs

            call prof_cdf_read(cdfid,prof_name(ista),0,'temperature'
     1                        ,2,t_sfc,status)
            i_qc_sfc = 0
            if(status.ne.0)then
                write(*,*)'bad t_sfc read ',status
                i_qc_sfc = 0
            endif

            write(6,*)
            write(6,*)prof_name(ista),' t_sfc ',t_sfc,i_qc_sfc

            call prof_cdf_read(cdfid,prof_name(ista),0,'relhumidity'
     1                        ,2,rh_sfc,status)
            call prof_cdf_read(cdfid,prof_name(ista),0,'windspeedsfc'
     1                        ,2,sp_sfc,status)
            call prof_cdf_read(cdfid,prof_name(ista),0,'winddirsfc'
     1                        ,2,di_sfc,status)

            call prof_cdf_read(cdfid,prof_name(ista),0,'staelev'
     1                        ,2,elev,status)
            if(status.ne.0)then
                write(*,*)'bad elev read ',status
                go to 900
            endif
            write(6,*)
            write(6,*)prof_name(ista),' elev  ',elev

            n_levels_tot = n_levels

!           i4time_ob = i4time_sys - lag_time
            i4time_ob = i4_timeobs - lag_time 

            call make_fnam_lp(i4time_ob,a9time_ob,istatus)

            write(6,401)ista,n_levels_tot,rlat,rlon,elev
     1                 ,prof_name(ista)(1:5),a9time_ob,'rass    '
            write(1,401)ista,n_levels_tot,rlat,rlon,elev
     1                 ,prof_name(ista)(1:5),a9time_ob,'rass    '
401         format(i12,i12,f11.3,f15.3,f15.0,5x,a5,3x,a9,1x,a8)
c
c           get the array of rass virtual temperatures for the profiler station.
c
            call prof_cdf_read(cdfid,prof_name(ista),0,'virtualtemp',2,
     $                     temp,status)
            if(status.ne.0)then
                write(*,*)'bad virtualtemp read ',status
                go to 900
            endif
c
c           get the associated quality control flags.
c
            call prof_cdf_read(cdfid,prof_name(ista),0,'qualitycode'
     $                         ,0,i_qc_flag,status)
            if(status.ne.0)then
                write(*,*)'bad qualitycode read ',status
                go to 900
            endif
c
c           get the associated levels
c
            call prof_cdf_read(cdfid,prof_name(ista),0,'levels',2,
     $                     rlevels,status)
            if(status.ne.0)then
                write(*,*)'bad levels read ',status
                go to 900
            endif
c
            write(6,*)prof_name(ista),' virtualtemp'

!           write surface temperature (and other data) as first level
!           write(1,*)elev,t_sfc,' ',i_qc_sfc,rh_sfc,di_sfc,sp_sfc,prs

            rms = 1.0

            do i = 1, n_levels
                iqc_flag  = i_qc_flag(i)
                iqc_flag2 = i_qc_flag(i)

                if(iqc_flag.eq.good)then
                        j = 1
                        iqc = 1
                else if(iqc_flag.eq.bad)then
                        j = 2
                        iqc = 0
                else if(iqc_flag.eq.missing)then
                        j = 3
                        iqc = 0
                else ! cover other cases too
                        iqc = 0
                endif

                height_msl = rlevels(i) + elev

                if(temp(i) .gt. r_missing_data)then
                    temp_out = r_missing_data
                else
                    temp_out = temp(i)
                endif

                if(temp_out .lt. 100. .or. temp_out .gt. 400.)then
!                   write(6,*)' setting iqc to 0'
                    iqc = 0
                endif

                write(1,*)height_msl,temp_out,' ',iqc,rms
                write(6,*)height_msl,temp_out,' ',iqc,rms
     1                                           ,iqc_flag,iqc_flag2     
            enddo

c
c           get the station name for white sands.  this shows how to get a
c           character string:  use the number of characters as the number of
c           array elements.
c
c           call prof_cdf_read(cdfid,prof_name(ista),0,'staname',1,
c    $                     staname,status)
c           if(status.ne.0)then
c               write(*,*)'bad staname read ',status
c               go to 900
c           endif
c           write(6,*)
c           write(6,*)prof_name(ista),' staname ',staname
c
c           get submode (should be 1 character, a.)
c
c           call prof_cdf_read(cdfid,prof_name(ista),0,'submode',1,
c    $                     submode,status)
c           if(status.ne.0)then
c               write(*,*)'bad submode read ',status
c               go to 900
c           endif
c           write(6,*)
c           write(6,*)prof_name(ista),' submode ',submode

          else !
            write(6,*)prof_name(ista),' is outside of domain'

          endif ! l_in_box


 900    enddo                    ! ista


        close(1)

c
c       close the netcdf file.  this isn't necessary, but here's a sample of
c       the call.  a program could have up to 4 profiler netcdf files open at
c       a time.  the cdfid's are what indicate which file is which.
c
        call prof_cdf_close(cdfid,status)
        if(status.ne.0)then
            write(*,*)'bad close ',status
        endif
c
        return
        end




