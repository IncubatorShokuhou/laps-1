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



        subroutine ingest_pro(i4time_sys,nx_l,ny_l,lun_out,istatus)

c       michael barth           12-aug-1993
c       steve albers               nov-1993         reworked for laps ingest
c                                  oct-1994         improve qc
c       steve albers               sep-1996         wfo database compatability
!       ken dritz                3-jul-1997  added nx_l, ny_l as dummy
!                                            arguments.
!       ken dritz                3-jul-1997  changed include of lapsparms.for
!                                            to laps_static_parameters.inc.
!       ken dritz                3-jul-1997  added call to get_r_missing_data.
c
c       this file shows examples of how the use prof_cdf subroutines to read
c       wpdn 60-minute rass data in netcdf files.
c
c       note: profiler winds are written out in knots, and are sorted by height

        integer cdfid,status,i,j,max_levels,max_levels_out,n_profilers
     1         ,file_n_prof
        parameter (max_levels = 72)
        parameter (max_levels_out = 72)
        parameter (n_profilers = 200) ! accomodates rsa

        real u(max_levels),v(max_levels),prs
        real ht_out(max_levels_out),di_out(max_levels_out)
     1                             ,sp_out(max_levels_out)

        character*1 c1_qc_flag(max_levels)        ! for /public
        integer i4_qc_flag(max_levels)          ! for wfo

        real level(max_levels)
        integer good,bad,missing, start(2), count(2), stanamlen
        integer start_time(1), count_time(1)
        parameter (good = 0)
        parameter (bad = 12)
        parameter (missing = -1)
        character*1 qc_char(3)
        data qc_char/'g','b','m'/
        integer byte_to_i4

        character*200 fnam_in
        character*180 dir_in
        character*255 c_filespec

        include 'lapsparms.for'

        integer max_files

        parameter(max_files = max_ingest_files)
        character*255 c_filenames(max_files)

        integer wsmr_wmo_id
        data wsmr_wmo_id/0/
        integer error_code
        data error_code/1/
        logical l_in_box
        data l_in_box/.true./
        character*8 c8_project
        character*5 c5_data_interval
        character*1 c1_char

        integer varid
        integer n_levels 
        include 'netcdf.inc'
        character*(maxncnam) dimname 
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

        character*13 filename13,outfile,asc13_tim,fname9_to_wfo_fname13       
        character*9 asc9_tim,a9time_ob,a9time_infile

        character*31    ext
        integer       len_dir_in

        character*40 c_vars_req
        character*180 c_values_req

        character*6 prof_name(n_profilers)
        character*9 a9_timeobs
        integer timeobs

        real lat(nx_l,ny_l),lon(nx_l,ny_l)
        real topo(nx_l,ny_l)

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting r_missing_data'
           return
        endif

        call get_laps_cycle_time(laps_cycle_time,istatus)
        if (istatus .ne. 1) then
           write(6,*)'error getting laps_cycle_time'
           return
        else
           write(6,*)'laps_cycle_time = ',laps_cycle_time
        endif

        r_mspkt = .518

        call get_latlon_perimeter(nx_l,ny_l,1.0
     1                           ,lat,lon,topo
     1                           ,rnorth,south,east,west,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading laps perimeter'
            return
        endif

        call prof_cdf_set_error(error_code,status)
        if(status.ne.0)then
                write(6,*)'bad set_error ',status
                return
        endif

        outfile = filename13(i4time_sys,'pro')
        asc9_tim = outfile(1:9)
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

        c_vars_req = 'path_to_raw_profiler'
        call get_static_info(c_vars_req,c_values_req,1,istatus)
        if(istatus .eq. 1)then
            write(6,*)c_vars_req(1:30),' = ',c_values_req
            dir_in = c_values_req
        else
            write(6,*)' error getting ',c_vars_req
            return
        endif

        call s_len(dir_in,len_dir_in) ! slash should be on the end of dir_in

!       determine whether we are using /public or wfo advanced filenames...
        if(.false.)then !! wni-bls ... changed to false to force file name determination
            call get_c8_project(c8_project,istatus)
            if (istatus .ne. 1) then
               write(6,*)'error getting c8_project'
               return
            else
               write(6,*)'c8_project = ',c8_project
            endif

        else
c           determine file format by looking at the file name convention
            call get_file_names(dir_in(1:len_dir_in),numoffiles_ret
     1                         ,c_filenames,max_files,istatus)

!           note that the gfn call may not work unless we also call
!           'filter_nonnumeric_fnames'.
            call filter_non_numeric_fnames(c_filenames,
     1                   numoffiles_ret,
     1                   numoffiles,
     1                   max_files,
     1                   istatus)

            if(istatus .ne. 1 .or. numoffiles .eq. 0)then
                write(6,*)' error calling get_file_names'
                istatus = 0
                return
            endif
            ipos = len_dir_in + 9
            c1_char = c_filenames(1)(ipos:ipos)
            if(c1_char .eq. '_')then
                c8_project = 'wfo' 
            else
                c8_project = 'nimbus'
            endif
            write(6,*)' 9th character of filename: ',c1_char
            write(6,*)' setting c8_project parameter: ',c8_project

        endif

        if(c8_project(1:6) .eq. 'nimbus')then
            write(6,*)' assumming /public filename format'
        else
            write(6,*)' assumming wfo filename format'
        endif

!       do we want hourly or 6 minute profiler data?
        if(c8_project(1:6) .eq. 'nimbus')then
            if(laps_cycle_time .le. 1800)then
!           if(.true.)then
                c5_data_interval = '0006o'
                write(6,*)' using 6 minute data'
                i4time_desired = (i4time_sys / 360) * 360
            else
                c5_data_interval = '0100o'
                write(6,*)' using hourly data'
                i4time_desired = i4time_sys
            endif
            c_filespec = dir_in(1:len_dir_in)//'*'//c5_data_interval

        else ! wfo
            c_filespec = dir_in(1:len_dir_in)//'*'
            i4time_desired = i4time_sys

        endif

        call make_fnam_lp(i4time_desired,a9time_infile,istatus)

        write(6,*)c_filespec(1:80)

c       wait for the data
        i4_check_interval = 10
        i4_thresh_age = 3600
        i4time_now = i4time_now_gg()
        i4_total_wait = min(300,i4time_desired+25*60 - i4time_now)

        open(31,file='zzzz', status = 'old', err=10)
        read(31,*,err=10)i4_check_interval
        read(31,*,err=10)i4_total_wait
        read(31,*,err=10)i4_thresh_age
 10     continue
        close(31)

        if(i4_total_wait .gt. 0)then
            call wait_for_data(c_filespec,i4time_desired
     1               ,i4_check_interval,i4_total_wait
     1               ,i4_thresh_age       ! only loop through the waiting
                                          ! if data is younger than this
                                          ! threshold
     1               ,istatus)

            if(istatus .ne. 1)then
                write(6,*)' no recent data'
                return        ! normal action
!               continue      ! do this for testing on the wfo
            endif
        endif

c       read in the raw profiler data
        if(c8_project(1:6) .eq. 'nimbus')then
            fnam_in = dir_in(1:len_dir_in)//a9time_infile
     1                                    //c5_data_interval
        else ! wfo
!           convert from asc9_tim to asc13_tim
!           asc13_tim = '19960903_2200'                      ! hardwired for testing.
            asc13_tim = fname9_to_wfo_fname13(a9time_infile) ! john smart's routine
            fnam_in = dir_in(1:len_dir_in)//asc13_tim
        endif

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
            write(6,*)' warning: bad open ',status
            return
        endif

! added to read from file and set lag_time lw 8-27-98
c 	read global attribute avgtimeperiod from input file and set lag_time
        call prof_i4_avg_wdw(i4_avg_wdw_sec,cdfid,istatus)
        if(istatus .eq. 1)then
            lag_time = i4_avg_wdw_sec/2
        else
            write(6,*)' ingest_sub_pro: '
     1               ,'error obtaining i4_avg_wdw_sec'
            return
        endif

!       get the number of levels from the netcdf file
        status = nf_inq_dimid(cdfid,'level',varid)
        status = nf_inq_dim(cdfid,varid,dimname,n_levels)

        write(6,*)' # of levels = ',n_levels

        if(n_levels .gt. max_levels)then
            write(6,*)' error: too many levels in the data'
            istatus = 0
            return
        endif

!       get the number of profilers from the netcdf file
        status = nf_inq_dimid(cdfid,'recnum',varid)
        status = nf_inq_dim(cdfid,varid,dimname,file_n_prof)
        write(6,*)' # of profilers = ',file_n_prof
        if (file_n_prof .gt. n_profilers) then
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

        do ista = 1,file_n_prof

!       fill profiler name into prof_name(ista)
          status = nf_inq_varid(cdfid,'staname',varid)
          start(2) = ista
         
          status = nf_get_vara_text(cdfid, varid, start, count,
     1                prof_name(ista))
          prof_name(ista)(6:6) = ' '
c
c         get the surface pressure for platteville.  this time we'll use the
c         5-character site name (plus a terminating blank) to select the station.
c
          call prof_cdf_read(cdfid,prof_name(ista),0,'pressure',2,prs
     1                    ,status)
          if(status.ne.0)then
            if(status .eq. -3)then
                write(6,*)prof_name(ista),' not found'
                goto 900
            else
                write(6,*)' warning: bad pressure read ',status
                return
            endif
          endif
          write(6,*)
          write(6,*)prof_name(ista),' pressure ',prs

!         get the latitude
          call prof_cdf_read(cdfid,prof_name(ista),0,'stalat',2,rlat
     1                      ,status)
          if(status.ne.0)then
              write(6,*)' warning: bad lat read ',status
              return
          endif
          write(6,*)
          write(6,*)prof_name(ista),' lat ',rlat

!         get the longitude
          call prof_cdf_read(cdfid,prof_name(ista),0,'stalon',2,rlon
     1                      ,status)
          if(status.ne.0)then
              write(6,*)' warning: bad lon read ',status
              return
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

            call prof_cdf_read
     1          (cdfid,prof_name(ista),0,'staelev',2,elev,status)
            if(status.ne.0)then
                write(6,*)' warning: bad elev read ',status
                return
            endif
            write(6,*)
            write(6,*)prof_name(ista),' elev ',elev

            call prof_sfcob_read(cdfid,prof_name(ista)          ! i
     1                          ,r_missing_data                 ! i
     1                          ,di_sfc                         ! o
     1                          ,sp_sfc                         ! o
     1                          ,p_sfc_hpa                      ! o
     1                          ,t_sfc_k                        ! o
     1                          ,rh_sfc_pct                     ! o
     1                          ,status)                        ! o

!           test whether di and sp lie within valid range
            if(abs(sp_sfc) .gt. 500.)status = 1
            if(abs(di_sfc) .gt. 500.)status = 1

            if(status .eq. 0)then
                n_good_sfc = 1
            else
                n_good_sfc = 0
            endif
c
c           get the array of profiler winds for the profiler station
c           at platteville.  for this call, we'll use the wmo
c           identifier to tell the subroutine what station we want.
c
            call prof_cdf_read(cdfid,prof_name(ista),0,'ucomponent',2,
     1                                             u,status)
            if(status.ne.0)then
                write(6,*)' warning: bad ucomponent read ',status
                return
            endif

            call prof_cdf_read(cdfid,prof_name(ista),0,'vcomponent',2,
     1                              v,status)
            if(status.ne.0)then
                write(6,*)' warning: bad vcomponent read ',status
                return
            endif
c
c           get the associated quality control flags.
c
            call prof_cdf_read(cdfid,prof_name(ista),0
     $                 ,'uvqualitycode',0,i4_qc_flag,status)
            if(status.ne.0)then
                write(6,*)' warning: bad qualitycode read ',status       
                return
            endif
c
c           get the associated levels for this profiler (also works with
c           global data statement on nimbus)
c
            call prof_cdf_read(cdfid,prof_name(ista),0,'levels',2,
     $                     level,status)
            if(status.ne.0)then
                write(6,*)' warning: bad level read ',status
                return
            endif
c
            write(6,*)prof_name(ista)

            n_good_levels = 0

            do i = n_levels, 1, -1

                iqc_flag = i4_qc_flag(i)

                if(iqc_flag.eq.good)then
                        j = 1
                        iqc = 1
                else if(iqc_flag.eq.bad)then
                        j = 2
                        iqc = 0
                        write(6,*)prof_name(ista),' ',i,iqc_flag,' bad'
                else if(iqc_flag.eq.missing)then
                        j = 3
                        iqc = 0
                        write(6,*)prof_name(ista),' ',i,iqc_flag,' msg'
                else
                        iqc = 0
                        write(6,*)prof_name(ista),' ',i,iqc_flag,' qcd'
                endif
                height_msl = level(i) + elev
                write(6,*)i,height_msl,u(i),v(i),' ',iqc,iqc_flag ! ,iqc_flag2

                mode_flag = 1

                if(n_good_levels .ge. 1)then
                    if(height_msl .ge. ht_out(n_good_levels))then
                        mode_flag = 0 ! this is a low mode level repeating the
                                      ! wind from the high mode
                    endif
                endif


                if(  (.not.
     1    (u(i) .gt. r_missing_data .or. v(i) .gt. r_missing_data)  )
     1                          .and.
     1                      iqc .eq. 1
     1                          .and.
     1                           mode_flag .eq. 1
     1                                                          )then

                    n_good_levels = n_good_levels + 1
                    ht_out(n_good_levels) = height_msl
                    call uv_to_disp(u(i),v(i)
     1         ,di_out(n_good_levels),sp_out(n_good_levels))
                endif
            enddo ! i

            i4time_ob = i4_timeobs - lag_time ! i4time_sys - lag_time 

            call make_fnam_lp(i4time_ob,a9time_ob,istatus)
c
c           open an output file if needed
            ext = 'pro'
            call open_ext(lun_out,i4time_sys,ext,istatus)
            if(istatus .ne. 1)then
                write(6,*)' error opening product file',ext
                return
            endif

            write(6,401)wsmr_wmo_id,n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(ista),a9time_ob,'profiler'       
            write(1,401)wsmr_wmo_id,n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(ista),a9time_ob,'profiler'       
401         format(i12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)

            rms = 1.0

            if(n_good_sfc .eq. 1)then
!               write surface data as first level
                write(1,301)elev,di_sfc,sp_sfc,rms ! /r_mspkt
     1                          ,p_sfc_hpa 
     1                          ,t_sfc_k   
     1                          ,rh_sfc_pct
                write(6,301)elev,di_sfc,sp_sfc,rms ! /r_mspkt
     1                          ,p_sfc_hpa 
     1                          ,t_sfc_k   
     1                          ,rh_sfc_pct
            endif

            do i = n_good_levels, 1, -1
                write(1,301,err=303)ht_out(i),di_out(i),sp_out(i),rms 
                write(6,301,err=303)ht_out(i),di_out(i),sp_out(i),rms 
301             format(1x,f6.0,f6.0,2f6.1,3f7.1)
303             continue
            enddo ! i

        else !
             write(6,*)prof_name(ista),' is outside of domain'

        endif ! l_in_box

        write(6,*)prof_name(ista)

900     enddo ! stations


c       close the netcdf file.  this isn't necessary, but here's a sample of
c       the call.  a program could have up to 4 profiler netcdf files open at
c       a time.  the cdfid's are what indicate which file is which.
c
        call prof_cdf_close(cdfid,status)
        if(status.ne.0)then
                write(6,*)' warning: bad close ',status
        endif
c
        return
        end

        subroutine prof_sfcob_read(cdfid,prof_name              ! i
     1                          ,r_missing_data                 ! i
     1                          ,di_sfc                         ! o
     1                          ,sp_sfc                         ! o
     1                          ,p_sfc_hpa                      ! o
     1                          ,t_sfc_k                        ! o
     1                          ,rh_sfc_pct                     ! o
     1                          ,istatus)                       ! o

        integer istatus   ! a value of zero indicates good dir & speed

        integer cdfid
        character*8 c8_project
        character*(*) prof_name

        di_sfc = r_missing_data
        sp_sfc = r_missing_data

        call get_sfc_badflag(sfc_badflag,istatus)
        if(istatus .ne. 1)return

        t_sfc_k = sfc_badflag
        rh_sfc_pct = sfc_badflag
        p_sfc_hpa = sfc_badflag

        call get_c8_project(c8_project,istatus)
        if(istatus .ne. 1)return

        call prof_cdf_read(cdfid,prof_name,0
     1                         ,'windspeedsfc',2,sp_sfc,istat_sp)

        if(c8_project(1:6) .eq. 'nimbus')then
            call prof_cdf_read(cdfid,prof_name,0
     1                         ,'winddirsfc',2,di_sfc,istat_di)
        else
            call prof_cdf_read(cdfid,prof_name,0
     1                         ,'winddirsfc',0,i4_di_sfc,istat_di)

            di_sfc = i4_di_sfc
        endif

        call prof_cdf_read(cdfid,prof_name,0
     1                         ,'pressure',2,p_sfc_hpa,istatus2)

        call prof_cdf_read(cdfid,prof_name,0
     1                         ,'temperature',2,t_sfc_k,istatus2)

        call prof_cdf_read(cdfid,prof_name,0
     1                         ,'relhumidity',2,rh_sfc_pct,istatus2)

        istatus = istat_di * istat_sp

        if(abs(p_sfc_hpa) .gt. 5000.)p_sfc_hpa  = sfc_badflag
        if(abs(t_sfc_k) .gt. 500.)   t_sfc_k    = sfc_badflag
        if(abs(rh_sfc_pct) .gt. 500.)rh_sfc_pct = sfc_badflag

        return
        end
