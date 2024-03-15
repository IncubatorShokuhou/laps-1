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



        subroutine ingest_blppro(i4time_sys,nx_l,ny_l,lun_out,istatus)       

c       michael barth           12-aug-1993
c       steve albers               nov-1993         reworked for laps ingest
c                                  oct-1994         improve qc
c       steve albers                   1996         bl profilers
!       ken dritz                3-jul-1997  added nx_l, ny_l as dummy
!                                            arguments.
!       ken dritz                3-jul-1997  changed include of lapsparms.for
!                                            to laps_static_parameters.inc.
!       ken dritz                3-jul-1997  added call to get_r_missing_data.
c
c       this file shows examples of how the use prof_cdf subroutines to read
c       wpdn 60-minute rass data in netcdf files.
c
c       note: profiler winds are written out in knots

        integer cdfid,status,i,j,max_levels_out,n_profilers,file_n_prof       
        parameter (n_profilers = 1000)

	parameter (max_modes = 3)
	parameter (max_levels = 50)
        parameter (max_levels_out = max_modes * max_levels)
        character*2 nmodes_short
	integer nmodes
!       equivalence(nmodes,nmodes_short)

        real u(max_modes,max_levels)
        real v(max_modes,max_levels), prs

        integer qc_flag(max_modes,max_levels)
        character*4 c4_qc
        real level(max_modes,max_levels)

        character*2  ngates_short(max_modes)
        character tmpgates    (max_modes*4)
        integer   ngates      (max_modes)


!       equivalence(ngates,tmpgates)

        real ht_out(max_levels_out)
        real di_out(max_levels_out)
        real sp_out(max_levels_out)

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
        integer wsmr_wmo_id
        data wsmr_wmo_id/74533/
        integer error_code
        data error_code/1/
        logical l_in_box
        data l_in_box/.true./

        integer varid
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

        character*13 filename13,outfile
        character*9 asc9_tim,a9time_ob

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
            write(6,*)' warning: bad set_error ',status
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

!       dir_in = path_to_raw_blpprofiler

        c_vars_req = 'path_to_raw_blpprofiler'
        call get_static_info(c_vars_req,c_values_req,1,istatus)
        if(istatus .eq. 1)then
            write(6,*)c_vars_req(1:30),' = ',c_values_req
            dir_in = c_values_req
        else
            write(6,*)' error getting ',c_vars_req
            return
        endif


        call s_len(dir_in,len_dir_in)

c       wait for the data

!       c_filespec = dir_in(1:len_dir_in)//'*0100o'
        c_filespec = dir_in(1:len_dir_in)//'*0??0o'

        write(6,*)c_filespec(1:80)

        i4time_desired = i4time_sys

        i4time_now = i4time_now_gg()
        i4_hour = (i4time_now/3600) * 3600
        minutes_now = (i4time_now - i4_hour) / 60
        i4time_stop_waiting = i4time_sys + 26 * 60
        i4_wait_period = i4time_stop_waiting - i4time_now

        i4_check_interval = 10
        i4_total_wait = min(300,i4_wait_period)
        i4_thresh_age = 3600

        if(minutes_now .ge. 19 .and. minutes_now .lt. 26)then
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
        endif ! minutes_now

c       read in the raw profiler data

        fnam_in = dir_in(1:len_dir_in)//asc9_tim//'0100o'
!       fnam_in = dir_in(1:len_dir_in)//asc9_tim//'0??0o'
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

        call prof_i4_avg_wdw(i4_avg_wdw_sec,cdfid,istatus)
!       this is hard-wired until 'prof_i4_avg_wdw' can be made more
!       general for blp profilers.
            
!       istatus = 1
!       i4_avg_wdw_sec = 3600

        if(istatus .eq. 1)then
            lag_time = i4_avg_wdw_sec/2
        else
            write(6,*)' ingest_sub_blppro: '
     1                   ,'error obtaining i4_avg_wdw_sec'
            return
        endif

!       get the number of profilers from the netcdf file
        status = nf_inq_dimid(cdfid,'recnum',varid)
        status = nf_inq_dim(cdfid,varid,dimname,file_n_prof)
        write(6,*)' # of profilers = ',file_n_prof
        if (file_n_prof .gt. n_profilers) then
          write(6,*)' error: too many profilers to process'
          istatus = 0
          return
        endif

!	get the number of characters in the station name
        status = nf_inq_dimid(cdfid,'stanamlen',varid)
        status = nf_inq_dim(cdfid,varid,dimname,stanamlen)

        start(1) = 1
        count(2) = 1
        count(1) = stanamlen

!       do ista = 1,n_profilers

        do ista = 1,file_n_prof

!       fill profiler name into prof_name(ista)
          status = nf_inq_varid(cdfid,'staname',varid)
          start(2) = ista
          status = nf_get_vara_text(cdfid, varid, start, count,
     1                prof_name(ista))
          prof_name(ista)(6:6) = ' '
          

          write(6,*)
          write(6,*)' looping for profiler ',prof_name(ista)
c
          if(.false.)then
              nmodes = 0
	      call prof_cdf_read(cdfid,prof_name(ista),0
     1                          ,'nummodesused',0,nmodes_short,
     $			         status)

              if(status .ne. 0)then
                  write(6,*)prof_name(ista)
     1             ,'  warning: bad status reading nummodesused',status       
                  go to 900
              else
                  write(6,*)'nummodesused =',nmodes
              endif

	      call prof_cdf_read(cdfid,prof_name(ista),0
     1                          ,'numlevelsused',0,ngates_short,status)       
	      do i = 1, max_modes
                  ngates(i) = 0
                  do j=1,2
                      tmpgates((i-1)*4+j+2) = ngates_short(2*(i-1)+j)
                  enddo
	      enddo

          else ! this may help for big/little endian conversions
              nmodes = 0
	      call prof_cdf_read(cdfid,prof_name(ista),0
     1                          ,'nummodesused',0,nmodes,
     $			         status)

!             call short_to_i4(nmodes_short,nmodes)

              if(status .ne. 0)then
                  write(6,*)prof_name(ista)
     1             ,'  warning: bad status reading nummodesused',status       
                  go to 900
              else
                  write(6,*)'nummodesused =',nmodes
              endif

	      call prof_cdf_read(cdfid,prof_name(ista),0
     1                          ,'numlevelsused',0,ngates,status)       
!             do i = 1, max_modes
!                 ngates(i) = 0
!                 call short_to_i4(ngates_short(i),ngates(i))
!             enddo

          endif

          if(status .ne. 0)then
            write(6,*)' warning: bad status reading numlevelsused'
            go to 900
          else
            write(6,*)' ngates array = ',ngates
          endif

c
c       get the surface pressure.  this time we'll use the
c       5-character site name (plus a terminating blank) to select the station.
c
          call prof_cdf_read(cdfid,prof_name(ista),0,'pressure',2,prs
     1                      ,status)
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

            call prof_cdf_read(cdfid,prof_name(ista),0
     1                         ,'windspeedsfc',2,sp_sfc,istatus1)
            call prof_cdf_read(cdfid,prof_name(ista),0
     1                         ,'winddirsfc',2,di_sfc,istatus2)

            if(abs(sp_sfc) .gt. 500.)istatus1 = 1
            if(abs(di_sfc) .gt. 500.)istatus2 = 1

            if(istatus1 .eq. 0 .and. istatus2 .eq. 0)then
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
            call prof_cdf_read(cdfid,prof_name(ista),0,'uvqualitycode'
     $                     ,0,qc_flag,status)
            if(status.ne.0)then
                write(6,*)' warning: bad qualitycode read ',status
                return
            endif
c
c           get the associated levels for this profiler
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
!           do i = 1, max_levels

            do im = 2,1,-1

              write(6,*)
              write(6,*)' checking mode ',im,' for ',prof_name(ista)

              nlevels = ngates(im)

!             check for used levels exceeding levels as stated in netcdf file
              if(nlevels .gt. max_levels)then
                  write(6,*)' warning: too many levels in blp-profiler'
     1                     ,nlevels,max_levels
                  nlevels = max_levels
              endif

              do i = 1, nlevels
                iqc_flag = qc_flag(im,i)

                if(iqc_flag.eq.good)then
                        j = 1
                        iqc = 1
                        c4_qc = 'good'
                else if(iqc_flag.eq.bad)then
                        j = 2
                        iqc = 0
                        c4_qc = 'bad'
                else if(iqc_flag.eq.missing)then
                        j = 3
                        iqc = 0
                        c4_qc = 'msg'
                else
                        iqc = 0
                        c4_qc = 'qcd'
                endif

                height_msl = level(im,i) + elev

                mode_flag = 1

!               use mode 1 (100 m) only if it's above the highest good level 
!               from mode 2
                if(n_good_levels .ge. 1 .and. im .eq. 1)then
                    if(height_msl .le. ht_out(n_good_levels))then
                        mode_flag = 0 
                    endif
                endif

                write(6,251,err=252)i,height_msl,u(im,i),v(im,i)
     1                             ,iqc_flag,iqc,mode_flag,c4_qc
 251            format(i3,f8.0,2f8.1,1x,2i3,i3,1x,a4)
 252            continue

                if(  (.not.
     1                    (u(im,i) .gt. r_missing_data      .or. 
     1                     v(im,i) .gt. r_missing_data)  )
     1                          .and.
     1                      iqc .eq. 1
     1                          .and.
     1                           mode_flag .eq. 1
     1                                                          )then

                    n_good_levels = n_good_levels + 1
                    ht_out(n_good_levels) = height_msl
                    call uv_to_disp(u(im,i),v(im,i)
     1                ,di_out(n_good_levels),sp_out(n_good_levels))
                endif
              enddo ! i
            enddo ! im

            write(6,*)

            i4time_ob = i4_timeobs - lag_time ! i4time_sys - lag_time

            call make_fnam_lp(i4time_ob,a9time_ob,istatus)
c
c           open an output file.
c
            ext = 'pro'

            call s_len(ext,len_ext)
            call open_ext(lun_out,i4time_sys,ext(1:len_ext),istatus)
            if(istatus .ne. 1)then
                write(6,*)' error opening product file',ext(1:len_ext)
                return
            else
                write(6,*)' opened/reopened extension ',ext(1:len_ext)       
            endif

            write(*,401)wsmr_wmo_id,n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(ista),a9time_ob,'profiler'
            write(1,401)wsmr_wmo_id,n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(ista),a9time_ob,'profiler'
401         format(i12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)

            rms = 1.0

            if(n_good_sfc .eq. 1)then
!               write surface winds as first level
                write(1,301)elev,di_sfc,sp_sfc,rms ! /r_mspkt
                write(6,301)elev,di_sfc,sp_sfc,rms ! /r_mspkt
            endif

            do i = 1,n_good_levels
!           do i = n_good_levels, 1, -1
                write(1,301,err=303)ht_out(i),di_out(i),sp_out(i),rms
                write(6,301,err=303)ht_out(i),di_out(i),sp_out(i),rms
301             format(1x,f6.0,f6.0,2f6.1)
303             continue
            enddo ! i

          else !
             write(6,*)prof_name(ista),' is outside of domain'

          endif ! l_in_box

          write(6,*)prof_name(ista)

900     enddo ! ista

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


        subroutine short_to_i4(shortvar,i4var)

        character*1 shortvar(2) 
        character*1 shortvar_loc(4)
        character*1 i4_to_byte
        integer i4var, i4var_loc
        equivalence (shortvar_loc, i4var_loc)

        shortvar_loc(1) = shortvar(1)
        shortvar_loc(2) = shortvar(2)
        shortvar_loc(3) = i4_to_byte(0)
        shortvar_loc(4) = i4_to_byte(0)

        iarg1 = i4var_loc
        call in_to_im(2,4,i4var_loc,1)                 
        iarg2 = i4var_loc

        if(abs(iarg1) .le. abs(iarg2))then
            i4var = iarg1
        else
            i4var = iarg2
        endif

        write(6,*)' short_to_i4:',iarg1,iarg2,i4var

        return
        end        

