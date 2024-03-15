
        subroutine ingest_madis_map(i4time_sys,nx_l,ny_l,ext,istatus)

        integer status,max_profiles,max_levels

        parameter (max_subdirs = 1)

        character*200 fnam_in
        character*180 dir_in
        character*255 c_filespec
        logical l_exist

        character*13 a13_time,filename13,cvt_i4time_wfo_fname13,outfile       
        character*9 asc9_tim,a9time_ob

        character*3     ext
        integer       len_dir_in

        character*40 c_vars_req
        character*180 c_values_req

        character*9 a9_timeobs
        integer timeobs

        real lat(nx_l,ny_l),lon(nx_l,ny_l)
        real topo(nx_l,ny_l)

        write(6,*)' start ingest_madis_map, ext = ',ext

        lun_out = 1

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

        outfile = filename13(i4time_sys,ext)
        asc9_tim = outfile(1:9)

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

        i4_prof_window = 1800 ! could be reset to laps_cycle_time

        i4time_raw = i4time_sys

        n_good_obs = 0

c       read in the raw profiler/rass data

!       the assumption is that the stated observation time (end of averaging 
!       period) is for the 30 minute period ending at the file name....

 500    a13_time = cvt_i4time_wfo_fname13(i4time_raw)
        fnam_in = dir_in(1:len_dir_in)//a13_time
        call s_len(fnam_in,len_fnam_in)
        write(6,*)' file = ',fnam_in(1:len_fnam_in)

        inquire(file=fnam_in(1:len_fnam_in),exist=l_exist)

        if(l_exist)then

            istatus = 0

            call read_ldad_prof(i4time_sys,i4_prof_window              ! i
     1                                    ,nx_l,ny_l                   ! i
     1                                    ,ext                         ! i
     1                                    ,lun_out                     ! i
     1                                    ,fnam_in(1:len_fnam_in)      ! i
     1                                    ,n_good_obs                  ! i/o
     1                                    ,istatus)                    ! o

            if(istatus.ne.1)then
                write(6,*)' warning: bad status in ingest_madis_map'
                goto980
            endif

            write(6,*)' n_good_obs = ',n_good_obs

 980        continue

        else
                write(6,*)' warning: cannot find file '
     1                   ,fnam_in(1:len_fnam_in)

        endif ! file exists

        i4time_raw = i4time_raw - 1800

        if(abs(i4time_raw-i4time_sys) .le. 1800 .and. 
     1                     n_good_obs .eq. 0          )then
            write(6,*)' looping back for earlier file time'
            goto 500
        endif

        write(6,*)

        return
        end


