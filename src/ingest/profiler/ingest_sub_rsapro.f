
        subroutine ingest_rsapro(i4time_sys,nx_l,ny_l,lun_out,istatus)       

        integer cdfid,status,max_profiles,max_levels,file_n_prof       

!       parameter (max_profiles = 1000)
!       parameter (max_levels = 300)
        parameter (max_subdirs = 3)

!       real ht_out(max_levels)
!       real di_out(max_levels)
!       real sp_out(max_levels)

!       character*6 prof_name(max_profiles)

        character*255 prof_subdirs(max_subdirs)

!       integer i4_mid_window_pr(max_profiles)
!       integer wmo_id(max_profiles)

!       real lat_pr(max_profiles)
!       real lon_pr(max_profiles)
!       real elev_m_pr(max_profiles)
!       real n_lvls_pr(max_profiles)
!       real ht_m_pr(max_profiles,max_levels)
!       real dir_dg_pr(max_profiles,max_levels)
!       real spd_ms_pr(max_profiles,max_levels)
!       real u_std_ms_pr(max_profiles,max_levels)
!       real v_std_ms_pr(max_profiles,max_levels)

        character*200 fnam_in
        character*180 dir_in
        character*255 c_filespec
        logical l_exist

        character*13 a13_time,filename13,cvt_i4time_wfo_fname13,outfile       
        character*9 asc9_tim,a9time_ob

        character*31    ext
        integer       len_dir_in

        character*40 c_vars_req
        character*180 c_values_req

        character*9 a9_timeobs
        integer timeobs

        real lat(nx_l,ny_l),lon(nx_l,ny_l)
        real topo(nx_l,ny_l)

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

        outfile = filename13(i4time_sys,'pro')
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

        prof_subdirs(1) = '50mhz'
        prof_subdirs(2) = '915mhz'
        prof_subdirs(3) = 'minisodar'

        ext = 'pro'

        i4_prof_window = 1800 ! could be reset to laps_cycle_time

        do idir = 1,max_subdirs
            call s_len(prof_subdirs(idir),len_subdir)
 
c           read in the raw profiler data
            a13_time = cvt_i4time_wfo_fname13(i4time_sys)
            fnam_in = dir_in(1:len_dir_in)
     1                //prof_subdirs(idir)(1:len_subdir)
     1                //'/netcdf/'//a13_time
            call s_len(fnam_in,len_fnam_in)
            write(6,*)' file = ',fnam_in(1:len_fnam_in)

            inquire(file=fnam_in(1:len_fnam_in),exist=l_exist)

            if(l_exist)then

                istatus = 0

                n_good_obs = 0

                if(idir .eq. 1 .or. idir .eq. 2)then
                    call read_ldad_prof(i4time_sys,i4_prof_window      ! i
     1                                    ,nx_l,ny_l                   ! i
     1                                    ,ext                         ! i
     1                                    ,lun_out                     ! i
     1                                    ,fnam_in(1:len_fnam_in)      ! i
     1                                    ,n_good_obs                  ! i/o
     1                                    ,istatus)                    ! o

                    write(6,*)' n_good_obs = ',n_good_obs

                else ! idir .eq. 3
                    i4time_earliest = i4time_sys - i4_prof_window
                    i4time_latest   = i4time_sys + i4_prof_window

                    call get_sodar_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l        ! i
     +                   ,i4time_earliest,i4time_latest                ! i
     +                   ,fnam_in(1:len_fnam_in)                       ! i
     +                   ,lun_out                                      ! i
     +                   ,istatus)                                     ! o

                endif ! idir

                if(istatus.ne.1)then
                    write(6,*)' warning: bad status on '
     1                       ,prof_subdirs(idir),istatus           
                endif

            else
                write(6,*)' warning: cannot find file '
     1                   ,fnam_in(1:len_fnam_in)

            endif ! file exists

            write(6,*)

        enddo ! idir

        return
        end


