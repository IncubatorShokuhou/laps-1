
        subroutine get_interval_precip(c_pcp,ext_in
     1                                ,i4time_start,i4time_end       
     1                                ,laps_cycle_time,i4_model_init
     1                                ,r_missing_data
     1                                ,ni,nj,nf,pcp_2d,istatus)

!       works for analysis, forecast, or stage iv                        
!       special exception for 'nam-nh'  

        real pcp_2d(ni,nj,nf)
        real pcp_buf_2d(ni,nj,nf)
        character*3 var_2d(nf)
        character*(*)ext_in
        character*31 ext
        character*1 c_pcp ! used only if a single field is chosen
        character*10  units_2d(nf)
        character*125 comment_2d(nf)
        character*9 a9_start,a9_end
        character*9 a9_init,a9_valid
        character*10 c_model
        character*255 directory
        character*4 lvl_coord_2d

        ext = ext_in(1:31)

        write(6,*)' get_interval_precip ext_in/ext: '
     1                     ,trim(ext_in),' ',trim(ext)

        call make_fnam_lp(i4time_start,a9_start,istatus)   
        call make_fnam_lp(i4time_end  ,a9_end  ,istatus)   

        pcp_2d = 0.
        c_model = ' '

        if(trim(ext) .eq. 'st4')then
            var_2d(1) = 'ppt'       
            ipcp_cycle_time = 3600         
        elseif(nf .eq. 2)then
            var_2d(1) = 'r01'
            var_2d(2) = 's01'
            ipcp_cycle_time = laps_cycle_time
        elseif(trim(ext) .eq. 'l1s')then
            var_2d(1) = c_pcp(1:1)//'01'
            ipcp_cycle_time = laps_cycle_time 
        else ! i.e. fsf
            c_model = ext
            ext = 'fsf'
            var_2d(1) = c_pcp(1:1)//'01'
            ipcp_cycle_time = laps_cycle_time ! (e.g. from model_fcst_intvl)
        endif

        i4_loop_start = i4time_start+ipcp_cycle_time
        i4_loop_end   = i4time_end
        i4_loop_int   = ipcp_cycle_time

!       exception for 'nam-nh' (precip convention depends on file time)
        if(trim(c_model) .eq. 'nam-nh')then
            i4_fcst_end = i4time_end - i4_model_init
!           6h fcst ending at 00/12z for initialization at 00/12z
            if(mod(i4_fcst_end,43200) .eq. 0 .and. 
     1         mod(i4_model_init,43200) .eq. 0 .and.
     1         i4time_end-i4time_start .eq. 21600)then
                write(6,*)' nam-nh exception 1'
                i4_loop_start = i4time_start
                i4_loop_end   = i4time_end
                i4_loop_int   = ipcp_cycle_time
            else
                write(6,*)' nam-nh situation pass through..'
                write(6,*)' i4_model_init = ',i4_model_init
                write(6,*)' i4_fcst_end   = ',i4_fcst_end
                write(6,*)' i4time_end    = ',i4time_end
                write(6,*)' i4time_start  = ',i4time_start
!               istatus = 0
!               return
            endif
        endif

        write(6,*)' get_interval_precip ext_in/ext/var: '
     1                     ,trim(ext_in),' ',trim(ext),' ',var_2d(1)
        write(6,*)' get_interval_precip times: ',a9_start,' ',a9_end

        i4time_interval = i4time_end - i4time_start
c    abdel       
        if (ipcp_cycle_time.eq.0) i4time_interval = 0
	    
        if(i4time_interval .ne. 
     1       (i4time_interval/ipcp_cycle_time) *ipcp_cycle_time)then
            write(6,*)' error: interval is not an integral # of cycles'
            pcp_2d = r_missing_data
            istatus = 0
            return
        endif

        do i4time = i4_loop_start,i4_loop_end,i4_loop_int

!           if(nf .le. 2)then
            if(.false.)then
                call get_laps_multi_2d(i4time,ext,var_2d
     1                                ,units_2d,comment_2d,ni,nj,nf
     1                                ,pcp_buf_2d,istatus_file)
            elseif(trim(ext) .eq. 'fsf')then
                i4_valid = i4_time
                call get_directory(ext,directory,len_dir)
                call s_len(c_model,len_model)
                directory=
     1              directory(1:len_dir)//c_model(1:len_model)//'/'

                level = 0
                call make_fnam_lp(i4_model_init,a9_init,istatus)
                call make_fnam_lp(i4time,a9_valid,istatus)
                write(6,*)' call read_laps for pcp '
     1                        ,i4_model_init,i4time
     1                        ,trim(ext),' ',a9_init,' ',a9_valid
                call read_laps(i4_model_init,i4time,directory,
     1                        ext,ni,nj,1,1,       
     1                        var_2d,level,lvl_coord_2d,
     1                        units_2d,comment_2d,
     1                        pcp_buf_2d,istatus_file)

            else ! this can be turned on for testing if needed
                write(6,*)' call get_laps_2d for pcp'
                call get_laps_2d(i4time,ext,var_2d
     1                          ,units_2d,comment_2d,ni,nj
     1                          ,pcp_buf_2d,istatus_file)
                write(6,*)' istatus_file from get_laps_2d:',istatus_file

            endif
            if(istatus_file .ne. 1 .and. istatus_file .ne. -1)then
                write(6,*)' error: missing precip file: ',istatus_file 
                if(trim(ext) .eq. 'fsf')then
                    write(6,*)' directory = ',trim(directory)
                    write(6,*)' ext = ',trim(ext)
                endif
                pcp_2d = r_missing_data
                istatus = 0
                return
            endif                

            do i = 1,nf
                call array_range(pcp_buf_2d,ni,nj,rmin,rmax
     1                          ,r_missing_data)
                write(6,*)' get_interval_precip range: ',nf,rmin,rmax 
            enddo ! i

            do i = 1,nf
                if(trim(c_model) .eq. 'nam-nh')then
                    if(mod(i4_fcst_end,43200)   .eq. 0  .and. 
     1                 mod(i4_model_init,43200) .eq. 0  .and.
     1                 i4time_end-i4time_start  .eq. 21600 .and.
     1                 i4time .eq. i4time_start                )then
                        write(6,*)' nam-nh exception 2: negate precip'
                        pcp_buf_2d(:,:,i) = -pcp_buf_2d(:,:,i)
                    endif
                endif
                call add_miss(pcp_2d(1,1,i),pcp_buf_2d(1,1,i)
     1                       ,pcp_2d(1,1,i),ni,nj)
            enddo ! i

        enddo ! i4time

        do i = 1,nf
            call array_range(pcp_2d(1,1,i),ni,nj,rmin,rmax
     1                                          ,r_missing_data)
            write(6,*)' get_interval_precip total range: ',i,rmin,rmax
        enddo ! i

        return
        end
