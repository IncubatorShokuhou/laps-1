
        subroutine read_pro_metadata(lun,i4time_prof,ext              ! i
     1                         ,max_pr,max_pr_levels                  ! i
     1                         ,n_profiles                            ! o
     1                         ,nlevels_obs_pr                        ! i/o
     1                         ,lat_pr,lon_pr,elev_pr                 ! o
     1                         ,c5_name,i4time_ob_pr,obstype          ! o
     1                         ,ob_pr_ht_obs                          ! o
     1                         ,istatus)                              ! o

cdoc    returns wind profile metadata from the pro file. 
cdoc    jacket for read_pro_data.

        real lat_pr(max_pr)
        real lon_pr(max_pr)
        real elev_pr(max_pr)

        real ob_pr_ht_obs(max_pr,max_pr_levels)
        real ob_pr_u_obs(max_pr,max_pr_levels)
        real ob_pr_v_obs(max_pr,max_pr_levels)
        real ob_pr_rms_obs(max_pr,max_pr_levels)

        real sfc_t(max_pr), sfc_p(max_pr), sfc_rh(max_pr)
        real sfc_u(max_pr), sfc_v(max_pr) 

        integer i4time_ob_pr(max_pr)
        integer nlevels_obs_pr(max_pr)

        character*5 c5_name(max_pr)
        character*8 obstype(max_pr)
        character ext*(*)

        call read_pro_data(     lun,i4time_prof,ext                   ! i
     1                         ,max_pr,max_pr_levels                  ! i
     1                         ,n_profiles                            ! o
     1                         ,nlevels_obs_pr                        ! i/o
     1                         ,lat_pr,lon_pr,elev_pr                 ! o
     1                         ,c5_name,i4time_ob_pr,obstype          ! o
     1                         ,ob_pr_ht_obs                          ! o
     1                         ,ob_pr_u_obs,ob_pr_v_obs               ! o
     1                         ,ob_pr_rms_obs                         ! o
     1                         ,sfc_t,sfc_p,sfc_rh,sfc_u,sfc_v        ! o
     1                         ,istatus)                              ! o

        return
        end 


        subroutine read_pro_data(lun,i4time_prof,ext                  ! i
     1                         ,max_pr,max_pr_levels                  ! i
     1                         ,n_profiles                            ! o
     1                         ,nlevels_obs_pr                        ! i/o
     1                         ,lat_pr,lon_pr,elev_pr                 ! o
     1                         ,c5_name,i4time_ob_pr,obstype          ! o
     1                         ,ob_pr_ht_obs                          ! o
     1                         ,ob_pr_u_obs,ob_pr_v_obs               ! o
     1                         ,ob_pr_rms_obs                         ! o
     1                         ,sfc_t,sfc_p,sfc_rh,sfc_u,sfc_v        ! o
     1                         ,istatus)                              ! o

cdoc    returns wind profile data from the pro file

!       profile stuff
        real lat_pr(max_pr)
        real lon_pr(max_pr)
        real elev_pr(max_pr)

        real ob_pr_ht_obs(max_pr,max_pr_levels)
        real ob_pr_di_obs(max_pr_levels)
        real ob_pr_sp_obs(max_pr_levels)
        real ob_pr_u_obs(max_pr,max_pr_levels) ! includes sfc wind when valid
        real ob_pr_v_obs(max_pr,max_pr_levels) ! includes sfc wind when valid
        real ob_pr_rms_obs(max_pr,max_pr_levels)

!       note that surface data access is still under development and may not
!       be reliable yet
        real sfc_t(max_pr), sfc_p(max_pr), sfc_rh(max_pr)
        real sfc_u(max_pr), sfc_v(max_pr) 

        integer i4time_ob_pr(max_pr)
        integer nlevels_obs_pr(max_pr)  ! includes sfc wind when valid
                                        ! optionally can be initialized to 0.
                                        ! by the calling program even though 
                                        ! this isn't necessary for the proper
                                        ! operation of this routine
        character*5 c5_name(max_pr)
        character*8 obstype(max_pr)
        character ext*(*)
        character*9 a9time_ob
        character*132 line

        logical l_sfc, l_goodwind

        istatus = 0

        call get_sfc_badflag(badflag,istat_badflag)

        rcycles_pr = 0

        call open_lapsprd_file_read(lun,i4time_prof,ext,istatus)
        if(istatus .ne. 1)go to 420

        do i_pr = 1,max_pr

            read(lun,401,err=430,end=450)
     1           ista,nlevels_in,lat_pr(i_pr),lon_pr(i_pr)
     1          ,elev_pr(i_pr),c5_name(i_pr),a9time_ob,obstype(i_pr)       
401         format(i12,i12,f11.0,f15.0,f15.0,5x,a5,1x,3x,a9,1x,a8)

            call i4time_fname_lp(a9time_ob,i4time_ob_pr(i_pr),istatus)

            write(6,407)i_pr,ista,nlevels_in,lat_pr(i_pr)
     1                 ,lon_pr(i_pr)
     1                 ,elev_pr(i_pr),rcycles_pr,c5_name(i_pr)
     1                 ,a9time_ob,obstype(i_pr)
407         format(/' profile #',i3,i6,i5,2f8.2,e10.3,f8.2,1x,a6,3x,a9
     1                          ,1x,a8)

            level_out = 0

            do level = 1,nlevels_in
                read(lun,402,end=430)line
 402            format(a)

                read(line,*,err=410,end=410)                 ! rd sfc if there
     1                              ob_pr_ht_obs_in 
     1                             ,ob_pr_di_obs_in
     1                             ,ob_pr_sp_obs_in
     1                             ,sfc_p_in,sfc_t_in,sfc_rh_in

                sfc_p(i_pr) = sfc_p_in
                sfc_t(i_pr) = sfc_t_in
                sfc_rh(i_pr) = sfc_rh_in
                
                ob_pr_rms_obs_in = 1.0 ! hardwired until available in pro file

                l_sfc = .true.

                if(level .gt. 1)then
                    write(6,*)' error: sfc data exists past level 1'       
                    istatus = 0
                    return
                endif

                if( abs(ob_pr_ht_obs_in - elev_pr(i_pr)) .gt. 1.)then
                    write(6,*)' error: elevation disagrees with sfc'     
     1                       ,ob_pr_ht_obs_in,elev_pr(i_pr)
                    istatus = 0
                    return
                endif                    

                goto415

 410            read(line,*,err=430)ob_pr_ht_obs_in          ! rd upr lvl only
     1                             ,ob_pr_di_obs_in
     1                             ,ob_pr_sp_obs_in
     1                             ,ob_pr_rms_obs_in
                l_sfc = .false.

                if( abs(ob_pr_ht_obs_in - elev_pr(i_pr)) .lt. 1.)then
                    write(6,*)' note: elevation agrees with sfc'
     1                       ,ob_pr_ht_obs_in,elev_pr(i_pr)

                    if(level .gt. 1)then
                        write(6,*)' error: elevation agrees with sfc'
     1                       ,ob_pr_ht_obs_in,elev_pr(i_pr)
                        istatus = 0
                        return
                    else
                        write(6,*)' note: no p, t, rh present at sfc'       
                        sfc_p(i_pr) = badflag
                        sfc_t(i_pr) = badflag
                        sfc_rh(i_pr) = badflag
                    endif

                endif ! gate elevation matches the surface

 415            if(    ob_pr_di_obs_in .ne. badflag 
     1           .and. ob_pr_sp_obs_in .ne. badflag )then    ! good wind
                    level_out = level_out + 1

                    if(level_out .gt. max_pr_levels)then
                        write(6,*)
     1                      ' error: too many profiler (.pro) levels '
     1                      ,i_pr,level_out,max_pr_levels
                        goto430
                    endif

                    ob_pr_ht_obs(i_pr,level_out)  = ob_pr_ht_obs_in 
                    ob_pr_di_obs(level_out)       = ob_pr_di_obs_in 
                    ob_pr_sp_obs(level_out)       = ob_pr_sp_obs_in 
                    ob_pr_rms_obs(i_pr,level_out) = ob_pr_rms_obs_in 

                    nlevels_obs_pr(i_pr) = level_out

                    call disp_to_uv(ob_pr_di_obs(level_out),
     1                              ob_pr_sp_obs(level_out),
     1                              ob_pr_u_obs(i_pr,level_out),
     1                              ob_pr_v_obs(i_pr,level_out))

                    if(l_sfc)then
                        sfc_u(i_pr) = ob_pr_u_obs(i_pr,level_out)
                        sfc_v(i_pr) = ob_pr_v_obs(i_pr,level_out)
                    endif

                else                                         ! bad wind
                    if(l_sfc)then
                        sfc_u(i_pr) = badflag
                        sfc_v(i_pr) = badflag
                    else
                        write(6,*)' error: non-sfc wind should be good'       
                        istatus = 0
                        return
                    endif

                endif                

311             format(1x,2i4,5f8.1)
312             continue
            enddo ! level
        enddo ! profiler site

        write(6,*)' error: # of profilers has reached'
     1           ,' dimensions of max_pr ',max_pr

        n_profiles=max_pr
        go to 500
c
c       profiler reading error handling
c
  420   write(6,*) ' warning, could not open pro file'
        n_profiles=0
        go to 500

  430   write(6,*) ' error during read of pro file'
        write(6,*) ' while reading profiler number ',i_pr
        n_profiles=i_pr-1
        go to 500
c
  450   continue ! used for end of file
        n_profiles=i_pr-1

        istatus = 1

        close(lun)

  500   continue

        return
        end








