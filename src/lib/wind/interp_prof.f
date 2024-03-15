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

        subroutine interp_prof(ob_pr_ht_obs,ob_pr_u_obs,ob_pr_v_obs, ! i
     1                             u_diff       , v_diff,            ! i
     1                             u_interp     , v_interp,          ! o
     1                             di_interp    , sp_interp,         ! o
     1                             i_pr,ht,level,nlevels_obs_pr,     ! i
     1                             lat_pr,lon_pr,i_ob,j_ob,          ! i
     1                             r_missing_data,                   ! i
     1                             heights_3d,ni,nj,nk,              ! i
     1                             max_pr,max_pr_levels,             ! i
     1                             n_vel_grids,istatus)              ! i/o

!************************arrays.for******************************************

!       profiler observations

        integer nlevels_obs_pr(max_pr)
        real ob_pr_ht_obs(max_pr,max_pr_levels)
        real ob_pr_u_obs(max_pr,max_pr_levels)
        real ob_pr_v_obs(max_pr,max_pr_levels)

!**************************** end arrays.for ********************************
        real heights_3d(ni,nj,nk)

        u_interp = r_missing_data
        v_interp = r_missing_data
        di_interp = r_missing_data
        sp_interp = r_missing_data

!  ***  interpolate the profiler observations to the input height *******
        do i_obs = 1,nlevels_obs_pr(i_pr)

          if(i_obs .gt. 1)then

            if(ob_pr_ht_obs(i_pr,i_obs-1) .le. ht .and.
     1         ob_pr_ht_obs(i_pr,i_obs  ) .ge. ht       )then

                h_lower  = ob_pr_ht_obs(i_pr,i_obs-1)
                h_upper =  ob_pr_ht_obs(i_pr,i_obs  )

                height_diff = h_upper - h_lower

                fracl = (h_upper - ht) / height_diff
                frach = 1.0 - fracl

                u_interp = ob_pr_u_obs(i_pr,i_obs) * frach
     1                    +       ob_pr_u_obs(i_pr,i_obs-1) * fracl

                v_interp = ob_pr_v_obs(i_pr,i_obs) * frach
     1                    +       ob_pr_v_obs(i_pr,i_obs-1) * fracl

!               correct for the time lag
!               u_diff      = du/dt * [t(anal) - t(ob)]
                u_interp = u_interp + u_diff

!               v_diff      = dv/dt * [t(anal) - t(ob)]
                v_interp = v_interp + v_diff

!               calculate direction and speed
                call uv_to_disp( u_interp,
     1                           v_interp,
     1                           di_interp,
     1                           sp_interp)

             endif
          endif

        enddo ! level

        if(.true.)return


!       lower tail

        if( float(level)
     1    .lt. height_to_zcoord2(ob_pr_ht_obs(i_pr,1),heights_3d
     1                          ,ni,nj,nk,i_ob,j_ob,istatus)
     1                          .and.
     1  (height_to_zcoord2(ob_pr_ht_obs(i_pr,1),heights_3d
     1                          ,ni,nj,nk,i_ob,j_ob,istatus)
     1       - float(level)) .le. 0.5)then

                u_interp  = ob_pr_u_obs(i_pr,1)
                v_interp  = ob_pr_v_obs(i_pr,1)

!               correct for the time lag
                u_interp = u_interp + u_diff
                v_interp = v_interp + v_diff

!               calculate direction and speed
                call uv_to_disp( u_interp,
     1                           v_interp,
     1                           di_interp,
     1                           sp_interp)

        endif

        if(istatus .ne. 1)then
            write(6,*)
     1        ' warning, out of bounds in height_to_zcoord2/interp_prof'       
        endif

!       upper tail

        if( float(level) .gt.
     1  height_to_zcoord2(ob_pr_ht_obs(i_pr,nlevels_obs_pr(i_pr))
     1                  ,heights_3d,ni,nj,nk,i_ob,j_ob,istatus)
     1                           .and.
     1    (float(level) -
     1     height_to_zcoord2(ob_pr_ht_obs(i_pr,nlevels_obs_pr(i_pr))
     1                  ,heights_3d,ni,nj,nk,i_ob,j_ob,istatus))
     1                                          .le. 0.5)then

                u_interp  = ob_pr_u_obs(i_pr,nlevels_obs_pr(i_pr))
                v_interp  = ob_pr_v_obs(i_pr,nlevels_obs_pr(i_pr))

!               correct for the time lag
                u_interp = u_interp + u_diff
                v_interp = v_interp + v_diff

!               calculate direction and speed
                call uv_to_disp( u_interp,
     1                           v_interp,
     1                           di_interp,
     1                           sp_interp)

        endif

        if(istatus .ne. 1)then
            write(6,*)' warning, out of bounds in height_to_zcoord2/inte
     1rp_prof'
        endif

        return

        end

