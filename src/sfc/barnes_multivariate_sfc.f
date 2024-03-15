
        subroutine barnes_multivariate_sfc_jacket(c_field,obs,mxstn   ! i
     1                                           ,tb                  ! i
     1                                           ,badflag,ni,nj       ! i
     1                                           ,rms_thresh_norm     ! i
     1                                           ,bad_mult_land       ! i
     1                                           ,bad_mult_water      ! i
     1                                           ,topo,ldf            ! i
     1                                           ,wt_bkg_a            ! i
     1                                           ,t_2d,istatus)       ! o

!       this subroutine can be a substitute to the call to 'spline' and 
!       interfaces more directly with 'barnes_multivariate_sfc'. it also
!       makes more complete use the data structures for storing the obs.

        character*(*) c_field

!       input observation data structure (all variables)
        include 'sfcob.inc'
        type (sfcob) obs(mxstn)

!       second observation data structure (just one variable)
        include 'barnesob.inc'
        type (barnesob_qc) obs_barnes(mxstn)

        real tb(ni,nj)                                ! background field
        real t_2d(ni,nj)                              ! analyzed field
        real topo(ni,nj),ldf(ni,nj)                   ! topo & landfrac
        real ob_diff(mxstn)
        real ob_full(mxstn)
        real ob_bkg(mxstn)
        real to_2d_dum(ni,nj)
        real wt_bkg_a(ni,nj)                         

        write(6,*)' subroutine barnes_multivariate_sfc_jacket for...'
     1           ,c_field      

!       transfer from 'obs' data structure to 'obs_barnes' data structure
        if(c_field .eq. 'tgd' .or. .true.)then
            rinst_err = 1.5 ! set eventually according to expected accuracy
        endif

        call get_sfcob_field(obs,mxstn,'tgd',ob_full,istatus)

        call get_systime_i4(i4time_sys,istatus) ! temporary while no time wt

        n_obs = 0
        do i = 1,mxstn
            n_obs = n_obs + 1

            ista = obs(i)%i
            jsta = obs(i)%j
            obs_barnes(i)%i = ista
            obs_barnes(i)%j = jsta

            obs_barnes(i)%k = 1
            obs_barnes(i)%weight = 1. / rinst_err**2
            obs_barnes(i)%vert_rad_rat = 1.
            obs_barnes(i)%i4time = i4time_sys ! no time weight for now

            obs_barnes(i)%ldf = obs(i)%ldf
c
c th: 29 november 2002 begin hack.
c
            if ((c_field .eq. 'p') .or. (c_field .eq. 'msl')) then
               obs_barnes(i)%mask_sea = 0
            else
               obs_barnes(i)%mask_sea = 1
            endif
c
c th: end hack.
c

            if(obs(i)%sfct_f .ne. badflag)then
                obs_barnes(i)%qc = .true.
            else
                obs_barnes(i)%qc = .false.
            endif

!           reject obs (for now) if they are outside of the domain
            if(ista .lt. 1 .or. ista .gt. ni 
     1    .or. jsta .lt. 1 .or. jsta .gt. nj )then
                obs_barnes(i)%qc = .false.
            endif

!           calculate difference of ob from the background
            if(obs_barnes(i)%qc)then
                ob_bkg(i)  = tb(ista,jsta)
                ob_diff(i) = ob_full(i) - ob_bkg(i)
                obs_barnes(i)%value(1) = ob_diff(i)
            endif

        enddo ! i

!       calculate rms (stdev) of the ob-background differences
        cnt = 0
        sumsq = 0.
        do i = 1,mxstn
            if(obs_barnes(i)%qc)then
                cnt = cnt + 1.
                sumsq = sumsq + ob_diff(i)**2
            endif
        enddo ! i

        if(cnt .gt. 0.)then
            std = sqrt(sumsq / cnt)
        else
            std = 0.
        endif

        bad_land = bad_mult_land * std
	print *,' std dev: ',std,', bad value: ',bad_land
     1                          ,', ratio: ',bad_mult_land      

        bad_water = bad_mult_water * std
	print *,' std dev: ',std,', bad value: ',bad_water
     1                          ,', ratio: ',bad_mult_water      

!       eliminate bad data from the data structure
        do i = 1,mxstn
            if(obs_barnes(i)%qc)then

                bad =        obs_barnes(i)%ldf  * bad_land + 
     1                (1.0 - obs_barnes(i)%ldf) * bad_water

                if(ob_diff(i) .ge. bad)then
                    obs_barnes(i)%qc = .false.
                    write(6,1099) obs_barnes(i)%i, obs_barnes(i)%j
     1                          , ob_full(i), ob_diff(i)
1099	            format(1x,'bad data at i,j ',2i5,': value ',e12.4       
     1                    ,', diff ',e12.4)
                endif
            endif
        enddo ! i

        call get_fnorm_max(ni,nj,r0_norm,r0_value_min,fnorm_max)   
        n_fnorm = int(fnorm_max) + 1

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

!       if(.true.)return ! temporary

!       note that we do not need to fill to_2d_dum when l_struct = .true.
!       note that we don't need to fill l_boundary when l_struct = .true.
!                                                  and  l_barnes_wide = .false.

        call barnes_multivariate_sfc(to_2d_dum,ni,nj             ! inputs
     1                   ,r_missing_data                         ! input
     1                   ,rms_thresh_norm                        ! input
     1                   ,rinst_err                              ! input
     1                   ,bad                                    ! input
     1                   ,wt_bkg_a                               ! input
     1                   ,n_fnorm                                ! input
     1                   ,l_boundary,.false.,.true.              ! input
     1                   ,mxstn,obs_barnes                       ! input
     1                   ,topo,ldf                               ! input
     1                   ,t_2d                                   ! output
     1                   ,istatus)                               ! output

        write(6,*)' adding incremental analysis to background'       
        call add(t_2d,tb,t_2d,ni,nj)

        return
        end

        subroutine barnes_multivariate_sfc(to_2d_in,ni,nj        ! inputs
     1                   ,r_missing_data                         ! input
     1                   ,rms_thresh_norm                        ! input
     1                   ,rinst_err                              ! input
     1                   ,bad                                    ! input
     1                   ,wt_bkg_a                               ! input
     1                   ,n_fnorm                                ! input
     1                   ,l_boundary,l_barnes_wide,l_struct_in   ! input
     1                   ,n_obs,obs_barnes                       ! input
     1                   ,topo,ldf                               ! input
     1                   ,t_2d                                   ! output
     1                   ,istatus)                               ! output

        include 'barnesob.inc'
        type (barnesob_qc) obs_barnes(n_obs)
        type (barnesob)    obs_barnes_valid(n_obs)

        integer nk
        parameter (nk=1)

        real t_2d(ni,nj)                              ! analyzed field
        real to_2d_in(ni,nj)                          ! observations
        real wt_bkg_a(ni,nj)                         
        real topo(ni,nj),ldf(ni,nj)                   ! topo & landfrac

        real wt_2d(ni,nj)
        integer n_obs_lvl

        logical l_analyze(nk), l_use_ob, l_boundary(ni,nj)

        logical l_barnes_wide          ! use 'barnes_wide' routine on boundary?
!       data l_barnes_wide /.true./    

        logical l_struct_in,l_struct_out,l_not_struct_out ! use data structures?
        data l_struct_out /.true./     ! use data structures?

        logical limit_analysis_increment
        data limit_analysis_increment /.true./

        integer  n_fnorm
        dimension fnorm(0:n_fnorm)

        write(6,*)' subroutine barnes_univariate_sfc'

        i4_elapsed = ishow_timer()

        call get_grid_spacing_cen(grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)return

        call get_systime_i4(i4time_sys,istatus) ! temporary while no time wt

        l_analyze = .false.             ! array set
        wt_2d = 0.                      ! array set

        write(6,*)' setting observation to array weights'

        sumsq_inst = 0.
        n_obs_valid = 0
        ibound = 0
        iskip_bound = 10 ! how far to skip between boundary obs that get used.
                         ! 1 means use all boundary obs, 2 means use every
                         ! other boundary ob, 3 means use every 3rd boundary ob.

        boundary_err = rinst_err
            
        if(l_struct_in)then

!         determine if point would represent a boundary ob that gets skipped
          if(l_barnes_wide)then
            do iob = 1,n_obs
                l_use_ob = .true.
                if(obs_barnes(iob)%qc)then
                    i = obs_barnes(iob)%i
                    j = obs_barnes(iob)%j              
                    if(l_boundary(i,j))then
                        ibound = ibound + 1
                        if(ibound .ne. (ibound/iskip_bound)*iskip_bound
     1                                                            )then
                            l_use_ob = .false.
                        endif
                    endif
                endif
            enddo ! iob
          endif ! l_barnes_wide

!         place valid obs from obs_barnes into obs_barnes_valid
          do iob = 1,n_obs
            if(obs_barnes(iob)%qc)then
                n_obs_valid = n_obs_valid + 1
                obs_barnes_valid(n_obs_valid)%i = obs_barnes(iob)%i
                obs_barnes_valid(n_obs_valid)%j = obs_barnes(iob)%j
                obs_barnes_valid(n_obs_valid)%k = obs_barnes(iob)%k
                obs_barnes_valid(n_obs_valid)%weight 
     1                                    = obs_barnes(iob)%weight
                obs_barnes_valid(n_obs_valid)%vert_rad_rat 
     1                                    = obs_barnes(iob)%vert_rad_rat       
                obs_barnes_valid(n_obs_valid)%i4time
     1                                    = obs_barnes(iob)%i4time
                obs_barnes_valid(n_obs_valid)%value(1) 
     1                                    = obs_barnes(iob)%value(1)
                sumsq_inst = sumsq_inst 
     1                     + 1. / obs_barnes_valid(n_obs_valid)%weight

                obs_barnes_valid(n_obs_valid)%elev 
     1                                        = obs_barnes(iob)%elev
                obs_barnes_valid(n_obs_valid)%ldf = obs_barnes(iob)%ldf
c
c th: 29 november 2002 begin hack
c
                obs_barnes_valid(n_obs_valid)%mask_sea =
     1             obs_barnes(iob)%mask_sea
c
c th: end hack.
c

            endif ! valid ob
          enddo ! iob

        else ! .not. l_struct_in
          do i = 1,ni
          do j = 1,nj
            l_use_ob = .true.

!           determine if point would represent a boundary ob that gets skipped
            if(l_barnes_wide)then
                if(l_boundary(i,j))then
                    ibound = ibound + 1
                    if(ibound .ne. (ibound/iskip_bound)*iskip_bound)then
                        l_use_ob = .false.
                    endif
                endif
            endif

            if(to_2d_in(i,j) .ne. 0.0 .and. l_use_ob)then
                wt_2d(i,j) = 1.0 / rinst_err**2 ! set differently for boundary?
                sumsq_inst = sumsq_inst + 1./wt_2d(i,j)
                n_obs_valid = n_obs_valid + 1

!               fill data structure element
                obs_barnes_valid(n_obs_valid)%i = i
                obs_barnes_valid(n_obs_valid)%j = j
                obs_barnes_valid(n_obs_valid)%k = 1
                obs_barnes_valid(n_obs_valid)%weight=1.0 / rinst_err**2
                obs_barnes_valid(n_obs_valid)%vert_rad_rat = 1.0
                obs_barnes_valid(n_obs_valid)%i4time = i4time_sys ! no time wt
                obs_barnes_valid(n_obs_valid)%value(1) = to_2d_in(i,j)
                obs_barnes_valid(n_obs_valid)%elev = topo(i,j)
                obs_barnes_valid(n_obs_valid)%ldf = ldf(i,j)
c
c th: 29 november 2002 begin hack.
c
                obs_barnes_valid(n_obs_valid)%mask_sea =
     1             obs_barnes(1)%mask_sea 
c
c th: end hack.
c

            endif

          enddo ! j
          enddo ! i
        endif

        if(n_obs_valid .gt. 0)then
            rms_inst = sqrt(sumsq_inst/float(n_obs_valid))
        else
            rms_inst = 0.
        endif

        write(6,*)' n_obs_valid,rms_inst = ',n_obs_valid,rms_inst       

!       set the rms threshold for iteration cutoff, based on instrument error
        rms_thresh = rms_thresh_norm * rms_inst

        write(6,*)'rms_thresh_norm,rms_thresh'
     1            ,rms_thresh_norm,rms_thresh      

        l_not_struct_out = .not. l_struct_out
        n_var = 1
        rep_pres_intvl = 5000. ! hardwire should work for a 2-d analysis

        r0_barnes_max_m = 140000.
        barnes_conv_rate = 0.5

        call barnes_multivariate(
     1                      t_2d                                  ! outputs
     1                     ,n_var,n_obs_valid,obs_barnes_valid    ! input
     1                     ,ni,nj,nk,grid_spacing_cen_m           ! inputs
     1                     ,rep_pres_intvl                        ! input
     1                     ,wt_bkg_a                              ! input
     1                     ,i4_loop_total                         ! output
     1                     ,wt_2d                                 ! not used
     1                     ,fnorm,n_fnorm                         ! inputs
     1                     ,l_analyze,l_not_struct_out,rms_thresh ! input
     1                     ,r0_barnes_max_m,barnes_conv_rate      ! input
     1                     ,topo,ldf,ni,nj                        ! input
     1                     ,n_obs_lvl,istatus)                    ! outputs

        call stats(t_2d,ni,nj)

        write(6,*)' limit of allowed increment should be ',bad

        if(limit_analysis_increment)then
            n_limited = 0
            do i = 1,ni
            do j = 1,nj
                if(t_2d(i,j) .lt. -bad)then
                    t_2d(i,j) = -bad
                    n_limited = n_limited + 1
                endif

                if(t_2d(i,j) .gt. +bad)then
                    t_2d(i,j) = +bad
                    n_limited = n_limited + 1
                endif
            enddo ! j
            enddo ! i

            write(6,*)' number of points with increment limited = '
     1               ,n_limited

        endif

        return
        end


        subroutine get_sfcob_field(obs,mxstn,c_field,ob_1d,istatus)

        character*(*) c_field
        real ob_1d(mxstn)

!       input observation data structure (all variables)
        include 'sfcob.inc'
        type (sfcob) obs(mxstn)

        if(c_field .eq. 'sfct_f' .or. c_field .eq. 'tgd')then
            do i = 1,mxstn
                ob_1d(i) = obs(i)%sfct_f
            enddo ! i

        elseif(c_field .eq. 't_f')then
            do i = 1,mxstn
                ob_1d(i) = obs(i)%t_f
            enddo ! i

        elseif(c_field .eq. 't_ea_f')then
            do i = 1,mxstn
                ob_1d(i) = obs(i)%t_ea_f
            enddo ! i

        elseif(c_field .eq. 'td_f')then
            do i = 1,mxstn
                ob_1d(i) = obs(i)%td_f
            enddo ! i

        elseif(c_field .eq. 'td_ea_f')then
            do i = 1,mxstn
                ob_1d(i) = obs(i)%td_ea_f
            enddo ! i

        else
            write(6,*)' error, unknown field in get_sfcob_field'
            istatus = 0
            return
        endif

        istatus = 1

        return
        end
