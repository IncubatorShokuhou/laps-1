
        subroutine precip_barnes_jacket(          c_field             ! i
     1                                           ,ilaps,jlaps         ! i
     1                                           ,pcp_gauge           ! i
     1                                           ,maxsta              ! i
     1                                           ,pcp_bkg_in          ! i
     1                                           ,badflag,ni,nj       ! i
     1                                           ,topo,ldf            ! i
     1                                           ,wt_bkg_a            ! i
     1                                           ,pcp_2d_in,istatus)  ! o

        include 'constants.inc'
        character*(*) c_field

        include 'barnesob.inc' ! observation data structure to be passed onward
        type (barnesob_qc) obs_barnes(maxsta)

        integer ilaps(maxsta),jlaps(maxsta)

        real pcp_bkg_in(ni,nj)                        ! background field
        real pcp_2d_in(ni,nj)                         ! analyzed field (in)
        real topo(ni,nj),ldf(ni,nj)                   ! topo & landfrac
        real pcp_gauge(maxsta)
        real ob_bkg(maxsta)
        real ob_full(maxsta)
        real ob_diff(maxsta)
        real wt_bkg_a(ni,nj)                         

        logical l_boundary(ni,nj)

        write(6,*)' subroutine precip_barnes_jacket for...'
     1           ,c_field      

        call get_systime_i4(i4time_sys,istatus) ! temporary while no time wt

        n_obs = 0
        do i = 1,maxsta
            n_obs = n_obs + 1

            ista = ilaps(i)
            jsta = jlaps(i)
            obs_barnes(i)%i = ista
            obs_barnes(i)%j = jsta

            ob_full(i) = pcp_gauge(i) ! inches
            if(pcp_gauge(i) .ge. 0.)then
                obs_barnes(i)%qc = .true.
            else
                obs_barnes(i)%qc = .false.
            endif

            obs_barnes(i)%k = 1
            obs_barnes(i)%weight = 1. ! / rinst_err**2
            obs_barnes(i)%i4time = i4time_sys ! no time weight for now

!           reject obs (for now) if they are outside of the domain
            if(ista .lt. 1 .or. ista .gt. ni 
     1    .or. jsta .lt. 1 .or. jsta .gt. nj )then
                obs_barnes(i)%qc = .false.
            endif

!           calculate difference of ob from the background
            if(obs_barnes(i)%qc)then
                ob_bkg(i)  = pcp_bkg_in(ista,jsta)
                ob_diff(i) = ob_full(i) - ob_bkg(i)
                obs_barnes(i)%value(1) = ob_diff(i)
            endif

        enddo ! i

!       calculate rms (stdev) of the ob-background differences
        cnt = 0
        sumsq = 0.
        do i = 1,maxsta
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


        call get_fnorm_max(ni,nj,r0_norm,r0_value_min,fnorm_max)   
        n_fnorm = int(fnorm_max) + 1

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

!       if(.true.)return ! temporary

!       note that we don't need to fill l_boundary when l_struct = .true.
!                                                  and  l_barnes_wide = .false.

        call precip_barnes_sfc(ni,nj                             ! inputs
     1                   ,r_missing_data                         ! input
     1                   ,rinst_err                              ! input
     1                   ,bad                                    ! input
     1                   ,wt_bkg_a                               ! input
     1                   ,n_fnorm                                ! input
     1                   ,l_boundary,.false.,.true.              ! input
     1                   ,maxsta,obs_barnes                      ! input
     1                   ,topo,ldf                               ! input
     1                   ,pcp_2d_in                              ! output
     1                   ,istatus)                               ! output

        write(6,*)' adding incremental analysis to background'       

        call array_range(pcp_2d_in,ni,nj,pmin,pmax,r_missing_data)
        write(6,*)' pcp_2d_in range',pmin,pmax

        call array_range(pcp_bkg_in,ni,nj,pmin,pmax,r_missing_data)
        write(6,*)' pcp_bkg_in range',pmin,pmax
        
        call add(pcp_2d_in,pcp_bkg_in,pcp_2d_in,ni,nj)

        write(6,*)' returning from precip_barnes_jacket'
        write(6,*)

        return
        end

        subroutine precip_barnes_sfc(ni,nj                       ! inputs
     1                   ,r_missing_data                         ! input
     1                   ,rinst_err                              ! input
     1                   ,bad                                    ! input
     1                   ,wt_bkg_a                               ! input
     1                   ,n_fnorm                                ! input
     1                   ,l_boundary,l_barnes_wide,l_struct_in   ! input
     1                   ,n_obs,obs_barnes                       ! input
     1                   ,topo,ldf                               ! input
     1                   ,pcp_2d_in                              ! output
     1                   ,istatus)                               ! output

        include 'barnesob.inc'
        type (barnesob_qc) obs_barnes(n_obs)
        type (barnesob)    obs_barnes_valid(n_obs)

        integer nk
        parameter (nk=1)

        real pcp_2d_in(ni,nj)                         ! analyzed field
        real wt_bkg_a(ni,nj)                         
        real topo(ni,nj),ldf(ni,nj)                   ! topo & landfrac

        real wt_2d(ni,nj)
        integer n_obs_lvl

        logical l_analyze(nk), l_use_ob, l_boundary(ni,nj)

        logical l_barnes_wide          ! use 'barnes_wide' routine on boundary?

        logical l_struct_in,l_struct_out,l_not_struct_out ! use data structures?
        data l_struct_out /.true./     ! use data structures?

        integer  n_fnorm
        dimension fnorm(0:n_fnorm)

        write(6,*)' subroutine precip_barnes_sfc'

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

        if(l_struct_in)then

!         place valid obs from obs_barnes into obs_barnes_valid
          do iob = 1,n_obs
            if(obs_barnes(iob)%qc)then
                n_obs_valid = n_obs_valid + 1
                obs_barnes_valid(n_obs_valid)%i = obs_barnes(iob)%i
                obs_barnes_valid(n_obs_valid)%j = obs_barnes(iob)%j
                obs_barnes_valid(n_obs_valid)%k = obs_barnes(iob)%k
                obs_barnes_valid(n_obs_valid)%weight 
     1                                        = obs_barnes(iob)%weight
                obs_barnes_valid(n_obs_valid)%vert_rad_rat = 1.
                obs_barnes_valid(n_obs_valid)%i4time
     1                                        = obs_barnes(iob)%i4time
                obs_barnes_valid(n_obs_valid)%value(1) 
     1                                        = obs_barnes(iob)%value(1)
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

        endif

        if(n_obs_valid .gt. 0)then
            rms_inst = sqrt(sumsq_inst/float(n_obs_valid))
        else
            rms_inst = 0.
        endif

        write(6,*)' n_obs_valid,rms_inst = ',n_obs_valid,rms_inst       

!       set the rms threshold for iteration cutoff, based on instrument error
        rms_thresh = 0.01 ! rms_inst (inches)

        write(6,*)'rms_thresh',rms_thresh      

        l_not_struct_out = .not. l_struct_out
        n_var = 1
        rep_pres_intvl = 5000. ! hardwire should work for a 2-d analysis

        r0_barnes_max_m = 140000.
        barnes_conv_rate = 0.5

!       perform gauge analysis of incremental precip in inches
        call barnes_multivariate(
     1                      pcp_2d_in                             ! outputs
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

        call stats(pcp_2d_in,ni,nj)

        return
        end


