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

        subroutine get_modelfg(cf_modelfg,t_modelfg                      ! o
     1                        ,model_q_3d                                ! o
     1                        ,model_ref_3d                              ! o
     1                        ,default_clear_cover,r_missing_data        ! i
     1                        ,temp_3d,heights_3d,cld_hts                ! i
     1                        ,i4time_needed,ilaps_cycle_time            ! i
     1                        ,ni,nj,klaps,kcloud                        ! i
     1                        ,istatus)                                  ! o

!       obtain model first guess cloud cover and temperature fields
!       this should probably be free of very small scale horizontal structures
!       for best results when combining with the satellite data

!                steve albers   use model first guess info to guess cloud cover
!       1994     steve albers   read in sh from model background (instead of
!                                                                         rh)
!       1995 dec steve albers   qc check to prevent cf_modelfg > 1.0
!       1999     steve albers   added in lwc/ice model first guess
!                               simple usage for the moment.

        real heights_3d(ni,nj,klaps)       ! input
        real temp_3d(ni,nj,klaps)          ! input
        real cf_modelfg(ni,nj,kcloud)      ! output
        real t_modelfg(ni,nj,kcloud)       ! output
        real cld_hts(kcloud)               ! input
        real model_q_3d(ni,nj,klaps)       ! output
        real model_ref_3d(ni,nj,klaps)     ! output

        real model_lwc_3d(ni,nj,klaps)     ! local (units are mixing ratio)
        real model_ice_3d(ni,nj,klaps)     ! local (units are mixing ratio)

        real make_rh,lwc_modelfg,ice_modelfg

        real icethresh, liqthresh

        character*31 ext
        character*3 var_2d

        write(6,*)
        write(6,*)' getting first guess cloud cover'

!       mode_fg = 1 ! prefer cloud derived from model rh 
        mode_fg = 2 ! use maximum of first guess lwc and rh implied cloud
!       mode_fg = 3 ! always use first guess lwc
!       mode_fg = 4 ! set first guess to default_clear_cover

!       initialize model first guess cover field with default value
        do k = 1,kcloud
        do j = 1,nj
        do i = 1,ni
            cf_modelfg(i,j,k) = default_clear_cover
            t_modelfg(i,j,k) = 0.
        enddo
        enddo
        enddo

        i_hum_high = 0
        i_hum_low = 0
        i_hum_ok = 0

        i_condensate = 0

        istat_sh  = 0
        istat_lwc = 0
        istat_ice = 0

        write(6,*)' getting model lwc background'

!       get model first guess lwc
        var_2d = 'lwc'
        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,klaps
     1                     ,model_lwc_3d,istat_lwc)

        if(istat_lwc .ne. 1)then
            write(6,*)' no first guess available for ',var_2d
        endif

        write(6,*)' getting model ice background'

!       get model first guess ice
        var_2d = 'ice'
        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,klaps
     1                     ,model_ice_3d,istat_ice)

        if(istat_ice .ne. 1)then
            write(6,*)' no first guess available for ',var_2d
        endif

        write(6,*)' getting model sh background'

!       get model first guess sh
        var_2d = 'sh'
        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,klaps
     1                     ,model_q_3d,istat_sh)

        if(istat_sh .ne. 1)then
            write(6,*)' no first guess available for ',var_2d
            write(6,*)' cf_modelfg set to ',default_clear_cover
            write(6,*)' returning from get_modelfg with istatus = 0'
            istatus = 0
            return
        endif

        write(6,*)' getting model ref background'

!       get model first guess ref
        var_2d = 'ref'
        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,klaps
     1                     ,model_ref_3d,istat_ref)

        if(istat_ref .ne. 1)then
            write(6,*)' no first guess available for ',var_2d
            model_ref_3d = r_missing_data
        endif

!       check status of lwc/ice from model first guess
        if(istat_lwc .ne. 1 .or. istat_ice .ne. 1)then
            continue

        else ! we are using model lwc/ice fields
            call get_grid_spacing_cen(grid_spacing_cen_m,istatus)
            if(istatus .ne. 1)return

            if(grid_spacing_cen_m .le. 10000.)then  ! grid spacing in meters
                icethresh  = 0.000005 ! these are in kg kg{-1} (mixing ratio)
                snowthresh = 0.000003
                liqthresh  = 0.000003
            else
                icethresh  = 0.000005
                snowthresh = 0.000025
                liqthresh  = 0.000025
            endif

        endif ! status for model data

!       remap to cloud height grid and convert to cloud cover
        t_ref = -10. ! colder than this ice saturation is assumed

        write(6,*)' istatus: sh/lwc/ice ',istat_sh,istat_lwc,istat_ice

        do k = 1,kcloud
        do j = 1,nj
          istat_z = 0
          do i = 1,ni

!           find the model pressure at this location in the cloud height grid
            if(i-1 .eq. (i-1)/10*10)then ! update every 10th grid point
                z_laps = height_to_zcoord2(cld_hts(k),heights_3d
     1                  ,ni,nj,klaps,i,j,istat_z)

                if(istat_z .ne. 1)then
!                   we're ok if cloud height grid is above pressure grid
                    if(cld_hts(k) .le. heights_3d(i,j,klaps))then
                        write(6,*)' error: bad status from '
     1                           ,'height_to_zcoord2'
                        istatus = 0
                        return
                    endif

                else
                    z_laps = max(1.,min(z_laps,float(klaps)-.001))
                    iz_laps = int(z_laps)
                    frac = z_laps - iz_laps

                    p_modelfg = pressure_of_level(iz_laps) * (1. - frac)       
     1                        + pressure_of_level(iz_laps+1)  * frac

                    p_modelfg_mb = p_modelfg * .01

                endif ! istat_z .ne. 1

            endif

            if(istat_z .ne. 1)then
                i_grid_high = i_grid_high + 1
                cf_modelfg(i,j,k) = default_clear_cover
                go to 1000
            endif

            if(istat_lwc .eq. 1 .and. istat_ice .eq. 1)then
                lwc_modelfg =  model_lwc_3d(i,j,iz_laps)   * (1. - frac)
     1                      +  model_lwc_3d(i,j,iz_laps+1)       * frac

                ice_modelfg =  model_ice_3d(i,j,iz_laps)   * (1. - frac)
     1                      +  model_ice_3d(i,j,iz_laps+1)       * frac

                if(lwc_modelfg .ge. liqthresh .or. 
     1             ice_modelfg .ge. icethresh        )then
                    cf_modelfg(i,j,k) = 1.
                    i_condensate = i_condensate + 1
                else
                    cf_modelfg(i,j,k) = default_clear_cover
                endif

            endif

            if(istat_sh .eq. 1)then
              if(.not. (istat_lwc .eq. 1 .and. istat_ice .eq. 1) 
     1                          .or.   mode_fg .eq. 2
     1                                                           )then
!               find the model temp at this location in the cloud height grid
                t_modelfg(i,j,k) =  temp_3d(i,j,iz_laps)   * (1. - frac)      
     1                           +  temp_3d(i,j,iz_laps+1) * frac

                t_modelfg_c = t_modelfg(i,j,k) - 273.15

!               find the model sh at this location in the cloud height grid
                q_modelfg =  model_q_3d(i,j,iz_laps)    * (1. - frac)       
     1                    +  model_q_3d(i,j,iz_laps+1)  * frac

                q_modelfg_gkg = q_modelfg * 1000.

                rh_modelfg = make_rh(p_modelfg_mb              ! fractional rh
     1                     ,t_modelfg_c,q_modelfg_gkg,t_ref)

!               qc the rh
                rh_qc = rh_modelfg                             ! fractional rh

                if(rh_qc .gt. 1.0)then
                    rh_qc = 1.0
                    i_hum_high = i_hum_high + 1
                elseif(rh_qc .lt. 0.0)then
                    rh_qc = 0.0
                    i_hum_low = i_hum_low + 1
                else
                    i_hum_ok = i_hum_ok + 1
                endif

                if(cld_hts(k) .gt. 11000.)rh_qc = .01   ! set upper lvls to dry
                                                        ! counters model (ruc) 
                                                        ! moist bias

                if(mode_fg .eq. 1)then
                    cf_modelfg(i,j,k) = rh_to_cldcv(rh_qc)    ! fractional_rh
                elseif(mode_fg .eq. 2)then
                    cf_modelfg(i,j,k) = 
     1              max(cf_modelfg(i,j,k),rh_to_cldcv(rh_qc)) ! fractional_rh
                elseif(mode_fg .eq. 3)then
                    continue
                elseif(mode_fg .eq. 4)then
                    cf_modelfg(i,j,k) = default_clear_cover
                else
                    write(6,*)' error: mode_fg not defined'
                    stop
                endif

              endif

            endif

 1000     enddo ! i
        enddo ! j
        enddo ! k (cloud height array level)

        write(6,*)' # rh values qced  high/low/ok'
     1            ,i_hum_high,i_hum_low,i_hum_ok
        write(6,*)' # cloud height grids above pressure grid'
     1            ,i_grid_high
        write(6,*)' # points set to cloud based on condensate = '
     1           ,i_condensate
               
        if(mode_fg .eq. 2)then
            write(6,*)
     1          ' used maximum of first guess lwc and rh implied cloud'
        endif

        write(6,*)' returning from get_modelfg with istatus = 1'
        istatus = 1

        return
        end
