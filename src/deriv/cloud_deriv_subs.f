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

        subroutine insert_thin_lwc_ice(clouds_3d,clouds_3d_pres
     1        ,heights_3d
     1        ,temp_3d,cldalb_in,cld_hts,ni,nj,nk,kcld,idb,jdb
     1        ,thresh_thin_lwc_ice
     1        ,pres_3d,slwc,cice,istatus)

        use mem_namelist, only: r_missing_data
        use cloud_rad

        include 'laps_cloud.inc'

!       insert lwc and ice in areas of thin cloud layers
        real clouds_3d(ni,nj,kcld)       ! input cloud fraction   
        real clouds_3d_pres(ni,nj,nk)    ! input
        real cldalb_in(ni,nj)            ! input
        real heights_3d(ni,nj,nk)        ! input
        real temp_3d(ni,nj,nk)           ! input
        real slwc(ni,nj,nk)              ! input/output (g/m**3)
        real cice(ni,nj,nk)              ! input/output (g/m**3)
        real pres_3d(ni,nj,nk)           ! input
        real pressures_pa(nk)            ! local
        real k_to_c

        integer max_layers
        parameter (max_layers=100)

        real a(max_layers)       ! cloud fractions of layers
        integer ik(max_layers)   ! height level representative of cloud layers
        integer ktop(max_layers) ! height level representative of cloud layers
        integer kbot(max_layers) ! height level representative of cloud layers
        integer lyr_indx(kcloud) ! layer index for each cloud analysis level

        write(6,*)' subroutine insert_thin_lwc_ice'

        zero = 0.

        iwrite = 0

        write(6,*)
     1  '  i   j   k thick    tau     cvr    t_c     slwc     cice' 
     1  ,'   ilyr nlyr'

!       convert from cloud cover to discreet cloud layer indices (cvr to a)
        do i = 1,ni
        do j = 1,nj
          if(i .eq. idb .and. j .eq. jdb)then
            idebug = 1
          else
            idebug = 0
          endif
           
          do k = 1,nk
             pressures_pa(k) = pres_3d(i,j,k)
          enddo ! k

          nlyr = 0

          do k = kcld,1,-1
             if(k .eq. kcld)then
                 clouds_3d_u = 0.
             else
                 clouds_3d_u = clouds_3d(i,j,k+1)
             endif

             if(clouds_3d(i,j,k)   .ge. thresh_thin_lwc_ice .and.
     1          clouds_3d_u        .lt. thresh_thin_lwc_ice)then ! top of layer
                nlyr = nlyr + 1
                a(nlyr) = clouds_3d(i,j,k)
                ik(nlyr) = k
                ktop(nlyr) = k

             else
                if(nlyr .ge. 1)then
                    if(clouds_3d(i,j,k) .gt. a(nlyr))then ! max within layer
                        a(nlyr) = clouds_3d(i,j,k)
                        ik(nlyr) = k
                    endif
                endif
             endif

             if(clouds_3d(i,j,k) .ge. thresh_thin_lwc_ice)then   ! still within layer
                 lyr_indx(k) = nlyr
                 kbot(nlyr) = k
             else                               ! below layer
                 lyr_indx(k) = 0
             endif

          enddo ! k

!         set up layer thickness, mean svp arrays

!         find thickness of layer
          do ilyr = 1,nlyr
              if(kbot(ilyr) .eq. 1)then
                  cloud_base = cld_hts(1)
              else
                  cloud_base = 0.5 * ( cld_hts(kbot(ilyr)  )
     1                              +  cld_hts(kbot(ilyr)-1)  )
              endif

              if(ktop(ilyr) .eq. kcloud)then
                  cloud_top = cld_hts(kcloud)
              else
                  cloud_top =  0.5 * ( cld_hts(ktop(ilyr)  )
     1                              +  cld_hts(ktop(ilyr)+1)  )
              endif

              thickness = max(cloud_top - cloud_base,1.)

!             compute the condensate concentration within the layer, the maximum
!             cloud cover for the layer is assummed to represent average of
!             the whole layer.

!             improvements are being considered based on this reference:
!             http://curry.eas.gatech.edu/map/pdf/mcfheyms97.pdf

              if(a(ilyr) .lt. 0.999)then
                  tau = -log(1. - a(ilyr))                        ! optical depth
                  tau = min(tau,30.)      
              else
                  tau = 30.
              endif

              scattering_coeff = tau / thickness                  ! meters**-1

!             convert from implied albedo to column optical depth
!             should maximum albedo of whole cloud column be used instead of
!             just this layer?
!             albarg = min(a(ilyr),0.930) ! bounded cloud fraction

              if(cldalb_in(i,j) .ne. r_missing_data)then          ! use vis sat
                  albarg = min(cldalb_in(i,j),0.930)
              else
                  albarg = min(maxval(a(1:nlyr)),.930)
              endif

              call albedo_to_clouds2(albarg
     1                             ,cloud_trans_l,cloud_trans_i
     1                             ,cloud_od_l,cloud_od_i
     1                             ,cloud_opac_l,cloud_opac_i)

              tau_l = cloud_od_l / float(nlyr)
              tau_i = cloud_od_i / float(nlyr)

              if(idebug .eq. 1)then
                  write(6,11)ilyr,a(ilyr),albarg,cloud_od_l,cloud_od_i
11                format(' ctr ilyr,a,albarg/cloud_od_l/cloud_od_i'
     1                             ,i3,4f9.4)
              endif

              scattering_coeff_l = tau_l / thickness              ! meters**-1
              scattering_coeff_i = tau_i / thickness              ! meters**-1

!             note tau = 3./2. * lwp / reff
!                  lwp = (tau / 1.5) * reff
!                  density = (alpha / 1.5) * reff
!             is scattering efficiency (q) being accounted for and is density 
!             being overestimated by a factor of ~2? consider adding new
!             q parameters from module 'cloud_rad'.

!             ice clouds
              rmvd = .000060                                      ! meters     (60 microns)
              rma =         3.14 * (rmvd/2.) ** 2                 ! meters**2  (pi r**2)
              rmv = 4./3. * 3.14 * (rmvd/2.) ** 3                 ! meters**3  (4/3 pi r**3)
   
              rnumber_density = scattering_coeff / rma            ! meters**-3
              h20_density = 1e6                                   ! g/m**3
              density_ave_i = rnumber_density * rmv * h20_density ! g/m**3
              density_ave_i = (scattering_coeff_i/1.5) * reff_cice
     1                                                 * h20_density

!             liquid clouds
              rmvd = .000020                                      ! meters     (20 microns)
              rma =         3.14 * (rmvd/2.) ** 2                 ! meters**2  (pi r**2)
              rmv = 4./3. * 3.14 * (rmvd/2.) ** 3                 ! meters**3  (4/3 pi r**3)
   
              rnumber_density = scattering_coeff / rma            ! meters**-3
              h20_density = 1e6                                   ! g/m**3
              density_ave_l = rnumber_density * rmv * h20_density ! g/m**3
              density_ave_l = (scattering_coeff_l/1.5) * reff_clwc
     1                                                 * h20_density

!             determine laps pressure levels of cloud base and top
!             qc cases where cloud layer on the height grid is near or above 
!                                                the top laps pressure level

              if(cloud_base .gt. heights_3d(i,j,nk) - 2.)then ! skip layer
                  write(6,*)' warning: skip layer '
     1                     ,i,j,cloud_base,cloud_top
                  go to 100            
              endif

              k_1d_base = int(height_to_zcoord3(cloud_base,heights_3d
     1                    ,pressures_pa,ni,nj,nk,i,j,istatus_base)) + 1   

              if(cloud_top  .gt. heights_3d(i,j,nk))then
                  cloud_top_before = cloud_top
                  cloud_top  =   heights_3d(i,j,nk)
                  write(6,*)' warning: clip layer ',i,j,cloud_base
     1                     ,cloud_top_before,cloud_top      
                  k_1d_top = nk ! avoids roundoff error at very top of domain
                  istatus_top = 1

              else
                  k_1d_top  = int(height_to_zcoord3(cloud_top 
     1                    ,heights_3d,pressures_pa,ni,nj,nk
     1                    ,i,j,istatus_top))

              endif

!             we have now defined a cloud base and top
              if(       k_1d_base .lt. 1 .or. k_1d_base .gt. nk  
     1             .or. k_1d_top  .lt. 1 .or. k_1d_top  .gt. nk
     1             .or. istatus_base .ne. 1 .or. istatus_top .ne. 1)then       
                  write(6,*)' error: bad return from height_to_zcoord3'       
                  write(6,*)k_1d_base,istatus_base,cloud_base
     1                     ,heights_3d(i,j,1)
                  write(6,*)k_1d_top,istatus_top,cloud_top
     1                     ,heights_3d(i,j,nk)     
                  istatus = 0
                  return
              endif

!             insert density value into lwc/ice field
              do k = k_1d_base,k_1d_top

!                 test if no condensate is currently at the grid point from the
!                 modified s-f method that was only applied to thick clouds
                  if((slwc(i,j,k) .eq. zero .and. 
     1                cice(i,j,k) .eq. zero       )  
     1                          .or. ! use density value from cloud
     1                               ! optical depth when available, instead of s-f
     1               cldalb_in(i,j) .ne. r_missing_data)then

!                   density is normalized according to cloud fraction within layer
                    density_l = density_ave_l * clouds_3d_pres(i,j,k) 
     1                                      / a(ilyr)

                    density_i = density_ave_i * clouds_3d_pres(i,j,k) 
     1                                      / a(ilyr)

                    temp1_c = -10.
                    temp2_c = -30.

                    temp1_k = c_to_k(temp1_c)
                    temp2_k = c_to_k(temp2_c)

                    if(temp_3d(i,j,k) .le. temp2_k)then          ! ice
                        cice(i,j,k) = density_i
                        tau_eff = tau_i 
                    elseif(temp_3d(i,j,k) .ge. temp1_k)then      ! lwc
                        slwc(i,j,k) = density_l
                        tau_eff = tau_l 
                    else                                        ! mixed
                        frac = (temp_3d(i,j,k) - temp2_k) 
     1                               / (temp1_k - temp2_k)
                        slwc(i,j,k) = density_l * frac
                        cice(i,j,k) = density_i * (1. - frac)
                        tau_eff = tau_l * frac + tau_i * (1. - frac)
                    endif


                    if(iwrite .lt. 30. .or. idebug .eq. 1)then
                        iwrite = iwrite + 1
                        write(6,1)i,j,k,thickness,tau_eff               
     1                           ,clouds_3d_pres(i,j,k)
     1                           ,k_to_c(temp_3d(i,j,k))
     1                           ,slwc(i,j,k),cice(i,j,k)
     1                           ,ilyr,nlyr
1                       format(3i4,f7.0,2f7.3,f7.1,2f9.4,2i5)
                    endif

                  else ! write out some thick cloud values
                    if(iwrite .lt. 30. .or. idebug .eq. 1)then
                        iwrite = iwrite + 1
                        write(6,2)i,j,k,thickness,tau             
     1                           ,clouds_3d_pres(i,j,k)
     1                           ,k_to_c(temp_3d(i,j,k))
     1                           ,slwc(i,j,k),cice(i,j,k)
     1                           ,ilyr,nlyr
2                       format(3i4,f7.0,f8.3,f7.3,f7.1,2f9.4,2i5,' k')
                    endif

                  endif

              enddo ! k

 100      enddo ! ilyr

        enddo ! j
        enddo ! i

        istatus = 1

        return
        end

        subroutine sao_drizzle_correction(pcp_type_2d,ni,nj
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k
     1              ,cloud_ceiling,r_missing_data)

!       adds drizzle to sfc precip type using sao wx string data

        real cloud_ceiling(ni,nj)
        real pcp_type_2d(ni,nj)
        real ri_s(maxstns), rj_s(maxstns)
        real lat(ni,nj),lon(ni,nj)
        real t_sfc_k(ni,nj)

!       declarations for lso file stuff
        real lat_s(maxstns), lon_s(maxstns)
        character wx_s(maxstns)*8, obstype(maxstns)*8

        logical l_valid_sta(maxstns)

        write(6,*)' adding drizzle information from saos'

        n_valid_sta = 0

!       set up i's and j's of the saos
        do ista = 1,n_obs_pos_b
            call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista)
     1                   ,lat,lon,ni,nj,ri_s(ista),rj_s(ista),istatus)

            l_valid_sta(ista) = .false.

!           does this station report precip?
            if(wx_s(ista)(1:7) .ne. 'unknown')then
                l_valid_sta(ista) = .true.
                n_valid_sta = n_valid_sta + 1
            endif
        enddo

        write(6,*)' nobs,n_valid = ',n_obs_pos_b,n_valid_sta

        n_low_ceil = 0
        n_no_pcp = 0
        n_l = 0
        n_zl = 0

        do i = 1,ni
        do j = 1,nj
            if(cloud_ceiling(i,j) .ne. r_missing_data
     1          .and. cloud_ceiling(i,j) .lt. 200.)then
                n_low_ceil = n_low_ceil + 1
                if(pcp_type_2d(i,j) .eq. 0.)then
                    n_no_pcp = n_no_pcp + 1

!                   determine nearest station to the grid point
!                   can be inside or outside the domain
                    distsq_min = 999999.
                    ista_nearest = 0

                    do ista = 1,n_obs_pos_b

!                       does this station report precip?
                        if(l_valid_sta(ista))then
!                           calculate distance of stations (grid points **2)
                            distsq = (float(i) - ri_s(ista))**2
     1                             + (float(j) - rj_s(ista))**2
                            if(distsq .lt. distsq_min)then ! nearest station
                                distsq_min = distsq
                                ista_nearest = ista
                            endif
                        endif

                    enddo ! ista

                    if(ista_nearest .ne. 0)then ! we found a reporting station
!                       parse nearest station
                        i_l = 0

                        call parse_wx_pcp(wx_s(ista_nearest),'l'
     1                                              ,ipresent,istatus)       
                        if(ipresent .eq. 1)i_l = 1

                        call parse_wx_pcp(wx_s(ista_nearest),'f'
     1                                              ,ipresent,istatus)       
                        if(ipresent .eq. 1)i_l = 1

                        if(i_l .eq. 1)then ! l or zl present
                            if(t_sfc_k(i,j) .gt. 273.15)then
                                n_l = n_l + 1
                                pcp_type_2d(i,j) = 6.0 ! drizzle
                            else
                                n_zl = n_zl + 1
                                pcp_type_2d(i,j) = 7.0 ! freezing drizzle
                            endif
                        endif
                    endif

                endif
            endif
        enddo
        enddo

        write(6,*)' # of grid points with ceiling < 200m = ',n_low_ceil
        write(6,*)' # of grid points also with no precip = ',n_no_pcp
        write(6,*)' # of grid points with drizzle =        ',n_l
        write(6,*)' # of grid points with frzing drizzle = ',n_zl

        return
        end


        subroutine sao_snow_correction(pcp_type_2d,ni,nj
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k,twet_snow
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)

!       adds snow to sfc precip type using sao wx string data
!       1995 nov 29 s. albers - improve use of sao's to add snow in ptt field.
!                               cloud ceiling threshold replaced with
!                               thresholds on cloud cover and sfc dewpoint
!                               depression.

        real cvr_max(ni,nj)
        real pcp_type_2d(ni,nj)
        real t_sfc_k(ni,nj)
        real td_sfc_k(ni,nj)
        real dbz_low_2d(ni,nj)
        real ri_s(maxstns), rj_s(maxstns)
        real lat(ni,nj),lon(ni,nj)

        integer ipresent_s(maxstns)

!       declarations for lso file stuff
        real lat_s(maxstns), lon_s(maxstns)
        character wx_s(maxstns)*8, obstype(maxstns)*8

        logical l_valid_sta(maxstns)

        write(6,*)' adding snow information from saos'

        n_valid_sta = 0
        ipresent_s = 0

!       set up i's and j's of the saos
        do ista = 1,n_obs_pos_b
            call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista)
     1                   ,lat,lon,ni,nj,ri_s(ista),rj_s(ista),istatus)

            l_valid_sta(ista) = .false.

            i_sta = nint(ri_s(ista))
            j_sta = nint(rj_s(ista))

            if(    i_sta .ge. 1 .and. i_sta .le. ni
     1       .and. j_sta .ge. 1 .and. j_sta .le. nj)then

!               does this station report precip and snow?
                if(wx_s(ista)(1:7) .ne. 'unknown')then
                    l_valid_sta(ista) = .true.
                    n_valid_sta = n_valid_sta + 1
                    call parse_wx_pcp(wx_s(ista),'s',ipresent_s(ista)
     1                               ,istatus)    

                endif
            endif
        enddo

        write(6,*)' nobs,n_valid = ',n_obs_pos_b,n_valid_sta

        n_cvr= 0
        n_moist = 0
        n_no_pcp = 0
        n_cold = 0
        n_s = 0

        do i = 1,ni
        do j = 1,nj
            if(cvr_max(i,j) .gt. 0.5)then ! ample cloud cover
              n_cvr = n_cvr + 1

              if( (t_sfc_k(i,j) - td_sfc_k(i,j)) .lt. 10.)then ! moist air
                n_moist = n_moist + 1

                if(pcp_type_2d(i,j) .eq. 0.)then ! no snow diagnosed at grid pt
                  n_no_pcp = n_no_pcp + 1

!                 approx wet bulb temp
                  tw_c = 0.5 * (t_sfc_k(i,j)
     1                       + td_sfc_k(i,j)) - 273.15
                  if(tw_c .le. twet_snow)then ! cold enough for snow
                    n_cold = n_cold + 1

!                   determine nearest station to the grid point
!                   must be inside the domain
                    distsq_min = 999999.
                    ista_nearest = 0
                    dbz_at_sta = 999.

                    do ista = 1,n_obs_pos_b

!                       valid station is in domain and reports precip
                        if(l_valid_sta(ista))then

!                           calculate distance of stations (grid points **2)
                            distsq = (float(i) - ri_s(ista))**2
     1                             + (float(j) - rj_s(ista))**2
                            if(distsq .lt. distsq_min)then ! nearest station
                                distsq_min = distsq
                                ista_nearest = ista
                                i_sta = nint(ri_s(ista))
                                j_sta = nint(rj_s(ista))
                                dbz_at_sta = dbz_low_2d(i_sta,j_sta)
                            endif
                        endif ! valid station

                    enddo ! ista

                    if(ista_nearest .ne. 0)then ! we found a reporting station
!                       parse nearest station
                        i_s = 0

                        if(ipresent_s(ista_nearest) .eq. 1)i_s = 1

                        if(i_s .eq. 1)then ! s present
                            if(dbz_at_sta .le. 0.)then ! station has no echo
                                n_s = n_s + 1
                                pcp_type_2d(i,j) = 2.0 ! snow
                            endif
                        endif ! s present
                    endif ! we found a reporting station

                  endif ! cold enough for snow
                endif ! no snow diagnosed
              endif ! moist air
            endif ! ample cloud cover
        enddo
        enddo

        write(6,*)' # of grid points with cvr > 0.5 = ',n_cvr
        write(6,*)' # of grid points also with t-td < 10 = ',n_moist
        write(6,*)' # of grid points also with no precip = ',n_no_pcp
        write(6,*)' # of grid points also cold enough = ',n_cold
        write(6,*)' # of grid points with snow added =     ',n_s

        return
        end


        subroutine sao_rain_correction(pcp_type_2d,ni,nj
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k,twet_snow
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)

!       adds rain to sfc precip type using sao wx string data
!       1996 oct 31 s. albers - improve use of sao's to add rain/freezing 
!                               rain in ptt field. cloud ceiling threshold 
!                               replaced with thresholds on cloud cover 

        real cvr_max(ni,nj)
        real pcp_type_2d(ni,nj)
        real t_sfc_k(ni,nj)
        real td_sfc_k(ni,nj)
        real dbz_low_2d(ni,nj)
        real ri_s(maxstns), rj_s(maxstns)
        real lat(ni,nj),lon(ni,nj)

!       declarations for lso file stuff
        real lat_s(maxstns), lon_s(maxstns)
        character wx_s(maxstns)*8, obstype(maxstns)*8

        logical l_valid_sta(maxstns)

        write(6,*)' adding rain information from stations'

        n_valid_sta = 0

!       set up i's and j's of the saos
        do ista = 1,n_obs_pos_b
            call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista)
     1                   ,lat,lon,ni,nj,ri_s(ista),rj_s(ista),istatus)

            l_valid_sta(ista) = .false.

            i_sta = nint(ri_s(ista))
            j_sta = nint(rj_s(ista))

            if(    i_sta .ge. 1 .and. i_sta .le. ni
     1       .and. j_sta .ge. 1 .and. j_sta .le. nj)then

!               does this station report precip?
                if(wx_s(ista)(1:7) .ne. 'unknown')then
                    l_valid_sta(ista) = .true.
                    n_valid_sta = n_valid_sta + 1
                endif
            endif
        enddo

        write(6,*)' nobs,n_valid = ',n_obs_pos_b,n_valid_sta

        n_cvr= 0
        n_moist = 0
        n_no_pcp = 0
        n_r = 0
        n_zr = 0
        n_r_sta = 0
        n_zr_sta = 0
        n_sta_noecho = 0

        do i = 1,ni
        do j = 1,nj
            if(cvr_max(i,j) .gt. 0.5)then ! ample cloud cover
                n_cvr = n_cvr + 1
                if(pcp_type_2d(i,j) .eq. 0.)then ! no rain diagnosed at grid pt
                    n_no_pcp = n_no_pcp + 1

!                   determine nearest station to the grid point
!                   must be inside the domain
                    distsq_min = 999999.
                    ista_nearest = 0
                    dbz_at_sta = 999.

                    do ista = 1,n_obs_pos_b
!                       valid station is in domain and reports precip
                        if(l_valid_sta(ista))then

!                           calculate distance of stations (grid points **2)
                            distsq = (float(i) - ri_s(ista))**2
     1                             + (float(j) - rj_s(ista))**2
                            if(distsq .lt. distsq_min)then ! nearest station
                                distsq_min = distsq
                                ista_nearest = ista
                                i_sta = nint(ri_s(ista))
                                j_sta = nint(rj_s(ista))
                                dbz_at_sta = dbz_low_2d(i_sta,j_sta)
                            endif
                      endif ! valid station
                    enddo ! ista

                    if(ista_nearest .ne. 0)then ! we found a reporting station
!                       parse nearest station
                        call parse_wx_pcp(wx_s(ista_nearest),'r'
     1                                             ,ipresent_r,istatus)       

                        call parse_wx_pcp(wx_s(ista_nearest),'z'
     1                                             ,ipresent_zr,istatus)       

                        if(ipresent_r .eq. 1)then ! r/zr present
                            n_r_sta = n_r_sta + 1

                            if(dbz_at_sta .le. 0.)then ! station has no echo
                                n_sta_noecho = n_sta_noecho + 1

!                               approx wet bulb temp
                                tw_c = 0.5 * (t_sfc_k(i,j)
     1                                     + td_sfc_k(i,j)) - 273.15

                                if(ipresent_zr .eq. 0)then     ! stn has rain
                                    if(tw_c .gt. twet_snow)then
                                        n_r = n_r + 1
                                        pcp_type_2d(i,j) = 1.0 ! rain
                                    else ! wet bulb below zero (indeterminate)
                                    endif
                                elseif(ipresent_zr .eq. 1)then ! stn has zr    
                                    n_zr_sta = n_zr_sta + 1
                                    if(t_sfc_k(i,j) .gt. 273.15)then
                                        n_r = n_r + 1
                                        pcp_type_2d(i,j) = 1.0 ! rain
                                    else ! dry bulb below zero
                                        n_zr = n_zr + 1
                                        pcp_type_2d(i,j) = 3.0 ! frz rain
                                    endif
                                endif
                            endif
                        endif ! r present
                    endif ! we found a reporting station

                endif ! no rain diagnosed
            endif ! ample cloud cover
        enddo
        enddo

        write(6,*)
     1      ' # of grid points with cvr > 0.5 =           ',n_cvr
        write(6,*)
     1      ' # of grid points also with no precip =      ',n_no_pcp
        write(6,*)
     1      ' # of grid points also with rain/zr at sta = ',n_r_sta
        write(6,*)
     1      ' # of grid points also with no echo at sta = ',n_sta_noecho      
        write(6,*)
     1      ' # of grid points with rain added =          ',n_r
        write(6,*)
     1      ' # of grid points with zr at sta =           ',n_zr_sta
        write(6,*)
     1      ' # of grid points with zr added =            ',n_zr

        return
        end



        subroutine sao_precip_correction(pcp_type_2d,ni,nj
     1              ,n_obs_pos_b,obstype,wx_s,lat_s,lon_s,maxstns
     1              ,ri_s,rj_s,lat,lon
     1              ,t_sfc_k,td_sfc_k
     1              ,dbz_low_2d
     1              ,cvr_max,r_missing_data)

!       increase radar reflectivity threshold for precip if the radar has
!       echo over the nearest station, but that station says no precip.
!       this routine sets the precip type at such grid points to "no precip".
!       the highest allowable reflectivity threshold for precip is 10 dbz.

        real cvr_max(ni,nj)
        real pcp_type_2d(ni,nj)
        real t_sfc_k(ni,nj)
        real td_sfc_k(ni,nj)
        real dbz_low_2d(ni,nj)
        real ri_s(maxstns), rj_s(maxstns)
        real lat(ni,nj),lon(ni,nj)

!       declarations for lso file stuff
        real lat_s(maxstns), lon_s(maxstns)
        character wx_s(maxstns)*8, obstype(maxstns)*8

!       this routine blanks the precip type if the neighboring station has no
!       precip and the reflectivity at the grid point is less than the
!       reflectivity at the station.

        write(6,*)' adjusting radar reflectivity threshold around stns'

!       set up i's and j's of the stns
        do ista = 1,n_obs_pos_b
            call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista)
     1                   ,lat,lon,ni,nj,ri_s(ista),rj_s(ista),istatus)
        enddo

        n_keep_precip = 0
        n_subtract = 0
        n_pcp = 0
        n_no_pcp_sta = 0
        n_no_pcp_sta_echo = 0

        do i = 1,ni
        do j = 1,nj
                if(pcp_type_2d(i,j) .gt. 0.)then ! precip diagnosed at grid pt
                    n_pcp = n_pcp + 1

!                   determine nearest station to the grid point
!                   must be inside the domain
                    distsq_min = 999999.
                    ista_nearest = 0
                    dbz_at_sta = 999.

                    do ista = 1,n_obs_pos_b
                      i_sta = nint(ri_s(ista))
                      j_sta = nint(rj_s(ista))

                      if(    i_sta .ge. 1 .and. i_sta .le. ni
     1                 .and. j_sta .ge. 1 .and. j_sta .le. nj)then

!                       does this station report precip?
                        if(      
!    1                           obstype(ista)(1:4) .ne. 'meso'
!    1                     .and. obstype(ista)(1:4) .ne. 'cdot'
!    1                     .and. obstype(ista)(7:8) .ne. '1a'
     1                           wx_s(ista)(1:7)    .ne. 'unknown'
     1                                                           )then
!                           calculate distance of stations (grid points **2)
                            distsq = (float(i) - ri_s(ista))**2
     1                             + (float(j) - rj_s(ista))**2
                            if(distsq .lt. distsq_min)then ! nearest station
                                distsq_min = distsq
                                ista_nearest = ista
                                dbz_at_sta = dbz_low_2d(i_sta,j_sta)
                            endif
                        endif
                      endif ! station in domain
                    enddo ! ista

                    if(ista_nearest .ne. 0)then ! we found a reporting station

!                       parse nearest station
                        i_precip_sta = 0

!                       snow
                        call parse_wx_pcp(wx_s(ista_nearest),'s'
     1                                              ,ipresent,istatus)       
                        if(ipresent .eq. 1)i_precip_sta = 1

!                       rain
                        call parse_wx_pcp(wx_s(ista_nearest),'r'
     1                                              ,ipresent,istatus)       
                        if(ipresent .eq. 1)i_precip_sta = 1

!                       drizzle
                        call parse_wx_pcp(wx_s(ista_nearest),'l'
     1                                              ,ipresent,istatus)       
                        if(ipresent .eq. 1)i_precip_sta = 1

!                       frz drizzle
                        call parse_wx_pcp(wx_s(ista_nearest),'f'
     1                                              ,ipresent,istatus)       
                        if(ipresent .eq. 1)i_precip_sta = 1

!                       ice pellets
                        call parse_wx_pcp(wx_s(ista_nearest),'i'
     1                                              ,ipresent,istatus)       
                        if(ipresent .eq. 1)i_precip_sta = 1

                        if(i_precip_sta .eq. 0)then ! precip absent at station

                            n_no_pcp_sta = n_no_pcp_sta + 1

                            if(dbz_at_sta .ge. 0.)then ! station has echo
                                n_no_pcp_sta_echo =
     1                          n_no_pcp_sta_echo + 1
                                if(dbz_low_2d(i,j) .le. dbz_at_sta .and.
     1                             dbz_low_2d(i,j) .le. 10.        )then
                                    n_subtract = n_subtract + 1
                                    pcp_type_2d(i,j) = 0.0 ! no precip
                                else
                                    n_keep_precip = n_keep_precip + 1
                                endif
                            endif
                        endif ! precip present
                    endif ! we found a reporting station

                endif ! precip diagnosed
        enddo
        enddo

        write(6,*)' # of grid points with precip=                  '
     1                                            ,n_pcp
        write(6,*)' # of grid points also with no precip at stn=   '
     1                                            ,n_no_pcp_sta
        write(6,*)' # of grid points also with radar echo at stn=  '
     1                                            ,n_no_pcp_sta_echo
        write(6,*)' # of grid points with precip subtracted=       '
     1                                            ,n_subtract
        write(6,*)' # of grid points with precip kept=             '
     1                                            ,n_keep_precip

        return
        end

        subroutine parse_wx_pcp(c8_wx,c1_pcp,ipresent,istatus)

        character*8 c8_wx
        character*2 c2_pcp
        character*1 c1_pcp
 
        ipresent = 0

!       snow
        if(c1_pcp .eq. 's')then
            call parse_wx_string(c8_wx,'sn',iparse,istatus)           
            if(iparse .eq. 1)ipresent = 1

            call parse_wx_string(c8_wx,'sg',iparse,istatus)           
            if(iparse .eq. 1)ipresent = 1

            call parse_wx_string(c8_wx,'ic',iparse,istatus)           
            if(iparse .eq. 1)ipresent = 1
        endif

!       rain
        if(c1_pcp .eq. 'r')then
            call parse_wx_string(c8_wx,'ra',iparse,istatus)           
            if(iparse .eq. 1)ipresent = 1
        endif

!       freezing rain
        if(c1_pcp .eq. 'z')then
            call parse_wx_string(c8_wx,'fz',iparse_fz,istatus)           
            call parse_wx_string(c8_wx,'ra',iparse_ra,istatus)           
            if(iparse_fz .eq. 1 .and. iparse_ra .eq. 1)ipresent = 1
        endif

!       drizzle
        if(c1_pcp .eq. 'l')then
            call parse_wx_string(c8_wx,'dz',iparse,istatus)           
            if(iparse .eq. 1)ipresent = 1
        endif

!       freezing drizzle
        if(c1_pcp .eq. 'f')then
            call parse_wx_string(c8_wx,'fz',iparse_fz,istatus)           
            call parse_wx_string(c8_wx,'dz',iparse_dz,istatus)           
            if(iparse_fz .eq. 1 .and. iparse_dz .eq. 1)ipresent = 1
        endif

!       ice pellets
        if(c1_pcp .eq. 'i')then
            call parse_wx_string(c8_wx,'pe',iparse,istatus)           
            if(iparse .eq. 1)ipresent = 1
            call parse_wx_string(c8_wx,'pl',iparse,istatus)           
            if(iparse .eq. 1)ipresent = 1
        endif

        istatus = 1

        return
        end



        subroutine parse_wx_string(c8_wx,c2_pcp,iparse,istatus)

        character*8 c8_wx
        character*2 c2_pcp
 
        iparse = 0

        do i = 1,7
            if(c8_wx(i:i+1) .eq. c2_pcp)then
                iparse = 1
            endif
        enddo ! i

        istatus = 1

        return
        end


        subroutine put_laps_3d_multi(i4time,ext,var_3d
     1                                     ,units_3d,comment_3d
     1                                     ,array1,array2
     1                                     ,array3,array4
     1                                     ,array5,array6
     1                                     ,nx_l_1,ny_l_1,nz_l_1
     1                                     ,nx_l_2,ny_l_2,nz_l_2
     1                                     ,nx_l_3,ny_l_3,nz_l_3
     1                                     ,nx_l_4,ny_l_4,nz_l_4
     1                                     ,nx_l_5,ny_l_5,nz_l_5
     1                                     ,nx_l_6,ny_l_6,nz_l_6
     1                                     ,n_3d_fields,istatus)       


        real array1(nx_l_1,ny_l_1,nz_l_1)
        real array2(nx_l_2,ny_l_2,nz_l_2)
        real array3(nx_l_3,ny_l_3,nz_l_3)
        real array4(nx_l_4,ny_l_4,nz_l_4)
        real array5(nx_l_5,ny_l_5,nz_l_5)
        real array6(nx_l_6,ny_l_6,nz_l_6)
        character*125 comment_3d(n_3d_fields)
        character*10 units_3d(n_3d_fields)
        character*3 var_3d(n_3d_fields)
        character*(*) ext

        write(6,*)' subroutine put_laps_3d_multi...'

        l = 1
        call put_laps_multi_3d(i4time,ext,var_3d(l),units_3d(l),
     1     comment_3d(l),array1,nx_l_1,ny_l_1,nz_l_1,1,istatus)       
        if(istatus .ne. 1)return
        if(l .eq. n_3d_fields)return

        l = 2
        call put_laps_multi_3d_append(i4time,ext,var_3d(l),units_3d(l),       
     1     comment_3d(l),array2,nx_l_2,ny_l_2,nz_l_2,1,istatus)       
        if(istatus .ne. 1)return
        if(l .eq. n_3d_fields)return

        l = 3
        call put_laps_multi_3d_append(i4time,ext,var_3d(l),units_3d(l),       
     1     comment_3d(l),array3,nx_l_3,ny_l_3,nz_l_3,1,istatus)       
        if(istatus .ne. 1)return
        if(l .eq. n_3d_fields)return

        l = 4
        call put_laps_multi_3d_append(i4time,ext,var_3d(l),units_3d(l),       
     1     comment_3d(l),array4,nx_l_4,ny_l_4,nz_l_4,1,istatus)       
        if(istatus .ne. 1)return
        if(l .eq. n_3d_fields)return

        l = 5
        call put_laps_multi_3d_append(i4time,ext,var_3d(l),units_3d(l),       
     1     comment_3d(l),array5,nx_l_5,ny_l_5,nz_l_5,1,istatus)       
        if(istatus .ne. 1)return
        if(l .eq. n_3d_fields)return

        l = 6
        call put_laps_multi_3d_append(i4time,ext,var_3d(l),units_3d(l),       
     1     comment_3d(l),array6,nx_l_6,ny_l_6,nz_l_6,1,istatus)       
        if(istatus .ne. 1)return
        if(l .eq. n_3d_fields)return

        write(6,*)' error: n_3d_fields exceeds limit ',n_3d_fields

        return
        end


        subroutine put_laps_multi_3d_append(i4time,ext,var_2d,units_2d,
     1                          comment_2d,field_3d,ni,nj,nk,nf,istatus)

        logical ltest_vertical_grid

        character*150 directory
        character*(*) ext

        character*125 comment_3d(nk*nf),comment_2d(nf)
        character*10 units_3d(nk*nf),units_2d(nf)
        character*3 var_3d(nk*nf),var_2d(nf)
        integer lvl_3d(nk*nf)
        character*4 lvl_coord_3d(nk*nf)

        real field_3d(ni,nj,nk,nf)

        istatus = 0

        call get_directory(ext,directory,len_dir)

        do l = 1,nf
            write(6,11)directory,ext(1:5),var_2d(l)
11          format(' writing 3d ',a50,1x,a5,1x,a3)
        enddo ! l

        do l = 1,nf
          do k = 1,nk

            iscript_3d = (l-1) * nk + k

            units_3d(iscript_3d)   = units_2d(l)
            comment_3d(iscript_3d) = comment_2d(l)
            if(ltest_vertical_grid('height'))then
                lvl_3d(iscript_3d) = zcoord_of_level(k)/10
                lvl_coord_3d(iscript_3d) = 'msl'
            elseif(ltest_vertical_grid('pressure'))then
                lvl_3d(iscript_3d) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(iscript_3d) = 'hpa'
            else
                write(6,*)' error, vertical grid not supported,'
     1                   ,' this routine supports pressure or height'
                istatus = 0
                return
            endif

            var_3d(iscript_3d) = var_2d(l)

          enddo ! k
        enddo ! l

        call write_laps_multi(i4time,directory,ext,ni,nj,
     1  nk*nf,nk*nf,var_3d,lvl_3d,lvl_coord_3d,units_3d,
     1                     comment_3d,field_3d,istatus)

        if(istatus .ne. 1)return

        istatus = 1

        return
        end
