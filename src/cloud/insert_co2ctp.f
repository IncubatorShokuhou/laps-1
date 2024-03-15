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
c
c

        subroutine insert_co2ctp(i4time,cld_hts,heights_3d                ! i
     1            ,imax,jmax,kmax,kcloud,r_missing_data                   ! i
     1            ,l_use_co2,latency_co2                                  ! i 
     1            ,default_clear_cover                                    ! i
     1            ,lat,lon,ix_low,ix_high,iy_low,iy_high                  ! i
     1            ,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd       ! i/o
     1            ,istatus)                                               ! o

!       input
        real lat(imax,jmax),lon(imax,jmax)
        real cld_hts(kcloud)
        real heights_3d(imax,jmax,kmax)

!       arrays for cloud soundings
        real cld_snd(max_cld_snd,kcloud)
        real wt_snd(max_cld_snd,kcloud)
        integer i_snd(max_cld_snd)
        integer j_snd(max_cld_snd)

        logical lstat_co2_a(imax,jmax)
        logical l_use_co2

!       local
        real pct_pa(imax,jmax)
        real lca(imax,jmax)
        real cldtop_co2_pa_a(imax,jmax)
        real cloud_frac_co2_a(imax,jmax)

        character*31 ext
        character var*3,units*10
        character*125 comment

        write(6,*)' subroutine insert_c02ctp...'

        write(6,*)
     1    ' number of cumulative cloud soundings before co2-slice = '
     1                                                     ,n_cld_snd

!       obtain nesdis cloud-top pressure
        if(l_use_co2)then
            i4_co2_window = latency_co2

            write(6,*)' getting nesdis cloud-top pressure'
            ext = 'ctp'
            var = 'pct'
            ilevel = 0
            call get_laps_2dgrid(i4time,i4_co2_window,i4time_nearest
     1                      ,ext,var       
     1                      ,units,comment,imax,jmax,pct_pa,ilevel
     1                      ,istat_pct)       
            if(abs(istat_pct) .ne. 1)then
                write(6,*)' note: cannot read nesdis cloud-top pressure'       
            endif

!           obtain nesdis cloud-fraction
            write(6,*)' getting nesdis cloud-fraction'
            ext = 'ctp'
            var = 'lca'
            ilevel = 0
            call get_laps_2dgrid(i4time,i4_co2_window,i4time_nearest
     1                          ,ext,var
     1                          ,units,comment,imax,jmax,lca,ilevel
     1                          ,istat_lca)       
            if(abs(istat_lca) .ne. 1)then
                write(6,*)' note: cannot read nesdis cloud-fraction'
            endif

        endif ! l_use_co2

!       place co2-slicing cloud-top pressure in sparse array

        cloud_frac_co2_a = r_missing_data
        cldtop_co2_pa_a  = r_missing_data
        lstat_co2_a = .false.
        icount = 0

        if(abs(istat_pct) .eq. 1 .and. abs(istat_lca) .eq. 1 
     1                           .and. l_use_co2                )then
            write(6,*)' extracting co2-slicing info from nesdis data'
            do j = 1,jmax
            do i = 1,imax
!               test for partial cloudiness
                if(lca(i,j) .gt. 0. .and. lca(i,j) .lt. 1.0)then ! use co2 data
                    if(pct_pa(i,j) .lt. 100000.)then
                        icount = icount + 1
                        cloud_frac_co2_a(i,j) = lca(i,j)
                        cldtop_co2_pa_a(i,j) = pct_pa(i,j) 
                        lstat_co2_a(i,j) = .true.
                    endif
                endif
            enddo ! i
            enddo ! j

            istat_co2 = 0 ! for now
        else
            istat_co2 = 0
        endif

!       add to cloud sounding arrays

        n_cloud = 0 ! # of co2 slicing soundings

        do j = 1,jmax
        do i = 1,imax
            if(lstat_co2_a(i,j))then ! create this cloud sounding
                call pressure_to_height(cldtop_co2_pa_a(i,j),heights_3d       
     1                                 ,imax,jmax,kmax,i,j,ctop_m
     1                                 ,istatus)       

                cbase_m = ctop_m - cld_thk(ctop_m)
                cbuf_low = cbase_m - 500.

                cover_rpt = cloud_frac_co2_a(i,j)

                if(n_cloud .le. 10)then
                    write(6,*)' working with co2 cloud top:'
     1                       ,nint(cldtop_co2_pa_a(i,j))
     1                       ,nint(cbase_m),nint(ctop_m),cover_rpt
     1                       ,' at ',i,j
                endif

                n_cloud = n_cloud + 1

                n_cld_snd = n_cld_snd + 1 ! # of cumulative cloud soundings

                do k=1,kcloud
                    cover = cover_rpt

!                   fill in cloud layer
                    if(cld_hts(k) .ge. cbase_m .and.
     1                 cld_hts(k) .lt. ctop_m                  )then
                        call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                              ,n_cld_snd,max_cld_snd
     1                              ,kcloud,i,j,k,cover,1.)
                        if(n_cloud .le. 10)then
                            write(6,*)' fill in k = ',k,cover,cld_hts(k)
                        endif
                    endif

                    cover = default_clear_cover

!                   fill in clear buffer below cloud layer
                    if(cld_hts(k) .lt. cbase_m .and.
     1                 cld_hts(k) .ge. cbuf_low               )then      
                        call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                              ,n_cld_snd,max_cld_snd
     1                              ,kcloud,i,j,k,cover,1.)
                        if(n_cloud .le. 10)then
                            write(6,*)' fill in k = ',k,cover,cld_hts(k)       
                        endif
                    endif

!                   fill in clear sky above cloud layer
                    if(cld_hts(k) .ge. ctop_m                 )then      
                        call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                              ,n_cld_snd,max_cld_snd
     1                              ,kcloud,i,j,k,cover,1.)
                        if(n_cloud .le. 10)then
                            write(6,*)' fill in k = ',k,cover,cld_hts(k)       
                        endif
                    endif

                enddo ! k cld_hts

            endif ! co2 measurement at this grid point
        enddo ! i
        enddo ! j

        write(6,*)' number of valid co2-slicing soundings = ',n_cloud

        write(6,*)
     1    ' number of cumulative cloud soundings after co2-slice = '
     1                                                     ,n_cld_snd

        return
        end
