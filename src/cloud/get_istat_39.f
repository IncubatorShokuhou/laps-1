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

        subroutine get_istat_39(t39_k,tb8_k,solar_alt,r_missing_data    ! i
     1                         ,rlaps_land_frac,ni,nj                   ! i
     1                         ,static_albedo                           ! i
     1                         ,istat_39_a)                             ! o

!       this routine returns the status of stratus cloud detection for each 
!       grid point (containing liquid water).
!       -1 = cloud is not present (of any kind)
!        0 = indeterminate cloud existence
!       +1 = cloud is present (and is stratus)

!       we think this works even with thin cloud scenarios that could cause
!       an underestimation of cloud-top temperature when using 'tb8_k'.

        real    t39_k(ni,nj)
        real    tb8_k(ni,nj)
        real    solar_alt(ni,nj)
        real    rlaps_land_frac(ni,nj)
        real    static_albedo(ni,nj)          ! static albedo database
        integer istat_39_a(ni,nj)
        integer istat_39_buff(ni,nj)
        integer icount(-1:+1)

        real k_to_c

        write(6,*)' subroutine get_istat_39...'

        do ic = -1,+1
            icount(ic) = 0
        enddo ! ic

        do i = 1,ni
        do j = 1,nj
            if(tb8_k(i,j) .ne. r_missing_data .and.
     1         t39_k(i,j) .ne. r_missing_data      )then
                t39_c = k_to_c(t39_k(i,j))
                tb8_c = k_to_c(tb8_k(i,j))

                if(tb8_c          .gt. -10. ! cloud-top "definitely" warm
     1                                      ! enough for possibility of liquid
     1                                      ! water
     1                     .and.
     1             solar_alt(i,j) .lt. 0.        )then

                    t39_diff = t39_k(i,j) - tb8_k(i,j)

                    istat_39_a(i,j) = 0                       ! indeterminate

                    if(
     1                  (t39_diff .le. -2.5 .and. 
     1                                    rlaps_land_frac(i,j) .gt. 0.5)      
     1                                 .or.
     1                  (t39_diff .le. -2.0 .and. 
     1                                    rlaps_land_frac(i,j) .le. 0.5)      
     1                                             )then ! sufficient diff     
                        istat_39_a(i,j) = +1

                        if(icount(1) .le. 20)then             ! log output
                            write(6,11)i,j,t39_c,tb8_c,t39_c-tb8_c
 11                         format(' stratus present t39/tb8/df: '
     1                            ,2i5,3f8.1)
                        endif 
                    endif

                    if(t39_diff .gt. +0.0 .and. t39_diff .lt. +2.5
     1                          .and. tb8_c .gt. 0.    
     1                          .and. rlaps_land_frac(i,j) .gt. 0.5)then ! land
                        istat_39_a(i,j) = 0 ! -1              ! warm - no cloud?
                    endif

                    if(t39_diff .gt. +0.0 .and. t39_diff .lt. +2.5
     1                          .and. tb8_c .gt. 0.    
     1                          .and. rlaps_land_frac(i,j) .lt. 0.5)then ! sea
                        istat_39_a(i,j) = 0 ! -1              ! warm - no cloud?
                    endif

                else ! sun above horizon or possibly cold cloud top: ambiguous
                    istat_39_a(i,j) = 0 ! indeterminate

                endif

            else ! we have missing data
                istat_39_a(i,j) = 0     ! indeterminate

            endif

!           screen points according to static albedo
            if(static_albedo(i,j) .eq. r_missing_data .or.
     1         static_albedo(i,j) .gt. 0.18                  )then
                istat_39_a(i,j) = 0     ! indeterminate
            endif

            icount(istat_39_a(i,j)) = icount(istat_39_a(i,j)) + 1

        enddo 
        enddo

        write(6,12)icount(-1),icount(0),icount(+1)
 12     format('  3.9u cloud mask status vals [-1,0,+1] =',3i8)       

!       filter out isolated points
        istat_39_buff = istat_39_a

        do i = 2,ni-1
        do j = 2,nj-1
            if(istat_39_buff(i,j) .eq. 1)then
!               count neighbors
                n_neighbors = 0
                do ii = i-1,i+1
                do jj = j-1,j+1
                    if(istat_39_buff(ii,jj) .eq. 1)then
                        n_neighbors = n_neighbors + 1
                    endif
                enddo ! jj
                enddo ! ii

                if(n_neighbors .gt. 2)then
!                   keep if at least 2 neighbors (plus central pt) are there
                    istat_39_a(i,j) = istat_39_buff(i,j) ! 1
                else
!                   throw out if isolated point or just one neighbor
                    istat_39_a(i,j) = 0
                endif
            else
                istat_39_a(i,j) = istat_39_buff(i,j)
            endif
        enddo ! j
        enddo ! i

        return
        end 


        subroutine get_istat_39_lwc(t39_k,tb8_k,solar_alt
     1                             ,r_missing_data,ni,nj,istat_39_lwc_a)

!       this routine returns the status of lwc detection given that a cloud 
!       has been shown to be present by other means.
!       -1 = lwc is present at cloud-top
!       0  = indeterminate
!       +1 = lwc is absent at cloud-top (cloud-ice at cloud-top)

!       we think this works even with thin cloud scenarios that could cause
!       an underestimation of cloud-top temperature when using 'tb8_k'.

        real    t39_k(ni,nj)
        real    tb8_k(ni,nj)
        real    solar_alt(ni,nj)
        integer istat_39_lwc_a(ni,nj)
        integer icount(-1:+1)

        real k_to_c

        do ic = -1,+1
            icount(ic) = 0
        enddo ! ic

        do i = 1,ni
        do j = 1,nj
            if(tb8_k(i,j)  .ne. r_missing_data .and.
     1         t39_k(i,j)  .ne. r_missing_data      )then
                t39_c = k_to_c(t39_k(i,j))
                tb8_c = k_to_c(tb8_k(i,j))

                if(tb8_c .gt. 0.)then ! "definitely" warm cloud top
                    istat_39_lwc_a(i,j) = +1

                elseif(solar_alt(i,j) .lt. 0.        )then
                    if(t39_k(i,j) - tb8_k(i,j) .lt. -3.)then ! sufficient diff
                        istat_39_lwc_a(i,j) = +1
                    else                                     ! insuff diff
                        istat_39_lwc_a(i,j) = -1
                    endif

                else ! sun above horizon and possibly cold cloud top: ambiguous
                    istat_39_lwc_a(i,j) = 0

                endif

            else ! we have missing data
                istat_39_lwc_a(i,j) = 0

            endif

            icount(istat_39_lwc_a(i,j)) = 
     1      icount(istat_39_lwc_a(i,j)) + 1       

        enddo 
        enddo

        write(6,12)icount(-1),icount(0),icount(+1)
 12     format('  3.9u cloud lwc status vals [-1,0,+1] =',3i8)       

        return
        end 

