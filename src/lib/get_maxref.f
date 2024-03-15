cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis

        subroutine get_max_ref(grid_ra_ref,imax,jmax,kmax,radar_array)

!       this routine converts 3d reflectivity volume data to a 2d reflectivity
!       field. this is done by vertically projecting the maximum reflectivity
!       onto a horizontal surface.

!       steve albers            1990
!       steve albers            1994     test for missing data values

        real grid_ra_ref(imax,jmax,kmax)      ! input 3d array
        real radar_array(imax,jmax)           ! output 2d array

        common /laps_diag/ no_laps_diag

        if(no_laps_diag .eq. 0)
     1   write(6,*)' projecting maximum reflectivity onto horizontal sur
     1face'

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in get_max_ref, stop'
            stop
        endif

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in get_max_ref, stop'
            stop
        endif

!       initialize radar array
        do j = 1,jmax
        do i = 1,imax
            radar_array(i,j) = ref_base
        enddo
        enddo

        do k = 1,kmax
        do j = 1,jmax
        do i = 1,imax
            if(grid_ra_ref(i,j,k) .ne. r_missing_data)
     1       radar_array(i,j) = max(radar_array(i,j),grid_ra_ref(i,j,k))
        enddo
        enddo
        enddo

        return
        end


        subroutine get_max_reflect(grid_ra_ref,imax,jmax,kmax
     1                            ,r_missing_out,radar_array)

!       this routine converts 3d reflectivity volume data to a 2d reflectivity
!       field. this is done by vertically projecting the maximum reflectivity
!       onto a horizontal surface.

!       steve albers            fsl

        real grid_ra_ref(imax,jmax,kmax)      ! input 3d array
        real r_missing_out                    ! input - flag value for 
                                                !         missing output
        real radar_array(imax,jmax)           ! output 2d array

        common /laps_diag/ no_laps_diag

        write(6,*)
     1  ' projecting maximum reflectivity onto horizontal surface'

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in get_max_reflect, stop'
            stop
        endif

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in get_max_reflect, stop'
            stop
        endif

!       initialize radar array
        radar_array = r_missing_data

        do k = 1,kmax
        do j = 1,jmax
        do i = 1,imax
            if(grid_ra_ref(i,j,k) .ne. r_missing_data)then
                if(radar_array(i,j) .eq. r_missing_data)then
                    radar_array(i,j) = grid_ra_ref(i,j,k)
                else
                    radar_array(i,j) = max(radar_array(i,j)
     1                                    ,grid_ra_ref(i,j,k))
                endif
            endif
        enddo
        enddo
        enddo

!       update radar array with proper missing output flag if needed 
        do j = 1,jmax
        do i = 1,imax
            if(r_missing_out    .eq. ref_base)then
                if(     radar_array(i,j) .eq. r_missing_data 
     1             .or. radar_array(i,j) .eq. -101.          ! qc flags
     1             .or. radar_array(i,j) .eq. -102.    )then       
                    radar_array(i,j) = ref_base
                endif

            else ! r_missing_out = r_missing_data
                if(     radar_array(i,j) .eq. -101.          ! qc flags
     1             .or. radar_array(i,j) .eq. -102.    )then       
                    radar_array(i,j) = ref_base
                endif

            endif
        enddo
        enddo

        return
        end

