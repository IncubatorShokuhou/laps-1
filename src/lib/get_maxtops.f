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

        subroutine get_maxtops(grid_ra_ref,heights_3d
     1                        ,imax,jmax,kmax,rmax_tops_m)

        real grid_ra_ref(imax,jmax,kmax)             ! i
        real heights_3d(imax,jmax,kmax)              ! i
        real rmax_tops_m(imax,jmax)                  ! o

        highest_top = 0.
        rf_thr = 0.
        do j = 1,jmax
        do i = 1,imax

!           initialize point
            rmax_tops_m(i,j) = 0.

!           find level of max tops
            k_grid = 0.
            do k = 1,kmax-1
                if(grid_ra_ref(i,j,k)   .ge. rf_thr .and.
     1             grid_ra_ref(i,j,k+1) .lt. rf_thr)then
                    k_grid = k
                    frac_k = (rf_thr - grid_ra_ref(i,j,k))
     1             / (grid_ra_ref(i,j,k+1) - grid_ra_ref(i,j,k))
                endif
            enddo

!           convert from grid level to height
            if(k_grid .gt. 0.)then
                height_low  = heights_3d(i,j,k_grid)
                height_high = heights_3d(i,j,k_grid+1)
                rmax_tops_m(i,j) = height_low  * (1. - frac_k)
     1                           + height_high * frac_k
            endif

            highest_top = max(highest_top,rmax_tops_m(i,j))

        enddo ! j
        enddo ! i


        return
        end

