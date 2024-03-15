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
        subroutine get_low_ref(z_3d_in,pres_sfc_pa,ni,nj,nk
     1                        ,radar_2d_out)

!       steve albers            1990

        real z_3d_in(ni,nj,nk)
        real radar_2d_out(ni,nj)
        real pres_sfc_pa(ni,nj)

        write(6,*)' converting from 3d z field to low level z field'

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in get_low_ref, stop'
            stop
        endif

        n_low_pts = 0

        do j = 1,nj
        do i = 1,ni

!           k_topo = max(int(height_to_zcoord(topo(i,j),istatus)),1)

            k_topo = int(zcoord_of_pressure(pres_sfc_pa(i,j)))
            radar_2d_out(i,j) = z_3d_in(i,j,k_topo+1)
            if(radar_2d_out(i,j) .ne. ref_base)then
                n_low_pts = n_low_pts + 1
            endif

        enddo
        enddo

        write(6,*)' n_low_pts  = ',n_low_pts

        return
        end
