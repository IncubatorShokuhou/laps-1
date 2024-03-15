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
        subroutine zs(precip_rate,temp_col_max,ni,nj,s_2d_out)

!           1994    steve albers
!       jan 1998    steve albers     remove lapsparms.inc and other cleanup

        real precip_rate(ni,nj)
        real temp_col_max(ni,nj) ! deg k
        real s_2d_out(ni,nj)

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting ref_base in zs'
            stop
        endif

        n_snow_pts = 0
!       n_warm_pts = 0

        do j = 1,nj
        do i = 1,ni

!           if(dbz .eq. ref_base)then
!               s_2d_out(i,j) = +1e-30

!           else
                ratio = snow_to_rain_ratio(temp_col_max(i,j))
                n_snow_pts = n_snow_pts + 1
                s_2d_out(i,j) = precip_rate(i,j) * ratio

                if(n_snow_pts .eq. (n_snow_pts/200) * 200   )then
                    write(6,*)i,j,temp_col_max(i,j)-273.15, ratio
                endif

!           endif

        enddo
        enddo

        write(6,*)' n_snow_pts = ',n_snow_pts

        return
        end

        function snow_to_rain_ratio(temp_col_max)

!       note that temp_col_max is in deg k

        temp_col_c = temp_col_max - 273.15   ! convert from k to c

        if(temp_col_c .ge. 0.0)then          !  t > 0c, use 10

            snow_to_rain_ratio = 10.

        elseif(temp_col_c .ge. -3.0)then     !  0c > t >  -3c, ramp 10 - 15

            frac = temp_col_c/ (-3.0)
            snow_to_rain_ratio = 10. * (1. - frac) + 15. * frac

        elseif(temp_col_c .ge. -10.0)then    !  -3c > t > -10c, ramp 15 - 25

            frac = (temp_col_c - (-3.0)) / (-7.0)
            snow_to_rain_ratio = 15. * (1. - frac) + 25. * frac

        elseif(temp_col_c .ge. -18.0)then    ! -10c > t > -18c, use 25

            snow_to_rain_ratio = 25.

        elseif(temp_col_c .ge. -22.0)then    ! -18c > t > -22c, ramp 25 - 15

            frac = (temp_col_c - (-18.0)) / (-4.0)
            snow_to_rain_ratio = 25. * (1. - frac) + 15. * frac

        else                                 !        t < -22c, use 15

            snow_to_rain_ratio = 15.

        endif

        return
        end
