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
        subroutine cloud_bogus_w (dx, cloud_type, height, nk       ! i
     1                           ,vv_to_height_ratio_cu            ! i
     1                           ,vv_to_height_ratio_sc            ! i
     1                           ,vv_for_st                        ! i
     1                           ,l_deep_vv                        ! i
     1                           ,w)                               ! o

!original version october 1990.

!modified may 1991 when we realized that a grid box is a lot bigger than
!any updraft.  we reduced the maximum vv in the parabolic profiles in cumulus
!clouds by a fairly large amount (30 m/s in 10 km cu to 5 m/s), and the vv max
!for stratocumulus by a smaller amount (50 cm/s in 4-km sc to 20 cm/s).

!  modified june 2002 - once again reduced the maximum cloud vv magnitude
!                       and made it dependent on grid spacing.  also
!                       changed parabolic vv profile for cumuliform clouds
!                       to only go down to the cloud base, rather than 
!                       1/3rd of the cloud depth below base to try and
!                       and improve elevated convection cases.

!can be used with either regular laps analysis grid or the cloud analysis grid.
        implicit none
        integer nk, cloud_type(nk)
        real dx, height(nk), w(nk) ! dx is m and w is m/s

!the following specifies the maximum vv in two cloud types as functions
!of cloud depth.  make parabolic vv profile, except for stratiform clouds,
!which get a constant value.  the values are tuned to give values that
!an nwp model would typically produce.  
        real vv_to_height_ratio_cu
        real vv_to_height_ratio_sc
        real vv_for_st

        real ratio, vv, parabolic_vv_profile

        integer k, k1, kbase, ktop
        real zbase, ztop

        logical l_deep_vv

!   cloud type      /'  ','st','sc','cu','ns','ac','as','cs','ci','cc','cb'/
!   integer value     0     1    2    3    4    5    6    7    8    9   10

!zero out return vector.
        do k = 1, nk
         w(k) = 0.
        end do

!put in the vv's for cumuliform clouds (cu or cb) first.
        ratio = vv_to_height_ratio_cu / dx
        do k = 1, nk
         if (cloud_type(k) .eq. 3  .or.  cloud_type(k) .eq. 10) then
          kbase = k
          go to 10
         end if
        end do
        go to 100

10      do k = kbase, nk
         if(l_deep_vv)then ! we are using adan's change
          if (cloud_type(k) .ne. 0) then ! change to the cloudtop by adan
           ktop = k
          else
           go to 20
          end if
         else ! older strategy with shallower parabolic profiles
          if (cloud_type(k) .eq. 3  .or.  cloud_type(k) .eq. 10) then
           ktop = k
          else
           go to 20
          end if
         endif ! l_deep_vv
        end do

20      k1 = k          ! save our place in the column
        zbase = height(kbase)
        ztop  = height(ktop)
        do k = 1, nk
         vv = parabolic_vv_profile (zbase, ztop, ratio, height(k))
         if (vv .gt. 0.) then
          w(k) = vv
         else
          w(k) = 0. ! could wipe out vv from another layer
         end if
        end do
        k1 = k1 + 1
        if (k1 .ge. nk) go to 100

!try for another level of cu.
        do k = k1, nk
         if (cloud_type(k) .eq. 3  .or.  cloud_type(k) .eq. 10) then
          kbase = k
          go to 10
         end if
        end do

!now do the stratocumulus or similar clouds (sc, ac, cc, ns).
100     ratio = vv_to_height_ratio_sc/dx
        do k = 1, nk
         if (cloud_type(k) .eq. 2  .or.  cloud_type(k) .eq. 4  .or.
     1     cloud_type(k) .eq. 5  .or.  cloud_type(k) .eq. 9) then
          kbase = k
          go to 110
         end if
        end do
        go to 200

110     do k = kbase, nk
         if (cloud_type(k) .eq. 2  .or.  cloud_type(k) .eq. 4  .or.
     1     cloud_type(k) .eq. 5  .or.  cloud_type(k) .eq. 9) then
          ktop = k
         else
          go to 120
         end if
        end do

120     k1 = k          ! save our place in the column
        zbase = height(kbase)
        ztop  = height(ktop)
        do k = 1, nk
         vv = parabolic_vv_profile (zbase, ztop, ratio, height(k))
         if (vv .gt. w(k)) w(k) = vv
        end do
        k1 = k1 + 1
        if (k1 .ge. nk) go to 200       ! try for stratiform clouds

!try for another level of sc.
        do k = k1, nk
         if (cloud_type(k) .eq. 2  .or.  cloud_type(k) .eq. 4  .or.
     1     cloud_type(k) .eq. 5  .or.  cloud_type(k) .eq. 9) then
          kbase = k
          go to 110
         end if
        end do

!make sure there is non-zero vv wherever there are clouds of any kind.
!also, return missing-data value for any non-bogussed vv value.
200     do k = 1, nk
         if (cloud_type(k).ne.0 .and. w(k).lt.vv_for_st) w(k) = vv_for_s
     1t
         if (w(k) .eq. 0.) w(k) = 1e37
        end do

        return
        end

!-------------------------------------------------------------------
        real function parabolic_vv_profile (zbase, ztop, ratio, z)
!the vertical velocity is zero at cloud top, peaks one third of the way up
!from the base, and extends below the base by one third of the cloud depth.

!  june 2002 - no longer extending profile to below cloud base.

        implicit none
        real zbase, ztop, ratio, z
        real depth, vvmax, vvspan, halfspan, height_vvmax, x

        depth = ztop - zbase
        if (depth .le. 0.) then
         parabolic_vv_profile = 0.
         return
        end if

        vvmax = ratio * depth
        vvspan = depth * 1.1    
        halfspan = vvspan / 2.
        height_vvmax = ztop - halfspan
        x = -vvmax/(halfspan*halfspan)

        parabolic_vv_profile = x * (z-height_vvmax)**2 + vvmax

        return
        end


! the variant below is called from 'get_radar_deriv.f/radar_bogus_w'

!-------------------------------------------------------------------
        real function parabolic_vv_profile1 (zbase, ztop, ratio, z)
!the vertical velocity is zero at cloud top, peaks one third of the way up
!from the base, and extends below the base by one third of the cloud depth.

!  june 2002 - no longer extending profile to below cloud base.

        implicit none
        real zbase, ztop, ratio, z
        real depth, vvmax, vvspan, halfspan, height_vvmax, x

        depth = ztop - zbase
        if (depth .le. 0.) then
         parabolic_vv_profile1 = 0.
         return
        end if

        vvmax = ratio * depth
        vvspan = depth
        halfspan = vvspan / 2.
        height_vvmax = ztop - halfspan
        x = -vvmax/(halfspan*halfspan)

        parabolic_vv_profile1 = x * (z-height_vvmax)**2 + vvmax

        return
        end
