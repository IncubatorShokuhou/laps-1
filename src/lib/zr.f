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
        subroutine zr(z_2d_in
!       1       ,temp_sfc_k,td_sfc_k,pres_sta_pa,tw_sfc_k
     1                                          ,ni,nj,r_2d_out)

!           1991                                        steve albers, j. smart
!       apr 1994 rate_max parameter disabled            sa
!       jan 1998 remove lapsparms.inc                   sa

        real z_2d_in(ni,nj)
!       real temp_sfc_k(ni,nj)
!       real td_sfc_k(ni,nj)
!       real pres_sta_pa(ni,nj)
!       real tw_sfc_k(ni,nj)
        real r_2d_out(ni,nj)

        real a,b,rate_max
        parameter (a = 200.)        ! z-r relationship
        parameter (b = 1.6)         ! z-r relationship
!       parameter (rate_max = 10.0) ! mm/hr; equiv to 39 dbz with z=200r**1.6
        parameter (rate_max = 1000.0) ! disabled

        write(6,*)
     1   ' converting from 2d z to rainfall rate field'

        call get_ref_base_useable(ref_base_useable,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting ref_base_useable in zr'
            stop
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting r_missing_data in zr'
            stop
        endif

        n_leq_pts = 0
!       n_warm_pts = 0

        aterm = alog10(1./a)
        bterm = 1./b
        cterm = .001 / 3600  ! (m/s) / (mm/hr)

        do j = 1,nj
        do i = 1,ni
            dbz = z_2d_in(i,j)

            if(dbz .eq. r_missing_data)then
                r_2d_out(i,j) = r_missing_data

            elseif(dbz .lt. ref_base_useable)then
!               r_2d_out(i,j) = +1e-30
                r_2d_out(i,j) = 0.

            else
                n_leq_pts = n_leq_pts + 1

              ! generate r (mm/hr) in z=a*r**b
                r_mm_hr = 10.**(bterm*(aterm + dbz/10.))
                r_2d_out(i,j) = min(r_mm_hr,rate_max) * cterm

!               if(tw_sfc_k(i,j) .gt. 273.65)then ! wet bulb > 0.5 deg c
!                   r_2d_out(i,j) = 0.
!                   n_warm_pts = n_warm_pts + 1
!               endif

            endif

        enddo
        enddo

        write(6,*)' n_leq_pts = ',n_leq_pts

        return
        end
