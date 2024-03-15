
        subroutine get_precip_inc(i4time_beg, i4time_end     ! input
1       , imax, jmax, kmax                             ! input
1       , max_radar_files                            ! input
1       , lat, lon, topo                               ! input
1       , ilaps_cycle_time, grid_spacing_cen_m        ! input
1       , radarext_3d_accum                          ! input
1       , snow_accum, precip_accum, frac_sum           ! outputs
1       , istatus)                                   ! output

        use mem_namelist, only: l_accum_fg, l_accum_radar, l_accum_gauge

!       calculate and return incremental precip for the analysis time cycle

!       input
        real lat(imax, jmax)
        real lon(imax, jmax)
        real topo(imax, jmax)
        real ldf(imax, jmax)

        character*31 radarext_3d_accum

!       output
        real snow_accum(imax, jmax)   ! m
        real precip_accum(imax, jmax) ! m

!       local
        real pcp_bkg_m(imax, jmax)      ! background field for gauge analysis
        real closest_radar(imax, jmax)  ! m
        real temp_col_max(imax, jmax)
        character var_req*4

        call get_r_missing_data(r_missing_data, istatus)
        if (istatus .ne. 1) return

        if (l_accum_radar .eqv. .true.) then
           call get_precip_radar(i4time_beg, i4time_end       ! input
1          , imax, jmax, kmax                             ! input
1          , max_radar_files                            ! input
1          , lat, lon, topo                               ! input
1          , ilaps_cycle_time, grid_spacing_cen_m        ! input
1          , radarext_3d_accum                          ! input
1          , snow_accum, precip_accum, frac_sum           ! output
1          , closest_radar                              ! output
1          , istat_radar)                               ! output
        else
           write (6, *) ' skipping call to get_precip_radar'
           istat_radar = 0
           precip_accum = r_missing_data
        end if

        if (istat_radar .ne. 1) then
           write (6, *) ' istat_radar != 1, not using radar data'
           precip_accum = r_missing_data
        else
           write (6, *) ' radar precip good status'
           write (6, *) ' range of radar ', minval(precip_accum)
1          , maxval(precip_accum)
        end if

!       read precip first guess (lgb/fsf). try 'r01' variable
        if (l_accum_fg .eqv. .true.) then
           var_req = 'r01'
           call get_modelfg_2d(i4time_end, var_req, imax, jmax, pcp_bkg_m
1          , istat_bkg)

           if (istat_bkg .ne. 1) then
              var_req = 'pcp'
              call get_modelfg_2d(i4time_end, var_req, imax, jmax, pcp_bkg_m
1             , istat_bkg)
           end if

        else
           write (6, *) ' skipping call to get_modelfg_2d'
           istat_bkg = 0
        end if

        if (istat_bkg .ne. 1) then
!           write(6,*)' no model first guess precip, using zero field'
!           pcp_bkg_m = 0.
           write (6, *) ' no model first guess precip, set to missing'
           pcp_bkg_m = r_missing_data
        else
           write (6, *) ' first guess precip good status'
           write (6, *) ' range of fg ', minval(pcp_bkg_m)
1          , maxval(pcp_bkg_m)
           write (6, *) ' setting floor of zero'
           pcp_bkg_m = max(pcp_bkg_m, 0.)
           write (6, *) ' new range of fg ', minval(pcp_bkg_m)
1          , maxval(pcp_bkg_m)
        end if

!       compare to gauge values
        if (ilaps_cycle_time .le. 3600 .and.
1       (l_accum_gauge .eqv. .true.)) then
        if (ilaps_cycle_time .lt. 3600) write (6, 51) ilaps_cycle_time/60
51      format('  warning: assuming gauge data read in as 1hr'
1       , ' data is actually every ', i2, ' minutes')

        call get_maxstns(maxsta, istatus)
        if (istatus .ne. 1) return

!           this is set to one to disable the land/water weighting for the
!           precip analysis.
        ldf(:, :) = 1.0

        call blend_gauge_data(i4time_end, imax, jmax, maxsta  ! i
1       , r_missing_data               ! i
1       , lat, lon, topo, ldf             ! i
1       , pcp_bkg_m                    ! i
1       , ilaps_cycle_time             ! i
1       , closest_radar, istat_radar    ! i
1       , precip_accum)                ! i/o
        else
        write (6, '(" skipping gauge processing")')
     end if

!       derive snow accum from overall precip accum field
     write (6, '(/" derive overall snow accum")')
     n_snow_pts = 0
     n_pcp_pts = 0

     do j = 1, jmax
     do i = 1, imax

        if (precip_accum(i, j) .ne. 0. .and.
1       precip_accum(i, j) .ne. r_missing_data) then

        n_pcp_pts = n_pcp_pts + 1
        if (n_pcp_pts .eq. 1) then
           i4_tol = 0
           call get_tcol_max(i4time_end, i4_tol, imax, jmax, kmax, ! i
1          temp_col_max, istatus)             ! o
        end if

        n_snow_pts = n_snow_pts + 1
        ratio = snow_to_rain_ratio(temp_col_max(i, j))
        snow_accum(i, j) = precip_accum(i, j)*ratio

        if (n_snow_pts .le. 50) then
           write (6, '(" snow accum",2i5,2f9.2)') i, j
1          , temp_col_max(i, j) - 273.15, ratio
        end if
        end if

     end do
     end do

     write (6, '(" npts pcp/snow = ",2i5)') n_pcp_pts, n_snow_pts

     return
  end

  subroutine get_precip_radar(i4time_beg, i4time_end   ! input
1 , imax, jmax, kmax                             ! input
1 , max_radar_files                            ! input
1 , lat, lon, topo                               ! input
1 , ilaps_cycle_time, grid_spacing_cen_m        ! input
1 , radarext_3d_accum                          ! input
1 , snow_accum, precip_accum, frac_sum           ! outputs
1 , closest_radar                              ! output
1 , istatus)                                   ! output

!       steve albers 1991
!       steve albers 1995 dec  modify radar call to read_radar_ref
!       returns accumulated snow and liquid precip
!       this routine uses 3d precip type (using temp data from the
!       closest cycle) as a cutoff for snow/no snow. this is calculated and
!       applied every time we are leaving the valid time window for a
!       particular set of environmental data. typically, the first half of the
!       cycle time comes from the previous cycle's temp data, the second half of
!       the period uses current data. a mask is also used so that the 3d
!       precip type calculations are performed only where there are accumulating
!       echoes. in other words, this code is more complex that it would
!       otherwise be so that real-time speed is optimized.

!       input
  real lat(imax, jmax)
  real lon(imax, jmax)
  real topo(imax, jmax)

!       output
  real snow_accum(imax, jmax)   ! m
  real precip_accum(imax, jmax) ! m
  real closest_radar(imax, jmax) ! m

!       local
  real precip_rateave(imax, jmax)
  real snow_rateave(imax, jmax)
  real snow_accum_pd(imax, jmax)
  real snow_rate(imax, jmax) ! m/s
  real precip_rate(imax, jmax) ! m/s
  real dbz_2d(imax, jmax)
  real t_sfc_k(imax, jmax)
  real td_sfc_k(imax, jmax)
  real pres_sfc_pa(imax, jmax)
  real tw_sfc_k(imax, jmax)
  real temp_3d(imax, jmax, kmax)
  real height_3d(imax, jmax, kmax)
  real temp_col_max(imax, jmax)
  real rh_3d(imax, jmax, kmax)
  real pres_3d(imax, jmax, kmax)
  logical l_mask_pcp(imax, jmax)
  logical l_mask_rdr(imax, jmax)
  integer i2_pcp_type_2d(imax, jmax)
  integer i2_cldpcp_type_3d(imax, jmax, kmax)
  integer ipcp_1d(kmax)

  real grid_ra_ref(imax, jmax, kmax)

  character*9 asc_tim_9, asc_tim_9_beg, asc_tim_9_end
  integer i4time_file(max_radar_files)
  real frac(max_radar_files)

  character*4 radar_name ! local

  character*3 var_2d
  character*3 ext_local
  character*31 ext, radarext_3d_accum
  character*10 units_2d
  character*125 comment_2d

  integer iarg

  logical l_first_sfc_update_completed

  character*80 grid_fnam_common
  common/grid_fnam_cmn/grid_fnam_common

  data mode_radar/1/

  twet_snow = +1.3

  max_radar_gap = float(ilaps_cycle_time) + 1200. ! * 1.33334

  im = imax/2 ! 15
  jm = jmax/2 ! 56

  call make_fnam_lp(i4time_beg, asc_tim_9_beg, istatus)
  call make_fnam_lp(i4time_end, asc_tim_9_end, istatus)

  call get_r_missing_data(r_missing_data, istatus)
  if (istatus .ne. 1) return

  write (6, *) ' radar accumulation from ', asc_tim_9_beg,
1 ' to ', asc_tim_9_end

!       get file times and fractional time periods
  i4time_now = i4time_now_gg()
  i_wait = (i4time_end + (35*60) - i4time_now)/60
  write (6, *) ' number of potential wait cycles = ', i_wait

50 continue

  if (radarext_3d_accum(1:3) .eq. 'xxx') then ! automatically select type(s)
     ext_local = 'vrc'
     call get_fracs(i4time_beg, i4time_end, max_radar_gap        ! i
1    , i_nbr_files_ret                            ! i
1    , ext_local                                  ! i
1    , max_radar_files                            ! i
1    , nscans_vrc                                 ! o
1    , i4time_file                                ! 0
1    , frac, frac_sum, istatus)                     ! o

     write (6, *) ' ext_local / nscans = ', ext_local, ' ', nscans_vrc

     ext_local = 'vrz'
     call get_fracs(i4time_beg, i4time_end, max_radar_gap        ! i
1    , i_nbr_files_ret                            ! i
1    , ext_local                                  ! i
1    , max_radar_files                            ! i
1    , nscans_vrz                                 ! o
1    , i4time_file                                ! 0
1    , frac, frac_sum, istatus)                     ! o

     write (6, *) ' ext_local / nscans = ', ext_local, ' ', nscans_vrz

     if (nscans_vrc .gt. nscans_vrz) then
        ext_local = 'vrc'
     elseif (nscans_vrc .lt. nscans_vrz) then
        ext_local = 'vrz'
     else ! equal case
        ext_local = 'vrc' ! or 'vrz'
     end if

     write (6, *)
     write (6, *) ' auto select: nscans_vrc/nscans_vrz/ext_local:'
1    , nscans_vrc, nscans_vrz, ' ', ext_local
     write (6, *)

  else
     ext_local = radarext_3d_accum(1:3)

  end if

  call get_fracs(i4time_beg, i4time_end, max_radar_gap            ! i
1 , i_nbr_files_ret                                ! i
1 , ext_local                                      ! i
1 , max_radar_files                                ! i
1 , nscans                                         ! o
1 , i4time_file                                    ! 0
1 , frac, frac_sum, istatus)                         ! o

  write (6, *) ' ext_local / nscans = ', ext_local, ' ', nscans

  if (istatus .ne. 1) then
     if (i_wait .gt. 0 .and. frac_sum .ge. 0.30) then
        write (6, *) ' waiting 1 min for possible radar data'
1       , frac_sum
        call snooze_gg(60.0, istat_snooze)
        i_wait = i_wait - 1
        goto50

     else ! will not wait for radar data
        write (6, *) ' warning: insufficient data for accum'
        write (6, *) ' no accum output being generated'
        return

     end if
  end if

!       initialize both period and total precip/snow

!       snow is accumulated in "periods" corresponding to when the ancillary
!       laps data is read in. each time new ancillary data is read in, the
!       "period" snow accumulation is added to the overall snow total for the
!       window.

!       this is done so that precip type (which determines snow/no snow) only
!       has to be calculated only once for every time ancillary data is read in.
!       the precip type is calculated over a masked area corresponding to where
!       precip (potentially snow) has occurred over the "period".

!       this approximation may lead to slightly inaccurate accumulations
!       over the period if the precip type is changing over the period.
!       note that precip type is a function of temperature, rh, and even
!       reflectivity - all of which can change in various ways over the
!       laps cycle. the approximation is there to increase the speed of
!       execution. when computers get better, it may be desirable to
!       redesign this so the precip type (and snow/rain ratio) gets calculated
!       for each radar scan. of course this can also be taken care of by
!       decreasing the overall laps analysis cycle time.

  do j = 1, jmax
  do i = 1, imax
     snow_rateave(i, j) = 0.
     snow_accum_pd(i, j) = 0.
     precip_rateave(i, j) = 0.
     l_mask_pcp(i, j) = .false.
     l_mask_rdr(i, j) = .false.
  end do ! i
  end do ! j

  i4_interval = (i4time_end - i4time_beg)
  i4time_mid = i4time_beg + i4_interval/2
  i4_tol = max(float(ilaps_cycle_time)*0.6, 20.)
  i4time_temp = 0

  l_first_sfc_update_completed = .false.

!       loop through all the radar scan times
  do ifile = 1, i_nbr_files_ret
     if (frac(ifile) .gt. 0.) then ! this scan is needed for the window at hand

        i4time_radar = i4time_file(ifile)

        call make_fnam_lp(i4time_radar, asc_tim_9, istatus)

!           write(6,101)asc_tim_9,frac(ifile)
!101        format(' time, frac = ',a9,2x,f6.3)

        write (6, *)

!           determine whether we need to update the sfc data to match the radar
!           and if we should add the period snow accum onto the total
        if (abs(i4time_radar - i4time_temp) .gt. ilaps_cycle_time/2
1       ) then

!               determine whether to add to snow accum total based on precip
!               type using ancillary laps data valid for the interval just
!               ended. note that the first time around, we are just
!               initializing and don't need to add to the snow accum, we just
!               read in the initial ancillary laps data.

        if (l_first_sfc_update_completed) then
           do j = 1, jmax
           do i = 1, imax
              if (snow_accum_pd(i, j) .gt. 1e-10) then
                 l_mask_pcp(i, j) = .true.
              end if
           end do
           end do

           write (6, *) ' compute 3d precip type over masked area'

!                   note that the reflectivities here are valid at the end
!                   of the accumulating period. this means that the wb
!                   threshold is calculated using reflectivities that may not
!                   be representative for the entire accumulation subperiod
!                   (typically one half of the laps cycle time).

           call cpt_pcp_type_3d(temp_3d, rh_3d, pres_3d
1          , grid_ra_ref, l_mask_pcp, grid_spacing_cen_m
1          , imax, jmax, kmax, twet_snow, i2_cldpcp_type_3d
1          , istatus)
           if (istatus .ne. 1) then
              return
           end if

           i4_elapsed = ishow_timer()

           write (6, *) ' compute sfc precip type'
           call get_sfc_preciptype(pres_sfc_pa, t_sfc_k, td_sfc_k
1          , i2_cldpcp_type_3d, twet_snow, dbz_2d
1          , i2_pcp_type_2d, imax, jmax, kmax)

           i4_elapsed = ishow_timer()

           write (6, *)
1          ' adding in accumulation for intervening period'

           n_pcp_pts = 0
           n_snw_pts = 0
           n_nopcp_pts = 0
           n_zr_pts = 0

           do j = 1, jmax
           do i = 1, imax

              if (l_mask_pcp(i, j)) then
                 n_pcp_pts = n_pcp_pts + 1
                 iarg = i2_pcp_type_2d(i, j)/16

                 if (iarg .eq. 2) then     ! precip type is snow
                    r_pcp_type = iarg
                    r_pcp_type_thresh = iarg

                    if (.true.) then  ! testcode
                       call nowrad_virga_correction(
1                      r_pcp_type,
1                      r_pcp_type_thresh,
1                      t_sfc_k(i, j),
1                      td_sfc_k(i, j),
1                      istatus_3dref)
                    end if

                    if (r_pcp_type_thresh .eq. 2.0) then
                       snow_rateave(i, j) =
1                      snow_rateave(i, j) + snow_accum_pd(i, j)
                       n_snw_pts = n_snw_pts + 1
                    else
                       write (6, *)
1                      ' nowrad_virga thresh out #1', i, j
                    end if

                 elseif (iarg .eq. 0) then ! precip type is no precip
                    n_nopcp_pts = n_nopcp_pts + 1
                    if (n_nopcp_pts .le. 20)
1                   write (6, *) ' no2dpcptype pt at', i, j

                 elseif (iarg .eq. 3) then ! precip type is freezing rain
                    n_zr_pts = n_zr_pts + 1

                 end if

              end if
           end do ! i
           end do ! j

!                   this is useful for debugging
           do k = 1, kmax
              iarg = i2_cldpcp_type_3d(im, jm, k)/16
              ipcp_1d(k) = iarg
           end do ! k
           write (6, 101) l_mask_pcp(im, jm), snow_accum_pd(im, jm)
1          , snow_rateave(im, jm)
1          , t_sfc_k(im, jm), td_sfc_k(im, jm)
1          , i2_pcp_type_2d(im, jm)
1          , (ipcp_1d(k), k=1, min(kmax, 10))
101        format(1x, 'accum/tw/type2,3d', l2, 2e11.3, 2f6.1, i3, 1x,
1          10i2)

           write (6, *) ' # of points snow/precip/zr = '
1          , n_snw_pts, n_pcp_pts, n_zr_pts

           if (n_nopcp_pts .gt. 0)
1          write (6, *) ' warning: n_nopcp_pts = ', n_nopcp_pts

        end if ! l_first_sfc_update_completed

!               initialize
        do j = 1, jmax
        do i = 1, imax
           snow_accum_pd(i, j) = 0.
           l_mask_pcp(i, j) = .false.
        end do ! i
        end do ! j

!               read in surface data that is time matched to the radar data
        write (6, *)
        write (6, *) ' updating surface and 3d temp information'

!               read in surface temp data
        var_2d = 't'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_radar, i4_tol, i4time_temp
1       , ext, var_2d, units_2d, comment_2d
1       , imax, jmax, t_sfc_k, 0, istatus)
        if (istatus .ne. 1) then
           write (6, *) ' warning: laps sfc temp not available'
           frac_sum = -1.0 ! turns off the wait loop for more radar
           return
        end if

!               read in surface dewpoint data
        var_2d = 'td'
        ext = 'lsx'
        call get_laps_2d(i4time_temp, ext, var_2d
1       , units_2d, comment_2d, imax, jmax, td_sfc_k, istatus)
        if (istatus .ne. 1) then
           write (6, *) ' warning: laps sfc dewpoint'
1          , ' not available'
           frac_sum = -1.0 ! turns off the wait loop for more radar
           return
        end if

!               read in surface pressure data
        var_2d = 'ps'
        ext = 'lsx'
        call get_laps_2d(i4time_temp, ext, var_2d
1       , units_2d, comment_2d, imax, jmax, pres_sfc_pa, istatus)
        if (istatus .ne. 1) then
           write (6, *) ' warning: laps sfc pressure '
1          , 'not available', istatus
           frac_sum = -1.0 ! turns off the wait loop for more radar
           write (6, *) ' range of pres_sfc_pa'
1          , minval(pres_sfc_pa), maxval(pres_sfc_pa)
           do i = 1, imax
           do j = 1, jmax
              if (pres_sfc_pa(i, j) .eq. r_missing_data) then
                 write (6, *) i, j, pres_sfc_pa(i, j), topo(i, j)
              end if
           end do ! j
           end do ! i
           return
        end if

        i4_elapsed = ishow_timer()

        write (6, *) ' getting lt1 height/temp'
        var_2d = 'ht'
        ext = 'lt1'

        call get_laps_3d(i4time_temp
1       , imax, jmax, kmax, ext, var_2d
1       , units_2d, comment_2d, height_3d, istatus)
        if (istatus .ne. 1) then
           write (6, *) ' warning: laps 3d height not available'
           frac_sum = -1.0 ! turns off the wait loop for more radar
           return
        end if

        var_2d = 't3'
        ext = 'lt1'

        call get_laps_3d(i4time_temp
1       , imax, jmax, kmax, ext, var_2d
1       , units_2d, comment_2d, temp_3d, istatus)
        if (istatus .ne. 1) then
           write (6, *) ' warning: laps 3d temp not available'
           frac_sum = -1.0 ! turns off the wait loop for more radar
           return
        end if

!               calculate column max temperatures
        do j = 1, jmax
        do i = 1, imax
           temp_col_max(i, j) = t_sfc_k(i, j)
           k_sfc = int(zcoord_of_pressure(pres_sfc_pa(i, j)))
           do k = k_sfc + 1, kmax
              temp_col_max(i, j) =
1             max(temp_col_max(i, j), temp_3d(i, j, k))
           end do ! k
        end do ! i
        end do ! j

        write (6, *) ' getting lh3 file (or equivalent)'
        var_2d = 'rhl'
        ext = 'lh3'
        call get_laps_3dgrid(i4time_temp, ilaps_cycle_time ! *2
1       , i4time_rh
1       , imax, jmax, kmax, ext, var_2d
1       , units_2d, comment_2d, rh_3d, istatus)
        if (istatus .ne. 1) then
           write (6, *) ' warning: laps 3d rh not available'
           frac_sum = -1.0 ! turns off the wait loop for more radar
           return
        end if

        l_first_sfc_update_completed = .true.

        write (6, *) ' ancillary laps data update completed'
        write (6, *)

        i4_elapsed = ishow_timer()

     else ! l_first_sfc_update_completed = .true.
        write (6, *) ' no ancillary laps data update needed'

     end if ! we need to update the laps information

!           get laps reflectivities at the surface (or immediately above it)
     write (6, *)

     i4_elapsed = ishow_timer()

!           repeat read in radar data with low level reflectivities filled in
     call get_ref_base(ref_base, istatus)
     if (istatus .ne. 1) return

     i4_tol_radar = 1200

     ext = '                               '
     ext(1:3) = ext_local

     if (ext(1:3) .ne. 'lmt') then

        call read_radar_3dref_new(i4time_radar,             ! i
1       i4_tol_radar, i4_ret,                               ! i/o
1       .true., r_missing_data, imax, jmax, kmax,              ! i
1       ext, lat, lon, topo,
1       .true., .false.,
1       height_3d,
1       grid_ra_ref,
1       closest_radar,                                     ! o
1       rlat_radar, rlon_radar, rheight_radar, radar_name,
1       n_ref, istatus_2dref, istatus_3dref)

        if (istatus_2dref .eq. 0) then
           write (6,
1          '(" error in reading radar data in precip accum")')
           frac_sum = -1.0 ! turns off the wait loop for more radar
           istatus = 0
           return
        end if

        i4_elapsed = ishow_timer()

!               for now, we can change the 'r_missing_data' values to 'ref_base'
        do i = 1, imax
        do j = 1, jmax
        do k = 1, kmax
           if (grid_ra_ref(i, j, k) .eq. r_missing_data) then
              grid_ra_ref(i, j, k) = ref_base
           else
              l_mask_rdr(i, j) = .true.
           end if
        end do ! k
        end do ! j
        end do ! i

        write (6, *) ' call get_low_ref'

        call get_low_ref(grid_ra_ref, pres_sfc_pa, imax, jmax, kmax
1       , dbz_2d)
     else ! lmt/llr
!               read analyzed low level reflectivity
        var_2d = 'llr'
        ext = 'lmt'
        call get_laps_2d(i4time_radar, ext, var_2d
1       , units_2d, comment_2d, imax, jmax, dbz_2d, istatus)
        l_mask_rdr(:, :) = .true.
     end if

     write (6, *) ' incrementing precip accumulation '
1    , 'rate for this scan (call zr)'

     call zr(dbz_2d, imax, jmax, precip_rate)

     do j = 1, jmax
     do i = 1, imax
        if (precip_rate(i, j) .ne. r_missing_data .and.
1       l_mask_rdr(i, j) .and.
1       precip_rateave(i, j) .ne. r_missing_data) then
        precip_rateave(i, j) = precip_rateave(i, j)
1       +precip_rate(i, j)*frac(ifile)
        else
        precip_rateave(i, j) = r_missing_data
        end if
     end do ! i
     end do ! j

     i4_elapsed = ishow_timer()

     write (6, *) ' incrementing snow accumulation_pd for this scan'

     call zs(precip_rate, temp_col_max, imax, jmax, snow_rate)

     do j = 1, jmax
     do i = 1, imax
        if (snow_rate(i, j) .ne. r_missing_data .and.
1       l_mask_rdr(i, j) .and.
1       snow_accum_pd(i, j) .ne. r_missing_data) then
        snow_accum_pd(i, j) = snow_accum_pd(i, j)
1       +snow_rate(i, j)*frac(ifile)
        else
        snow_accum_pd(i, j) = r_missing_data
        end if
     end do ! i
     end do ! j

     write (6, 202) dbz_2d(im, jm), snow_rate(im, jm)
1    , snow_accum_pd(im, jm)
202  format(1x, 'dbz/rate/accum_pd', 3e12.4)

  end if ! frac > 0 (this radar file is within the accumulating window

!         write(6,*)' cycle to next file'

end do ! ifile (loop through all radar file times)

!       add in snow accumulation from final period; first define the mask
do j = 1, jmax
do i = 1, imax
   if (snow_accum_pd(i, j) .gt. 1e-10 .and.
1  snow_accum_pd(i, j) .ne. r_missing_data) then
   l_mask_pcp(i, j) = .true.
   end if
end do
end do

write (6, *) ' compute 3d precip type over masked area'

call cpt_pcp_type_3d(temp_3d, rh_3d, pres_3d, grid_ra_ref
1 , l_mask_pcp, grid_spacing_cen_m
1 , imax, jmax, kmax, twet_snow, i2_cldpcp_type_3d, istatus)
if (istatus .ne. 1) then
   return
end if

write (6, *) ' compute sfc precip type'
call get_sfc_preciptype(pres_sfc_pa, t_sfc_k, td_sfc_k
1 , i2_cldpcp_type_3d, twet_snow, dbz_2d, i2_pcp_type_2d
1 , imax, jmax, kmax)

write (6, *) ' adding in accumulation for the last period'

n_pcp_pts = 0
n_snw_pts = 0
n_nopcp_pts = 0
n_zr_pts = 0

do j = 1, jmax
do i = 1, imax

   if (l_mask_pcp(i, j)) then
      n_pcp_pts = n_pcp_pts + 1
      iarg = i2_pcp_type_2d(i, j)/16

      if (iarg .eq. 2) then     ! precip type is snow

         r_pcp_type = iarg
         r_pcp_type_thresh = iarg

         if (.true.) then ! testcode
            call nowrad_virga_correction(
1           r_pcp_type,
1           r_pcp_type_thresh,
1           t_sfc_k(i, j),
1           td_sfc_k(i, j),
1           istatus_3dref)
         else
            write (6, *) ' nowrad_virga thresh out #2', i, j
         end if

         if (r_pcp_type_thresh .eq. 2.0) then
            snow_rateave(i, j) = snow_rateave(i, j)
1           +snow_accum_pd(i, j)
            n_snw_pts = n_snw_pts + 1
         end if

      elseif (iarg .eq. 0) then ! precip type is no precip
         n_nopcp_pts = n_nopcp_pts + 1
         if (n_nopcp_pts .le. 20) write (6, *)
1        ' no2dpcptype pt at ', i, j
      elseif (iarg .eq. 3) then ! precip type is freezing rain
         n_zr_pts = n_zr_pts + 1
      end if

   end if

end do ! i
end do ! j

!       this is useful for debugging
do k = 1, kmax
   iarg = i2_cldpcp_type_3d(im, jm, k)/16
   ipcp_1d(k) = iarg
end do ! k
write (6, 101) l_mask_pcp(im, jm), snow_accum_pd(im, jm)
1 , snow_rateave(im, jm)
1 , t_sfc_k(im, jm), td_sfc_k(im, jm)
1 , i2_pcp_type_2d(im, jm)
1 , (ipcp_1d(k), k=1, min(kmax, 10))

write (6, *) ' # of points snow/precip/zr = ', n_snw_pts, n_pcp_pts
1 , n_zr_pts

if (n_nopcp_pts .gt. 0) write (6, *) ' warning: n_nopcp_pts = '
1 , n_nopcp_pts

write (6, *) ' converting from average rate to actual accumulation'

!       convert from time averaged rate to accumulation
do j = 1, jmax
do i = 1, imax
   if (snow_rateave(i, j) .ne. r_missing_data) then
      snow_accum(i, j) = snow_rateave(i, j)*i4_interval
   else
      snow_accum(i, j) = r_missing_data
   end if

   if (precip_rateave(i, j) .ne. r_missing_data) then
      precip_accum(i, j) = precip_rateave(i, j)*i4_interval
   else
      precip_accum(i, j) = r_missing_data
   end if
end do ! i
end do ! j

write (6, *) ' radar snow_accum(im,jm) = ', snow_accum(im, jm)
write (6, *) ' radar precip_accum(im,jm) = ', precip_accum(im, jm)

return
end

subroutine get_fracs(i4time_beg, i4time_end, max_radar_gap        ! i
1 , i_nbr_files_ret                            ! i
1 , ext                                        ! i
1 , max_radar_files                            ! i
1 , nscans                                     ! o
1 , i4time_file, frac, frac_sum, istatus)         ! o

!       steve albers    1991    this routine calculates linear combination
!                               coefficients for radar scans which can be used
!                               to arrive at an integrated precipitation rate
!                               over a specified time window

integer max_radar_files

character*255 c_filespec
character*9 asc_tim_9
character c_fnames(max_radar_files)*80
character*3 ext

real frac(max_radar_files)
integer i4time_file(max_radar_files)

nscans = 0 ! initialize

call get_filespec(ext, 2, c_filespec, istatus)

call get_file_names(c_filespec,
1 i_nbr_files_ret,
1 c_fnames,
1 max_radar_files,
1 i_status)

if (i_nbr_files_ret .gt. 0) then
   call get_directory_length(c_fnames(1), lenf)
else ! error condition
   write (6, *) ' warning: no radar data available for'
1  , ' snow/precip accumulation ', ext
   istatus = 0
   return
end if

do i = 1, i_nbr_files_ret
   asc_tim_9 = c_fnames(i) (lenf + 1:lenf + 9)
   call i4time_fname_lp(asc_tim_9, i4time_file(i), istatus)
   if (istatus .ne. 1) then
      write (6, *) ' bad return from i4time_fname_lp,',
1     ' called from get_precip_accum: ', asc_tim_9
      return
   end if
end do

i4_interval = i4time_end - i4time_beg
nscans = 0

do i = 1, i_nbr_files_ret
!           write(6,301)i4time_beg,i4time_file(i),i4time_end
!301        format(1x,3i12)
   frac(i) = 0.

   if (i4time_file(i) .ge. i4time_beg .and.
1  i4time_file(i) .lt. i4time_end) then
   nscans = nscans + 1
   end if
end do

do i = 1, i_nbr_files_ret - 1
   ibeg = i
   iend = i + 1

   interval_between_files = i4time_file(iend)
1  -i4time_file(ibeg)

   frac_between_files = float(interval_between_files)
1  /float(i4_interval)

   if (i4time_file(ibeg) .ge. i4time_beg .and.
1  i4time_file(iend) .le. i4time_end) then ! fully within pd
   frac(ibeg) = frac(ibeg) + 0.5*frac_between_files
   frac(iend) = frac(iend) + 0.5*frac_between_files

   if (interval_between_files .gt. max_radar_gap) then
      write (6, *) ' warning: gap in radar files (min) >'
1     , max_radar_gap/60
      istatus = 0
      return
   end if

   end if

   if (i4time_file(ibeg) .lt. i4time_beg .and.
1  i4time_file(iend) .gt. i4time_beg .and.
1  i4time_file(iend) .le. i4time_end) then ! straddle beginning
   a = i4time_beg - i4time_file(ibeg)
   b = i4time_file(iend) - i4time_beg
   partial_frac = b/(a + b)
   frac(ibeg) = frac(ibeg) + (0.5*partial_frac)
1  *frac_between_files*partial_frac
   frac(iend) = frac(iend) + (1.0 - 0.5*partial_frac)
1  *frac_between_files*partial_frac

   if (interval_between_files .gt. max_radar_gap) then
      write (6, *) ' warning: gap in radar files (min) >'
1     , max_radar_gap/60
      istatus = 0
      return
   end if

end if

if (i4time_file(ibeg) .lt. i4time_end .and.
1 i4time_file(ibeg) .ge. i4time_beg .and.
1 i4time_file(iend) .gt. i4time_end) then ! straddle end
a = i4time_end - i4time_file(ibeg)
b = i4time_file(iend) - i4time_end
partial_frac = a/(a + b)
frac(ibeg) = frac(ibeg) + (1.0 - 0.5*partial_frac)
1 *frac_between_files*partial_frac
frac(iend) = frac(iend) + (0.5*partial_frac)
1 *frac_between_files*partial_frac

if (interval_between_files .gt. max_radar_gap) then
   write (6, *) ' warning: gap in radar files (min) >'
1  , max_radar_gap/60
   istatus = 0
   return
end if

end if

if (i4time_file(ibeg) .lt. i4time_beg .and.
1 i4time_file(iend) .gt. i4time_end) then ! brackets the pd
i4time_mid = i4time_beg + (i4time_end - i4time_beg)/2
frac_mid = float(i4time_mid - i4time_file(ibeg))
1 /float(i4time_file(iend) - i4time_file(ibeg))
frac(ibeg) = 1.0 - frac_mid
frac(iend) = frac_mid

if (interval_between_files .gt. max_radar_gap) then
   write (6, *) ' warning: gap in radar files (min) >'
1  , max_radar_gap/60
   istatus = 0
   return
end if

end if

end do ! i

frac_sum = 0.
do ifile = 1, i_nbr_files_ret
   if (frac(ifile) .gt. 0.) then
      call make_fnam_lp(i4time_file(ifile), asc_tim_9, istat_fnam)
      write (6, 101) asc_tim_9, frac(ifile)
101   format('        filetime, frac = ', a9, 2x, f8.5)
      frac_sum = frac_sum + frac(ifile)
   end if
end do

if (abs(frac_sum - 1.0) .gt. 1e-5) then
!           note: we can here potentially wait for more radar data in driver
   write (6, *) ' warning: fractions do not add up to 1.0'
1  , frac_sum
   istatus = 0
end if

if (i_nbr_files_ret .gt. 0) then
   if (i4time_file(1) .gt. i4time_beg) then
      write (6, *)
1     ' radar files begin after start of accumulation window'
      frac_sum = -1.0 ! turns off the wait loop for more radar
   end if
end if

999 if (istatus .ne. 1) then
   write (6, *) ' insufficient files within time window'
   return
end if

istatus = 1
return
end

subroutine get_tcol_max(i4time_radar, i4_tol, imax, jmax, kmax,  ! i
1 temp_col_max, istatus)                ! o

real t_sfc_k(imax, jmax)
real td_sfc_k(imax, jmax)
real pres_sfc_pa(imax, jmax)
real temp_3d(imax, jmax, kmax)
real temp_col_max(imax, jmax)

character*3 var_2d
character*3 ext_local
character*31 ext
character*10 units_2d
character*125 comment_2d

!       read in surface data that is time matched to the radar data
write (6, *)
write (6, *) ' updating surface and 3d temp information'

!       read in surface temp data
var_2d = 't'
ext = 'lsx'
call get_laps_2dgrid(i4time_radar, i4_tol, i4time_temp
1 , ext, var_2d, units_2d, comment_2d
1 , imax, jmax, t_sfc_k, 0, istatus)
if (istatus .ne. 1) then
   write (6, '(" warning: laps sfc temperature not available")')
   frac_sum = -1.0 ! turns off the wait loop for more radar
   return
end if

!       read in surface dewpoint data
var_2d = 'td'
ext = 'lsx'
call get_laps_2d(i4time_temp, ext, var_2d
1 , units_2d, comment_2d, imax, jmax, td_sfc_k, istatus)
if (istatus .ne. 1) then
   write (6, '(" warning: laps sfc dewpoint not available")')
   frac_sum = -1.0 ! turns off the wait loop for more radar
   return
end if

!       read in surface pressure data
var_2d = 'ps'
ext = 'lsx'
call get_laps_2d(i4time_temp, ext, var_2d
1 , units_2d, comment_2d, imax, jmax, pres_sfc_pa, istatus)
if (istatus .ne. 1) then
   write (6, '(" warning: laps sfc pressure not available"
1  , i3) ')istatus
   frac_sum = -1.0 ! turns off the wait loop for more radar
   write (6, *) ' range of pres_sfc_pa'
1  , minval(pres_sfc_pa), maxval(pres_sfc_pa)
   do i = 1, imax
   do j = 1, jmax
      if (pres_sfc_pa(i, j) .eq. r_missing_data) then
         write (6, '(2i4,f9.3)') i, j, pres_sfc_pa(i, j)
      end if
   end do ! j
   end do ! i
   return
end if

i4_elapsed = ishow_timer()

var_2d = 't3'
ext = 'lt1'

call get_laps_3d(i4time_temp
1 , imax, jmax, kmax, ext, var_2d
1 , units_2d, comment_2d, temp_3d, istatus)
if (istatus .ne. 1) then
   write (6, *) ' warning: laps 3d temp not available'
   frac_sum = -1.0 ! turns off the wait loop for more radar
   return
end if

!       calculate column max temperatures
do j = 1, jmax
do i = 1, imax
   temp_col_max(i, j) = t_sfc_k(i, j)
   k_sfc = int(zcoord_of_pressure(pres_sfc_pa(i, j)))
   do k = k_sfc + 1, kmax
      temp_col_max(i, j) =
1     max(temp_col_max(i, j), temp_3d(i, j, k))
   end do ! k
end do ! i
end do ! j

return
end

