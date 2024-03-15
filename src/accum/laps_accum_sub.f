
        subroutine laps_accum_sub(i4time,
1       nx_l, ny_l, nz_l, max_radar_files,
1       i_diag,
1       n_prods,
1       iprod_number,
1       j_status)

!       1991            steve albers
!       1996 may        steve albers  call 'get_laps_cycle_time'
!                                     use 'ilaps_cycle_time' to calculate
!                                     'precip_reset_thresh'.
!       1997 mar        steve albers  added calls to 'move' to replace
!                                     equivalencing.
!       1997 jun        ken dritz     made nx_l, ny_l, nz_l dummy arguments,
!                                     making non-dummy arrays dimensioned
!                                     therewith dynamic (automatic).
!       1997 jun        ken dritz     changed include of lapsparms.for to
!                                     include of laps_static_parameters.inc.

        use mem_namelist, only: precip_cycle_time

        integer ss_normal, rtsys_bad_prod, rtsys_no_data
1       , rtsys_abort_prod
        parameter(ss_normal=1, ! success
1       rtsys_bad_prod = 2, ! inappropriate data, insufficient data
1       rtsys_no_data = 3, ! no data
1       rtsys_abort_prod = 4) ! failed to make a prod

        include 'laps_static_parameters.inc'

        integer i4time, i_diag, n_prods

        real lat(nx_l, ny_l), lon(nx_l, ny_l)
        real topo(nx_l, ny_l)

        real snow_2d(nx_l, ny_l)
        real snow_2d_tot(nx_l, ny_l)
        real snow_2d_depth(nx_l, ny_l)
        real snow_cover(nx_l, ny_l)

        real precip_2d(nx_l, ny_l)
        real precip_2d_tot(nx_l, ny_l)

        integer nfields
        parameter(nfields=5)

        real field_2d(nx_l, ny_l, nfields)

        real precip_type_2d(nx_l, ny_l)

        character*9 filename, filename_start

        c character*13 filename13

        character*31 radarext_3d_accum

        c character var*3, comment*125, directory*150, ext*31, units*10
        character var*3, directory*150, ext*31, units*10

        character*125 comment_s, comment_r

        c character*80 c80_domain_file

        integer j_status(20), iprod_number(20)

        character*40 c_vars_req
        character*180 c_values_req

        logical l_reset

        istat = init_timer()

        write (6, *)
1       ' welcome to the laps gridded snow/precip accum analysis'

        c_vars_req = 'radarext_3d_accum'

        call get_static_info(c_vars_req, c_values_req, 1, istatus)

        if (istatus .eq. 1) then
           radarext_3d_accum = c_values_req(1:3)
        else
           write (6, *) ' error getting radarext_3d_accum'
           return
        end if

        write (6, *) ' radarext_3d_accum = ', radarext_3d_accum

        c read in laps lat/lon and topo
        call get_laps_domain(nx_l, ny_l, laps_domain_file, lat, lon, topo
1       , istatus)
        if (istatus .eq. 0) then
           write (6, *) ' error getting laps domain'
           return
        end if

        icen = nx_l/2 + 1
        jcen = ny_l/2 + 1
        call get_grid_spacing_actual_xy(lat(icen, jcen), lon(icen, jcen)
1       , grid_spacing_actual_mx
1       , grid_spacing_actual_my
1       , istatus)
        if (istatus .ne. 1) then
           return
        end if

        grid_spacing_cen_m = grid_spacing_actual_my
        write (6, *) ' actual grid spacing in domain center = '
1       , grid_spacing_cen_m

        ilaps_cycle_time = precip_cycle_time
        write (6, *) ' ilaps_cycle_time = ', ilaps_cycle_time

        n_prods = 1

        do i = 1, n_prods
           j_status(i) = rtsys_abort_prod
        end do

        n_l1s = 1

        ext = 'l1s'

        if (ext(1:3) .ne. 'l1s') then
           iprod_number(n_l1s) = 36352 ! l1s
        else
           iprod_number(n_l1s) = 27280 ! s1s
        end if

!       get incremental accumulations from radar data & update storm totals
        i4time_now = i4time_now_gg()
        i_wait = (i4time + (35*60) - i4time_now)/60
        write (6, *) ' number of potential wait cycles = ', i_wait

        istatus_inc = 1

        i4time_end = i4time
        i4time_beg = i4time_end - ilaps_cycle_time

        minutes = ilaps_cycle_time/60

        write (6, *) ' getting snow/precip accumulation over ', minutes
1       , ' min'
50      call get_precip_inc(i4time_beg, i4time_end, nx_l, ny_l, nz_l   ! i
1       , max_radar_files                                   ! i
1       , lat, lon, topo                                      ! i
1       , ilaps_cycle_time, grid_spacing_cen_m               ! i
1       , radarext_3d_accum                                 ! i
1       , snow_2d, precip_2d, frac_sum                        ! o
1       , istatus_inc)                                      ! o

        if (istatus_inc .ne. 1) then
           write (6, *) ' no incremental precip was generated'
           return
        end if

        comment_r = '           null comment'
        comment_s = '           null comment'

        var = 'sdp'
        call get_laps_2d_prior(i4time, ilaps_cycle_time, ext, var, units
1       , comment_r, nx_l, ny_l, snow_2d_depth, istatus_sdp)
        if (istatus_sdp .eq. 1) then
           write (6, *) ' success reading prior snow depth'
        else
           write (6, *) ' failed reading prior snow depth'
           snow_2d_depth = 0.
        end if

!       get storm total from previous analysis cycle
        if (istatus_inc .eq. 1) then
           write (6, *) ' getting previous storm total accumulations'
           ext = 'l1s'
           var = 'sto'
           call get_laps_2d(i4time - ilaps_cycle_time, ext, var, units
1          , comment_s, nx_l, ny_l, snow_2d_tot, istatus_tot)

           if (istatus_tot .eq. 1) then
              var = 'rto'
              call get_laps_2d(i4time - ilaps_cycle_time, ext, var, units
1             , comment_r, nx_l, ny_l, precip_2d_tot, istatus_tot)
           end if

           if (istatus_tot .eq. 1) then
           end if

           if (istatus_tot .ne. 1) then
              write (6, *) ' warning: could not read previous storm'
1             , ' total accumulations'
           end if
        end if

!       set rate threshold based on length of storm total so far
        rate_thresh_mph = .0001                ! meters per hour
        filename_start = comment_r(1:9)

        if (istatus_tot .eq. 1) then
           call i4time_fname_lp(filename_start, i4time_start_tot
1          , istatus)
        else
           istatus = 0
        end if

        if (istatus .ne. 1) then
           write (6, *) ' could not get start time for storm total'
1          , comment_r
           i4_prev_total = 0
           istatus_tot = 0

        else ! valid storm total start time in comment
           i4_prev_total = i4time_beg - i4time_start_tot
           if (i4_prev_total .gt. 48*3600) then
              rate_thresh_mph = .100
           end if
        end if

        write (6, *) i4_prev_total/3600, ' hours for previous storm total,'
1       , ' rate_thresh_mph = ', rate_thresh_mph

!       decide whether to reset storm total based on insignificant precip over
!       current cycle time
        i_suff_pcp = 0

        precip_reset_thresh =
1       rate_thresh_mph*float(ilaps_cycle_time)/3600.

        if (istatus_inc .eq. 1 .and. istatus_tot .eq. 1) then
           precip_max = 0.
           do j = 1, ny_l
           do i = 1, nx_l
              precip_max = max(precip_max, precip_2d(i, j))
              if (precip_2d(i, j) .gt. precip_reset_thresh) then
                 i_suff_pcp = 1
              end if
           end do ! i
           end do ! j

           write (6, *) ' max precip over domain = ', precip_max
           write (6, *) ' threshold precip over domain = '
1          , precip_reset_thresh
           if (i_suff_pcp .eq. 0) then
              write (6, *) ' resetting due to insufficient precip'
           else
              write (6, *)
1             ' we have sufficient precip to continue storm total'
           end if

        end if

        i4_beg_gmt = i4time_beg - ((i4time_beg/86400)*86400)

!       reset at 1200 gmt
        if (i4_beg_gmt .eq. 12*3600) then
           l_reset = .true.
        else
           l_reset = .false.
        end if

!       if(i_suff_pcp .eq. 1)then ! old method
!           l_reset = .false.
!       else
!           l_reset = .true.
!       endif

        write (6, *) ' l_reset = ', l_reset

!       nominal melt rate of 1mm/hr
        rate_melt = .001/3600.
        snow_melt = rate_melt*float(laps_cycle_time)

        if (istatus_inc .eq. 1) then
!           add current hour precip/snow accumulation to storm total
           if (istatus_tot .eq. 1 .and. (.not. l_reset)) then
              write (6, *) ' adding latest increment for new storm total'
1             , ' accumulation & depth'
              call add_miss(snow_2d, snow_2d_tot, snow_2d_tot
1             , nx_l, ny_l)
              call add_miss(precip_2d, precip_2d_tot, precip_2d_tot
1             , nx_l, ny_l)
              call add_miss(snow_2d, snow_2d_depth, snow_2d_depth
1             , nx_l, ny_l)

              call add_miss(snow_2d, snow_2d_depth, snow_2d_depth
1             , nx_l, ny_l)
              snow_2d_depth = snow_2d_depth - snow_melt

           else ! reset storm total
              write (6, *) ' resetting storm total accumulations to '
1             , minutes, ' min values'

!               put the new reset time in the file header
              call make_fnam_lp(i4time - ilaps_cycle_time, filename
1             , istatus)
              comment_s = filename//
1             ' time that storm total snow begins at.'
              comment_r = filename//
1             ' time that storm total precip begins at.'

              call move(snow_2d, snow_2d_tot, nx_l, ny_l)
              call move(precip_2d, precip_2d_tot, nx_l, ny_l)
              call move(snow_2d, snow_2d_depth, nx_l, ny_l)

           end if ! valid continuation of storm total

!           compare snow depth with analyzed snow cover field
           ext = 'lm2'
           var = 'sc'
           call get_laps_2d(i4time, ext, var, units
1          , comment_s, nx_l, ny_l, snow_cover, istatus_cover)
           if (istatus_cover .eq. 1) then
              write (6, *) ' read in snow cover ok for comparison'
           end if

           do i = 1, nx_l
           do j = 1, ny_l
              if (snow_cover(i, j) .ne. r_missing_data) then
                 if (snow_cover(i, j) .gt. 0.3 .and.
1                snow_2d_depth(i, j) .eq. 0.) then
                 snow_2d_depth(i, j) = 3.
                 elseif (snow_cover(i, j) .eq. 0.3 .and.
1                snow_2d_depth(i, j) .gt. 0.) then
                 snow_2d_depth(i, j) = 0.
              end if
              end if
           end do ! j
           end do ! i

           write (6, *) ' writing incr / storm total accumulations'
           write (6, *) comment_r(1:80)
           write (6, *) comment_s(1:80)
           ext = 'l1s'
           call get_directory(ext, directory, len_dir)
           units = 'm'

           call move(snow_2d, field_2d(1, 1, 1), nx_l, ny_l)
           call move(snow_2d_tot, field_2d(1, 1, 2), nx_l, ny_l)
           call move(precip_2d, field_2d(1, 1, 3), nx_l, ny_l)
           call move(precip_2d_tot, field_2d(1, 1, 4), nx_l, ny_l)
           call move(snow_2d_depth, field_2d(1, 1, 5), nx_l, ny_l)

           call put_precip_2d(i4time, directory, ext, var, units
1          , comment_s, comment_r, nx_l, ny_l, field_2d, ilaps_cycle_time
1          , nfields, istatus)
           if (istatus .eq. 1) then
              j_status(n_l1s) = ss_normal
           else
              write (6, *)
1             ' warning: bad status returned from put_precip_2d'
1             , istatus
           end if

           else
           write (6, *) ' not writing incr / storm total accumulations'
           j_status(n_l1s) = rtsys_no_data

           end if

           i4_elapsed = ishow_timer()

999        continue

           return
        end

        subroutine put_precip_2d(i4time, directory, ext, var, units,
1       comment_s, comment_r, imax, jmax, field_2dsnow
1       , ilaps_cycle_time
1       , nfields, istatus)

        character*(*) directory
        character*31 ext

        character*125 comment_s, comment_r, comment_2d(nfields)
        character*10 units, units_2d(nfields)
        character*3 var, var_2d(nfields)
        integer lvl, lvl_2d(nfields)
        character*4 lvl_coord, lvl_coord_2d(nfields)

        real field_2dsnow(imax, jmax, nfields)

        lend = len(directory)
        write (6, 11) directory(1:lend), ext(1:5)
11      format(' writing 2d snow/precip ', a, 1x, a5, 1x, a3)

        lvl = 0
        lvl_coord = 'msl'

        var_2d(1) = 's01'
        var_2d(2) = 'sto'
        var_2d(3) = 'r01'
        var_2d(4) = 'rto'
        var_2d(5) = 'sdp'

        minutes = ilaps_cycle_time/60

        write (comment_2d(1), 21) minutes
21      format('laps', i3, ' minute snow accumulation')

        comment_2d(2) = comment_s

        write (comment_2d(3), 22) minutes
22      format('laps', i3, ' minute precip accumulation')

        comment_2d(4) = comment_r

        comment_2d(5) = 'laps snow depth'

        do k = 1, nfields
           lvl_2d(k) = lvl
           lvl_coord_2d(k) = lvl_coord
           units_2d(k) = 'm'
           write (6, *) ' range of ', var_2d(k), ' '
1          , minval(field_2dsnow(:, :, k))
1          , maxval(field_2dsnow(:, :, k))
        end do

        call write_laps_data(i4time, directory, ext, imax, jmax,
1       nfields, nfields, var_2d, lvl_2d, lvl_coord_2d, units_2d,
1       comment_2d, field_2dsnow, istatus)
        if (istatus .ne. 1) then
           write (6, *) ' bad status returned from write_laps_data'
        end if

        return
     end

     subroutine get_laps_2d_prior(i4time, ilaps_cycle_time, ext, var
1    , units, comment, nx_l, ny_l, field_2d, istatus)

     character*31 ext
     character*125 comment
     character*10 units
     character*3 var

     real field_2d(nx_l, ny_l)

     do iloop = 1, 50
        call get_laps_2d(i4time - ilaps_cycle_time*iloop, ext, var, units
1       , comment, nx_l, ny_l, field_2d, istatus)
        if (istatus .eq. 1) then
           write (6, *) ' success in get_laps_2d_prior'
1          , iloop, iloop*ilaps_cycle_time
           return
        end if
     end do

     return
  end
