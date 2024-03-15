        use mem_namelist, only: laps_cycle_time, precip_cycle_time

!       1997 jun        ken dritz     added calls to get_grid_dim_xy and
!                                     get_laps_dimensions to get values of
!                                     nx_l, ny_l, nz_l.
!       1997 jun        ken dritz     now pass nx_l, ny_l, nz_l to laps_accum.

        c integer j_status(20), iprod_number(20), i4time_array(20)
        integer j_status(20), iprod_number(20)
        character*9 a9_time

        call get_systime(i4time, a9_time, istatus)
        if (istatus .ne. 1) go to 999

        call get_grid_dim_xy(nx_l, ny_l, istatus)
        if (istatus .ne. 1) then
           write (6, *) 'error getting horizontal domain dimensions'
           go to 999
        end if

        call get_laps_dimensions(nz_l, istatus)
        if (istatus .ne. 1) then
           write (6, *) 'error getting vertical domain dimension'
           go to 999
        end if

        call get_max_radar_files(max_radar_files, istatus)
        if (istatus .ne. 1) then
           write (6, *) 'error getting max_radar_files'
           go to 999
        end if

        if (i4time .eq. (i4time/precip_cycle_time)*
1       precip_cycle_time) then

        call laps_accum_sub(i4time,
1       nx_l, ny_l, nz_l, max_radar_files,
1       i_diag,
1       n_prods,
1       iprod_number,
1       j_status)
        else
        write (6, '(" skipping call to laps_accum_sub")')
     end if

999  continue

  end

