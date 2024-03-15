

        subroutine wind_post_process(i4time_sys
     1                              ,uanl,vanl                            ! i
     1                              ,wanl                                 ! o
     1                              ,nx_l,ny_l,nz_l                       ! i
     1                              ,n_3d_fields                          ! i
     1                              ,heights_3d                           ! i
     1                              ,uanl_sfcitrp,vanl_sfcitrp            ! i
     1                              ,topo,lat,lon,grid_spacing_m          ! i
     1                              ,r_missing_data,l_grid_north_out      ! i
     1                              ,istat_lw3)

        real uanl(nx_l,ny_l,nz_l),vanl(nx_l,ny_l,nz_l) ! wrt true north ! i
        real heights_3d(nx_l,ny_l,nz_l)                                 ! i
        real wanl(nx_l,ny_l,nz_l)                                       ! o
        real uanl_sfcitrp(nx_l,ny_l),vanl_sfcitrp(nx_l,ny_l)            ! o

        real lat(nx_l,ny_l)
        real lon(nx_l,ny_l)
        real topo(nx_l,ny_l)

        real rk_terrain(nx_l,ny_l)

        logical l_grid_north_out

csms$ignore begin
        write(6,*)' subroutine wind_post_process...'

!  **** generate interpolated sfc analysis ****

        write(6,*)' generating interpolated laps surface wind'

        i_sfc_bad = 0

        do j = 1,ny_l
        do i = 1,nx_l

!           interpolate from three dimensional grid to terrain surface
            zlow = height_to_zcoord2(topo(i,j),heights_3d,nx_l,ny_l,nz_l
     1                                                  ,i,j,istatus)
            if(istatus .ne. 1)then
                write(6,*)' lapswind_anal: error in height_to_zcoord2'
     1                   ,' in sfc wind interpolation',istatus
                write(6,*)i,j,zlow,topo(i,j),
     1                    (heights_3d(i,j,k),k=1,nz_l)
                return
            endif

            rk_terrain(i,j) = zlow

            klow = max(zlow,1.)
            khigh = klow + 1
            fraclow = float(khigh) - zlow
            frachigh = 1.0 - fraclow

            if( uanl(i,j,klow)  .eq. r_missing_data
     1     .or. vanl(i,j,klow)  .eq. r_missing_data
     1     .or. uanl(i,j,khigh) .eq. r_missing_data
     1     .or. vanl(i,j,khigh) .eq. r_missing_data        )then

                write(6,3333)i,j
3333            format(' warning: cannot interpolate to sfc at ',2i3)
                i_sfc_bad = 1
                uanl_sfcitrp(i,j) = r_missing_data
                vanl_sfcitrp(i,j) = r_missing_data

            else
                uanl_sfcitrp(i,j) = uanl(i,j,klow ) * fraclow
     1                            + uanl(i,j,khigh) * frachigh

                vanl_sfcitrp(i,j) = vanl(i,j,klow ) * fraclow
     1                            + vanl(i,j,khigh) * frachigh

            endif

        enddo ! j
        enddo ! i

        i4_elapsed = ishow_timer()

        write(6,*)' computing omega'
        call vert_wind(uanl,vanl,uanl_sfcitrp,vanl_sfcitrp                ! i
     1                ,nx_l,ny_l,nz_l                                     ! i
     1                ,wanl                                               ! o
     1                ,topo,lat,lon,grid_spacing_m                        ! i
     1                ,rk_terrain,r_missing_data,l_grid_north_out         ! i
     1                ,istatus)                                           ! o

        if(istatus .ne. 1)then
            write(6,*)' error: bad data detected by vert_wind'
            write(6,*)
     1      ' check for missing data on one or more levels of uanl/vanl'
            return
        endif

        i4_elapsed = ishow_timer()

        return
        end


        subroutine write_wind_output(i4time_sys,ext,var_3d
     1                              ,uanl,vanl                            ! i
     1                              ,wanl                                 ! i
     1                              ,uanl_sfcitrp,vanl_sfcitrp            ! i
     1                              ,nx_l,ny_l,nz_l                       ! i
     1                              ,n_3d_fields                          ! i
     1                              ,r_missing_data                       ! i
     1                              ,istat_lw3)

        use mem_wind, only: num_wind_obs       

!       stuff for 3d winds
        real uanl(nx_l,ny_l,nz_l),vanl(nx_l,ny_l,nz_l) ! wrt true north ! i
        real wanl(nx_l,ny_l,nz_l)                                       ! i

        character*125 comment_3d(n_3d_fields)
        character*10 units_3d(n_3d_fields)
        character*3 var_3d(n_3d_fields)
        character*3 ext

!       stuff for sfc winds
        real uanl_sfcitrp(nx_l,ny_l),vanl_sfcitrp(nx_l,ny_l)            ! i
        real out_sfc_3d(nx_l,ny_l,2)

        character*125 comment_a(2)
        character*10 units_a(2)
        character*3 var_a(2)

        write(6,*)' subroutine write_wind_output...'

!       header information for 3d wind
        ext = 'lw3'

        units_3d(1) = 'm/s'
        units_3d(2) = 'm/s'
        units_3d(3) = 'pa/s'

        var_3d(1) = 'u3'
        var_3d(2) = 'v3'
        var_3d(3) = 'om'

        do i = 1,3
            write(comment_3d(i),1)num_wind_obs
 1          format('3dwind num_wind_obs = ',i9)
        enddo ! i

        i4_elapsed = ishow_timer()

        write(6,*)' calling write routine for all grids ',ext(1:3)
     1                                  ,i4time_sys

        call check_nan3(uanl,nx_l,ny_l,nz_l,istat_lw3)
        if(istat_lw3 .ne. 1)then
            write(6,*)' error: nan detected in uanl field'
            return
        endif

        call check_nan3(vanl,nx_l,ny_l,nz_l,istat_lw3)
        if(istat_lw3 .ne. 1)then
            write(6,*)' error: nan detected in vanl field'
            return
        endif

        call check_nan3(wanl,nx_l,ny_l,nz_l,istat_lw3)
        if(istat_lw3 .ne. 1)then
            write(6,*)' error: nan detected in wanl field'
            return
        endif

        call put_laps_multi_3d_jacket(i4time_sys,ext,var_3d
     1                               ,units_3d,comment_3d
     1                               ,uanl,vanl,wanl
     1                               ,nx_l,ny_l,nz_l,n_3d_fields
     1                               ,istat_lw3)
        if(istat_lw3 .eq. 1)then
            write(6,*)' success writing out lw3 file'
        else
            write(6,*)' error writing out lw3 file'
        endif

!       write out derived winds file (sfc wind)
        call move(uanl_sfcitrp,out_sfc_3d(1,1,1),nx_l,ny_l)
        call move(vanl_sfcitrp,out_sfc_3d(1,1,2),nx_l,ny_l)

        ext = 'lwm'

        var_a(1) = 'su'
        var_a(2) = 'sv'

        do i = 1,2
            units_a(i) = 'm/s'
            comment_a(i) = 'sfcwind'
        enddo

        call put_laps_multi_2d(i4time_sys,ext,var_a
     1      ,units_a,comment_a,out_sfc_3d,nx_l,ny_l,2,istat_lwm)

        if(istat_lwm .eq. 1)then
            write(6,*)' success in writing out lwm file (su,sv)'
        else
            write(6,*)' error writing out lwm file (su,sv)'
        endif

        return
        end



        subroutine put_laps_multi_3d_jacket(i4time_sys,ext,var_3d
     1                                     ,units_3d,comment_3d
     1                                     ,uanl,vanl,wanl
     1                                     ,nx_l,ny_l,nz_l
     1                                     ,n_3d_fields,istat_lw3)

        real uanl(nx_l,ny_l,nz_l),vanl(nx_l,ny_l,nz_l) ! wrt true north
        real wanl(nx_l,ny_l,nz_l)

        real outarray_4d(nx_l,ny_l,nz_l,n_3d_fields)
        character*125 comment_3d(n_3d_fields)
        character*10 units_3d(n_3d_fields)
        character*3 var_3d(n_3d_fields)
        character*3 ext

csms$ignore begin
        write(6,*)' subroutine put_laps_multi_3d_jacket...'

        call move_3d(uanl,outarray_4d(1,1,1,1),nx_l,ny_l,nz_l)
        call move_3d(vanl,outarray_4d(1,1,1,2),nx_l,ny_l,nz_l)
        call move_3d(wanl,outarray_4d(1,1,1,3),nx_l,ny_l,nz_l)

        call put_laps_multi_3d(i4time_sys,ext,var_3d,units_3d,
     1     comment_3d,outarray_4d,nx_l,ny_l,nz_l,n_3d_fields,istat_lw3)

csms$ignore end
        return
        end
