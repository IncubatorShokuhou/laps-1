
        subroutine hsect_img(i4time_ref,lun,nx_l,ny_l,nz_l
     1                      ,r_missing_data,laps_cycle_time
     1                      ,plot_parms,namelist_parms,ifield_found)

        use ppm

        real field_2d(nx_l,ny_l)
        integer isky_rgb_cyl(3,nx_l,ny_l)
        character*31 ext
        character*4 var_2d
        character*20 units_2d

        write(6,*)' subroutine hsect_img'

!       read in cloud albedo from lil file
        var_2d = 'cla'
        vel = 0
        ext = 'lil'
        call get_laps_2dgrid(i4time_ref,86400,i4time_cloud,
     1                       ext,var_2d,units_2d,comment_2d,
     1                       nx_l,ny_l,field_2d,0,istatus)

        if(istatus .ne. 1)then
            write(6,*)' skipping image write in hsect_img'
        endif

        do ic = 1,3
            isky_rgb_cyl(ic,:,:) = nint(field_2d(:,:) * 255.)
        enddo ! ic

!       write as png
        write(6,*)' write cloud albedo ppm file '
        call writeppm3matrix(
     1             isky_rgb_cyl(1,:,:),isky_rgb_cyl(2,:,:)
     1            ,isky_rgb_cyl(3,:,:),'cloud_albedo')


        return
        end
