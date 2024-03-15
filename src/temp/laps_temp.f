

        subroutine laps_temp(i4time_needed)

!       1997 jun        ken dritz     made nx_l, ny_l, nz_l dummy arguments,
!                                     making non-dummy arrays dimensioned
!                                     therewith dynamic (automatic).
!       1997 jun        ken dritz     changed include to 
!                                     laps_static_parameters.inc.
!       1997 jun        ken dritz     added call to get_laps_cycle_time.
!       1999 jan        steve albers  added stability stuff

        use mem_namelist, only: nx_l,ny_l,nz_l=>nk_laps
     1                         ,laps_cycle_time,c6_maproj
     1                         ,grid_spacing_m_parm=>grid_spacing_m
     1                         ,iwrite_output

        use mem_grid, only: lat,lon,topo
 
        use mem_temp, only: temp_3d,heights_3d,pres_3d_pa

        include 'laps_static_parameters.inc'

        integer j_status(20),iprod_number(20)

!  ************ declarations **************************************************

        real output_4d(nx_l,ny_l,nz_l,2)

!       real temp_3d(nx_l,ny_l,nz_l)
!       real heights_3d(nx_l,ny_l,nz_l)
!       real pres_3d_pa(nx_l,ny_l,nz_l)
        real pres_3d_mb(nx_l,ny_l,nz_l)
        real temp_sfc_k(nx_l,ny_l)
        real pres_sfc_pa(nx_l,ny_l), pres_sfc_mb(nx_l,ny_l)
        real pres_msl_pa(nx_l,ny_l)
        real pbl_top_pa(nx_l,ny_l), pbl_top_mb(nx_l,ny_l)
        real pbl_depth_m(nx_l,ny_l)

        character*31 ext

        character*10  units_2d
        character*125 comment_2d
        character*3 var_2d

        parameter (max_fields=2)
        real field_array(nx_l,ny_l,max_fields)
        character*125 comment_a(max_fields)
        character*10 units_a(max_fields)
        character*3 var_a(max_fields)

!       obtain grid spacing at the center
!       test for a conformal map projection
        if(c6_maproj .ne. 'latlon' .and. c6_maproj .ne. 'icshdr')then 
            icen = nx_l/2 + 1
            jcen = ny_l/2 + 1
            call get_grid_spacing_actual(lat(icen,jcen),lon(icen,jcen)       
     1                                  ,grid_spacing_m,istatus)
            if(istatus .ne. 1)then
                write(6,*)' error return from get_grid_spacing_actual'       
                return
            endif

        else
            write(6,*)' non-conformal map projection: ',c6_maproj
            write(6,*)' set grid_spacing_cen_m to parameter value'
            grid_spacing_m = grid_spacing_m_parm

        endif

!       read in surface temp data
        var_2d = 't'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_needed,laps_cycle_time/2
     1                      ,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l       
     1                      ,temp_sfc_k,0,istatus)

        if(istatus .ne. 1)then
            write(6,*)' laps sfc temp not available'
            write(6,*)' not calling put_temp_anal'
            go to 999
        endif

!       read in surface pressure data
        var_2d = 'ps'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_needed,laps_cycle_time/2
     1                      ,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,pres_sfc_pa,0,istatus)

        if(istatus .ne. 1)then
            write(6,*)' laps sfc pres not available'
            write(6,*)' not calling put_temp_anal'
            go to 999
        endif

!       read in msl pressure data
        var_2d = 'msl'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_needed,laps_cycle_time/2
     1                      ,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,pres_msl_pa,0,istatus)

        if(istatus .ne. 1)then
            write(6,*)' laps msl pres not available'
            write(6,*)' not calling put_temp_anal'
            go to 999
        endif

!  ************ updated argument list ****************************************

        call put_temp_anal(i4time_needed
     1          ,nx_l,ny_l,nz_l                  ! input
     1          ,heights_3d                      ! output
     1          ,lat,lon,topo                    ! input
     1          ,temp_sfc_k                      ! input
     1          ,pres_sfc_pa                     ! input
     1          ,pres_msl_pa                     ! input
     1          ,laps_cycle_time                 ! input
     1          ,grid_spacing_m                  ! input
     1          ,comment_2d                      ! output
     1          ,temp_3d,pres_3d_pa,istatus)     ! output

        if(iwrite_output .ge. 0 .and. istatus .eq. 1)then
            call write_temp_anal(i4time_needed,nx_l,ny_l,nz_l,temp_3d       
     1                  ,heights_3d,comment_2d,istatus)
        endif

        i4_elapsed = ishow_timer()


!  ******************** pbl section ******************************************

        if(.true. .and. istatus .eq. 1)then
            write(6,*)' start pbl section'
            pres_3d_mb  = pres_3d_pa / 100.
            pres_sfc_mb = pres_sfc_pa / 100.

            call ghbry (i4time_needed,pres_3d_mb,pres_sfc_mb          ! i
     1                 ,temp_sfc_k,temp_3d                            ! i
     1                 ,pbl_top_mb                                    ! o
     1                 ,nx_l,ny_l,nz_l                                ! i
     1                 ,istatus)                                      ! o
            if(istatus .ne. 1)then
                write(6,*)' error: on pbl istatus returned from ghbry'       
                return
            endif

            i4_elapsed = ishow_timer()

            pbl_top_pa = pbl_top_mb * 100.

!           convert to pbl height agl

!           note that the 'pres_to_ht' call uses a linear interpolation that
!           can be upgraded (within the routine) to logp interpolation

!           the 'pressure_to_height' call uses log interpolation but will
!           become an unusable routine if we switch away from a constant
!           pressure vertical grid.

            do i = 1,nx_l
            do j = 1,ny_l
                call pres_to_ht(pbl_top_pa(i,j),pres_3d_pa,heights_3d
     1                         ,nx_l,ny_l,nz_l,i,j,pbl_top_m,istatus)       

!               call pressure_to_height(pbl_top_pa(i,j),heights_3d
!    1                                 ,nx_l,ny_l,nz_l,i,j
!    1                                 ,pbl_top_m,istatus)       

                pbl_depth_m(i,j) = max(pbl_top_m - topo(i,j),0.)
            enddo ! j
            enddo ! i

!           write pbl file
            if(iwrite_output .ge. 0)then
                call move(pbl_top_pa ,field_array(1,1,1),nx_l,ny_l)
                call move(pbl_depth_m,field_array(1,1,2),nx_l,ny_l)

                ext = 'pbl'
                var_a(1) = 'ptp'
                var_a(2) = 'pdm'
                units_a(1) = 'pa'
                units_a(2) = 'm'
                comment_a(1) = 'pbl top pressure'
                comment_a(2) = 'pbl depth'
                call put_laps_multi_2d(i4time_needed,ext,var_a,units_a
     1                                ,comment_a,field_array,nx_l,ny_l
     1                                ,2,istatus)
            endif    

        else
            write(6,*)' no pbl calculation done for pbl file'

        endif

! ************* notification stuff *********************************************

        i4_elapsed = ishow_timer()

        iprod_number(1) = 28261 ! lt1
        n_prods = 1

999     if(istatus .eq. 1)then
            j_status(1) = 1 ! normal product
        else
            j_status(1) = 4 ! no product
        endif

! ****************************************************************************

        return
        end

