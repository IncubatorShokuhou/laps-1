c
c
        subroutine insert_sat(i4time,cldcv,
     1  cldcv_sao,cld_hts,rlat,rlon,                                    ! i
     1  pct_req_lvd_s8a,default_clear_cover,                            ! i
     1  tb8_k,istat_tb8,                                                ! i
     1  sst_k,istat_sst,                                                ! i
     1  istat_39_a, l_use_39,                                           ! i
     1  di_dh_ir,dj_dh_ir,                                              ! i
     1  di_dh_vis,dj_dh_vis,                                            ! i
     1  i_fill_seams,idb,jdb,                                           ! i
     1  offset_ir_i,offset_ir_j,                                        ! i
     1  istat_39_add_a,                                                 ! o
     1  tb8_cold_k,                                                     ! o
     1  grid_spacing_m,surface_sao_buffer,                              ! i
     1  cloud_frac_vis_a,istat_vis_potl_a,                              ! i
     1  istat_vis_added_a,                                              ! o
     1  solar_alt,solar_ha,solar_dec,                                   ! i
     1  lstat_co2_a, cloud_frac_co2_a, cldtop_co2_pa_a,                 ! i
     1  rlaps_land_frac,                                                ! i
     1  topo,heights_3d,temp_3d,t_sfc_k,td_sfc_k,pres_sfc_pa,           ! i
     1  t_modelfg,sh_modelfg,pres_3d,                                   ! i
     1  cvr_snow,imax,jmax,kcld,klaps,r_missing_data,sfc_albedo,        ! i
     1  t_gnd_k,                                                        ! o
     1  cldtop_co2_m,cldtop_tb8_m,cldtop_m,ht_sao_top,                  ! o
     1  cldht_prlx_top,cldht_prlx_unsm,                                 ! o
     1  istatus)                                                        ! o
c
c*************************************************************************
c
c       routine to process satellite data for clouds and to modify the
c       3d cloud cover array.
c       currently band 8 (11.2 micron) brightness temps are used with a
c       dummy call for the co2 slicing method.
c       cloud building for vis is under construction
c
c
c       draft summary of general goals:
c
c       set minimum cloud top from tb8 (100% cover)
c
c       allow higher cloud top based on vis (thin cloud is higher)
c
c       allow higher cloud top based on rh (via pre-existing cloud, determined 
c       from rh of first guess).
c
c       iterative corrective step of cloud layers vs tb8 is done at the end.
c       this could be moved up earlier? perhaps this and other measures could
c       help avoid inserting new cloud tops too low in the column. a higher
c       top could be used more consistent with rh?
c
c*************************************************************************
c
!       prevents clearing out using satellite (hence letting saos dominate)
!       below this altitude (m agl)
        real surface_sao_buffer

!       prevents adding cloud using satellite (hence letting saos dominate)
!       below this altitude (m agl)
        real surface_ir_buffer
        parameter (surface_ir_buffer = 5000.) ! 3000.

!       default thickness of clouds inserted by satellite
!       real thk_lyr
!       parameter (thk_lyr            = 1500.)

!       cloud cover threshold in evaluating presence of sao/pirep layers
        real thr_sao_cvr
        parameter (thr_sao_cvr = 0.1)

!       threshold for ir cloud detection (sfc temp - ir temp)
        real thresh_ir_diff1
        parameter (thresh_ir_diff1 = 7.)

!       second threshold for ir cloud detection (sfc temp - ir temp)
        real thresh_ir_diff2
        parameter (thresh_ir_diff2 = 21.)

        character*3 lvd_ext
        data lvd_ext /'lvd'/

!       input/output
        real cldcv(imax,jmax,kcld)       ! 3d cloud cover array

        include 'laps_cloud.inc'

!       input
!       real cld_hts(kcloud)
        real cldcv_sao(imax,jmax,kcld)       ! 3d cloud cover array
        real rlat(imax,jmax),rlon(imax,jmax)
        real tb8_k(imax,jmax)
        real di_dh_ir(imax,jmax)           ! parallax offset for sat data
        real dj_dh_ir(imax,jmax)           ! parallax offset for sat data
        real di_dh_vis(imax,jmax)          ! parallax offset for sat data
        real dj_dh_vis(imax,jmax)          ! parallax offset for sat data
        real offset_ir_i(imax,jmax)        ! sat i minus actual i
        real offset_ir_j(imax,jmax)        ! sat j minus actual j
        real sst_k(imax,jmax)
        real tb8_cold_k(imax,jmax)
        real topo(imax,jmax)
        real rlaps_land_frac(imax,jmax)
        real cloud_frac_vis_a(imax,jmax)   ! used for cloud building with vis
        real cloud_frac_vis_s(imax,jmax)   ! smoothed cloud_frac_vis_a
        real solar_alt(imax,jmax)
        real solar_ha(imax,jmax)
        real temp_3d(imax,jmax,klaps)
        real t_modelfg(imax,jmax,klaps)
	real sh_modelfg(imax,jmax,klaps)
        real pres_3d(imax,jmax,klaps)
        real t_sfc_k(imax,jmax)
        real td_sfc_k(imax,jmax)
        real cvr_snow(imax,jmax)
        real sfc_albedo(imax,jmax)
        real pres_sfc_pa(imax,jmax)
        real heights_3d(imax,jmax,klaps)
        real cloud_frac_co2_a(imax,jmax)
        real cldtop_co2_pa_a(imax,jmax)
        real cldht_prlx_top(imax,jmax)
        real cldht_prlx_unsm(imax,jmax)

        integer istat_39_a(imax,jmax)
        integer istat_39_add_a(imax,jmax)
        integer istat_vis_potl_a(imax,jmax)  ! pot'l cloud building with vis
                                             ! (image space)
        integer istat_vis_added_a(imax,jmax) ! actual cloud building with vis
                                             ! (gridpoint space?)
        integer i_fill_seams(imax,jmax)

!       output
        real t_gnd_k(imax,jmax)
        real cldtop_m(imax,jmax)
        real cldtop_co2_m(imax,jmax)
        real cldtop_tb8_m(imax,jmax)
        real ht_sao_top(imax,jmax)

!       local
        real p(1),t(1),td(1),lcl_agl,tlcl_pbe,plcl_pbe ! used for lcl (abdel)
	real lcl_2d(imax,jmax)

        real zcoords_1d(klaps)
        real cldcv_1d(kcloud)
        real laps_p(klaps)

        character*31 ext
        character var*3,comment*125,units*10

        logical  l_tb8,l_cloud_present,l_use_39,l_poss_extrap
     1          ,l_no_sao_vis

        logical lstat_co2_a(imax,jmax)
        logical l_clear_ir /.true./
        logical l_add_ir /.true./
        logical l_correct_cldtop /.true./

        integer idebug_a(imax,jmax)        ! l
        integer k_terrain(imax,jmax)

        real k_to_c

!       control search box for sao analyzed data
        integer idelt(3)
        integer jdelt(3)
!       data idelt/-5,0,+5/
!       data jdelt/0,+5,-5/

        integer nidelt,njdelt
        data nidelt/3/,njdelt/3/

        write(6,*)' subroutine insert_sat...'

!       initialize
        ht_sao_top = r_missing_data
        idebug_tb8 = 0

!       calculate lcl
        do j=1,jmax
        do i=1,imax
            p(1)=pres_sfc_pa(i,j)/100.
            t(1)=t_sfc_k(i,j)-273.15
            td(1)=td_sfc_k(i,j)-273.15

            call lcl_fast(p(1),t(1),td(1),lcl_agl,tlcl_pbe,plcl_pbe)       
            lcl_2d(i,j)=lcl_agl+topo(i,j)
        enddo
        enddo 

        r_missing_ht = 1e30

        idelt_max = nint(50000. / grid_spacing_m)

        idelt(1) = -idelt_max
        idelt(2) = 0
        idelt(3) = +idelt_max

        jdelt(1) = 0
        jdelt(2) = +idelt_max
        jdelt(3) = -idelt_max

        iwrite = 0
        iwrite_lcl = 0

        call filter_2dx_array(cloud_frac_vis_a,r_missing_data,imax,jmax       
     1                       ,cloud_frac_vis_s)              

        do k = 1,klaps
            zcoords_1d(k) = zcoord_of_level(k)
            laps_p(k)     = pressure_of_level(k)
        enddo ! k

        if(istat_tb8 .ne. 1)then
            if(pct_req_lvd_s8a .gt. 0.)then
                write(6,*)
     1              ' error reading tb8_k - return from insert_sat'            
                istatus = 0
                return

            else ! pct_req_lvd_s8a .eq. 0
                write(6,*)' no lvd files, however pct_req_lvd_s8a = 0.'       
                istatus = 1
                return

            endif

        endif

!       calculate valid percentage
        icount = 0
        do j = 1,jmax
        do i = 1,imax
            if(tb8_k(i,j) .eq. r_missing_data)then
                icount = icount + 1
            endif
        enddo
        enddo

        i4_elapsed = ishow_timer()

        valid_pct = (1.- float(icount)/float(imax*jmax) ) * 100.

        if(valid_pct .lt. pct_req_lvd_s8a)then
            write(6,51)pct_req_lvd_s8a,valid_pct
 51         format(' warning: insufficient lvd s8a data ',2f8.2)
            istatus = 0
            return

        else
            write(6,52)pct_req_lvd_s8a,valid_pct
 52         format(' pct lvd s8a data required/actual...',2f8.2)
            
        endif

!       calculate cold filtered temperatures (over ~10km box - mode_sao=2)
        i_delt = max(1,nint(1000./grid_spacing_m)) ! for tb8_cold_k
        do j = 1,jmax
            jl = max(1   ,j-i_delt)
            jh = min(jmax,j+i_delt)

            do i = 1,imax
                il = max(1   ,i-i_delt)
                ih = min(imax,i+i_delt)
                tb8_cold_k(i,j) = tb8_k(i,j)

                do jj = jl,jh,1
                do ii = il,ih,1
                  if(tb8_k(ii,jj) .ne. r_missing_data)then
                    tb8_cold_k(i,j) = min(tb8_cold_k(i,j),tb8_k(ii,jj))
                  endif
                enddo ! ii
                enddo ! jj
            enddo ! i
        enddo ! j


!       calculate ground temperature
        if(istat_sst .ne. 1)then
            write(6,*)' note: sst_k not available'
        endif

        do j = 1,jmax
        do i = 1,imax
          if(rlaps_land_frac(i,j) .ge. 0.5)then
            t_gnd_k(i,j) = t_ground_k(t_sfc_k(i,j),solar_alt(i,j)
     1                               ,solar_ha(i,j),solar_dec,rlat(i,j)       
     1                               ,cvr_snow(i,j),r_missing_data,i,j
     1                               ,imax,jmax)
          else ! water environment
            if(istat_sst .eq. 1 
     1                   .and. sst_k(i,j) .ne. r_missing_data)then     
                t_gnd_k(i,j) = sst_k(i,j)
            else
                t_gnd_k(i,j) = t_sfc_k(i,j)
            endif
          endif
        enddo
        enddo

!       check number of points thrown out as a function of offset to check
!       for ir navigation errors
        thresh1 = 5.
!       call correlation(t_gnd_k,tb8_k,thresh1,imax,jmax)
        call correlation(t_gnd_k,tb8_k,thresh_ir_diff1,imax,jmax)
        thresh3 = 15.
!       call correlation(t_gnd_k,tb8_k,thresh3,imax,jmax)

!       find terrain locations
        do j=1,jmax
        do i=1,imax
             k_terrain(i,j) = 1 ! minimum value allowed
        enddo ! i
        enddo ! j

        do kl=1,klaps
        do j=1,jmax
        do i=1,imax
             if(heights_3d(i,j,kl) .lt. topo(i,j))k_terrain(i,j) = kl
        enddo ! i
        enddo ! j
        enddo ! k

        i4_elapsed = ishow_timer()

        write(6,*)
        write(6,*)'  cldtp_t  i  j  k   frac_k   cldtop_temp_f   ',       
     1            'cldtop_m     sfc t (f)   topo'

        init_co2 = 0 ! initializes co2 slicing routine
        n_valid_co2 = 0
        n_missing_co2 = 0
        n_no_sao1 = 0
        n_no_sao2 = 0
        n_no_sao3 = 0
        n_no_sao_vis = 0

        nskip_max = 4 ! 'see barnes_r5'

        if(imax .lt. 1000)then
            iskd = 20 ! 4
            jskd = 20
        else
            iskd = 50
            jskd = 50
        endif

        write(6,81)
     1  (heights_3d(imax/2,jmax/2,k),cldcv(imax/2,jmax/2,k),k=kcld,1,-1)
81      format(' cldcv section 0 (fg/sao):'                
     1                  /'    ht      cvr',50(/f8.1,f8.3,'   ctr0'))

        mode_prlx = 3 ! 2 (fixed) or 3 (variable)
        cldht_prlx_fixed = 3000.

        if(mode_prlx .eq. 3)then ! extra call to cloud_top to support region growing

          do j=1,jmax
          do i=1,imax

           jp10 = j+10
!          if(j .eq. (j/jskd)*jskd .and. i .eq. (i/iskd)*iskd)then
           if(i .eq. idb .and. j .eq. jdb)then
             idebug_a(i,j) = 1
             write(6,91)i,j,rlat(i,j),rlon(i,j)
91           format(/' debugging at lat/lon ',2i6,2f8.2)
           else
             idebug_a(i,j) = 0
           endif
           idebug = idebug_a(i,j)

           if(tb8_k(i,j) .ne. r_missing_data)then

!           compare brightness temp to surface temperature
            if(idebug .eq. 1)then
              tb8_c = k_to_c(tb8_k(i,j))
              t_gnd_c = k_to_c(t_gnd_k(i,j))
              write(6,111,err=112)i,j,tb8_c,t_gnd_c
     1                               ,tb8_c-t_gnd_c
111           format(1x,2i4,' ctr 1st cloud_top call: tb8/sfc/diff'
     1              ,7x,f8.1,4x,f8.1,8x,f8.1)
112         endif

!           calculate cloud top height from band 8 and/or co2 slicing method
            call cloud_top( init_co2,i4time,tb8_k(i,j),idebug
     1     ,cloud_frac_co2_a(i,j)                                        ! i
     1     ,t_gnd_k(i,j),sfc_albedo,pres_sfc_pa
     1     ,thresh_ir_diff1,topo(i,j),r_missing_data
     1     ,i,j,imax,jmax,klaps,heights_3d,temp_3d,k_terrain(i,j),laps_p      
     1     ,istat_39_a(i,j), l_use_39                                    ! i
     1     ,istat_39_add_a(i,j),istat_vis_added_a(i,j)                   ! o
     1     ,cloud_frac_vis_a(i,j),istat_vis_potl_a(i,j)                  ! i
     1     ,di_dh_vis(i,j),dj_dh_vis(i,j)                                ! i
     1     ,cloud_frac_vis_s(i,j)                                        ! i
     1     ,lstat_co2_a(i,j)                                             ! i
     1     ,t_modelfg,sh_modelfg                                         ! i
     1     ,pres_3d                                                      ! i
     1     ,n_valid_co2,n_missing_co2,cldtop_co2_m(i,j),istat_co2        ! o
     1     ,cldtop_tb8_m(i,j),l_tb8                                      ! o
     1     ,cldtop_m(i,j),l_cloud_present                                ! o
     1     ,sat_cover)                                                   ! o

           endif ! tb8_k

!          this can be set up for mode_prlx
           if(cldtop_tb8_m(i,j) .ne. r_missing_data)then
!          if(cldtop_m(i,j) .ne. r_missing_data)then
              cldht_prlx_top(i,j) = cldtop_tb8_m(i,j)
!             cldht_prlx_top(i,j) = cldtop_m(i,j)
!             cldht_prlx_top(i,j) = cldht_prlx_fixed
           else
!             cldht_prlx_top(i,j) = cldht_prlx_fixed
              cldht_prlx_top(i,j) = r_missing_data
           endif
!          cldht_prlx_top(i,j) =
!    1                     min(max(cldtop_tb8_m(i,j),3000.),3000.)

          enddo ! i
          enddo ! j

!         we might also insert more terrain heights in areas of missing data
          if(.true.)then

           cldht_prlx_unsm(:,:) = cldht_prlx_top(:,:)

           call region_grow(cldht_prlx_top,imax,jmax
     1                     ,r_missing_data,1,30)

           where(cldht_prlx_top(:,:) .eq. r_missing_data)
            cldht_prlx_top(:,:) = cldht_prlx_fixed
           endwhere

           n_cross_in = 41 ! 40km box on a side
           i_l = 1    + n_cross_in/2
           i_h = imax - n_cross_in/2
           j_l = 1    + n_cross_in/2
           j_h = jmax - n_cross_in/2
           call smooth_cross_laps(imax,jmax,i_l,i_h,j_l,j_h
     1                           ,cldht_prlx_top,n_cross_in)

           cldht_prlx_top(:,:) = cldht_prlx_top(:,:) * 0.7 ! experiment

         else ! uniform parallax array
           cldht_prlx_top(:,:) = cldht_prlx_fixed

         endif ! true
             
        endif ! true

!       loop in gridpoint space
        do j=1,jmax 
        do i=1,imax

!        location of satellite cloud tops above grid point
!        testing this earlier calculation of it/jt
         if(mode_prlx .eq. 3)then
             cldht_prlx = cldht_prlx_top(i,j)
         elseif(mode_prlx .eq. 2)then
             cldht_prlx = cldht_prlx_fixed
         else
             cldht_prlx = cld_hts(k)
         endif 
               
         it = i + nint(di_dh_ir(i,j)*cldht_prlx)
         jt = j + nint(dj_dh_ir(i,j)*cldht_prlx)
         it = max(min(it,imax),1)
         jt = max(min(jt,jmax),1)

!        compare to sfc_albedo using parallax and topo?
!        find ig,jg given i,j
         ig = i + nint(di_dh_vis(i,j) * (cldtop_tb8_m(i,j)-topo(i,j)))
         jg = j + nint(dj_dh_vis(i,j) * (cldtop_tb8_m(i,j)-topo(i,j)))
         ig = max(min(ig,imax),1)
         jg = max(min(jg,jmax),1)

         jp10 = j+10
!        if(j .eq. (j/jskd)*jskd .and. i .eq. (i/iskd)*iskd)then
         if(i .eq. idb .and. j .eq. jdb)then
             idebug_a(i,j) = 1
             write(6,91)i,j,rlat(i,j),rlon(i,j)
         else
             idebug_a(i,j) = 0
         endif
         idebug = idebug_a(i,j)

         if(imax-i .le. nskip_max .or. jmax-j .le. nskip_max)then
             l_poss_extrap = .true. ! extrapolation edge effects possible
         else
             l_poss_extrap = .false.
         endif

         if(tb8_k(it,jt) .ne. r_missing_data)then

!         compare brightness temp to surface temperature
          if(idebug .eq. 1)then
              tb8_c = k_to_c(tb8_k(it,jt))
              t_gnd_c = k_to_c(t_gnd_k(ig,jg))
              write(6,1131,err=1132)i,j,tb8_c,t_gnd_c
     1                                 ,tb8_c-t_gnd_c
     1                                 ,istat_vis_potl_a(it,jt)
1131          format(1x,2i4,' ctr 2nd cloud_top call: tb8/sfc/diff/vs'
     1              ,7x,f8.1,4x,f8.1,8x,f8.1,i3)
1132      endif

!         calculate cloud top height from band 8 and/or co2 slicing method
          call cloud_top( init_co2,i4time,tb8_k(it,jt),idebug
     1     ,cloud_frac_co2_a(i,j)                                        ! i
     1     ,t_gnd_k(ig,jg),sfc_albedo,pres_sfc_pa                        ! i
     1     ,thresh_ir_diff1,topo(i,j),r_missing_data                     ! i
     1     ,i,j,imax,jmax,klaps,heights_3d,temp_3d,k_terrain(i,j),laps_p ! i  
     1     ,istat_39_a(i,j), l_use_39                                    ! i
     1     ,istat_39_add_a(i,j),istat_vis_added_a(i,j)                   ! o
     1     ,cloud_frac_vis_a(it,jt),istat_vis_potl_a(it,jt)              ! i
     1     ,di_dh_vis(i,j),dj_dh_vis(i,j)                                ! i
     1     ,cloud_frac_vis_s(it,jt)                                      ! i
     1     ,lstat_co2_a(i,j)                                             ! i
     1     ,t_modelfg,sh_modelfg                                         ! i
     1     ,pres_3d                                                      ! i
     1     ,n_valid_co2,n_missing_co2,cldtop_co2_m(i,j),istat_co2        ! o
     1     ,cldtop_tb8_m(i,j),l_tb8                                      ! o
     1     ,cldtop_m(i,j),l_cloud_present                                ! o
     1     ,sat_cover)                                                   ! o

          if(istat_vis_potl_a(it,jt) .eq. 1)then ! vis potl added
              iwrite = iwrite + 1
              if(idebug .eq. 1)then
                  write(6,113)i,j,istat_vis_added_a(i,j),l_tb8
     1                       ,l_cloud_present
     1                       ,cldtop_m(i,j),cldtop_tb8_m(i,j)
     1                       ,cloud_frac_vis_a(it,jt),sat_cover
     1                       ,tb8_k(it,jt),t_gnd_k(ig,jg)
 113              format(' vis potl added ',2i5,i2,2l2,6e10.3)
              endif
          endif

          if(idebug .eq. 1)write(6,'(" tb8/tgnd check 1 ",2f9.3)')
     1                     tb8_k(it,jt),t_gnd_k(ig,jg)

          if(istat_39_a(i,j) .eq. 1)then ! 3.9u potl added
              iwrite = iwrite + 1
              if((iwrite .lt. 1000 .and. iwrite .eq. (iwrite/10)*10)
     1                              .or. 
     1                    (iwrite .eq. (iwrite/100)*100)
     1                              .or. 
     1                       (.not. l_cloud_present)
     1                              .or. 
     1                       (idebug .eq. 1)
     1                                                          )then
                  write(6,114)i,j,istat_39_add_a(i,j),l_tb8
     1                       ,l_cloud_present
     1                       ,cldtop_m(i,j),cldtop_tb8_m(i,j)
     1                       ,tb8_k(it,jt),t_gnd_k(ig,jg)
 114              format(' 3.9u potl added ',2i5,i2,2l2,4e10.3)
              endif
          endif

          if(idebug .eq. 1)write(6,'(" tb8/tgnd check 2 ",2f9.3)')
     1                     tb8_k(it,jt),t_gnd_k(ig,jg)

          if(lstat_co2_a(i,j))then ! using co2 slicing method

!           clear out those levels higher than what the satellite is showing
!           with the co2 slicing method

            do k=kcld,1,-1
              if(cldcv(i,j,k) .gt. .04)then ! efficiency test
                if(cldtop_co2_m(i,j) .ne. r_missing_data)then
                  if(cld_hts(k) .gt. cldtop_co2_m(i,j))then
                    cldcv(i,j,k) = default_clear_cover ! include surface buffer?
                  endif
                endif
              endif
            enddo ! k

          endif

          if(idebug .eq. 1)write(6,'(" tb8/tgnd check 2a ",2f9.3)')
     1                     tb8_k(it,jt),t_gnd_k(ig,jg)

          if(l_clear_ir)then ! use band 8 (11.2 mm)

!           modify (clear) those levels where the satellite shows warmer than
!           the calculated brightness temp from the analysis of sao/pireps

            do k=kcld,1,-1

!             if(idebug .eq. 1)write(6,'(" tb8/tgnd check 2b ",2f9.3)')
!    1                         tb8_k(it,jt),t_gnd_k(ig,jg)

              if(mode_prlx .eq. 3)then
                  cldht_prlx = cldht_prlx_top(i,j)
              elseif(mode_prlx .eq. 2)then
                  cldht_prlx = cldht_prlx_fixed
              else
                  cldht_prlx = cld_hts(k)
              endif 
                  
              it = i + nint(di_dh_ir(i,j)*cldht_prlx)
              jt = j + nint(dj_dh_ir(i,j)*cldht_prlx)
              it = max(min(it,imax),1)
              jt = max(min(jt,jmax),1)
              if(i_fill_seams(i,j) .ne. 0)then
                  itn = min(i,it)
                  itx = max(i,it)
!                 itn = it
!                 itx = it
              else ! consider gradients
                  if(mode_prlx .eq. 3)then
                    ip = min(i+1,imax)
                    itp = ip - nint(di_dh_ir(ip,j)*cldht_prlx_top(ip,j))
                    itp = min(max(itp,1),imax)
                  else
                    itp = it
                  endif

!                 if(itp-it .ge. 2 .or. .true.)then
                  if(itp-it .ge. 2)then
                    itn = it
                    itx = min(it+1,imax)
                  else
                    itn = it
                    itx = it
                  endif

                  if(mode_prlx .eq. 3)then
                    jp = min(j+1,jmax)
                    jtp = jp - nint(di_dh_ir(i,jp)*cldht_prlx_top(i,jp))
                    jtp = min(max(jtp,1),jmax)
                  else
                    jtp = jt
                  endif

!                 if(jtp-jt .ge. 2 .or. .true.)then
                  if(jtp-jt .ge. 2)then
                    jtn = jt
                    jtx = min(jt+1,jmax)
                  else
                    jtn = jt
                    jtx = jt
                  endif
              endif

!             if(idebug .eq. 1)write(6,'(" tb8/tgnd check 2c ",2f9.3)')
!    1                         tb8_k(it,jt),t_gnd_k(ig,jg)

!             if(cldcv(it,jt,k) .gt. .04)then ! efficiency test
!             if(maxval(cldcv(itn:itx,jtn:jtx,k)) .gt. .04)then ! efficiency test
              if(       cldcv(i,j,k) .gt. .04)then ! efficiency test

                z_temp = height_to_zcoord2(cld_hts(k),heights_3d,
     1                                     imax,jmax,klaps,i,j,istatus1)       

                if(istatus1 .ne. 1)then
                    if(cld_hts(k) .gt. heights_3d(i,j,klaps))then
!                       we cannot perform the cloud clearing when the cloud
!                       height is above the top of the pressure grid
                        go to 500
                    else
                        write(6,*)
     1                   ' error status returned from height_to_zcoord2'       
     1                  ,k,klaps,cld_hts(k),heights_3d(i,j,klaps)
                        istatus = 0
                        return
                    endif
                endif
                z_temp = max(1.,min(z_temp,float(klaps)-.001))
                iz_temp = int(z_temp)
                frac = z_temp - iz_temp
                temp_grid_point = temp_3d(i,j,iz_temp)    * (1. - frac)  
     1                          + temp_3d(i,j,iz_temp+1)  * frac

                if(cldcv(i,j,k) .ne. r_missing_data)then
                  if(cldcv(i,j,k) .gt. 1.070 
     1                    .and. .not. l_poss_extrap)then ! excessively over 1.0
                      write(6,*)
     1                    ' error in insert_sat, interior cldcv > 1.070'     
     1                    ,i,j,k,cldcv(i,j,k)
                      istatus = 0
                      return
                  elseif(cldcv(i,j,k) .gt. 1.0)then ! slightly over 1.0
                      write(6,*)' warning in insert_sat, cldcv > 1.0'       
     1                         ,i,j,k,cldcv(i,j,k)
                      write(6,*)' resetting cloud cover to 1.0'
                      cldcv(i,j,k) = 1.0
                  endif

                  call get_tb8_fwd(temp_grid_point,t_gnd_k(ig,jg)       ! i   
     1                            ,cldcv(i,j,k)                         ! i
     1                            ,tb8_calculated,istatus)              ! o

                  if(istatus .ne. 1)then
                      write(6,*)' error in get_tb8_fwd'
                      write(6,*)'i,j,k',i,j,k
                      istatus = 0
                      return
                  endif

                  tb8_calculated_test = temp_grid_point * cldcv(i,j,k)       
     1                              + t_gnd_k(ig,jg) * (1.-cldcv(i,j,k))
                else
                  tb8_calculated = t_gnd_k(ig,jg)
                endif

!               set thresh_tb8 to +8 in snowy situations where more ir clearing
!               is needed and to +0 in dark land situations where the vis can
!               do more clearing. same is done if no visible data are present
                if(cvr_snow(i,j) .gt. 0. .and. 
     1             cvr_snow(i,j) .ne. r_missing_data)then
!                   thresh_tb8_clr =  0.
!                   tb8_calculated = t_gnd_k(i,j)
                    thresh_tb8_clr  = -8. * (cvr_snow(i,j)**0.3)
                    thresh_tb8_clr2 = -8.          
                elseif(cloud_frac_vis_a(it,jt) .eq. r_missing_data)then ! low sun/night
                    thresh_tb8_clr  = -8.
                    thresh_tb8_clr2 = -8.          
                else ! no snow case
!                   thresh_tb8_clr =  0.    ! low aerosol case?
                    thresh_tb8_clr =  -8.   ! high aerosol case?
                    thresh_tb8_clr2 = 99999.
                endif

!               if(idebug .eq. 1)
!    1               write(6,'(" tb8/tgnd check 2d ",2f9.3)')
!    1                           tb8_k(it,jt),t_gnd_k(ig,jg)

!               test if clouds analyzed by sao/pirep should have been
!               detected by the satellite (sat tb8 may countermand previous analysis)
                iclr = 0
                if(tb8_k(it,jt) - tb8_calculated .gt. thresh_tb8_clr
     1         .or.tb8_k(it,jt) - t_gnd_k(ig,jg) .gt. thresh_tb8_clr2
     1                                                             )then ! sat warmer than calc

!                 only touch points above surface buffer    
                  if(cld_hts(k) - topo(i,j) .gt. surface_sao_buffer)then
!                 if(cld_hts(k) - topo(i,j) .gt. 1200.)then
!                   if(thresh_tb8_clr .gt. 0.)then ! clear clouds when snow present
                    if(.true.)then ! key line for artifacts                                              
                      cldcv(i,j,k)=default_clear_cover
                      iclr = 1
                    else ! .false.
!                     does satellite still imply at least some cloud?
                      if(tb8_k(it,jt) - t_gnd_k(ig,jg)  .lt. -8.0)then ! some cloud
                        if(cldcv(i,j,k) .gt. 0.9)then ! lower top of solid cld
                            cldcv(i,j,k)=default_clear_cover
                            iclr = 2
                        else                          ! cover < 0.9, correct it
                            call get_band8_cover(tb8_k(it,jt)      
     1                                   ,t_gnd_k(ig,jg),temp_grid_point
     1                                   ,idebug,cldcv(i,j,k),istatus)
                            if(cldcv(i,j,k) .gt. 1.0 .or.
     1                         cldcv(i,j,k) .lt. 0.0       )then
                                write(6,*)' error--cover out of bounds'
                                istatus = 0
                                return
                            endif
                        endif
                      else ! band 8 nearly matches ground, clear it
!                       ensure that "black (or grey) stratus" is not present
                        temp_thresh = 
     1                            min(t_gnd_k(ig,jg),t_sfc_k(i,j)-10.0)       
                        if(temp_grid_point .lt. temp_thresh)then
!                       if(temp_grid_point .lt. t_gnd_k(ig,jg))then
                            cldcv(i,j,k)=default_clear_cover ! not in inversion, 
                                                             ! clear it out
                            iclr = 3
                        endif
                      endif ! ir signature present
                    endif ! .true.
                  endif ! analysis has cloud above surface buffer

                endif ! sat is warmer than sao grid point (calc tb8)

!               if(idebug .eq. 1)
!    1               write(6,'(" tb8/tgnd check 2e ",2f9.3)')
!    1                           tb8_k(it,jt),t_gnd_k(ig,jg)

                if(idebug .ge. 1)then
                  write(6,141)i,j,k,t_gnd_k(ig,jg)
     1                       ,tb8_k(it,jt),tb8_calculated
     1                       ,tb8_k(it,jt)-tb8_calculated,thresh_tb8_clr
     1                       ,cvr_snow(i,j),iclr
 141              format(' clr: ijk/gnd/tb8/calc/diff/thr/sncv/iclr = ' 
     1                  ,3i4,6f8.2,i2)
                endif

              endif ! current cloud cover is significant (> .04)

 500        enddo ! k (for clearing clouds)

          endif ! l_clear_ir

          if(idebug .eq. 1)write(6,'(" tb8/tgnd check 3 ",2f9.3)')
     1                     tb8_k(it,jt),t_gnd_k(ig,jg)

!         test if we're confident that a cloud is present and we know where
!         the cloud top is.
          if(l_cloud_present .and. l_add_ir) then ! insert satellite clouds

!           set initial satellite cloud base
!           use lcl as a constraint when cloud base is just below it     
            thk_lyr = cld_thk(cldtop_m(i,j))
            buf_lcl = 0.

            if(cldtop_m(i,j) .gt. (lcl_2d(i,j) + thk_lyr) )then
!               cloud top and base are above lcl
                htbase=cldtop_m(i,j) - thk_lyr
            elseif(cldtop_m(i,j) .gt. (lcl_2d(i,j) + buf_lcl) )then
!               cloud top is above lcl with a buffer
	        htbase=max(lcl_2d(i,j),cldtop_m(i,j) - thk_lyr)
            else 
!               cloud top is near/below lcl - consider relocation later
                htbase=cldtop_m(i,j) - thk_lyr
                iwrite_lcl = iwrite_lcl + 1
                if(iwrite_lcl .le. 100)then
                    write(6,151)i,j,cldtop_m(i,j),lcl_2d(i,j),htbase
 151                format(' note: sat cloud top near or below lcl '
     1                     ,2i6,3f8.1)       
                endif
            endif

!           initialize lowest sao cloud base & highest sao/co2 top
            ht_sao_base = r_missing_ht
            ht_sao_top(i,j) = r_missing_ht

!           locate lowest sao cloud base
            cldcv_below = 0.
            do k=1,kcld
              if(cldcv_below      .le. thr_sao_cvr .and.
     1           cldcv_sao(i,j,k) .gt. thr_sao_cvr)then
                    ht_sao_base = cld_hts(k)
                    goto181
              endif
              cldcv_below = cldcv_sao(i,j,k)
            enddo ! k

 181        i_sao = i
            j_sao = j

!           locate highest sao/co2 top
            cldcv_above = 0.
            do k=kcld,2,-1
              if(cldcv_above      .le. thr_sao_cvr .and.
     1           cldcv_sao(i,j,k) .gt. thr_sao_cvr)then
                    ht_sao_top(i,j) = cld_hts(k)
                    goto186
              endif
              cldcv_above = cldcv_sao(i,j,k)
            enddo ! k

 186        continue

!           if(.false.)then ! search for nearby sao cloud layers
            if(ht_sao_base .eq. r_missing_ht)then 
                                          ! search for nearby sao cloud layers
                                          ! because the satellite says cloud
                                          ! and the sao doesn't
              n_no_sao1 = n_no_sao1 + 1

              do jdelt_index = 1,njdelt
              jj = j + jdelt(jdelt_index)

              do idelt_index = 1,nidelt
              ii = i + idelt(idelt_index)
                if(ii .ge. 1 .and. ii .le. imax .and.
     1             jj .ge. 1 .and. jj .le. jmax        )then ! in bounds

!                   locate lowest neighboring sao cloud base
                    cldcv_below = 0. 
                    do k=1,kcld
                      if(cldcv_below        .le. thr_sao_cvr .and.
     1                   cldcv_sao(ii,jj,k) .gt. thr_sao_cvr)then
                          ht_sao_base = cld_hts(k)
                          goto191
                      endif
                      cldcv_below = cldcv_sao(ii,jj,k)
                    enddo ! k

 191                if(ht_sao_base .ne. r_missing_ht)then
                        i_sao = ii
                        j_sao = jj
                        goto201
                    endif

                endif  ! in bounds
              enddo ! ii
              enddo ! jj
            endif

201         continue

            l_no_sao_vis = .false.
            mode_sao = 0

            if(idebug .eq. 1)write(6,'(" tb8/tgnd check 4 ",2f9.3)')
     1                       tb8_k(it,jt),t_gnd_k(ig,jg)

            if(        ( istat_vis_added_a(i,j) .eq. 1   ! vis/3.9 sat present 
     1              .or. istat_39_add_a(i,j)    .eq. 1 )
     1                                .and.              ! and
     1                 ( ht_sao_base .eq. r_missing_ht   ! no obvious sao base
     1              .or. ht_sao_base .gt. cldtop_m(i,j) )
     1                                                      )then 
              n_no_sao_vis = n_no_sao_vis + 1
              l_no_sao_vis = .true.
              mode_sao = 1

!             calculate/utilize cloud top based on vis cover and tb8 temp
              cover=sat_cover
              htbase = max( topo(i,j), cldtop_m(i,j)-1000. )

              if(cldtop_m(i,j) .eq. 0.           .or. 
     1           abs(cldtop_m(i,j)) .gt. 100000.      )then
                  write(6,*)' warning: cldtop_m(i,j) = ',cldtop_m(i,j)
     1                     ,i,j,mode_sao     
              endif

            elseif(ht_sao_base .eq. r_missing_ht .and.
     1             l_correct_cldtop .eqv. .true.      )then 
                                                ! non-vis sat with no sao cloud
              n_no_sao2 = n_no_sao2 + 1
              mode_sao = 2
              cover=sat_cover
              htbase = ht_sao_base

              if(tb8_k(it,jt) - t_gnd_k(ig,jg)
     1                                    .lt. -thresh_ir_diff2)then 
                  buffer = 2100.             ! we more likely have a cloud
              else                            
                  buffer = surface_ir_buffer ! weed out ir tops w/higher buffer
              endif

!             calculate new cloud top and cover (based on filtered tb8/11u)
!             this gives a cloud edge with more uniform height for an isolated
!             cloud when the edge has a "soft" appearance in the imagery.
!             this helps account for thin clouds by assuming that a thicker
!             cloud nearby will exist with the same cloud top height. this
!             technique might be compared with clustering and a histogram 
!             slope approach using 11u and 6.7u ir.

              cldtop_m_avg = cldtop_m(i,j)

!             compare brightness temp to surface temperature
              if(i .eq. idb .and. j .eq. jdb)idebug = 1
              if(idebug .eq. 1)then
                 write(6,206,err=207)i,j,it,jt,k_to_c(tb8_cold_k(it,jt))        
     1                                   ,k_to_c(t_gnd_k(ig,jg))
206              format(1x,4i4,' ctr 3rd cloud_top call: tb8_cold/sfc'
     1                         ,6x,f8.1,8x,f8.1)
207           endif

              call cloud_top(init_co2,i4time,tb8_cold_k(it,jt),idebug
     1            ,cloud_frac_co2_a(i,j)                                ! i
     1            ,t_gnd_k(ig,jg),sfc_albedo,pres_sfc_pa
     1            ,thresh_ir_diff1,topo(i,j),r_missing_data
     1            ,i,j,imax,jmax,klaps,heights_3d,temp_3d
     1            ,k_terrain(i,j),laps_p
     1            ,istat_39_a(i,j), l_use_39                            ! i
     1            ,istat_39_add_dum,istat_vis_added_dum                 ! o
     1            ,cloud_frac_vis_a(it,jt),istat_vis_potl_a(it,jt)      ! i
     1            ,di_dh_vis(i,j),dj_dh_vis(i,j)                        ! i
     1            ,cloud_frac_vis_s(it,jt)                              ! i
     1            ,lstat_co2_a(i,j)                                     ! i
     1            ,t_modelfg,sh_modelfg                                 ! i
     1            ,pres_3d                                              ! i
     1            ,n_valid_co2,n_missing_co2,cldtop_co2_m(i,j),istat_co2! o
     1            ,cldtop_tb8_m(i,j),l_tb8                              ! o
     1            ,cldtop_m(i,j),l_cloud_present                        ! o
     1            ,sat_cover)                                           ! o

!             calculate the cover (opacity) given the brightness temperature,
!             ground temperature, and assumed ambient cloud-top temperature.
              call get_band8_cover(tb8_k(it,jt),t_gnd_k(ig,jg)
     1                          ,tb8_cold_k(it,jt),idebug,cover,istatus)       
              if(istatus .ne. 1)write(6,*)' bad band8_cover status #1'      
              thk_lyr = cld_thk(cldtop_m(i,j))
              htbase = max( topo(i,j) + buffer , cldtop_m(i,j)-thk_lyr )

              if(htbase .gt. cldtop_m(i,j))then
                  n_no_sao3 = n_no_sao3 + 1
                  mode_sao = 3
              endif

              if(idebug .eq. 1)then
                  write(6,211,err=212)i,j,mode_sao
     1                               ,tb8_k(it,jt),tb8_cold_k(it,jt)
     1                               ,cldtop_m_avg,cldtop_m(i,j)
!       1                            ,cldtop_m_avg,cldtop_m_cold
     1                               ,cover
211               format(1x,2i4,' mode_sao =',i2,': avg/cold'
     1                  ,2f7.1,2x,2f7.0,f9.3)
212               continue
              endif

            elseif(htbase .lt. lcl_2d(i,j) .and. ! .false. .and.      
     1             ( (lcl_2d(i,j)+thk_lyr) .gt. ht_sao_top(i,j) .or.
     1               ht_sao_top(i,j) .eq. r_missing_ht    )
     1                                                        )then        
                                             ! cloud base is below lcl
                                             ! & lcl+thk is higher than ht_sao_top

              if(idebug .eq. 1)write(6,'(" tb8/tgnd check 5 ",2f9.3)')
     1                         tb8_k(it,jt),t_gnd_k(ig,jg)

              mode_sao = 4
              cover=sat_cover
              htbase = lcl_2d(i,j)                     
              thk_lyr = cld_thk(htbase)
              cldtop_old = cldtop_m(i,j)
              cldtop_m(i,j) = htbase + thk_lyr

!             find a thinner value for cloud cover consistent with the new
!             higher cloud top and the known brightness temperature.
!             this works ok if tb8_k is warmer than t at the assumed cloud top
              if(.true.)then ! should this depend on co2?

                  call correct_cover(cover,cover_new,cldtop_old
     1                              ,cldtop_m(i,j)
     1                              ,temp_3d,tb8_k(it,jt),t_gnd_k(ig,jg)
     1                              ,heights_3d
     1                              ,imax,jmax,klaps,i,j,idebug,istatus)       
                  if(istatus .ne. 1)then
                      if(idebug .eq. 1)then
                          write(6,*)' correct_cover: tb8_k < t_cld #4'
                          thk_lyr = cld_thk(cldtop_m(i,j))
                          write(6,*)cldtop_old,cldtop_m(i,j)
     1                             ,htbase,thk_lyr,cover
                          write(6,*)(heights_3d(i,j,k),k=1,klaps)
!                         return
                      endif
                  endif
                  if(idebug .eq. 1)then
                      write(6,261)i,j,it,jt,cover,cover_new
     1                           ,cldtop_old,cldtop_m(i,j)
     1                           ,tb8_k(it,jt),t_gnd_k(ig,jg)
261                   format(' mode_sao = 4: i/j/cvrs/tops/ts '
     1                                                ,4i4,6f8.2)
                  endif
                  cover = cover_new
              endif ! .true.

            elseif(ht_sao_top(i,j) .gt. cldtop_m(i,j) .and. 
     1             ht_sao_top(i,j) .ne. r_missing_ht .and. .true.)then 
                                                    ! satellite top below 
                                                    ! sao top
              mode_sao = 5
              cover=sat_cover
              cldtop_old = cldtop_m(i,j)
              cldtop_m(i,j) = ht_sao_top(i,j)
              htbase = ht_sao_base
              thk_lyr = cld_thk(ht_sao_top(i,j))
              htbase = ht_sao_top(i,j) - thk_lyr
              htbase = max(htbase,lcl_2d(i,j))          

!             find a thinner value for cloud cover consistent with the new
!             higher cloud top and the known brightness temperature.
!             this works ok if tb8_k is warmer than t at the assumed cloud top
              if(.true.)then ! should this depend on co2?

                  call correct_cover(cover,cover_new,cldtop_old
     1                              ,cldtop_m(i,j)
     1                              ,temp_3d,tb8_k(it,jt),t_gnd_k(ig,jg)
     1                              ,heights_3d
     1                              ,imax,jmax,klaps,i,j,idebug,istatus)       
                  if(istatus .ne. 1)then
                      if(idebug .eq. 1)then
                          write(6,*)' correct_cover: tb8_k < t_cld #5'
                          thk_lyr = cld_thk(cldtop_m(i,j))
                          write(6,*)cldtop_old,cldtop_m(i,j)
     1                             ,htbase,thk_lyr,cover
                          write(6,*)(heights_3d(i,j,k),k=1,klaps)
!                         return
                      endif
                  endif
                  cover = cover_new
              endif ! .true.

            else ! normal use of satellite data
              mode_sao = 6
              cover=sat_cover

              if(cldtop_m(i,j) .eq. 0.           .or. 
     1           abs(cldtop_m(i,j)) .gt. 100000.      )then
                  write(6,'(" warning: cldtop_m(i,j) = ",f9.3,2i5,i3)')
     1                                 cldtop_m(i,j),i,j,mode_sao     
              endif

              if(htbase .eq. 0.           .or. 
     1           abs(htbase) .gt. 100000.      )then
                  write(6,*)' warning: htbase = ',htbase,i,j     
              endif

!             locate sao cloud base below satellite cloud top, modify
!             satellite cloud base. highest sao ceiling within default thickness
!             range of satellite layer is used.
              do k=kcld-1,1,-1

                if(       cld_hts(k) .ge. htbase
     1             .and.  cld_hts(k) .le. cldtop_m(i,j)      )then

                  if(cldcv_sao(i_sao,j_sao,k)   .le. thr_sao_cvr .and.
     1               cldcv_sao(i_sao,j_sao,k+1) .gt. thr_sao_cvr
     1                                       )then ! we have an sao base
                    htbase = cld_hts(k+1)

!                   if sao (hence satellite) base is above the satellite
!                   cloud top, lower the satellite base by one grid level
                    if(htbase .gt. cldtop_m(i,j))then
                        htbase = cld_hts(k)
                    endif

                    goto301

                  endif ! we have an sao base

                endif ! in satellite layer

              enddo ! k

            endif ! satellite top below sao ceiling

301         continue

            if(idebug .eq. 1)then
                write(6,302,err=303)i,j,mode_sao,htbase,cover
     1                             ,cldtop_m(i,j),ht_sao_top(i,j)
     1                             ,l_no_sao_vis,istat_vis_potl_a(it,jt)
302             format(1x,2i4,' mode_sao =',i2,f7.0,f7.2,2f7.0,l2,i2)       
303             continue
            endif
            
!           test for questionable cloud layer
            if(htbase .ge. cldtop_m(i,j) .and. mode_sao .ne. 3)then       
                write(6,*)' warning: htbase>cldtp = '
     1                   ,htbase,cover,i,j,mode_sao
     1                   ,cldtop_m(i,j),l_no_sao_vis
     1                   ,istat_vis_potl_a(it,jt)
            endif

            if(htbase .lt. lcl_2d(i,j))then
!               consider relocating cloud above the lcl
                if(idebug .eq. 1)then
                    write(6,305)i,j,htbase,lcl_2d(i,j),mode_sao
     1                         ,cldtop_m(i,j),ht_sao_top(i,j)
 305                format(' warning: sat cloud base below lcl     '
     1                    ,2i6,2f8.1,i4,2f8.1)       
                endif
            endif

!           test for unreasonable cloud layer
            ierr = 0
            if(htbase      .lt. topo(i,j)     .or. 
     1         abs(htbase) .gt. 100000.       .or.
     1         cover       .le. 0.01              )then
                write(6,322,err=323)i,j,mode_sao,htbase,cover
     1                   ,cldtop_m(i,j),ht_sao_top(i,j)
     1                   ,l_no_sao_vis,istat_vis_potl_a(it,jt)
322             format(1x,2i4,' warning: htbase/cover = '
     1                ,i2,f7.0,f7.2,2f7.0,l2,i2)       
                if(cover .lt. 0.)then
                    write(6,*)' error, cover < 0., reset to 0.', cover       
                    cover = 0.
                    ierr = 1
                endif
323             continue
            endif

!           add satellite cloud to array
            if(ierr .eq. 0)then
              do k=kcld,1,-1
                if(mode_prlx .eq. 3)then
                    cldht_prlx = cldht_prlx_top(i,j)
                elseif(mode_prlx .eq. 2)then
                    cldht_prlx = cldht_prlx_fixed
                else
                    cldht_prlx = cld_hts(k)
                endif

                if(cld_hts(k) .ge. htbase  .and.
     1             cld_hts(k) .le. cldtop_m(i,j) )then ! in satellite layer
                   if(istat_vis_added_a(i,j) .eq. 1)then
                       it = i + nint(di_dh_vis(i,j) * cldht_prlx)
                       jt = j + nint(dj_dh_vis(i,j) * cldht_prlx)
                   else
                       it = i + nint(di_dh_ir(i,j) * cldht_prlx)
                       jt = j + nint(dj_dh_ir(i,j) * cldht_prlx)
                   endif
                   it = max(min(it,imax),1)
                   jt = max(min(jt,jmax),1)

                   if(i_fill_seams(i,j) .ne. 0)then
                       itn = min(i,it)
                       itx = max(i,it)
!                      itn = it
!                      itx = it
                   else ! consider gradients
                       if(mode_prlx .eq. 3)then
                         ip = min(i+1,imax)
                         itp = ip -
     1                         nint(di_dh_ir(ip,j)*cldht_prlx_top(ip,j))
                         itp = min(max(itp,1),imax)
                       else
                         itp = it
                       endif

!                      if(itp-it .ge. 2 .or. .true.)then
                       if(itp-it .ge. 2)then
                         itn = it
                         itx = min(it+1,imax)
                       else
                         itn = it
                         itx = it
                       endif

                       if(mode_prlx .eq. 3)then
                         jp = min(j+1,jmax)
                         jtp = jp -
     1                         nint(di_dh_ir(i,jp)*cldht_prlx_top(i,jp))
                         jtp = min(max(jtp,1),jmax)
                       else
                         jtp = jt
                       endif

!                      if(jtp-jt .ge. 2 .or. .true.)then
                       if(jtp-jt .ge. 2)then
                         jtn = jt
                         jtx = min(jt+1,jmax)
                       else
                         jtn = jt
                         jtx = jt
                       endif
                   endif

!                  cldcv(itn:itx,jtn:jtx,k)=cover
                   cldcv(i,j,k)=cover
!                  if(idebug .eq. 1)then
                   if(i .eq. idb .and. j .eq. jdb)then
!                  if(it.eq.660 .and. jt.gt.50 .and. jt.lt.120)then
                       write(6,331)i,j,it,jt,k,cover,cldht_prlx
     1                            ,istat_vis_added_a(i,j)
331                    format(' ctr added ijk/cvr/htprlx:',5i5,2f8.2,i2)       
                   endif
                endif
              enddo
            endif ! ierr = 0 (unreasonable cloud that was below zero cover)

          else
            if(idebug .eq. 1)then
                write(6,*)' cloud not present for adding at i,j = ',i,j
            endif

          endif ! l_cloud_present (cloudy) & l_add_ir

         elseif(idebug .eq. 1)then
            write(6,'("ctr - tb8(it,jt) is missing: ",4i5)')i,j,it,jt
         endif ! tb8_k(it,jt) .ne. r_missing_data

        enddo ! imax
        enddo ! jmax

        write(6,332)
     1  (heights_3d(idb,jdb,k),cldcv(idb,jdb,k),k
     1                                           ,k=kcld,1,-1)
332     format(' cldcv section 1 (after subtract/add):'
     1             /'    ht      cvr    k',50(/f8.1,f8.3,i4,'   ctr1'))

!       write stats on co2 and band 8 (11.2mm) methods
        write(6,*)' n_valid_co2 = '  ,n_valid_co2
     1           ,' n_missing_co2 = ',n_missing_co2
        write(6,*)' n_no_sao (1/2/3/vis) = '
     1             ,n_no_sao1,n_no_sao2,n_no_sao3,n_no_sao_vis

        i4_elapsed = ishow_timer()

        call compare_radiation(kcld,temp_3d,klaps,imax,jmax
     1       ,cldcv,cldcv_1d,cld_hts,t_sfc_k,t_gnd_k,tb8_k
     1       ,idebug_a,r_missing_data,cvr_snow,heights_3d,nlyr,istatus)
        if(istatus .ne. 1)then
            write(6,*)
     1           ' warning: bad status returned from compare_radiation'       
        endif

!       visible addition stats
        n_vis_add_potl = 0
        n_vis_added = 0
        do i = 1,imax
        do j = 1,jmax
            n_vis_add_potl = n_vis_add_potl + istat_vis_potl_a(it,jt)
            n_vis_added    = n_vis_added    + istat_vis_added_a(i,j)
        enddo ! j
        enddo ! i
        write(6,*)' n_vis_add_potl = ',n_vis_add_potl
     1           ,' n_vis_added = ',n_vis_added

        write(6,*)'section pt vis potl/added'
     1        ,istat_vis_potl_a(idb,jdb)
     1        ,istat_vis_added_a(idb,jdb)

        write(6,342)
     1  (heights_3d(idb,jdb,k),cldcv(idb,jdb,k),k=kcld,1,-1)
342     format(' cldcv section 1a (comp rad):'
     1                  /'    ht      cvr',50(/f8.1,f8.3))

        istatus = 1
        return
        end

        subroutine cloud_top( init_co2,i4time,tb8_k,idebug             ! i
     1  ,cloud_frac_co2                                                ! i
     1  ,t_gnd_k,sfc_albedo,pres_sfc_pa,thresh_ir_diff1,topo           ! i
     1  ,r_missing_data                                                ! i
     1  ,i,j,imax,jmax,klaps,heights_3d,temp_3d,k_terrain,laps_p       ! i
     1  ,istat_39, l_use_39                                            ! i
     1  ,istat_39_add,istat_vis_added                                  ! o
     1  ,cloud_frac_vis_a,istat_vis_potl,di_dh_vis,dj_dh_vis           ! i
     1  ,cloud_frac_vis_s                                              ! i
     1  ,lstat_co2                                                     ! i
     1  ,t_modelfg,sh_modelfg                                          ! i
     1  ,pres_3d                                                       ! i
     1  ,n_valid_co2,n_missing_co2,cldtop_co2_m,istat_co2              ! o
     1  ,cldtop_tb8_m,l_tb8                                            ! o
     1  ,cldtop_m,l_cloud_present                                      ! o
     1  ,sat_cover)                                                    ! o


!       this routine computes the cloud top height given a band 8 brightness
!       temperature and 3d fields of temp and height. the co2 method is also
!       employed to yield a cloud top pressure and height. this is not yet
!       fully integrated into the final cloud analysis.

        real zeros
        parameter (zeros = 1.e-30)

!       argument list
        integer init_co2                      ! input
        integer i4time                        ! input
        real tb8_k                            ! input
        real cloud_frac_vis_a                 ! input (vis cloud building)
        real cloud_frac_vis_s                 ! input (vis cloud building)
        integer i,j,imax,jmax,klaps           ! input
        real t_gnd_k                          ! input
        real sfc_albedo(imax,jmax)            ! input
        real pres_sfc_pa(imax,jmax)           ! input
        real t_modelfg(imax,jmax,klaps)       ! input
        real sh_modelfg(imax,jmax,klaps)      ! input
        real pres_3d(imax,jmax,klaps)         ! input
        real thresh_ir_diff1                  ! input
        real topo                             ! input
        real r_missing_data                   ! input
        integer istat_vis_potl                ! input (vis cloud building)
                                              !       (sat image space)
        real heights_3d(imax,jmax,klaps)      ! input
        real temp_3d(imax,jmax,klaps)         ! input
        integer k_terrain                     ! input
        real laps_p(klaps)                    ! input
        integer n_valid_co2,n_missing_co2     ! input/output
        real cldtop_co2_m                     ! output
        logical lstat_co2,l_use_39            ! input
        integer istat_co2                     ! output (not presently used)
        real cldtop_tb8_m                     ! output
        logical l_tb8                         ! output
        real cldtop_m                         ! output
        logical l_cloud_present               ! output
        real sat_cover                        ! output

!       local
!       real dum_3d(imax,jmax,klaps)          ! local (dummy array for q)
        real arg,frac_k,temp_above,cldtop_temp_k
        integer kl

!       function calls
        real k_to_c, jcost_cldtop

!       used for adding visible
        visthr = 0.1 ! 0.0 (experiment)

!       call the co2 slicing method to get cloud tops
        if(init_co2 .eq. 0)istatus_co2 = 0
        istatus_co2 = -3

!       convert co2 cloud top from pressure to height
        if(istatus_co2 .eq. 0)then
            istat_co2 = 1
            cldtop_co2_pa = iipcld * 100.
            cldtop_co2_z = zcoord_of_logpressure(cldtop_co2_pa)

            iarg = int(cldtop_co2_z)
            frac = cldtop_co2_z - float(iarg)
            if(iarg .lt. klaps)then
                cldtop_co2_m = heights_3d(i,j,iarg) * (1.0 - frac) +
     1                 heights_3d(i,j,iarg+1) * frac
            else
                cldtop_co2_m = r_missing_data
                write(6,*)' warning, co2 pressure is out of bounds '
     1                   ,iipcld
            endif

        else
            istat_co2 = 0
            cldtop_co2_m = r_missing_data

        endif

        if(cldtop_co2_m .ne. r_missing_data)then
            n_valid_co2 = n_valid_co2 + 1
        else
            n_missing_co2 = n_missing_co2 + 1
            cldtop_co2_m = r_missing_data 
        endif

!       this section finds the cloud top using band 8 data and temperatures
!       estimate whether tb8_k - t < threshold
        cldtop_tb8_m = r_missing_data ! zeros
        if(tb8_k - t_gnd_k .lt. -thresh_ir_diff1) then ! probably clds
            l_tb8 = .true.

        else ! no clouds according to satellite (band 8 - 11.2mm)
            l_tb8 = .false.

        endif

        cloud_frac_tb8 = 1.0

        if(               l_tb8 
     1                     .or. 
     1       (l_use_39 .and. istat_39 .eq. 1)
     1                     .or. 
     1       (istat_vis_potl .eq. 1 .and. 
     1        cloud_frac_vis_a .ne. r_missing_data)       
     1                                            )then ! get 11u cloud top
            cldtop_temp_k_before = tb8_k

!           correct cloud top temperature for thin clouds using vis data
            call correct_cldtop_t_rad(tb8_k,t_gnd_k                    ! i
     1                           ,cloud_frac_vis_s                     ! i
     1                           ,istat_vis_potl,idebug                ! i
     1                           ,cldtop_temp_k,istat_vis_corr)        ! o
            if(istatus .eq. -1)then
                write(6,*)
     1      '-1 status returned from correct_cldtop_t_rad at i,j = ',i,j
            endif

!           locate cloud top in 3-d temperature grid (using lowest crossing point)

            costmin_x = 9999. ! initial cost value from crossing points

	    if (cldtop_temp_k .le. temp_3d(i,j,klaps))then ! restrict to top
                cldtop_temp_k = temp_3d(i,j,klaps)
            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! original version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            temp_above = temp_3d(i,j,klaps)
            do kl = klaps-1,k_terrain,-1

                if( (temp_3d(i,j,kl) - cldtop_temp_k) *
     1              (temp_above      - cldtop_temp_k) .lt. 0.)then ! crossing pt

                    frac_k = (cldtop_temp_k - temp_3d(i,j,kl))
     1                    /  (temp_above    - temp_3d(i,j,kl))

                    arg = heights_3d(i,j,kl) + frac_k *
     1                   (heights_3d(i,j,kl+1) - heights_3d(i,j,kl))

                    if(arg .ge. topo)then
                        cldtop_tb8_m = arg
                    endif

                    qamb = sh_modelfg(i,j,kl) + frac_k *
     1                    (sh_modelfg(i,j,kl+1) - sh_modelfg(i,j,kl))

                    p_pa = pres_3d(i,j,kl) + frac_k *
     1                    (pres_3d(i,j,kl+1) - pres_3d(i,j,kl))

                    costmin_x = jcost_cldtop(cldtop_temp_k,tb8_k,qamb
     1                                                 ,p_pa,idebug)      

                    if(idebug .eq. 1)then
                        if(arg .lt. topo)then
                            write(6,*)
     1                      ' cloud top below ground - not used'
                        endif

                        write(6,111,err=121)cldtop_tb8_m,cldtop_co2_m
111                     format(1x,'cloud_top tb8/co2:'f10.0,1x,f10.0)

121                     write(6,122,err=123)i,j,kl,frac_k
     1                           ,k_to_c(cldtop_temp_k_before)
     1                           ,k_to_c(cldtop_temp_k)
     1                           ,arg,topo,k_to_c(temp_3d(i,j,kl))
     1                           ,k_to_c(temp_above)
     1                           ,costmin_x
122                     format(1x,2i4,' cldtp_t',i4,f8.3,2f8.1,f11.1
     1                           ,f21.1,2f6.1,f9.2)
123                 endif

                endif

                temp_above = temp_3d(i,j,kl)

            enddo ! kl

!           evaluate cost function at all the laps levels
            costmin_lvl = 9999.       

            if(idebug .eq. 1)then
                write(6,*)
     1      '          p      tb8     del_t   del_q  del_rh     jcost'
            endif

            do kl = klaps,k_terrain,-1
                costlvl = jcost_cldtop(temp_3d(i,j,kl),tb8_k
     1                                ,sh_modelfg(i,j,kl)
     1                                ,pres_3d(i,j,kl),idebug)

                if(costlvl .lt. costmin_lvl)then
                    costmin_lvl = costlvl
                    kl_min = kl
                endif
            enddo ! kl

            if(costmin_lvl .lt. costmin_x)then

                if(idebug .eq. 1)then
                    write(6,*)i,j,
     1            'cloud_top: found lower cost function, call corr-cvr '
                endif

!               calculate pot'l new cloud top and cover                             
                cldtop_new_potl_m = heights_3d(i,j,kl_min)

                if(temp_3d(i,j,kl_min) .lt. tb8_k)then
                    call correct_cover(cloud_frac_tb8
     1                                ,cloud_frac_tb8_potl
     1                                ,cldtop_tb8_m,cldtop_new_potl_m          
     1                                ,temp_3d,tb8_k,t_gnd_k
     1                                ,heights_3d,imax,jmax
     1                                ,klaps,i,j,idebug,istatus)
                else
                    cloud_frac_tb8_potl = cloud_frac_tb8
                endif

                if((cloud_frac_tb8_potl .ge. 0. .or.
     1              cloud_frac_tb8_potl .le. 1.) .and.
     1                                             istatus .eq. 1)then 
!                   set new cloud top and cover
                    cldtop_tb8_m = cldtop_new_potl_m            
                    cloud_frac_tb8 = cloud_frac_tb8_potl
                elseif(idebug .eq. 1)then
                    write(6,131)cloud_frac_tb8_potl,istatus
131                 format(' warning in cloud_top: cloud_frac_tb8_potl '
     1                    ,'out of bounds or istatus=0 ',f9.3,i2)
                    idebug = 1
                endif

                if(idebug .eq. 1)then
                    write(6,132)i,j,kl_min                        
     1                     ,costmin_lvl,tb8_k,temp_3d(i,j,kl_min)
     1                     ,cloud_frac_tb8_potl
     1                     ,cldtop_new_potl_m              
132                 format(' cloud_top: found lower cost function at: '
     1                    ,2i6,i4,f12.4,2f8.1,f6.2,f8.0)
                endif

            else ! no lower cost function
                if(idebug .eq. 1)then
                    write(6,*)' cloud_top: no lower cost function '          
     1                        ,costmin_lvl,costmin_x                 
                endif
            endif

        endif ! we will want to use a sat/model determined cloud top

        istat_39_add = 0
        istat_vis_added = 0
        mode_top = 0
        ig = i
        jg = j

!       set variables depending on whether in band 8 or co2 mode
        if(lstat_co2)then ! using co2 method
            mode_top = 1
            if(cldtop_co2_m .ne. r_missing_data)then
                l_cloud_present = .true.
            else
                l_cloud_present = .false.
            endif

            cldtop_m = cldtop_co2_m
            sat_cover = cloud_frac_co2

        elseif( (.not. l_tb8) .and. (l_use_39 .and. istat_39 .eq. 1) 
     1                        .and. cldtop_tb8_m .ne. r_missing_data 
     1                                                            ) then      

!           band 8 (11mm) threshold says no but 3.9 micron says yes
!           we did still get a valid band 8 derived cloud top

            mode_top = 2
            l_cloud_present = .true.
            cldtop_m = cldtop_tb8_m
            sat_cover = cloud_frac_tb8
            istat_39_add = 1

        elseif( (.not. l_tb8) .and. istat_vis_potl .eq. 1 
     1                .and. cloud_frac_vis_a .ne. r_missing_data
     1                .and. cldtop_tb8_m     .ne. r_missing_data
     1                                                            ) then      

!           band 8 (11mm) threshold says no but visible says yes
!           we did still get a valid band 8 derived cloud top

!           compare to sfc_albedo using parallax and topo?
!           find ig,jg given i,j
            ig = i - nint(di_dh_vis * (cldtop_tb8_m-topo))
            jg = j - nint(dj_dh_vis * (cldtop_tb8_m-topo))
            ig = max(min(ig,imax),1)
            jg = max(min(jg,jmax),1)

!           test only works if 'sfc_albedo' is set to non-missing value
!           if(cloud_frac_vis_a .gt. sfc_albedo(ig,jg) + visthr)then       
            if(cloud_frac_vis_a .gt.                     visthr)then       
                mode_top = 3
                l_cloud_present = .true.
                cldtop_m = cldtop_tb8_m
                sat_cover = cloud_frac_vis_a
                istat_vis_added = 1
                if(idebug .eq. 1)then
                    write(6,*)i,j,' visible is being added',ig,jg
                endif
            else
                mode_top = 4
                l_cloud_present = .false.
                sat_cover = 0.
            endif

        elseif( l_tb8 .and. istat_vis_potl .eq. 1 
     1                .and. cloud_frac_vis_a .ne. r_missing_data
     1                .and. cldtop_tb8_m     .ne. r_missing_data
     1                                                            ) then      
            l_cloud_present = l_tb8
            cldtop_m = cldtop_tb8_m
            if(cloud_frac_vis_a .gt. visthr)then       
                mode_top = 5              ! possibly use vis?
                sat_cover = cloud_frac_vis_a
                istat_vis_added = 1
                if(idebug .eq. 1)then
                    write(6,*)i,j,' visible is being added',ig,jg
                endif
            else
                mode_top = 6              ! possibly use vis?
                sat_cover = cloud_frac_tb8
            endif

        else                          ! using band 8 (11.2mm) data only
            mode_top = 7
            l_cloud_present = l_tb8
            cldtop_m = cldtop_tb8_m
            sat_cover = cloud_frac_tb8

        endif

        ierr = 0
        if(l_cloud_present)then
          if(cldtop_m .eq. r_missing_data .or.
     1       sat_cover .lt. 0.            .or.
     1       sat_cover .gt. 1.                 )then
            write(6,'(" error in cloud_top: i,j,cldtp_m,sat_cvr",2i5
     1                         ,f9.3,f9.4)')i,j,cldtop_m,sat_cover
            ierr = 1
            if(sat_cover .gt. 1.)then
                sat_cover = 1.
            endif
            if(sat_cover .gt. 0.)then
                sat_cover = 0.
            endif
          endif
        endif

        if(idebug .eq. 1 .or. ierr .eq. 1)then
            write(6,201)i,j,l_tb8,l_cloud_present,mode_top
     1             ,istat_vis_potl,istat_vis_added,istat_39_add
     1             ,cldtop_m,cldtop_tb8_m
     1             ,tb8_k,t_gnd_k,cloud_frac_vis_a
     1             ,cloud_frac_vis_s,sfc_albedo(ig,jg)
     1             ,cldtop_temp_k,sat_cover,cloud_frac_tb8
 201        format(' cloud_top info: ',2i5,2l2,4i2,2f9.0,2f9.2      
     1               ,' vis ',2f7.4,2f9.2,2f7.3)
        endif

        return
        end


        subroutine correlation(t,tb8_k,thresh,ni,nj)

        integer ibox
        parameter (ibox = 3)

        real t(ni,nj),tb8_k(ni,nj)
        integer i_corr_array(-ibox:ibox,-ibox:ibox)

        write(6,*)
        write(6,*)' checking navigation of ir satellite, thresh = '
     1                                          ,thresh

        thresh_minus = -thresh
        min_count = 999999

        do joff = -ibox,ibox
        do ioff = -ibox,ibox
            icount = 0

            do j = 1,nj
              jj = max(min(j+joff,nj),1)

              do i = 1,ni

                ii = max(min(i+ioff,ni),1)

                if(tb8_k(ii,jj) - t(i,j) .lt. thresh_minus)then
                    icount = icount+1 ! clouds indicated
                endif

              enddo ! i

            enddo ! j

            if(icount .lt. min_count)then
                min_count = icount
                ioff_min = ioff
                joff_min = joff
            endif

            i_corr_array(ioff,joff) = icount

        enddo ! ioff
        enddo ! joff

        do j = -ibox,ibox
            write(6,101)(i_corr_array(i,j),i=-ibox,ibox)
101         format(10i8)
        enddo

        write(6,201)min_count,ioff_min,joff_min
201     format('  minimum of',i6,' points flagged with an offset of',2i4
     1)

        return
        end

        subroutine correct_cover(cover_in,cover_new_f,cldtop_old
     1             ,cldtop_new,temp_3d,tb8_k,t_gnd_k,heights_3d
     1             ,imax,jmax,klaps,i,j,idebug,istatus)

!       find a thinner value for cloud cover consistent with the new
!       higher cloud top and the known brightness temperature.
!       the cover in this routine should represent visible opacity as this
!       is most in line with ir emissivity

        real temp_3d(imax,jmax,klaps)
        real heights_3d(imax,jmax,klaps)

        integer iwrite
        data iwrite /0/
        save iwrite

        logical lwrite

        istatus = 1

!       find temperature of old cloud top
        z_temp = height_to_zcoord2(cldtop_old,heights_3d,imax,jmax,klaps
     1                                          ,i,j,istatus)
        if(istatus .ne. 1)then
            return
        endif

        z_temp = max(1.,min(z_temp,float(klaps)-.001))
        iz_temp = int(z_temp)
        frac = z_temp - iz_temp
        temp_old = temp_3d(i,j,iz_temp)    * (1. - frac)
     1          +  temp_3d(i,j,iz_temp+1)  * frac

!       find temperature of new cloud top
        z_temp = height_to_zcoord2(cldtop_new,heights_3d,imax,jmax,klaps
     1                                          ,i,j,istatus)
        if(istatus .ne. 1)then
            cldtop_new = cldtop_old
            cover_new_f = cover_in
            return
        endif

        z_temp = max(1.,min(z_temp,float(klaps)-.001))
        iz_temp = int(z_temp)
        frac = z_temp - iz_temp
        temp_new = temp_3d(i,j,iz_temp)    * (1. - frac)
     1          +  temp_3d(i,j,iz_temp+1)  * frac

!       this one utilizes a linear approximation to the sigma t**4 relationship
        cover_old = min(cover_in,1.0)
        cover_new = cover_old *
     1  (tb8_k - t_gnd_k) / (temp_new - t_gnd_k)

        if(idebug .eq. 1)then
            lwrite = .true.
        else
            lwrite = .false.
        endif

!       this one utilizes the sigma t**4 relationship
        call get_band8_cover(tb8_k,t_gnd_k,temp_new,idebug
     1                      ,cover_new_f,istatus)
        if(istatus .ne. 1)then
            if(lwrite)then
                write(6,*)' bad band8_cover status #2 ',cover_new_f      
            endif
        endif

        if(cover_new_f .lt. -0.0001)then
            if(lwrite)then
              write(6,*)' error in correct_cover: corrected cover << 0.'      
     1                  ,i,j,cover_new_f
            endif
            istatus = 0
            cldtop_new = cldtop_old
            cover_new_f = cover_in
        endif

        if(cover_new_f .gt. 1.0)then
            if(lwrite)then
              write(6,*)' error in correct_cover: corrected cover > 1. '       
     1                  ,i,j,cover_new_f
            endif
            istatus = 0
            cldtop_new = cldtop_old
            cover_new_f = cover_in
        endif

        if(lwrite)then
            write(6,1,err=2)i,j,tb8_k,t_gnd_k,temp_old,temp_new
     1             ,cldtop_old,cldtop_new,cover_in,cover_new,cover_new_f
1           format(1x,'corr-cvr t/top/cvr',2i5,4f7.1,2f8.0,3f8.2)
2           continue
        endif

        return
        end


        subroutine get_band8_cover(tb8_k,t_gnd_k,t_cld,idebug   ! i
     1                            ,band8_cover,istatus)         ! o

        r_sfc = temp_to_rad(t_gnd_k)
        r_sat = temp_to_rad(tb8_k)
        r_cld = temp_to_rad(t_cld)

        band8_cover = (r_sat - r_sfc) / (r_cld - r_sfc)

        istatus = 1

        if(band8_cover .gt. 1.0)then
            if(idebug .eq. 1)then
              write(6,*)
     1             ' warning: get_band8_cover resetting down to 1.0'      
              write(6,11)tb8_k,t_gnd_k,t_cld,band8_cover
 11           format(' tb8_k,t_gnd_k,t_cld,band8_cover:',4f9.3,f8.4)
            endif
            band8_cover = 1.0 
            istatus = 0
        elseif(band8_cover .lt. 0.0)then
            if(idebug .eq. 1)then
              write(6,*)' warning: get_band8_cover resetting up to 0.0'       
              write(6,11)tb8_k,t_gnd_k,t_cld,band8_cover
            endif
            band8_cover = 0.0 
            istatus = 0
        endif

        return
        end

        subroutine compare_radiation(kcld,temp_3d,klaps,imax,jmax
     1      ,cldcv,cldcv_1d,cld_hts,t_sfc_k,t_gnd_k,tb8_k
     1      ,idebug_a,r_missing_data,cvr_snow,heights_3d,nlyr,istatus)

        include 'laps_cloud.inc'

        real temp_3d(imax,jmax,klaps)
        real heights_3d(imax,jmax,klaps)
        real cvr_snow(imax,jmax)
        real tb8_k(imax,jmax),t_gnd_k(imax,jmax),t_sfc_k(imax,jmax)
        real a(100),f(100)
        integer ilyr(kcloud) ! dimension needs to be changed to kcloud
        real a_new(100),f_new(100)
        integer ilyr_new(kcloud) ! dimension needs to be changed to kcloud
        real cldcv(imax,jmax,kcld)
        real cldcv_1d(kcld) ! ,cld_hts(kcld)
        logical l_correct,l_output
        integer idebug_a(imax,jmax)        

!       this routine compares the analyzed clouds to the 11.2mm radiation
!       and determines adjusted cloud fractions of the cloud layers to yield
!       a better fit.

        write(6,*)' comparing radiation'

        istatus=1

        corr_thr = 4.
        cover_step = 0.05
        iter_max = 10

        iwrite = 0
        n_clear = 0
        n_clear_ns = 0
        n_clear_sn = 0
        n_cloudy = 0
        n_total = 0
        n_correct = 0
        n_iter = 0
        tdiff_sumsq = 0.
        tdiff_cld_sumsq = 0.
        tdiff_corr_sumsq = 0.
        tdiff_corr_sumsq2 = 0.
        tdiff_corr_sumsq3 = 0.
        tdiff_corr_cld_sumsq = 0.
        tb8_g_clr_sum   = 0.
        tb8_a_clr_sum   = 0.
        tb8_g_clr_sumsq = 0.
        tb8_a_clr_sumsq = 0.
        tb8_g_clr_ns_sum   = 0.
        tb8_a_clr_ns_sum   = 0.
        tb8_g_clr_ns_sumsq = 0.
        tb8_a_clr_ns_sumsq = 0.
        tb8_g_clr_sn_sum   = 0.
        tb8_a_clr_sn_sum   = 0.
        tb8_g_clr_sn_sumsq = 0.
        tb8_a_clr_sn_sumsq = 0.

        idebug_tb8 = 0

        do j = 1,jmax
        do i = 1,imax

          if(idebug_a(i,j) .eq. 1)then
            l_output = .true.
          else
            l_output = .false.
          endif

          if(tb8_k(i,j) .ne. r_missing_data)then
            do k = 1,kcld
                cldcv_1d(k) = cldcv(i,j,k)
            enddo

            call cvr_to_tb8_effective(kcld,temp_3d,klaps,i,j,imax,jmax
     1      ,a,f,ilyr,cldcv_1d,cld_hts,t_gnd_k(i,j),heights_3d
     1      ,t_effective,nlyr,idebug_tb8,istatus)
            tdiff = tb8_k(i,j)-t_effective ! band 8 - calculated
            frac_clouds = 1.-f(nlyr)

            if(nlyr .eq. 1)then
                tb8_g_clr_sum   = tb8_g_clr_sum   + tdiff
                tb8_g_clr_sumsq = tb8_g_clr_sumsq + tdiff**2
                tb8_a_clr_sum   = tb8_a_clr_sum   + tb8_k(i,j) 
     1                                            - t_sfc_k(i,j)
                tb8_a_clr_sumsq = tb8_a_clr_sumsq +
     1                                  (tb8_k(i,j) - t_sfc_k(i,j))**2
                n_clear = n_clear + 1

                if(cvr_snow(i,j) .ne. r_missing_data)then
                    if(cvr_snow(i,j) .lt. 0.5)then
                        tb8_g_clr_ns_sum   = tb8_g_clr_ns_sum   + tdiff
                        tb8_g_clr_ns_sumsq = tb8_g_clr_ns_sumsq + tdiff*
     1*2
                        tb8_a_clr_ns_sum   = tb8_a_clr_ns_sum
     1                             + tb8_k(i,j) - t_sfc_k(i,j)
                        tb8_a_clr_ns_sumsq = tb8_a_clr_ns_sumsq +
     1                              (tb8_k(i,j) - t_sfc_k(i,j))**2
                        n_clear_ns = n_clear_ns + 1
                    else
                        tb8_g_clr_sn_sum   = tb8_g_clr_sn_sum   + tdiff
                        tb8_g_clr_sn_sumsq = tb8_g_clr_sn_sumsq + tdiff*
     1*2
                        tb8_a_clr_sn_sum   = tb8_a_clr_sn_sum
     1                             + tb8_k(i,j) - t_sfc_k(i,j)
                        tb8_a_clr_sn_sumsq = tb8_a_clr_sn_sumsq +
     1                              (tb8_k(i,j) - t_sfc_k(i,j))**2
                        n_clear_sn = n_clear_sn + 1
                    endif
                endif
            else
                n_cloudy = n_cloudy + 1
                tdiff_cld_sumsq = tdiff_cld_sumsq + tdiff**2

            endif

            n_total = n_total + 1
            tdiff_sumsq = tdiff_sumsq + tdiff**2

!           apply corrections?
            if(nlyr .ge. 2 .and. abs(tdiff) .gt. corr_thr
     1                   .and. t_gnd_k(i,j) - tb8_k(i,j) .gt. 8.
     1             .and. frac_clouds .gt. 0.4             )then
                l_correct = .true.
                n_correct = n_correct + 1
                tdiff_corr_sumsq3 = tdiff_corr_sumsq3 + tdiff**2
            else
                l_correct = .false.
            endif

            if(abs(tdiff) .gt. corr_thr .and. l_output)then
                write(6,901,err=902)i,j,nlyr-1,tb8_k(i,j),t_effective
     1               ,tdiff,frac_clouds,l_correct,(a(l),l=nlyr-1,1,-1)       
901             format(1x,2i4,i3,2f7.0,f7.1,f7.3,l2,' *',10f6.2)
902         endif

!           this corrective section is turned off since it is not yet parallax
!           corrected
            iter = 0
            delta_cover = 0.

905         if(l_correct .and. .false.)then
                iter = iter + 1
                tdiff_ref = tdiff
                delta_cover_ref = delta_cover

                if(tdiff .lt. 0.)then
                    delta_cover = delta_cover + cover_step
                else
                    delta_cover = delta_cover - cover_step
                endif

                call apply_correction(cldcv,i,j,imax,jmax
     1                       ,cldcv_1d,1,1,1,1,kcld,f,ilyr,delta_cover)

                call cvr_to_tb8_effective(kcld,temp_3d,klaps,i,j
     1              ,imax,jmax,a_new,f_new,ilyr_new,cldcv_1d,cld_hts       
     1              ,t_gnd_k(i,j),heights_3d
     1              ,t_effective,nlyr_new,idebug_tb8,istatus)

                tdiff = tb8_k(i,j)-t_effective
                frac_clouds = 1.-f_new(nlyr_new)

!               continue to apply corrections?
                if(nlyr_new .ge. 2 .and. frac_clouds .gt. 0.4)then
                    i_correct = 1
                else
                    i_correct = 0
                endif

                if(iwrite .le. 100 .and. l_output)then
                    iwrite = iwrite + 1
                    write(6,911,err=912)i,j,nlyr_new-1,tb8_k(i,j)
     1               ,t_effective
     1               ,tdiff,frac_clouds,i_correct
     1               ,(a_new(l),l=nlyr_new-1,1,-1)
911                 format(1x,2i4,i3,2f7.0,f7.1,f7.3,i2,'  ',10f6.2)
912             endif

                if(i_correct .eq. 1 .and. iter .lt. iter_max
     1          .and. abs(tdiff) .lt. abs(tdiff_ref)
     1          .and.     tdiff * tdiff_ref .gt. 0.
     1                          )goto905 ! loop back & increment cover

                n_iter = n_iter + iter

!               final iteration
                if(.false. .or. tdiff*tdiff_ref .ge. 0.)then
                  ! select best of last two increments
                    if(abs(tdiff) .ge. abs(tdiff_ref))then
                        delta_cover = delta_cover_ref
                        tdiff = tdiff_ref
                    endif

                else           ! do one newton iteration
                    frac = tdiff / (tdiff - tdiff_ref)
                    delta_cover = delta_cover
     1                  + frac * (delta_cover_ref - delta_cover)

                    call apply_correction(cldcv,i,j,imax,jmax
     1                        ,cldcv_1d,1,1,1,1,kcld,f,ilyr,delta_cover)

                    call cvr_to_tb8_effective(kcld,temp_3d,klaps,i,j
     1                                       ,imax,jmax
     1                                       ,a_new,f_new,ilyr_new
     1                                       ,cldcv_1d,cld_hts
     1                                       ,t_gnd_k(i,j),heights_3d       
     1                                       ,t_effective,nlyr_new
     1                                       ,idebug_tb8,istatus)
                    tdiff = tb8_k(i,j)-t_effective
                    n_iter = n_iter + 1
                endif

!               if(iwrite .le. 100 .and. l_output)write(6,913)delta_cover,tdiff
                if(iter .gt. 0 .and. l_output)write(6,913)delta_cover
     1                                                   ,tdiff
913             format(70x,'best delta_cover/tdiff = ',f6.2,f6.1)

!               apply correction to 3d cloud cover field
                call apply_correction(cldcv,i,j,imax,jmax
     1                   ,cldcv,i,j,imax,jmax,kcld,f,ilyr,delta_cover)

            endif ! corrections were made

            tdiff_corr_sumsq     = tdiff_corr_sumsq + tdiff**2
            if(nlyr .gt. 1)then ! cloudy (at least before adjustment)
                tdiff_corr_cld_sumsq = tdiff_corr_cld_sumsq + tdiff**2
            endif

            if(l_correct)then
                tdiff_corr_sumsq2 = tdiff_corr_sumsq2 + tdiff**2
            endif

          endif ! tb8_k(i,j) is valid
        enddo ! i
        enddo ! j

!       write out statistics on consistency of band 8 data and cloud cover
        if(n_clear .gt. 0)then
            write(6,951,err=960)n_clear, tb8_g_clr_sum  /float(n_clear)
     1                 ,sqrt    (tb8_g_clr_sumsq/float(n_clear))
951         format(/
     1      ' mean/rms band 8/gnd  temp residual in clear skies =  '
     1                                                ,i5,2f9.3)
960         write(6,961,err=970)n_clear, tb8_a_clr_sum  /float(n_clear)
     1                 ,sqrt    (tb8_a_clr_sumsq/float(n_clear))
961         format(' mean/rms band 8/air  temp residual in clear skies =
     1  '
     1                                                ,i5,2f9.3)
970     endif

        if(n_clear_ns .gt. 0)then
            write(6,956,err=965)n_clear_ns
     1                         ,tb8_g_clr_ns_sum  /float(n_clear_ns)
     1                       ,sqrt(tb8_g_clr_ns_sumsq/float(n_clear_ns))
956         format(/
     1      ' mean/rms band 8/gnd temp resid in clear/nsnow skies ='
     1                                                ,i5,2f9.3)
965         write(6,966,err=975)n_clear_ns
     1                         ,tb8_a_clr_ns_sum  /float(n_clear_ns)
     1                       ,sqrt(tb8_a_clr_ns_sumsq/float(n_clear_ns))
966         format(
     1      ' mean/rms band 8/air temp resid in clear/nsnow skies ='
     1                                                ,i5,2f9.3)
975     endif

        if(n_clear_sn .gt. 0)then
            write(6,976,err=985)n_clear_sn
     1                         ,tb8_g_clr_sn_sum  /float(n_clear_sn)
     1                       ,sqrt(tb8_g_clr_sn_sumsq/float(n_clear_sn))
976         format(/
     1      ' mean/rms band 8/gnd temp resid in clear/snow skies = '
     1                                                ,i5,2f9.3)
985         write(6,986,err=995)n_clear_sn
     1                         ,tb8_a_clr_sn_sum  /float(n_clear_sn)
     1                       ,sqrt(tb8_a_clr_sn_sumsq/float(n_clear_sn))
986         format(
     1      ' mean/rms band 8/air temp resid in clear/snow skies = '
     1                                                ,i5,2f9.3)
995     endif

        if(n_total .gt. 0)then
            write(6,971,err=980)n_total,sqrt(tdiff_sumsq/float(n_total))
971         format(/
     1      ' rms  band 8/teff residual (bfr corr) in all  skies = '
     1                                                ,i5, f9.3)
980         write(6,981,err=990)n_total
     1                         ,sqrt(tdiff_corr_sumsq/float(n_total))
981         format(
     1      ' rms  band 8/teff residual (aft corr) in all  skies = '
     1                                                ,i5, f9.3)
990     endif

        if(n_cloudy .gt. 0)then
            write(6,991,err=1000)n_cloudy
     1                          ,sqrt(tdiff_cld_sumsq/float(n_cloudy))
991         format(
     1      ' rms  band 8/teff residual (bfr corr) in cldy skies = '
     1                                                ,i5, f9.3)
1000        write(6,1001,err=1010)
     1              n_cloudy,sqrt(tdiff_corr_cld_sumsq/float(n_cloudy))
1001        format(
     1      ' rms  band 8/teff residual (aft corr) in cldy skies = '
     1                                                ,i5, f9.3)
1010    endif

        if(n_correct .gt. 0)then
            write(6,1011,err=1020)
     1          n_correct,sqrt(tdiff_corr_sumsq3/float(n_correct))
1011        format(/
     1      ' rms  band 8/teff resid (bfr corr - corrected pts)  = '
     1                                                ,i5, f9.3)
1020        write(6,1021,err=1030)
     1          n_correct,sqrt(tdiff_corr_sumsq2/float(n_correct))
1021        format(
     1      ' rms  band 8/teff resid (aft corr - corrected pts)  = '
     1                                                ,i5, f9.3)
            write(6,*)' total/average # of iterations = ',n_iter,
     1                          float(n_iter)/float(n_correct)
            write(6,*)
1030    endif

        return
        end


        subroutine apply_correction(cldcv_in,i_in,j_in,imax_in,jmax_in
     1                             ,cldcv_out,i_out,j_out,imax_out
     1                             ,jmax_out,kcld,f,ilyr,delta_cover)

        real cldcv_in(imax_in,jmax_in,kcld)
        real cldcv_out(imax_out,jmax_out,kcld)
        real f(100)
        integer ilyr(kcld)    ! dimension needs to be changed to kcloud

!       apply correction to 3d cloud cover field
        do k = 1,kcld
            if(ilyr(k) .gt. 0)then
                if(f(ilyr(k)) .gt. 0.)then
                    cldcv_out(i_out,j_out,k) =
     1          min(cldcv_in(i_in,j_in,k),1.0) + delta_cover
                    cldcv_out(i_out,j_out,k) =
     1          max(min(cldcv_out(i_out,j_out,k),1.0),0.0)
                endif
            endif
        enddo

        return
        end

        subroutine correct_cldtop_t(tb8_k,t_gnd_k,cloud_frac_vis         ! i
     1                             ,istat_vis_potl                       ! i
     1                             ,cldtop_temp_k,istatus)               ! o

        cldtop_temp_k = tb8_k ! default value
        istatus = 0

!       correct the cloud top temperature for thin clouds using vis data
        if(cloud_frac_vis .ge. 0. .and. cloud_frac_vis .le. 1.
     1                            .and. istat_vis_potl .eq. 1)then
            diff = tb8_k - t_gnd_k
            if(diff .lt. 0.)then            
                diff2 = diff * 1./cloud_frac_vis
                cldtop_temp_k = t_gnd_k + diff2
                istatus = 1                
            endif
        endif

        return
        end

        subroutine correct_cldtop_t_rad(tb8_k,t_gnd_k,cloud_frac_vis      ! i
     1                                 ,istat_vis_potl,idebug             ! i
     1                                 ,cldtop_temp_k,istatus)            ! o

        use cloud_rad

        data iwrite /0/
        save iwrite

        cldtop_temp_k = tb8_k ! default value
        istatus = 0

!       correct the cloud top temperature for thin clouds using vis data
        if(cloud_frac_vis .ge. 0. .and. cloud_frac_vis .le. 1.
     1                            .and. istat_vis_potl .eq. 1)then
            tb8_rad = temp_to_rad(tb8_k)
            t_gnd_rad = temp_to_rad(t_gnd_k)

            diff = tb8_rad - t_gnd_rad
            if(diff .lt. 0.)then            
                call albedo_to_clouds2(cloud_frac_vis
     1                                ,cloud_trans_l,cloud_trans_i
     1                                ,cloud_od_l,cloud_od_i
     1                                ,cloud_opac_l,cloud_opac_i)
                diff2 = diff * 1./cloud_opac_i
                cldtop_temp_rad = t_gnd_rad + diff2

                if (cldtop_temp_rad .le. 0.0) then
                   if(iwrite .le. 100 .or. idebug .eq. 1)then
                       write(6,*) 'warning: cldtop_temp_rad <= 0'
                       write(6,'(" tb8_k   ",f9.3,"   tb8_rad  ",f9.3)')
     1                             tb8_k,             tb8_rad
                       write(6,'(" t_gnd_k ",f9.3,"   t_gnd_rad",f9.3)')
     1                             t_gnd_k,           t_gnd_rad
                       write(6,*) 'cloud_frac_vis ',cloud_frac_vis
                       write(6,*) 'cloud_opac_i ',cloud_opac_i
                       iwrite = iwrite + 1
                   endif
                   istatus = -1
                   return
                endif

                cldtop_temp_k = rad_to_temp(cldtop_temp_rad)
                istatus = 1                
            endif
        endif

        return
        end

        subroutine get_tb8_fwd(t_ambient_k,t_gnd_k,cldcv               ! i
     1                        ,tb8_fwd_k,istatus)                      ! o

!       forward model for single cloud layer

        tb8_calculated_rad = 
     1                temp_to_rad(t_ambient_k)     * cldcv +
     1                temp_to_rad(t_gnd_k)         * (1.-cldcv)       

        if(tb8_calculated_rad .lt. 0.)then
            write(6,*)' error, tb8_calculated_rad < 0 '
     1                ,tb8_calculated_rad
            write(6,*)'cldcv',cldcv
            write(6,*)'t_ambient_k,t_gnd_k'
     1                ,t_ambient_k,t_gnd_k
            istatus = 0
            return
        endif

        tb8_fwd_k = rad_to_temp(tb8_calculated_rad)

        istatus = 1

        return
        end

        subroutine filter_2dx_array(array_in,r_missing_data,ni,nj
     1                             ,array_out)              

        real array_in(ni,nj)
        real array_out(ni,nj)

        do i = 1,ni
        do j = 1,nj       
            if(array_in(i,j) .eq. r_missing_data)then
                array_out(i,j) = 0.
            else
                array_out(i,j) = array_in(i,j)
            endif
        enddo ! j
        enddo ! i

        call filter_2dx(array_out,ni,nj,1,+0.5)
        call filter_2dx(array_out,ni,nj,1,-0.5)

        return
        end

        function jcost_cldtop(tamb_k,tb8_k,qamb,p_pa,idebug)

        real jcost_cldtop, k_to_c, make_ssh

        delta_t = tb8_k - tamb_k

        p_mb = p_pa / 100.

        tamb_c = k_to_c(tamb_k)

!       qsat_amb = ssh(p_mb,tamb_c) / 1000.
        qsat_amb = make_ssh(p_mb,tamb_c,1.0,-5.) / 1000.

        delta_q = max(qsat_amb - qamb,0.) ! set supersaturated to 0.

        rh_amb = qamb / qsat_amb ! w.r.t. liquid or ice

        delta_rh = max(1.0 - rh_amb,0.) ! set supersaturated to 0.

!       favor warm brightness temp with still some allowance for cold
!       brightness temps in the boundary layer
        if(delta_t .gt. 0.)then
            t_coeff = 0.04
        else
            t_coeff = 2.
        endif

        jcost_cldtop = abs(delta_t) * t_coeff
     1               + abs(delta_q) * 100.
     1               + abs(delta_rh) * 10.

        if(idebug .eq. 1)then
            write(6,1)p_mb,tb8_k,delta_t,delta_q*1000.,delta_rh*100.
     1               ,jcost_cldtop
1           format(1x,'jcost',2f8.1,f8.1,f8.1,f8.1,f12.4)
        endif

        return
        end
