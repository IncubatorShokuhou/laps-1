
        subroutine calc_allsky(i4time_solar,exposure ! ,clwc_3d,cice_3d
!    1                     ,heights_3d                              ! i
!    1                     ,rain_3d,snow_3d                         ! i
!    1                     ,pres_3d,aod_3d                          ! i
     1                     ,topo_sfc,topo,swi_2d,pw_2d              ! i
     1                     ,topo_albedo_2d,land_frac,snow_cover     ! i
     1                     ,htagl                                   ! i
     1                     ,aod_ref                                 ! i
     1                     ,nx_l,ny_l,nz_l,newloc                   ! i
     1                     ,ri_obs,rj_obs                           ! i
     1                     ,alt_a_roll,azi_a_roll                   ! i
     1                     ,sol_alt_2d,sol_azi_2d                   ! i
     1                     ,solar_alt,solar_az                      ! i
     1                     ,solar_lat,solar_lon,r_au                ! i
     1                     ,alt_norm                                ! i
     1                     ,moon_alt_2d,moon_azi_2d,alm,azm         ! i
     1                     ,moon_mag,moon_mag_thr,elgms             ! i
     1                     ,l_solar_eclipse,eobsc,emag              ! i
     1                     ,rlat,rlon,lat,lon                       ! i
     1                     ,minalt,maxalt,minazi,maxazi,nsp         ! i
     1                     ,ni_cyl,nj_cyl,elong_a                   ! o
     1                     ,alt_scale,azi_scale                     ! i
     1                     ,grid_spacing_m,r_missing_data           ! i
     1                     ,l_binary,l_terrain_following            ! i
     1                     ,mode_cloud_mask,camera_cloud_mask       ! i
     1                     ,cloud_od,dist_2_topo                    ! o
     1                     ,sky_rgb_cyl,istatus)                    ! o

        use mem_allsky
        use mem_namelist, only: earth_radius, ssa
        use mem_namelist, only: fcterm, angstrom_exp_a

        include 'trigd.inc'
        include 'wa.inc'

        real ext_o(nc)               ! od per airmass
        data ext_o /.037,.029,.001/  ! interp from gflash.bas (300du)

        parameter (pi = 3.14159265)
        parameter (rpd = pi / 180.)

!       statement functions
        trans(od) = exp(-min(od,80.))
        opac(od) = 1.0 - trans(od)

        rayleigh_pf(theta) =
     1      (1.00 * ( (1. + cosd(theta)**2) / (4./3.) ) ) + (0.00 * 1.0) 

        addlogs(x,y) = log10(10.**x + 10.**y)
        angdist(p1,p2,dlon) = acosd(sind(p1) * sind(p2)
     1                      + cosd(p1) * cosd(p2) * cosd(dlon))

!       input arrays
!       real clwc_3d(nx_l,ny_l,nz_l)      ! control variable
!       real cice_3d(nx_l,ny_l,nz_l)      ! control variable
!       real heights_3d(nx_l,ny_l,nz_l)   ! control variable
!       real rain_3d(nx_l,ny_l,nz_l)      ! control variable
!       real snow_3d(nx_l,ny_l,nz_l)      ! control variable
!       real pres_3d(nx_l,ny_l,nz_l)      ! control variable
!       real aod_3d(nx_l,ny_l,nz_l)       
        real topo(nx_l,ny_l)
        real land_frac(nx_l,ny_l)
        real snow_cover(nx_l,ny_l)
        real du_2d(nx_l,ny_l)
        real aod_2d(nx_l,ny_l)
        real swi_2d(nx_l,ny_l)
        real pw_2d(nx_l,ny_l)
        real topo_albedo_2d(nc,nx_l,ny_l)
        real alt_a_roll(minalt:maxalt,minazi:maxazi)
        real azi_a_roll(minalt:maxalt,minazi:maxazi)
        real sol_alt_2d(nx_l,ny_l)
        real sol_azi_2d(nx_l,ny_l)
        real alt_norm(nx_l,ny_l)   ! solar alt w.r.t. terrain normal
        real eobsc(nx_l,ny_l)
        real moon_alt_2d(nx_l,ny_l)
        real moon_azi_2d(nx_l,ny_l)
        real lat(nx_l,ny_l)
        real lon(nx_l,ny_l)
        real elong_a(minalt:maxalt,minazi:maxazi)
        integer camera_cloud_mask(minalt:maxalt,minazi:maxazi)
        integer sim_cloud_mask(minalt:maxalt,minazi:maxazi)
        integer diff_cloud_mask(minalt:maxalt,minazi:maxazi)

!       output arrays
        real sky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi) ! observed variable
        real sky_sprad(0:2,minalt:maxalt,minazi:maxazi) 
        real sky_reflectance(0:2,minalt:maxalt,minazi:maxazi) 
        real cloud_od(minalt:maxalt,minazi:maxazi)
        real*8 dist_2_topo(minalt:maxalt,minazi:maxazi)

!       local arrays (e.g. outputs from get_cloud_rays)
!       real transm_3d(nx_l,ny_l,nz_l)
!       real transm_4d(nx_l,ny_l,nz_l,nc) 

        real r_cloud_3d(minalt:maxalt,minazi:maxazi)
        real cloud_od_sp(minalt:maxalt,minazi:maxazi,nsp)
        real cloud_od_sp_w(minalt:maxalt,minazi:maxazi,nsp)
        real airmass_2_cloud_3d(minalt:maxalt,minazi:maxazi)
        real airmass_2_topo_3d(minalt:maxalt,minazi:maxazi)
        real topo_swi(minalt:maxalt,minazi:maxazi)
        real topo_albedo(nc,minalt:maxalt,minazi:maxazi)
        real topo_lf(minalt:maxalt,minazi:maxazi)
        real topo_sc(minalt:maxalt,minazi:maxazi)
        real du_a(minalt:maxalt,minazi:maxazi)
        real aod_a(minalt:maxalt,minazi:maxazi)
        real topo_ri(minalt:maxalt,minazi:maxazi)
        real topo_rj(minalt:maxalt,minazi:maxazi)
        real*8 topo_lat(minalt:maxalt,minazi:maxazi)
        real*8 topo_lon(minalt:maxalt,minazi:maxazi)
        real topo_lat_r4(minalt:maxalt,minazi:maxazi)
        real topo_lon_r4(minalt:maxalt,minazi:maxazi)
        real trace_ri(minalt:maxalt,minazi:maxazi)
        real trace_rj(minalt:maxalt,minazi:maxazi)
        real topo_solalt(minalt:maxalt,minazi:maxazi)
        real topo_solazi(minalt:maxalt,minazi:maxazi)
        real trace_solalt(minalt:maxalt,minazi:maxazi)
        real eobsc_sky(minalt:maxalt,minazi:maxazi)
        real gtic(nc,minalt:maxalt,minazi:maxazi)
        real dtic(nc,minalt:maxalt,minazi:maxazi)
        real btic(nc,minalt:maxalt,minazi:maxazi)
        real emic(nc,minalt:maxalt,minazi:maxazi)
        real aod_2_cloud(minalt:maxalt,minazi:maxazi)
        real aod_2_topo(minalt:maxalt,minazi:maxazi)
        real aod_ill(minalt:maxalt,minazi:maxazi)
        real aod_ill_dir(minalt:maxalt,minazi:maxazi)
        real aod_tot(minalt:maxalt,minazi:maxazi)
        real r_cloud_rad(minalt:maxalt,minazi:maxazi)
        real cloud_rad_c(nc,minalt:maxalt,minazi:maxazi)
        real cloud_rad_w(minalt:maxalt,minazi:maxazi)
        real cloud_sfc_c(nc,minalt:maxalt,minazi:maxazi)  
        real cloud_rad_c_nt(nc,minalt:maxalt,minazi:maxazi)
        real clear_rad_c(nc,minalt:maxalt,minazi:maxazi)
        real clear_radf_c(nc,minalt:maxalt,minazi:maxazi)
        real clear_rad_c_nt(nc,minalt:maxalt,minazi:maxazi)

!       other local arrays
        real blog_moon_roll(minalt:maxalt,minazi:maxazi)
        real blog_sun_roll(minalt:maxalt,minazi:maxazi)
        real glow_stars(nc,minalt:maxalt,minazi:maxazi)

        real moon_mag,moon_mag_thr
        real*8 bi_coeff(2,2),ritopo,rjtopo,fi,fj
        logical l_solar_eclipse, l_binary, l_phase
        logical l_terrain_following

        integer mode_cloud_mask ! 1 is ignore cloud mask
                                ! 2 is display cloud mask differences
                                ! 3 perform cloud mask clearing
                                ! 4 no action
                                ! 5 no action

        write(6,*)' subroutine calc_allsky...'

!       if(mode_cloud_mask .ge. 4)goto 60 ! for testing

!       set up rayleigh scattering
        write(6,*)' rayleigh pf = ',rayleigh_pf(170.)
        write(6,*)'       wa       od       tr      refl    rmag'

!       bodhaine et. al "on rayleigh optical depth calculations", jtech
!       http://web.gps.caltech.edu/~vijay/papers/rayleigh_scattering/bodhaine-etal-99.pdf
        od_550 = .097069
!       od_550 = .140

        do ic = 1,nc
           ext_g(ic) = od_550 * (0.55 / wa(ic))**4.0
           tr = trans(ext_g(ic))
           rmag = -log10(tr) * 2.5
           refl = opac(ext_g(ic)) * 0.5 * rayleigh_pf(170.)
           write(6,1)ic,wa(ic),ext_g(ic),tr,refl,rmag
1          format(i3,f9.3,4f9.4)
        enddo ! ic

        iobs = nint(ri_obs)
        jobs = nint(rj_obs)

        if(iobs .lt. 1 .or. iobs .gt. nx_l .or. 
     1     jobs .lt. 1 .or. jobs .gt. ny_l)then
          write(6,*)' error, viewpoint outside of domain',iobs,jobs
          return
        endif

        isound = iobs
        jsound = jobs

        write(6,*)' albedo_sfc - calc_als',topo_albedo_2d(:,iobs,jobs)

        swi_obs = swi_2d(iobs,jobs) ! initialize
        write(6,*)' swi_2d at observer location 1 = ',swi_obs,'w/m^2'

        pw_obs = pw_2d(iobs,jobs) ! initialize
        write(6,*)' pw_2d at observer location = ',pw_obs

        eobsl = eobsc(iobs,jobs)
        write(6,*)' eobsl = ',eobsl
        eobsc_sky = 0. ! initialize

        clear_rad_c(:,:,:) = 0. ! initialize

        write(6,*)' range of clwc_3d is',minval(clwc_3d),maxval(clwc_3d)
        write(6,*)' range of cice_3d is',minval(cice_3d),maxval(cice_3d)
        write(6,*)' range of rain_3d is',minval(rain_3d),maxval(rain_3d)
        write(6,*)' range of snow_3d is',minval(snow_3d),maxval(snow_3d)
        write(6,*)' max top of cice_3d is',maxval(cice_3d(:,:,nz_l))
        write(6,*)' max top of clwc_3d is',maxval(clwc_3d(:,:,nz_l))
        write(6,*)' ssa values are ',ssa(:)

        if(angstrom_exp_a .eq. -99.)then
          angstrom_exp_a = 2.4 - (fcterm * 15.)
          if(angstrom_exp_a .lt. 0.25)then
            angstrom_exp_a = 0.25
            write(6,*)' using floor angstrom_exp_a:',angstrom_exp_a
          else
            write(6,*)' using fcterm angstrom_exp_a:',angstrom_exp_a
     1                                               ,fcterm
          endif
        else
          write(6,*)' using namelist angstrom_exp_a:',angstrom_exp_a
        endif

        write(6,*)' call get_cloud_rays...'

!         get line of sight from isound/jsound
          call get_cloud_rays(i4time_solar,clwc_3d,cice_3d
     1                     ,heights_3d                              ! i
     1                     ,rain_3d,snow_3d                         ! i
     1                     ,pres_3d,topo_sfc,topo                   ! i
     1                     ,topo_albedo_2d                          ! i
     1                     ,swi_2d                                  ! i
     1                     ,topo_swi,topo_albedo                    ! o
     1                     ,gtic,dtic,btic,emic                     ! o
     1                     ,topo_ri,topo_rj                         ! o
     1                     ,trace_ri,trace_rj                       ! o
     1                     ,swi_obs                                 ! o
!    1                     ,ghi_2d,dhi_2d                           ! o
     1                     ,aod_vrt,aod_2_cloud,aod_2_topo          ! o
     1                     ,dist_2_topo                             ! o
     1                     ,aod_ill,aod_ill_dir                     ! o
     1                     ,aod_tot,transm_obs,obs_glow_gnd         ! o
     1                     ,transm_3d,transm_4d                     ! o
     1                     ,r_cloud_3d,cloud_od,cloud_od_sp         ! o
     1                     ,cloud_od_sp_w                           ! o
     1                     ,r_cloud_rad,cloud_rad_c,cloud_rad_w     ! o
     1                     ,cloud_rad_c_nt                          ! o
     1                     ,clear_rad_c,clear_radf_c,patm           ! o
     1                     ,clear_rad_c_nt                          ! o
     1                     ,airmass_2_cloud_3d,airmass_2_topo_3d    ! o
     1                     ,htmsl,horz_dep,twi_0                    ! o
!    1                     ,solalt_limb_true                        ! o
     1                     ,htagl                                   ! i
     1                     ,aod_ref,ext_g                           ! i
     1                     ,nx_l,ny_l,nz_l,isound,jsound,newloc     ! i
     1                     ,ri_obs,rj_obs                           ! i
     1                     ,alt_a_roll,azi_a_roll                   ! i
     1                     ,sol_alt_2d,sol_azi_2d                   ! i
     1                     ,alt_norm                                ! i
     1                     ,moon_alt_2d,moon_azi_2d                 ! i
     1                     ,moon_mag,moon_mag_thr                   ! i
     1                     ,l_solar_eclipse,eobsc,rlat,rlon,lat,lon ! i
     1                     ,minalt,maxalt,minazi,maxazi             ! i
     1                     ,alt_scale,azi_scale                     ! i
     1                     ,l_binary,l_terrain_following            ! i
     1                     ,grid_spacing_m,r_missing_data           ! i
     1                     ,istatus)                                ! o

          if(istatus .ne. 1)then
             write(6,*)' bad status back from get_cloud_rays',istatus  
             return
          endif

          write(6,*)' return from get_cloud_rays: ',a9time
     1             ,' aod_vrt is ',aod_vrt

          write(6,*)' observer htagl/htmsl ',htagl,htmsl
          write(6,*)' swi_2d at observer location 2 = ',swi_obs,'w/m^2'

          solalt_limb_true = solar_alt + horz_dep
          write(6,*)' solalt_true = ',solar_alt
          write(6,*)' solalt_limb_true = ',solalt_limb_true

!         moon glow in cylindrical coordinates (add color info)?                   
          blog_moon_roll = 0.
!         if(moon_mag .lt. moon_mag_thr .and.
          if(.true.                     .and.
     1       alm      .gt. 0.           .and.
     1       l_solar_eclipse .eqv. .false.    )then
            write(6,*)' moon glow being calculated: ',alm,azm
            diam_deg = 0.5
            l_phase = .true.
            call great_circle(alm,azm,solar_alt,solar_az,gcdist,va)
            rill = (1. - cosd(elgms)) / 2.
            write(6,*)' vertex angle / rill is ',va,rill
            call get_glow_obj(i4time,alt_a_roll,azi_a_roll
     1                       ,minalt,maxalt,minazi,maxazi 
     1                       ,alt_scale,azi_scale
     1                       ,htmsl,patm
     1                       ,alm,azm,moon_mag,.false.
     1                       ,l_phase,rill,va
     1                       ,dum1,dum2,dum3
     1                       ,diam_deg,horz_dep,blog_moon_roll)

            write(6,*)' range of blog_moon_roll is',
     1          minval(blog_moon_roll),maxval(blog_moon_roll)
          endif

!         sun glow in cylindrical coordinates, treated as round?
          blog_sun_roll = 0.
          if(l_solar_eclipse .eqv. .true.)then
              if(eobsl .ge. 1.0)then ! totality
                  s_mag = -12.5
              else                   ! partial
                  s_mag = -26.74 - (log10(1.0-eobsl))*2.5
              endif
              diam_deg = 8.0    
              write(6,*)' corona glow being calculated: '
              l_phase = .false.
              call get_glow_obj(i4time,alt_a_roll,azi_a_roll
     1                     ,minalt,maxalt,minazi,maxazi 
     1                     ,alt_scale,azi_scale
     1                     ,htmsl,patm
     1                     ,solar_alt,solar_az,s_mag,l_solar_eclipse
     1                     ,l_phase,rill,va
     1                     ,alm,azm,emag ! used for solar eclipse
     1                     ,diam_deg,horz_dep,blog_sun_roll)
              write(6,31)minval(blog_sun_roll)
     1                  ,maxval(blog_sun_roll),eobsl,s_mag
 31           format('  range of blog_sun_roll (corona) is',4f10.4)
          endif

          if(emag .lt. 1.0)then
            s_mag = -26.74      ! at mean distance
            diam_deg = 0.533239 ! at mean distance
            write(6,*)' sun glow being calculated: '
            l_phase = .false.
            call get_glow_obj(i4time,alt_a_roll,azi_a_roll
     1                     ,minalt,maxalt,minazi,maxazi 
     1                     ,alt_scale,azi_scale
     1                     ,htmsl,patm
     1                     ,solar_alt,solar_az,s_mag,l_solar_eclipse
     1                     ,l_phase,rill,va
     1                     ,alm,azm,emag ! used for solar eclipse
     1                     ,diam_deg,horz_dep,blog_sun_roll)
            write(6,*)' range of blog_sun_roll is',
     1          minval(blog_sun_roll),maxval(blog_sun_roll),diam_deg
          else
            write(6,*)' total eclipse, sun glow not calculated'
          endif

!         if(solar_alt .ge. 0.)then
!         if(.true.)then
          i4_elapsed = ishow_timer()
!         write(6,*)' call get_skyglow_cyl for sun or moon...'

!         get all sky for cyl   
          ni_cyl = maxalt - minalt + 1
          nj_cyl = maxazi - minazi + 1

          i4_elapsed = ishow_timer()

          if(solar_alt .gt. -3.)then
                call get_idx(solar_alt,minalt,alt_scale,ialt_sun)
                call get_idx(solar_az ,minazi,azi_scale,jazi_sun)
                ialt_sun = ialt_sun - minalt + 1
                jazi_sun = jazi_sun - minazi + 1
                write(6,*)' solar_alt,minalt,alt_scale = '
     1                     ,solar_alt,minalt,alt_scale
                write(6,*)' ialt_sun,jazi_sun = ',ialt_sun,jazi_sun
          endif

          i4_elapsed = ishow_timer()

          do j = minazi,maxazi
          do i = minalt,maxalt

!            initialize
             trace_solalt(i,j) = solar_alt
             topo_solalt(i,j) = solar_alt
             topo_solazi(i,j) = solar_az
             topo_lf(i,j) = 1.
             topo_sc(i,j) = 0.

             itrace = nint(trace_ri(i,j))
             jtrace = nint(trace_rj(i,j))

             if(alt_a_roll(i,j) .lt. 0.)then
               if(htmsl .gt. 7000.)then
                 if(itrace .ge. 1 .and. itrace .le. nx_l .and.
     1              jtrace .ge. 1 .and. jtrace .le. ny_l)then
                   trace_solalt(i,j) = sol_alt_2d(itrace,jtrace)
                   eobsc_sky(i,j) = eobsc(itrace,jtrace)
                   itrace_updated = 1
                 else
                   itrace_updated = 0
                 endif

!                else ! outside domain - estimate solar altitude at sea level

                 if(htmsl .gt. 100000. .or. itrace_updated .eq. 0)then

!                  estimate solar altitude at sea level
                   htmsl = htagl + topo_sfc
                   if(azi_a_roll(i,j) .eq. 90.)then
                      iverbose = 1
                   else
                      iverbose = 0
                   endif

                   if(dist_2_topo(i,j) .gt. 0.)then ! ray hits sfc
                     call latlon_ray(rlat,rlon,htmsl,alt_a_roll(i,j)
     1               ,azi_a_roll(i,j),dist_2_topo(i,j),rlat_sfc,rlon_sfc
     1               ,iverbose,istatus)
                     solzen_sfc = angdist(rlat_sfc,solar_lat 
     1                                   ,rlon_sfc-solar_lon)
                     trace_solalt(i,j) = 90. - solzen_sfc
                     if(iverbose .eq. 1)then
                       write(6,41)i,j,alt_a_roll(i,j),azi_a_roll(i,j)
     1                          ,dist_2_topo(i,j),rlat_sfc,rlon_sfc
     1                          ,trace_solalt(i,j)
41                    format(
     1                    ' setting trace_solalt hi or outside domain: '
     1                      ,2i6,2f9.3,f13.1,3f9.3)
                     endif

                   else                             ! ray misses sfc
                     gcdist_km = (horz_dep * rpd * earth_radius) / 1000.
                     call razm_lat_lon_gm(rlat,rlon,gcdist_km 
     1                     ,azi_a_roll(i,j),rlat_limb,rlon_limb,istatus)
                     solzen_sfc = angdist(rlat_limb,solar_lat 
     1                                   ,rlon_limb-solar_lon)
                     trace_solalt(i,j) = 90. - solzen_sfc
                     if(iverbose .eq. 1)then
                       write(6,42)i,j,alt_a_roll(i,j),azi_a_roll(i,j)
     1                          ,gcdist_km,rlat_limb,rlon_limb
     1                          ,trace_solalt(i,j)
42                     format(' setting trace_solalt off limb: '
     1                      ,2i6,2f9.3,f10.1,3f9.3)
                     endif

                   endif ! dist_2_topo > 0
                 endif ! htmsl > 100000. or not updated
               endif ! htmsl > 7000.
             endif ! alt_a_roll < 0.

             ritopo = topo_ri(i,j)
             rjtopo = topo_rj(i,j)
             itopo = nint(ritopo)
             jtopo = nint(rjtopo)

             if(itopo .ge. 1 .and. itopo .le. nx_l .and.
     1          jtopo .ge. 1 .and. jtopo .le. ny_l)then
                    topo_solalt(i,j) = sol_alt_2d(itopo,jtopo)
                    topo_solazi(i,j) = sol_azi_2d(itopo,jtopo)
                    eobsc_sky(i,j) = eobsc(itopo,jtopo)
                    if(.true.)then ! bilin interp
                        i1 = max(min(int(ritopo),nx_l-1),1)
                        fi = ritopo - i1
                        i2 = i1+1
                        j1 = max(min(int(rjtopo),ny_l-1),1)
                        fj = rjtopo - j1
                        j2 = j1+1

                        bi_coeff(1,1) = (1d0-fi) * (1d0-fj)
                        bi_coeff(2,1) = fi       * (1d0-fj)
                        bi_coeff(1,2) = (1d0-fi) *      fj 
                        bi_coeff(2,2) = fi       *      fj 
                        topo_lf(i,j) = nint(sum(bi_coeff(:,:) 
     1                               * land_frac(i1:i2,j1:j2)))
                        topo_lat(i,j) = sum(bi_coeff(:,:) 
     1                                * lat(i1:i2,j1:j2))
                        topo_lon(i,j) = sum(bi_coeff(:,:) 
     1                                * lon(i1:i2,j1:j2))
                    else
                        topo_lf(i,j) = land_frac(itopo,jtopo) 
                        topo_lat(i,j) = lat(itopo,jtopo) 
                        topo_lon(i,j) = lon(itopo,jtopo) 
                    endif

                    if(snow_cover(itopo,jtopo) .ne. r_missing_data)then
                        topo_sc(i,j) = snow_cover(itopo,jtopo)
                    else
                        topo_sc(i,j) = 0.
                    endif
                    du_a(i,j) = du_2d(itopo,jtopo) ! bilin interp?
                    aod_a(i,j) = aod_2d(itopo,jtopo) ! bilin interp?
             endif

             if(alt_a_roll(i,j) .eq. -90. .and. j .eq. 1)then
                write(6,*)' nadir info'
                write(6,*)' i/j/lat/lon/lf',itopo,jtopo
     1                       ,topo_lat(i,j),topo_lon(i,j),topo_lf(i,j)
     1                       ,topo_albedo(:,i,j)
             endif

          enddo ! i
          enddo ! j

          topo_lat_r4(:,:) = topo_lat(:,:)
          topo_lon_r4(:,:) = topo_lon(:,:)

          if(.true.)then
             call drape_topo_albedo(
     1                 topo_lat,topo_lon                            ! i
     1                ,nc,minalt,maxalt,minazi,maxazi               ! i
     1                ,grid_spacing_m,r_missing_data                ! i
     1                ,topo_albedo)                                 ! i/o

             i4_elapsed = ishow_timer()
          endif

          if(.false.)then
             call add_topo_mask(topo_ri,topo_rj,lat,lon,nx_l,ny_l   ! i
     1                          ,nc,minalt,maxalt,minazi,maxazi     ! i
     1                          ,grid_spacing_m                     ! i
     1                          ,topo_albedo)                       ! i/o
          endif

          if(l_binary .eqv. .false.)then
              write(6,*)' call get_sky_rgb with cyl data'
              if(htagl .ge. 1000e3)then                    ! high alt
                corr1_in = 9.1             ! for high scattering angle
              elseif(htagl .eq. 300.)then                  ! bao
                if(solar_alt .lt. 30.)then
                  corr1_in = 9.2                            ! low sun
                elseif(solar_alt .gt. 60.)then
                  corr1_in = 8.9                            ! high sun
                else
                  corr1_in = 9.2 - (solar_alt-30.)*(0.3/30.)! med sun
                endif
              else
                corr1_in = 9.1
                if(solar_alt .lt. 0.)corr1_in = 9.26 ! volcanic value
              endif
              corr1_in = corr1_in ! - log10(exposure)

              write(6,*)' corr1 in calc_allsky ',corr1_in

!             this can be more accurate by using surface pressure
              patm_sfc = ztopsa(topo(isound,jsound)) / 1013.25
              patm_sfc = max(patm_sfc,patm)

              write(6,*)' patm/patm_sfc in calc_allsky ',patm,patm_sfc
              write(6,*)' range of clear_rad_c 3 =',
     1                   minval(clear_rad_c(3,:,:)),
     1                   maxval(clear_rad_c(3,:,:))

              write(6,*)' call get_sky_rgb with cyl data '
     1                   ,l_solar_eclipse
              call get_sky_rgb(r_cloud_3d    ! cloud opacity
     1                    ,cloud_od          ! cloud optical depth
     1                    ,cloud_od_sp,nsp   ! cloud species optical depth
     1                    ,cloud_od_sp_w     ! front weighted species
     1                    ,r_cloud_rad       ! cloud solar transmittance
     1                    ,cloud_rad_c       ! cloud solar transmittance / color
     1                    ,cloud_rad_c_nt    ! cloud night color
     1                    ,cloud_rad_w       ! cloud solar transmittance * rad
     1                    ,cloud_sfc_c       ! cld rad from sfc lighting (sru)
     1                    ,clear_rad_c       ! clear sky illumination by sun     
     1                    ,l_solar_eclipse,i4time_solar,rlat,rlon,eobsl
     1                    ,clear_radf_c      ! clear sky frac illumination by sun     
     1                    ,patm,patm_sfc,htmsl
     1                    ,clear_rad_c_nt    ! night sky brightness
!    1                    ,blog_v_roll       ! skyglow
     1                    ,blog_sun_roll     ! sunglow
     1                    ,blog_moon_roll    ! moonglow
     1                    ,glow_stars        ! starglow
     1                    ,ext_g,ext_o,wa,nc,day_int0
     1                    ,aod_vrt,aod_ref 
     1                    ,transm_obs,obs_glow_gnd ! observer illumination
     1                    ,ialt_sun,jazi_sun ! sun location
     1                    ,airmass_2_cloud_3d      
     1                    ,airmass_2_topo_3d      
     1                    ,swi_obs           ! sw at ground below observer 
     1                    ,topo_swi,topo_albedo,gtic,dtic,btic,emic
     1                    ,topo_albedo_2d(:,isound,jsound)
     1                    ,topo_lat_r4,topo_lon_r4,topo_lf,topo_sc       ! i
     1                    ,aod_2_cloud,aod_2_topo,aod_ill,aod_ill_dir
     1                    ,aod_tot
     1                    ,dist_2_topo,topo_solalt,topo_solazi
     1                    ,trace_solalt,eobsc_sky
     1                    ,alt_a_roll,azi_a_roll                         ! i   
     1                    ,ni_cyl,nj_cyl,alt_scale,azi_scale  
     1                    ,solar_alt,solar_az                            ! i
     1                    ,solar_lat,solar_lon,r_au                      ! i
     1                    ,minalt,maxalt,minazi,maxazi                   ! i
     1                    ,twi_0,horz_dep
     1                    ,solalt_limb_true
     1                    ,alm,azm,moon_mag  ! moon alt/az/mag
     1                    ,corr1_in,exposure
     1                    ,elong_a                                       ! o
     1                    ,sky_rgb_cyl,sky_sprad,sky_reflectance)        ! o   

!             add bounds and scaling to rgb values. final image should have
!             either a max of 128 or a green average of 64, whichever represents
!             a more conservative adjustment.              
              skymax = maxval(sky_rgb_cyl)
              if(skymax .lt. 128. .and. skymax .gt. 0.)then
                final_scaling = 128./skymax
              else
                final_scaling = 1.0               
              endif

              skyave = sum(sky_rgb_cyl(2,:,:)) / float(ni_cyl*nj_cyl)
              if(skyave .lt. 64. .and. skyave .gt. 0.)then
                final_scaling = min(final_scaling,64./skyave)
              else
                final_scaling = 1.0               
              endif

              write(6,51)skymax,skyave,final_scaling
51            format('  max/ave/final_scaling = ',3f9.3)             

              do j = minazi,maxazi
              do i = minalt,maxalt

                sky_rgb_cyl(:,i,j) = sky_rgb_cyl(:,i,j)*final_scaling

                if(.false.)then ! preserve colors in bright saturated areas
                  colmax = maxval(sky_rgb_cyl(:,i,j))
                  if(colmax .gt. 255.)then
                    col_ratio = 255. / colmax
                    sky_rgb_cyl(:,i,j) = sky_rgb_cyl(:,i,j) * col_ratio                                        
                  endif                   
                else ! simple clip at 255.
                    sky_rgb_cyl(:,i,j) = 
     1                  max(min(sky_rgb_cyl(:,i,j),255.),0.)       
                endif

                if(.false.)then ! grayscale reflectance image
                    sky_rgb_cyl(:,i,j) = sky_reflectance(1,i,j) * 1000.
                    sky_rgb_cyl(:,i,j) = min(sky_rgb_cyl(:,i,j),255.)
                endif

              enddo ! i
              enddo ! j

          else
              continue ! use cloud_od to drive categorical output

          endif ! l_binary

60        continue

!         cloud mask section
          if(mode_cloud_mask .eq. 1)then
              write(6,*)' skipping camera_cloud_mask processing'
          elseif(mode_cloud_mask .eq. 2 .or. mode_cloud_mask .eq. 3)then
              write(6,*)' showing cloud mask differences'
      
              write(6,*)' camera cloud mask...'
              do ialt = maxalt,minalt,-10
                 write(6,61)ialt,(camera_cloud_mask(ialt,ia)
     1                          ,ia=minazi,maxazi,6)
61               format(1x,i3,1x,300i1)
              enddo ! ialt
      
              do j = minazi,maxazi
              do i = minalt,maxalt
                  if(cloud_od(i,j) .gt. 0.3)then
                      sim_cloud_mask(i,j) = 2
                  else
                      sim_cloud_mask(i,j) = 1
                  endif
              enddo ! i
              enddo ! j
              write(6,*)' simulated cloud mask...'
              do ialt = maxalt,minalt,-10
                 write(6,61)ialt,(sim_cloud_mask(ialt,ia)
     1                          ,ia=minazi,maxazi,6)
              enddo ! ialt
      
              write(6,*)' difference cloud mask...'
              nclear_mask = 0
              do j = minazi,maxazi
              do i = minalt,maxalt
                  if(camera_cloud_mask(i,j) .eq. 2 .and.
     1               sim_cloud_mask(i,j)    .eq. 1        )then
                      diff_cloud_mask(i,j) = 2 ! add cloud (potentially)
                  elseif(camera_cloud_mask(i,j) .eq. 1 .and.
     1                   sim_cloud_mask(i,j)    .eq. 2        )then
                      diff_cloud_mask(i,j) = 1 ! clearing
                      nclear_mask = nclear_mask + 1
                  else
                      diff_cloud_mask(i,j) = 0 ! no change
                  endif
!                 if(j .eq. (9*maxazi)/10)then
!                     write(6,*)i,j,camera_cloud_mask(i,j)
!    1                         ,sim_cloud_mask(i,j),diff_cloud_mask(i,j)
!                 endif
              enddo ! i
              enddo ! j
              write(6,*)' nclear_mask = ',nclear_mask
              do ialt = maxalt,minalt,-10
                 write(6,61)ialt,(diff_cloud_mask(ialt,ia)
     1                          ,ia=minazi,maxazi,6)
              enddo ! ialt

              if(mode_cloud_mask .eq. 3)then
                 write(6,*)' performing camera_cloud_mask clearing'
                 call clear_3d_mask(
     1                 minalt,maxalt,minazi,maxazi                  ! i
     1                ,alt_scale,azi_scale                          ! i
     1                ,htmsl,ri_obs,rj_obs                          ! i
     1                ,lat,lon,nx_l,ny_l,nz_l,grid_spacing_m        ! i
     1                ,diff_cloud_mask                              ! i
     1                ,sky_rgb_cyl)                                 ! i/o
              endif

          else ! mode_cloud_mask = 4 or 5 (correlation)
              continue

          endif

          write(6,*)' end of calc_allsky...'

          return
          end

          subroutine latlon_ray(rlat,rlon,htmsl,altray,aziray,tpdist
     1                         ,tlat,tlon,iverbose,istatus)

!         find intersection lat,lon of a light ray with the earth's surface

          use mem_namelist, only: r_missing_data,earth_radius

          include 'trigd.inc'
          
          real*8 tpdist,dslant1_h
          real*8 dsdst,dradius_start,daltray,dcurvat2

          dcurvat2(dsdst,dradius_start,daltray) = 
     1        sqrt(dsdst**2 + dradius_start**2 
     1   - (2d0*dsdst*dradius_start*cosd(90d0+daltray))) - dradius_start   

!         dradius_start = earth_radius + htmsl

!         sinarc = sind(90d0+dble(altray))*tpdist/dble(earth_radius)

!         dz1_h = dcurvat2(dslant1_h,dradius_start,dble(altray))   
          gc_deg = asind(cosd(altray) * tpdist / (earth_radius)) 
!    1           *     tpdist / (dradius_start+dz1_h)) 
!    1           *     tpdist / (earth_radius)) 
          
          ycost = cosd(aziray)

!         obtain latlon from spherical geometry
          tlat = asind(ycost*sind(gc_deg)*cosd(rlat) 
     1         + sind(rlat)*cosd(gc_deg))

          cosdlon = (cosd(gc_deg) - sind(rlat)*sind(tlat)) 
     1            / (cosd(rlat)*cosd(tlat))
          if(abs(cosdlon).gt.1.)cosdlon=sign(1.,cosdlon)
          dlon = acosd(cosdlon)
                  
          if(aziray.ge..0.and.aziray.le.180.)then ! east
             tlon=rlon+dlon
          else ! west
             tlon=rlon-dlon
          endif

          if(iverbose .eq. 1)then
            write(6,11)rlat,rlon,htmsl,altray,aziray,tpdist
     1                ,gc_deg,dradius_start,dlon,tlat,tlon
11          format('   latlon_ray ',2f9.3,f13.1,2f9.3,f13.1,f9.3,
     1                                    f13.1,3f9.3)
          endif

          istatus = 1
          
          return
          end

          subroutine clear_3d_mask(
     1                  minalt,maxalt,minazi,maxazi             ! i
     1                 ,alt_scale,azi_scale                     ! i
     1                 ,htmsl,ri_obs,rj_obs                     ! i
     1                 ,lat,lon,nx_l,ny_l,nz_l,grid_spacing_m   ! i
     1                 ,diff_cloud_mask                         ! i
     1                 ,sky_rgb_cyl)                            ! i/o

          use mem_allsky ! 3d model grids

          real lat(nx_l,ny_l)
          real lon(nx_l,ny_l)
          real sky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi) ! observed variable

          integer diff_cloud_mask(minalt:maxalt,minazi:maxazi)

!         color scheme for grid points
!         yellow  - no clearing of gridpoint having cloud
!!        green   - no clearing of gridpoint having no cloud         
!         green   - potential gridpoints for cloud addition
!         magenta - clearing of gridpoint having cloud
!         orange  - clearing of gridpoint having no cloud

          logical l_draw_grid /.true./

!         observer location
          ht_obs = htmsl
          call bilinear_laps(ri_obs,rj_obs,nx_l,ny_l,lat,rlat_obs)
          call bilinear_laps(ri_obs,rj_obs,nx_l,ny_l,lon,rlon_obs)

          write(6,*)' clear_3d_mask: observer lat/lon',rlat_obs,rlon_obs

          nclear = 0
          nmag = 0; nred = 0; ngrn = 0; nyel = 0

          do i = 1,nx_l
          do j = 1,ny_l      
          do k = 1,nz_l

!           determine altitude / azimuth of this grid point
            ht_grid = heights_3d(i,j,k)
            rlat_grid = lat(i,j)
            rlon_grid = lon(i,j)

!           approximate for now that radar travels similarly to light
            call latlon_to_radar(rlat_grid,rlon_grid,ht_grid    ! i
     1                          ,azimuth,slant_range,elev       ! o
     1                          ,rlat_obs,rlon_obs,ht_obs)      ! i
            
!           determine approximate alt/az range covered by the grid box for
!           option that prevents aliasing            

!           hgrid_alt =
!           hgrid_azi =
!           vgrid_alt =             

            ialt = nint(elev    / alt_scale)
            jazi = nint(azimuth / azi_scale)
           
            if(ialt .ge. minalt .and. ialt .le. maxalt .and.
     1         jazi .ge. minazi .and. jazi .le. maxazi      )then

              cond_max = max(clwc_3d(i,j,k),cice_3d(i,j,k)
     1                      ,rain_3d(i,j,k),snow_3d(i,j,k))
               
              if(diff_cloud_mask(ialt,jazi) .eq. 1)then ! clear this grid point
                nclear = nclear + 1
                write(6,11)nclear,i,j,k,rlat_grid,rlon_grid
     1                    ,elev,azimuth,slant_range,ialt,jazi
11              format('   clearing',4i5,2f9.2,f10.2,f10.4,f10.0,2i5)          

                if(l_draw_grid)then
                  if(cond_max .gt. 0.)then 
                    sky_rgb_cyl(0,ialt,jazi) = 255. ! magenta
                    sky_rgb_cyl(1,ialt,jazi) = 0.
                    sky_rgb_cyl(2,ialt,jazi) = 255.
                    nmag = nmag + 1
                  else
                    sky_rgb_cyl(0,ialt,jazi) = 180. ! red / gray
                    sky_rgb_cyl(1,ialt,jazi) = 140.
                    sky_rgb_cyl(2,ialt,jazi) = 140.
                    nred = nred + 1
                  endif
                endif
               
                clwc_3d(i,j,k) = 0.
                cice_3d(i,j,k) = 0.
                rain_3d(i,j,k) = 0.
                snow_3d(i,j,k) = 0.

              else ! no clearing
                if(i .eq. nint(ri_obs) .and. k .eq. nz_l/2)then ! informational
                  write(6,12)nclear,i,j,k,rlat_grid,rlon_grid
     1                      ,elev,azimuth,slant_range,ialt,jazi
     1                      ,diff_cloud_mask(ialt,jazi)
12                format('no clearing',4i5,2f9.2,f10.2,f10.4,f10.0,2i5
     1                                ,i2)
                endif

                if(l_draw_grid)then
                  if(cond_max .gt. 1e-7)then 
                    rint = 200. + 30. * log10(cond_max / 1e-5)
                    rint = min(max(rint,0.),255.)
                    sky_rgb_cyl(0,ialt,jazi) = rint ! yellow
                    sky_rgb_cyl(1,ialt,jazi) = rint
                    sky_rgb_cyl(2,ialt,jazi) = 0.
                  else
!                   sky_rgb_cyl(0,ialt,jazi) = 120.
!                   sky_rgb_cyl(1,ialt,jazi) = 205. ! green
!                   sky_rgb_cyl(2,ialt,jazi) = 120.
                  endif
                endif

              endif
              if(diff_cloud_mask(ialt,jazi) .eq. 2)then ! potential addition
                if(l_draw_grid)then
                  sky_rgb_cyl(0,ialt,jazi) = 0.
                  sky_rgb_cyl(1,ialt,jazi) = 255. ! green
                  sky_rgb_cyl(2,ialt,jazi) = 0.
                endif
              endif
            endif ! in camera domain

          enddo ! k
          enddo ! j
          enddo ! i

          write(6,*)' clear_3d_mask: nclear = ',nclear
          if(l_draw_grid)then
            write(6,*)' n red/mag = ',nred,nmag
          endif

          return
          end

        subroutine add_topo_mask(topo_ri,topo_rj,lat,lon,ni,nj          ! i
     1                          ,nc,minalt,maxalt,minazi,maxazi         ! i
     1                          ,grid_spacing_m                         ! i
     1                          ,topo_albedo)                           ! i/o

!       add high resolution vector albedo information to 'topo_albedo'

        real topo_ri(minalt:maxalt,minazi:maxazi)
        real topo_rj(minalt:maxalt,minazi:maxazi)
        real topo_albedo(nc,minalt:maxalt,minazi:maxazi)

        parameter (n_rect = 2)

        real   lat(ni,nj)
        real   lon(ni,nj)
        real   lat_rect(n_rect)
        real   lon_rect(n_rect)

        real*8 width_rect(n_rect)
        real*8 height_rect(n_rect)
        real*8 azi_rect(n_rect)
        real*8 albedo_rect(n_rect)

        real*8 ri_rect,rj_rect,ridelt,rjdelt

        write(6,*)' adding topo mask'

        lat_rect(1) = 39.849
        lon_rect(1) = -104.674
        width_rect(1) = 3000. / grid_spacing_m
        height_rect(1) = 150. / grid_spacing_m
        albedo_rect(1) = 1.

        lat_rect(2) = 39.849
        lon_rect(2) = -104.674
        width_rect(2) = 2800. / grid_spacing_m
        height_rect(2) = 150. / grid_spacing_m
        albedo_rect(2) = 0.

        do m = 1,n_rect
          call latlon_to_rlapsgrid(lat_rect(m),lon_rect(m),lat,lon,ni,nj
     1                            ,r4i_rect,r4j_rect,istatus)
          ri_rect = r4i_rect
          rj_rect = r4j_rect

          do ialt = minalt,maxalt
          do jazi = minazi,maxazi
            ridelt = topo_ri(ialt,jazi) - ri_rect
            rjdelt = topo_rj(ialt,jazi) - rj_rect

!           these can be rotated by azimuth, presently is n-s
            if(abs(ridelt) .lt. height_rect(m)/2.  .and.
     1         abs(rjdelt) .lt. width_rect(m)/2.        )then  
              topo_albedo(:,ialt,jazi) = albedo_rect(m)
!             if(ialt .eq. -15 .or. ialt .eq. -25)then
              if(jazi .eq. 40)then
                 write(6,1)m,ialt,jazi,topo_ri(ialt,jazi)
     1                                ,topo_rj(ialt,jazi)
1                format(i2,2i5,2f9.2,' runway added')               
              endif
            else
!             if(ialt .eq. -15 .or. ialt .eq. -25)then
              if(jazi .eq. 40)then
                 write(6,2)m,ialt,jazi,topo_ri(ialt,jazi)
     1                                ,topo_rj(ialt,jazi)
2                format(i2,2i5,2f9.2,' runway not added')               
              endif
              continue
            endif

          enddo ! jazi
          enddo ! ialt
        enddo ! m
       
        return
        end

        subroutine drape_topo_albedo(
     1                 topo_lat,topo_lon                             ! i
     1                ,nc,minalt,maxalt,minazi,maxazi                ! i
     1                ,grid_spacing_m,r_missing_data                 ! i
     1                ,topo_albedo)                                  ! i/o

!       add high resolution vector albedo information to 'topo_albedo'

        use ppm

        real topo_albedo(nc,minalt:maxalt,minazi:maxazi)

        real nlattile

        real*8 topo_lat(minalt:maxalt,minazi:maxazi)
        real*8 topo_lon(minalt:maxalt,minazi:maxazi)

        real*8 arglat,arglon,pix_we,pix_sn,rlat_start,rlon_start
     1        ,pix_latlon_we,pix_latlon_sn

        integer, allocatable :: img(:,:,:)                   

        character*255 directory, file_dc, file
        character*20 adum,ctype,cropname
        integer u,u_out
        logical l_there_dc
        integer counts(3)
        real bi_coeff(2,2)

        write(6,*)' drape topo albedo'

        call get_directory('static',directory,len_dir)

        do itile = 1,100
        
          write(cropname,1)itile
1         format('vhires_crop',i2.2,'.ppm')        
          file_dc=trim(directory)//trim(cropname)

          inquire(file=trim(file_dc),exist=l_there_dc)
          write(6,*)' file being inquired is ',trim(file_dc),' '
     1                                        ,l_there_dc      
          if(l_there_dc)then                   ! naip/descartes data
            pix_latlon_we = 1. / 865.954              
            pix_latlon_sn = 1. / 1130.26
            file = trim(file_dc)
            perimeter = 0.05
          else
            write(6,*)
     1        ' naip data not present - returning from drape_topo'
            return
          endif

          write(6,*)' open for reading ',trim(file)

!         read section of descartes image in ppm format
          u = 11
          open(u,file=trim(file),status='old',err=999)
          read(u,*)   
          read(u,*)iwidth,iheight
          rewind(u)
          write(6,*)' dynamic dims ',nc,iwidth,iheight
          allocate(img(nc,iwidth,iheight))
          call read_ppm (u,img,nc,iwidth,iheight)
          read(u,*,end=999)adum,adum,pixdgw
          read(u,*,end=999)adum,adum,pixdgh
          read(u,*,end=999)adum,adum,nlattile
          read(u,*,end=999)adum,adum,wlontile
          read(u,*,end=999)adum,adum,ctype
          write(6,*)' comment info ',pixdgw,pixdgh,nlattile,wlontile
     1                              ,ctype
          close(u)

          pix_latlon_sn = 1. / pixdgh
          pix_latlon_we = 1. / pixdgw

          drape_res_m = 110000. * pix_latlon_sn
          write(6,*)' drape data resolution (m) is ',drape_res_m

          rlat_start = nlattile ! 40.4176  ! rnorth                         
          rlon_start = wlontile ! -105.806 ! west                          
          rlat_end = rlat_start - float(iheight) * pix_latlon_sn
          rlon_end = rlon_start + float(iwidth)  * pix_latlon_we

          write(6,*)' input image lat range ',rlat_start,rlat_end
          write(6,*)' input image lon range ',rlon_start,rlon_end

!         for each ray topo lat/lon, interpolate from image albedo arrays
          write(6,*)
          write(6,*)' sample drape points'
          write(6,*)'                     ialt   jazi      arglat    '
     1             ,'arglon  pix_we   pix_sn'
        
          do ialt = minalt,maxalt
          do jazi = minazi,maxazi
            arglat = topo_lat(ialt,jazi)
            arglon = topo_lon(ialt,jazi)
            if(arglat .ne. r_missing_data .and.
     1         arglon .ne. r_missing_data       )then
              pix_we = (arglon - rlon_start) / pix_latlon_we
              pix_sn = (rlat_start - arglat) / pix_latlon_sn

              i1 = int(pix_we)
              fi = pix_we - i1
              i2 = i1+1
              j1 = int(pix_sn)
              fj = pix_sn - j1
              j2 = j1+1
            
              in = nint(pix_we)
              jn = nint(pix_sn)

!             nearest neighbor for now
              if(i1 .ge. 1 .and. i1 .le. iwidth-1 .and.
     1           j1 .ge. 1 .and. j1 .le. iheight-1      )then

                if(.false.)then
                  counts(:) = img(:,in,jn)
                else ! bilinear interpolation
                  bi_coeff(1,1) = (1.-fi) * (1.-fj)
                  bi_coeff(2,1) = fi      * (1.-fj)
                  bi_coeff(1,2) = (1.-fi) *     fj 
                  bi_coeff(2,2) = fi      *     fj 
                  do ic = 1,3
                    counts(ic) = sum(bi_coeff(:,:)
     1                         * img(ic,i1:i2,j1:j2))
                  enddo ! ic
                endif

                do ic = 1,3
                  if(counts(ic) .le. 179.)then
                    topo_albedo(ic,ialt,jazi) = .0010 * counts(ic)
     1                                      + 6.92e-11 * counts(ic)**4.0
                  elseif(counts(ic) .le. 255.)then
                    topo_albedo(ic,ialt,jazi)
     1                         = 0.25 + (counts(ic) - 179.)    * .00258
     1                                + (counts(ic) - 179.)**2 * .000092
                  else ! bad / default value
                    topo_albedo(ic,ialt,jazi) = 0.2
                  endif
                enddo ! ic

                if(trim(ctype) .eq. 'usda')then
                  topo_albedo(:,ialt,jazi) = topo_albedo(:,ialt,jazi)
     1                                     * 0.7
                endif

                if(jazi .eq. minazi)then
                  write(6,11)ialt,jazi,arglat,arglon
     1           ,pix_we,pix_sn,in,jn,counts(:),topo_albedo(:,ialt,jazi)
11                format(' sample drape point',2i7,2f12.6,2f9.2,2i6,2x
     1                                        ,3i6,3f9.3)
                endif ! printing point

              else ! outside image
                if(jazi .eq. minazi)then
!                 write(6,21)ialt,jazi,pix_we,pix_sn,arglat,arglon
21                format(' outside image ',2i7,2f13.3,2f10.4)
                endif              

              endif ! inside image
            endif ! valid lat/lon
          enddo ! jazi
          enddo ! ialt

          deallocate(img)

          write(6,*)' end of tile ',itile

          i4_elapsed = ishow_timer()

        enddo ! itile

!       normal end
        istatus = 1
        goto 9999

!       error condition
999     istatus = 0
        write(6,*)' error in drape_topo_albedo'

9999    if(allocated(img))deallocate(img)

        return
        end
      
        subroutine get_camsite(rlat,rlon,site)

        character*10 site

        if(rlat .gt. 39.9)then
           site = 'dsrc'
        else
           site = 'nrel'
        endif

        return
        end

        subroutine diffimg(img1,img2,nc,ni,nj,a_t,b_t,avecorr
     1                    ,isun,jsun,idb,fname)

        use ppm

        integer img1(nc,ni,nj),img2(nc,ni,nj),imgdiff(nc,ni,nj)
        real a_t(nc),b_t(nc)
        real tmp1(ni,nj),tmp2(ni,nj)

        character*(*)fname

        write(6,*)
        write(6,*)' subroutine diffimg...'
        write(6,*)' cam-sim difference image ',trim(fname),ni,nj
        write(6,*)' a_t is ',a_t
        write(6,*)' b_t is ',b_t

!       determine difference wrt regression line
        do ic = 1,nc
            tmp1(:,:) = img1(ic,:,:)       ! sim
            tmp2(:,:) = img2(ic,:,:)       ! cam
            if(avecorr .lt. 0.5)then
                imgdiff(ic,:,:) = 128 +        ! wrt scaled images
     1            nint( tmp2(:,:) - tmp1(:,:) ) 
            else
                imgdiff(ic,:,:) = 128 +        ! wrt regression line
     1            nint( tmp2(:,:) - (a_t(ic) * tmp1(:,:) + b_t(ic)) )
            endif
!           imgdiff(ic,:,:) = img1(ic,:,:) ! test simulated image
!           imgdiff(ic,:,:) = img2(ic,:,:) ! test camera image

            imgdiff(ic,:,:) = min(max(imgdiff(ic,:,:),0),255)
        enddo ! ic

        call writeppm3matrix(imgdiff(1,:,:),imgdiff(2,:,:)
     1                      ,imgdiff(3,:,:)
     1                      ,trim(fname))

        if(idb .ge. 1)then
            write(6,*)' top sim  = ',img1(:,ni,1)
            write(6,*)' top cam  = ',img2(:,ni,1)
            write(6,*)' top diff = ',imgdiff(:,ni,1)
            write(6,*)' solar sim  = ',img1(:,isun,jsun)
            write(6,*)' solar cam  = ',img2(:,isun,jsun)
            write(6,*)' solar diff = ',imgdiff(:,isun,jsun)
        endif

        return
        end

        subroutine compare_camera(iloop,rlat,rlon,nc                  ! i
     1                           ,minalt,maxalt,minazi,maxazi         ! i
     1                           ,alt_scale,azi_scale,camera_rgb      ! i
     1                           ,i4time_solar,isun,jsun,idb          ! i
     1                           ,alt_a_roll,elong_a,r_missing_data   ! i
     1                           ,sky_rgb_cyl                         ! i/o
     1                           ,correlation,a_t,b_t                 ! o
     1                           ,istatus)                            ! o

        include 'trigd.inc'

        real camera_rgb(nc,minalt:maxalt,minazi:maxazi)
        real camera_rgbf(nc,minalt:maxalt,minazi:maxazi)
        real wt_a(minalt:maxalt,minazi:maxazi)
        real alt_a_roll(minalt:maxalt,minazi:maxazi)
        real elong_a(minalt:maxalt,minazi:maxazi)
        real sky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi) ! observed variable
        real correlation(nc),xbar_a(nc),ybar_a(nc),a_t(nc),b_t(nc)
        character*10 site
        

        if(iloop .eq. 1 .or. .true.)then
            write(6,*)' calling get_camera_image: ',rlat,rlon
            mode_cam = 2
            call get_camsite(rlat,rlon,site)
            if(trim(site) .eq. 'dsrc')then
                  i4time_camera = i4time_solar - 60
            else
                  i4time_camera = i4time_solar
            endif
            call get_camera_image(minalt,maxalt,minazi,maxazi,nc,         ! i
     1                                alt_scale,azi_scale,                ! i
     1                                i4time_camera,mode_cam,             ! i
     1                                site,                               ! i
     1                                camera_rgb,istatus)                 ! o
            if(istatus .eq. 0)then
                  write(6,*)' return from calc_allsky sans camera image'
                  return
            endif

            i4_elapsed = ishow_timer()
            wt_a(:,:) = 1.0

!           subsample flagged array to account to weight cylindrical projection             
            iskip_max = 6
            do ialt = maxalt,minalt,-1
                  if(ialt .eq. maxalt)then
                    iskip = iskip_max
                  else
                    iskip = nint(1./cosd(alt_a_roll(ialt,minazi)))
                    iskip = min(iskip,iskip_max)
                  endif
                  do jazi = minazi,maxazi
                    if(jazi .eq. (jazi/iskip)*iskip)then
                      camera_rgbf(:,ialt,jazi) = camera_rgb(:,ialt,jazi)
                    else
!                     camera_rgbf(:,ialt,jazi) = r_missing_data
                      wt_a(ialt,jazi) = 0.0
                    endif
                    if(elong_a(ialt,jazi) .lt. 3.0)then
!                     camera_rgbf(:,ialt,jazi) = r_missing_data
                      wt_a(ialt,jazi) = 0.0
                    endif
                  enddo ! jazi
            enddo ! ialt
        else
            write(6,*)' skip cam img read - use saved array'
        endif ! iloop

        cam_checksum = sum(min(camera_rgbf,255.))
        write(6,*)' camera_rgbf checksum = ',cam_checksum
        if(cam_checksum .gt. 0. .and. cam_checksum .lt. 1.0)then
            write(6,*)' warning: 0 < cam_checksum < 1'
            istatus = 0
            return
        endif

!       we have a choice of doing the weighting by setting pixels
!       to missing, or by passing in variable weights
        write(6,*)' performing correlation calculation (stats_2d)'
        do ic = 1,nc
            call stats_2d(maxalt-minalt+1,maxazi-minazi+1
     1                       ,camera_rgbf(ic,:,:),sky_rgb_cyl(ic-1,:,:)
     1                       ,wt_a,a_t(ic),b_t(ic),xbar_a(ic),ybar_a(ic)
     1                       ,correlation(ic),bias,std
     1                       ,r_missing_data,istatus)
        enddo

        corr_scaling = sum(ybar_a(:)) / sum(xbar_a(:))
        write(6,*)' correlation scaling (sim/cam) is '
     1                  ,corr_scaling
        if(corr_scaling .lt. 1.0)then
            write(6,*)' brightening simulated image to match camera'
            sky_rgb_cyl(:,:,:) =
     1                 min(sky_rgb_cyl(:,:,:)/corr_scaling,255.)
            a_t(:) = a_t(:) / corr_scaling
            b_t(:) = b_t(:) / corr_scaling
        else
            write(6,*)' simulated image is brighter than camera'
        endif

        if(idb .eq. 1)then
           write(6,*)' solar sim = ',sky_rgb_cyl(:,isun,jsun)
           write(6,*)' solar cam = ',camera_rgbf(:,isun,jsun)
        endif

        write(6,*)' a_t is ',a_t
        write(6,*)' b_t is ',b_t

        return
        end
