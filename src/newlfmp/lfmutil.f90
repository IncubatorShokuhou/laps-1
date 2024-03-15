subroutine lfm_namelist(lfmprd_dir)

! read namelist variables.
! default values are set in lfmgrid

   use lfmgrid

   implicit none

   integer :: istatus

   character(len=*) :: lfmprd_dir
   character(len=256) :: nlfile

   namelist /lfmpost_nl/ out_cdf, out_grib, out_v5d &
      , make_micro, make_firewx, make_points &
      , verbose, realtime, write_to_lapsdir &
      , make_donefile &
      , precip_dt, c_m2z, advection_time &
      , n3d_pts_thr, p3d_pts_thr &
      , wrf_version

   nlfile = trim(lfmprd_dir)//'/../static/lfmpost.nl'
   open (1, file=trim(nlfile), form='formatted', status='old', iostat=istatus)
   if (istatus /= 0) then
      print *, 'error opening namelist file: ', trim(nlfile)
      stop
   end if
   read (1, nml=lfmpost_nl)
   close (1)

   return
end

!===============================================================================

subroutine get_native_dims(mtype, filename, nx, ny, nz, istatus)

   use lfmgrid, only: n3d_pts_thr, large_ngrid

   implicit none

   integer :: nx, ny, nz, istatus

   character(len=*) :: mtype, filename

! istatus: 1=good return, 0=error/bad return
   istatus = 1  ! assume good return

   select case (trim(mtype))
   case ('mm5')
      call get_mm5_dims(filename, nx, ny, nz)
   case ('wrf')
      call get_wrf_dims(filename, nx, ny, nz, istatus)
   case ('nmm')
      call get_nmm_dims(filename, nx, ny, nz)
   case ('st4')
      call get_st4_dims(filename, nx, ny, nz)
   end select

! set native grid size threshold based on an 8gb machine
   write (6, *) ' # of 3d native grid points, large_ngrid thresh = ', nx*ny*nz, n3d_pts_thr
   if (nx*ny*nz > n3d_pts_thr) then
      write (6, *) ' large native grid - process reduced set of 3-d fields'
      large_ngrid = .true.
!  large_pgrid = .true.
   end if

   return
end

!===============================================================================

subroutine get_laps_static(laps_data_root)

   use lfmgrid

   implicit none

   include 'netcdf.inc'

   integer, parameter :: nprmax = 150
   integer :: icode, ncid, nid, nk_laps, istatus, i

   real, dimension(nprmax) :: pressures

   character(len=256) :: stlaps
   character(len=132) :: gridtype
   character(len=*) :: laps_data_root

   namelist /pressures_nl/ pressures

   stlaps = trim(laps_data_root)//'/static/static.nest7grid'
   icode = nf_open(trim(stlaps), nf_nowrite, ncid)
   if (icode /= 0) then
      print *, 'laps static file not found: ', trim(stlaps)
      stop
   end if
   icode = nf_inq_dimid(ncid, 'x', nid)
   icode = nf_inq_dimlen(ncid, nid, lx)
   icode = nf_inq_dimid(ncid, 'y', nid)
   icode = nf_inq_dimlen(ncid, nid, ly)
   icode = nf_inq_varid(ncid, 'grid_type', nid)
   icode = nf_get_var_text(ncid, nid, gridtype)
   icode = nf_inq_varid(ncid, 'grid_spacing', nid)
   icode = nf_get_var_real(ncid, nid, grid_spacing)
   icode = nf_inq_varid(ncid, 'latin1', nid)
   icode = nf_get_var_real(ncid, nid, truelat1)
   icode = nf_inq_varid(ncid, 'latin2', nid)
   icode = nf_get_var_real(ncid, nid, truelat2)
   icode = nf_inq_varid(ncid, 'lov', nid)
   icode = nf_get_var_real(ncid, nid, stdlon)
   if (stdlon > 180.) stdlon = stdlon - 360.

   allocate (llat(lx, ly), llon(lx, ly))
   icode = nf_inq_varid(ncid, 'lat', nid)
   icode = nf_get_var_real(ncid, nid, llat)
   icode = nf_inq_varid(ncid, 'lon', nid)
   icode = nf_get_var_real(ncid, nid, llon)
   icode = nf_close(ncid)

   if (gridtype(1:5) == 'polar') then
      projection = 'polar stereographic'
   elseif (gridtype(1:14) == 'secant lambert') then
      projection = 'lambert conformal'
   elseif (gridtype(1:18) == 'tangential lambert') then
      projection = 'lambert conformal'
   elseif (gridtype(1:18) == 'mercator') then
      projection = 'mercator'
   else
      print *, 'unrecognized laps grid type: ', trim(gridtype)
      stop
   end if

   pressures = rmsg
   stlaps = trim(laps_data_root)//'/static/pressures.nl'
   open (2, file=trim(stlaps), status='old', err=900)
   read (2, pressures_nl, err=901)
   close (2)

! determine number of isobaric levels by checking pressure data.
!   (assume there is at least one level, and that data is ordered correctly)

   do lz = 1, nprmax
      if (pressures(lz + 1) > 200000.) exit
   end do

   allocate (lprs(lz), lprsl(lz))
   lprs(1:lz) = pressures(1:lz)
   lprsl(1:lz) = alog(lprs(1:lz))

   write (6, *) ' # of 3d laps grid points, large_pgrid thresh = ', lx*ly*lz, p3d_pts_thr
   if (lx*ly*lz > p3d_pts_thr) then
      write (6, *) ' large laps grid - process reduced set of 3-d fields'
      large_pgrid = .true.
!  large_ngrid = .true.
   end if

   return

900 continue
   print *, 'could not open laps namelist file: ', trim(stlaps)
   stop

901 continue
   print *, 'error reading laps namelist file: ', trim(stlaps)
   stop

end

!===============================================================================

subroutine set_laps_projection

   use lfmgrid

   implicit none

   select case (trim(projection))
   case ('lambert conformal')
      call map_set(proj_lc, llat(1, 1), llon(1, 1), grid_spacing &
                   , stdlon, truelat1, truelat2, lx, ly, proj)
   case ('polar stereographic')
      call map_set(proj_ps, llat(1, 1), llon(1, 1), grid_spacing &
                   , stdlon, truelat1, truelat2, lx, ly, proj)
   case ('mercator')
      call map_set(proj_merc, llat(1, 1), llon(1, 1), grid_spacing &
                   , stdlon, truelat1, truelat2, lx, ly, proj)
   end select

   return
end

!===============================================================================

subroutine fill_native_grid

   use lfmgrid

   implicit none

! fill native grids.

   select case (trim(mtype))
   case ('mm5')
      call fill_mm5_grid
   case ('wrf')
      call fill_wrf_grid
   case ('nmm')
      call fill_nmm_grid
   case ('st4')
      call fill_st4_grid
   end select

! set map utilities for use by horizontal interpolation.

   select case (trim(nprojection))
   case ('lambert conformal')
      call map_set(proj_lc, nlat(1, 1), nlon(1, 1), ngrid_spacingx &
                   , nstdlon, ntruelat1, ntruelat2, nx, ny, proj)
   case ('polar stereographic')
      call map_set(proj_ps, nlat(1, 1), nlon(1, 1), ngrid_spacingx &
                   , nstdlon, ntruelat1, ntruelat2, nx, ny, proj)
   case ('mercator')
      call map_set(proj_merc, nlat(1, 1), nlon(1, 1), ngrid_spacingx &
                   , nstdlon, ntruelat1, ntruelat2, nx, ny, proj)
   case ('rotated lat-lon')
   end select

   return
end

!===============================================================================

subroutine lfm_derived

   use lfmgrid
   use constants
   use cloud_rad ! cloud radiation and microphysics parameters

   implicit none

   integer :: i, j, k, k500, k700, k850, k1000
   real :: dz, dtdz, dewpt, relhum, potential_temp, eq_potential_temp, heatindex &
           , redp_lvl = 1500.
   real, allocatable, dimension(:, :) :: pcp_06, pcp_init, fallen_precip_type
   real, allocatable, dimension(:, :, :) :: htdsig, hrhsig, hthetasig, hthetaesig, hzsigf &
                                            , pcpmr_sig, rhodrysig, tvsig, rhomoistsig

!beka
   real :: dx(lx, ly), dy(lx, ly)
   real :: r_missing_data, aglhgt
   real :: stefan_boltzmann, b_olr, eff_emissivity, b_tau
   real :: bt_flux_equiv
   real :: coeff
   real :: a1, b1, c1, d1, e1, alpha(lx, ly)
   real :: const_lwp, const_iwp, const_rwp, const_swp, const_gwp
   real :: const_lwp_bks, const_iwp_bks, const_rwp_bks, const_swp_bks, const_gwp_bks
   real :: clwc2alpha, cice2alpha, rain2alpha, snow2alpha, pice2alpha

!beka

   character*256 directory
   real rspacing_dum
   character*125 comment_2d
   character*10 units_2d
   character*50 ext
   integer len_dir, istatus
   character*3 var_2d
   real :: ldf(lx, ly), lat(lx, ly), lon(lx, ly), avg(lx, ly), static_albedo(lx, ly)
   real :: windspeed(lx, ly), soil_moist(lx, ly), snow_cover(lx, ly)
   real :: intrain(lx, ly), intsnow(lx, ly), intgraupel(lx, ly), intliq700(lx, ly)
   real :: const_pcpadj

   integer ::        ismoist, isnow
   integer ::        status

   integer i4_elapsed, ishow_timer

!cj variables to compute omega, added 6/21/2007
   real :: tvprs

! variables for low level reflectivity
   integer :: klow, khigh
   real :: ht_1km, frack
   real :: umean(lx, ly), vmean(lx, ly), ustorm(lx, ly), vstorm(lx, ly), ushear(lx, ly), vshear(lx, ly)
   real :: array_buf(lx, ly), array_out(lx, ly)

   real :: ghi_ratio(lx, ly)
   real :: arg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!beka!!!!!!!!!!!!! obtaining ldf, lat and lon !!!!!!!!!!!!!!!!!!!!

   ext = 'nest7grid'

   call get_directory(ext, directory, len_dir)
   var_2d = 'ldf'
   call rd_laps_static(directory, ext, lx, ly, 1, var_2d, &
                       units_2d, comment_2d, ldf, &
                       rspacing_dum, istatus)
   var_2d = 'avg'
   call rd_laps_static(directory, ext, lx, ly, 1, var_2d, &
                       units_2d, comment_2d, avg, &
                       rspacing_dum, istatus)

   var_2d = 'lat'
   call rd_laps_static(directory, ext, lx, ly, 1, var_2d, &
                       units_2d, comment_2d, lat, &
                       rspacing_dum, istatus)

   var_2d = 'lon'
   call rd_laps_static(directory, ext, lx, ly, 1, var_2d, &
                       units_2d, comment_2d, lon, &
                       rspacing_dum, istatus)

   call get_static_field_interp('albedo', laps_valtime, lx, ly &
                                , static_albedo, istatus)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! extrapolate a surface temp if not available.

   if (.not. large_ngrid) then
      if (maxval(tsfc) > 1000.) then
         if (verbose) then
            print *, ' '
            print *, 'tsfc not available...will extrapolate from lowest level...'
         end if
         do j = 1, ly
         do i = 1, lx

! compute dz between lowest sigma layer and 2m level

            dz = hzsig(i, j, 1) - (zsfc(i, j) + 2.)

! compute the lapse rate of temp for the two lowest sigma levels.
! use of the hypsometric equation did not work well for these thin layers.

            dtdz = (htsig(i, j, 2) - htsig(i, j, 1)) &
                   /(hzsig(i, j, 2) - hzsig(i, j, 1))
            tsfc(i, j) = htsig(i, j, 1) - dtdz*dz
         end do
         end do
      end if
   end if

! fill other missing surface fields from lowest model level.

   print *, 'min/max mrsfc (model sfc grids) =', minval(mrsfc), maxval(mrsfc)
   print *, 'substitute lowest sigma level for surface mixing ratio'

!if (maxval(mrsfc) > 1000.)mrsfc(:,:)=hmrsig(:,:,1)
   if (.true.) mrsfc(:, :) = hmrsig(:, :, 1)

   print *, 'min/max mrsfc (model sfc grids) =', minval(mrsfc), maxval(mrsfc)

   if (l_process_uv) then
      if (maxval(usfc) > 1000.) usfc(:, :) = husig(:, :, 1)
      if (maxval(vsfc) > 1000.) vsfc(:, :) = hvsig(:, :, 1)
   end if
   if (l_process_w) then
      if (maxval(wsfc) > 1000.) wsfc(:, :) = hwsig(:, :, 1)
   end if

! other derived surface fields.
   where (mrsfc < zero_thresh) mrsfc = zero_thresh

   if (maxval(mrsfc) < 1000.) then
      do j = 1, ly
      do i = 1, lx
         rhsfc(i, j) = min(relhum(tsfc(i, j), mrsfc(i, j), psfc(i, j)), 1.)
         tdsfc(i, j) = dewpt(tsfc(i, j), rhsfc(i, j))
         thetasfc(i, j) = potential_temp(tsfc(i, j), psfc(i, j))
         thetaesfc(i, j) = eq_potential_temp(tsfc(i, j), psfc(i, j), mrsfc(i, j), rhsfc(i, j))
         rhsfc(i, j) = rhsfc(i, j)*100.
      end do
      end do
   end if

! generate reduced pressure for laps usage.

   if (.not. large_pgrid) then
      if (verbose) then
         print *, ' '
         print *, 'reducing pressure to ', redp_lvl, ' meters'
      end if
      call interp_press_to_z(lprs, zprs, redp_lvl, redp, lx, ly, lz)

! use same routine to interpolate sea-level pressure.  this method
! is used in lieu of original reduction routine, because it will
! keep the msl field consistent with the height field, which has been
! properly reduced.  it also produces a bit smoother field over the mountains.

      if (verbose) print *, 'reducing pressure to msl...'
      call interp_press_to_z(lprs, zprs, 0., pmsl, lx, ly, lz)
   end if

! microphysical surface variables.

   if (make_micro) then
      where (hmrsig < zero_thresh) hmrsig = zero_thresh
      where (hcldliqmr_sig < zero_thresh) hcldliqmr_sig = 0.0
      where (hcldicemr_sig < zero_thresh) hcldicemr_sig = 0.0
      where (hsnowmr_sig < zero_thresh) hsnowmr_sig = 0.0
      where (hrainmr_sig < zero_thresh) hrainmr_sig = 0.0
      where (hgraupelmr_sig < zero_thresh) hgraupelmr_sig = 0.0
      clwmrsfc(:, :) = hcldliqmr_sig(:, :, 1)
      icemrsfc(:, :) = hcldicemr_sig(:, :, 1)
      snowmrsfc(:, :) = hsnowmr_sig(:, :, 1)
      rainmrsfc(:, :) = hrainmr_sig(:, :, 1)
      graupmrsfc(:, :) = hgraupelmr_sig(:, :, 1)

!if(.not. large_ngrid)then
      if (.true.) then

! generate cloud fields.

         cldbase = rmsg
         cldtop = rmsg
         ceiling = rmsg  ! unlimited ceiling
         if (maxval(hcldliqmr_sig) < rmsg .and. maxval(hcldicemr_sig) < rmsg .and. &
             maxval(hsnowmr_sig) < rmsg) then
            cldamt = 0.
            call lfmclouds(lx, ly, nz, ngrid_spacingx, hcldliqmr_sig, hcldicemr_sig &
                           , hsnowmr_sig, hzsig, zsfc, cldbase, cldtop, ceiling, cldamt)
         else
            cldamt = rmsg
         end if

! generate integrated hydrometeor fields

         allocate (pcpmr_sig(lx, ly, nz), rhodrysig(lx, ly, nz)) ! for possible use
         pcpmr_sig = 0.
         if (maxval(hsnowmr_sig) < rmsg) pcpmr_sig = pcpmr_sig + hsnowmr_sig
         if (maxval(hrainmr_sig) < rmsg) pcpmr_sig = pcpmr_sig + hrainmr_sig
         if (maxval(hgraupelmr_sig) < rmsg) pcpmr_sig = pcpmr_sig + hgraupelmr_sig
         where (pcpmr_sig < zero_thresh) pcpmr_sig = 0.
         rhodrysig = hpsig/(r*htsig)
         call lfm_integrated_liquid(lx, ly, nz, hcldliqmr_sig, hcldicemr_sig, hmrsig, &
                                    hrainmr_sig, hsnowmr_sig, hgraupelmr_sig, rhodrysig, hzsig, hpsig, zsfc, &
                                    intliqwater, intliq700, intcldice, totpcpwater, intrain, intsnow, intgraupel, rmsg)
      end if ! large_ngrid

      if (.true.) then ! supercedes binary cldamt calculated earlier
         ! units of integrated vertical quantities (e.g. intliqwater) is meters
         ! for liquid we can augment this noting that tau = (3 * intliqwater * rholiq) / (2 * rho * radius)
         ! cldamt should be related to opacity and optical depth
         ! we might use a new module for this derived from these files
         ! in  'src/lib/cloud/get_cloud_rad.f90':
         ! and 'src/lib/modules/module_cloud_rad.f90'

         const_lwp = (1.5*rholiq)/(rholiq*reff_clwc)
         const_lwp_bks = const_lwp*bksct_eff_clwc

         const_iwp = (1.5*rholiq)/(rholiq*reff_cice)
         const_iwp_bks = const_iwp*bksct_eff_cice

         const_rwp = (1.5*rholiq)/(rholiq*reff_rain)
         const_rwp_bks = const_rwp*bksct_eff_rain

         const_swp = (1.5*rholiq)/(rhosnow*reff_snow)
         const_swp_bks = const_swp*bksct_eff_snow

         const_gwp = (1.5*rholiq)/(rhograupel*reff_graupel)
         const_gwp_bks = const_gwp*bksct_eff_graupel

         write (6, *) ' simvis constants ', const_lwp_bks, const_iwp_bks, const_rwp_bks, const_swp_bks, const_gwp_bks

         do j = 1, ly
            do i = 1, lx
               if (intliqwater(i, j) < rmsg .and. intcldice(i, j) < rmsg) then
!        cldamt(i,j) = (intliqwater(i,j) * 100.) ** 0.3333333 ! this might work empirically if small integrated liq/ice
!        cldamt(i,j) = min(cldamt(i,j),1.0)                   ! values have smaller mean diameters, so the dependence
!                                                             ! is a weaker one
!        liquid water is assumed to have smaller particles contributing to a larger constant
!        cloud amount is opacity of cloud liquid and cloud ice hydrometeors
                  cldamt(i, j) = 1.-(exp(-(const_lwp*intliqwater(i, j) + const_iwp*intcldice(i, j))))

!        rain, snow, and graupel are added for the cldalb & simvis computation
                  b_tau = const_lwp_bks*intliqwater(i, j) + const_iwp_bks*intcldice(i, j) + &
                          const_rwp_bks*intrain(i, j) + const_swp_bks*intsnow(i, j) + const_gwp_bks*intgraupel(i, j)
                  cldalb(i, j) = b_tau/(b_tau + 1.)
                  simvis(i, j) = cldalb(i, j) + (1.-cldalb(i, j))**2*(static_albedo(i, j)/(1.-cldalb(i, j)*static_albedo(i, j)))
               end if
            end do ! i
         end do ! j
      end if ! true

      if (.not. large_ngrid) then
         if (verbose) then
            print *, ' '
            print *, 'called integrated_liquid.'
            print *, 'min/max pcpmr_sig =', minval(pcpmr_sig), maxval(pcpmr_sig)
            print *, 'min/max mrsig =', minval(hmrsig), maxval(hmrsig)
            print *, 'min/max rhodrysig =', minval(rhodrysig), maxval(rhodrysig)
            print *, 'min/max zsig =', minval(hzsig), maxval(hzsig)
         end if
         deallocate (pcpmr_sig, rhodrysig)
      end if ! large_ngrid

! generate derived reflectivity fields.
      if (.true.) then ! always run this part
         if (verbose) then
            print *, ' '
            print *, 'calling lfm_reflectivity.'
            print *, 'min/max hrainmr_sig =', minval(hrainmr_sig), maxval(hrainmr_sig)
            print *, 'min/max hcldicemr_sig =', minval(hcldicemr_sig), maxval(hcldicemr_sig)
            print *, 'min/max hsnowmr_sig =', minval(hsnowmr_sig), maxval(hsnowmr_sig)
            print *, 'min/max hgraupelmr_sig =', minval(hgraupelmr_sig), maxval(hgraupelmr_sig)
         end if
         if (minval(hrainmr_sig) < rmsg .and. minval(hcldicemr_sig) < rmsg .and. &
             minval(hsnowmr_sig) < rmsg .and. minval(hgraupelmr_sig) < rmsg) then
            allocate (tvsig(lx, ly, nz), rhomoistsig(lx, ly, nz))
            tvsig = htsig*(1.+0.61*hmrsig)
            rhomoistsig = hpsig/(r*tvsig)
            call lfm_reflectivity(lx, ly, nz, rhomoistsig, hzsig &
                                  , hrainmr_sig, hcldicemr_sig, hsnowmr_sig, hgraupelmr_sig &
                                  , hrefl_sig, hzdr_sig, hldr_sig, htsig, hpsig, hmrsig &
                                  , max_refl, echo_tops)

!     determine low level reflectivity
            if (.false.) then ! use lowest sigma level (sfc)
               refl_sfc(:, :) = hrefl_sig(:, :, 1)
            elseif (large_ngrid) then ! assume a certain sigma level is 1km agl
               refl_sfc(:, :) = hrefl_sig(:, :, 8)
            else ! interpolate to 1km agl
               do j = 1, ly
               do i = 1, lx
                  ht_1km = zsfc(i, j) + 1000.
                  do k = 1, nz - 1
                     if (hzsig(i, j, k) .le. ht_1km .and. hzsig(i, j, k + 1) .ge. ht_1km) then
                        frack = (ht_1km - hzsig(i, j, k))/(hzsig(i, j, k + 1) - hzsig(i, j, k))
                        klow = k
                        khigh = klow + 1
                        refl_sfc(i, j) = hrefl_sig(i, j, klow)*(1.-frack) + hrefl_sig(i, j, khigh)*frack
                     end if
                  end do ! k
               end do ! i
               end do ! j
            end if

!     run advection routine on reflectivity
            if (fcsttime .le. i4_adv_pcp) then
               i4_elapsed = ishow_timer()

               ext = 'lmr'
               var_2d = 'r'
               call get_laps_2d(laps_reftime, ext, var_2d, units_2d, comment_2d, lx, ly, max_refl, istatus)
               if (istatus .eq. 1) then
                  write (6, *) ' obtained initial time max_refl field from lmr analysis lmr file'
               else
                  write (6, *) ' warning: could not find initial time max_refl field from lmr analysis lmr file'
               end if

               call mean_wind_bunkers(husig, hvsig, hzsig(:, :, 1), lx, ly, nz &     ! i
                                      , hzsig &     ! i
                                      , umean, vmean, ushear, vshear, ustorm, vstorm, status)    ! o
               write (6, *) ' advect surface reflectivity ', fcsttime, i4_adv_pcp
               call advect(umean, vmean, refl_sfc, array_buf, ngrid_spacingx &
                           , lx, ly, array_out, float(fcsttime), 1.0, lon, rmsg)
               refl_sfc = array_out
               write (6, *) ' advect max reflectivity ', fcsttime, i4_adv_pcp
               call advect(umean, vmean, max_refl, array_buf, ngrid_spacingx &
                           , lx, ly, array_out, float(fcsttime), 1.0, lon, rmsg)
               max_refl = array_out
            else
               write (6, *) ' using model output reflectivity (no advection)', fcsttime, i4_adv_pcp
            end if

            if (fcsttime .le. i4_adv_cld) then
               i4_elapsed = ishow_timer()
               write (6, *) ' advect 3d clouds ', fcsttime, i4_adv_cld
               do k = 1, nz
                  call advect(husig(:, :, k), hvsig(:, :, k), hcldliqmr_sig(:, :, k), array_buf, ngrid_spacingx &
                              , lx, ly, array_out, float(fcsttime), 1.0, lon, rmsg)
                  hcldliqmr_sig(:, :, k) = array_out(:, :)

                  call advect(husig(:, :, k), hvsig(:, :, k), hcldicemr_sig(:, :, k), array_buf, ngrid_spacingx &
                              , lx, ly, array_out, float(fcsttime), 1.0, lon, rmsg)
                  hcldicemr_sig(:, :, k) = array_out(:, :)
               end do ! k
            else
               write (6, *) ' using model output clouds (no advection)', fcsttime, i4_adv_cld
            end if

            i4_elapsed = ishow_timer()

!     conversion of microphysical cloud mixing ratios to concentrations
            hcldliqmr_sig(:, :, :) = hcldliqmr_sig(:, :, :)*rhomoistsig(:, :, :)
            hcldicemr_sig(:, :, :) = hcldicemr_sig(:, :, :)*rhomoistsig(:, :, :)

            deallocate (tvsig, rhomoistsig)
            if (verbose) then
               print *, ' '
               print *, 'after lfm_reflectivity.'
               print *, 'min/max hrefl_sig =', minval(hrefl_sig), maxval(hrefl_sig)
               print *, 'min/max hzdr_sig =', minval(hzdr_sig), maxval(hzdr_sig)
               print *, 'min/max hldr_sig =', minval(hldr_sig), maxval(hldr_sig)
               print *, 'min/max refl_sfc =', minval(refl_sfc), maxval(refl_sfc)
            end if
         else
            hrefl_sig = rmsg
            hzdr_sig = rmsg
            hldr_sig = rmsg
            max_refl = rmsg
            echo_tops = rmsg
         end if

! do some qc on cloud species fields.

         where (cldliqmr_prs < zero_thresh) cldliqmr_prs = 0.
         where (cldicemr_prs < zero_thresh) cldicemr_prs = 0.
         where (snowmr_prs < zero_thresh) snowmr_prs = 0.
         where (rainmr_prs < zero_thresh) rainmr_prs = 0.
         where (graupelmr_prs < zero_thresh) graupelmr_prs = 0.

! generate surface and ua precip type.

         if (verbose) then
            print *, ' '
            print *, 'calling lfm_ua_pcptype.'
            print *, 'min/max hrainmr_sig =', minval(hrainmr_sig), maxval(hrainmr_sig)
            print *, 'min/max hsnowmr_sig =', minval(hsnowmr_sig), maxval(hsnowmr_sig)
            print *, 'min/max hgraupelmr_sig =', minval(hgraupelmr_sig), maxval(hgraupelmr_sig)
         end if
         if (minval(hrainmr_sig) < rmsg .and. &
             minval(hsnowmr_sig) < rmsg .and. minval(hgraupelmr_sig) < rmsg) then

            call lfm_sfc_pcptype(lx, ly, nz, hrainmr_sig, hsnowmr_sig, hgraupelmr_sig &
                                 , pcptype_sfc)
            if (verbose) then
               print *, ' '
               print *, 'after lfm_sfc_pcptype.'
               print *, 'min/max pcptype_sfc =', minval(pcptype_sfc), maxval(pcptype_sfc)
            end if

            if (large_pgrid .eqv. .false.) then
               call lfm_ua_pcptype(lx, ly, lz, rainmr_prs, snowmr_prs, graupelmr_prs &
                                   , pcptype_prs)
               if (verbose) then
                  print *, ' '
                  print *, 'after lfm_ua_pcptype.'
                  print *, 'min/max pcptype_prs =', minval(pcptype_prs), maxval(pcptype_prs)
               end if
            end if

         else
            pcptype_sfc = rmsg
            pcptype_prs = rmsg
         end if

!  convert mr to concentration for microphysical species?
!  to do this we'd need rho on the pressure grid. otherwise we can convert
!  on the model sig grid using rhosig, and this is now done at the end of
!  'lfm_reflectivity' (and right after returning)

      end if ! true

   end if ! make_micro

   i4_elapsed = ishow_timer()

! precip fields.

   do k = 1, lz
      if (lprs(k) == 100000) k1000 = k
      if (lprs(k) == 85000) k850 = k
      if (lprs(k) == 70000) k700 = k
      if (lprs(k) == 50000) k500 = k
   end do

   allocate (pcp_init(lx, ly))
   allocate (pcp_06(lx, ly))
   call fill_precip(pcp_init, pcp_06, snow_tot)

! fcsttime is in seconds
! change nmm, nmm zeroes accumulated precipitation every interval period

   if (fcsttime > 0) then
      if (trim(mtype) == 'nmm') then
         pcp_inc = pcp_tot
         print *, 'pcp_inc from pcp_tot: ', trim(mtype)
      else
         pcp_inc = pcp_tot - pcp_init
         print *, 'pcp_inc from pcp_init: ', trim(mtype)
         print *, 'fcsttime: ', fcsttime
         print *, maxval(pcp_inc), maxval(pcp_tot), maxval(pcp_init)
      end if
      where (pcp_inc < .0001) pcp_inc = 0.
      where (pcp_tot < .0001) pcp_tot = 0.

      if (.not. large_ngrid) then
         print *, 'min/max snow tot (previous) = ', minval(snow_tot), maxval(snow_tot)

!     allocate(fallen_precip_type(lx,ly))
!     call wintprec(htsig,hzsig,zprs,psfc,tsfc,tdsfc,zsfc,pcp_inc,lx,ly  &
!               ,nz,lz,k700,k850,k1000,pcp_inc,fallen_precip_type)
!     call snowfall(htsig,tsfc,pcp_inc,fallen_precip_type,lx,ly,nz,snow_inc,snow_tot)
         call snowfall(htsig, tsfc, pcp_inc, pcptype_sfc, lx, ly, nz, snow_inc, snow_tot)
!     deallocate(fallen_precip_type)
      end if

   else ! no incremental precip at time zero
      pcp_inc = r_missing_data

   end if

   if (trim(mtype) == 'nmm') then
      print *, 'pcp_inc update: ', trim(mtype)
      pcp_06 = pcp_tot
      pcp_tot = 0.
      pcp_tot = pcp_init + pcp_inc
   else
      print *, 'skip pcp_inc update: ', trim(mtype)
!pcp_tot = pcp_tot + pcp_inc
   end if

! write intermediate precip file for future use.

   open (1, file=trim(filename0), status='unknown', form='unformatted')
   if (trim(mtype) == 'nmm') then
      write (1) pcp_tot, pcp_06, snow_tot
      print *, 'write intermediate pcp file for model: ', trim(mtype)
   else
      write (1) pcp_tot, snow_tot
   end if
   close (1)
   if (verbose) then
      print *, ' '
      print *, 'writing new precip data to: ', trim(filename0)
   end if
   deallocate (pcp_init)
   deallocate (pcp_06)

! adjust precip/snow fields with a constant (for now) bias correction
   const_pcpadj = 1.00 ! 0.65
   where (pcp_inc .ne. r_missing_data) pcp_inc = pcp_inc*const_pcpadj
   where (pcp_tot .ne. r_missing_data) pcp_tot = pcp_tot*const_pcpadj
   where (snow_inc .ne. r_missing_data) snow_inc = snow_inc*const_pcpadj
   where (snow_tot .ne. r_missing_data) snow_tot = snow_tot*const_pcpadj

! temporary variables needed to derive some fields.

   allocate (hrhsig(lx, ly, nz), hthetasig(lx, ly, nz), hthetaesig(lx, ly, nz), hzsigf(lx, ly, nz + 1))
   do k = 1, nz
   do j = 1, ly
   do i = 1, lx
      hrhsig(i, j, k) = min(1., relhum(htsig(i, j, k), hmrsig(i, j, k), hpsig(i, j, k)))
      hthetasig(i, j, k) = potential_temp(htsig(i, j, k), hpsig(i, j, k))
      hthetaesig(i, j, k) = eq_potential_temp(htsig(i, j, k), hpsig(i, j, k) &
                                              , hmrsig(i, j, k), hrhsig(i, j, k))
   end do
   end do
   end do
   do k = 2, nz
      hzsigf(:, :, k) = (hzsig(:, :, k - 1) + hzsig(:, :, k))*0.5
   end do
   hzsigf(:, :, 1) = zsfc
   hzsigf(:, :, nz + 1) = 2.*hzsig(:, :, nz) - hzsigf(:, :, nz)

! generate pbl height if not available from model.

   if (maxval(pblhgt) > 999999.) then
      if (verbose) then
         print *, ' '
         print *, 'generating pbl height using theta and surface temp.'
      end if

      do j = 1, ly
      do i = 1, lx
         thetasfc(i, j) = potential_temp(tsfc(i, j), psfc(i, j))
      end do
      end do
      call model_pblhgt(hthetasig, thetasfc, hpsig, hzsig, zsfc, lx, ly, nz, pblhgt)
   end if
   where (pblhgt < 0.) pblhgt = 0.

   i4_elapsed = ishow_timer()

   if (.not. large_ngrid) then

! helicity, cape, cin, li.

      write (6, *) ' calculating helicity, stability indices'
      call helicity(husig, hvsig, hzsig, usfc, vsfc, zsfc, lx, ly, nz, srhel)

      i4_elapsed = ishow_timer()

      call updraft_helicity(husig, hvsig, hwsig, hzsig, hzsigf, zsfc, llat, llon, lx, ly, nz, uhel)
      print *, 'min/max uhel =', minval(uhel), maxval(uhel)

      call bulk_shear(husig, hvsig, hzsig, hzsigf, zsfc, llat, llon, lx, ly, nz, bshr)

      i4_elapsed = ishow_timer()

      call capecin(hpsig*0.01, htsig, hthetaesig, hthetasig, hrhsig &
                   , hzsigf, tprs, liftedind, cape, cin, k500, lx, ly, nz, lz)

! supercell composite parameter
      do i = 1, lx
      do j = 1, ly
         if (bshr(i, j) .ge. 20.) then
            arg = 1.0
         elseif (bshr(i, j) .ge. 10.) then
            arg = bshr(i, j)/20.
         else
            arg = 0.0
         end if
         scp(i, j) = cape(i, j)/1000.*(srhel(i, j)/50.)*arg
      end do ! j
      end do ! i
   end if ! large_ngrid

   deallocate (hthetasig, hthetaesig, hzsigf)

   i4_elapsed = ishow_timer()

! height of wet-bulb = 0c.

   write (6, *) ' calculating wet bulb zero height'
   allocate (htdsig(lx, ly, nz))
   if (.not. large_ngrid) then
      htdsig = htsig/((-rvolv*alog(hrhsig)*htsig) + 1.0)

      call height_tw(hpsig, hzsig, htsig, htdsig, psfc, zsfc, tsfc, tdsfc, mrsfc, pmsl &
                     , 0.0, ztw0, lx, ly, nz)

      i4_elapsed = ishow_timer()

      call height_tw(hpsig, hzsig, htsig, htdsig, psfc, zsfc, tsfc, tdsfc, mrsfc, pmsl &
                     , 1.3, ztw1, lx, ly, nz)

      i4_elapsed = ishow_timer()

      aglhgt = 80.
      call wind_agl(husig, hvsig, hzsig, aglhgt, zsfc, lx, ly, nz &
                    , u80, v80, r_missing_data)

   end if ! large_ngrid

! visibility.

   if (.true.) then
! a1 = 75.; b1 = 37.5; c1 = 1.5; d1 = 5.357; e1 = 1.2

      clwc2alpha = 1.5/(rholiq*reff_clwc)
      cice2alpha = 1.5/(rholiq*reff_cice)
      rain2alpha = 1.5/(rholiq*reff_rain)
      snow2alpha = 1.5/(rhosnow*reff_snow)
      pice2alpha = 1.5/(rhograupel*reff_graupel)

      a1 = clwc2alpha
      b1 = cice2alpha
      c1 = rain2alpha
      d1 = snow2alpha
      e1 = pice2alpha

! note that hydrometeor mixing ratios have previously been converted to content
      alpha(:, :) = a1*hcldliqmr_sig(:, :, 1) + b1*hcldicemr_sig(:, :, 1) &
                    + c1*hrainmr_sig(:, :, 1) + d1*hsnowmr_sig(:, :, 1) + e1*hgraupelmr_sig(:, :, 1)
! meteorological visibility, tau ~= 2.8, with lower bound added on alpha
      visibility(:, :) = 2.8/max(alpha(:, :), .00004)
   else
      visibility = 6000000.*(tsfc - tdsfc)/(rhsfc**1.75)  ! in meters
      where (visibility > 99990.) visibility = 99990.
   end if

! compute heat index if temp is above 80f (300k).

   do j = 1, ly
   do i = 1, lx
      if (tsfc(i, j) >= 300.) then
         heatind(i, j) = heatindex(tsfc(i, j), rhsfc(i, j))
      else
         heatind(i, j) = tsfc(i, j)
      end if
   end do
   end do

! fire weather indices.

   if (make_firewx .and. .not. large_ngrid) then
      call ventilation(husig, hvsig, hzsig, pblhgt, zsfc, lx, ly, nz, upbl, vpbl, vnt_index)

! mid-level haines index.

      call haines_layer(hpsig, htsig, htdsig, ham_index, lx, ly, nz, 850., 700.)

! high-level haines index.

      call haines_layer(hpsig, htsig, htdsig, hah_index, lx, ly, nz, 700., 500.)

! fosberg fwi.

      call fosberg_fwi(tsfc, rhsfc, usfc, vsfc, lx, ly, fwi_index)

! kelsch fwx.
      do i = 1, lx
      do j = 1, ly
         windspeed(i, j) = sqrt(usfc(i, j)**2 + vsfc(i, j)**2)
      end do ! j
      end do ! i
      ismoist = 0
      isnow = 0
      call lp_fire_danger(lx, ly, rhsfc, tsfc, windspeed, &
                          soil_moist, snow_cover, zsfc, ldf, &
                          ismoist, isnow, fwx_index, istatus)

   end if

   if (allocated(hrhsig)) deallocate (hrhsig)
   if (allocated(htdsig)) deallocate (htdsig)

!cj compute omega, added 6/21/2007
   if (.not. large_pgrid) then
      do k = 1, lz
      do j = 1, ly
      do i = 1, lx
         tvprs = tprs(i, j, k)*(1.+0.61*shprs(i, j, k))
!  yhl  omprs(i,j,k) = -(lprsl(k)*wprs(i,j,k)*grav)/(r*tvprs) * 100. ! pa/mb
         omprs(i, j, k) = -(lprs(k)*wprs(i, j, k)*grav)/(r*tvprs)  ! changed by steve. a, huiling yuan 9/11/2009
      end do
      end do
      end do
   end if

   if (fcsttime .le. i4_adv_cld) then
      call mean_wind_bunkers(husig, hvsig, hzsig(:, :, 1), lx, ly, nz &     ! i
                             , hzsig &     ! i
                             , umean, vmean, ushear, vshear, ustorm, vstorm, status)    ! o
      ext = 'lcv'
      var_2d = 'swi'
      call get_laps_2d(laps_reftime, ext, var_2d, units_2d, comment_2d, lx, ly, swdown, istatus)
      if (istatus .eq. 1) then
         write (6, *) ' obtained initial time swdown field from swi analysis lcv file'
      else
         write (6, *) ' warning: could not find initial time swdown field from swi analysis lcv file'
      end if
      write (6, *) ' advect solar radiation ', fcsttime, i4_adv_cld
      call advect(umean, vmean, swdown, array_buf, ngrid_spacingx &
                  , lx, ly, array_out, float(fcsttime), 1.0, lon, rmsg)
      write (6, *) ' range of advection input is ', minval(swdown), maxval(swdown)
      write (6, *) ' range of initial advection output is ', minval(array_out), maxval(array_out)

      write (6, *) ' calling get_ghi_ratio for times ', laps_reftime, laps_valtime
      call get_ghi_ratio(laps_reftime, laps_valtime, lat, lon, lx, ly, ghi_ratio)
      write (6, *) ' range of ghi_ratio is ', minval(ghi_ratio), maxval(ghi_ratio)

      swdown = array_out*ghi_ratio
      write (6, *) ' range of final advection output is ', minval(swdown), maxval(swdown)
   else
      if (fcsttime .eq. 0) then ! obtain initial time swdown field from swi analysis lcv files
         ext = 'lcv'
         var_2d = 'swi'
         call get_laps_2d(laps_reftime, ext, var_2d, units_2d, comment_2d, lx, ly, swdown, istatus)
         if (istatus .eq. 1) then
            write (6, *) ' obtained initial time swdown field from swi analysis lcv file'
         else
            write (6, *) ' warning: could not find initial time swdown field from swi analysis lcv file'
         end if
      else
         if (.true.) then ! post-process swdown field
            do j = 1, ly
            do i = 1, lx
               coeff = 0.93 + (cldamt(i, j)*0.5)
               swdown(i, j) = swdown(i, j)*coeff
            end do ! i
            end do ! j
         end if
      end if
   end if

   if (verbose .and. .not. large_ngrid) then
      print *, ' '
      print *, 'min/max tsfc      = ', minval(tsfc), maxval(tsfc)
      print *, 'min/max rhsfc     = ', minval(rhsfc), maxval(rhsfc)
      print *, 'min/max tdsfc     = ', minval(tdsfc), maxval(tdsfc)
      print *, 'min/max thetasfc  = ', minval(thetasfc), maxval(thetasfc)
      print *, 'min/max thetaesfc = ', minval(thetaesfc), maxval(thetaesfc)
      print *, 'min/max redp      = ', minval(redp), maxval(redp)
      print *, 'min/max pmsl      = ', minval(pmsl), maxval(pmsl)
      print *, 'min/max ztw0      = ', minval(ztw0), maxval(ztw0)
      print *, 'min/max ztw1      = ', minval(ztw1), maxval(ztw1)
      print *, 'min/max pblhgt    = ', minval(pblhgt), maxval(pblhgt)
      print *, 'min/max cldbase   = ', minval(cldbase), maxval(cldbase)
      print *, 'min/max cldtop    = ', minval(cldtop), maxval(cldtop)
      print *, 'min/max ceiling   = ', minval(ceiling), maxval(ceiling)
      print *, 'min/max cldamt    = ', minval(cldamt), maxval(cldamt)
      print *, 'min/max intliqwat = ', minval(intliqwater), maxval(intliqwater)
      print *, 'min/max intcldice = ', minval(intcldice), maxval(intcldice)
      print *, 'min/max totpcpwat = ', minval(totpcpwater), maxval(totpcpwater)
      print *, 'min/max intrain   = ', minval(intrain), maxval(intrain)
      print *, 'min/max intsnow   = ', minval(intsnow), maxval(intsnow)
      print *, 'min/max intgraupel= ', minval(intgraupel), maxval(intgraupel)
      print *, 'min/max simvis    = ', minval(simvis), maxval(simvis)
      print *, 'min/max max refl  = ', minval(max_refl), maxval(max_refl)
      print *, 'min/max echo tops = ', minval(echo_tops), maxval(echo_tops)
      print *, 'min/max refl sfc  = ', minval(refl_sfc), maxval(refl_sfc)
      print *, 'min/max pcptypsfc = ', minval(pcptype_sfc), maxval(pcptype_sfc)
      print *, 'min/max pcp inc   = ', minval(pcp_inc), maxval(pcp_inc)
      print *, 'min/max snow inc  = ', minval(snow_inc), maxval(snow_inc)
      print *, 'min/max pcp tot   = ', minval(pcp_tot), maxval(pcp_tot)
      print *, 'min/max snow tot  = ', minval(snow_tot), maxval(snow_tot)
      print *, 'min/max helicity  = ', minval(srhel), maxval(srhel)
      print *, 'min/max cape      = ', minval(cape), maxval(cape)
      print *, 'min/max cin       = ', minval(cin), maxval(cin)
      print *, 'min/max li        = ', minval(liftedind), maxval(liftedind)
      print *, 'min/max vis       = ', minval(visibility), maxval(visibility)
      print *, 'min/max heat ind  = ', minval(heatind), maxval(heatind)
      print *, 'min/max lw out    = ', minval(lwout), maxval(lwout)
      print *, 'min/max sw out    = ', minval(swout), maxval(swout)
      print *, 'min/max lw down   = ', minval(lwdown), maxval(lwdown)
      print *, 'min/max sw down   = ', minval(swdown), maxval(swdown)
      print *, 'min/max sh flux   = ', minval(shflux), maxval(shflux)
      print *, 'min/max lh flux   = ', minval(lhflux), maxval(lhflux)
      if (make_firewx) then
         print *, 'min/max upbl      = ', minval(upbl), maxval(upbl)
         print *, 'min/max vpbl      = ', minval(vpbl), maxval(vpbl)
         print *, 'min/max vent indx = ', minval(vnt_index), maxval(vnt_index)
         print *, 'min/max md haines = ', minval(ham_index), maxval(ham_index)
         print *, 'min/max up haines = ', minval(hah_index), maxval(hah_index)
         print *, 'min/max fosberg   = ', minval(fwi_index), maxval(fwi_index)
         print *, 'min/max laps/kelsch = ', minval(fwx_index), maxval(fwx_index)
      end if
   end if

   i4_elapsed = ishow_timer()

!beka moisture flux

   call get_grid_spacing_array(lat, lon, lx, ly, dx, dy)

!       print *,'min/max zsig =',minval(hzsig),maxval(hzsig)
!       write(6,*)'lx,ly,lz',lx,ly,lz

!       do i = 1,lx
!       do j = 1,ly
!           if(hzsig(i,j,1) .gt. hzsig(i,j,lz))then
!               write(6,*)'hzsig qc issue ',i,j
!           endif
!       enddo ! j
!       enddo ! i

   if (.not. large_ngrid) then

      call up_mflux(lx, ly, nz, zsfc, ldf, dx, dy &
                    , husig, hvsig, totpcpwater, upflux &
                    , hzsig, rmsg)

!       write(*,*)'upflux/totpcpwater fields'

!       write(*,*)upflux,totpcpwater

      write (*, *) 'beka beka beka'

   end if

   if (fcsttime .eq. 0) then ! obtain initial time bt11u field from s8a analysis lcv files
      ext = 'lcv'
      var_2d = 's8a'
      call get_laps_2d(laps_reftime, ext, var_2d, units_2d, comment_2d, lx, ly, bt11u, istatus)
      if (istatus .eq. 1) then
         write (6, *) ' obtained initial time bt11u field from s8a analysis lcv file'
      else
         write (6, *) ' warning: could not find initial time swdown field from s8a analysis lcv file'
      end if

   else ! calculate brightness temperature
      stefan_boltzmann = 5.67e-8
      if (.false.) then ! first method
         eff_emissivity = 0.6
         bt11u(:, :) = ((lwout(:, :)/eff_emissivity)/stefan_boltzmann)**0.25
      else          ! second method: adapted from ohring, george, arnold gruber, robert ellingson, 1984:
         ! satellite determinations of the relationship between total longwave radiation flux
         ! and infrared window radiance. j. climate appl. meteor., 23, 416-425.
         ! based on a look at graphs on this and another paper a quadratic fit is being used
         ! using these three points: 200,200 240,255 and 280,320
         do j = 1, ly
         do i = 1, lx
            bt_flux_equiv = (lwout(i, j)/stefan_boltzmann)**0.25
            bt11u(i, j) = 255.+1.5*(bt_flux_equiv - 240.) + .003125*(bt_flux_equiv - 240.)**2
         end do ! i
         end do ! j
      end if

   end if

   return
end

!===============================================================================

subroutine fill_precip(pcp_init, pcp_06, snow_init)

   use lfmgrid

   implicit none

   integer :: nc, itry
   real, dimension(lx, ly) :: pcp_init, pcp_06, snow_init
   character(len=10) :: atime
   logical :: back, there

! read previous accumulated precip data from an intermediate file that
!  is created in the model run directory.

   back = .true.
   nc = index(filename, '/', back)

   if (fcsttime > 0) then
      if (nc > 0) filename0 = filename(1:nc)//'/'//domnum_fstr
      write (atime, '(i6.6)') max(0, fcsttime - precip_dt)
      filename0 = trim(filename0)//'_'//trim(adjustl(atime))//'.pcp'
      itry = 0

80    inquire (file=trim(filename0), exist=there)
      if (there) then
90    if (verbose) then
         print *, ' '
         print *, 'reading previous precip data from: ', trim(filename0)
      end if
      open (1, file=trim(filename0), status='old', form='unformatted')
      if (trim(mtype) == 'nmm') then
         read (1, end=101, err=101) pcp_init, pcp_06, snow_init
         print *, 'reading intermediate file for model: ', trim(mtype)
      else
         read (1, end=101, err=101) pcp_init, snow_init
         pcp_06 = 0.
      end if
      goto 102

101   print *, '  error reading previous precip'

      itry = itry + 1
      if (itry .le. 5) then
         print *, '  trying once again after waiting'
         call sleep(60)
         goto 90
      end if

      print *, '  prior precip set to zero.'
      pcp_init = 0.
      snow_init = 0.

102   close (1)

      else ! not there
      print *, 'could not find previous precip file: ', trim(filename0)

      itry = itry + 1
      if (itry .le. 5) then
         print *, '  trying once again after waiting'
         call sleep(60)
         goto 80
      end if

      print *, '  prior precip set to zero.'
      pcp_init = 0.
      snow_init = 0.
      end if

   else ! fcsttime = 0
      pcp_init = 0.
      snow_init = 0.
      pcp_inc = 0.
      snow_inc = 0.
      pcp_tot = 0.
      snow_tot = 0.

   end if

! create new intermediate filename for current accumulated precip data.

   if (nc > 0) filename0 = filename(1:nc)//'/'//domnum_fstr
   write (atime, '(i6.6)') fcsttime
   filename0 = trim(filename0)//'_'//trim(adjustl(atime))//'.pcp'

   return
end

!===============================================================================

subroutine model_pblhgt(theta, thsfc, psig, zsig, topo, nx, ny, nz, pblhgt)

!  subroutine to estimate height agl in meters of pbl from native
!  coordinate model data.  adapted from the laps routine for
!  terrain-following model coordinates.

   implicit none

   integer :: nx, ny, nz, i, j, k, ktop
   real :: thresh_k, topwgt, botwgt
   real, dimension(nx, ny, nz) :: theta(nx, ny, nz), psig, zsig
   real, dimension(nx, ny) :: thsfc, topo, pblhgt
   logical :: found_pbl_top

   do j = 1, ny
   do i = 1, nx

! compute threshold value that theta needs to exceed
! to be above pbl.  we use surface theta plus an
! additional 3k for slop.

      thresh_k = thsfc(i, j) + 3.0

! now begin at the bottom and work our way up until we find
! the first level with a theta exceeding the threshold.

      found_pbl_top = .false.
      do k = 1, nz

         if (theta(i, j, k) >= thresh_k) then
            ktop = k
            found_pbl_top = .true.
            exit
         end if
      end do

! if we did not find a good pbl, set pbl to first level
! and print out some diagnostics.

      if (.not. found_pbl_top) then
         print *, 'pbl height not found at i/j = ', i, j
         print *, 'surface theta = ', thsfc(i, j)
         print *, 'theta in the column:'
         print *, 'pressure height  theta'
         print *, '-------- ------- --------'
         do k = 1, nz
            print '(f8.0,f8.0,f8.2)', psig(i, j, k), zsig(i, j, k), theta(i, j, k)
         end do
         ktop = 1
         pblhgt(i, j) = zsig(i, j, 1) - topo(i, j)
      else

! we found the top k-level bounding the pbl so interpolate
! to the actual level.

         if (ktop == 1) then
            pblhgt(i, j) = zsig(i, j, 1) - topo(i, j)
         else
! interpolate to get height at thresh_k.
            botwgt = ((theta(i, j, ktop) - thresh_k) &
                      /(theta(i, j, ktop) - theta(i, j, ktop - 1)))
            topwgt = 1.0 - botwgt
            pblhgt(i, j) = botwgt*zsig(i, j, ktop - 1) &
                           + topwgt*zsig(i, j, ktop) - topo(i, j)
         end if
      end if
   end do
   end do

   return
end

!===============================================================================

subroutine lfmclouds(nx, ny, nz, grid_spacing, cldliqmr, cldicemr, snowmr &
                     , height, topo, cldbase, cldtop, ceiling, cldamt)

! this routine is used to compute cloud ceiling (agl), cloud base
! height (msl), cloud top height (msl), and coverage (fraction) given
! mixing ratios of the various microphysical species.

! adapted from the usaf weather agency mm5 post processer.
!  brent shaw, noaa forecast systems lab, dec 2000

   implicit none

   integer :: nx, ny, nz, i, j, k
   real :: grid_spacing, icethresh, liqthresh, snowthresh
   real, dimension(nx, ny) :: topo, cldbase, cldtop, ceiling, cldamt
   real, dimension(nx, ny, nz) :: cldliqmr, cldicemr, snowmr, height

! set thresholds for cloud determination based on grid
! resolution.  kind of hokey, but will get us by for now.

   if (grid_spacing <= 10000) then
      icethresh = 0.000005
      snowthresh = 0.000003
      liqthresh = 0.000003
   else
      icethresh = 0.000005
      snowthresh = 0.000025
      liqthresh = 0.000025
   end if

! loop through using these thresholds to determine cloudiness.

   do j = 1, ny
   do i = 1, nx
      do k = 1, nz
         if ((cldliqmr(i, j, k) >= liqthresh) .or. &
             (cldicemr(i, j, k) >= icethresh) .or. &
             (snowmr(i, j, k) >= snowthresh)) then
            cldbase(i, j) = height(i, j, k)
            ceiling(i, j) = height(i, j, k) - topo(i, j)

! for now all we can do is use cldamt as a yes/no.
! we should look at coming up with a fraction function.

            cldamt(i, j) = 1.0
            exit
         end if
      end do
      do k = nz, 1, -1
         if ((cldliqmr(i, j, k) >= liqthresh) .or. &
             (cldicemr(i, j, k) >= icethresh) .or. &
             (snowmr(i, j, k) >= snowthresh)) then
            cldtop(i, j) = height(i, j, k)
            exit
         end if
      end do
   end do
   end do

   return
end

!===============================================================================

subroutine lfm_integrated_liquid(nx, ny, nz, cond_mr, cice_mr, vapor_mr, rain_mr, snow_mr, graupel_mr, rho, height, pressure, topo &
                                 , intliqwater, intliq700, intcldice, totpcpwater, intrain, intsnow, intgraupel, rmsg)

! computes integrated liquid water and total precip. water in a column.
!  adapted from usaf weather agency mm5 post processor
!  brent shaw, noaa forecast systems lab, dec 2000

   implicit none

   integer :: nx, ny, nz, i, j, k
   real :: height_top, height_bot, dz, dz700, rmsg
   real, dimension(nx, ny) :: newtopo, intliqwater, totpcpwater, intcldice  ! m
   real, dimension(nx, ny) :: topo, intrain, intsnow, intgraupel         ! m
   real, dimension(nx, ny) :: intliq700                               ! m
   real, dimension(nx, ny, nz) :: cond_mr, vapor_mr, cice_mr             ! kg/kg
   real, dimension(nx, ny, nz) :: rain_mr, snow_mr, graupel_mr           ! kg/kg
   real, dimension(nx, ny, nz) :: rho                                  ! kg/m**3
   real, dimension(nx, ny, nz) :: height                               ! m
   real, dimension(nx, ny, nz) :: pressure                             ! pa
   real rho_h2o                                                      ! kg/m**3

   rho_h2o = 1000.                                                   ! kg/m**3

   do j = 1, ny
   do i = 1, nx
      intliqwater(i, j) = 0.0
      intliq700(i, j) = 0.0
      totpcpwater(i, j) = 0.0
      intcldice(i, j) = 0.0
      intrain(i, j) = 0.0
      intsnow(i, j) = 0.0
      intgraupel(i, j) = 0.0
      do k = 1, nz
! compute layer thickness
         if (k == 1) then
            height_bot = topo(i, j)
            height_top = 0.5*(height(i, j, 1) + height(i, j, 2))
         elseif (k == nz) then
            height_bot = 0.5*(height(i, j, nz - 1) + height(i, j, nz))
            height_top = 2*height(i, j, nz) - height_bot
         else
            height_bot = 0.5*(height(i, j, k - 1) + height(i, j, k))
            height_top = 0.5*(height(i, j, k) + height(i, j, k + 1))
         end if

         dz = height_top - height_bot
         if (pressure(i, j, k) .lt. 70000.) then
            dz700 = dz
         else
            dz700 = 0.
         end if

         if (cond_mr(i, j, k) < rmsg) intliqwater(i, j) = intliqwater(i, j) + cond_mr(i, j, k)*(rho(i, j, k)/rho_h2o)*dz  ! meters
         if (cond_mr(i, j, k) < rmsg) intliq700(i, j) = intliq700(i, j) + cond_mr(i, j, k)*(rho(i, j, k)/rho_h2o)*dz700 ! meters
         if (cice_mr(i, j, k) < rmsg) intcldice(i, j) = intcldice(i, j) + cice_mr(i, j, k)*(rho(i, j, k)/rho_h2o)*dz  ! meters
         if (rain_mr(i, j, k) < rmsg) intrain(i, j) = intrain(i, j) + rain_mr(i, j, k)*(rho(i, j, k)/rho_h2o)*dz  ! meters
         if (snow_mr(i, j, k) < rmsg) intsnow(i, j) = intsnow(i, j) + snow_mr(i, j, k)*(rho(i, j, k)/rho_h2o)*dz  ! meters
         if (graupel_mr(i, j, k) < rmsg) intgraupel(i, j) = intgraupel(i, j) + graupel_mr(i, j, k)*(rho(i, j, k)/rho_h2o)*dz  ! meters
         totpcpwater(i, j) = totpcpwater(i, j) + vapor_mr(i, j, k)*(rho(i, j, k)/rho_h2o)*dz  ! meters
      end do
   end do
   end do

   return
end

!===============================================================================

subroutine lfm_reflectivity(lx, ly, nz, rho, hgt, rainmr, icemr, snowmr, graupelmr &
                            , refl, zdr, ldr, tmp, p, qv, max_refl, echo_tops)

   use lfmgrid, only: c_m2z, large_ngrid
   use src_versuch, only: versuch

! subroutine to compute estimated radar reflectivity from
! the precipitation mixing ratios.  will also return
! column max reflectivity and echo tops.  the estimation
! is done using formulae from kessler (1969) and
! rogers and yau (1989).

! adapted from usaf weather agency routine.
!  brent shaw, noaa forecast system lab, dec 2000

   implicit none

   integer :: lx, ly, nz, i, j, k
   integer :: in0r, in0s, in0g, iliqskin
   real, parameter :: svnfrth = 7.0/4.0, max_top_thresh = 5.0
   real :: w
   real, dimension(lx, ly) :: max_refl, echo_tops
   real, dimension(lx, ly, nz) :: rho, rainmr, icemr, snowmr, graupelmr, refl
   real, dimension(lx, ly, *)  :: hgt
   real, dimension(lx, ly, nz) :: zdr, ldr
   real, dimension(lx, ly, nz) :: tmp, p, qv

! local arrays

   real, dimension(lx, ly, nz) :: zhhx, ahhx, zhvx, zvhx, zvvx    ! versuch outputs
   real, dimension(lx, ly, nz) :: delahvx, kkdp, rhohvx         ! versuch outputs

   real :: zvvxmax, zhhxmax

   max_refl = 0.0
   echo_tops = 1.0e37
   if (c_m2z .ne. 'wrf') then
      refl = 0.0
   else
      refl = max(refl, -10.)
      where (refl .lt. 0. .and. refl .gt. -10.) refl = 0.
   end if
   zdr = 0.0
   ldr = 0.0

   print *, 'c_m2z = ', c_m2z

   if (c_m2z == 'rip') then
      in0r = 0
      in0s = 0
      in0g = 0
      iliqskin = 0
!call dbzcalc(qv,rainmr,snowmr,graupelmr,tmp,p,refl, &
!             ly,lx,nz,in0r,in0s,in0g,iliqskin)

      do j = 1, ly
      do i = 1, lx
!  compute the max value in the column
         max_refl(i, j) = maxval(refl(i, j, :))
      end do ! i
      end do ! j

   else
      do j = 1, ly
      do i = 1, lx
         do k = 1, nz

! compute the basic reflectivity using rams reflectivity algorithm.

            if (c_m2z == 'rams') then

               w = 264083.11*(rainmr(i, j, k) &
                              + 0.2*(snowmr(i, j, k)) &
                              + 2.0*graupelmr(i, j, k))
               w = max(1., w)
               refl(i, j, k) = 17.8*alog10(w)

               if (refl(i, j, k) .eq. 0.) then
                  refl(i, j, k) = -10.0
               end if

! compute the basic reflectivity using kessler reflectivity algorithm.

            elseif (c_m2z == 'kessler') then

               refl(i, j, k) = 17300.0* &
                               (rho(i, j, k)*1000.0* &
                                max(0.0, rainmr(i, j, k)))**svnfrth

               ! add the ice component
               refl(i, j, k) = refl(i, j, k) + &
                               38000.0*(rho(i, j, k)*1000.0* &
                                        !                   max(0.0,icemr(i,j,k)+snowmr(i,j,k)+graupelmr(i,j,k)))**2.2
                                        max(0.0, snowmr(i, j, k) + graupelmr(i, j, k)))**2.2

               ! convert to dbz (with a lower limit of zero, given the max application)
               refl(i, j, k) = 10.*alog10(max(refl(i, j, k), 1.0))

               ! set zero reflectivity to -10
               if (refl(i, j, k) .eq. 0.) then
                  refl(i, j, k) = -10.0
               end if

! computing the reflectivity by using the synpolrad code

            elseif (c_m2z == 'synp') then

               call versuch(rainmr(i, j, k), snowmr(i, j, k), graupelmr(i, j, k), tmp(i, j, k) &   ! i
                            , zhhx(i, j, k) &   ! o
                            , ahhx(i, j, k) &   ! o
                            , zhvx(i, j, k), zvhx(i, j, k), zvvx(i, j, k) &   ! o
                            , p(i, j, k) &   ! i
                            , qv(i, j, k) &   ! i
                            , delahvx(i, j, k), kkdp(i, j, k), rhohvx(i, j, k))                   ! o

!       qr = specific rain   (kg/kg)
!       qs = specific snow (kg/kg)
!       qg = specific graupel (kg/kg)
!       t  = temperature (k)
!       p  = pressure (pa) - druck
!       qv = specific water vapor (kg/kg) - qd

!       intent(out):
!       zhh = reflectivity-hh (z), not in dbz
!       zhv
!       zvh
!       zvv
!       kappa = extinction coeff (km-1) - ahhx
!       delahv   = specific differential attenuationa (db km-1)
!       kkdp   = specific differential phase (degree km-1)
!       rhohv = hv-vh cross correlation coeffcient (n/a)

!       convert to dbz
               zvvxmax = max(zvvx(i, j, k), 1.0)
               zhhxmax = max(zhhx(i, j, k), 1.0)
               refl(i, j, k) = 10.*alog10(zhhxmax)
               zdr(i, j, k) = 10.*alog10(zhhx(i, j, k)/zvvxmax)
               ldr(i, j, k) = 10.*alog10(zvhx(i, j, k)/zhhxmax)

            end if ! c_m2z (reflectivity algorithm)

! since we are going from the ground up, we can
! check threshold and set echo top.

            if (refl(i, j, k) >= max_top_thresh) echo_tops(i, j) = hgt(i, j, k)
         end do

! compute the max value in the column

         max_refl(i, j) = maxval(refl(i, j, :))
      end do
      end do

   end if

! convert microphysical precipitation mixing ratios to concentrations - multiplying by rho
   rainmr(:, :, :) = rainmr(:, :, :)*rho(:, :, :)
   snowmr(:, :, :) = snowmr(:, :, :)*rho(:, :, :)
! icemr(:,:,:)     = icemr(:,:,:) * rho(:,:,:)
   graupelmr(:, :, :) = graupelmr(:, :, :)*rho(:, :, :)

   return
end

!===============================================================================

subroutine lfm_ua_pcptype(nx, ny, nz, rainmr_prs, snowmr_prs, graupelmr_prs &
                          , pcptype_prs)

   implicit none

   real, parameter :: nonecode = 0.
   real, parameter :: raincode = 1.
   real, parameter :: snowcode = 2.
   real, parameter :: zraincode = 3.
   real, parameter :: sleetcode = 4.
   real, parameter :: hailcode = 5.
   real, parameter :: drizzlecode = 6.
   real, parameter :: zdrizzlecode = 7.
   real, parameter :: rainsnowcode = 8.
   real, parameter :: rainicecode = 9.

   integer :: nx, ny, nz, i, j, k
   real, dimension(nx, ny, nz) :: rainmr_prs, snowmr_prs, graupelmr_prs, pcptype_prs

   pcptype_prs = nonecode

   do k = 1, nz
   do j = 1, ny
   do i = 1, nx
      if (rainmr_prs(i, j, k) > 0.) then
         if (snowmr_prs(i, j, k) > 0.) then
            if (graupelmr_prs(i, j, k) > snowmr_prs(i, j, k)) then
               pcptype_prs(i, j, k) = rainicecode
            else
               pcptype_prs(i, j, k) = rainsnowcode
            end if
         else
            if (graupelmr_prs(i, j, k) > 0.) then
               pcptype_prs(i, j, k) = rainicecode
            else
               pcptype_prs(i, j, k) = raincode
            end if
         end if
      else
         if (snowmr_prs(i, j, k) > 0.) then
            if (graupelmr_prs(i, j, k) > snowmr_prs(i, j, k)) then
               pcptype_prs(i, j, k) = sleetcode
            else
               pcptype_prs(i, j, k) = snowcode
            end if
         else
            if (graupelmr_prs(i, j, k) > 0) then
               pcptype_prs(i, j, k) = sleetcode
            else
               pcptype_prs(i, j, k) = nonecode
            end if
         end if
      end if
   end do
   end do
   end do

   return
end

!===============================================================================

subroutine lfm_sfc_pcptype(nx, ny, nz, rainmr_sig, snowmr_sig, graupelmr_sig &
                           , pcptype_sfc)

   implicit none

   real, parameter :: nonecode = 0.
   real, parameter :: raincode = 1.
   real, parameter :: snowcode = 2.
   real, parameter :: zraincode = 3.
   real, parameter :: sleetcode = 4.
   real, parameter :: hailcode = 5.
   real, parameter :: drizzlecode = 6.
   real, parameter :: zdrizzlecode = 7.
   real, parameter :: rainsnowcode = 8.
   real, parameter :: rainicecode = 9.

   integer :: nx, ny, nz, i, j, k
   real, dimension(nx, ny) :: pcptype_sfc
   real, dimension(nx, ny, nz) :: rainmr_sig, snowmr_sig, graupelmr_sig

   pcptype_sfc = nonecode

   do j = 1, ny
   do i = 1, nx
      if (rainmr_sig(i, j, 1) > 0.) then
         if (snowmr_sig(i, j, 1) > 0.) then
            if (graupelmr_sig(i, j, 1) > snowmr_sig(i, j, 1)) then
               pcptype_sfc(i, j) = rainicecode
            else
               pcptype_sfc(i, j) = rainsnowcode
            end if
         else
            if (graupelmr_sig(i, j, 1) > 0.) then
               pcptype_sfc(i, j) = rainicecode
            else
               pcptype_sfc(i, j) = raincode
            end if
         end if
      else
         if (snowmr_sig(i, j, 1) > 0.) then
            if (graupelmr_sig(i, j, 1) > snowmr_sig(i, j, 1)) then
               pcptype_sfc(i, j) = sleetcode
            else
               pcptype_sfc(i, j) = snowcode
            end if
         else
            if (graupelmr_sig(i, j, 1) > 0.) then
               pcptype_sfc(i, j) = sleetcode
            else
               pcptype_sfc(i, j) = nonecode
            end if
         end if
      end if
   end do
   end do

   return
end

!===============================================================================

subroutine interp_press_to_z(press_levels, heights, new_level, press_out &
                             , nx, ny, nz)

! given a 1d array of pressure levels and a 3d array of heights (m) at
! those levels, this routine interpolates the pressure to the desired
! new_level.
! pressures are in pa, heights are in m!

   implicit none

   integer :: nx, ny, nz, i, j, k

   real :: press_levels(nz), heights(nx, ny, nz), new_level, press_out(nx, ny) &
           , ptop, pbot, ztop, zbot

   do j = 1, ny
   do i = 1, nx
   do k = 2, nz
      if (heights(i, j, k) >= new_level) then
         ztop = heights(i, j, k)
         zbot = heights(i, j, k - 1)
         ptop = press_levels(k)*0.01    ! convert to mb
         pbot = press_levels(k - 1)*0.01
         press_out(i, j) = exp((new_level*alog(pbot/ptop) &
                                - ztop*alog(pbot) &
                                + zbot*alog(ptop)) &
                               /(zbot - ztop))*100.  ! convert back to pa
         exit
      end if
   end do
   end do
   end do

   return
end

!===============================================================================

subroutine helicity(usig, vsig, zsig, usf, vsf, terrain, imax, jmax, ksig, srelhel)

! name: helicity
!
! purpose
! =======
! calculate storm relative helicity.
!
! remarks
! =======
! helicity is equal to the vertical integral of:
!
! (v - c) x dv/dz
!
! where v is the model predicted wind and c is the
! estimated storm motion vector.
!
! updates
! =======
! ??? 97   initial version......................................dnxm
! dec 98   modified call to wdir function.  an error in wdir was
!          corrected............................................dnxm
! jan 99   changed variable names usigcrs and vsigcrs to usig and
!          and vsig as winds are now at the cross point.........dnxm
! jan 01   adapted for use at fsl by b. shaw, noaa/fsl

   use constants

   implicit none

   integer :: i, imax, j, jmax, k, k3km, k10km, ksig
   real :: dz, stormdir, stormspd, stormu, stormv, sumdz, sumhel, sumu, sumv, wdir, wspd, dudz, dvdz, udif, vdif, zdif
   real, dimension(imax, jmax) :: srelhel, terrain, usf, vsf
   real, dimension(imax, jmax, ksig) :: usig, vsig, zsig

! determine indices of 3 and 10 km layers (above ground level).

   do j = 1, jmax
   do i = 1, imax
      do k = 1, ksig

         if (((zsig(i, j, k) - terrain(i, j)) >= 3000.)) then
            k3km = k
            exit
         end if
      end do

      do k = ksig, k3km, -1
         if (((zsig(i, j, k) - terrain(i, j)) <= 10000.)) then
            k10km = k
            exit
         end if
      end do

! estimate storm motion vector.  speed is 75 percent of the mean
!   wind between 3 and 10 km.  direction is 30 degrees to the
!   right of the mean wind.

      sumu = 0.
      sumv = 0.
      sumdz = 0.

      do k = k3km, k10km
         dz = zsig(i, j, k) - zsig(i, j, k - 1)
         sumdz = sumdz + dz
         sumu = sumu + (0.5*dz*(usig(i, j, k) + usig(i, j, k - 1)))
         sumv = sumv + (0.5*dz*(vsig(i, j, k) + vsig(i, j, k - 1)))
      end do
      stormu = sumu/sumdz
      stormv = sumv/sumdz
      stormspd = wspd(stormu, stormv)*0.75

! when calling wdir, send in a cone factor of zero
!   so that direction is grid relative.

      stormdir = wdir(stormu, stormv, 0., 0., 0.) + 30.
      if (stormdir > 360.) stormdir = stormdir - 360.

      stormu = -stormspd*sin(stormdir*deg2rad)
      stormv = -stormspd*cos(stormdir*deg2rad)

      sumhel = 0.

! calculate helicity.  integrate between the ground and 3 km,
!   a depth that is frequently used in the literature.

      do k = 1, k3km
         if (k .eq. 1) then
            udif = 0.5*(usig(i, j, k) + usf(i, j)) - stormu
            vdif = 0.5*(vsig(i, j, k) + vsf(i, j)) - stormv
            dz = zsig(i, j, k) - terrain(i, j)
            dudz = (usig(i, j, k) - usf(i, j))/dz
            dvdz = (vsig(i, j, k) - vsf(i, j))/dz
         else
            udif = 0.5*(usig(i, j, k) + usig(i, j, k - 1)) - stormu
            vdif = 0.5*(vsig(i, j, k) + vsig(i, j, k - 1)) - stormv
            dz = zsig(i, j, k) - zsig(i, j, k - 1)
            dudz = (usig(i, j, k) - usig(i, j, k - 1))/dz
            dvdz = (vsig(i, j, k) - vsig(i, j, k - 1))/dz
         end if
         sumhel = sumhel + (udif*dvdz - vdif*dudz)*dz
      end do
      srelhel(i, j) = -sumhel
   end do
   end do

   return
end

!===============================================================================

subroutine updraft_helicity(usig, vsig, wsig, zsig, zsigf, terrain, lat, lon, imax, jmax, ksig, uhel)

! updraft_helicity
!
! calculate the helicity of updrafts within the 2-5km agl layer of the updraft
!
! updraft helicity is equal to the vertical integral bounded by 2-5km of
! w x zeta
! where w is the vertical velocity and zeta is the relative vorticity
! a centered-difference scheme is used.
!
! updates
! =======
! 09 07   initial version.........................................cja

   use constants

   implicit none

   integer :: i, imax, j, jmax, k, k2km, k5km, ksig
   real :: dz, dx, dy, dvdx, dudy, sumhel
   real, dimension(imax, jmax) :: uhel, terrain, lat, lon
   real, dimension(imax, jmax, ksig) :: wsig, usig, vsig, zsig
   real, dimension(imax, jmax, ksig + 1) :: zsigf

! determine indices of 2 and 5 km layers (above ground level).

   do j = 2, jmax - 1
   do i = 2, imax - 1
      do k = 1, ksig

         if (((zsig(i, j, k) - terrain(i, j)) >= 2000.)) then
            k2km = k
            exit
         end if
      end do

      do k = ksig, k2km, -1
         if (((zsig(i, j, k) - terrain(i, j)) <= 5000.)) then
            k5km = k
            exit
         end if
      end do

! calculate updraft helicity.  integrate between 2 and 5 km.

      sumhel = 0.
      do k = k2km, k5km
         dz = zsigf(i, j, k) - zsigf(i, j, k - 1)
         dx = (lon(i + 1, j) - lon(i - 1, j))*cos(lat(i, j)*deg2rad)*eradius*deg2rad
         dvdx = (vsig(i + 1, j, k) - vsig(i - 1, j, k))/dx
         dy = (lat(i, j + 1) - lat(i, j - 1))*eradius*deg2rad
         dudy = (usig(i, j + 1, k) - usig(i, j - 1, k))/dy
         sumhel = sumhel + wsig(i, j, k)*(dvdx - dudy)*dz
         if (k == k5km) then
!      print *, i,j,k,dx,dy
!      print *, dz,dvdx,dudy,wsig(i,j,k),sumhel
         end if
      end do
      uhel(i, j) = sumhel
   end do
   end do

   return
end

!===============================================================================

subroutine bulk_shear(usig, vsig, zsig, zsigf, terrain, lat, lon, imax, jmax, ksig, bshr)

! bulk_shear - sfc to 6km agl

   use constants

   implicit none

   integer :: i, imax, j, jmax, k, k0km, k6km, ksig
   real    :: usigl, usigh, vsigl, vsigh, ushr, vshr
   real, dimension(imax, jmax) :: bshr, terrain, lat, lon
   real, dimension(imax, jmax, ksig) :: usig, vsig, zsig
   real, dimension(imax, jmax, ksig + 1) :: zsigf

! determine indices of 0 and 6 km layers (above ground level).

   do j = 1, jmax
   do i = 1, imax
      do k = 1, ksig
         if (((zsig(i, j, k) - terrain(i, j)) >= 0.)) then
            k0km = k
            exit
         end if
      end do

      do k = ksig, k0km, -1
         if (((zsig(i, j, k) - terrain(i, j)) <= 6000.)) then
            k6km = k
            exit
         end if
      end do

! calculate bulk shear.  difference between 0 and 6 km.

      usigl = usig(i, j, k0km)
      usigh = usig(i, j, k6km)
      vsigl = vsig(i, j, k0km)
      vsigh = vsig(i, j, k6km)

      ushr = usigh - usigl
      vshr = vsigh - vsigl

      bshr(i, j) = sqrt(ushr**2 + vshr**2)

   end do
   end do

   return
end

!===============================================================================

subroutine ventilation(usig, vsig, zsig, pblhgt, topo, nx, ny, nz &
                       , upbl, vpbl, vent_ind)

   implicit none

   integer :: nx, ny, nz, i, j, k, nbl
   real :: usum, vsum, umean, vmean, spmean
   real, dimension(nx, ny) :: pblhgt, topo, upbl, vpbl, vent_ind
   real, dimension(nx, ny, nz) :: usig, vsig, zsig

   do j = 1, ny
   do i = 1, nx
      if (pblhgt(i, j) > 0.) then

! compute mean wind within the pbl.

         nbl = 0
         usum = 0.
         vsum = 0.
         umean = 0.
         vmean = 0.
         do k = 1, nz
            if (zsig(i, j, k) - topo(i, j) <= pblhgt(i, j)) then
               nbl = nbl + 1
               usum = usum + usig(i, j, k)
               vsum = vsum + vsig(i, j, k)
            else
               exit
            end if
         end do
         if (nbl > 0) then
            umean = usum/float(nbl)
            vmean = vsum/float(nbl)

! compute mean wind speed for this layer.

            spmean = sqrt(umean**2 + vmean**2)

! multiply mean speed by pbl depth to get index.

            vent_ind(i, j) = pblhgt(i, j)*spmean
            upbl(i, j) = umean
            vpbl(i, j) = vmean
         else

! pbl height is lower than the lowest model level...use
! lowest model wind.

            spmean = sqrt(usig(i, j, 1)**2 + vsig(i, j, 1)**2)
            vent_ind(i, j) = pblhgt(i, j)*spmean
            upbl(i, j) = usig(i, j, 1)
            vpbl(i, j) = vsig(i, j, 1)
         end if
      else
         if (pblhgt(i, j) < 0.) &
            print *, 'warning:  pbl height <0 in ventilation index:', pblhgt(i, j)
         vent_ind(i, j) = 0.
         upbl(i, j) = 0.
         vpbl(i, j) = 0.
      end if
   end do
   end do

   return
end

!===============================================================================

subroutine wind_agl(usig, vsig, zsig, aglhgt, topo, lx, ly, nz &
                    , uagl, vagl, r_missing_data)

   implicit none

   integer :: lx, ly, nz, i, j, k, nbl
   logical :: l_debug
   real :: aglhgt, refhgt, usum, vsum, umean, vmean, spmean, r_missing_data, zlow, zhigh
   real, dimension(lx, ly) :: topo, uagl, vagl
   real, dimension(lx, ly, nz) :: usig, vsig, zsig

   write (6, *) ' subroutine wind_agl...'

   do j = 1, ly
   do i = 1, lx

! compute mean wind within the pbl.
      if (i .eq. lx/2 .and. j .eq. ly/2) then
         l_debug = .true.
      else
         l_debug = .false.
      end if

      nbl = 0
      usum = 0.
      vsum = 0.
      umean = 0.
      vmean = 0.
      refhgt = topo(i, j) + aglhgt
      do k = 1, nz - 1
         if (zsig(i, j, k) <= refhgt .and. zsig(i, j, k + 1) >= refhgt) then
            zlow = refhgt - zsig(i, j, k)
            zhigh = zsig(i, j, k + 1) - refhgt
            nbl = nbl + 1

            usum = 0.
            usum = usum + usig(i, j, k + 1)*(zlow/(zlow + zhigh))
            usum = usum + usig(i, j, k)*(zhigh/(zlow + zhigh))

            vsum = 0.
            vsum = vsum + vsig(i, j, k + 1)*(zlow/(zlow + zhigh))
            vsum = vsum + vsig(i, j, k)*(zhigh/(zlow + zhigh))

            exit
         end if
      end do

      if (nbl > 0) then
         if (l_debug) then
            write (6, *) 'topo,refhgt,zsig1,zsig2', topo(i, j), refhgt, zsig(i, j, k), zsig(i, j, k + 1)
            write (6, *) 'zlow,zhigh,usig(i,j,k),usig(i,j,kp1),usum = ', zlow, zhigh, usig(i, j, k), usig(i, j, k + 1), usum
         end if
         uagl(i, j) = usum
         vagl(i, j) = vsum
      else
         uagl(i, j) = r_missing_data
         vagl(i, j) = r_missing_data
      end if

   end do
   end do

   write (6, *) aglhgt, ' meter wind u range = ', minval(uagl), maxval(uagl)
   write (6, *) aglhgt, ' meter wind v range = ', minval(vagl), maxval(vagl)

   return
end

!===============================================================================

subroutine haines_layer(p3d, t3d_k, td3d_k, haines2d, nx, ny, nz &
                        , pmbbot, pmbtop)

! computes haines index for layer bounded by top and bottom pressure
!  levels pmbtop and pmbbot (mb).

   implicit none

   integer :: nx, ny, nz, i, j, k, kk, km1
   real :: p3d(nx, ny, nz)     ! 3d pressure in pa
   real :: t3d_k(nx, ny, nz)  ! 3d temp in k
   real :: td3d_k(nx, ny, nz)  ! 3d dewpoint in k
   real :: haines2d(nx, ny)   ! 2d haines index
   real :: pmbbot, pmbtop     ! bounding mb levels
   real :: ppabot, ppatop, factor
   real :: tmkt, tmkb, tdkb, deltat, dpdep, factor1, factor2

   ppabot = pmbbot*100.
   ppatop = pmbtop*100.

   do j = 1, ny
   do i = 1, nx
      if (p3d(i, j, 1) >= ppabot) then
         do k = 1, nz - 1

! find temperature at the top.

            if (p3d(i, j, k) > ppatop .and. p3d(i, j, k + 1) <= ppatop) then
               tmkt = t3d_k(i, j, k) + (t3d_k(i, j, k + 1) - t3d_k(i, j, k)) &
                      *(alog(ppatop/p3d(i, j, k))) &
                      /(alog(p3d(i, j, k + 1)/p3d(i, j, k)))
            end if

! find temp/dewpoint at the bottom of the layer.

            if (p3d(i, j, k) > ppabot .and. p3d(i, j, k + 1) <= ppabot) then
               factor = (alog(ppabot/p3d(i, j, k))) &
                        /(alog(p3d(i, j, k + 1)/p3d(i, j, k)))
               tmkb = t3d_k(i, j, k) + (t3d_k(i, j, k + 1) - t3d_k(i, j, k))*factor
               tdkb = td3d_k(i, j, k) + (td3d_k(i, j, k + 1) - td3d_k(i, j, k))*factor
            end if

         end do

         deltat = tmkb - tmkt
         dpdep = tmkb - tdkb

         if (nint(pmbbot) == 700) then  ! high haines
            if (deltat <= 17.5) then
               factor1 = 1.
            elseif (deltat > 17.5 .and. deltat <= 21.5) then
               factor1 = 2.
            else
               factor1 = 3.
            end if

            if (dpdep <= 14.5) then
               factor2 = 1.
            elseif (dpdep > 14.5 .and. dpdep <= 20.5) then
               factor2 = 2.
            else
               factor2 = 3.
            end if

         elseif (nint(pmbbot) == 850) then  ! mid-level haines
            if (deltat <= 5.5) then
               factor1 = 1.
            elseif (deltat > 5.5 .and. deltat <= 10.5) then
               factor1 = 2.
            else
               factor1 = 3.
            end if

            if (dpdep <= 5.5) then
               factor2 = 1.
            elseif (dpdep > 5.5 .and. dpdep <= 12.5) then
               factor2 = 2.
            else
               factor2 = 3.
            end if

         else
            print *, 'haines_layer needs 850 or 700 as bottom layer'
            print *, 'bottom level (mb) specified:', pmbbot
            stop 'bad_haines_layer'

         end if

         haines2d(i, j) = factor1 + factor2

      end if

   end do
   end do

   return
end

!===============================================================================

subroutine fosberg_fwi(t2k, rh2pct, u10, v10, nx, ny, fwi)

! computes the fosberg fire weather index

   implicit none

   integer :: nx, ny, i, j
   real :: uuu10, vvv10, m, n, t2f, rh2
   real :: t2k(nx, ny)      ! sfc temp (k)
   real :: rh2pct(nx, ny)   ! sfc rh (%)
   real :: u10(nx, ny)      ! 10 m u wind (m/s)
   real :: v10(nx, ny)      ! 10 m v wind (m/s)
   real :: fwi(nx, ny)      ! fosberg index

   do j = 1, ny
   do i = 1, nx

! convert temperature from k to f.

      t2f = 1.8*(t2k(i, j) - 273.15) + 32.0

! convert u/v from m/s to mph.

      uuu10 = u10(i, j)*2.237
      vvv10 = v10(i, j)*2.237
      rh2 = rh2pct(i, j)

      if (rh2 <= 10.5) then
         m = 0.03229 + (0.281073*rh2) - (0.000578*rh2*t2f)
      elseif (rh2 > 10.5 .and. rh2 <= 50.5) then
         m = 2.22749 + (0.160107*rh2) - (0.014784*t2f)
      elseif (rh2 > 50.5 .and. rh2 <= 100.) then
         m = 21.0606 + (0.005565*rh2**2) - (0.00035*rh2*t2f) - (0.483199*rh2)
      else
         m = 21.0606 + (0.005565*100**2) - (0.00035*100*t2f) - (0.483199*100)
      end if

      n = 1.-2.*(m/30.) + 1.5*(m/30.)**2 - 0.5*(m/30.)**3
      fwi(i, j) = (n*sqrt(1.+uuu10**2 + vvv10**2))/0.3002

   end do
   end do

   return
end

!===============================================================================

subroutine capecin(psig, tsig, thetaesig, thetasig &
                   , rhsig, zsigfull, tprs, li, posbuoyen, negbuoyen &
                   , k500, imax, jmax, ksig, kprs)

! purpose
! =======
! calculate positive buoyant energy (or convective available
! energy) and negative buoyant energy (or convective inhibition).
!
! also calculate lifted index.
!
! reference
! =========
! doswell and rasmussen (1994), wea and fcsting, p 625.
!
! updates
! =======
! 4 jan 01  - adapted from usaf weather agency routine
!             b. shaw, noaa/fsl

   use constants
   use lfmgrid, only: large_pgrid

   implicit none

   integer :: i, imax, j, jmax, k, k500, kmax, ksig, kprs
   real, parameter :: pref = 1000.
   real :: prslcl, thw, thwmax, tlcl, tmplcl, tparcel, wobf, deltaz, dtheta, thetaparcel, potential_temp, psave
   real, dimension(imax, jmax) :: li, negbuoyen, posbuoyen
   real, dimension(imax, jmax, ksig) :: psig, rhsig, thetasig, thetaesig, tsig
   real, dimension(imax, jmax, kprs) :: tprs
   real, dimension(imax, jmax, ksig + 1) :: zsigfull
   real, allocatable, dimension(:) :: buoy
   logical :: compute_cin

   posbuoyen = 0.
   negbuoyen = 0.
   allocate (buoy(ksig))
   do j = 1, jmax
   do i = 1, imax
      thwmax = -9999.0

      do k = 1, ksig

! pick the most unstable parcel in the lowest 50 mb as
!   indicated by the sigma level with the highest wet bulb
!   potential temperature.  store index in kmax.

         if (((psig(i, j, 1) - psig(i, j, k)) < 50.)) then
            thw = thetaesig(i, j, k) - wobf(thetaesig(i, j, k) - t0)
            if (thw > thwmax) then
               kmax = k
               thwmax = thw
            end if
         else
            exit
         end if

      end do

! calculate lifted index by lifting the most unstable parcel.

      if (.not. large_pgrid) then
         call the2t(thetaesig(i, j, kmax), 500., tparcel)
         li(i, j) = tprs(i, j, k500) - tparcel
      end if

! calculate the temperature and pressure of the lifting
!   condensation level.

      tmplcl = tlcl(tsig(i, j, kmax), rhsig(i, j, kmax))
      prslcl = psig(i, j, kmax)*(tmplcl/tsig(i, j, kmax))**cpor

! calculate the buoyancy.

      posbuoyen(i, j) = 0.
      negbuoyen(i, j) = 0.
      do k = kmax, ksig

! above the lcl, calculate virtual temperature of the
!   parcel as it moves along a moist adiabat.  below the
!   lcl, lift parcel along a dry adiabat.

         if (psig(i, j, k) <= prslcl) then
            call the2t(thetaesig(i, j, kmax), psig(i, j, k), tparcel)
         else
            tparcel = thetasig(i, j, kmax)/(pref/psig(i, j, k))**kappa
         end if

! compute the potential temperature of the parcel.

         thetaparcel = potential_temp(tparcel, psig(i, j, k)*100.)
         dtheta = thetaparcel - thetasig(i, j, k)
         deltaz = zsigfull(i, j, k + 1) - zsigfull(i, j, k)
         buoy(k) = deltaz*dtheta/thetasig(i, j, k)
      end do

! now loop through the column again, partitioning the buoyency
!   into positive (cape) and negative (cin) component.  we terminate
!   the contribution to cin when/if a layer of cape greater than
!   150 mb deep is found.

      compute_cin = .true.
      psave = -100.
      do k = kmax, ksig
         if (buoy(k) > 0.) then
            if (psave < 0) then
               psave = psig(i, j, k)
            else
               if ((psave - psig(i, j, k)) > 150.) then
                  compute_cin = .false.
               end if
            end if
            posbuoyen(i, j) = posbuoyen(i, j) + buoy(k)
         elseif (buoy(k) < 0.) then
            psave = -100.
            if (compute_cin) then
               negbuoyen(i, j) = negbuoyen(i, j) + buoy(k)
            end if
         end if
      end do
   end do
   end do

   posbuoyen = grav*posbuoyen
   negbuoyen = grav*negbuoyen

! cap the negative buoyancy to a maximum value of 700 j/kg

   where (negbuoyen < -700) negbuoyen = -700.
   deallocate (buoy)

   return
end

!===============================================================================

subroutine wintprec(tsig, zsig, zprs, psfc, tsfc, tdsfc, terrain, prcpinc, imax, jmax &
                    , ksig, kprs, k700, k850, k1000, prcpconv, preciptype)

! (winter) precipitation algorithm
!
! purpose:  this product identifies areas of precipitation and the
! type expected, based on mm5 data.  the process is essentially two-
! fold.  first, areas of precipitation are identified when the mm5
! precipitation array (prcpinc) exceeds 0.01 inch.
!
! second, thickness thresholds are used at the gridpoints for which
! it has been determined that precipitation is occurring and the
! surface pressure is greater than 850mb (i.e., non-mountainous
! regions).  the thickness thresholds utilized are based on
! meteorological research from the pertinent sources
! (e.g., mwr, w&f), and are as follows:
!
!                             (thick1)     (thick2)
!                           1000mb-850mb  850mb-700mb  <--thickness
!
! precipitation type:   rain   gt 1310      gt 1540
!
!              freezing rain   gt 1310      gt 1540 [sig1 t < 0co]
!
!                  ice/mixed   le 1310      gt 1540
!
!                       snow   le 1310      le 1540.
!
! over mountainous terrain, precipitation type is limited to either
! rain or snow.  this is consistent with climatic data presented in
! "a regional climatology of freezing precipitation for the
! contiguous united states" (10th conf. on applied climatology, 20-
! 24 oct 97).  where a precipitation occurrence has been determined,
! the temperatures in the lowest 1500 m are checked.  if all are
! below freezing, snow is forecasted; otherwise rain is predicted.
!
! modification:  added ability to predict regions where thunderstorm
! activity may occur.  prior to exiting the main loop, a check is
! made:  where rain is predicted and the convective component of the
! precip exceeds 0.01", forecast for thunderstorms.
!
! updates
! =======
! jan 2001 1998  initial version, adapted from usaf weather agency
!     brent shaw, noaa/forecast systems lab

   use constants

   implicit none

   integer :: i, imax, j, jmax, k, k1000, k850, k700, kprs, ksig, k1500
   real :: tsfcf, thickhigh, thicklow, tsig1, tsig2, tsig3, fahren, td_c, tw_c, pres_mb, twet_fast
   real, dimension(imax, jmax) :: prcpconv, prcpinc, preciptype, psfc, terrain, tsfc, tdsfc
   real, dimension(imax, jmax, ksig) :: tsig, zsig
   real, dimension(imax, jmax, kprs) :: zprs

   do j = 1, jmax
   do i = 1, imax

! the threshold for calculating precip type is 0.0001 meter
!  per time period.

      if (prcpinc(i, j) <= 0.0001) then
         preciptype(i, j) = 0
      else

! check the surface pressure to determine whether high or low-elevation
!  logic is used.  850 mb is the cutoff.

! low elevation grid point.

         if (psfc(i, j) > 85000.) then

! calculate thicknesses that will be used to determine precip type.

            thicklow = zprs(i, j, k850) - zprs(i, j, k1000)
            thickhigh = zprs(i, j, k700) - zprs(i, j, k850)

! rain, or if surface temperature is below freezing, freezing rain.

            if ((thicklow > 1310.) .and. (thickhigh > 1540.)) then
               if (tsfc(i, j) >= t0) then
                  preciptype(i, j) = 1
               else
                  preciptype(i, j) = 3
               end if

! ice/mixed.

            elseif (thicklow <= 1310. .and. thickhigh > 1540.) then
               preciptype(i, j) = 4

! rain or snow.

            elseif (thicklow <= 1310. .and. thickhigh <= 1540.) then
               tsfcf = fahren(tsfc(i, j) - t0)
               td_c = tdsfc(i, j) - t0
               pres_mb = psfc(i, j)/100.
               tw_c = twet_fast(tsfc(i, j) - t0, td_c, pres_mb)
!           if (tsfcf >= 37.) then
               if (tw_c >= +1.3) then
                  preciptype(i, j) = 1
               else
                  preciptype(i, j) = 5
               end if

! rain.

            else
               preciptype(i, j) = 1
            end if

! high terrain grid point.

         else

! find top of 1500 m agl layer.

            do k = 1, ksig
               if (zsig(i, j, k) - terrain(i, j) >= 1500.) then
                  k1500 = k
                  exit
               end if
            end do

! if the model top is ever lowered, the above code
!  could fail on the top of a high mountain.

            k1500 = max(k1500, 1)

! find temperature at the bottom, top and middle of
!  1500 m agl layer.  for middle layer, recycle variable k1500.

            tsig1 = tsig(i, j, 1)
            tsig3 = tsig(i, j, k1500)
            k1500 = max(1, nint(float(k1500)/2.))
            tsig2 = tsig(i, j, k1500)

! snow.
            td_c = tdsfc(i, j) - t0
            pres_mb = psfc(i, j)/100.
            tw_c = twet_fast(tsfc(i, j) - t0, td_c, pres_mb)

            !       if (tsig1 < t0 .and. tsig2 < t0 .and. tsig3 .lt. t0) then
            if (tw_c <= +1.3) then
               preciptype(i, j) = 5

! rain & check for thunderstorms.

            else
               preciptype(i, j) = 1
            end if
         end if

      end if

      if (preciptype(i, j) == 1 .and. prcpconv(i, j) >= 0.001) &
         preciptype(i, j) = 2         ! thunderstorm

   end do
   end do

   return
end

!===============================================================================

subroutine snowfall(tsig, tsfc, prcpinc, preciptype, imax, jmax, ksig, snowinc, snowtot)

! name: snow accumulation algorithm
!
! purpose
! =======
! this algorithm calculates incremental snow accumulation and
! total snow accumulation.
!
! method
! ======
! - if precip type is snow
!    - calculate column max temperature on dot
!    - calculate incremental precip on dot
!    - obtain snow_to_rain ratio from laps analysis subroutine
! - if precip type is not snow
!    - snow accumulation is zero
!
! variables
! =========
! name             type      i/o     definition
! ----             ----      ---     ----------
! imax             integer  input   grid dimension, i direction
! jmax             integer  input   grid dimension, j direction
! prcpinc          real     input   3hr precip accum
!                                   array(2-d)
! preciptype       integer  input   0 - is no precip
!                                     - is rain
!                                   2 - is snow
!                                     - is ice/mixed
! snowinc          real     output  3hr snow accum
!                                       array(2-d)
! snowtot          real     output  total snow accum
!                                      array(2-d)
! tsfc             real     input   surface temp array (2-d)
! tsfcf            real     local   surface temp in f
!
! updates
! =======
! 20 feb 98  initial version..................capt. john lewis/dnxt
! 10 nov 98  changed preciptype flag for snow
!            from 4 to 5; this corresponds
!            with changes made to the precip
!            type algorithm (wintprec.f)...capt david beberwyk/dnxt
!  5 jan 99  removed interpolations to dot grid as mmpost now
!            operates on the cross grid........................dnxm
!  4 jan 01  adapted by fsl for use with laps.. b. shaw, noaa/fsl
!    jan 12  use column max temperature and analysis library routine
!            to obtain snow_to_rain ratio (s. albers)

   use constants

   implicit none

   integer :: i, imax, j, jmax, ksig
   real :: fahren, tsfcf, tcolmax, snow_to_rain_ratio
   real, dimension(imax, jmax) :: prcpinc, preciptype, snowinc, snowtot, tsfc, tdsfc
   real, dimension(imax, jmax, ksig) :: tsig

   do j = 1, jmax
   do i = 1, imax

! check if precipitation type is snow.

      if (preciptype(i, j) == 2) then

         if (.true.) then ! liquid equivalent of snow depends on column max temp
            tcolmax = maxval(tsig(i, j, :))
            snowinc(i, j) = snow_to_rain_ratio(tcolmax)*prcpinc(i, j)

         else ! liquid equivalent of snow depends on surface temperature.
            tsfcf = fahren(tsfc(i, j) - t0)
            if (tsfcf >= 10.0) then
               snowinc(i, j) = 10.*prcpinc(i, j)
            else
               snowinc(i, j) = 15.*prcpinc(i, j)
            end if
         end if

! if precip type is not snow then snow accum is zero.

      else
         snowinc(i, j) = 0.
      end if

   end do
   end do

! update snow total.

   snowtot = snowtot + snowinc

   return
end

!===============================================================================

subroutine height_tw(pr, ht, tp, td, spr, sht, stp, std, smr, slp, threshold, ztw0, lx, ly, nz)

   implicit none

   integer :: lx, ly, nz, i, j, k, ku, kl
   real, parameter :: lapse = 0.0065
   real :: tw, twl, twu, rat, pr0, zbot, ztop, t0, rh0, td0, relhum, dewpt, threshold
   real, dimension(lx, ly, nz) :: pr, ht, tp, td
   real, dimension(lx, ly) :: spr, sht, stp, std, smr, slp, ztw0

   ku = -999
   kl = -999

   do j = 1, ly
   do i = 1, lx
      if (maxval(tp(i, j, :)) .lt. 273.15 + threshold) then ! subfreezing sigma column
         goto 2
      end if
      if (kl .ne. -999) then ! use the previously found bracketing levels for efficiency
         twu = tw(tp(i, j, ku), td(i, j, ku), pr(i, j, ku))
         do k = kl, kl
            twl = tw(tp(i, j, k), td(i, j, k), pr(i, j, k))
            if (twu <= 273.15 + threshold .and. twl > 273.15 + threshold) then
               rat = (twl - (273.15 + threshold))/(twl - twu)
               pr0 = pr(i, j, k) + rat*(pr(i, j, k + 1) - pr(i, j, k))
               if (pr0 .lt. 0.) then
                  write (6, *) ' error in height_tw: pr0 < 0. ', pr0
                  go to 1
               end if
               zbot = ht(i, j, k)
               ztop = ht(i, j, k + 1)
               rat = alog(pr0/pr(i, j, k))/alog(pr(i, j, k + 1)/pr(i, j, k))
               ztw0(i, j) = zbot + rat*(ztop - zbot)
               goto 1
            end if
         end do
      end if
      twu = tw(tp(i, j, nz), td(i, j, nz), pr(i, j, nz))
      do k = nz - 1, 1, -1
         twl = tw(tp(i, j, k), td(i, j, k), pr(i, j, k))
         if (twu <= 273.15 + threshold .and. twl > 273.15 + threshold) then
            rat = (twl - (273.15 + threshold))/(twl - twu)
            pr0 = pr(i, j, k) + rat*(pr(i, j, k + 1) - pr(i, j, k))
            if (pr0 .lt. 0.) then
               write (6, *) ' error in height_tw: pr0 < 0. ', pr0
               go to 1
            end if
            zbot = ht(i, j, k)
            ztop = ht(i, j, k + 1)
            rat = alog(pr0/pr(i, j, k))/alog(pr(i, j, k + 1)/pr(i, j, k))
            ztw0(i, j) = zbot + rat*(ztop - zbot)
            kl = k
            ku = k + 1
            goto 1
         end if
         twu = twl
      end do
2     continue ! check surface wet bulb
      twl = tw(stp(i, j), std(i, j), spr(i, j))
      if (twu <= 273.15 + threshold .and. twl > 273.15 + threshold) then
         rat = (twl - (273.15 + threshold))/(twl - twu)
         pr0 = spr(i, j) + rat*(pr(i, j, 1) - spr(i, j))
         zbot = sht(i, j)
         ztop = ht(i, j, 1)
         rat = alog(pr0/spr(i, j))/alog(pr(i, j, 1)/spr(i, j))
         ztw0(i, j) = zbot + rat*(ztop - zbot)
      else
! height of wet-bulb = 0c is below ground.
! calculate wet-bulb temp at height = 0 and interpolate.
         twu = twl
         t0 = stp(i, j) + lapse*sht(i, j)
         rh0 = min(relhum(t0, smr, slp(i, j)), 1.)
         td0 = dewpt(t0, rh0)
         twl = tw(t0, td0, slp(i, j))
         if (twu <= 273.15 + threshold .and. twl > 273.15 + threshold) then
            rat = (twl - (273.15 + threshold))/(twl - twu)
            pr0 = slp(i, j) + rat*(spr(i, j) - slp(i, j))
            ztop = sht(i, j)
            rat = alog(pr0/slp(i, j))/alog(spr(i, j)/slp(i, j))
            ztw0(i, j) = rat*ztop
         else
            ztw0(i, j) = 0.
         end if
      end if
1     continue
   end do
   end do

   return
end

subroutine get_ghi_ratio(i4time1, i4time2, lat, lon, ni, nj, ghi_ratio)

! ratio of ghi at time 2 to ghi at time 1

   include 'trigd.inc'

   real lat(ni, nj)
   real lon(ni, nj)

   real solalt1(ni, nj)
   real solalt2(ni, nj)

   real ghi_ratio(ni, nj)

   call get_solalt_2d(lat, lon, i4time1, ni, nj, solalt1)
   call get_solalt_2d(lat, lon, i4time2, ni, nj, solalt2)

   write (6, *) ' range of solalt1 is ', minval(solalt1), maxval(solalt1)
   write (6, *) ' range of solalt2 is ', minval(solalt2), maxval(solalt2)

   do i = 1, ni
   do j = 1, nj

      floor = max(.01*(1.+solalt1(i, j)/3.), 0.)
      zenfrac1 = max(sind(solalt1(i, j)), floor)

      floor = max(.01*(1.+solalt2(i, j)/3.), 0.)
      zenfrac2 = max(sind(solalt2(i, j)), floor)

      if (zenfrac1 .gt. 0) then
         ratio = zenfrac2/zenfrac1
      else
         ratio = 1.
      end if

!   if(i .eq. ni/2)then
!       write(6,*)'solalt1/2,zenfrac1,zenfrac2,ratio',i,j,solalt1(i,j),solalt2(i,j),zenfrac1,zenfrac2,ratio
!   endif

      ghi_ratio(i, j) = min(max(ratio, .1), 10.)

   end do ! j
   end do ! i

   return
end

