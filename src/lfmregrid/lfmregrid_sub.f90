
!gfortran subroutine lfmregrid_sub(nx_bg,ny_bg,nzbg,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,fname_in,nx_l,ny_l,nz_l,gproj,bgmodel,cmodel,laps_data_root,mtype,laps_reftime,laps_valtime,l_process_grib,l_process_cdf,l_grib_fua,l_grib_fsf,l_cdf_fua,l_cdf_fsf)
!gfortran modifications begin
subroutine lfmregrid_sub(nx_bg, ny_bg, nzbg, nzbg_tp, nzbg_sh, nzbg_uv, nzbg_ww, &
                         fname_in, nx_l, ny_l, nz_l, gproj, bgmodel, cmodel, &
                         laps_data_root, mtype, laps_reftime, laps_valtime, &
                         l_process_grib, l_process_cdf, l_grib_fua, l_grib_fsf, &
                         l_cdf_fua, l_cdf_fsf)
!gfortran modifications end

!use mem_namelist
   use storage_module, only: get_plvls
   integer, parameter :: maxbglvl = 52 ! dimension is maxbglvl in 'degrib_nav' routine

   character*256 fname_in, fullname_in, laps_data_root, fname_bg, bgpath
   character*132 cmodel
   character*4 af_bg
   character*2 gproj
   character*256 syscmd
   character*30 mtype

   integer bgmodel

   integer ct, n2df, n3df
   parameter(n2df=11)
   parameter(n3df=4)

   real sgrid(nx_l, ny_l, n2df)
   real pgrid(nx_l, ny_l, n3df*nz_l)
   character*132 com2d(n2df)
   character*132 com3d(n3df*nz_l)
   character*3 name2d(n2df)
   character*3 name3d(n3df*nz_l)
   integer lvls3d(n3df*nz_l)  ! mb for pressure grid

   real pres_1d(nz_l) ! pa
   real pr1d_pa(nz_l) ! pa
   real pr1d_mb(nz_l) ! mb

! *** background model grid data.
!
   real prbght(nx_bg, ny_bg, nzbg) !pressure (mb) ht and temp
   real prbgsh(nx_bg, ny_bg, nzbg) !pressure (mb) q
   real prbguv(nx_bg, ny_bg, nzbg) !pressure (mb) u- v-components
   real prbgww(nx_bg, ny_bg, nzbg) !pressure (mb) omega
   real pr(nzbg)

   real htbg(nx_bg, ny_bg, nzbg)   !height (m)
   real tpbg(nx_bg, ny_bg, nzbg)   !temperature (k)
   real shbg(nx_bg, ny_bg, nzbg)   !specific humidity (kg/kg)
   real tdbg(nx_bg, ny_bg, nzbg)   !dewpoint (k)
   real uwbg(nx_bg, ny_bg, nzbg)   !u-wind (m/s)
   real vwbg(nx_bg, ny_bg, nzbg)   !v-wind (m/s)
   real wwbg(nx_bg, ny_bg, nzbg)   !w-wind (pa/s)
   real lwbg(nx_bg, ny_bg, nzbg)   !cloud liquid
   real icbg(nx_bg, ny_bg, nzbg)   !cloud ice
   real plvl_grib(maxbglvl) ! dimension is maxbglvl in 'degrib_nav' routine

   real mslpbg(nx_bg, ny_bg)         !mslp  (mb)
   real htbg_sfc(nx_bg, ny_bg)
   real prbg_sfc(nx_bg, ny_bg)
   real shbg_sfc(nx_bg, ny_bg)       !specific humidity (kg/kg)
   real uwbg_sfc(nx_bg, ny_bg)
   real vwbg_sfc(nx_bg, ny_bg)
   real tdbg_sfc(nx_bg, ny_bg)
   real tpbg_sfc(nx_bg, ny_bg)
   real t_at_sfc(nx_bg, ny_bg)
   real pcpbg(nx_bg, ny_bg)          !precip at surface, acpc (k/m^2)

!     local input model variables for the time being
   real r01(nx_bg, ny_bg)
   real lmr(nx_bg, ny_bg)
   real llr(nx_bg, ny_bg)
   real swi(nx_bg, ny_bg)
   real s8a(nx_bg, ny_bg)
   real tpw(nx_bg, ny_bg)

   real htvi(nx_bg, ny_bg, nz_l)
   real tpvi(nx_bg, ny_bg, nz_l)
   real shvi(nx_bg, ny_bg, nz_l)
   real uwvi(nx_bg, ny_bg, nz_l)
   real vwvi(nx_bg, ny_bg, nz_l)
   real wwvi(nx_bg, ny_bg, nz_l)

   real lat(nx_l, ny_l)         ! laps lat
   real lon(nx_l, ny_l)         ! laps lon
   real topo(nx_l, ny_l)        ! laps lon
   real grx(nx_l, ny_l)         ! hinterp factor (background grid x index)
   real gry(nx_l, ny_l)         ! hinterp factor (background grid y index)
   real lmr_laps(nx_l, ny_l)
   real tsf_laps(nx_l, ny_l)
   real dsf_laps(nx_l, ny_l)
   real psf_laps(nx_l, ny_l)
   real usf_laps(nx_l, ny_l)
   real vsf_laps(nx_l, ny_l)
   real swi_laps(nx_l, ny_l)
   real s8a_laps(nx_l, ny_l)
   real tpw_laps(nx_l, ny_l)
   real rto_laps(nx_l, ny_l)
   real r01_laps(nx_l, ny_l)

   real ssh, k_to_c

   logical wrapped, l_process_grib, l_process_cdf
   logical l_grib_fua, l_grib_fsf, l_cdf_fua, l_cdf_fsf

   write (6, *)
   write (6, *) ' subroutine lfmregrid_sub...'
   write (6, *) ' nx_bg/ny_bg = ', nx_bg, ny_bg
   write (6, *) ' cmodel = ', cmodel

   write (6, *) ' l_grib_fua,l_grib_fsf,l_cdf_fua,l_cdf_fsf:', l_grib_fua, l_grib_fsf, l_cdf_fua, l_cdf_fsf

   call get_laps_domain(nx_l, ny_l, 'nest7grid' &
                        , lat, lon, topo, istatus)
   if (istatus .lt. 1) then
      print *, 'error reading lat, lon, topo data from get_laps_domain'
      stop
   end if

   call get_pres_1d(i4_valtime, nz_l, pr1d_pa, istatus)
   if (istatus .ne. 1) then
      print *, 'error returned from get_pres_1d'
      print *, 'check pressures.nl or nk_laps in nest7grid.parms'
      stop
   end if
   pr1d_mb(:) = pr1d_pa(:)/100.  ! pa to mb

   if (l_process_grib .eqv. .true.) then
!   bgmodel = 13
!   if(mtype .eq. 'nam')then
!       cmodel = 'nam'
!   else
!       cmodel = 'hrrr'
!   endif
!   write(6,*)' calling get_bkgd_mdl_info for cmodel: ',trim(cmodel)
!   call get_bkgd_mdl_info(bgmodel,cmodel,fname_in  &
!     ,nx_bg,ny_bg,nzbg,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww &
!     ,gproj,dlat,dlon,centrallat,centrallon,dxbg,dybg &
!     ,lat0,lat1,lon0,sw,ne,cgrddef,istatus)

      if (.true.) then
         call get_basename(fname_in, fname_bg)
         call get_directory_length(fname_in, lend)
         bgpath = fname_in(1:lend)
         af_bg = fname_bg(10:13)
         write (6, *) ' calling read_bgdata: dims are ', nx_bg, ny_bg, nzbg, nzbg_tp, nzbg_sh, nzbg_uv, nzbg_ww
         call read_bgdata(nx_bg, ny_bg, &
                          nzbg, nzbg_tp, nzbg_sh, nzbg_uv, nzbg_ww, 'lapsb' &
                          , bgpath, fname_bg, af_bg, fname_in, cmodel, bgmodel &
                          , prbght, prbgsh, prbguv, prbgww &
                          , htbg, tpbg, uwbg, vwbg, shbg, wwbg &
                          , htbg_sfc, prbg_sfc, shbg_sfc, tdbg_sfc, tpbg_sfc &
                          , t_at_sfc, uwbg_sfc, vwbg_sfc, mslpbg, pcpbg, lmr, tpw, swi, istatus)

         write (6, *) ' returned from read_bgdata'

         istatus = ishow_timer()

         nx_pr = 1
         ny_pr = 1

         call get_plvls(plvl_grib, 100, nlvl_grib)
         write (6, *) ' grib plvls info: ', nlvl_grib, plvl_grib(1:nlvl_grib)
         nzbg_ht = 0
         do k = 1, nlvl_grib
            if (plvl_grib(k) .lt. 150000) then
               nzbg_ht = nzbg_ht + 1
            end if
         end do ! k

         write (6, *) ' reset nzbg_ht, etc. to ', nzbg_ht

         nzbg_tp = nzbg_ht
         nzbg_sh = nzbg_ht
         nzbg_uv = nzbg_ht
         nzbg_ww = nzbg_ht

         ixmin = 1
         ixmax = nx_bg
         iymin = 1
         iymax = ny_bg

         nz_laps = nz_l
         write (6, *) ' calling vinterp: dims are ', nx_bg, ny_bg, nzbg_ht, nz_laps
         call vinterp(nz_laps, nx_bg, ny_bg, nx_pr, ny_pr &
                      , ixmin, ixmax, iymin, iymax &
                      , nzbg_ht, nzbg_tp, nzbg_sh, nzbg_uv, nzbg_ww &
                      , pr1d_mb, prbght, prbgsh, prbguv, prbgww &
                      , htbg, tpbg, shbg, uwbg, vwbg, wwbg &
                      , htvi, tpvi, shvi, uwvi, vwvi, wwvi)

         nx_laps = nx_l
         ny_laps = ny_l

!       if((l_cdf_fua .eqv. .true.) .or. (l_grib_fua .eqv. .true.))then
         if (.false.) then

            write (6, *) ' calling init_hinterp'
            wrapped = .false.
            bgmodel = 0
            call init_hinterp(nx_bg, ny_bg, nx_l, ny_l, gproj, &
                              lat, lon, grx, gry, bgmodel, cmodel, wrapped)

            print *, 'use bilinear_laps_3d for t ', trim(cmodel)
            call bilinear_laps_3df(grx, gry, nx_bg, ny_bg &
                                   , nx_laps, ny_laps, nz_laps, tpvi, tp)

            istatus = ishow_timer()
            print *, 'use bilinear_laps_3d for q ', trim(cmodel)
            call bilinear_laps_3df(grx, gry, nx_bg, ny_bg &
                                   , nx_laps, ny_laps, nz_laps, shvi, sh)

         end if

         if (istatus .eq. 1) then
            write (6, *) ' calling write_lga'
            call write_lga(nx_laps, ny_laps, nz_laps, laps_reftime, &
                           bgvalid, cmodel, missingflag, pr1d_mb, ht, tp, sh, uw, vw, ww, istatus)
            if (istatus .ne. 1) then
               print *, 'error writing lga - returning to main'
            end if

         end if
      end if
   end if

   if (l_process_cdf .eqv. .true.) then
!   fullname_in = trim(fname_in)//'.fsf'

      write (6, *) ' calling read_fuafsf_cdf'

!   initialize variables
      tpw = r_missing_data

      call read_fuafsf_cdf(fname_in &
                           , nx_bg, ny_bg, nzbg &
                           , htbg, pr, wwbg, shbg, tpbg, uwbg, vwbg &
                           , uwbg_sfc, vwbg_sfc, tpbg_sfc, tdbg_sfc &
                           , prbg_sfc, mslpbg, htbg_sfc &
                           , r01 &
                           , pcpbg &
                           , lmr, llr, s8a, swi, tpw &
                           , istatus)
      if (istatus .ne. 1 .and. istatus .ne. -1) then
         print *, 'error returned: read_fuafsf_cdf'
         return
      end if

   end if

! if(l_process_cdf .eqv. .true.)then
   if (.true.) then
      wrapped = .false.
      bgmodel = 0

      write (6, *) ' calling init_hinterp'

      call init_hinterp(nx_bg, ny_bg, nx_l, ny_l, gproj, &
                        lat, lon, grx, gry, bgmodel, cmodel, wrapped)

      if ((l_grib_fua .eqv. .true.) .or. (l_cdf_fua .eqv. .true.)) then
         call get_pres_1d(i4_valtime, nz_l, pres_1d, istatus)
         lz = nz_l
         ct = 1

         print *, 'use bilinear_laps_3df for t3 starting at pgrid level', ct
         call bilinear_laps_3df(grx, gry, nx_bg, ny_bg &
                                , nx_l, ny_l, nz_l, tpbg, pgrid(1, 1, ct))
         name3d(ct:ct + lz - 1) = 't3 '; com3d(ct:ct + lz - 1) = 'temperature'; lvls3d(ct:ct + lz - 1) = nint(pres_1d(lz:1:-1)/100.); ct = ct + lz

         if (trim(cmodel) .eq. 'hrrr') then
!         convert sh from dpt into actual sh
            itesth = 1; jtesth = 1
            itestl = 506; jtestl = 128
            if (nx_l .gt. itestl .and. ny_l .gt. jtestl) then ! dfw laps grid location (hwt domain)
               itesth = nint(grx(506, 128))
               jtesth = nint(gry(506, 128))
               write (6, *) ' convert sh from dpt to actual sh at lat/lon ', lat(itestl, jtestl), lon(itestl, jtestl)
               write (6, *) ' laps gridpt', itestl, jtestl, ' hrrr gridpt', itesth, jtesth
            end if

            do k = 1, nz_l
               kflip = (nz_l + 1) - k
               p_mb = pres_1d(kflip)/100.
               do i = 1, nx_bg
               do j = 1, ny_bg
                  sh_orig = shbg(i, j, k)
                  td_c = k_to_c(shbg(i, j, k))
                  shbg(i, j, k) = ssh(p_mb, td_c)/1000.
                  if (i .eq. itesth .and. j .eq. jtesth) then
                     write (6, 11) p_mb, tpbg(i, j, k), sh_orig, td_c, shbg(i, j, k)
                  end if
11                format(' p,t,sh_orig,td_c,sh: ', 4f10.3, f10.4)
               end do ! j
               end do ! i
            end do ! k
         end if

         istatus = ishow_timer()

         print *, 'use bilinear_laps_3df for sh starting at pgrid level', ct
         write (6, *) shbg(1, 1, :)
         call bilinear_laps_3df(grx, gry, nx_bg, ny_bg &
                                , nx_l, ny_l, nz_l, shbg, pgrid(1, 1, ct))
         write (6, *) pgrid(1, 1, ct:ct + lz - 1)
         name3d(ct:ct + lz - 1) = 'sh '; com3d(ct:ct + lz - 1) = 'specific humidity'; lvls3d(ct:ct + lz - 1) = nint(pres_1d(lz:1:-1)/100.); ct = ct + lz

         istatus = ishow_timer()

         print *, 'use bilinear_laps_3df for u3 starting at pgrid level', ct
         call bilinear_laps_3df(grx, gry, nx_bg, ny_bg &
                                , nx_l, ny_l, nz_l, uwbg, pgrid(1, 1, ct))
         name3d(ct:ct + lz - 1) = 'u3 '; com3d(ct:ct + lz - 1) = 'u-component wind'; lvls3d(ct:ct + lz - 1) = nint(pres_1d(lz:1:-1)/100.); ct = ct + lz

         istatus = ishow_timer()

         print *, 'use bilinear_laps_3df for v3 starting at pgrid level', ct
         call bilinear_laps_3df(grx, gry, nx_bg, ny_bg &
                                , nx_l, ny_l, nz_l, vwbg, pgrid(1, 1, ct))
         name3d(ct:ct + lz - 1) = 'v3 '; com3d(ct:ct + lz - 1) = 'v-component wind'; lvls3d(ct:ct + lz - 1) = nint(pres_1d(lz:1:-1)/100.); ct = ct + lz

         istatus = ishow_timer()

      else
         write (6, *) ' skip 3d interpolation - fua processing set to false'

      end if

      if (.false.) then

         write (6, *) ' calling hinterp_field_2d for lmr'
         call hinterp_field_2d(nx_bg, ny_bg, nx_l, ny_l, 1, &
                               grx, gry, lmr, lmr_laps, wrapped)

         write (6, *) ' calling hinterp_field_2d for tsf'
         call hinterp_field_2d(nx_bg, ny_bg, nx_l, ny_l, 1, &
                               grx, gry, tpbg_sfc, tsf_laps, wrapped)

         write (6, *) ' calling hinterp_field_2d for dsf'
         call hinterp_field_2d(nx_bg, ny_bg, nx_l, ny_l, 1, &
                               grx, gry, tdbg_sfc, dsf_laps, wrapped)

         write (6, *) ' calling hinterp_field_2d for psf'
         call hinterp_field_2d(nx_bg, ny_bg, nx_l, ny_l, 1, &
                               grx, gry, prbg_sfc, psf_laps, wrapped)

         write (6, *) ' calling hinterp_field_2d for usf'
         call hinterp_field_2d(nx_bg, ny_bg, nx_l, ny_l, 1, &
                               grx, gry, uwbg_sfc, usf_laps, wrapped)

         write (6, *) ' calling hinterp_field_2d for vsf'
         call hinterp_field_2d(nx_bg, ny_bg, nx_l, ny_l, 1, &
                               grx, gry, vwbg_sfc, vsf_laps, wrapped)

         write (6, *) ' calling hinterp_field_2d for swi'
         call hinterp_field_2d(nx_bg, ny_bg, nx_l, ny_l, 1, &
                               grx, gry, swi, swi_laps, wrapped)

         write (6, *) ' calling hinterp_field_2d for rto'
         call hinterp_field_2d(nx_bg, ny_bg, nx_l, ny_l, 1, &
                               grx, gry, pcpbg, rto_laps, wrapped)

         write (6, *) ' calling hinterp_field_2d for r01'
         call hinterp_field_2d(nx_bg, ny_bg, nx_l, ny_l, 1, &
                               grx, gry, pcpbg, r01_laps, wrapped)

         write (6, *) ' calling hinterp_field_2d for tpw'
         call hinterp_field_2d(nx_bg, ny_bg, nx_l, ny_l, 1, &
                               grx, gry, tpw, tpw_laps, wrapped)

      else

         write (6, *) ' calling bilinear_laps_2d for lmr'
         call bilinear_laps_2d(grx, gry, nx_bg, ny_bg, nx_l, ny_l, lmr, lmr_laps)

         write (6, *) ' calling bilinear_laps_2d for tsf'
         call bilinear_laps_2d(grx, gry, nx_bg, ny_bg, nx_l, ny_l, tpbg_sfc, tsf_laps)

         write (6, *) ' calling bilinear_laps_2d for dsf'
         call bilinear_laps_2d(grx, gry, nx_bg, ny_bg, nx_l, ny_l, tdbg_sfc, dsf_laps)

         write (6, *) ' calling bilinear_laps_2d for psf'
         call bilinear_laps_2d(grx, gry, nx_bg, ny_bg, nx_l, ny_l, prbg_sfc, psf_laps)

         write (6, *) ' calling bilinear_laps_2d for usf'
         call bilinear_laps_2d(grx, gry, nx_bg, ny_bg, nx_l, ny_l, uwbg_sfc, usf_laps)

         write (6, *) ' calling bilinear_laps_2d for vsf'
         call bilinear_laps_2d(grx, gry, nx_bg, ny_bg, nx_l, ny_l, vwbg_sfc, vsf_laps)

         write (6, *) ' calling bilinear_laps_2d for swi'
         call bilinear_laps_2d(grx, gry, nx_bg, ny_bg, nx_l, ny_l, swi, swi_laps)

         write (6, *) ' calling bilinear_laps_2d for rto'
         call bilinear_laps_2d(grx, gry, nx_bg, ny_bg, nx_l, ny_l, pcpbg, rto_laps)

         write (6, *) ' calling bilinear_laps_2d for r01'
         call bilinear_laps_2d(grx, gry, nx_bg, ny_bg, nx_l, ny_l, pcpbg, r01_laps)

         write (6, *) ' calling bilinear_laps_2d for tpw'
         call bilinear_laps_2d(grx, gry, nx_bg, ny_bg, nx_l, ny_l, tpw, tpw_laps)

      end if

      s8a_laps = r_missing_data

      where (rto_laps(:, :) .ne. r_missing_data)
         rto_laps(:, :) = rto_laps(:, :)/1000.
      end where
      where (r01_laps(:, :) .ne. r_missing_data)
         r01_laps(:, :) = r01_laps(:, :)/1000.
      end where
      where (tpw_laps(:, :) .ne. r_missing_data)
         tpw_laps(:, :) = tpw_laps(:, :)/1000.
      end where

      write (6, *) ' lmr ranges: ', minval(lmr), maxval(lmr), minval(lmr_laps), maxval(lmr_laps)
      write (6, *) ' tsf ranges: ', minval(tpbg_sfc), maxval(tpbg_sfc), minval(tsf_laps), maxval(tsf_laps)
      write (6, *) ' dsf ranges: ', minval(tdbg_sfc), maxval(tdbg_sfc), minval(dsf_laps), maxval(dsf_laps)
      write (6, *) ' psf ranges: ', minval(prbg_sfc), maxval(prbg_sfc), minval(psf_laps), maxval(psf_laps)
      write (6, *) ' usf ranges: ', minval(uwbg_sfc), maxval(uwbg_sfc), minval(usf_laps), maxval(usf_laps)
      write (6, *) ' vsf ranges: ', minval(vwbg_sfc), maxval(vwbg_sfc), minval(vsf_laps), maxval(vsf_laps)
      write (6, *) ' swi ranges: ', minval(swi), maxval(swi), minval(swi_laps), maxval(swi_laps)
      write (6, *) ' s8a ranges: ', minval(s8a), maxval(s8a), minval(s8a_laps), maxval(s8a_laps)
      write (6, *) ' rto ranges: ', minval(pcpbg), maxval(pcpbg), minval(rto_laps), maxval(rto_laps)
      write (6, *) ' r01 ranges: ', minval(r01), maxval(r01), minval(r01_laps), maxval(r01_laps)
      write (6, *) ' tpw ranges: ', minval(tpw), maxval(tpw), minval(tpw_laps), maxval(tpw_laps)

      write (6, *) ' setup output arrays'
      ct = 1
      sgrid(1:nx_l, 1:ny_l, ct) = lmr_laps; name2d(ct) = 'lmr'; com2d(ct) = 'composite reflectivity'; ct = ct + 1
      sgrid(1:nx_l, 1:ny_l, ct) = tsf_laps; name2d(ct) = 'tsf'; com2d(ct) = 'sfc temperature'; ct = ct + 1
      sgrid(1:nx_l, 1:ny_l, ct) = dsf_laps; name2d(ct) = 'dsf'; com2d(ct) = 'sfc dewpoint temperature'; ct = ct + 1
      sgrid(1:nx_l, 1:ny_l, ct) = psf_laps; name2d(ct) = 'psf'; com2d(ct) = 'surface pressure'; ct = ct + 1
      sgrid(1:nx_l, 1:ny_l, ct) = usf_laps; name2d(ct) = 'usf'; com2d(ct) = 'sfc u-component wind'; ct = ct + 1
      sgrid(1:nx_l, 1:ny_l, ct) = vsf_laps; name2d(ct) = 'vsf'; com2d(ct) = 'sfc v-component wind'; ct = ct + 1
      sgrid(1:nx_l, 1:ny_l, ct) = swi_laps; name2d(ct) = 'swi'; com2d(ct) = 'incoming sw radiation'; ct = ct + 1
      sgrid(1:nx_l, 1:ny_l, ct) = s8a_laps; name2d(ct) = 's8a'; com2d(ct) = '11u brightness temperature'; ct = ct + 1
      sgrid(1:nx_l, 1:ny_l, ct) = rto_laps; name2d(ct) = 'rto'; com2d(ct) = 'run total precip accum'; ct = ct + 1
      sgrid(1:nx_l, 1:ny_l, ct) = r01_laps; name2d(ct) = 'r01'; com2d(ct) = 'incremental precip accum'; ct = ct + 1
      sgrid(1:nx_l, 1:ny_l, ct) = tpw_laps; name2d(ct) = 'tpw'; com2d(ct) = 'precipitable water depth'; ct = ct + 1

      if (ct - 1 .ne. n2df) then
         write (6, *) ' error in lfmregrid_sub: value of ct is inconsistent with n2df ', ct, n2df
      end if

!   if(n3df .eq. 0)nz_l = 1
!   syscmd = 'rm -f '//trim(fullname_in)
!   call system (syscmd)

      write (6, *)
      write (6, *) ' calling output_laps_rg: nx_l,ny_l,nz_l is ', nx_l, ny_l, nz_l
      write (6, *) ' laps_data_root is: ', laps_data_root
      write (6, *) ' pgrid is ', pgrid(1, 1, :)
!gfortran    call output_laps_rg(laps_data_root,mtype,domnum_in,laps_reftime,laps_valtime,pgrid,sgrid,name2d,name3d,com2d,com3d,lvls3d,n2df,n3df,pres_1d,nx_l,ny_l,nz_l)
!gfortran modifications begin
      call output_laps_rg(laps_data_root, mtype, domnum_in, laps_reftime, &
                          laps_valtime, pgrid, sgrid, name2d, name3d, com2d, com3d, &
                          lvls3d, n2df, n3df, pres_1d, nx_l, ny_l, nz_l)
!gfortran modifications end

   end if ! .true. (formerly l_process_cdf)

   write (6, *) ' end of lfmregrid_sub'

! stop ! prevent seg fault upon return

   return

end

subroutine get_basename(fullname, basename)

   character*(*) fullname, basename

   call get_directory_length(fullname, lend)
   call s_len(fullname, lenfn)

   basename = fullname(lend + 1:lenfn)

   write (6, *) ' basename = ', trim(basename)

   return
end
