program lfmregrid

   use mem_namelist

   character*256 fname_in, fullname_in, laps_data_root
   character*132 cmodel
   character*150 static_dir, filename

   character*16 atime, afcst
   character*9 a9time
   character*30 projname
   character*2 gproj
   character*1 cgrddef
   character*30 mtype

   integer fcsttime, fcsthr, fcstmn, bgmodel

   real la1, lo1, la1in, la2in, lov, la2, lo2
   real lat0, lat1, lon0
   real sw(2), ne(2)

   logical l_process_grib, l_process_cdf, l_parse
   character*1 a_process_grib, a_process_cdf

   logical l_grib_fua, l_grib_fsf, l_cdf_fua, l_cdf_fsf
   character*1 a_grib_fua, a_grib_fsf, a_cdf_fua, a_cdf_fsf

!
! *** common block variables for polar stereographic grid.
!
   integer nx_ps, ny_ps, nz_ps    !no. of ps domain grid points
   real lat0_ps, lon0_ps, rota &  !pol ste. std lat, lon and rotation
      , sw_ps(2), ne_ps(2)     !sw lat, lon, ne lat, lon
   common/psgrid/nx_ps, ny_ps, nz_ps, lat0_ps, lon0_ps &
      , rota, sw_ps, ne_ps

!
! *** common block variables for lambert-conformal grid.
!
   integer nx_lc, ny_lc, nz_lc
   real lat1_lc, lat2_lc, lon0_lc, sw_lc(2), ne_lc(2)
   common/lcgrid/nx_lc, ny_lc, nz_lc, lat1_lc, lat2_lc &
      , lon0_lc, sw_lc, ne_lc

   istatus = init_timer()

   call getarg(1, fname_in)       ! input file name (without the .f?? extension)
   call getarg(2, a9time)         ! ascii 9 character time for initialization
   call getarg(3, afcst)          ! hhmm of the forecast
   call getarg(4, a_grib_fua)
   call getarg(5, a_grib_fsf)
   call getarg(6, a_cdf_fua)
   call getarg(7, a_cdf_fsf)
   call getarg(8, mtype)
   call getarg(9, laps_data_root)

   call cv_asc_i4time(a9time, laps_reftime)
   read (afcst, '(i2)', err=900) fcsthr
   read (afcst, '(2x,i2)', err=900) fcstmn
   fcsttime = fcsthr*3600 + fcstmn*60
   laps_valtime = laps_reftime + fcsttime
   read (a_grib_fua, *) l_grib_fua
   read (a_grib_fsf, *) l_grib_fsf
   read (a_cdf_fua, *) l_cdf_fua
   read (a_cdf_fsf, *) l_cdf_fsf

   l_process_grib = (l_grib_fua .or. l_grib_fsf)
   l_process_cdf = (l_cdf_fua .or. l_cdf_fsf)

   write (6, *) ' a9time/laps_reftime ', a9time, ' ', laps_reftime
   write (6, *) ' l_process_grib, l_process_cdf:', l_process_grib, l_process_cdf
   write (6, *) ' l_grib_fua,l_grib_fsf,l_cdf_fua,l_cdf_fsf:', l_grib_fua, l_grib_fsf, l_cdf_fua, l_cdf_fsf
   write (6, *) ' fname_in: ', trim(fname_in)
   write (6, *) ' mtype: ', trim(mtype)
   write (6, *)

   if (l_process_cdf) then

      if (l_cdf_fua .eqv. .true.) then
         fullname_in = trim(fname_in)//'.fua'
      else
         fullname_in = trim(fname_in)//'.fsf'
      end if

!   obtain input file horizontal dimensions, note that z will be missing for fsf case
      write (6, *) ' get dimensions from input file: ', trim(fullname_in)
      call getdims_lapsprd(fullname_in, nxbg, nybg, nzbg, istatus)
      if (istatus .ne. 1) then
         goto 900
      end if

      write (6, *) ' input nxbg,nybg,nzbg = ', nxbg, nybg, nzbg

      if (.false.) then
         write (6, *) ' get nav info from input file... '
         call read_lapsprd_attr(fullname_in, &
                                dxbg, dybg, la1, lo1, la1in, la2in, lov, &
                                projname, la2, lo2, istatus)
         if (istatus .ne. 1) then
            print *, 'error returned: read_lapsprd_attr'
            goto 900
         end if
      else
         write (6, *) ' assume hard wired nav info from input file'
         dxbg = 3000.
         dybg = 3000.
         la1 = 21.13812
         lo1 = -122.7195
         la1in = 38.500
         la2in = 38.500
         lov = 262.500
         projname = 'tangential lambert conformal'
         la2 = 47.84364
         lo2 = -60.90137
      end if

      if (lo1 .gt. 180) lo1 = lo1 - 360
      if (lo2 .gt. 180) lo2 = lo2 - 360
      if (lov .gt. 180) lov = lov - 360
      nzbg_ht = nzbg
      nzbg_tp = nzbg
      nzbg_sh = nzbg
      nzbg_uv = nzbg
      nzbg_ww = nzbg
      sw(1) = la1
      sw(2) = lo1
      ne(1) = la2
      ne(2) = lo2
      lon0 = lov
      lat0 = la1in
      cenlat = la1in
      cenlon = lov
      lat1 = la2in

!     specify whether the standard lat is southern or northern boundary
      cgrddef = 'n'

      call s_len2(projname, l)

      if (projname(1:l) .eq. 'polar' .or. &
          projname(1:l) .eq. 'polar stereographic') then
         gproj = 'ps'
      elseif (projname(1:l) .eq. 'mercator') then
         gproj = 'mc'
      elseif (projname(1:l) .eq. 'lambert') then
         gproj = 'lc'
      elseif (projname(1:l) .eq. 'secant lambert conformal') then
         gproj = 'lc'
      elseif (projname(1:l) .eq. 'tangential lambert conformal') then
         gproj = 'lc'
      else
         print *, 'error: unable to determine gproj setting ', projname
      end if

      write (6, *) ' lat0/lat1/lon0 = ', lat0, lat1, lon0
      write (6, *) ' sw = ', sw
      write (6, *) ' ne = ', ne

      call init_gridconv_cmn(gproj, nxbg, nybg, nzbg_ht &
                             , dlat, dlon, cenlat, cenlon, lat0, lat1, lon0 &
                             , sw(1), sw(2), ne(1), ne(2), cgrddef, istatus)

      write (6, *)
      if (gproj .eq. 'ps') then
         write (6, *) ' psgrid common block variables'
         write (6, *) ' nx_ps,ny_ps,nz_ps ', nx_ps, ny_ps, nz_ps
         write (6, *) ' lat0_ps,lon0_ps,rota,sw_ps,ne_ps ', lat0_ps, lon0_ps, rota, sw_ps, ne_ps
      elseif (gproj .eq. 'lc') then
         write (6, *) ' lcgrid common block variables'
         write (6, *) ' nx_lc,ny_lc,nz_lc ', nx_lc, ny_lc, nz_lc
         write (6, *) ' lat1_lc,lat2_lc,lon0_lc,sw_lc,ne_lc ', lat1_lc, lat2_lc, lon0_lc, sw_lc, ne_lc
      end if

   end if ! l_process_cdf

! read global laps parameters into module memory structure
   call get_directory('static', static_dir, len_dir)
   filename = static_dir(1:len_dir)//'/nest7grid.parms'
   call read_namelist_laps('lapsparms', filename)

   if ((l_grib_fua .eqv. .true.) .or. (l_cdf_fua .eqv. .true.)) then
      nz_l = nk_laps
   else
      nz_l = 1
   end if

   write (6, *) ' grid dims nx_l,ny_l,nz_l = ', nx_l, ny_l, nz_l

   if (mtype .eq. 'wrf-hrrr') then
      cmodel = 'hrrr'
   else
      call upcase(mtype, cmodel)
   end if

   if (l_process_grib .eqv. .true.) then
      bgmodel = 13
      write (6, *) ' calling get_bkgd_mdl_info for cmodel: ', trim(cmodel)
      call get_bkgd_mdl_info(bgmodel, cmodel, fname_in &
                             , nxbg, nybg, nzbg, nzbg_tp, nzbg_sh, nzbg_uv, nzbg_ww &
                             , gproj, dlat, dlon, centrallat, centrallon, dxbg, dybg &
                             , lat0, lat1, lon0, sw, ne, cgrddef, istatus)
      if (istatus .ne. 1) then
         write (6, *) ' error: bad status return from get_bkgd_mdl_info'
         goto 900
      end if

      call init_gridconv_cmn(gproj, nxbg, nybg, nzbg &
                             , dlat, dlon, centrallat, centrallon, lat0, lat1, lon0 &
                             , sw(1), sw(2), ne(1), ne(2), cgrddef, istatus)
      if (istatus .ne. 1) then
         write (6, *) ' error: bad status return from init_gridconv_cmn'
         goto 900
      end if
   end if

   call lfmregrid_sub(nxbg, nybg, nzbg, nzbg_tp, nzbg_sh, nzbg_uv, nzbg_ww, fname_in, nx_l, ny_l, nz_l, gproj &
                      , bgmodel, cmodel, laps_data_root, mtype, laps_reftime, laps_valtime &
                      , l_process_grib, l_process_cdf, l_grib_fua, l_grib_fsf, l_cdf_fua, l_cdf_fsf)

   write (6, *) ' returned from lfmregrid_sub - program end'

900 continue

end
