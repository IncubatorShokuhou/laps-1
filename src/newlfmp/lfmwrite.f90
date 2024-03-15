subroutine output_laps(lfmprd_dir, laps_data_root, domnum_in, laps_reftime, laps_valtime)

! creates the laps *.fua and *.fsf files for model output.

! use of this routine requires a valid mm5_data_root that contains
! a subdirectory in mm5_data_root/mm5prd/dxx/fua (and fsf), where the
! xx is the two digit domain number.  the fsf.cdl and fua.cdl files
! in laps_data_root/cdl should also have dimensions equal to lx/ly.

! history
! =======
! initial version:  brent shaw, noaa/fsl 5 jan 01
! modified to output files to domain specific directories in
! mm5_data_root instead of laps_data_root.  2 nov 01

! note that the only variable that is supported in the netcdf files but
! not produced by this routine is the fire index.

   use lfmgrid, only: domnum, verbose, lx, ly, lz, lprs, mtype, write_to_lapsdir, &
                      nvar3dout, name3d, lvltype3d, units3d, com3d, pgrid, realtime, make_donefile, &
                      lvls2d, lvls3d, nvar2dout, name2d, lvltype2d, &
                      units2d, com2d, sgrid

   implicit none

   integer :: domnum_in, laps_reftime, laps_valtime, istatus, fnlen, extended
   real, allocatable, dimension(:) :: cdl_levels
   character(len=*) :: lfmprd_dir, laps_data_root
   character(len=256) :: output_file, donefile
   character(len=255) :: output_dir, cdl_dir
   character(len=9) :: gtime
   character(len=5) :: fcst_hhmm
   character(len=2) :: domnum_str
   logical :: outdir_defined

!beka
   integer istat, i4_elapsed, ishow_timer, init_timer

! lw assign domnum_in, from calling args, to domnum, declared in lfmgrid module
   domnum = domnum_in

   if (verbose) then
      print *, ' '
      print *, 'outputting laps format (fua/fsf) files...'
   end if

   write (domnum_str, '(i2.2)') domnum
   allocate (cdl_levels(lz))
   cdl_levels = lprs(lz:1:-1)*0.01

! lets make the fua file first (contains 3d variables).
   if (trim(mtype) /= 'st4') then

! write out the 3d stuff using laps library routine

      write (6, *) ' write_to_lapsdir = ', write_to_lapsdir

      if (.not. write_to_lapsdir) then
         output_dir = trim(lfmprd_dir)//'/d'//domnum_str//'/fua/'
      else
         output_dir = trim(laps_data_root)//'/lapsprd/fua/'//trim(mtype)//'/'
      end if
      cdl_dir = trim(laps_data_root)//'/cdl/'

! build the output file name so we can create a "donefile" if
! running in realtime mode.

      call make_fnam_lp(laps_reftime, gtime, istatus)

      call make_fcst_time(laps_valtime, laps_reftime, fcst_hhmm, istatus)

      call cvt_fname_v3(output_dir, gtime, fcst_hhmm, 'fua', 3, output_file &
                        , fnlen, istatus)
!beka

      i4_elapsed = ishow_timer()

      print *, ' '
      print *, 'writing 3d fields to netcdf: ', trim(output_file)
      inquire (file=trim(output_dir), exist=outdir_defined)
      if (outdir_defined .eqv. .false.) then
         write (6, *) ' error: output directory does not exist'
         stop
      end if
      lvls3d = lvls3d*0.01
      call write_laps_lfm(laps_reftime, laps_valtime, trim(output_dir), trim(cdl_dir) &
                          , 'fua' &
                          , lx, ly, lz*nvar3dout, lz*nvar3dout, name3d, lvls3d &
                          , lvltype3d, units3d, com3d, lz, cdl_levels &
                          , pgrid, istatus)

      if (istatus /= 1) then
         print *, 'error writing laps 3d (fua) netcdf file.'
      elseif (verbose) then
         print *, 'done writing 3d netcdf file.'
      end if

!beka
      i4_elapsed = ishow_timer()

      if (realtime .and. istatus == 1 .and. make_donefile) then
         donefile = trim(output_file)//'.done'
         open (unit=2, file=donefile, status='unknown')
         close (2)
      end if
      deallocate (cdl_levels)

   end if

! do 2d variables.

   lvls2d(:) = 0.

   if (.not. write_to_lapsdir) then
      output_dir = trim(lfmprd_dir)//'/d'//domnum_str//'/fsf/'
   else
      output_dir = trim(laps_data_root)//'/lapsprd/fsf/'//trim(mtype)//'/'
   end if

   cdl_dir = trim(laps_data_root)//'/cdl/'
   call make_fnam_lp(laps_reftime, gtime, istatus)

   call make_fcst_time(laps_valtime, laps_reftime, fcst_hhmm, istatus)

   call cvt_fname_v3(output_dir, gtime, fcst_hhmm, 'fsf', 3, output_file &
                     , fnlen, istatus)

   print *, 'writing 2d fields to netcdf: ', trim(output_file)
   inquire (file=trim(output_dir), exist=outdir_defined)
   if (outdir_defined .eqv. .false.) then
      write (6, *) ' error: output directory does not exist'
      stop
   end if

   call write_laps_lfm(laps_reftime, laps_valtime, trim(output_dir), cdl_dir &
                       , 'fsf' &
                       , lx, ly, nvar2dout, nvar2dout, name2d, lvls2d, lvltype2d &
                       , units2d, com2d, 1, 0. &
                       , sgrid, istatus)

   if (istatus /= 1) then
      print *, 'error writing laps 2d (fsf) netcdf file.'
   elseif (verbose) then
      print *, 'done writing 2d netcdf file.'
   end if

   if (realtime .and. istatus == 1 .and. make_donefile) then
      donefile = trim(output_file)//'.done'
      open (unit=2, file=donefile, status='unknown')
      close (2)
   end if

   return
end

!===============================================================================

subroutine grib_sfc_vars(laps_reftime, laps_valtime)

   use lfmgrid, only: precip_dt, lx, ly, nvar2dout, gribit, sgrid, verbose, rmsg, &
                      com2d, timerange, table_version, center_id, subcenter_id, &
                      process_id, param, leveltype, level1, level2, scalep10, &
                      igds, funit, startbyte, nbytes
   use grib

   implicit none

   integer :: laps_reftime, laps_valtime, istatus &
              , id(27), yyyyr, mmr, ddr, hhr, minr, itype &
              , timeunit, timeperiod1, timeperiod2 &
              , fcsttime_now, fcsttime_prev, itot &
              , n, shape(1)

   real, pointer, dimension(:, :) :: fld2d
   real, allocatable, dimension(:) :: fld

   character(len=24) :: atime
   character(len=3) :: amonth, amonths(12)

   data amonths/'jan', 'feb', 'mar', 'apr', 'may', 'jun' &
      , 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/

! compute year, month, day of month, hour, and minute from laps_reftime.

   call cv_i4tim_asc_lp(laps_reftime, atime, istatus)
   read (atime, '(i2.2,x,a3,x,i4.4,x,i2.2,x,i2.2)') ddr, amonth, yyyyr, hhr, minr

   do mmr = 1, 12
      if (amonth == amonths(mmr)) exit
   end do

! determine appropriate timeunit.

   if (mod(precip_dt, 3600) == 0) then
! time unit in hours.
      timeunit = 1
      fcsttime_now = (laps_valtime - laps_reftime)/3600
      if (fcsttime_now > 0) then
         fcsttime_prev = fcsttime_now - (precip_dt/3600)
      else
         fcsttime_prev = 0
      end if
   else
! time unit in minutes.
      timeunit = 0
      fcsttime_now = (laps_valtime - laps_reftime)/60
      if (fcsttime_now > 0) then
         fcsttime_prev = fcsttime_now - (precip_dt/60)
      else
         fcsttime_prev = 0
      end if
   end if
   itype = 0

! grib up each variable.

   allocate (fld(lx*ly))
   shape(1) = lx*ly
   do n = 1, nvar2dout
      if (gribit(n)) then
         fld2d => sgrid(1:lx, 1:ly, n)
         if (minval(fld2d) >= rmsg) cycle
         fld = reshape(fld2d, shape)
         if (verbose) write (6, '('' gribbing '',a25,'', min/max = '',2(f10.4,1x))') &
            com2d(n) (1:25), minval(fld), maxval(fld)
         if (timerange(n) == 0) then
            timeperiod1 = fcsttime_now
            timeperiod2 = 0
         elseif (timerange(n) == 4) then
            timeperiod1 = fcsttime_prev
            timeperiod2 = fcsttime_now
         end if
         call make_id(table_version, center_id, subcenter_id, process_id &
                      , param(n), leveltype(n), level1(n), level2(n), yyyyr, mmr, ddr &
                      , hhr, minr, timeunit, timerange(n), timeperiod1, timeperiod2 &
                      , scalep10(n), id)
         call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
         nbytes = nbytes + itot
         startbyte = nbytes + 1
      end if
   end do
   deallocate (fld)

   return
end

!===============================================================================

subroutine grib_ua_vars(laps_reftime, laps_valtime)

   use lfmgrid, only: precip_dt, lx, ly, lz, lprs, nvar3dout, gribitua, pgrid, &
                      rmsg, verbose, com3d, table_version, center_id, subcenter_id, &
                      process_id, paramua, scalep10, igds, funit, startbyte, nbytes
   use grib

   implicit none

   integer :: laps_reftime, laps_valtime, istatus &
              , id(27), yyyyr, mmr, ddr, hhr, minr &
              , itype, lvltype, lvl1, lvl2, tmrange &
              , timeunit, timeperiod1, timeperiod2 &
              , fcsttime_now, fcsttime_prev, itot &
              , k, n, shape(1)

   real, pointer, dimension(:, :) :: fld2d
   real, allocatable, dimension(:) :: fld

   character(len=24) :: atime
   character(len=3) :: amonth, amonths(12)

   data amonths/'jan', 'feb', 'mar', 'apr', 'may', 'jun' &
      , 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/

! compute year, month, day of month, hour, and minute from laps_reftime

   call cv_i4tim_asc_lp(laps_reftime, atime, istatus)
   read (atime, '(i2.2,x,a3,x,i4.4,x,i2.2,x,i2.2)') ddr, amonth, yyyyr, hhr, minr

   do mmr = 1, 12
      if (amonth == amonths(mmr)) exit
   end do

! determine appropriate timeunit.

   if (mod(precip_dt, 3600) == 0) then
! time unit in hours.
      timeunit = 1
      fcsttime_now = (laps_valtime - laps_reftime)/3600
      if (fcsttime_now > 0) then
         fcsttime_prev = fcsttime_now - (precip_dt/3600)
      else
         fcsttime_prev = 0
      end if
   else
! time unit in minutes.
      timeunit = 0
      fcsttime_now = (laps_valtime - laps_reftime)/60
      if (fcsttime_now > 0) then
         fcsttime_prev = fcsttime_now - (precip_dt/60)
      else
         fcsttime_prev = 0
      end if
   end if

! grib up each variable at each level.

   itype = 0
   lvltype = 100
   lvl2 = 0
   tmrange = 0
   timeperiod1 = fcsttime_now
   timeperiod2 = 0

   allocate (fld(lx*ly))
   shape(1) = lx*ly
   do k = 1, lz
      lvl1 = nint(lprs(k))*0.01
      if (lvl1 <= 1000) then
         do n = 1, nvar3dout
            if (gribitua(n)) then
               fld2d => pgrid(1:lx, 1:ly, (n - 1)*lz + k)
               if (minval(fld2d) >= rmsg) cycle
               fld = reshape(fld2d, shape)
               if (verbose) write (6, '('' gribbing '',a25,'', level='',i4,''mb, min/max = '',2(f10.4,1x))') &
                  com3d((n - 1)*lz + k) (1:25), lvl1, minval(fld), maxval(fld)
               call make_id(table_version, center_id, subcenter_id, process_id &
                            , paramua(n), lvltype, lvl1, lvl2, yyyyr, mmr, ddr &
                            , hhr, minr, timeunit, tmrange, timeperiod1, timeperiod2 &
                            , scalep10(n), id)
               call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
               nbytes = nbytes + itot
               startbyte = nbytes + 1
            end if
         end do
      end if
   end do
   deallocate (fld)

   return
end
