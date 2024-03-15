!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! program:         laps2grib
!!
!! purpose:         converts a subset of laps analysis grids to a grib2 file
!!
!! author:        brent shaw, weathernews inc.
!!
!! history:
!!                07 dec 2006        brent shaw        initial version
!!                25 nov 2011         paula mccaslin        modified to add accumulated varaibles
!!                15 dec 2011         paula mccaslin        modified to add model varaibles (from fsf, fua)
!!                28 may 2014            craig hartsough        modified to add multiple output
!!                                              option
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program laps2grib

  !! requires some modules in laps_src_root/lib/modules
   use grib2
  !! local modules
   use laps2grib_config
   use lapsdata
   include 'lapsparms.cmn'     ! laps_cycle_time, model_cycle_time, model_fcst_intvl

   implicit none

   integer                     :: istatus, year, month, day, hour, minute, second
   integer                     :: filestatus, g2lun
   integer                     :: ncid, varid, ncstat
   integer                     :: k, n
   integer, parameter          :: data_type = 0
   integer                     :: reftime_sig !0=analysis, 1=start of fcst
   integer, parameter          :: process_type = 0
   integer, parameter          :: bg_process_id = 255
   integer, parameter          :: cutoff_hr = 0
   integer, parameter          :: cutoff_min = 0
   integer                     :: time_unit_indicator
   integer                     :: time_range
   integer, parameter          :: levtype_iso = 100 !fixed level type
   integer, parameter          :: levscale_iso = 0
   integer, parameter          :: g2miss = 255
   integer, parameter          :: newrec = 1
   integer, parameter          :: pack_method = 2
   integer                     :: miss_mgmt, inomiss
   integer                     :: val_time, ref_time, accum_time
   integer                     :: eyear, emonth, eday, ehour, eminute, esecond
   integer                     :: etime_value, etime_unit
   integer                     :: igrid

   real, allocatable           :: lapsdata2d(:, :), lapsdata3d(:, :, :), slab(:, :)
   real                        :: r_missing, rhhmm
   character(len=32), parameter :: vtag = "laps2grib v1.2, 01 dec 2011, wni" !"laps2grib v1.0, 08 dec 2006, wni"
   character(len=256)          :: laps_data_root, lapsfile
   character(len=256)          :: g2file, g2file_tmp
   character(len=512)          :: syscmd

   integer                     :: i, iargc, max_args, ihhmm
   integer                     :: modeltime_passed_in
   character(len=100)          :: vtab, forecast_id
   character(len=100)          :: vtab2
   character(len=5)            :: hhmm
   character(len=14)           :: file_a9time, modeltime
   logical                     :: dir_exists

   ! print banner
   print *, "=================================================================="
   print *, "********** ", vtag, " **********"
   print *, ""
   print *, " usage:       laps2grib.exe [vtab]"
   print *, " model usage: laps2grib.exe vtab [hh]hmm forecast_id [modeltime]"
   print *, "        (e.g. laps2grib.exe wfr2grib.vtab 1200 wrf-hrrr 130451800)"
   print *, "        (optional modeltime has the format yyjjjhhmm)"
   print *, "=================================================================="

   vtab = 'laps2grib.vtab'
   vtab2 = 'laps2grib.vtab2'
   reftime_sig = 0
   forecast_id = ''
   hhmm = ''
   max_args = iargc()

   if (max_args .eq. 1) then
      ! expect vtable name
      call getarg(1, vtab)

   else if ((max_args .eq. 3) .or. (max_args .eq. 4)) then
      reftime_sig = 1
      ! expect vtable name
      call getarg(1, vtab)

      ! expect hhmm for forecast time
      call getarg(2, hhmm)
      read (hhmm, *, iostat=ihhmm) rhhmm
      if (rhhmm .eq. 0. .and. ihhmm .ne. 0) then
         print *, "the usage requires hhmm for the forecast time, this is not an expected number: ", hhmm
         stop
      else if (rhhmm .eq. 0.) then
         hhmm = '0000'
      else if (len_trim(hhmm) .lt. 3) then
         print *, "the usage requires hhmm for the forecast time, the string is too short: ", hhmm
         stop
      else if (len_trim(hhmm) .le. 3 .and. index(hhmm, ' ') .ne. 0) then
         print *, "--> index value", index(hhmm, ' ')
         hhmm = '0'//hhmm
      end if

      ! expect ensemble forecast_id name, .e.g. mean, or wrf-hrrr
      call getarg(3, forecast_id)
      forecast_id = '/'//forecast_id
      modeltime_passed_in = 0
      if (max_args .eq. 4) then
         call getarg(4, modeltime)
         if (len_trim(modeltime) .ne. 9) then
            print *, 'modeltime passed in is not the correct length: ', trim(modeltime)
            print *, 'defaulting to time in file modeltime.dat'
         else
            modeltime_passed_in = 1
         end if
      end if
   else if (max_args .ne. 0) then
      stop "check the usage statement above..."
   end if

   ! get the laps_data_root
   call getenv('laps_data_root', laps_data_root)
   if (len_trim(laps_data_root) .lt. 1) then
      print *, "laps_data_root not set!"
      stop 'no laps_data_root'
   end if

   print *, "-- laps_data_root=", trim(laps_data_root)

   ! read the laps2grib namelist
   call read_laps2grib_nl(laps_data_root)

   if (.not. lrun_laps2grib) then
      write (6, *) 'lrun_laps2grib is set to false - stopping program'
      stop
   end if

   ! set up projection
   call get_laps_proj(laps_data_root)
   call get_laps_plevels(laps_data_root)

   ! get laps analysis valid time
   if (len_trim(hhmm) .eq. 0) then
      call get_laps_analtime
   else
      call get_laps_modeltime(modeltime, modeltime_passed_in)
   end if
   call cv_i4tim_int_lp(i4time, year, month, day, hour, minute, second, istatus)
   year = year + 1900
   file_a9time = a9time
   if (hhmm .ne. '') file_a9time = a9time//hhmm
   print *, "-- using timestamp: ", file_a9time

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! run the remainder of the program once for each desired output file  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do igrid = 1, num_grids
      if (igrid .gt. 1) then
         vtab = vtab2
         output_path = output_path2
      end if

      ! configure the variable list
      call get_data_config(laps_data_root, vtab)

      ! get the rmissing value
      call get_r_missing_data(r_missing, istatus)
      print *, "-- r_missing: ", r_missing

      ! check the grib dir
      g2file = trim(output_path)//trim(forecast_id)
      inquire (file=g2file, exist=dir_exists)
      if (.not. dir_exists) then
         print *, "dir does not exist: ", trim(g2file)
         stop
      end if

      ! open the grib file
      g2file = trim(output_path)//trim(forecast_id)//'/'//trim(file_a9time)//'.gr2'
      g2file_tmp = trim(output_path)//trim(forecast_id)//'/.'//trim(file_a9time)//'.gr2'
      call init_grib2_file(g2file_tmp, laps_proj, center_id, subcenter_id, &
                           reftime_sig, year, month, day, hour, minute, second, &
                           prod_status, data_type, g2lun, istatus)

      ! process any 3d isobaric variables
      if (n_iso .gt. 0) then
         allocate (lapsdata3d(nx, ny, nz))
         allocate (slab(nx, ny))
         print *, "-- processing 3d isobaric fields"
         loop3d: do n = 1, n_iso

            if (meta_iso3d(n)%qbal .eq. 0) then
               lapsfile = trim(laps_data_root)//'/lapsprd/'//meta_iso3d(n)%ext//trim(forecast_id)//'/'//trim(file_a9time)//'.'// &
                          meta_iso3d(n)%ext
            else
               lapsfile = trim(laps_data_root)//'/lapsprd/balance/'//meta_iso3d(n)%ext//'/'// &
                          trim(file_a9time)//'.'//meta_iso3d(n)%ext
            end if

            filestatus = 0
            call ncread_3d(lapsfile, nx, ny, nz, meta_iso3d(n)%var, lapsdata3d, val_time, ref_time, istatus)
            if (istatus .eq. 1) then
               filestatus = 1  ! ncread_3d returned data...ok to write out grib2
               call calc_val_ref_times(val_time, ref_time, time_unit_indicator, time_range, istatus)
               if (istatus .ne. 1) then
                  filestatus = 0  ! could not calc val_ref times...do not write out grib2
                  print *, "could not calc val_ref times: ", meta_iso3d(n)%ext, " var ", meta_iso3d(n)%var
               end if
            else
               print *, "could not read file: ", meta_iso3d(n)%ext, " var ", meta_iso3d(n)%var
            end if

            if (filestatus .eq. 1) then
               press_loop: do k = 1, nz
                  if ((plevels_pa(k) .ge. meta_iso3d(n)%pmin) .and. &
                      (plevels_pa(k) .le. meta_iso3d(n)%pmax)) then

                     slab = lapsdata3d(:, :, k)
                     if (maxval(slab) .eq. r_missing) then
                        inomiss = 0
                        miss_mgmt = 1
                        if (meta_iso3d(n)%conv_fac .ne. 1.) then
                           where (slab .ne. r_missing) slab = slab*meta_iso3d(n)%conv_fac
                        end if

                     else
                        inomiss = 1
                        miss_mgmt = 0
                        if (meta_iso3d(n)%conv_fac .ne. 1.) then
                           slab = slab*meta_iso3d(n)%conv_fac
                        end if
                     end if
                     print '("grib: ",a3,1x,a3,1x,i5,"mb",i3,i3,i4,f12.4,f12.4)', meta_iso3d(n)%ext, &
                        meta_iso3d(n)%var, nint(plevels_pa(k)*0.01), &
                        meta_iso3d(n)%discipline, meta_iso3d(n)%cat_id, meta_iso3d(n)%parm_id, &
                        minval(slab), maxval(slab)
                     call write_grib2_template0(g2lun, meta_iso3d(n)%discipline, &
                                                meta_iso3d(n)%cat_id, meta_iso3d(n)%parm_id, process_type, &
                                                bg_process_id, process_id, cutoff_hr, cutoff_min, time_unit_indicator, &
                                                time_range, levtype_iso, levscale_iso, nint(plevels_pa(k)), g2miss, g2miss, &
                                                g2miss, pack_method, meta_iso3d(n)%scale_fac, miss_mgmt, &
                                                nx, ny, newrec, inomiss, r_missing, r_missing, slab)
                  else
                     print '("skip: ",a3,1x,a3,1x,i5,"mb")', meta_iso3d(n)%ext, meta_iso3d(n)%var, &
                        nint(plevels_pa(k)*0.01)
                  end if
               end do press_loop

            else
               print *, "problem getting: ", trim(forecast_id), " ext .", meta_iso3d(n)%ext, " var ", meta_iso3d(n)%var
            end if
         end do loop3d
         deallocate (lapsdata3d)
         deallocate (slab)
      end if

      ! process any 2d variables
      if (n_2d .gt. 0) then
         allocate (lapsdata2d(nx, ny))
         print *, "-- processing 2d fields"
         loop2d: do n = 1, n_2d
            if (meta_2d(n)%ext .ne. 'n7g') then
               lapsfile = trim(laps_data_root)//'/lapsprd/'//meta_2d(n)%ext//trim(forecast_id)//'/'//trim(file_a9time)//'.'// &
                          meta_2d(n)%ext

            else
               lapsfile = trim(laps_data_root)//'/static/static.nest7grid'
            end if

            filestatus = 0
            call ncread_2d(lapsfile, nx, ny, meta_2d(n)%var, lapsdata2d, val_time, ref_time, istatus)
            if (istatus .eq. 1) then
               filestatus = 1  ! ncread_2d returned data...ok to write out grib2
               call calc_val_ref_times(val_time, ref_time, time_unit_indicator, time_range, istatus)
               if (istatus .ne. 1) then
                  filestatus = 0  ! could not calc val_ref times...do not write out grib2
                  print *, "could not calc val_ref times: ", meta_2d(n)%ext, " var ", meta_2d(n)%var
               end if
            else
               print *, "could not read file: ", meta_2d(n)%ext, " var ", meta_2d(n)%var
            end if

            if (filestatus .eq. 1) then
               if (maxval(lapsdata2d) .eq. r_missing) then
                  inomiss = 0
                  miss_mgmt = 1
                  if (meta_2d(n)%conv_fac .ne. 1.) then
                     where (lapsdata2d .ne. r_missing) lapsdata2d = lapsdata2d*meta_2d(n)%conv_fac
                  end if

               else
                  inomiss = 1
                  miss_mgmt = 0
                  if (meta_2d(n)%conv_fac .ne. 1.) then
                     lapsdata2d = lapsdata2d*meta_2d(n)%conv_fac
                  end if
               end if

               print '("grib: ",a3,1x,a3,1x,f12.4,f12.4)', &
                  meta_2d(n)%ext, meta_2d(n)%var, minval(lapsdata2d), maxval(lapsdata2d)
               call write_grib2_template0(g2lun, meta_2d(n)%discipline, &
                                          meta_2d(n)%cat_id, meta_2d(n)%parm_id, process_type, &
                                          bg_process_id, process_id, cutoff_hr, cutoff_min, time_unit_indicator, &
                                          time_range, meta_2d(n)%lev1_type, meta_2d(n)%lev1_scale, meta_2d(n)%lev1_value, &
                                          meta_2d(n)%lev2_type, meta_2d(n)%lev2_scale, meta_2d(n)%lev2_value, &
                                          pack_method, meta_2d(n)%scale_fac, miss_mgmt, &
                                          nx, ny, newrec, inomiss, r_missing, r_missing, lapsdata2d)

            else
               print *, "problem getting: ", trim(forecast_id), " ext .", meta_2d(n)%ext, " var ", meta_2d(n)%var
            end if
         end do loop2d
         deallocate (lapsdata2d)
      end if

      ! process any accumulated 2d variables
      if (n_accum2d .gt. 0) then
         allocate (lapsdata2d(nx, ny))
         print *, "-- processing accumulated 2d fields"
         loop_a2d: do n = 1, n_accum2d
            lapsfile = trim(laps_data_root)//'/lapsprd/'//meta_accum2d(n)%ext//trim(forecast_id)//'/'//trim(file_a9time)//'.'// &
                       meta_accum2d(n)%ext

            filestatus = 0
            call ncread_2d(lapsfile, nx, ny, meta_accum2d(n)%var, lapsdata2d, val_time, ref_time, istatus)
            print *, "ncread_2d(", lapsfile, nx, ny, n, val_time, ref_time, istatus, ")"
            if (istatus .eq. 1) then
               filestatus = 1  ! ncread_2d returned data...ok to write out grib2
               call calc_val_ref_times(val_time, ref_time, time_unit_indicator, time_range, istatus)
               if (istatus .ne. 1) then
                  filestatus = 0  ! could not calc val_ref times...do not write out grib2
                  print *, "could not calc val_ref times: ", meta_accum2d(n)%ext, " var ", meta_accum2d(n)%var
               end if
            else
               print *, "could not read file: ", meta_accum2d(n)%ext, " var ", meta_accum2d(n)%var
            end if

            if (filestatus .eq. 1) then
               if (maxval(lapsdata2d) .eq. r_missing) then
                  inomiss = 0
                  miss_mgmt = 1
                  if (meta_accum2d(n)%conv_fac .ne. 1.) then
                     where (lapsdata2d .ne. r_missing) lapsdata2d = lapsdata2d*meta_accum2d(n)%conv_fac
                  end if

               else
                  inomiss = 1
                  miss_mgmt = 0
                  if (meta_accum2d(n)%conv_fac .ne. 1.) then
                     lapsdata2d = lapsdata2d*meta_accum2d(n)%conv_fac
                  end if
               end if

               call calc_accum_time(lapsfile, meta_accum2d(n)%ext, meta_accum2d(n)%var, val_time, ref_time, accum_time, &
                                    time_unit_indicator, time_range, etime_unit, etime_value, istatus)
               if (istatus .eq. 0) print *, "problem with calc_accum_time"

               print *, 'after calc_accum_time ', ref_time, ', ', ref_time, ' ', &
                  accum_time
               print *, 'etime: ', etime_unit, ', ', etime_value

               ! octet 49-50, set to current runtime
               meta_accum2d(n)%etime_unit = etime_unit
               meta_accum2d(n)%etime_value = etime_value

               ! cv_i4tim_int_lp expects laps time  klh - 20 mar 2014
               val_time = val_time + 315619200
               ! --- valid date, time to describe data for grib2.
               call cv_i4tim_int_lp(val_time, eyear, emonth, eday, ehour, eminute, esecond, istatus)
               !print *, 'cv_i4tim_int_lp(', val_time,eyear,emonth,eday,ehour,&
               !eminute,esecond,istatus, ')'
               ! add 10 to eyear to have time in grib2 files to account for 10 yr difference between i4time and unixtime
               ! octet 35-40, set to current runtime
               !meta_accum2d(n)%eyear = eyear + 1900 + 10
               ! this 10 year offset is not needed when valid laps time is passed to cv_itim_int_lp   klh -- 20 mar 2014

               meta_accum2d(n)%eyear = eyear + 1900
               meta_accum2d(n)%emon = emonth
               meta_accum2d(n)%eday = eday
               meta_accum2d(n)%ehour = ehour
               meta_accum2d(n)%emin = 0

               print '("grib: ",a3,1x,a3,1x,f12.4,f12.4)', &
                  meta_accum2d(n)%ext, meta_accum2d(n)%var, minval(lapsdata2d), maxval(lapsdata2d)
               call write_grib2_template8(g2lun, meta_accum2d(n)%discipline, &
                                          meta_accum2d(n)%cat_id, meta_accum2d(n)%parm_id, process_type, &
                                          bg_process_id, process_id, cutoff_hr, cutoff_min, time_unit_indicator, &
                                    time_range, meta_accum2d(n)%lev1_type, meta_accum2d(n)%lev1_scale, meta_accum2d(n)%lev1_value, &
                                          meta_accum2d(n)%lev2_type, meta_accum2d(n)%lev2_scale, meta_accum2d(n)%lev2_value, &
                                          meta_accum2d(n)%eyear, meta_accum2d(n)%emon, meta_accum2d(n)%eday, &
                                          meta_accum2d(n)%ehour, meta_accum2d(n)%emin, meta_accum2d(n)%esec, &
                                          meta_accum2d(n)%ntimes, meta_accum2d(n)%ntimes_miss, &
                                          meta_accum2d(n)%stattype, meta_accum2d(n)%periodtype, &
                                          meta_accum2d(n)%etime_unit, meta_accum2d(n)%etime_value, &
                                          pack_method, meta_accum2d(n)%scale_fac, miss_mgmt, &
                                          nx, ny, newrec, inomiss, r_missing, r_missing, lapsdata2d)
            else
               print *, "problem getting: ", trim(forecast_id), " ext .", meta_accum2d(n)%ext, " var ", meta_accum2d(n)%var
            end if
         end do loop_a2d
         deallocate (lapsdata2d)
      end if

      ! close the file and rename it to the final output name
      print *, "-- closing file"
      call close_grib2_file(g2lun)
      syscmd = 'mv '//trim(g2file_tmp)//' '//trim(g2file)
      print *, trim(syscmd)
      call system(syscmd)

   end do  ! igrid

   print *, "======================================================"
   print *, "**********      laps2grib completed         **********"
   print *, "======================================================"
end program laps2grib

!!!!!!!!!!!!!!!!!!!!!!!
subroutine ncread_3d(ncfile, nx, ny, nz, ncvar, data3d, nc_val, nc_ref, istatus)

   implicit none
   character(len=*), intent(in)  :: ncfile
   integer, intent(in)           :: nx, ny, nz
   character(len=3), intent(in)  :: ncvar
   real, intent(out)              :: data3d(nx, ny, nz)
   integer, intent(out)          :: istatus

   integer  :: ncid, vid, ncstat
   integer  :: nc_val, nc_ref
   logical  :: file_exists
   include 'netcdf.inc'

   istatus = 1
   inquire (file=ncfile, exist=file_exists)
   if (.not. file_exists) then
      print *, "file not found: ", trim(ncfile)
      istatus = 0
      return
   end if
   ncstat = nf_open(ncfile, nf_nowrite, ncid)
   if (ncstat .ne. nf_noerr) then
      print *, "error opening: ", trim(ncfile)
      istatus = 0
      return
   end if
   ncstat = nf_inq_varid(ncid, ncvar, vid)
   if (ncstat .eq. nf_noerr) then
      ncstat = nf_get_var_real(ncid, vid, data3d)
      if (ncstat .ne. nf_noerr) then
         print *, "problem getting data for ", trim(ncvar)
         istatus = 0
      end if
   else
      print *, "could not find var id for ", trim(ncvar)
      istatus = 0
   end if

   ! looking valtime value
   ncstat = nf_inq_varid(ncid, 'valtime', vid)
   if (ncstat .eq. nf_noerr) then
      ncstat = nf_get_var_int(ncid, vid, nc_val)
      if (ncstat .ne. nf_noerr) then
         print *, "could not find value for var id for valtime", vid
         istatus = 0
      end if
   else
      print *, "could not find var id for valtime"
   end if

   ! looking reftime value
   ncstat = nf_inq_varid(ncid, 'reftime', vid)
   if (ncstat .eq. nf_noerr) then
      ncstat = nf_get_var_int(ncid, vid, nc_ref)
      if (ncstat .ne. nf_noerr) then
         print *, "could not find value for var id for reftime", vid
         istatus = 0
      end if
   else
      print *, "could not find var id for reftime"
   end if

   ncstat = nf_close(ncid)
   return
end subroutine ncread_3d
!!!!!!!!!!!!!!!!!!!!!!!
subroutine ncread_2d(ncfile, nx, ny, ncvar, data2d, nc_val, nc_ref, istatus)

   implicit none
   character(len=*), intent(in)  :: ncfile
   integer, intent(in)           :: nx, ny
   character(len=3), intent(in)  :: ncvar
   real, intent(out)              :: data2d(nx, ny)
   integer, intent(out)          :: istatus

   integer  :: ncid, vid, ncstat
   integer  :: nc_val, nc_ref
   logical  :: file_exists
   include 'netcdf.inc'

   istatus = 1
   inquire (file=ncfile, exist=file_exists)
   if (.not. file_exists) then
      print *, "file not found: ", trim(ncfile)
      istatus = 0
      return
   end if
   ncstat = nf_open(ncfile, nf_nowrite, ncid)
   if (ncstat .ne. nf_noerr) then
      print *, "error opening: ", trim(ncfile)
      istatus = 0
      return
   end if
   ncstat = nf_inq_varid(ncid, ncvar, vid)
   if (ncstat .eq. nf_noerr) then
      ncstat = nf_get_var_real(ncid, vid, data2d)
      if (ncstat .ne. nf_noerr) then
         print *, "problem getting data for ", trim(ncvar)
         istatus = 0
      end if
   else
      print *, "could not find var id for ", trim(ncvar)
      istatus = 0
   end if

   ! looking valtime value
   ncstat = nf_inq_varid(ncid, 'valtime', vid)
   if (ncstat .eq. nf_noerr) then
      ncstat = nf_get_var_int(ncid, vid, nc_val)
      if (ncstat .ne. nf_noerr) then
         print *, "could not find value for var id for valtime", vid
         istatus = 0
      end if
   else
      print *, "could not find var id for valtime"
   end if

   ! looking reftime value
   ncstat = nf_inq_varid(ncid, 'reftime', vid)
   if (ncstat .eq. nf_noerr) then
      ncstat = nf_get_var_int(ncid, vid, nc_ref)
      if (ncstat .ne. nf_noerr) then
         print *, "could not find value for var id for reftime", vid
         istatus = 0
      end if
   else
      print *, "could not find var id for reftime"
   end if

   ncstat = nf_close(ncid)
   return
end subroutine ncread_2d

!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_val_ref_times(nc_val, nc_ref, time_unit_indicator, time_range, istatus)

   include 'lapsparms.cmn' ! laps_cycle_time, model_cycle_time, model_fcst_intvl
   implicit none
   integer, intent(out)    :: istatus

   integer  :: nc_val, nc_ref, dif
   integer  :: remainder, my_cycle_time
   integer  :: time_unit_indicator, time_range

   istatus = 1
   ! --- set accum time arg with the use of current time intervals
   dif = nc_val - nc_ref
   my_cycle_time = laps_cycle_time
   if (dif .ne. 0) my_cycle_time = dif

   remainder = mod(my_cycle_time, 3600)
   if (remainder .eq. 0) then
      ! working with time interval of hours
      time_unit_indicator = 1
      time_range = my_cycle_time/3600 !determine no. hrs in laps cycle
      if (dif .eq. 0) time_range = 0
   else
      ! working with time interval of minutes
      remainder = mod(my_cycle_time, 60)
      if (remainder .eq. 0) then
         time_unit_indicator = 0
         time_range = my_cycle_time/60 !determine no. mins in laps cycle
      else
         print *, "problem setting time_unit and time_range"
         istatus = 0
      end if
   end if

   !print *," ---->> valtime", nc_val, " reftime", nc_ref
   return

end subroutine calc_val_ref_times

!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_accum_time(ncfile, ncext, ncvar, nc_val, nc_ref, accum_i4time, time_unit_indicator, startime_accum, &
                           etim_unit, etim_value, istatus)

   include 'lapsparms.cmn'     ! laps_cycle_time, model_cycle_time, model_fcst_intvl
   implicit none
   character(len=*), intent(in)  :: ncfile
   character(len=3), intent(in)  :: ncvar, ncext
   integer, intent(out)          :: istatus
   integer  :: nc_val, nc_ref, dif
   integer  :: accum_i4time
   character(len=100) :: str

   integer  :: eyear, emonth, eday, ehour, eminute, esecond
   integer  :: etim_value, etim_unit, remainder, my_cycle_time
   integer  :: time_unit_indicator, startime_accum
   integer  :: ncid, vid, ncstat
   logical  :: file_exists
   include 'netcdf.inc'

   istatus = 1
   if (ncext .ne. 'fsf') then

      inquire (file=ncfile, exist=file_exists)
      if (.not. file_exists) then
         print *, "file not found: ", trim(ncfile)
         istatus = 0
         return
      end if
      ncstat = nf_open(ncfile, nf_nowrite, ncid)
      if (ncstat .ne. nf_noerr) then
         print *, "error opening: ", trim(ncfile)
         istatus = 0
         return
      end if

      ! looking for accumulation start time (a9time) in sto_comment, or rto_comment
      ncstat = nf_inq_varid(ncid, trim(ncvar)//'_comment', vid)
      if (ncstat .eq. nf_noerr) then
         ncstat = nf_get_var_text(ncid, vid, str)
         if (ncstat .ne. nf_noerr) then
            print *, "could not find value for nc _comment: ", vid
            istatus = 0
         end if
      else
         print *, "could not find var id for ", trim(ncvar)//'_comment'
      end if
      ncstat = nf_close(ncid)

      ! convert a9time to i4time
      read (str, '(a9)') str
      call i4time_fname_lp(str, accum_i4time, istatus)
      call cv_i4tim_int_lp(accum_i4time, eyear, emonth, eday, ehour, eminute, esecond, istatus)

      ! --- set accum time arg with the use of current time intervals
      accum_i4time = accum_i4time - 315619200
      dif = nc_val - accum_i4time

   else
      !do not need var id for fsf, fua
      dif = nc_val - nc_ref
   end if

   my_cycle_time = laps_cycle_time
   if (dif .ne. 0) my_cycle_time = dif

   remainder = mod(my_cycle_time, 3600)
   if (remainder .eq. 0) then
      ! working with time interval of hours
      etim_value = my_cycle_time/3600
      etim_unit = 1 !units of hrs
      time_unit_indicator = 1
      startime_accum = ehour
      if (dif .eq. 0) startime_accum = 0
   else
      ! working with time interval of minutes
      remainder = mod(my_cycle_time, 60)
      if (remainder .eq. 0) then
         etim_value = my_cycle_time/60
         etim_unit = 0 !units of mins
         time_unit_indicator = 0
         startime_accum = eminute
         if (dif .eq. 0) startime_accum = 0
      else
         print *, "problem setting accum vars: etim_value and etim_unit"
         istatus = 0
      end if
   end if

   return

end subroutine calc_accum_time
