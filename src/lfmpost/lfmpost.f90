!dis
!dis    open source license/disclaimer, forecast systems laboratory
!dis    noaa/oar/fsl, 325 broadway boulder, co 80305
!dis
!dis    this software is distributed under the open source definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    in particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - all modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - if significant modifications or enhancements are made to this
!dis    software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    this software and its documentation are in the public domain
!dis    and are furnished "as is."  the authors, the united states
!dis    government, its instrumentalities, officers, employees, and
!dis    agents make no warranty, express or implied, as to the usefulness
!dis    of the software and documentation for any purpose.  they assume
!dis    no responsibility (1) for the use of the software and
!dis    documentation; or (2) to provide technical support to users.
!dis
!dis

program lfmpost

! f90 program to post-process local forecast model output files into
!  laps-format (netcdf fsf/fua), vis5d, grib, etc.  currently supports
!  mm5v3 and wrfv1.2.1 for input model data.
!
!  this program assumes you are processing model using a grid configuration
!  identical to the grid in the specified laps_data_root!
!
! brent shaw, sep 2002
!   base on original mm5post program by same author
!
!

   use setup
   use mm5v3_io
   use wrf_netcdf
   use postproc_lfm
   use constants
   use time_utils
   use vis5d
   use map_utils
   use grib

   implicit none
   integer                     :: lun_data
   integer                     :: t, k
   integer                     :: time_index
   integer                     :: status
   logical                     :: file_ready
   integer                     :: funit, nbytes
   character(len=255)          :: gribdir
   integer                     :: gdir_len, ddir_len, gfile_len
   character(len=255)          :: gribfile, gribdone
   character(len=24)          :: current_time
   integer                     :: laps_reftime
   integer                     :: laps_valtime, laps_valtime_prev
   integer                     :: period_sec
   integer                     :: pass
   real                        :: smth
   real, allocatable           :: rhodry(:, :, :)
   real, allocatable           :: cldliqcon_prs(:, :, :)
   real, allocatable           :: cldicecon_prs(:, :, :)
   real, allocatable           :: raincon_prs(:, :, :)
   real, allocatable           :: snowcon_prs(:, :, :)
   real, allocatable           :: graupelcon_prs(:, :, :)
   ! declarations added for make_points
   integer                     :: year
   integer                     :: month
   integer                     :: day
   integer                     :: hour
   integer                     :: minute
   integer                     :: second
   integer                     :: igds(18)
   integer                     :: ip
   integer                     :: jp
   integer                     :: dir_pt
   integer                     :: spd_pt
   integer                     :: ceiling_pt
   integer                     :: pt
   real                        :: t_pt
   real                        :: td_pt
   integer                     :: rh_pt
   real                        :: u_pt
   real                        :: v_pt
   real                        :: dir_pt_r
   real                        :: spd_pt_r
   real                        :: cld_pt
   real                        :: vis_pt
   real                        :: pcp_pt
   real                        :: snow_pt
   real                        :: vnt_pt
   character(len=16)           :: date_str
   character(len=8)            :: wx_pt
   real, external              :: bint
   real                        :: ceiling_pt_r
   real                        :: u_pt_g
   real                        :: v_pt_g
   real                        :: up_pt_g, vp_pt_g, vp_pt, up_pt, dirp_pt_r, spdp_pt_r, pbl_pt
   integer                     :: dirp_pt, spdp_pt
   integer                     :: t_v5d
   character(len=255)          :: flagfile
   character(len=9)            :: a9time
   character(len=4)            :: fcst_hhmm
   integer                     :: istatus
   integer                     :: startb
   integer                     :: flagunit

   call setup_lfmpost

   call mm5_to_lapstime(cycle_date, laps_reftime)

   ! start the loop over all times

   call s_len(lfmprd_dir, ddir_len)

   write (gribdir, '(a,"/d",i2.2,"/grib")') &
      lfmprd_dir(1:ddir_len), domain_num
   call s_len(gribdir, gdir_len)

   laps_valtime_prev = 0
   time_loop: do t = 1, num_times_to_proc, time_index_inc

      time_index = t + start_time_index - 1
      current_time = times_to_proc(t)
      ! compute laps reference time
      call mm5_to_lapstime(current_time, laps_valtime)
      ! open this times file for reading
      if (mtype .eq. 'mm5') then
         call make_data_file_name(mm5_data_root, domain_num_str, split_output, &
                                  time_index - 1, data_file, file_num3)
         if (split_output) then
            if (realtime) then
               flagfile = trim(data_file)//'.done'
               inquire (file=flagfile, exist=file_ready)
               if (.not. file_ready) then
                  call io_wait(data_file, max_wait_sec)
               end if
            end if
            call open_mm5v3(data_file, lun_data, status)
         else if (t .eq. 1) then
            call open_mm5v3(data_file, lun_data, status)
         end if

      elseif (mtype(1:3) .eq. 'wrf') then
         if (mtype .eq. 'wrf') then
            call make_wrf_file_name(lfmprd_dir, domain_num, sim_tstep(t), data_file)
         elseif (mtype .eq. 'wrf2') then
            call make_wrf2_file_name(lfmprd_dir, domain_num, times_to_proc(t), &
                                     data_file)
         end if
         inquire (file=data_file, exist=file_ready)
         if (.not. file_ready) then
            call wrfio_wait(data_file, max_wait_sec)
         else
            if (realtime) call sleep(60)
         end if
         call open_wrfnc(data_file, lun_data, status)
      end if
      if (status .ne. 0) then
         print *, 'stopping due to error opening ', trim(data_file)
         call abort
      end if

      ! allocate the variables.  allocation/deallocation of most of the arrays
      ! is done inside the time loop in case a user is running this program
      ! concurrently with the model (using split output) on a shared system.

      allocate (psfc(nx, ny))
      allocate (tsfc(nx, ny))
      allocate (thetasfc(nx, ny))
      allocate (thetaesfc(nx, ny))
      allocate (rhsfc(nx, ny))
      allocate (tdsfc(nx, ny))
      allocate (thetaprs(nx, ny, kprs))
      allocate (zprs(nx, ny, kprs))
      allocate (rhprs(nx, ny, kprs))
      allocate (tprs(nx, ny, kprs))
      allocate (tdprs(nx, ny, kprs))
      allocate (shprs(nx, ny, kprs))
      allocate (tkeprs(nx, ny, kprs))
      allocate (redp(nx, ny))
      allocate (pmsl(nx, ny))
      allocate (usfc(nx, ny))
      allocate (uprs(nx, ny, kprs))
      allocate (vsfc(nx, ny))
      allocate (vprs(nx, ny, kprs))
      allocate (upbl(nx, ny))
      allocate (vpbl(nx, ny))
      allocate (wsfc(nx, ny))
      allocate (clwmrsfc(nx, ny))
      allocate (icemrsfc(nx, ny))
      allocate (snowmrsfc(nx, ny))
      allocate (rainmrsfc(nx, ny))
      allocate (graupmrsfc(nx, ny))
      allocate (wprs(nx, ny, kprs))
      allocate (omprs(nx, ny, kprs))
      allocate (cldbase(nx, ny))
      allocate (cldtop(nx, ny))
      allocate (cldamt(nx, ny))
      allocate (ceiling(nx, ny))
      allocate (intliqwater(nx, ny))
      allocate (totpcpwater(nx, ny))
      allocate (max_refl(nx, ny))
      allocate (echo_tops(nx, ny))
      allocate (cldliqmr_prs(nx, ny, kprs))
      allocate (cldicemr_prs(nx, ny, kprs))
      allocate (rainmr_prs(nx, ny, kprs))
      allocate (snowmr_prs(nx, ny, kprs))
      allocate (graupelmr_prs(nx, ny, kprs))
      allocate (refl_prs(nx, ny, kprs))
      allocate (refl_sfc(nx, ny))
      allocate (pcptype_sfc(nx, ny))
      allocate (pcptype_prs(nx, ny, kprs))
      if (.not. allocated(pcp_init)) allocate (pcp_init(nx, ny))
      allocate (pcp_inc(nx, ny))
      allocate (con_pcp_inc(nx, ny))
      if (.not. allocated(con_pcp_init)) allocate (con_pcp_init(nx, ny))
      if (.not. allocated(pcp_tot)) allocate (pcp_tot(nx, ny))
      if (.not. allocated(con_pcp_tot)) allocate (con_pcp_tot(nx, ny))
      if (.not. allocated(snow_init)) allocate (snow_init(nx, ny))
      allocate (snow_inc(nx, ny))
      if (.not. allocated(snow_tot)) allocate (snow_tot(nx, ny))
      if (.not. allocated(srhel)) allocate (srhel(nx, ny))
      allocate (cape(nx, ny))
      allocate (cin(nx, ny))
      allocate (liftedind(nx, ny))
      allocate (visibility(nx, ny))
      allocate (heatind(nx, ny))
      allocate (lwout(nx, ny))
      allocate (swout(nx, ny))
      allocate (lwdown(nx, ny))
      allocate (swdown(nx, ny))
      allocate (albedo(nx, ny))
      allocate (pblhgt(nx, ny))
      allocate (shflux(nx, ny))
      allocate (lhflux(nx, ny))
      allocate (ground_t(nx, ny))
      allocate (vnt_index(nx, ny))
      allocate (ham_index(nx, ny))
      allocate (hah_index(nx, ny))
      allocate (fwi_index(nx, ny))

      if (make_v5d(domain_num)) then
         allocate (abs_vort(nx, ny, kprs))
         allocate (thick_10_5(nx, ny))
         allocate (snowcover(nx, ny))
      end if
      print '("processing time ",a24," file = ",a," unit # = ",i4)', &
         current_time, trim(data_file), lun_data
      call process_one_lfm(lun_data, current_time)
      if (mtype .eq. 'mm5') then
         if (split_output) close (lun_data)
      elseif (mtype(1:3) .eq. 'wrf') then
         if (split_output) call close_wrfnc(lun_data)
      end if

      ! if smoothing is desired, call the smooth routine.  two pass
      ! for 2dx filter on surface fields, 20 pass for 4dx on ua

      if (do_smoothing) then
         print *, 'performing extra smoothing on selected variables...'
         ! extra smoothing for upper-air heights, etc.
         do pass = 1, 10
            call smooth(pmsl, nx, ny, 1, 1.)
            call smooth(redp, nx, ny, 1, 1.)
            call smooth(zprs, nx, ny, kprs, 1.)
            call smooth(srhel, nx, ny, 1, 1.)
            if (make_v5d(domain_num)) then
               call smooth(abs_vort, nx, ny, kprs, 1.)
               call smooth(thick_10_5, nx, ny, 1, 1.)
            end if
         end do
      end if
      ! call output routines.  this is where you can insert a call
      ! to a routine to write whatever format you want

      if (make_laps(domain_num)) then
         print *, 'making laps output files (fua/fsf)...'
         ! for laps, some of the variables need to be manipulated a bit

         where (refl_prs .le. 0.) refl_prs = -10.
         where (refl_sfc .le. 0.) refl_sfc = -10.
         where (max_refl .le. 0.) max_refl = -10.

         ! convert intliqwater and totpcpwater to meters (mks)
         intliqwater = intliqwater*0.001
         totpcpwater = totpcpwater*0.001

         ! convert cloud species from mixing ratios to mass per volume. since
         ! they are mixing ratios, we need to simply multiply by the density
         ! of *dry* air
         allocate (rhodry(nx, ny, kprs))
         do k = 1, kprs
            rhodry(:, :, k) = prslvl(k)/(r*tprs(:, :, k))
         end do
         allocate (cldliqcon_prs(nx, ny, kprs))
         allocate (cldicecon_prs(nx, ny, kprs))
         allocate (raincon_prs(nx, ny, kprs))
         allocate (snowcon_prs(nx, ny, kprs))
         allocate (graupelcon_prs(nx, ny, kprs))
         cldliqcon_prs = cldliqmr_prs*rhodry
         cldicecon_prs = cldicemr_prs*rhodry
         raincon_prs = rainmr_prs*rhodry
         snowcon_prs = snowmr_prs*rhodry
         graupelcon_prs = graupelmr_prs*rhodry
         deallocate (rhodry)

         call output_laps_format(zprs, uprs, vprs, wprs, omprs, tprs, shprs, rhprs, &
                                 cldliqcon_prs, cldicecon_prs, raincon_prs, &
                                 snowcon_prs, graupelcon_prs, tkeprs, &
                                 refl_prs, pcptype_prs, usfc, vsfc, wsfc, &
                                 tsfc, tdsfc, rhsfc, cldbase, cldtop, pmsl, redp, psfc, &
                                 intliqwater, totpcpwater, pcp_inc, pcp_tot, &
                                 snow_inc, snow_tot, thetasfc, thetaesfc, cape, cin, &
                                 cldamt, ceiling, echo_tops, max_refl, refl_sfc, &
                                 pcptype_sfc, srhel, liftedind, heatind, visibility, &
                                 terdot, lwout, swout, &
                                 shflux, lhflux, pblhgt, ground_t, &
                                 upbl, vpbl, vnt_index, ham_index, hah_index, fwi_index, &
                                 prslvl*0.01, lfmprd_dir, laps_data_root, domain_num, &
                                 laps_reftime, laps_valtime, nx, ny, kprs, realtime, &
                                 write_to_lapsdir(domain_num), &
                                 model_name, make_donefile)

         deallocate (cldliqcon_prs)
         deallocate (cldicecon_prs)
         deallocate (raincon_prs)
         deallocate (snowcon_prs)
         deallocate (graupelcon_prs)
      end if

      ! make grib output as required

      if ((gribsfc(domain_num)) .or. (gribua(domain_num))) then
         if (t .eq. 1) then
            call make_igds(proj, igds)
         end if
         if (laps_valtime_prev .gt. 0) then
            period_sec = laps_valtime - laps_valtime_prev
         else
            period_sec = 0
         end if
         call make_fnam_lp(laps_reftime, a9time, istatus)
         call make_fcst_time(laps_valtime, laps_reftime, fcst_hhmm, istatus)
         if (fcst_hhmm(3:4) .eq. "00") then
            gribfile = gribdir(1:gdir_len)//'/'//a9time//'00'//fcst_hhmm(1:2) &
                       //'.grib'
         else
            gribfile = gribdir(1:gdir_len)//'/'//a9time//'00'//fcst_hhmm//'.grib'
         end if

         call s_len(gribfile, gfile_len)
         print *, 'opening ', gribfile(1:gfile_len)
         call open_grib_c(gribfile, funit)
         print *, 'opened funit ', funit, ' for ', trim(gribfile)
         startb = 1
         nbytes = 0
         if (gribsfc(domain_num)) then
            call grib_sfc_vars(table_version, center_id, subcenter_id, &
                               process_id(domain_num), laps_reftime, laps_valtime, &
                               period_sec, igds, nx, ny, &
                               tsfc, tdsfc, rhsfc, usfc, vsfc, &
                               wsfc, pmsl, psfc, totpcpwater, pcp_inc, pcp_tot, snow_inc, &
                               snow_tot, thetasfc, thetaesfc, cape, cin, srhel, &
                               liftedind, terdot, lwout, swout, lwdown, swdown, &
                               shflux, lhflux, pblhgt, ground_t, &
                               clwmrsfc, icemrsfc, rainmrsfc, snowmrsfc, graupmrsfc, &
                               cldamt, cldbase, cldtop, visibility, &
                               ceiling, echo_tops, max_refl, refl_sfc, &
                               pcptype_sfc, funit, startb, nbytes)
            startb = startb + nbytes
         end if
         if (gribua(domain_num)) then

            call grib_ua_vars(table_version, center_id, subcenter_id, &
                              process_id(domain_num), laps_reftime, laps_valtime, &
                              period_sec, igds, nx, ny, kprs, &
                              prslvl*0.01, zprs, uprs, vprs, wprs, omprs, tprs, shprs, rhprs, &
                              cldliqmr_prs, cldicemr_prs, rainmr_prs, snowmr_prs, graupelmr_prs, &
                              pcptype_prs, refl_prs, tkeprs, funit, startb, nbytes)

            startb = startb + nbytes
         end if
         print *, 'number of bytes written in grib: ', nbytes

         call close_grib_c(funit)
         gribdone = gribfile(1:gfile_len)//'.done'
         call get_file_unit(flagunit)
         open (unit=flagunit, file=gribdone, status='unknown')
         close (flagunit)
      end if
      ! make vis5d output if requested in namelist

      if (make_v5d(domain_num)) then

         print *, 'making vis5d output...'
         ! if this is the first output time, then the vis5d file
         ! needs to be initialized
         if (t .eq. 1) then

            print '(a)', 'initializing vis5d data file...'
            call v5dinit(v5d_compress)

         end if
         where (refl_prs .lt. 0) refl_prs = 0.
         where (refl_sfc .lt. 0) refl_sfc = 0.
         where (max_refl .lt. 0) max_refl = 0.
         ! call v5dout to fill data for this time period
         if (time_index_inc .ne. 1) then
            t_v5d = t/time_index_inc + 1
         else
            t_v5d = t
         end if
         print '(a,i3)', 'populating vis5d file for time period ', t_v5d
         call v5dout(t_v5d, zprs, tprs, tdprs, rhprs, uprs, vprs, wprs, omprs, &
                     abs_vort, shprs, cldliqmr_prs, cldicemr_prs, rainmr_prs, &
                     snowmr_prs, graupelmr_prs, refl_prs, pcptype_prs, tkeprs, &
                     thick_10_5, tsfc, tdsfc, rhsfc, usfc, vsfc, wsfc, &
                     pmsl, psfc, cldbase, cldtop, cldamt, ceiling, &
                     intliqwater, totpcpwater, pcp_inc, pcp_tot, &
                     con_pcp_inc, con_pcp_tot, &
                     snow_inc, snow_tot, pcptype_sfc, thetasfc, thetaesfc, &
                     cape, cin, liftedind, srhel, refl_sfc, max_refl, echo_tops, &
                     heatind, visibility, snowcover, lwout, swout, lwdown, swdown, &
                     shflux, lhflux, &
                     pblhgt, upbl, vpbl, ground_t)
         if (t .eq. num_times_to_proc) call v5dend

      end if

      if (make_points(domain_num)) then
         call split_date_char(point_times(t), year, month, day, hour, minute, second)
         write (date_str, '(i2.2,"/",i2.2,"/",i4.4,1x,i2.2,":",i2.2)') month, &
            day, year, hour, minute
         point_loop: do pt = 1, num_points
            ip = nint(point_rec(pt)%i)
            jp = nint(point_rec(pt)%j)
            t_pt = bint(point_rec(pt)%i, point_rec(pt)%j, tsfc, nx, ny)
            td_pt = bint(point_rec(pt)%i, point_rec(pt)%j, tdsfc, nx, ny)
            ! convert t/td units
            if (point_temp_units .eq. 'f') then
               t_pt = (t_pt - 273.15)*9./5.+32.
               td_pt = (td_pt - 273.15)*9./5.+32.
            elseif (point_temp_units .eq. 'c') then
               t_pt = t_pt - 273.15
               td_pt = td_pt - 273.15
            end if
            rh_pt = nint(bint(point_rec(pt)%i, point_rec(pt)%j, rhsfc, nx, ny))

            ! do winds...now includes the mean pbl wind (variables that
            ! have an extra p in their name are pbl values
            u_pt_g = bint(point_rec(pt)%i, point_rec(pt)%j, usfc, nx, ny)
            v_pt_g = bint(point_rec(pt)%i, point_rec(pt)%j, vsfc, nx, ny)
            up_pt_g = bint(point_rec(pt)%i, point_rec(pt)%j, upbl, nx, ny)
            vp_pt_g = bint(point_rec(pt)%i, point_rec(pt)%j, vpbl, nx, ny)
            call gridwind_to_truewind(londot(ip, jp), proj, u_pt_g, v_pt_g, &
                                      u_pt, v_pt)
            call gridwind_to_truewind(londot(ip, jp), proj, up_pt_g, vp_pt_g, &
                                      up_pt, vp_pt)
            call uv_to_disp(u_pt, v_pt, dir_pt_r, spd_pt_r)
            call uv_to_disp(up_pt, vp_pt, dirp_pt_r, spdp_pt_r)
            dir_pt = nint(dir_pt_r/10.)*10  ! integer to nearest 10 degrees
            dirp_pt = nint(dirp_pt_r/10.)*10.

            ! convert wind speed
            if (point_windspd_units .eq. 'kts') then
               spd_pt = nint(spd_pt_r*1.9425) ! convert to knots
               spdp_pt = nint(spdp_pt_r*1.9425)
            elseif (point_windspd_units .eq. 'mph') then
               spd_pt = nint(spd_pt_r*2.2369)
               spdp_pt = nint(spdp_pt_r*2.2369)
            elseif (point_windspd_units .eq. 'm/s') then
               spd_pt = nint(spd_pt_r)
               spdp_pt = nint(spdp_pt_r)
            end if
            vis_pt = bint(point_rec(pt)%i, point_rec(pt)%j, visibility, nx, ny) &
                     *0.00062317
            cld_pt = bint(point_rec(pt)%i, point_rec(pt)%j, cldamt, nx, ny)
            ceiling_pt_r = bint(point_rec(pt)%i, point_rec(pt)%j, ceiling, nx, ny)
            if (ceiling_pt_r .lt. 50000.) then
               ceiling_pt = nint(ceiling_pt_r*3.2808/100.)  ! hundreds of feet
            else
               ceiling_pt = 999
            end if
            select case (int(pcptype_sfc(ip, jp)))
            case (0)
               if (vis_pt .ge. 7.) then
                  if (cld_pt .ge. 0.75) then
                     wx_pt = 'cloudy  '
                  else if ((cld_pt .lt. 0.75) .and. (cld_pt .ge. 0.25)) then
                     wx_pt = 'pt cldy '
                  else
                     wx_pt = 'clear   '
                  end if
               else if ((vis_pt .lt. 7.) .and. (vis_pt .ge. 3.)) then
                  wx_pt = 'haze    '
               else
                  wx_pt = 'fog     '
               end if

            case (1)
               wx_pt = 'rain    '
            case (2)
               wx_pt = 'snow    '
            case (3)
               wx_pt = 'frzrain '
            case (4)
               wx_pt = 'sleet   '
            case (5)
               wx_pt = 'hail    '
            case (6)
               wx_pt = 'drizzle '
            case (7)
               wx_pt = 'mixed   '
            case (8)
               wx_pt = 'mixed   '
            case (9)
               wx_pt = 'rain/ice'
            case default
               wx_pt = 'unknown '
            end select
            pcp_pt = bint(point_rec(pt)%i, point_rec(pt)%j, pcp_inc, nx, ny)*39.37
            snow_pt = bint(point_rec(pt)%i, point_rec(pt)%j, snow_inc, nx, ny)*39.37

            vnt_pt = bint(point_rec(pt)%i, point_rec(pt)%j, vnt_index, nx, ny)
            pbl_pt = bint(point_rec(pt)%i, point_rec(pt)%j, pblhgt, nx, ny)*3.28
            if (point_vent_units .eq. 'kt-ft') vnt_pt = vnt_pt*6.3774
            write (point_rec(pt)%output_unit, &
                   '(a,3i4,1x,i3.3,"/",i2.2,1x,i3.3,1x,f4.1,1x,a,1x,f5.2,1x,f4.1,1x,i6,1x,i5,1x,i3.3,"/",i2.2,1x,i2,1x,i2,1x,i3)') &
               date_str, nint(t_pt), nint(td_pt), rh_pt, dir_pt, spd_pt, ceiling_pt, vis_pt, wx_pt, &
               pcp_pt, snow_pt, nint(vnt_pt), nint(pbl_pt), dirp_pt, spdp_pt, nint(ham_index(ip, jp)), &
               nint(hah_index(ip, jp)), nint(fwi_index(ip, jp))

            if (t_pt .gt. point_rec(pt)%hi_temp) then
               point_rec(pt)%hi_temp = t_pt
               point_rec(pt)%hi_temp_time = date_str
            end if
            if (t_pt .lt. point_rec(pt)%lo_temp) then
               point_rec(pt)%lo_temp = t_pt
               point_rec(pt)%lo_temp_time = date_str
            end if
            point_rec(pt)%total_pcp = point_rec(pt)%total_pcp + pcp_pt
            point_rec(pt)%total_snow = point_rec(pt)%total_snow + snow_pt
            point_rec(pt)%avg_temp = point_rec(pt)%avg_temp + t_pt
            point_rec(pt)%avg_dewpt = point_rec(pt)%avg_dewpt + td_pt
            if (t .eq. num_times_to_proc) then
               ! finish summary and close the file
               write (point_rec(pt)%output_unit, '(80x)')
               write (point_rec(pt)%output_unit, '(80x)')
               write (point_rec(pt)%output_unit, '("summary information for period")')
               write (point_rec(pt)%output_unit, &
                      '("--------------------------------------------------")')
               write (point_rec(pt)%output_unit, &
                      '("high temperature:  ",f6.1,2x,"at ",a)') point_rec(pt)%hi_temp, &
                  point_rec(pt)%hi_temp_time
               write (point_rec(pt)%output_unit, &
                      '("low temperature:   ",f6.1,2x,"at ",a)') point_rec(pt)%lo_temp, &
                  point_rec(pt)%lo_temp_time
               point_rec(pt)%avg_temp = point_rec(pt)%avg_temp/ &
                                        float(num_times_to_proc/time_index_inc)
               write (point_rec(pt)%output_unit, &
                      '("avg temperature:   ",f6.1)') point_rec(pt)%avg_temp
               point_rec(pt)%avg_dewpt = point_rec(pt)%avg_dewpt/ &
                                         float(num_times_to_proc/time_index_inc)
               write (point_rec(pt)%output_unit, &
                      '("avg dewpoint:      ",f6.1)') point_rec(pt)%avg_dewpt
               write (point_rec(pt)%output_unit, &
                      '("total precip:      ",f6.2)') point_rec(pt)%total_pcp
               write (point_rec(pt)%output_unit, &
                      '("total snow:        ",f6.2)') point_rec(pt)%total_snow
               close (point_rec(pt)%output_unit)
            end if
         end do point_loop
      end if
      ! deallocate all variables except pcp/snow init/total
      deallocate (psfc)
      deallocate (tsfc)
      deallocate (thetasfc)
      deallocate (thetaesfc)
      deallocate (rhsfc)
      deallocate (tdsfc)
      deallocate (thetaprs)
      deallocate (zprs)
      deallocate (rhprs)
      deallocate (tprs)
      deallocate (tdprs)
      deallocate (shprs)
      deallocate (redp)
      deallocate (pmsl)
      deallocate (usfc)
      deallocate (uprs)
      deallocate (vsfc)
      deallocate (upbl)
      deallocate (vpbl)
      deallocate (vprs)
      deallocate (wsfc)
      deallocate (wprs)
      deallocate (tkeprs)
      deallocate (omprs)
      deallocate (cldbase)
      deallocate (cldtop)
      deallocate (cldamt)
      deallocate (ceiling)
      deallocate (intliqwater)
      deallocate (totpcpwater)
      deallocate (max_refl)
      deallocate (echo_tops)
      deallocate (cldliqmr_prs)
      deallocate (cldicemr_prs)
      deallocate (rainmr_prs)
      deallocate (snowmr_prs)
      deallocate (graupelmr_prs)
      deallocate (refl_prs)
      deallocate (refl_sfc)
      deallocate (pcptype_sfc)
      deallocate (pcptype_prs)
      deallocate (pcp_inc)
      deallocate (con_pcp_inc)
      deallocate (snow_inc)
      deallocate (cape)
      deallocate (cin)
      deallocate (liftedind)
      deallocate (visibility)
      deallocate (heatind)
      deallocate (lwout)
      deallocate (swout)
      deallocate (lwdown)
      deallocate (swdown)
      deallocate (albedo)
      deallocate (shflux)
      deallocate (lhflux)
      deallocate (pblhgt)
      deallocate (ground_t)
      deallocate (clwmrsfc)
      deallocate (icemrsfc)
      deallocate (rainmrsfc)
      deallocate (snowmrsfc)
      deallocate (graupmrsfc)
      deallocate (vnt_index)
      deallocate (ham_index)
      deallocate (hah_index)
      deallocate (fwi_index)

      if (make_v5d(domain_num)) then
         deallocate (abs_vort)
         deallocate (thick_10_5)
         deallocate (snowcover)
      end if
      if (realtime) then
         if (split_output) then
            ! close(89,status='delete')
         end if
      end if
      laps_valtime_prev = laps_valtime
   end do time_loop

   ! deallocate the precip arrays
   deallocate (pcp_init)
   deallocate (con_pcp_tot)
   deallocate (con_pcp_init)
   deallocate (pcp_tot)
   deallocate (snow_init)
   deallocate (snow_tot)
   deallocate (latdot)
   deallocate (londot)
   deallocate (terdot)
end program lfmpost

