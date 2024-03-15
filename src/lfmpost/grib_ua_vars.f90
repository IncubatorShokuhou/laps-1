  subroutine grib_ua_vars(table_version, center_id, subcenter_id, &
                          process_id, laps_reftime, laps_valtime, period_sec, igds, nx, ny, kprs, &
                          plvlmb, zprs, uprs, vprs, wprs, omprs, tprs, shprs, rhprs, &
                          cldliqmr_prs, cldicemr_prs, rainmr_prs, snowmr_prs, graupelmr_prs, &
                          pcptype_prs, refl_prs, tkeprs, funit, startb, nbytes)

     use grib
     implicit none
     integer, intent(in)          :: table_version
     integer, intent(in)          :: center_id
     integer, intent(in)          :: subcenter_id
     integer, intent(in)          :: process_id
     integer, intent(in)          :: laps_reftime
     integer, intent(in)          :: laps_valtime
     integer, intent(in)          :: period_sec
     integer, intent(in)          :: igds(18)
     integer, intent(in)          :: nx
     integer, intent(in)          :: ny
     integer, intent(in)          :: kprs
     real, intent(in)             :: plvlmb(kprs)
     real, intent(in)             :: zprs(nx, ny, kprs)
     real, intent(in)             :: uprs(nx, ny, kprs)
     real, intent(in)             :: vprs(nx, ny, kprs)
     real, intent(in)             :: wprs(nx, ny, kprs)
     real, intent(in)             :: omprs(nx, ny, kprs)
     real, intent(in)             :: tprs(nx, ny, kprs)
     real, intent(in)             :: shprs(nx, ny, kprs)
     real, intent(in)             :: rhprs(nx, ny, kprs)
     real, intent(in)             :: cldliqmr_prs(nx, ny, kprs)
     real, intent(in)             :: cldicemr_prs(nx, ny, kprs)
     real, intent(in)             :: rainmr_prs(nx, ny, kprs)
     real, intent(in)             :: snowmr_prs(nx, ny, kprs)
     real, intent(in)             :: graupelmr_prs(nx, ny, kprs)
     real, intent(in)             :: pcptype_prs(nx, ny, kprs)
     real, intent(in)             :: refl_prs(nx, ny, kprs)
     real, intent(in)             :: tkeprs(nx, ny, kprs)
     integer, intent(in)          :: funit
     integer, intent(in)          :: startb
     integer, intent(out)         :: nbytes

     integer                     :: i, j, k
     integer                     :: lvl
     integer                     :: itype
     integer                     :: istatus
     integer                     :: id(27)
     integer                     :: param
     integer                     :: leveltype
     integer                     :: level1
     integer                     :: level2
     integer                     :: yyyyr
     integer                     :: mmr
     integer                     :: ddr
     integer                     :: hhr
     integer                     :: minr
     integer                     :: timeunit
     integer                     :: timerange
     integer                     :: timeperiod1
     integer                     :: timeperiod2
     integer                     :: scalep10
     character(len=24)           :: atime
     real                        :: fld(nx*ny)
     character(len=3)            :: amonth
     character(len=3)            :: amonths(12)
     integer                     :: fcsttime_now
     integer                     :: fcsttime_prev
     integer                     :: itot
     integer                    :: startbyte
     data amonths/'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
        'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/
     ! compute year, month, day of month, hour, and minute from laps_reftime

     call cv_i4tim_asc_lp(laps_reftime, atime, istatus)
     read (atime, '(i2.2,x,a3,x,i4.4,x,i2.2,x,i2.2)') ddr, amonth, yyyyr, &
        hhr, minr
     do i = 1, 12
        if (amonth .eq. amonths(i)) then
           mmr = i
           exit
        end if
     end do

     ! determine appropriate timeunit

     if (mod(period_sec, 3600) .eq. 0) then
        ! time unit shoud be hours
        timeunit = 1
        fcsttime_now = (laps_valtime - laps_reftime)/3600
        if (fcsttime_now .gt. 0) then
           fcsttime_prev = fcsttime_now - (period_sec/3600)
        else
           fcsttime_prev = 0
        end if
     else
        ! time unit in minutes
        timeunit = 0
        fcsttime_now = (laps_valtime - laps_reftime)/60
        if (fcsttime_now .gt. 0) then
           fcsttime_prev = fcsttime_now - (period_sec/60)
        else
           fcsttime_prev = 0
        end if
     end if

     ! grib up each variable at each level...

     nbytes = startb - 1
     startbyte = startb + nbytes

     levelloop: do k = 1, kprs
        lvl = nint(plvlmb(k))

        if (lvl .le. 1000) then
           ! geopotential height
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = zprs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 7
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 0
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! temperature
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = tprs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 11
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 2
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! specific humidity
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = shprs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 51
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 8
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! relative humidity
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = rhprs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 52
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 1
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! u wind component
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = uprs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 33
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 1
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! v wind component
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = vprs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 34
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 1
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! w wind component
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = wprs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 40
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 3
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! omega
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = omprs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 39
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 3
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! cloud liquid
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = cldliqmr_prs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 153
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 6
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! cloud ice
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = cldicemr_prs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 178
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 6
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! rain
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = rainmr_prs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 170
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 6
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1
           ! snow
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = snowmr_prs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 171
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 6
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1

           ! graupel
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = graupelmr_prs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 179
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 6
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1

           ! precip type code
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = pcptype_prs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 136
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 0
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1

           ! radar reflectivity
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = refl_prs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 128
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 0
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1

           ! turbulent kinetic energy
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = tkeprs(i, j + 1, k)
              end do
           end do
           itype = 0
           param = 158
           leveltype = 100
           level1 = lvl
           level2 = 0
           timerange = 0
           timeperiod1 = fcsttime_now
           timeperiod2 = 0
           scalep10 = 3
           call make_id(table_version, center_id, subcenter_id, process_id, &
                        param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                        hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                        scalep10, id)
           call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
           nbytes = nbytes + itot
           startbyte = nbytes + 1

        end if ! only levels 1000mb and higher
     end do levelloop
     return
  end subroutine grib_ua_vars
