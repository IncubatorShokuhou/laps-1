  subroutine grib_sfc_vars(table_version, center_id, subcenter_id, &
                           process_id, laps_reftime, laps_valtime, period_sec, igds, &
                           nx, ny, tsfc, tdsfc, rhsfc, usfc, vsfc, &
                           wsfc, pmsl, psfc, totpcpwater, pcp_inc, pcp_tot, snow_inc, snow_tot, &
                           thetasfc, thetaesfc, cape, cin, srhel, liftedind, terdot, lwout, &
                           swout, lwdown, swdown, shflux, lhflux, pblhgt, ground_t, &
                           clwmrsfc, icemrsfc, rainmrsfc, snowmrsfc, graupmrsfc, &
                           cldamt, cldbase, cldtop, &
                           visibility, ceiling, echo_tops, max_refl, refl_sfc, &
                           pcptype_sfc, funit, startb, nbytes)

     use grib
     implicit none

     integer, intent(in)          :: table_version
     integer, intent(in)          :: center_id
     integer, intent(in)          :: process_id
     integer, intent(in)          :: subcenter_id
     integer, intent(in)          :: laps_reftime
     integer, intent(in)          :: laps_valtime
     integer, intent(in)          :: period_sec
     integer, intent(in)          :: igds(18)
     integer, intent(in)          :: nx
     integer, intent(in)          :: ny
     real, intent(in)             :: tsfc(nx, ny)
     real, intent(in)             :: tdsfc(nx, ny)
     real, intent(in)             :: rhsfc(nx, ny)
     real, intent(in)             :: usfc(nx, ny)
     real, intent(in)             :: vsfc(nx, ny)
     real, intent(in)             :: wsfc(nx, ny)
     real, intent(in)             :: pmsl(nx, ny)
     real, intent(in)             :: psfc(nx, ny)
     real, intent(in)             :: totpcpwater(nx, ny)
     real, intent(in)             :: pcp_inc(nx, ny)
     real, intent(in)             :: pcp_tot(nx, ny)
     real, intent(in)             :: snow_inc(nx, ny)
     real, intent(in)             :: snow_tot(nx, ny)
     real, intent(in)             :: thetasfc(nx, ny)
     real, intent(in)             :: thetaesfc(nx, ny)
     real, intent(in)             :: cape(nx, ny)
     real, intent(in)             :: cin(nx, ny)
     real, intent(in)             :: srhel(nx, ny)
     real, intent(in)             :: liftedind(nx, ny)
     real, intent(in)             :: terdot(nx, ny)
     real, intent(in)             :: lwout(nx, ny)
     real, intent(in)             :: swout(nx, ny)
     real, intent(in)             :: lwdown(nx, ny)
     real, intent(in)             :: swdown(nx, ny)
     real, intent(in)             :: shflux(nx, ny)
     real, intent(in)             :: lhflux(nx, ny)
     real, intent(in)             :: pblhgt(nx, ny)
     real, intent(in)             :: ground_t(nx, ny)
     real, intent(in)             :: clwmrsfc(nx, ny)
     real, intent(in)             :: icemrsfc(nx, ny)
     real, intent(in)             :: rainmrsfc(nx, ny)
     real, intent(in)             :: snowmrsfc(nx, ny)
     real, intent(in)             :: graupmrsfc(nx, ny)
     real, intent(in)             :: cldamt(nx, ny)
     real, intent(in)             :: cldbase(nx, ny)
     real, intent(in)             :: cldtop(nx, ny)
     real, intent(in)             :: visibility(nx, ny)
     real, intent(in)             :: ceiling(nx, ny)
     real, intent(in)             :: echo_tops(nx, ny)
     real, intent(in)             :: max_refl(nx, ny)
     real, intent(in)             :: refl_sfc(nx, ny)
     real, intent(in)             :: pcptype_sfc(nx, ny)

     integer, intent(in)          :: funit
     integer, intent(in)          :: startb
     integer, intent(out)         :: nbytes

     integer                     :: i, j
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
     integer                     :: startbyte
     real, parameter             :: rmissing = 1.e36
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

     ! grib up each variable...
     nbytes = startb - 1
     startbyte = nbytes + startb

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! surface temperature
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = tsfc(i, j + 1)
        end do
     end do
     print *, 'gribbing tsfc min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 11
     leveltype = 105
     level1 = 2
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! surface dewpoint
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = tdsfc(i, j + 1)
        end do
     end do
     print *, 'gribbing tdsfc min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 17
     leveltype = 105
     level1 = 2
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

    !!!!!!!!!!!!!!!!!!!!!!!
     ! surface rh
    !!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = rhsfc(i, j + 1)
        end do
     end do
     print *, 'gribbing sfc rh min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 52
     leveltype = 105
     level1 = 2
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! surface u wind
    !!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = usfc(i, j + 1)
        end do
     end do
     print *, 'gribbing sfc u min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 33
     leveltype = 105
     level1 = 10
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

    !!!!!!!!!!!!!!!!!!!!!!!!
     ! surface v wind
    !!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = vsfc(i, j + 1)
        end do
     end do
     print *, 'gribbing sfc v min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 34
     leveltype = 105
     level1 = 10
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

    !!!!!!!!!!!!!!!!!!!!!!!!
     ! surface w wind
    !!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = wsfc(i, j + 1)
        end do
     end do
     print *, 'gribbing sfc w min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 40
     leveltype = 105
     level1 = 10
     level2 = 0
     timerange = 0
     timeperiod1 = fcsttime_now
     timeperiod2 = 0
     scalep10 = 4
     call make_id(table_version, center_id, subcenter_id, process_id, &
                  param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                  hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                  scalep10, id)
     call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
     nbytes = nbytes + itot
     startbyte = nbytes + 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! sea level pressure
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = pmsl(i, j + 1)
        end do
     end do
     print *, 'gribbing mslp min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 2
     leveltype = 102
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! surface pressure
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = psfc(i, j + 1)
        end do
     end do
     print *, 'gribbing sfc p min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 1
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! total preciptable water
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = totpcpwater(i, j + 1)*1000.
        end do
     end do
     print *, 'gribbing total pw min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 54
     leveltype = 200
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! potential temperature
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = thetasfc(i, j + 1)
        end do
     end do
     print *, 'gribbing theta min/max =', minval(fld), maxval(fld)
     itype = 0
     param = 13
     leveltype = 105
     level1 = 2
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! equiv potential temperature
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = thetaesfc(i, j + 1)
        end do
     end do
     print *, 'gribbing thetae min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 14
     leveltype = 105
     level1 = 2
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! cape
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = cape(i, j + 1)
        end do
     end do
     print *, 'gribbing cape min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 157
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! cin
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = cin(i, j + 1)
        end do
     end do
     print *, 'gribbing cin min/max =', minval(fld), maxval(fld)
     itype = 0
     param = 156
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! storm relative helicity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = srhel(i, j + 1)
        end do
     end do
     print *, 'gribbing srh min/max =', minval(fld), maxval(fld)
     itype = 0
     param = 190
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! lifted index
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = liftedind(i, j + 1)
        end do
     end do
     print *, 'gribbing li min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 131
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! ground t
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = ground_t(i, j + 1)
        end do
     end do
     print *, 'gribbing tgd min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 11
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! cloud liquid
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = clwmrsfc(i, j + 1)
        end do
     end do
     print *, 'gribbing clw min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 153
     leveltype = 1
     level1 = 0
     level2 = 0
     timerange = 0
     timeperiod1 = fcsttime_now
     timeperiod2 = 0
     scalep10 = 4
     call make_id(table_version, center_id, subcenter_id, process_id, &
                  param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                  hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                  scalep10, id)
     call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
     nbytes = nbytes + itot
     startbyte = nbytes + 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! cloud ice
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = icemrsfc(i, j + 1)
        end do
     end do
     print *, 'gribbing ice min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 178
     leveltype = 1
     level1 = 0
     level2 = 0
     timerange = 0
     timeperiod1 = fcsttime_now
     timeperiod2 = 0
     scalep10 = 4
     call make_id(table_version, center_id, subcenter_id, process_id, &
                  param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                  hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                  scalep10, id)
     call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
     nbytes = nbytes + itot
     startbyte = nbytes + 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! rain mr
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = rainmrsfc(i, j + 1)
        end do
     end do
     print *, 'gribbing rain min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 170
     leveltype = 1
     level1 = 0
     level2 = 0
     timerange = 0
     timeperiod1 = fcsttime_now
     timeperiod2 = 0
     scalep10 = 4
     call make_id(table_version, center_id, subcenter_id, process_id, &
                  param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                  hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                  scalep10, id)
     call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
     nbytes = nbytes + itot
     startbyte = nbytes + 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!
     ! snow mr
    !!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = snowmrsfc(i, j + 1)
        end do
     end do
     print *, 'gribbing snow min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 171
     leveltype = 1
     level1 = 0
     level2 = 0
     timerange = 0
     timeperiod1 = fcsttime_now
     timeperiod2 = 0
     scalep10 = 4
     call make_id(table_version, center_id, subcenter_id, process_id, &
                  param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                  hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                  scalep10, id)
     call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
     nbytes = nbytes + itot
     startbyte = nbytes + 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! graupel mr
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = graupmrsfc(i, j + 1)
        end do
     end do
     print *, 'gribbing graupel min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 179
     leveltype = 1
     level1 = 0
     level2 = 0
     timerange = 0
     timeperiod1 = fcsttime_now
     timeperiod2 = 0
     scalep10 = 4
     call make_id(table_version, center_id, subcenter_id, process_id, &
                  param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                  hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                  scalep10, id)
     call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
     nbytes = nbytes + itot
     startbyte = nbytes + 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! total cloud cover
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = cldamt(i, j + 1)*100.
        end do
     end do
     print *, 'gribbing lcv = ', minval(fld), maxval(fld)
     itype = 0
     param = 71
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!
     ! cloud base
    !!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = cldbase(i, j + 1)
        end do
     end do
     print *, 'gribbing cloud base min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 138
     leveltype = 102
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! cloud top
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = cldtop(i, j + 1)
        end do
     end do
     print *, 'gribbing cloud top min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 139
     leveltype = 102
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! cloud ceiling
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = ceiling(i, j + 1)
        end do
     end do
     print *, 'gribbing cloud ceiling min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 137
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! surface radar reflectivity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = refl_sfc(i, j + 1)
        end do
     end do
     print *, 'gribbing sfc refl min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 128
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! column max radar reflectivity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = max_refl(i, j + 1)
        end do
     end do
     print *, 'gribbing col max refl min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 129
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! max echo tops
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = echo_tops(i, j + 1)
        end do
     end do
     print *, 'gribbing max echo top min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 130
     leveltype = 102
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! surface precip type
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = pcptype_sfc(i, j + 1)
        end do
     end do
     print *, 'gribbing ptype sfc min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 136
     leveltype = 1
     level1 = 0
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! visbility
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = visibility(i, j + 1)
        end do
     end do
     print *, 'gribbing vis min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 20
     leveltype = 1
     level1 = 0
     level2 = 0
     timerange = 0
     timeperiod1 = fcsttime_now
     timeperiod2 = 0
     scalep10 = -3
     call make_id(table_version, center_id, subcenter_id, process_id, &
                  param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                  hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                  scalep10, id)
     call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
     nbytes = nbytes + itot
     startbyte = nbytes + 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! some fields only present in true forecast fields
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (fcsttime_now .gt. 0) then

        ! incremental precip
        do j = 0, ny - 1
           do i = 1, nx
              fld(j*nx + i) = pcp_inc(i, j + 1)*1000.
           end do
        end do
        print *, 'gribbing pcpinc min/max =', minval(fld), maxval(fld)
        itype = 0
        param = 61
        leveltype = 1
        level1 = 0
        level2 = 0
        timerange = 4
        timeperiod1 = fcsttime_prev
        timeperiod2 = fcsttime_now
        scalep10 = 1
        call make_id(table_version, center_id, subcenter_id, process_id, &
                     param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                     hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                     scalep10, id)
        call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
        nbytes = nbytes + itot
        startbyte = nbytes + 1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! total accum precip
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j = 0, ny - 1
           do i = 1, nx
              fld(j*nx + i) = pcp_tot(i, j + 1)*1000.
           end do
        end do
        print *, 'gribbing totpcp min/max = ', minval(fld), maxval(fld)
        itype = 0
        param = 142
        leveltype = 1
        level1 = 0
        level2 = 0
        timerange = 4
        timeperiod1 = 0
        timeperiod2 = fcsttime_now
        scalep10 = 1
        call make_id(table_version, center_id, subcenter_id, process_id, &
                     param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                     hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                     scalep10, id)
        call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
        nbytes = nbytes + itot
        startbyte = nbytes + 1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! incremental snow
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j = 0, ny - 1
           do i = 1, nx
              fld(j*nx + i) = snow_inc(i, j + 1)
           end do
        end do
        print *, 'gribbing snowinc min/max = ', minval(fld), maxval(fld)
        itype = 0
        param = 66
        leveltype = 1
        level1 = 0
        level2 = 0
        timerange = 4
        timeperiod1 = fcsttime_prev
        timeperiod2 = fcsttime_now
        scalep10 = 4
        call make_id(table_version, center_id, subcenter_id, process_id, &
                     param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                     hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                     scalep10, id)
        call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
        nbytes = nbytes + itot
        startbyte = nbytes + 1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! total accum  snow
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j = 0, ny - 1
           do i = 1, nx
              fld(j*nx + i) = snow_tot(i, j + 1)
           end do
        end do
        print *, 'gribbing snow tot min/max =', minval(fld), maxval(fld)
        itype = 0
        param = 141
        leveltype = 1
        level1 = 0
        level2 = 0
        timerange = 4
        timeperiod1 = 0
        timeperiod2 = fcsttime_now
        scalep10 = 4
        call make_id(table_version, center_id, subcenter_id, process_id, &
                     param, leveltype, level1, level2, yyyyr, mmr, ddr, &
                     hhr, minr, timeunit, timerange, timeperiod1, timeperiod2, &
                     scalep10, id)
        call write_grib(itype, fld, id, igds, funit, startbyte, itot, istatus)
        nbytes = nbytes + itot
        startbyte = nbytes + 1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! outgoing longwave radiation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (maxval(lwout) .lt. rmissing) then
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = lwout(i, j + 1)
              end do
           end do
           print *, 'gribbing lwout min/max = ', minval(fld), maxval(fld)
           itype = 0
           param = 212
           leveltype = 8
           level1 = 0
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
        end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! outgoing shortwave radiation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (maxval(swout) .lt. rmissing) then
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = swout(i, j + 1)
              end do
           end do
           print *, 'gribbing swout min/max = ', minval(fld), maxval(fld)
           itype = 0
           param = 211
           leveltype = 8
           level1 = 0
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
        end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! incoming shortwave radiation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (maxval(swdown) .lt. rmissing) then
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = swdown(i, j + 1)
              end do
           end do
           print *, 'gribbing swdown min/max = ', minval(fld), maxval(fld)
           itype = 0
           param = 111
           leveltype = 1
           level1 = 0
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
        end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! incoming longwave radiation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (maxval(lwdown) .lt. rmissing) then
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = lwdown(i, j + 1)
              end do
           end do
           print *, 'gribbing lwdown min/max = ', minval(fld), maxval(fld)
           itype = 0
           param = 112
           leveltype = 1
           level1 = 0
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
        end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  sensible heat flux
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (maxval(shflux) .lt. rmissing) then
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = shflux(i, j + 1)
              end do
           end do
           print *, 'gribbing shf min/max =', minval(fld), maxval(fld)
           itype = 0
           param = 122
           leveltype = 1
           level1 = 0
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
        end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  latent heat flux
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (maxval(lhflux) .lt. rmissing) then
           do j = 0, ny - 1
              do i = 1, nx
                 fld(j*nx + i) = lhflux(i, j + 1)
              end do
           end do
           print *, 'gribbing lhf min/max = ', minval(fld), maxval(fld)
           itype = 0
           param = 121
           leveltype = 1
           level1 = 0
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
        end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  pbl height
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j = 0, ny - 1
           do i = 1, nx
              fld(j*nx + i) = pblhgt(i, j + 1)
           end do
        end do
        print *, 'gribbing pblhgt min/max =', minval(fld), maxval(fld)
        itype = 0
        param = 221
        leveltype = 1
        level1 = 0
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
     end if

     ! topography
     do j = 0, ny - 1
        do i = 1, nx
           fld(j*nx + i) = terdot(i, j + 1)
        end do
     end do
     print *, 'gribbing topo min/max = ', minval(fld), maxval(fld)
     itype = 0
     param = 7
     leveltype = 1
     level1 = 0
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
     return
  end subroutine grib_sfc_vars

