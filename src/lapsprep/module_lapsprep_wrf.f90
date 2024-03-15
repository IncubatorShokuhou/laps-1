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

module lapsprep_wrf

! purpose
! =======
! module to contain the various output routines needed for lapsprep
! to support initializition of the wrf model.
!
! subroutines contained
! =====================
! output_gribprep_format  - used to support wrf initializations
! output_gribprep_header  - writes the grib prep headers
! output_metgrid_format  - used to support wrf initializations using wps
! output_metgrid_header  - writes the metgrid headers
! remarks
! =======
!
!
! history
! =======
! 4 dec 2000 -- original -- brent shaw

   use setup
   use laps_static
   use date_pack
! use lapsprep_constants
   implicit none

   private
   integer, parameter :: gp_version = 4
   integer, parameter :: gp_version_wps = 5
   real, parameter    :: xfcst = 0.
   character(len=32), parameter :: source = &
                                   'laps analysis                   '
   character(len=8), parameter:: knownloc = 'swcorner'
   character(len=24) :: hdate
   integer            :: llflag
   character(len=9)  :: field
   character(len=25) :: units
   character(len=46) :: desc
   integer, parameter :: output_unit = 78
   real, parameter    :: slp_level = 201300.0

   public output_gribprep_format
   public output_metgrid_format
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine output_gribprep_format(p, t, ht, u, v, rh, slp, psfc, &
                                     lwc, rai, sno, ice, pic, snocov, tskin)

      !  subroutine of lapsprep that will build a file the
      !  wrfsi "gribprep" format that can be read by hinterp

      implicit none

      ! arguments

      real, intent(in)                   :: p(:)        ! pressure (hpa)
      real, intent(in)                   :: t(:, :, :)    ! temperature (k)
      real, intent(in)                   :: ht(:, :, :)   ! height (m)
      real, intent(in)                   :: u(:, :, :)    ! u-wind (m s{-1})
      real, intent(in)                   :: v(:, :, :)    ! v-wind (m s{-1})
      real, intent(in)                   :: rh(:, :, :)   ! relative humidity (%)
      real, intent(in)                   :: slp(:, :)    ! sea-level pressure (pa)
      real, intent(in)                   :: psfc(:, :)   ! surface pressure (pa)
      real, intent(in)                   :: lwc(:, :, :)  ! cloud liquid (kg/kg)
      real, intent(in)                   :: rai(:, :, :)  ! rain (kg/kg)
      real, intent(in)                   :: sno(:, :, :)  ! snow (kg/kg)
      real, intent(in)                   :: ice(:, :, :)  ! ice (kg/kg)
      real, intent(in)                   :: pic(:, :, :)  ! graupel (kg/kg)
      real, intent(in)                   :: snocov(:, :) ! snow cover (fract)
      real, intent(in)                   :: tskin(:, :)  ! skin temperature

      ! local variables

      integer            :: valid_mm, valid_dd
      character(len=256):: output_file_name
      real, allocatable  :: d2d(:, :)
      real, allocatable  :: p_pa(:)
      integer            :: k, yyyyddd
      integer            :: istatus
      real               :: r_missing_data

      ! allocate a scratch 2d array
      allocate (d2d(x, y))
      allocate (p_pa(z3 + 1))
      ! build the output file name

      output_prefix = trim(laps_data_root)//'/lapsprd/lapsprep/wrf/laps'
      yyyyddd = valid_yyyy*1000 + valid_jjj
      call wrf_date_to_ymd(yyyyddd, valid_yyyy, valid_mm, valid_dd)
      write (hdate, '(i4.4,"-",i2.2,"-",i2.2,"_",i2.2,":",i2.2,":00.0000")') &
         valid_yyyy, valid_mm, valid_dd, valid_hh, valid_min
      if (valid_min .eq. 0) then
         output_file_name = trim(output_prefix)//':'//hdate(1:13)
      else
         output_file_name = trim(output_prefix)//':'//hdate(1:16)
      end if
      !  open the file for sequential, unformatted output
      open (file=trim(output_file_name), &
            unit=output_unit, &
            form='unformatted', &
            status='unknown', &
            access='sequential')

      ! convert p levels from mb to pascals

      p_pa = p*100.

      ! set llflag based on grid type

      if (grid_type(1:8) .eq. 'mercator') then
         llflag = 1
      else if ((grid_type(1:24) .eq. 'secant lambert conformal') .or. &
               (grid_type(1:28) .eq. 'tangential lambert conformal')) then
         llflag = 3
      else if (grid_type(1:19) .eq. 'polar stereographic') then
         llflag = 5
      else
         print '(a,a,a)', 'unknown map projection: ', trim(grid_type), '.  i quit.'
         stop 'unknown_projection'
      end if

      print *, 'gribprep version =', gp_version
      print *, 'source = ', source
      print *, 'hdate = ', hdate
      print *, 'xfcst = ', xfcst
      print *, 'nx = ', x
      print *, 'ny = ', y
      print *, 'iproj = ', llflag
      print *, 'knownloc = ', knownloc
      print *, 'startlat = ', la1
      print *, 'startlon = ', lo1
      print *, 'dx = ', dx
      print *, 'dy = ', dy
      print *, 'xlonc = ', lov
      print *, 'truelat1 = ', latin1
      print *, 'truelat2 = ', latin2

      ! output temperature
      field = 't        '
      units = 'k                        '
      desc = 'temperature                                   '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      var_t: do k = 1, z3 + 1
         if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
            cycle var_t
         end if
         call write_gribprep_header(field, units, desc, p_pa(k))
         d2d = t(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f5.1,a,f5.1)', 'level (pa):', p_pa(k), ' min: ', &
            minval(d2d), ' max: ', maxval(d2d)
      end do var_t

      ! do u-component of wind
      field = 'u        '
      units = 'm s{-1}                  '
      desc = 'u-component of velocity, rotated to grid      '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      var_u: do k = 1, z3 + 1
         if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
            cycle var_u
         end if
         call write_gribprep_header(field, units, desc, p_pa(k))
         d2d = u(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f5.1,a,f5.1)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)
      end do var_u

      ! do v-component of wind
      field = 'v        '
      units = 'm s{-1}                  '
      desc = 'v-component of velocity, rotated to grid      '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      var_v: do k = 1, z3 + 1
         if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
            cycle var_v
         end if
         call write_gribprep_header(field, units, desc, p_pa(k))
         d2d = v(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f5.1,a,f5.1)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)
      end do var_v

      ! relative humidity
      field = 'rh       '
      units = '%                        '
      desc = 'relative humidity                             '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      var_rh: do k = 1, z3 + 1
         if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
            cycle var_rh
         end if
         call write_gribprep_header(field, units, desc, p_pa(k))
         d2d = rh(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f5.1,a,f5.1)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)
      end do var_rh

      ! do the heights
      field = 'hgt      '
      units = 'm                        '
      desc = 'geopotential height                           '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      var_ht: do k = 1, z3 + 1
         if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
            cycle var_ht
         end if
         call write_gribprep_header(field, units, desc, p_pa(k))
         d2d = ht(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f8.1,a,f8.1)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)
      end do var_ht

      ! terrain height
      field = 'soilhgt '
      units = 'm                        '
      desc = 'height of topography                          '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      call write_gribprep_header(field, units, desc, p_pa(z3 + 1))
      d2d = ht(:, :, z3 + 1)
      write (output_unit) d2d
      print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', p_pa(z3 + 1), ' min: ', minval(d2d), &
         ' max: ', maxval(d2d)

      ! skin temperature
      field = 'skintemp '
      units = 'k                        '
      desc = 'skin temperature                              '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      call write_gribprep_header(field, units, desc, p_pa(z3 + 1))
      write (output_unit) tskin
      print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', p_pa(z3 + 1), &
         ' min: ', minval(tskin), ' max: ', maxval(tskin)

      ! sea-level pressure field
      field = 'pmsl     '
      units = 'pa                       '
      desc = 'sea-level pressure                            '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      call write_gribprep_header(field, units, desc, slp_level)
      write (output_unit) slp
      print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', slp_level, ' min: ', minval(slp), &
         ' max: ', maxval(slp)

      ! surface pressure field
      field = 'psfc     '
      units = 'pa                       '
      desc = 'surface pressure                              '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      call write_gribprep_header(field, units, desc, p_pa(z3 + 1))
      write (output_unit) psfc
      print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', p_pa(z3 + 1), ' min: ', minval(psfc), &
         ' max: ', maxval(psfc)

      call get_r_missing_data(r_missing_data, istatus)
      print *, 'output_gribprep_format r_missing_data ', r_missing_data
      if (istatus .ne. 1) then
         write (6, *) ' bad status for r_missing_data'
         stop
      end if

      if (minval(snocov) .ge. 0) then
         ! water equivalent snow depth
         field = 'snowcovr '
         units = '(dimensionless)          '
         desc = 'snow cover flag                               '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         call write_gribprep_header(field, units, desc, p_pa(z3 + 1))

!   initialize output snow cover field to missing value
         d2d = -999.

!   convert from fraction to mask using namelist entry snow_thresh
         if (snow_thresh .le. 1.0) then
            where (snocov .ge. snow_thresh .and. snocov .ne. r_missing_data) d2d = 1.0
            where (snocov .lt. snow_thresh .and. snocov .ne. r_missing_data) d2d = 0.0
         end if

         write (output_unit) d2d
         print '(a,f9.1,a,f9.2,a,f9.2)', 'level (pa):', p_pa(z3 + 1), &
            ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)

      end if

      ! get cloud species if this is a hot start
      if (hotstart) then
         field = 'qliquid  '
         units = 'kg kg{-1}               '
         desc = 'cloud liquid water mixing ratio             '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_lwc: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_lwc
            end if
            call write_gribprep_header(field, units, desc, p_pa(k))
            d2d = lwc(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. from f8.6 to f9.6 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_lwc

         ! cloud ice
         field = 'qice     '
         units = 'kg kg{-1}               '
         desc = 'cloud ice mixing ratio                      '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_ice: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_ice
            end if
            call write_gribprep_header(field, units, desc, p_pa(k))
            d2d = ice(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. from f8.6 to f9.6 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_ice

         ! cloud rain
         field = 'qrain    '
         units = 'kg kg{-1}               '
         desc = 'rain water mixing ratio                     '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_rai: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_rai
            end if
            call write_gribprep_header(field, units, desc, p_pa(k))
            d2d = rai(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. from f8.6 to f9.6 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_rai

         ! snow
         field = 'qsnow    '
         units = 'kg kg{-1}               '
         desc = 'snow mixing ratio                           '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_sno: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_sno
            end if
            call write_gribprep_header(field, units, desc, p_pa(k))
            d2d = sno(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. from f8.6 to f9.6 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_sno

         ! graupel
         field = 'qgraupel '
         units = 'kg kg{-1}               '
         desc = 'graupel mixing ratio                        '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_pic: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_pic
            end if
            call write_gribprep_header(field, units, desc, p_pa(k))
            d2d = pic(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. from f8.6 to f9.6 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_pic

      end if

      close (output_unit)
      deallocate (d2d)
      deallocate (p_pa)
      return
   end subroutine output_gribprep_format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine write_gribprep_header(field, units, desc, level)

      ! writes the gribprep header given the filed, units, description, and level

      implicit none
      character(len=9), intent(in)  :: field
      character(len=25), intent(in)  :: units
      character(len=46), intent(in)  :: desc
      real, intent(in)              :: level

      write (output_unit) gp_version
      write (output_unit) hdate, xfcst, source, field, units, desc, level, x, y, llflag
      select case (llflag)
      case (1)
         write (output_unit) knownloc, la1, lo1, dx, dy, latin1
      case (3)
         write (output_unit) knownloc, la1, lo1, dx, dy, lov, latin1, latin2
      case (5)
         write (output_unit) knownloc, la1, lo1, dx, dy, lov, latin1
      end select

   end subroutine write_gribprep_header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine output_metgrid_format(p, t, ht, u, v, w, rh, slp, psfc, &
                                    lwc, rai, sno, ice, pic, snocov, tskin, soilt, soilm)

      !  subroutine of lapsprep that will build a file the
      !  wps format that can be read by metgrid

      implicit none

      ! arguments

      real, intent(in)                   :: p(:)        ! pressure (hpa)
      real, intent(in)                   :: t(:, :, :)    ! temperature (k)
      real, intent(in)                   :: ht(:, :, :)   ! height (m)
      real, intent(in)                   :: u(:, :, :)    ! u-wind (m s{-1})
      real, intent(in)                   :: v(:, :, :)    ! v-wind (m s{-1})
      real, intent(in)                   :: w(:, :, :)    ! w-wind (m s{-1})
      real, intent(in)                   :: rh(:, :, :)   ! relative humidity (%)
      real, intent(in)                   :: slp(:, :)    ! sea-level pressure (pa)
      real, intent(in)                   :: psfc(:, :)   ! surface pressure (pa)
      real, intent(in)                   :: lwc(:, :, :)  ! cloud liquid (kg/kg)
      real, intent(in)                   :: rai(:, :, :)  ! rain (kg/kg)
      real, intent(in)                   :: sno(:, :, :)  ! snow (kg/kg)
      real, intent(in)                   :: ice(:, :, :)  ! ice (kg/kg)
      real, intent(in)                   :: pic(:, :, :)  ! graupel (kg/kg)
      real, intent(in)                   :: soilt(:, :, :)! soil temp (k)
      real, intent(in)                   :: soilm(:, :, :)! soil moist (kg/kg)
      real, intent(in)                   :: snocov(:, :) ! snow cover (fract)
      real, intent(in)                   :: tskin(:, :)  ! skin temperature

      ! local variables

      integer            :: valid_mm, valid_dd
      character(len=256):: output_file_name
      real, allocatable  :: d2d(:, :)
      real, allocatable  :: p_pa(:)
      integer            :: k, yyyyddd
      integer            :: istatus
      real               :: r_missing_data

      ! yuanfu added laps soil moisture level depth information in centi-meters:
      integer :: lk
      real, parameter :: lsm3(3) = (/15.2, 30.5, 91.4/)
      real    :: alpha

      call get_r_missing_data(r_missing_data, istatus)
      print *, 'output_metgrid_format r_missing_data ', r_missing_data
      if (istatus .ne. 1) then
         write (6, *) ' bad status for r_missing_data'
         stop
      end if

      ! allocate a scratch 2d array
      allocate (d2d(x, y))
      allocate (p_pa(z3 + 1))
      ! build the output file name

      output_prefix = trim(laps_data_root)//'/lapsprd/lapsprep/wps/laps'
      yyyyddd = valid_yyyy*1000 + valid_jjj
      call wrf_date_to_ymd(yyyyddd, valid_yyyy, valid_mm, valid_dd)
      write (hdate, '(i4.4,"-",i2.2,"-",i2.2,"_",i2.2,":",i2.2,":00.0000")') &
         valid_yyyy, valid_mm, valid_dd, valid_hh, valid_min
      if (valid_min .eq. 0) then
         output_file_name = trim(output_prefix)//':'//hdate(1:13)
      else
         output_file_name = trim(output_prefix)//':'//hdate(1:16)
      end if
      !  open the file for sequential, unformatted output
      open (file=trim(output_file_name), &
            unit=output_unit, &
            form='unformatted', &
            status='unknown', &
            access='sequential')

      ! convert p levels from mb to pascals

      p_pa = p*100.

      ! set llflag based on grid type

      if (grid_type(1:8) .eq. 'mercator') then
         llflag = 1
      else if ((grid_type(1:24) .eq. 'secant lambert conformal') .or. &
               (grid_type(1:28) .eq. 'tangential lambert conformal')) then
         llflag = 3
      else if (grid_type(1:19) .eq. 'polar stereographic') then
         llflag = 5
      else
         print '(a,a,a)', 'unknown map projection: ', trim(grid_type), '.  i quit.'
         stop 'unknown_projection'
      end if

      print *, 'metgrid version =', gp_version_wps
      print *, 'source = ', source
      print *, 'hdate = ', hdate
      print *, 'xfcst = ', xfcst
      print *, 'nx = ', x
      print *, 'ny = ', y
      print *, 'iproj = ', llflag
      print *, 'knownloc = ', knownloc
      print *, 'startlat = ', la1
      print *, 'startlon = ', lo1
      print *, 'dx = ', dx
      print *, 'dy = ', dy
      print *, 'xlonc = ', lov
      print *, 'truelat1 = ', latin1
      print *, 'truelat2 = ', latin2

      ! output temperature
      field = 'tt       '
      units = 'k                        '
      desc = 'temperature                                   '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      var_t: do k = 1, z3 + 1
         if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
            cycle var_t
         end if
         call write_metgrid_header(field, units, desc, p_pa(k))
         d2d = t(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f5.1,a,f5.1)', 'level (pa):', p_pa(k), ' min: ', &
            minval(d2d), ' max: ', maxval(d2d)
      end do var_t

      ! do u-component of wind
      field = 'uu       '
      units = 'm s-1                    '
      desc = 'u                                             '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      var_u: do k = 1, z3 + 1
         if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
            cycle var_u
         end if
         call write_metgrid_header(field, units, desc, p_pa(k))
         d2d = u(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f5.1,a,f5.1)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)
      end do var_u

      ! do v-component of wind
      field = 'vv       '
      units = 'm s-1                    '
      desc = 'v                                             '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      var_v: do k = 1, z3 + 1
         if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
            cycle var_v
         end if
         call write_metgrid_header(field, units, desc, p_pa(k))
         d2d = v(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f5.1,a,f5.1)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)
      end do var_v

      ! do w-component of wind
      if (use_laps_vv) then
         field = 'vvel     '
         units = 'm s-1                   '
         desc = 'w                                             '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_vv: do k = 1, z3 + 1
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_vv
            end if
            call write_metgrid_header(field, units, desc, p_pa(k))
            d2d = w(:, :, k)
            write (output_unit) d2d
            print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_vv
      end if

      ! relative humidity
      field = 'rh       '
      units = '%                        '
      desc = 'relative humidity                             '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      var_rh: do k = 1, z3 + 1
         if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
            cycle var_rh
         end if
         call write_metgrid_header(field, units, desc, p_pa(k))
         d2d = rh(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f5.1,a,f5.1)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)
      end do var_rh

      ! do the heights
      field = 'hgt      '
      units = 'm                        '
      desc = 'geopotential height                           '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      var_ht: do k = 1, z3 + 1
         if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
            cycle var_ht
         end if
         call write_metgrid_header(field, units, desc, p_pa(k))
         d2d = ht(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f8.1,a,f8.1)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)
      end do var_ht

      ! terrain height
      field = 'soilhgt '
      units = 'm                        '
      desc = 'terrain field of source analysis              '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      call write_metgrid_header(field, units, desc, p_pa(z3 + 1))
      d2d = ht(:, :, z3 + 1)
      write (output_unit) d2d
      print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', p_pa(z3 + 1), ' min: ', minval(d2d), &
         ' max: ', maxval(d2d)

      ! surface fractional snow cover
      ! dividing by 6 translates a 0.3 threshold into .05 for wrf
      if (snow_thresh .le. 1.0) then
         field = 'snowh    '
         units = '                         '
         desc = 'snow cover (fraction)                         '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         call write_metgrid_header(field, units, desc, p_pa(z3 + 1))

!     initialize output snow cover field to missing value
         d2d = 0.

!     convert from fraction to mask using namelist entry snow_thresh
         where (snocov .ge. snow_thresh .and. snocov .ne. r_missing_data) d2d = 1.0
         where (snocov .lt. snow_thresh .and. snocov .ne. r_missing_data) d2d = 0.0

         write (output_unit) d2d
         print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', p_pa(z3 + 1), ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)
      end if

      ! skin temperature
      if (use_laps_skintemp) then
         field = 'skintemp '
         units = 'k                        '
         desc = 'skin temperature                              '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         call write_metgrid_header(field, units, desc, p_pa(z3 + 1))
         write (output_unit) tskin
         print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', p_pa(z3 + 1), &
            ' min: ', minval(tskin), ' max: ', maxval(tskin)
      end if

      ! soil temperature and moisture: by yuanfu xie 2015/05
      if (num_soil_layers .gt. 0) then
         field = 'soilt       '
         units = 'k                        '
         desc = 'soil temperature                                '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_soilt: do k = 1, num_soil_layers
            call write_metgrid_header(field, units, desc, soil_layer_depths(k))
            d2d = tskin  ! test of skin temperature for now
            write (output_unit) d2d
            print '(a,f9.1,a,f5.1,a,f5.1)', 'level (cm):', soil_layer_depths(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
            print *, 'yuanfu debugger: ', d2d(1, 1), tskin(1, 1), x, y
         end do var_soilt

         field = 'soilm       '
         units = 'kg/kg                    '
         desc = 'soil moisture                                    '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_soilm: do k = 1, num_soil_layers
            call write_metgrid_header(field, units, desc, soil_layer_depths(k))

            ! test soil moisture from laps lm1: by yuanfu xie
            d2d = 12.3  ! test
            print *, 'soil_layer_depths: ', soil_layer_depths(1:4), lsm3(1:3), maxval(soilm)
            if (soil_layer_depths(k) .le. lsm3(3)) then
               ! interpolated from lm1:
               if (soil_layer_depths(k) .le. lsm3(1)) then
                  alpha = 1.0
                  lk = 1
               else if (soil_layer_depths(k) .le. lsm3(2)) then
                  alpha = (lsm3(2) - soil_layer_depths(k))/(lsm3(2) - lsm3(1))
                  lk = 1
               else
                  alpha = (lsm3(3) - soil_layer_depths(k))/(lsm3(3) - lsm3(2))
                  lk = 2
               end if
               d2d = alpha*soilm(:, :, lk) + (1.0 - alpha)*soilm(:, :, lk + 1)
            else
               ! set to the deepest level of laps lm1
               d2d = soilm(:, :, 3)
            end if

            write (output_unit) d2d/100.0  ! temporarily change the scale to 0-1
            print '(a,f9.1,a,f5.1,a,f5.1)', 'level (cm):', soil_layer_depths(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
            print *, 'xie debugger: ', d2d(1, 1)
         end do var_soilm
      end if

      ! sea-level pressure field
      field = 'pmsl     '
      units = 'pa                       '
      desc = 'sea-level pressure                            '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      call write_metgrid_header(field, units, desc, slp_level)
      print *, 'xxxx : ', slp(1, 198)
      write (output_unit) slp
      print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', slp_level, ' min: ', minval(slp), &
         ' max: ', maxval(slp)

      ! surface pressure field
      field = 'psfc     '
      units = 'pa                       '
      desc = 'surface pressure                              '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      call write_metgrid_header(field, units, desc, p_pa(z3 + 1))
      print *, 'yyyy : ', psfc(1, 198)
      write (output_unit) psfc
      print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', p_pa(z3 + 1), ' min: ', minval(psfc), &
         ' max: ', maxval(psfc)

      ! get cloud species if this is a hot start
      if (hotstart) then
         field = 'qc       '
         units = 'kg kg-1                 '
         desc = 'cloud liquid water mixing ratio             '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_lwc: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_lwc
            end if
            call write_metgrid_header(field, units, desc, p_pa(k))
            d2d = lwc(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_lwc

         ! cloud ice
         field = 'qi       '
         units = 'kg kg-1                 '
         desc = 'cloud ice mixing ratio                      '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_ice: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_ice
            end if
            call write_metgrid_header(field, units, desc, p_pa(k))
            d2d = ice(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_ice

         ! cloud rain
         field = 'qr       '
         units = 'kg kg-1                 '
         desc = 'rain water mixing ratio                     '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_rai: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_rai
            end if
            call write_metgrid_header(field, units, desc, p_pa(k))
            d2d = rai(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_rai

         ! snow
         field = 'qs       '
         units = 'kg kg-1                 '
         desc = 'snow mixing ratio                           '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_sno: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_sno
            end if
            call write_metgrid_header(field, units, desc, p_pa(k))
            d2d = sno(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_sno

         ! graupel
         field = 'qg       '
         units = 'kg kg-1                 '
         desc = 'graupel mixing ratio                        '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_pic: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_pic
            end if
            call write_metgrid_header(field, units, desc, p_pa(k))
            d2d = pic(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_pic

      end if

      close (output_unit)
      deallocate (d2d)
      deallocate (p_pa)
      return
   end subroutine output_metgrid_format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine write_metgrid_header(field, units, desc, level)

      ! writes the gribprep header given the filed, units, description, and level

      implicit none
      character(len=9), intent(in)  :: field
      character(len=25), intent(in)  :: units
      character(len=46), intent(in)  :: desc
      real, intent(in)              :: level
      logical, parameter            :: is_wind_grid_rel = .true.
      real                          :: radius_of_earth_m, radius_of_earth_km
      integer                       :: istatus

      call get_earth_radius(radius_of_earth_m, istatus)
      radius_of_earth_km = radius_of_earth_m/1000.
      print *, 'write_metgrid_header radius_of_earth_km ', radius_of_earth_km
      if (istatus .ne. 1) then
         write (6, *) ' bad status for radius_of_earth_km'
         stop
      end if

      write (output_unit) gp_version_wps
      write (output_unit) hdate, xfcst, source, field, units, desc, level, x, y, llflag
      select case (llflag)
      case (1)
         write (output_unit) knownloc, la1, lo1, dx, dy, latin1, radius_of_earth_km
      case (3)
         write (output_unit) knownloc, la1, lo1, dx, dy, lov, latin1, latin2, radius_of_earth_km
      case (5)
         write (output_unit) knownloc, la1, lo1, dx, dy, lov, latin1, radius_of_earth_km
      end select
      write (output_unit) is_wind_grid_rel

   end subroutine write_metgrid_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lapsprep_wrf

