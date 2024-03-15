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

module lapsprep_mm5

! purpose
! =======
! module to contain the various output routines needed for lapsprep
! to support initializition of the mm5v3 model.
!
! subroutines contained
! =====================
! output_pregrid_format  - used to support mm5 initializations
! write_pregrid_header   - writes pregrid headers
!
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
   implicit none

   private
   integer, parameter :: pg_version = 3
   real, parameter    :: xfcst = 0.
   character(len=24) :: hdate
   integer            :: llflag
   character(len=9)  :: field
   character(len=25) :: units
   character(len=46) :: desc
   integer, parameter :: output_unit = 78
   real, parameter    :: slp_level = 201300.0

   public output_pregrid_format
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine output_pregrid_format(p, t, ht, u, v, rh, slp, &
                                    lwc, rai, sno, ice, pic, snocov, tskin)

      !  subroutine of lapsprep that will build a file in the
      !  mm5v3 pregrid format that can be read by regrid

      implicit none

      ! arguments

      real, intent(in)                   :: p(:)        ! pressure (hpa)
      real, intent(in)                   :: t(:, :, :)    ! temperature (k)
      real, intent(in)                   :: ht(:, :, :)   ! height (m)
      real, intent(in)                   :: u(:, :, :)    ! u-wind (m s{-1})
      real, intent(in)                   :: v(:, :, :)    ! v-wind (m s{-1})
      real, intent(in)                   :: rh(:, :, :)   ! relative humidity (%)
      real, intent(in)                   :: slp(:, :)    ! sea-level pressure (pa)
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
      character(len=256)  :: sstfile
      real, allocatable  :: d2d(:, :)
      real, allocatable  :: p_pa(:)
      integer            :: k, yyyyddd

      ! allocate a scratch 2d array
      allocate (d2d(x, y))
      allocate (p_pa(z3 + 1))
      ! build the output file name

      output_prefix = trim(laps_data_root)//'/lapsprd/lapsprep/mm5/laps'
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

      print *, 'pregrid version =', pg_version
      print *, 'hdate = ', hdate
      print *, 'xfcst = ', xfcst
      print *, 'nx = ', x
      print *, 'ny = ', y
      print *, 'iproj = ', llflag
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
         call write_pregrid_header(field, units, desc, p_pa(k))
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
         call write_pregrid_header(field, units, desc, p_pa(k))
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
         call write_pregrid_header(field, units, desc, p_pa(k))
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
         call write_pregrid_header(field, units, desc, p_pa(k))
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
         call write_pregrid_header(field, units, desc, p_pa(k))
         d2d = ht(:, :, k)
         write (output_unit) d2d
         print '(a,f9.1,a,f8.1,a,f8.1)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)
      end do var_ht

      ! terrain height
      field = 'hgt      '
      units = 'm                        '
      desc = 'height of topography                          '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      call write_pregrid_header(field, units, desc, p_pa(z3 + 1))
      d2d = ht(:, :, z3 + 1)
      write (output_unit) d2d
      print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', p_pa(z3 + 1), ' min: ', minval(d2d), &
         ' max: ', maxval(d2d)

      ! sea-level pressure field
      field = 'pmsl     '
      units = 'pa                       '
      desc = 'sea-level pressure                            '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      call write_pregrid_header(field, units, desc, slp_level)
      write (output_unit) slp
      print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', slp_level, ' min: ', minval(slp), &
         ' max: ', maxval(slp)

      ! snow cover
      if ((minval(snocov) .ge. 0.) .and. (maxval(snocov) .lt. 1.1)) then
         ! water equivalent snow depth
         field = 'snowcovr '
         units = '(dimensionless)          '
         desc = 'snow cover flag                               '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         call write_pregrid_header(field, units, desc, p_pa(z3 + 1))

         ! convert the snow cover fraction (ranges from 0 -> 1.) into a flag
         ! value using the snow_thresh parameter set in the namelist.
         d2d(:, :) = 0.
         where (snocov .ge. snow_thresh) d2d = 1.

         write (output_unit) d2d
         print '(a,f9.1,a,f9.2,a,f9.2)', 'level (pa):', p_pa(z3 + 1), &
            ' min: ', minval(d2d), &
            ' max: ', maxval(d2d)

      end if

      ! get cloud species if this is a hot start
      if (hotstart) then
         field = 'clw      '
         units = 'kg kg{-1}               '
         desc = 'cloud liquid water mixing ratio             '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_lwc: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_lwc
            end if
            call write_pregrid_header(field, units, desc, p_pa(k))

            d2d = lwc(:, :, k)

            write (output_unit) d2d
!hj: w>=d+3. from f8.6 to f9.6 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_lwc

         ! cloud ice
         field = 'ice      '
         units = 'kg kg{-1}               '
         desc = 'cloud ice mixing ratio                      '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_ice: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_ice
            end if
            call write_pregrid_header(field, units, desc, p_pa(k))
            d2d = ice(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. from f8.6 to f9.6 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_ice

         ! cloud rain
         field = 'rnw      '
         units = 'kg kg{-1}               '
         desc = 'rain water mixing ratio                     '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_rai: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_rai
            end if
            call write_pregrid_header(field, units, desc, p_pa(k))
            d2d = rai(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. from f8.6 to f9.6 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_rai

         ! snow
         field = 'snow     '
         units = 'kg kg{-1}               '
         desc = 'snow mixing ratio                           '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_sno: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_sno
            end if
            call write_pregrid_header(field, units, desc, p_pa(k))
            d2d = sno(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. from f8.6 to f9.6 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_sno

         ! graupel
         field = 'graupel  '
         units = 'kg kg{-1}               '
         desc = 'graupel mixing ratio                        '
         print *, 'field = ', field
         print *, 'units = ', units
         print *, 'desc =  ', desc
         var_pic: do k = 1, z3
            if ((p_pa(k) .gt. 100100) .and. (p_pa(k) .lt. 200000)) then
               cycle var_pic
            end if
            call write_pregrid_header(field, units, desc, p_pa(k))
            d2d = pic(:, :, k)
            write (output_unit) d2d
!hj: w>=d+3. from f8.6 to f9.6 10/14/2013
            print '(a,f9.1,a,f9.6,a,f9.6)', 'level (pa):', p_pa(k), ' min: ', minval(d2d), &
               ' max: ', maxval(d2d)
         end do var_pic

      end if

      close (output_unit)

      ! change for rsa to write skin temp out as a laps:tskin file
      ! skin temperature field
      sstfile = trim(output_prefix)//':tskin'
      open (unit=output_unit, file=sstfile, form='unformatted', access='sequential')
      field = 'skintemp '
      units = 'k                        '
      desc = 'skin temperature                              '
      print *, 'field = ', field
      print *, 'units = ', units
      print *, 'desc =  ', desc
      call write_pregrid_header(field, units, desc, p_pa(z3 + 1))
      write (output_unit) tskin
      print '(a,f9.1,a,f9.1,a,f9.1)', 'level (pa):', p_pa(z3 + 1), ' min: ', &
         minval(tskin), ' max: ', maxval(tskin)
      close (output_unit)

      deallocate (d2d)
      deallocate (p_pa)
      return
   end subroutine output_pregrid_format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine write_pregrid_header(field, units, desc, level)

      ! writes the gribprep header given the filed, units, description, and level

      implicit none
      character(len=9), intent(in)  :: field
      character(len=25), intent(in)  :: units
      character(len=46), intent(in)  :: desc
      real, intent(in)              :: level

      write (output_unit) pg_version
      write (output_unit) hdate, xfcst, field, units, desc, level, x, y, llflag
      select case (llflag)
      case (1)
         write (output_unit) la1, lo1, dx, dy, latin1
      case (3)
         write (output_unit) la1, lo1, dx, dy, lov, latin1, latin2
      case (5)
         write (output_unit) la1, lo1, dx, dy, lov, latin1
      end select

   end subroutine write_pregrid_header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lapsprep_mm5

