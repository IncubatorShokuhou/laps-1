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

module wfoprep_mm5

! purpose
! =======
! module to contain the various output routines needed for lapsprep
! to support initializition of the mm5v3 model.
!
! subroutines contained
! =====================
! output_mm5v3_basic     - outputs state variables on pressure surfaces
!                          plus mslp and topography
! output_mm5v3_sfc       - outputs the surface fields
! write_pregrid_header   - writes pregrid headers
!
! remarks
! =======
!
!
! history
! =======
! 1 oct 2001 -- original -- brent shaw

   use map_utils
   implicit none

   integer, parameter :: pg_version = 3
   integer, parameter :: output_unit = 78
   public output_mm5v3_basic, output_mm5v3_sfc
   integer, parameter, public   :: mm5mode_new = 1
   integer, parameter, public   :: mm5mode_append = 2
   logical, parameter           :: verbose = .true.
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine output_mm5v3_basic(i4time_cycle, i4time_valid, proj, &
                                 np_ht, np_t, np_u, np_v, np_rh, &
                                 ht_plevels, t_plevels, u_plevels, v_plevels, rh_plevels, &
                                 z, t, u, v, rh, mslp, topo, &
                                 ext_data_path, output_name, mode, istatus)

      implicit none
      integer, intent(in)         :: i4time_cycle
      integer, intent(in)         :: i4time_valid
      type(proj_info), intent(in)  :: proj
      integer, intent(in)         :: np_ht
      integer, intent(in)         :: np_t
      integer, intent(in)         :: np_u
      integer, intent(in)         :: np_v
      integer, intent(in)         :: np_rh
      real, intent(in)            :: ht_plevels(np_ht)
      real, intent(in)            :: t_plevels(np_t)
      real, intent(in)            :: u_plevels(np_u)
      real, intent(in)            :: v_plevels(np_v)
      real, intent(in)            :: rh_plevels(np_rh)
      real, intent(in)            :: z(:, :, :)
      real, intent(in)            :: t(:, :, :)
      real, intent(in)            :: u(:, :, :)
      real, intent(in)            :: v(:, :, :)
      real, intent(in)            :: rh(:, :, :)
      real, intent(in)            :: mslp(:, :)
      real, intent(in)            :: topo(:, :)
      character(len=256), intent(in):: ext_data_path
      character(len=32), intent(in):: output_name
      integer, intent(in)         :: mode
      integer, intent(out)        :: istatus

      integer                     :: k
      character(len=256)          :: outfile
      character(len=24)           :: atime
      character(len=3)            :: amonth
      character(len=2)            :: amonth_num
      integer                     :: m
      character(len=24) :: hdate
      character(len=9)  :: field
      character(len=25) :: units
      character(len=46) :: desc
      real               :: xfcst
      real, parameter    :: slp_level = 201300.0
      real, parameter    :: sfc_level = 200100.0
      istatus = 1

      call make_pregrid_filename(ext_data_path, output_name, &
                                 i4time_valid, outfile)

      ! compute xfcst
      xfcst = float(i4time_valid - i4time_cycle)/3600.

      ! compute hdate
      call make_hdate_from_i4time(i4time_valid, hdate)

      ! open the file, using the mode dependency
      print *, 'opening file: ', trim(outfile)
      if (mode .eq. mm5mode_new) then
         open (file=trim(outfile), &
               unit=output_unit, &
               form='unformatted', &
               status='replace', &
               access='sequential')
      else if (mode .eq. mm5mode_append) then
         open (file=trim(outfile), &
               unit=output_unit, &
               form='unformatted', &
               status='unknown', &
               access='sequential', &
               position='append')
      else
         print *, 'uknown open mode for mm5v3 output: ', mode
         istatus = 0
         return
      end if

      if (verbose) then
         print *, 'pregrid version =', pg_version
         print *, 'hdate = ', hdate
         print *, 'xfcst = ', xfcst
         print *, 'nx = ', proj%nx
         print *, 'ny = ', proj%ny
         print *, 'iproj = ', proj%code
         print *, 'startlat = ', proj%lat1
         print *, 'startlon = ', proj%lon1
         print *, 'dx = ', proj%dx*0.001
         print *, 'dy = ', proj%dx*0.001
         print *, 'xlonc = ', proj%stdlon
         print *, 'truelat1 = ', proj%truelat1
         print *, 'truelat2 = ', proj%truelat2
      end if

      ! output temperature
      field = 't        '
      units = 'k                        '
      desc = 'temperature                                   '
      var_t: do k = 1, np_t
         call write_pregrid_header(proj, field, units, desc, t_plevels(k), hdate, xfcst)
         write (output_unit) t(:, :, k)
         if (verbose) then
            print '(a,1x,a,a,f9.1,a,f5.1,a,f5.1)', field, units, &
               'level (pa):', t_plevels(k), ' min: ', &
               minval(t(:, :, k)), ' max: ', maxval(t(:, :, k))
         end if
      end do var_t

      ! do u-component of wind
      field = 'u        '
      units = 'm s{-1}                  '
      desc = 'u-component of velocity, rotated to grid      '
      var_u: do k = 1, np_u
         call write_pregrid_header(proj, field, units, desc, u_plevels(k), hdate, xfcst)
         write (output_unit) u(:, :, k)
         if (verbose) then
            print '(a,1x,a,a,f9.1,a,f5.1,a,f5.1)', field, units, &
               'level (pa):', u_plevels(k), &
               ' min: ', minval(u(:, :, k)), &
               ' max: ', maxval(u(:, :, k))
         end if
      end do var_u

      ! do v-component of wind
      field = 'v        '
      units = 'm s{-1}                  '
      desc = 'v-component of velocity, rotated to grid      '
      var_v: do k = 1, np_v
         call write_pregrid_header(proj, field, units, desc, v_plevels(k), hdate, xfcst)
         write (output_unit) v(:, :, k)
         if (verbose) then
            print '(a,1x,a,a,f9.1,a,f5.1,a,f5.1)', field, units, &
               'level (pa):', v_plevels(k), &
               ' min: ', minval(v(:, :, k)), &
               ' max: ', maxval(v(:, :, k))
         end if
      end do var_v

      ! relative humidity
      field = 'rh       '
      units = '%                        '
      desc = 'relative humidity                             '
      var_rh: do k = 1, np_rh
         call write_pregrid_header(proj, field, units, desc, rh_plevels(k), hdate, xfcst)
         write (output_unit) rh(:, :, k)
         if (verbose) then
            print '(a,1x,a,a,f9.1,a,f5.1,a,f5.1)', field, units, &
               'level (pa):', rh_plevels(k), &
               ' min: ', minval(rh(:, :, k)), &
               ' max: ', maxval(rh(:, :, k))
         end if
      end do var_rh

      ! do the heights
      field = 'hgt      '
      units = 'm                        '
      desc = 'geopotential height                           '
      var_ht: do k = 1, np_ht
         call write_pregrid_header(proj, field, units, desc, ht_plevels(k), hdate, xfcst)
         write (output_unit) z(:, :, k)
         if (verbose) then
            print '(a,x,a,a,f9.1,a,f8.1,a,f8.1)', field, units, &
               'level (pa):', ht_plevels(k), &
               ' min: ', minval(z(:, :, k)), &
               ' max: ', maxval(z(:, :, k))
         end if
      end do var_ht

      ! terrain height
      field = 'soilhgt  '
      units = 'm                        '
      desc = 'height of topography                          '
      call write_pregrid_header(proj, field, units, desc, sfc_level, hdate, xfcst)
      write (output_unit) topo
      if (verbose) then
         print '(a,x,a,a,f9.1,a,f9.1,a,f9.1)', field, units, &
            'level (pa):', sfc_level, &
            ' min: ', minval(topo), &
            ' max: ', maxval(topo)
      end if
      ! sea-level pressure field
      field = 'pmsl     '
      units = 'pa                       '
      desc = 'sea-level pressure                            '
      call write_pregrid_header(proj, field, units, desc, slp_level, hdate, xfcst)
      write (output_unit) mslp
      if (verbose) then
         print '(a,x,a,a,f9.1,a,f9.1,a,f9.1)', field, units, &
            'level (pa):', slp_level, &
            ' min: ', minval(mslp), &
            ' max: ', maxval(mslp)
      end if
      close (output_unit)
      return
   end subroutine output_mm5v3_basic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine output_mm5v3_sfc(i4time_cycle, i4time_valid, proj, &
                               tsf, usf, vsf, rhsf, &
                               ext_data_path, output_name, mode, istatus)

      implicit none

      integer, intent(in)                     :: i4time_cycle
      integer, intent(in)                     :: i4time_valid
      type(proj_info), intent(in)             :: proj
      real, intent(in)                        :: tsf(:, :)
      real, intent(in)                        :: usf(:, :)
      real, intent(in)                        :: vsf(:, :)
      real, intent(in)                        :: rhsf(:, :)
      character(len=256), intent(in)          :: ext_data_path
      character(len=32), intent(in)           :: output_name
      integer, intent(in)                     :: mode
      integer, intent(out)                    :: istatus

      character(len=256)                     :: outfile

      character(len=24) :: hdate
      character(len=9)  :: field
      character(len=25) :: units
      character(len=46) :: desc
      real               :: xfcst
      real, parameter    :: slp_level = 201300.0
      real, parameter    :: sfc_level = 200100.0
      istatus = 1

      call make_pregrid_filename(ext_data_path, output_name, &
                                 i4time_valid, outfile)
      ! compute xfcst
      xfcst = float(i4time_valid - i4time_cycle)/3600.

      ! compute hdate
      call make_hdate_from_i4time(i4time_valid, hdate)

      ! open the file, using the mode dependency
      print *, 'opening file: ', trim(outfile)
      if (mode .eq. mm5mode_new) then
         open (file=trim(outfile), &
               unit=output_unit, &
               form='unformatted', &
               status='replace', &
               access='sequential')
      else if (mode .eq. mm5mode_append) then
         open (file=trim(outfile), &
               unit=output_unit, &
               form='unformatted', &
               status='unknown', &
               access='sequential', &
               position='append')
      else
         print *, 'uknown open mode for mm5v3 output: ', mode
         istatus = 0
         return
      end if

      ! output temperature
      field = 't        '
      units = 'k                        '
      desc = 'temperature                                   '
      call write_pregrid_header(proj, field, units, desc, sfc_level, hdate, xfcst)
      write (output_unit) tsf
      if (verbose) then
         print '(a,x,a,a,f9.1,a,f5.1,a,f5.1)', field, units, &
            'level (pa):', sfc_level, ' min: ', &
            minval(tsf), ' max: ', maxval(tsf)
      end if

      ! output u-wind component
      field = 'u        '
      units = 'm s{-1}                  '
      desc = 'u-component of wind, rotated to grid          '
      call write_pregrid_header(proj, field, units, desc, sfc_level, hdate, xfcst)
      write (output_unit) usf
      if (verbose) then
         print '(a,x,a,a,f9.1,a,f5.1,a,f5.1)', field, units, &
            'level (pa):', sfc_level, ' min: ', &
            minval(usf), ' max: ', maxval(usf)
      end if

      ! output v-wind component
      field = 'v        '
      units = 'm s{-1}                  '
      desc = 'v-component of wind, rotated to grid          '
      call write_pregrid_header(proj, field, units, desc, sfc_level, hdate, xfcst)
      write (output_unit) vsf
      if (verbose) then
         print '(a,x,a,a,f9.1,a,f5.1,a,f5.1)', field, units, &
            'level (pa):', sfc_level, ' min: ', &
            minval(vsf), ' max: ', maxval(vsf)
      end if

      ! output relative humidity
      field = 'rh       '
      units = '%                        '
      desc = 'relative humidity                             '
      call write_pregrid_header(proj, field, units, desc, sfc_level, hdate, xfcst)
      write (output_unit) rhsf
      if (verbose) then
         print '(a,x,a,a,f9.1,a,f5.1,a,f5.1)', field, units, &
            'level (pa):', sfc_level, ' min: ', &
            minval(rhsf), ' max: ', maxval(rhsf)
      end if

      close (output_unit)
      return
   end subroutine output_mm5v3_sfc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine output_mm5v3_fsst(i4time_cycle, i4time_valid, proj, &
                                fsst, ext_data_path, output_name, mode, istatus)

      implicit none

      integer, intent(in)                     :: i4time_cycle
      integer, intent(in)                     :: i4time_valid
      type(proj_info), intent(in)             :: proj
      real, intent(in)                        :: fsst(:, :)
      character(len=256), intent(in)          :: ext_data_path
      character(len=32), intent(in)           :: output_name
      integer, intent(in)                     :: mode
      integer, intent(out)                    :: istatus

      character(len=256)                     :: outfile

      character(len=24) :: hdate
      character(len=9)  :: field
      character(len=25) :: units
      character(len=46) :: desc
      real               :: xfcst
      real, parameter    :: sfc_level = 200100.0
      istatus = 1

      outfile = trim(ext_data_path)//trim(output_name)// &
                ':fsst'

      ! compute xfcst
      xfcst = float(i4time_valid - i4time_cycle)/3600.

      ! compute hdate
      call make_hdate_from_i4time(i4time_valid, hdate)
      ! open the file, using the mode dependency

      print *, 'opening file: ', trim(outfile)
      if (mode .eq. mm5mode_new) then
         open (file=trim(outfile), &
               unit=output_unit, &
               form='unformatted', &
               status='replace', &
               access='sequential')
      else if (mode .eq. mm5mode_append) then
         open (file=trim(outfile), &
               unit=output_unit, &
               form='unformatted', &
               status='unknown', &
               access='sequential', &
               position='append')
      else
         print *, 'uknown open mode for mm5v3 output: ', mode
         istatus = 0
         return
      end if

      ! output temperature
      field = 'sst      '
      units = 'k                        '
      desc = 'sea-surface temperature (faked)               '
      call write_pregrid_header(proj, field, units, desc, sfc_level, hdate, xfcst)
      write (output_unit) fsst
      if (verbose) then
         print '(a,x,a,a,f9.1,a,f5.1,a,f5.1)', field, units, &
            'level (pa):', sfc_level, ' min: ', &
            minval(fsst), ' max: ', maxval(fsst)
      end if
      close (output_unit)
      return
   end subroutine output_mm5v3_fsst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine make_pregrid_filename(ext_data_path, output_name, i4time, filename)

      implicit none
      character(len=256), intent(in)           :: ext_data_path
      character(len=32), intent(in)           :: output_name
      integer, intent(in)                     :: i4time
      character(len=256), intent(out)          :: filename
      character(len=24)                       :: hdate

      call make_hdate_from_i4time(i4time, hdate)

      ! build the output file name

      filename = trim(ext_data_path)//trim(output_name)// &
                 ':'//hdate(1:13)

      return

   end subroutine make_pregrid_filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine make_hdate_from_i4time(i4time, hdate)

      implicit none
      integer, intent(in)                    :: i4time
      character(len=24), intent(out)          :: hdate

      character(len=3)                       :: amonth(12)
      character(len=24)                      :: atime
      integer                                :: istatus
      integer                                :: m
      integer                                :: dom
      character(len=2)                       :: amonth_num

      data amonth/'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
         'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/

      ! build the file name by converting the i4time_valid
      ! to yyyy-mm-dd_hh format

      call cv_i4tim_asc_lp(i4time, atime, istatus)
      if (istatus .ne. 1) then
         print *, 'problem converting i4time!'
         return
      end if
      findmonth: do m = 1, 12
         if (atime(4:6) .eq. amonth(m)) exit findmonth
      end do findmonth

      ! ensure we make day into a 2-digit value
      read (atime(1:2), '(i2)') dom
      write (atime(1:2), '(i2.2)') dom
      write (amonth_num, '(i2.2)') m
      hdate = atime(8:11)//'-'// &
              amonth_num//'-'// &
              atime(1:2)//'_'// &
              atime(13:23)//'00'
      return
   end subroutine make_hdate_from_i4time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_pregrid_header(proj, field, units, desc, level, hdate, xfcst)

      ! writes the gribprep header given the filed, units, description, and level

      implicit none
      type(proj_info), intent(in)    :: proj
      character(len=9), intent(in)  :: field
      character(len=25), intent(in)  :: units
      character(len=46), intent(in)  :: desc
      real, intent(in)              :: level
      character(len=24), intent(in)  :: hdate
      real, intent(in)               :: xfcst

      write (output_unit) pg_version
      write (output_unit) hdate, xfcst, field, units, desc, level, proj%nx, &
         proj%ny, proj%code
      select case (proj%code)
      case (0)
         write (output_unit) proj%lat1, proj%lon1, proj%dlat, proj%dlon
      case (1)
         write (output_unit) proj%lat1, proj%lon1, proj%dx*0.001, proj%dx*0.001, &
            proj%truelat1
      case (3)
         write (output_unit) proj%lat1, proj%lon1, proj%dx*0.001, proj%dx*0.001, &
            proj%stdlon, proj%truelat1, proj%truelat2
      case (5)
         write (output_unit) proj%lat1, proj%lon1, proj%dx*0.001, proj%dx*0.001, &
            proj%stdlon, proj%truelat1
      end select

   end subroutine write_pregrid_header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module wfoprep_mm5

