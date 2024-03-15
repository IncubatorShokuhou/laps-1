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

module lapsprep_rams

! purpose
! =======
! module to contain the various output routines needed for lapsprep.
!
! subroutines contained
! =====================
! output_ralph2_format  - used to support rams 4.x initializations
!
! remarks
! =======

!
! history
! =======
! 28 nov 2000 -- original -- brent shaw

   use setup
   use laps_static
   use date_pack

   private
   public output_ralph2_format
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine output_ralph2_format(p, u, v, t, ht, rh, slp, psfc, snocov, tskin)

      ! subroutine to output data in rams 4.x "ralph 2" format. the skin temp from
      ! laps (lsx/tgd) is used for sst.

      ! note that the p array has a bogus 2001 mb value as the last entry
      ! to designate the surface for mm5 applications.  the surface values
      ! of the state variables are contained in the last layer of their
      ! respective 3d arrays.

      ! the u/v winds are grid-relative, and for now ralph 2 expects
      ! true winds, so this routine must rotate them to true.

      ! surface and msl pressure must be converted to mb for rams.

      ! p levels must be written as integer values

      ! rh must be written as a fraction

      ! 3d variables are stored top down, but rams wants them
      ! bottom up.

      ! note that ralph only supports lat/lon or lambert projections.
      ! however, lambert can be used to fake a polar stereographic by
      ! setting both standard lats equal to the pole. however, if the laps
      ! grid is a mercator, then you are out of luck.

      implicit none
      real, intent(in)               :: p(:)      !pressure levels in mb
      real, intent(in)               :: u(:, :, :)  !u-component of wind wrt grid in m/s
      real, intent(in)               :: v(:, :, :)  !v-component of wind wrf grid in m/s
      real, intent(in)               :: t(:, :, :)  !temperature in k
      real, intent(in)               :: ht(:, :, :) !geopotential height in m
      real, intent(in)               :: rh(:, :, :) !relative humidity in %
      real, intent(in)               :: slp(:, :)  !msl pressure in pa
      real, intent(in)               :: psfc(:, :) !surface pressure in pa
      real, intent(in)               :: snocov(:, :)! snow cover (fract)
      real, intent(in)               :: tskin(:, :) ! skin temp
      ! local variables

      integer, parameter             :: marker = 999999
      integer, parameter             :: version = 2
      integer, parameter             :: vt_increment = 0
      integer, parameter             :: latlon_proj = 1
      integer, parameter             :: lambert_proj = 2
      integer, parameter             :: level_flag = 1
      integer, parameter             :: output_unit = 78
      integer                        :: proj_flag
      integer                        :: yyyyddd, valid_mm, valid_dd
      integer                        :: i, j, k
      integer, allocatable            :: p_int(:)
      real                           :: latin1_out
      real                           :: latin2_out
      real, allocatable              :: ut(:, :)
      real, allocatable              :: vt(:, :)
      character(len=256)            :: output_file_name
      character(len=15)             :: date_string

      ! build the output file name

      output_prefix = trim(laps_data_root)//'/lapsprd/lapsprep/rams/laps'
      yyyyddd = valid_yyyy*1000 + valid_jjj
      call wrf_date_to_ymd(yyyyddd, valid_yyyy, valid_mm, valid_dd)
      write (date_string, '(i4.4,"-",i2.2,"-",i2.2,"-",i2.2,i2.2)') &
         valid_yyyy, valid_mm, valid_dd, valid_hh, valid_min
      output_file_name = trim(output_prefix)//':'//date_string

      !  open the file for sequential, unformatted output
      open (file=trim(output_file_name), &
            unit=output_unit, &
            form='formatted', &
            status='unknown', &
            access='sequential')
      rewind (output_unit)
      ! write the header record as described in the ralph 2 documentation

      ! record 1 - header
      write (output_unit, '(i6.6,i4)') marker, version

      ! record 2 - time and dimension information

      write (output_unit, '(i4.4,2i4,i6,4i4)') valid_yyyy, valid_mm, &
         valid_dd, valid_hh*100 + valid_min, &
         vt_increment, z3, x, y

      ! record 3 - projection info
      if (grid_type(1:8) .eq. 'latlon') then
         proj_flag = 1
         latin1_out = latin1
         latin2_out = latin2
      else if (grid_type(1:24) .eq. 'secant lambert conformal') then
         proj_flag = 2
         latin1_out = latin1
         latin2_out = latin2
      else if (grid_type(1:19) .eq. 'polar stereographic') then
         if (abs(latin2) .ne. 90.) then
            print '(a)', 'this is a local stereographic, so i quit.'
            stop 'unsupported projection'
         else
            proj_flag = 3
            latin1_out = latin1
            latin2_out = latin2
         end if
      else
         print '(a,a,a)', 'ralph2 unsupported map projection: ', &
            trim(grid_type), '.  i quit.'
         stop 'unsupported projection'
      end if
      write (output_unit, '(i1,2f10.1,7f10.3)') proj_flag, dx*1000, dy*1000, &
         la1, lo1, la2, lo2, &
         latin1_out, lov, &
         latin2_out

      ! record 4 - vertical level specification
      allocate (p_int(z3 + 1))
      p_int = nint(p)
      write (output_unit, *) level_flag, p_int(z3:1:-1)
      deallocate (p_int)

      allocate (ut(x, y)) ! array for true u-winds
      allocate (vt(x, y)) ! array for true v-winds
      ! now write out the pressure level data, looping by level

      level_loop: do k = z3, 1, -1

         ! rotate the u and v winds to true if proj_flag eq 2
         if (proj_flag .eq. 2) then
            rotate_winds_j: do j = 1, y
               rotate_winds_i: do i = 1, x
                  call uvgrid_to_uvtrue(u(i, j, k), v(i, j, k), &
                                        ut(i, j), vt(i, j), &
                                        lons(i, j))
               end do rotate_winds_i
            end do rotate_winds_j
         else
            ut(:, :) = u(:, :, k)
            vt(:, :) = v(:, :, k)
         end if

         write (output_unit, '(8f10.3)') ((ut(i, j), i=1, x), j=1, y)
         write (output_unit, '(8f10.3)') ((vt(i, j), i=1, x), j=1, y)
         write (output_unit, '(8f10.3)') ((t(i, j, k), i=1, x), j=1, y)
         write (output_unit, '(8f10.2)') ((ht(i, j, k), i=1, x), j=1, y)
         write (output_unit, '(8f10.5)') ((rh(i, j, k)*.01, i=1, x), j=1, y)

      end do level_loop
      deallocate (ut)
      deallocate (vt)

      ! write out the surface fields

      write (output_unit, '(8f10.3)') ((slp(i, j)*.01, i=1, x), j=1, y)
      write (output_unit, '(8f10.3)') ((psfc(i, j)*.01, i=1, x), j=1, y)
      write (output_unit, '(8f10.3)') ((t(i, j, z3 + 1), i=1, x), j=1, y)
      write (output_unit, '(8f10.4)') ((snocov(i, j), i=1, x), j=1, y)
      write (output_unit, '(8f10.4)') ((tskin(i, j), i=1, x), j=1, y)
      ! close the file

      close (output_unit)
      return
   end subroutine output_ralph2_format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module lapsprep_rams
