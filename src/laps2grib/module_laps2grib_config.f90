!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module:  laps2grib_config
!!

module laps2grib_config

   ! modules from the laps shared module library
   use map_utils

   implicit none

   type(proj_info)           :: laps_proj
   integer                    :: nx, ny, nz
   integer                    :: i4time
   integer                    :: num_grids
   character(len=9)           :: a9time
   real, allocatable          :: plevels_pa(:)
   character(len=512)         :: output_path, output_path2
   integer                    :: center_id, subcenter_id
   integer                    :: process_id, prod_status
   logical                    :: lrun_laps2grib

contains

   subroutine get_laps_proj(laps_data_root)

      implicit none

      character(len=*), intent(in) :: laps_data_root
      integer                 :: ncid, dimid, varid, nfstat, proj_code
      character(len=132)      :: grid_type
      character(len=256)      :: static_file
      logical                 :: file_exists
      real                    :: dx_m, la1, lo1, latin1, latin2, lov

      include 'netcdf.inc'

      static_file = trim(laps_data_root)//'/static/static.nest7grid'
      inquire (file=static_file, exist=file_exists)
      if (.not. file_exists) then
         print *, "-- static file not found:", trim(static_file)
         stop "no static_file"
      end if

      ! open the file
      nfstat = nf_open(static_file, nf_nowrite, ncid)
      if (nfstat .ne. nf_noerr) then
         print *, "-- error returned from nf_open: ", nfstat
         stop "netcdf file open error"
      end if

      ! get the data dimensions
      nx = 0
      ny = 0
      nfstat = nf_inq_dimid(ncid, 'x', dimid)
      nfstat = nf_inq_dimlen(ncid, dimid, nx)
      nfstat = nf_inq_dimid(ncid, 'y', dimid)
      nfstat = nf_inq_dimlen(ncid, dimid, ny)
      if ((nx .lt. 1) .or. (ny .lt. 1)) then
         print *, "-- bad nx or ny: ", nx, ny
         stop "bad dimensions"
      end if

      ! get the laps grid type
      nfstat = nf_inq_varid(ncid, 'grid_type', varid)
      nfstat = nf_get_var_text(ncid, varid, grid_type)
      if (grid_type(1:4) .eq. 'merc') then
         proj_code = proj_merc
      elseif ((grid_type(1:4) .eq. 'seca') .or. &
              (grid_type(1:4) .eq. 'tang')) then
         proj_code = proj_lc
      elseif (grid_type(1:4) .eq. 'pola') then
         proj_code = proj_ps
      elseif (grid_type(1:4) .eq. 'latl') then
         proj_code = proj_latlon
      else
         print *, "-- unknown laps grid_type: ", trim(grid_type)
         stop "bad grid_type"
      end if

      ! get the grid spacing
      nfstat = nf_inq_varid(ncid, 'dx', varid)
      nfstat = nf_get_var_real(ncid, varid, dx_m)

      ! get the truelat1/truelat2/lov
      nfstat = nf_inq_varid(ncid, 'latin1', varid)
      nfstat = nf_get_var_real(ncid, varid, latin1)
      nfstat = nf_inq_varid(ncid, 'latin2', varid)
      nfstat = nf_get_var_real(ncid, varid, latin2)
      nfstat = nf_inq_varid(ncid, 'lov', varid)
      nfstat = nf_get_var_real(ncid, varid, lov)

      ! get lower left corner lat/lon
      nfstat = nf_inq_varid(ncid, 'la1', varid)
      nfstat = nf_get_var_real(ncid, varid, la1)
      nfstat = nf_inq_varid(ncid, 'lo1', varid)
      nfstat = nf_get_var_real(ncid, varid, lo1)

      ! set up the laps_proj structure
      if (lo1 .gt. 180.) lo1 = lo1 - 360.
      if (lov .gt. 180.) lov = lov - 360.
      call map_set(proj_code, la1, lo1, dx_m, lov, latin1, latin2, nx, ny, &
                   laps_proj)

      ! close the file
      nfstat = nf_close(ncid)
      return
   end subroutine get_laps_proj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_laps_plevels(laps_data_root)

      implicit none

      character(len=*), intent(in)      :: laps_data_root
      character(len=256)                :: nl_file
      integer, parameter                :: max_lev = 100
      integer                           :: k
      logical                           :: file_exists
      real             :: pressures(max_lev)
      namelist /pressures_nl/ pressures

      ! look for and open the pressures.nl file
      nl_file = trim(laps_data_root)//'/static/pressures.nl'
      inquire (file=nl_file, exist=file_exists)
      if (.not. file_exists) then
         print *, "-- no pressures.nl file found: ", trim(nl_file)
         stop "no pressures.nl"
      end if

      pressures(:) = 0.0
      open (file=nl_file, unit=10, form='formatted', status='old')
      read (10, nml=pressures_nl)
      close (10)

      ! search through the levels in reverse order, because laps
      ! stores the data from the top down and the levels are specified
      ! in pressures_nl from the bottom up

      do k = max_lev, 1, -1

         if (pressures(k) .gt. 0.) then

            if (.not. allocated(plevels_pa)) then
               nz = k
               allocate (plevels_pa(nz))
            end if
            plevels_pa(nz - k + 1) = pressures(k)
         end if
      end do
      print *, "-- number of pressure levels found: ", nz
      print *, plevels_pa
      return
   end subroutine get_laps_plevels

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_laps_analtime

      implicit none
      integer   :: istatus

      call get_systime(i4time, a9time, istatus)

      if (istatus .ne. 1) then
         print *, "-- error getting laps analysis time!"
         stop "get_laps_analtime"
      end if
      print *, "-- laps analysis time: ", a9time
      return
   end subroutine get_laps_analtime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_laps_modeltime(modeltime, modeltime_passed_in)
      implicit none
      character(len=*), intent(in)   :: modeltime
      integer   :: modeltime_passed_in, istatus

      if (modeltime_passed_in .eq. 1) then
         print *, "generating model i4time using ", trim(modeltime)
         call i4time_fname_lp(modeltime, i4time, istatus)
         if (istatus .ne. 1) then
            print *, "-- error converting ", trim(modeltime), " to i4time!"
            stop "stop in get_laps_modeltime"
         end if
         a9time = trim(modeltime)
      else
         call get_modeltime(i4time, a9time, istatus)
         if (istatus .ne. 1) then
            print *, "-- error getting laps modeltime.dat!"
            stop "stop in get_laps_modeltime"
         end if
      end if

      print *, "-- laps model time: ", a9time
   end subroutine get_laps_modeltime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_laps2grib_nl(laps_data_root)

      implicit none
      character(len=*), intent(in)   :: laps_data_root
      character(len=256)             :: nl_file
      logical                        :: file_exists
      namelist /laps2grib_nl/ output_path, output_path2, num_grids, &
         center_id, subcenter_id, &
         process_id, prod_status, lrun_laps2grib

      ! set defaults
      num_grids = 1
      output_path(:) = " "
      output_path2(:) = " "
      center_id = 59 ! noaa esrl gsd (aka fsl), for use with templates
      subcenter_id = 1
      process_id = 100
      prod_status = 0
      lrun_laps2grib = .false.

      nl_file = trim(laps_data_root)//'/static/laps2grib.nl'
      inquire (file=nl_file, exist=file_exists)
      if (file_exists) then
         open (file=nl_file, unit=10, form='formatted', status='old')
         read (10, nml=laps2grib_nl)
         close (10)
      else
         print *, "-- no laps2grib.nl file found: ", trim(nl_file)
         print *, "-- will use default values!"
      end if
      if (len_trim(output_path) .lt. 1) then
         output_path = trim(laps_data_root)//'/lapsprd/gr2'
         output_path2 = trim(laps_data_root)//'/lapsprd/grb'
      end if

      print *, "-- other configuration: "
      print *, "   center id:       ", center_id
      print *, "   subcenter id:    ", subcenter_id
      print *, "   process id:      ", process_id
      print *, "   prod status:     ", prod_status
      print *, "   output path:     ", trim(output_path)
      return
   end subroutine read_laps2grib_nl

end module laps2grib_config
