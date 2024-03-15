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

module setup

   use mm5v3_io
   use wrfsi_static
   use wrf_netcdf
   use time_utils
   use map_utils
   ! contains routines for reading model setup

   implicit none

   integer, parameter            :: max_domains = 10
   ! file/path names
   character(len=1)              :: domain_num_str
   character(len=2)              :: domain_num_str2
   character(len=255)            :: laps_data_root
   character(len=255)            :: mm5_data_root
   character(len=255)            :: moad_dataroot
   character(len=255)            :: lfm_data_root
   character(len=255)            :: lfmprd_dir
   character(len=255)            :: data_file
   character(len=255)            :: terrain_file
   character(len=10)             :: mtype
   ! run configuration
   integer                       :: num_domains
   integer                       :: domain_num
   integer                       :: kprs
   integer                       :: domarg
   real, allocatable             :: prslvl(:)
   real                          :: redp_lvl
   logical                       :: keep_fdda
   logical                       :: split_output
   logical                       :: realtime
   logical                       :: make_donefile
   integer                       :: max_wait_sec
   logical                       :: proc_by_file_num
   integer                       :: start_file_num
   integer                       :: stop_file_num
   integer                       :: file_num_inc
   logical                       :: file_num3
   logical                       :: make_laps(max_domains)
   logical                       :: write_to_lapsdir(max_domains)
   logical                       :: make_v5d(max_domains)
   logical                       :: make_points(max_domains)
   integer                       :: v5d_compress
   character(len=32)             :: model_name
   logical                       :: do_smoothing
   logical                       :: use_model_pbl
   logical                       :: gribsfc(max_domains)
   logical                       :: gribua(max_domains)
   integer                       :: table_version
   integer                       :: center_id
   integer                       :: subcenter_id
   integer                       :: process_id(max_domains)

   ! time information
   integer                       :: num_times_avail
   integer                       :: num_times_to_proc
   integer                       :: start_time_index
   integer                       :: stop_time_index
   integer                       :: time_index_inc
   real                          :: output_freq_min
   character(len=24)             :: cycle_date
   character(len=24), allocatable :: times_to_proc(:)
   integer, allocatable           :: sim_tstep(:)
   character(len=24), parameter  :: static_date = '1900-01-01_00:00:00.0000'

   ! new stuff related to formatting of point forecast files

   character(len=24), allocatable :: point_times(:)
   character(len=3)             :: point_tz_label
   integer                       :: point_tz_utcoffset
   character(len=3)             :: point_windspd_units
   character(len=1)             :: point_temp_units
   character(len=5)             :: point_vent_units

   ! model domain configuration info
   integer                       :: nx
   integer                       :: ny
   integer                       :: ksigh
   integer                       :: ksigf
   type(proj_info)               :: proj
   real                          :: grid_spacing
   real, allocatable             :: terdot(:, :)
   real, allocatable             :: latdot(:, :)
   real, allocatable             :: londot(:, :)
   real, allocatable             :: mapfac_d(:, :)
   real, allocatable             :: coriolis(:, :)
   real                          :: proj_cent_lat
   real                          :: proj_cent_lon
   real                          :: truelat1
   real                          :: truelat2
   real                          :: cone_factor
   real                          :: pole_point
   real                          :: corner_lats(4)
   real                          :: corner_lons(4)
   real, allocatable             :: sigmah(:)
   real, allocatable             :: sigmaf(:)
   character(len=32)             :: projection
   logical                       :: clwflag
   logical                       :: iceflag
   logical                       :: graupelflag
   real                          :: ptop
   real                          :: pmslbase
   real                          :: tmslbase
   real                          :: dtdlnpbase
   real                          :: tisobase

   ! stuff for point output
   type point_struct
      character(len=10)           :: id
      real                        :: lat
      real                        :: lon
      real                        :: i
      real                        :: j
      integer                     :: elevation
      real                        :: hi_temp
      character(len=16)           :: hi_temp_time
      real                        :: lo_temp
      character(len=16)           :: lo_temp_time
      real                        :: avg_temp
      real                        :: avg_dewpt
      real                        :: total_pcp
      real                        :: total_snow
      integer                     :: output_unit
      character(len=80)           :: customer
   end type point_struct
   type(point_struct), allocatable :: point_rec(:)
   integer                        :: num_points

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine setup_lfmpost

      implicit none
      integer  :: lun_data, lun_terrain, status
      logical   :: file_ready
      character(len=2) :: domarg_str
      ! main lfmpost.nl namelist is in laps_data_root, so we must
      ! have one set

      call getenv("laps_data_root", laps_data_root)

      ! get the 1st argument which tells us which model type is being
      ! run.  (mtype)

      call getarg(1, mtype)
      if (mtype .eq. 'lfmpost.ex') then
         !must be an hp!
         call getarg(2, mtype)
         call getarg(3, domarg_str)
      else
         call getarg(2, domarg_str)
      end if

      read (domarg_str, '(i2)') domarg
      print *, 'domain number provided as argument: ', domarg
      if (domarg .gt. 0) then
         domain_num = domarg
      else
         domain_num = 1
      end if
      if ((mtype .eq. 'wrf       ') .or. (mtype .eq. 'wrf       ')) then
         mtype = 'wrf       '
      elseif ((mtype .eq. 'wrf2      ') .or. (mtype .eq. 'wrf2      ')) then
         mtype = 'wrf2      '
      elseif ((mtype .eq. 'mm5       ') .or. (mtype .eq. 'mm5       ')) then
         mtype = 'mm5       '
      else
         print *, 'unrecognized model type provided as first arg: ', mtype
         print *, 'wrf and mm5 are supported.'
         stop
      end if

      if (mtype(1:3) .eq. "mm5") then
         call getenv("mm5_data_root", mm5_data_root)
         if (mm5_data_root(1:3) .eq. "   ") then
            print *, 'setup: need to set mm5_data_root environment variable.'
            stop
         end if

         lfmprd_dir = trim(mm5_data_root)//"/mm5prd"
         lfm_data_root = mm5_data_root
      elseif (mtype(1:3) .eq. "wrf") then
         call getenv("moad_dataroot", moad_dataroot)
         if (moad_dataroot(1:3) .eq. "   ") then
            print *, 'setup: need to set moad_dataroot environment variable.'
            stop
         end if

         lfmprd_dir = trim(moad_dataroot)//"/wrfprd"
         lfm_data_root = moad_dataroot
      end if

      ! lets make sure we get it right
      print *, 'model prd directory:  ', trim(lfmprd_dir)

      call read_namelist

      write (domain_num_str, '(i1)') domain_num
      write (domain_num_str2, '(i2.2)') domain_num

      if (mtype .eq. 'mm5') then
         if (proc_by_file_num) then
            call make_data_file_name(mm5_data_root, domain_num_str, split_output, &
                                     start_file_num, data_file, file_num3)
         else
            call make_data_file_name(mm5_data_root, domain_num_str, split_output, &
                                     0, data_file, file_num3)
         end if
         print *, 'initial data file=', trim(data_file)
         call make_terrain_file_name(mm5_data_root, domain_num_str, terrain_file)
         print *, 'terrain file =', trim(terrain_file)
         call open_mm5v3(terrain_file, lun_terrain, status)
         if (split_output) then
            inquire (file=data_file, exist=file_ready)
            if (.not. file_ready) call io_wait(data_file, max_wait_sec)
         end if
         call open_mm5v3(data_file, lun_data, status)
         call time_setup(lun_data)
         call model_setup(lun_data, lun_terrain)
         close (lun_data)
         close (lun_terrain)

      elseif (mtype(1:3) .eq. 'wrf') then
         ! set up wrf model
         print *, 'setting up times for ', mtype
         if (mtype .eq. 'wrf') then
            call make_wrf_file_name(lfmprd_dir, domain_num, 0, data_file)
            call wrf_time_setup
         elseif (mtype .eq. 'wrf2') then
            call wrf2_time_setup
            call make_wrf2_file_name(lfmprd_dir, domain_num, &
                                     times_to_proc(1), data_file)
         end if
         print *, 'initial wrf file: ', trim(data_file)
         inquire (file=data_file, exist=file_ready)
         if (.not. file_ready) then
            call wrfio_wait(data_file, max_wait_sec)
         else
            if (realtime) call sleep(60)
         end if
         call open_wrfnc(data_file, lun_data, status)
         call model_setup(lun_data, lun_terrain)
         close (lun_data)

      end if
      ! if we want to make points, then set this up now
      if (make_points(domain_num)) then
         call init_points(status)
         if (status .ne. 0) then
            print *, 'point output requested, but cannot be fulfilled.'
            print *, 'ensure mm5_data_root/static/lfmpost_points.txt file exists.'
            make_points(domain_num) = .false.
         else
            print '(a,i3)', 'points initialized: ', num_points
         end if
      end if
      print *, 'ended setup_lfmpost'
   end subroutine setup_lfmpost
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_namelist

      implicit none
      integer, parameter          :: max_levels = 200
      integer                     :: k, unit, nml_unit
      integer                     :: status
      real                        :: levels_mb(max_levels)
      logical                     :: used
      character(len=255)          :: namelist_file
      character(len=32)           :: lfm_name(max_domains)

      namelist /lfmpost_nl/ num_domains, keep_fdda, split_output, levels_mb, &
         redp_lvl, lfm_name, proc_by_file_num, start_file_num, stop_file_num, &
         file_num_inc, file_num3, make_laps, realtime, write_to_lapsdir, &
         make_donefile, make_v5d, v5d_compress, max_wait_sec, do_smoothing, &
         gribsfc, gribua, table_version, center_id, subcenter_id, process_id, &
         make_points, point_tz_utcoffset, point_tz_label, point_windspd_units, &
         point_temp_units, point_vent_units, use_model_pbl

      if (lfmprd_dir(1:3) .eq. "   ") then
         namelist_file = "lfmpost.nl"
      else
         namelist_file = trim(lfmprd_dir)//'/../static/lfmpost.nl'
      end if

      ! get a unit number
      call get_file_unit(nml_unit)
      open (unit=nml_unit, file=trim(namelist_file), form='formatted', &
            status='old', iostat=status)
      if (status .ne. 0) then
         print *, 'error opening namelist file (', trim(namelist_file), '):', &
            status
         call abort
      end if

      ! initialize some values
      num_domains = 1
      model_name = '                                '
      kprs = 1
      levels_mb(:) = -1.
      redp_lvl = 0.
      keep_fdda = .true.
      split_output = .true.
      max_wait_sec = 3600
      proc_by_file_num = .false.
      start_file_num = 0
      stop_file_num = 999
      file_num_inc = 1
      file_num3 = .false.
      make_v5d(:) = .false.
      make_laps(:) = .true.
      make_points(:) = .false.
      v5d_compress = 2
      realtime = .true.
      make_donefile = .true.
      do_smoothing = .true.
      ! default grib settings
      gribsfc(:) = .false.
      gribua(:) = .false.
      table_version = 2
      center_id = 59   ! fsl
      subcenter_id = 2 ! lapb
      process_id(:) = 255
      write_to_lapsdir(:) = .false.
      point_tz_utcoffset = 0
      point_tz_label = 'utc'
      point_temp_units = 'f'
      point_windspd_units = 'kts'
      point_vent_units = 'm^2/s'
      use_model_pbl = .true.
      read (unit=nml_unit, nml=lfmpost_nl)
      close (nml_unit)

      ! see if we have to have laps_data_root
      if (make_laps(domain_num)) then
         if (laps_data_root(1:3) .eq. "   ") then
            print *, 'no laps_data_root environment variable set!'
            print *, 'either set laps_data_root or set make_laps = .false.'
            stop
         end if
      end if

      ! check domain number
      if (domain_num .gt. num_domains) then
         print *, 'you requested domain number: ', domain_num
         print *, 'but lfmpost.nl indicates only ', num_domains, &
            ' domains are configured.'
         stop
      end if

      ! count up number of levels requested.  they must be in monotonically
      ! decreasing order (bottom to top atmospherically)

      check_levels: do k = 2, max_levels
         if ((levels_mb(k) .gt. 0) .and. (levels_mb(k) .lt. levels_mb(k - 1))) then
            kprs = kprs + 1
         else
            exit check_levels
         end if
      end do check_levels

      if ((kprs .lt. 2) .and. (levels_mb(1) .lt. 0)) then
         print *, 'no valid pressure levels set in levels_mb!'
         call abort
      end if
      allocate (prslvl(kprs))
      prslvl(1:kprs) = levels_mb(1:kprs)*100.  ! convert to pascals!
      print *, 'setup information'
      print *, '-----------------'
      print '(a,i4,a)', 'interpolating to ', kprs, ' pressure levels.'
      do k = 1, kprs
         print '(a,f9.1,a)', 'level: ', prslvl(k), 'pa'
      end do
      model_name = lfm_name(domain_num)

      ! check point settings
      if (make_points(domain_num)) then
         ! timezone stuff
         if ((point_tz_utcoffset .lt. -12) .or. &
             (point_tz_utcoffset .gt. 12)) then
            print *, 'bad point_tz_utcoffset specified for points: ', &
               point_tz_utcoffset
            stop
         end if

         ! temp units
         if ((point_temp_units .eq. 'f') .or. (point_temp_units .eq. 'f')) then
            point_temp_units = 'f'
         elseif ((point_temp_units .eq. 'c') .or. (point_temp_units .eq. 'c')) then
            point_temp_units = 'c'
         elseif ((point_temp_units .eq. 'k') .or. (point_temp_units .eq. 'k')) then
            point_temp_units = 'k'
         else
            print *, 'bad point_temp_units:', point_temp_units
            stop
         end if

         ! wind units
         if ((point_windspd_units(1:1) .eq. 'k') .or. &
             (point_windspd_units(1:1) .eq. 'k')) then
            point_windspd_units = 'kts'
         elseif ((point_windspd_units(1:3) .eq. 'mph') .or. &
                 (point_windspd_units(1:3) .eq. 'mph')) then
            point_windspd_units = 'mph'
         elseif ((point_windspd_units(1:3) .eq. 'm/s') .or. &
                 (point_windspd_units(1:3) .eq. 'm/s')) then
            point_windspd_units = 'm/s'
         else
            print *, 'bad point_windspd_units:', point_windspd_units
            stop
         end if

         ! ventilation index units
         if ((point_vent_units(1:5) .eq. 'm^2/s') .or. &
             (point_vent_units(1:5) .eq. 'm^2/s')) then
            point_vent_units = 'm^2/s'
         elseif ((point_vent_units(1:5) .eq. 'kt-ft') .or. &
                 (point_vent_units(1:5) .eq. 'kt-ft')) then
            point_vent_units = 'kt-ft'
         else
            print *, 'bad point_vent_units:', point_windspd_units
            stop
         end if
      end if

   end subroutine read_namelist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine time_setup(lun_data)

      implicit none
      integer, intent(in)          :: lun_data
      integer                     :: min_to_add
      integer                     :: num_to_skip
      integer                     :: status
      integer                     :: sim_start_year
      integer                     :: sim_start_month
      integer                     :: sim_start_day
      integer                     :: sim_start_hour
      integer                     :: sim_start_min
      integer                     :: sim_start_sec
      integer                     :: sim_start_frac
      integer                     :: t
      real                        :: sim_stop_min
      real                        :: dom_start_min
      real                        :: dom_stop_min
      real                        :: fdda_start_min
      real                        :: fdda_stop_min
      real                        :: tapfrq
      real                        :: buffrq
      logical                     :: fdda_on

      character(len=24)           :: initial_date, new_date
      character(len=255)          :: wrfnl

      call get_mm5_time_info(lun_data, sim_start_year, sim_start_month, &
                             sim_start_day, sim_start_hour, &
                             sim_start_min, sim_start_sec, sim_start_frac, &
                             sim_stop_min, dom_start_min, dom_stop_min, &
                             fdda_on, fdda_start_min, fdda_stop_min, &
                             tapfrq, buffrq, status)

      ! compute total number of times available in this domain
      num_times_avail = nint(dom_stop_min - dom_start_min)/nint(tapfrq) + 1
      print *, 'num_times_avail = ', num_times_avail
      write (initial_date, &
             '(i4.4,"-",i2.2,"-",i2.2,"_",i2.2,":",i2.2,":",i2.2,".",i4.4)') &
         sim_start_year, sim_start_month, sim_start_day, sim_start_hour, &
         sim_start_min, sim_start_sec, sim_start_frac
      print '(2a)', 'simulation start time: ', initial_date
      cycle_date = initial_date
      if ((fdda_on) .and. (.not. keep_fdda)) then
         min_to_add = fdda_stop_min - dom_start_min
         call geth_newdate(new_date, cycle_date(1:16), min_to_add)
         cycle_date = new_date
      end if
      if (dom_start_min .gt. 0.) then
         call geth_newdate(new_date, initial_date(1:16), nint(dom_start_min))
         cycle_date = new_date
      end if
      print '(2a)', 'model cycle time (excluding fdda): ', cycle_date
      print '(2a)', 'domain initial time: ', initial_date
      if ((.not. proc_by_file_num) .and. (split_output)) then
         if ((fdda_on) .and. (.not. keep_fdda) .and. &
             (dom_start_min .lt. fdda_stop_min)) then
            num_to_skip = nint(fdda_stop_min - dom_start_min)/nint(tapfrq)
            num_times_to_proc = num_times_avail - num_to_skip
            start_time_index = num_to_skip + 1
         else
            num_to_skip = 0
            num_times_to_proc = num_times_avail
            start_time_index = 1
         end if
         stop_time_index = num_times_avail
         time_index_inc = 1
      else
         num_to_skip = start_file_num
         stop_file_num = min(stop_file_num, num_times_avail - 1)
         stop_time_index = stop_file_num + 1
         num_times_to_proc = stop_file_num - start_file_num + 1
         start_time_index = start_file_num + 1
         time_index_inc = file_num_inc
      end if
      allocate (times_to_proc(num_times_to_proc))
      print *, 'total number of output times: ', num_times_to_proc
      do t = 1, num_times_avail, time_index_inc
         min_to_add = (t - 1)*nint(tapfrq)
         if ((proc_by_file_num) .and. (split_output)) min_to_add = (t - 1)*nint(buffrq)
         call geth_newdate(new_date, initial_date(1:16), min_to_add)
         if ((t .ge. start_time_index) .and. (t .le. stop_time_index)) then
            print *, 'will process date: ', new_date
            times_to_proc(t - num_to_skip) = new_date
         end if
      end do
      print *, ' '

   end subroutine time_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine model_setup(lun_data, lun_terrain)

      implicit none
      integer              :: k, km1
      integer, intent(in)  :: lun_data
      integer, intent(in)  :: lun_terrain
      integer              :: status
      real, allocatable    :: tempsigma(:)
      real                 :: lat1, lon1, stdlon, dx_m, dy_m
      character(len=255)   :: wrfinitfile
      integer              :: luninit
      logical              :: fileready

      if (mtype .eq. 'mm5') then
         call get_mm5_map(lun_terrain, projection, proj_cent_lat, proj_cent_lon, &
                          truelat1, truelat2, cone_factor, pole_point, &
                          grid_spacing, nx, ny, status)
         allocate (latdot(nx, ny))
         allocate (londot(nx, ny))
         allocate (terdot(nx, ny))
         call get_mm5_2d(lun_terrain, 'latitdot ', static_date, latdot, 'd   ', status)
         call get_mm5_2d(lun_terrain, 'longidot ', static_date, londot, 'd   ', status)
         call get_mm5_2d(lun_terrain, 'terrain  ', static_date, terdot, 'd   ', status)
         lat1 = latdot(1, 1)
         lon1 = londot(1, 1)
         stdlon = proj_cent_lon
      elseif (mtype(1:3) .eq. 'wrf') then

         ! use routines in module wrfsi_static to read from static.wrfsi
         call get_wrfsi_static_dims(moad_dataroot, domain_num, nx, ny)
         call get_wrfsi_static_proj(moad_dataroot, domain_num, &
                                    projection, lat1, lon1, dx_m, dy_m, stdlon, truelat1, truelat2)
         if (dx_m .ne. dy_m) then
            print *, 'wrf dx != dy...not ready for this!', dx_m, dy_m
            stop
         else
            grid_spacing = dx_m
         end if
         allocate (latdot(nx, ny))
         allocate (londot(nx, ny))
         allocate (terdot(nx, ny))

         ! get the lat/lon for the non-staggered grid
         call get_wrfsi_static_latlon(moad_dataroot, domain_num, 'n', &
                                      latdot, londot)

         ! get the terrain height
!     call get_wrfsi_static_2d(moad_dataroot, 'avg', terdot)
         wrfinitfile = trim(moad_dataroot)//'/wrfprd/wrfinput_d'// &
                       domain_num_str2
         inquire (file=wrfinitfile, exist=fileready)
         if (.not. fileready) then
            call wrfio_wait(wrfinitfile, max_wait_sec)
         else
            if (realtime) call sleep(60)
         end if
         call open_wrfnc(wrfinitfile, luninit, status)
         if (status .ne. 0) then
            print *, 'wrf input file not found: ', trim(wrfinitfile)
            print *, 'needed for terrain heights!'
            stop
         end if
         call get_wrfnc_2d(luninit, "hgt", "a", nx, ny, 1, terdot, status)
         call close_wrfnc(luninit)
      end if

      ! use the map_set routine to set up the projection information structure
      select case (projection)
      case ('lambert conformal               ')
         call map_set(proj_lc, lat1, lon1, grid_spacing, &
                      stdlon, truelat1, truelat2, nx, ny, proj)
      case ('polar stereographic             ')
         call map_set(proj_ps, lat1, lon1, grid_spacing, &
                      stdlon, truelat1, truelat2, nx, ny, proj)
      case ('mercator                        ')
         call map_set(proj_merc, lat1, lon1, grid_spacing, &
                      stdlon, truelat1, truelat2, nx, ny, proj)
      end select

      if (make_v5d(domain_num)) then
         allocate (mapfac_d(nx, ny))
         allocate (coriolis(nx, ny))
         if (mtype .eq. 'mm5') then
            call get_mm5_2d(lun_terrain, 'mapfacdt ', static_date, &
                            mapfac_d, 'd   ', status)
            call get_mm5_2d(lun_terrain, 'coriolis ', static_date, &
                            coriolis, 'd   ', status)
         elseif (mtype(1:3) .eq. 'wrf') then
            call get_wrfsi_static_2d(moad_dataroot, domain_num, 'mfl', mapfac_d)
            call get_wrfsi_static_2d(moad_dataroot, domain_num, 'cph', coriolis)
         end if

      end if
      print *, ' '
      print '(a,2f10.2)', 'min/max value of terrain: ', minval(terdot), &
         maxval(terdot)
      print *, ' '
      corner_lats(1) = latdot(1, 1)
      corner_lons(1) = londot(1, 1)
      corner_lats(2) = latdot(1, ny)
      corner_lons(2) = londot(1, ny)
      corner_lats(3) = latdot(nx, ny)
      corner_lons(3) = londot(nx, ny)
      corner_lats(4) = latdot(nx, 1)
      corner_lons(4) = londot(nx, 1)
      print *, ' '
      print *, 'corner points from dot point lat/lon arrays:'
      print *, '============================================'
      print *, ' '
      print '(f8.3,1x,f8.3,10x,f8.3,1x,f8.3)', &
         corner_lats(2), corner_lons(2), corner_lats(3), corner_lons(3)
      print *, '      (nw)----------------------(ne)'
      print *, '        |                        |'
      print *, '        |                        |'
      print *, '      (sw)----------------------(se)'
      print '(f8.3,1x,f8.3,10x,f8.3,1x,f8.3)', &
         corner_lats(1), corner_lons(1), corner_lats(4), corner_lons(4)
      print *, ' '
      !deallocate(latdot)
      !deallocate(londot)

      if (mtype .eq. 'mm5') then
         call get_mm5_misc(lun_data, ksigh, ptop, pmslbase, tmslbase, &
                           dtdlnpbase, tisobase, &
                           clwflag, iceflag, graupelflag)
         ksigf = ksigh + 1
         allocate (tempsigma(ksigh))
         allocate (sigmah(ksigh))
         allocate (sigmaf(ksigf))
         call get_mm5_1d(lun_data, 'sigmah   ', static_date, tempsigma, status)
         if (status .ne. 0) then
            print *, 'sigmah not in output file!'
            stop
         end if
         ! re-order sigma levels to be from ground up and compute sigmaf
         do k = 1, ksigh
            sigmah(k) = tempsigma(ksigh + 1 - k)
         end do
         sigmaf(1) = 1.0
         sigmaf(ksigf) = 0.0
         do k = 2, ksigh
            km1 = k - 1
            sigmaf(k) = sigmah(km1) - (sigmaf(km1) - sigmah(km1))
         end do
         deallocate (tempsigma)
         print '(a,i4)', 'number of half sigma levels: ', ksigh
         print '(a)', 'level    sigfull   sighalf'
         do k = 1, ksigh
            print '(i4,4x,f8.6,2x,f8.6)', k, sigmaf(k), sigmah(k)
         end do
         print '(i4,4x,f8.6)', ksigf, sigmaf(ksigf)
         print *, ' '
      elseif (mtype(1:3) .eq. 'wrf') then
         if (mtype .eq. 'wrf') then
            call get_wrf_misc(lun_data, ksigh, ksigf, ptop, clwflag, iceflag, &
                              graupelflag)
         elseif (mtype .eq. 'wrf2') then
            call get_wrf2_misc(lun_data, ksigh, ksigf, ptop, clwflag, iceflag, &
                               graupelflag)
         end if
         allocate (sigmah(ksigh))
         allocate (sigmaf(ksigf))
         call get_wrf_1d(lun_data, 'znu', sigmah, status)
         call get_wrf_1d(lun_data, 'znw', sigmaf, status)

         ! some things that wrf does not provide but are initialized
         ! here in a hardcoded way to maintain compatibility with
         ! original mm5post core code.
         pmslbase = 100000.  !not currently used
         tmslbase = 275.     !not currently used
         dtdlnpbase = 50.    ! only used for downward t extrapolation
         tisobase = 0.      !not currently used

      end if

      print '(a,i4)', 'number of half sigma levels: ', ksigh
      print '(a)', 'level    sigfull   sighalf'
      do k = 1, ksigh
         print '(i4,4x,f8.6,2x,f8.6)', k, sigmaf(k), sigmah(k)
      end do
      print '(i4,4x,f8.6)', ksigf, sigmaf(ksigf)
      print *, ' '

      return
   end subroutine model_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_points(status)

      ! reads the lfmpost_points.txt file, checks validity of each requested
      ! point, initializes their output files, etc.  returns non-zero status
      ! if point initialization cannot be completed.
      use map_utils
      implicit none
      integer, intent(out)        :: status
      integer, parameter       :: max_points = 1000
      type(point_struct), allocatable :: points_temp(:)
      type(point_struct)       :: point
      character(len=200)       :: pointfile
      integer                  :: pointunit
      logical                  :: pointfileexists
      logical                  :: lunused
      integer                  :: ic, jc
      integer                  :: year, year2, month, day, hour, minute, sec, jjj
      character(len=9)        :: cyclestr
      character(len=200)       :: outfile
      integer                  :: outunit
      character(len=10)       :: id
      real                     :: lat
      real                     :: lon
      character(len=80)        :: customer
      real                     :: elevation
      real, external           :: bint
      character(len=3)         :: domnum_str
      integer                  :: nestpt, i
      character(len=13)        :: odate, ndate
      character(len=24)        :: tempdate
      character(len=3)         :: model_title
      status = 0
      num_points = 0
      pointfile = trim(lfmprd_dir)//'/../static/lfmpost_points.txt'
      inquire (file=pointfile, exist=pointfileexists)
      if (.not. pointfileexists) then
         print *, 'point file missing: ', pointfile
         status = 1
      else
         ! get a logical unit number to use
         call get_file_unit(pointunit)
         open (file=pointfile, unit=pointunit, status='old', form='formatted', &
               access='sequential')
         allocate (points_temp(max_points))
         do while (num_points .lt. max_points)
            read (pointunit, '(i2,1x,a10,1x,f8.4,1x,f9.4,1x,a)', end=99) &
               nestpt, id, lat, lon, customer

            if (nestpt .eq. domain_num) then
               point%id = id
               point%lat = lat
               point%lon = lon
               point%customer = customer
               ! compute the i/j  for this point

               call latlon_to_ij(proj, point%lat, point%lon, point%i, point%j)
               if (abs(point%i - 1.) .lt. .001) point%i = 1.
               if (abs(point%i - nx) .lt. .001) point%i = nx
               if (abs(point%j - 1.) .lt. .001) point%j = 1.
               if (abs(point%j - ny) .lt. .001) point%j = ny

               if ((point%i .ge. 1.) .and. (point%i .le. nx) .and. &
                   (point%j .ge. 1.) .and. (point%j .le. ny)) then
                  print *, 'initializing point location ', point%id
                  ic = nint(point%i)
                  jc = nint(point%j)

                  elevation = bint(point%i, point%j, terdot, nx, ny)
                  point%elevation = nint(elevation*3.2808)
                  point%hi_temp = -999.9
                  point%lo_temp = 999.9
                  point%hi_temp_time = '00/00/0000 00:00'
                  point%lo_temp_time = '00/00/0000 00:00'
                  point%avg_temp = 0.
                  point%avg_dewpt = 0.
                  point%total_pcp = 0.
                  point%total_snow = 0.
                  ! create output file name, find a logical unit number, open the
                  ! file, and write out the header
                  call split_date_char(times_to_proc(1), year, &
                                       month, day, hour, minute, sec)
                  jjj = compute_day_of_year(year, month, day)
                  year2 = mod(year, 1000)
                  write (cyclestr, '(i2.2,i3.3,i2.2,i2.2)') year2, jjj, hour, minute
                  write (domnum_str, '("d",i2.2)') domain_num
                  outfile = trim(lfmprd_dir)//'/'//domnum_str// &
                            '/points/'//trim(point%id)//'_'//cyclestr// &
                            '_fcst.txt'
                  call get_file_unit(outunit)
                  point%output_unit = outunit
                  num_points = num_points + 1
                  points_temp(num_points) = point
                  open (file=outfile, unit=outunit, form='formatted', &
                        access='sequential')
                  write (outunit, '("****************************************", &
                                 &"******************************************************")')
                  write (outunit, &
                    '("location: ",a,2x,"lat: ",f8.4,2x,"lon: ",f9.4,2x, &
                   &"i: ",f7.2,2x,"j: ",f7.2)') point%id, point%lat, &
                         point%lon, point%i, point%j
                  if (mtype(1:3) .eq. 'mm5') then
                     model_title = 'mm5'
                  elseif (mtype(1:3) .eq. 'wrf') then
                     model_title = 'wrf'
                  else
                     model_title = 'unk'
                  end if
                  write (outunit, &
                     '(a3,1x,f4.1," km    forecast cycle: ",a,2x,"dom: ",&
                     & i2,2x,"model elevation: ",i4)') &
                       model_title, grid_spacing/1000., cyclestr, domain_num, point%elevation
                  write (outunit, '("****************************************", &
                                 &"******************************************************")')
                  write (outunit, &
                         '("date       time  tmp dpt rh  wind   cei vis  weather  precp snow vent   mixht pblwnd hm hh fbg")')
                  write (outunit, &
                         '(a3,8x,a3,3x,a1,3x,a1,3x,"%",3x,"dg@",a3,1x,"hft mile",10x,"in",4x,"in",3x,a5,2x,"ftagl",1x,"dg@",a3)') &
                     point_tz_label, point_tz_label, point_temp_units, &
                     point_temp_units, point_windspd_units, point_vent_units, point_windspd_units

                  write (outunit, &
                         '("---------- ----- --- --- --- ------ --- ---- -------- ----- ---- ------ ----- ------ -- -- ---")')

               else
                  print *, 'point location ', trim(point%id), ' outside of domain!'
               end if
            end if
         end do
99       if (num_points .eq. 0.) then
            status = 1
            print *, 'no valid points found to process'
         else
            allocate (point_rec(num_points))
            point_rec = points_temp(1:num_points)
            deallocate (points_temp)
            status = 0
         end if
      end if

      ! set up point times time zone offset
      allocate (point_times(num_times_to_proc))
      if (point_tz_utcoffset .eq. 0) then
         point_times = times_to_proc
      else
         do i = 1, num_times_to_proc
            tempdate = times_to_proc(i)
            odate = tempdate(1:13)
            call geth_newdate(ndate, odate, point_tz_utcoffset)
            tempdate(1:13) = ndate
            point_times(i) = tempdate
         end do
      end if
      return
   end subroutine init_points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine wrf_time_setup

      ! sets up the list of expected times and files to find/process for
      ! a run of the wrf model by reading the namelist.input file.

      implicit none

      ! here are the variables declared for being able to read the wrf 1.2
      ! namelist sections required
      integer  :: time_step_max, max_dom, dyn_opt, rk_ord, diff_opt
      integer  :: km_opt, damp_opt, isfflx, ifsnow, icloud, num_soil_layers
      integer  :: spec_bdy_width, spec_zone, relax_zone, tile_sz_x, tile_sz_y
      integer  :: numtiles, debug_level
      namelist /namelist_01/ time_step_max, max_dom, dyn_opt, rk_ord, diff_opt, &
         km_opt, damp_opt, isfflx, ifsnow, icloud, num_soil_layers, &
         spec_bdy_width, spec_zone, relax_zone, tile_sz_x, tile_sz_y, &
         numtiles, debug_level

      integer  :: grid_id, level, s_we, e_we, s_sn, e_sn, s_vert, e_vert
      integer  :: time_step_count_output, frames_per_outfile
      integer  :: time_step_count_restart, time_step_begin_restart
      integer  :: time_step_sound
      namelist /namelist_02/ grid_id, level, s_we, e_we, s_sn, e_sn, s_vert, e_vert, &
         time_step_count_output, frames_per_outfile, &
         time_step_count_restart, time_step_begin_restart, &
         time_step_sound

      real     :: dx, dy, dt, ztop, zdamp, dampcoef
      logical  :: non_hydrostatic
      real     :: smdiv, emdiv, epssm, khdif, kvdif, mix_cr_len
      real     :: radt, bldt, cudt, gmt
      integer  :: julyr, julday
      namelist /namelist_03/ dx, dy, dt, ztop, zdamp, dampcoef, &
         non_hydrostatic, &
         smdiv, emdiv, epssm, khdif, kvdif, mix_cr_len, &
         radt, bldt, cudt, gmt, julyr, julday

      integer  :: start_year, start_month, start_day, start_hour, start_second
      integer  :: start_minute, end_year, end_month, end_day, end_hour, &
                  end_minute, end_second
      integer  :: interval_seconds, real_data_init_type
      namelist /namelist_05/ start_year, start_month, start_day, start_hour, &
         start_minute, start_second, end_year, end_month, &
         end_day, end_hour, &
         end_minute, end_second, &
         interval_seconds, real_data_init_type

      character(len=255)  :: wrfnl
      logical             :: found_wrfnl, used
      integer             :: nml_unit, unit
      integer             :: t, sec_to_add, status
      character(len=24)   :: new_date

      found_wrfnl = .false.
      ! first, look for "namelist.input" in lfmprd_dir
      wrfnl = trim(lfmprd_dir)//'/namelist.input'
      inquire (file=wrfnl, exist=found_wrfnl)

      ! if not found, look in lfmprd_dir/../static for wrf.nl
      if (.not. found_wrfnl) then
         wrfnl = trim(lfmprd_dir)//'/../static/wrf.nl'
         inquire (file=wrfnl, exist=found_wrfnl)
      end if

      if (.not. found_wrfnl) then
         wrfnl = 'namelist.input'
         inquire (file=wrfnl, exist=found_wrfnl)
      end if
      if (.not. found_wrfnl) then
         print *, 'unable to find valid wrf namelist.input'
         stop
      end if

      ! presumably, we have a valid file, so open and read the
      ! four appropriate namelist sections.

      ! get a unit number
      call get_file_unit(nml_unit)
      open (unit=nml_unit, file=trim(wrfnl), form='formatted', &
            status='old', iostat=status)
      if (status .ne. 0) then
         print *, 'error opening namelist file (', trim(wrfnl), '):', &
            status
         call abort
      end if
      rewind (nml_unit)
      read (nml_unit, nml=namelist_01)
      read (nml_unit, nml=namelist_02)
      read (nml_unit, nml=namelist_03)
      read (nml_unit, nml=namelist_05)
      close (nml_unit)

      num_times_avail = time_step_max/time_step_count_output + 1

      ! for now, hardcode to process all available times
      num_times_to_proc = num_times_avail
      start_time_index = 1
      stop_time_index = num_times_avail
      time_index_inc = 1

      ! fill cycle date
      write (cycle_date, &
             '(i4.4,"-",i2.2,"-",i2.2,"_",i2.2,":",i2.2,":",i2.2,".0000")') &
         start_year, start_month, start_day, start_hour, &
         start_minute, start_second

      allocate (times_to_proc(num_times_to_proc))
      allocate (sim_tstep(num_times_to_proc))
      print *, 'total number of output times: ', num_times_to_proc
      do t = 1, num_times_avail, time_index_inc
         sec_to_add = (t - 1)*nint(dt)*time_step_count_output
         sim_tstep(t) = nint(sec_to_add/dt)
         call geth_newdate(new_date, cycle_date(1:19), sec_to_add)
         print *, 'will process date: ', new_date, sim_tstep(t)
         times_to_proc(t) = new_date
      end do
      print *, ' '
      output_freq_min = dt*float(time_step_count_output)/60.

   end subroutine wrf_time_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine wrf2_time_setup

      ! this subroutine sets up list of times to process for wrfv2.  it
      ! requires a special "wrftimes.nl" file to be in the
      ! moad_dataroot/wrfprd area.  this namelist is a subset of the
      ! time control items from the standard wrfv2 namelist.input file

      implicit none

      ! namelist entries
      integer, dimension(max_domains) :: start_year, end_year
      integer, dimension(max_domains) :: start_month, end_month
      integer, dimension(max_domains) :: start_day, end_day
      integer, dimension(max_domains) :: start_hour, end_hour
      integer, dimension(max_domains) :: start_minute, end_minute
      integer, dimension(max_domains) :: start_second, end_second
      integer, dimension(max_domains) :: history_interval
      integer, dimension(max_domains) :: frames_per_outfile

      ! others
      integer, parameter  :: max_times = 500
      character(len=24), allocatable  :: possible_times(:)
      character(len=24) :: hdate, hdate_beg, hdate_end, hdate_new
      character(len=255)              :: wrftimesnl
      integer  :: first_time_index, unit, status

      namelist /wrftimes/ &
         start_year, start_month, start_day, &
         start_hour, start_minute, start_second, &
         end_year, end_month, end_day, end_hour, end_minute, &
         end_second, history_interval, frames_per_outfile

      ! initialization of some namelist variables that might
      ! not be present.
      start_minute(:) = 0
      start_second(:) = 0
      end_minute(:) = 0
      end_second(:) = 0
      history_interval(:) = -1
      frames_per_outfile(:) = 1

      wrftimesnl = trim(moad_dataroot)//'/wrfprd/wrftimes.nl'

      ! get a unit number
      call get_file_unit(unit)
      open (unit=unit, file=trim(wrftimesnl), form='formatted', &
            status='old', iostat=status)
      if (status .ne. 0) then
         print *, 'error opening namelist file (', trim(wrftimesnl), '):', &
            status
         call abort
      end if
      read (unit, nml=wrftimes)

      if (frames_per_outfile(domain_num) .ne. 1) then
         print *, 'frames_per_outfile must = 1 for lfmpost!'
         stop
      end if

      allocate (possible_times(max_times))
      num_times_to_proc = 1
      write (hdate, 50) &
         start_year(domain_num), start_month(domain_num), start_day(domain_num), &
         start_hour(domain_num), start_minute(domain_num), &
         start_second(domain_num)
50    format(i4.4, "-", i2.2, "-", i2.2, "_", i2.2, ":", i2.2, ":", i2.2, ".0000")
      possible_times(1) = hdate
      write (hdate_end, 50) end_year(domain_num), end_month(domain_num), &
         end_day(domain_num), end_hour(domain_num), end_minute(domain_num), &
         end_second(domain_num)
      do while (hdate .le. hdate_end)
         call geth_newdate(hdate, hdate(1:16), history_interval(domain_num))
         if (hdate .le. hdate_end) then
            num_times_to_proc = num_times_to_proc + 1
            if (num_times_to_proc .gt. max_times) then
               print *, 'exceeded maximum input times!'
               stop
            end if
            possible_times(num_times_to_proc) = hdate
         end if
      end do
      allocate (times_to_proc(num_times_to_proc))
      times_to_proc(:) = possible_times(1:num_times_to_proc)
      deallocate (possible_times)

      ! for now, hardcode to process all available times
      num_times_avail = num_times_to_proc
      start_time_index = 1
      stop_time_index = num_times_to_proc
      time_index_inc = 1
      cycle_date = times_to_proc(1)
      return
   end subroutine wrf2_time_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module setup
