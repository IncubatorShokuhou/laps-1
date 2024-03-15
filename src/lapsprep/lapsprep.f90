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

program lapsprep
   !
   ! purpose
   ! =======
   ! prepares laps analysis data for ingest by various nwp model pre-processors.
   ! currently supports mm5v3 (outputs pregrid v3 format), wrf (outputs
   ! gribprep format), and rams 4.3 (outputs ralph2 format).
   !
   ! arguments
   ! =========
   !  laps_valid_time   - optional command line argument of form yyjjjhhmm
   !                      specifying time for which to build output
   !                      if not present, the program will use the latest
   !                      available analysis based on laps systime.dat
   !
   ! remarks
   ! =======
   !  1. you must set the laps_data_root environment variable before running.
   !  2. other program controls in lapsprep.nl
   !
   ! history
   ! =======
   ! 28 nov 2000 -- original -- brent shaw
   !    (based on lapsreader program originally developed by dave gill of
   !     ncar to support mm5)

   ! module declarations

   use lapsprep_constants
   use setup
   use laps_static
   use lapsprep_mm5
   use lapsprep_wrf
   use lapsprep_rams
   use lapsprep_netcdf
   use mem_namelist, only: read_namelist_laps, model_cycle_time, r_missing_data

   ! variable declarations

   implicit none

   ! declarations for use of netcdf library

   include "netcdf.inc"
   integer :: cdfid, rcode
   integer :: zid
   integer :: z
   integer, dimension(4) :: start, count
   integer, dimension(2) :: startc, countc
   integer :: vid
   character(len=132) :: dum

   ! arrays for data
   real, allocatable, dimension(:, :, :) :: u, v, w, t, rh, ht, &
                                            lwc, rai, sno, pic, ice, sh, mr, vv, &
                                            virtual_t, rho, lcp, mvd, &
                                            soilt, soilm
   real, allocatable, dimension(:, :)   :: slp, psfc, snocov, d2d, tskin, &
                                           tpw_before, tpw_after
   real, allocatable, dimension(:)     :: p
   real, parameter                       :: tiny = 1.0e-30

   ! miscellaneous local variables

   integer :: out_loop, loop, var_loop, i, j, k, kbot, istatus, len_dir
   integer :: moment_uphys
   logical :: file_present, in_prebal_list
   real    :: rhmod, lwcmod, shmod, icemod
   real    :: rhadj, sh_sfc
   real    :: lwc_limit
   real    :: hydrometeor_scale_pcp, hydrometeor_scale_cld
   character*200 :: static_dir, filename

   ! some stuff for jax to handle lga problem
   ! with constant mr above 300 mb
   logical :: jaxsbn
   real  :: weight_top, weight_bot, newsh
   real, external ::make_rh
   integer :: k300

   jaxsbn = .false.

   ! beginning of code

   call get_systime(i4time, laps_file_time, istatus)

   if (istatus .ne. 1) then
      write (*, *) 'bad istatus from get_systime'
      stop
   end if

   print *, 'laps_file_time = ', laps_file_time

!   read global parameters into module memory structure
   call get_directory('static', static_dir, len_dir)
   filename = static_dir(1:len_dir)//'/nest7grid.parms'
   call read_namelist_laps('lapsparms', filename)

   read (laps_file_time, '(i2.2,i3.3,i2.2,i2.2)') valid_yyyy, valid_jjj, &
      valid_hh, valid_min
   if (valid_yyyy .lt. 80) then
      valid_yyyy = 2000 + valid_yyyy
   else
      valid_yyyy = 1900 + valid_yyyy
   end if

   print '(2a)', 'running lapsprep using a9_time of: ', laps_file_time

   ! get the laps_data_root from the environment.

   call getenv('laps_data_root', laps_data_root)
   print '(2a)', 'laps_data_root=', laps_data_root

   !  get the namelist items (from the setup module).

   call read_namelist

   if ((model_cycle_time .ne. 0) .and. (lapsprep_always_write .eqv. .false.)) then
      print *, 'model_cycle_time = ', model_cycle_time
      if (i4time .ne. ((i4time/model_cycle_time)*model_cycle_time) .or. model_cycle_time .lt. 0) then
         write (6, *) ' not on the model run cycle: stopping...'
         stop
      end if
   end if

   moment_uphys = 1

   ! get the static information (projection, dimensions,etc.)

   print '(a)', 'getting horizontal grid specs from static file.'
   call get_horiz_grid_spec(laps_data_root)

   ! now that we have laps grid info, set up the hydrometeor scaling
   ! factor, which scales the concentrations of hydormeteors for this
   ! grid spacing.  we assume the values from laps are approprate on
   ! a grid with radar scaling (approx. 2km)

!beka    hydrometeor_scale = 2./dx  ! dx is in km
   if (hydrometeor_scale_factor_pcp .ge. 0.) then
      hydrometeor_scale_pcp = hydrometeor_scale_factor_pcp
   else
      hydrometeor_scale_pcp = -hydrometeor_scale_factor_pcp/dx  ! dx is in km
   end if

   if (hydrometeor_scale_factor_cld .ge. 0.) then
      hydrometeor_scale_cld = hydrometeor_scale_factor_cld
   else
      hydrometeor_scale_cld = -hydrometeor_scale_factor_cld/dx  ! dx is in km
   end if

   print *, 'beka', hydrometeor_scale_pcp, hydrometeor_scale_factor_pcp
   print *, 'cj', hydrometeor_scale_cld, hydrometeor_scale_factor_cld

   !  loop through each of the requested extensions for this date.  each of the
   !  extensions has a couple of the variables that we want.

   print '(a)', 'starting loop for each laps file'
   file_loop: do loop = 1, num_ext

      print *, 'looking for ', ext(loop)
      !  if this is a microphysical species but not doing
      !  a hotstart, then cycle over this file.

      if (((trim(ext(loop)) .eq. 'lwc') .or. (trim(ext(loop)) .eq. 'lcp')) .and. &
          (.not. hotstart)) then
         cycle file_loop
      end if

      !  determine based on extension whether to use prebalanced directories

      if (use_sfc_bal) then
         if ((trim(ext(loop)) .ne. 'lw3') .and. & ! list of balance files
             (trim(ext(loop)) .ne. 'lt1') .and. &
             (trim(ext(loop)) .ne. 'lq3') .and. &
             (trim(ext(loop)) .ne. 'lsx') .and. &
             (trim(ext(loop)) .ne. 'lh3')) then
            in_prebal_list = .true.
         else
            in_prebal_list = .false.
         end if
      else ! false
         if ((trim(ext(loop)) .ne. 'lw3') .and. & ! list of balance files
             (trim(ext(loop)) .ne. 'lt1') .and. &
             (trim(ext(loop)) .ne. 'lq3') .and. &
             !             (trim(ext(loop)) .ne. 'lsx' ).and. &
             (trim(ext(loop)) .ne. 'lh3')) then
            in_prebal_list = .true.
         else
            in_prebal_list = .false.
         end if
      end if

      !  build the input file name.   the input file.
      if (in_prebal_list) then
         input_laps_file = trim(laps_data_root)//'/lapsprd/'// &
                           trim(ext(loop))//'/'//laps_file_time//'.'// &
                           trim(ext(loop))
      else
         if (balance) then
            input_laps_file = trim(laps_data_root)//'/lapsprd/balance/'// &
                              trim(ext(loop))//'/'//laps_file_time//'.'// &
                              trim(ext(loop))
         else
            input_laps_file = trim(laps_data_root)//'/lapsprd/'// &
                              trim(ext(loop))//'/'//laps_file_time//'.'// &
                              trim(ext(loop))
         end if
      end if
      print *, 'opening: ', input_laps_file

      ! determine if the file exists

      inquire (file=trim(input_laps_file), exist=file_present)
      if (.not. file_present) then
         if ((ext(loop) .eq. 'lt1') .or. &
             (ext(loop) .eq. 'lw3') .or. &
             (ext(loop) .eq. 'lh3') .or. &
             (ext(loop) .eq. 'lsx')) then
            print '(a)', 'mandatory file not available:', input_laps_file
            stop 'not_enough_data'
         else if ((ext(loop) .eq. 'lq3') .or. &
                  (ext(loop) .eq. 'lwc')) then
            print '(a)', 'file not available, cannot do hotstart.'
            hotstart = .false.
            cycle file_loop
         else if (ext(loop) .eq. 'lm2') then
            print '(a)', 'file not available, cannot do snowcover.'
            cycle file_loop
         else
            print '(a)', 'file not available, but not mandatory.'
            cycle file_loop
         end if
      end if

      ! open the netcdf file and get the vertical dimension

      cdfid = ncopn(trim(input_laps_file), ncnowrit, rcode)

      zid = ncdid(cdfid, 'z', rcode)
      call ncdinq(cdfid, zid, dum, z, rcode)

      if ((ext(loop) .eq. 'lsx') .or. &
          (ext(loop) .eq. 'lm2')) then
         z2 = z
      else
         z3 = z
      end if

      if (loop .eq. 1) then

         ! allocate space for all of the variables in this data set.  note
         ! that some ofthe 3d variables are allocated by z+1 instead of just z, as
         ! we are going to put the surface values into these arrays as well.

         allocate (u(x, y, z3 + 1))
         allocate (v(x, y, z3 + 1))
         allocate (vv(x, y, z3 + 1))
         allocate (w(x, y, z3 + 1))
         allocate (t(x, y, z3 + 1))
         allocate (rh(x, y, z3 + 1))
         allocate (ht(x, y, z3 + 1))
         allocate (slp(x, y))
         allocate (psfc(x, y))
         allocate (d2d(x, y))
         allocate (tskin(x, y))
         allocate (tpw_before(x, y))
         allocate (tpw_after(x, y))
         allocate (p(z3 + 1))
         ! the following variables are not "mandatory"
         allocate (lwc(x, y, z3))
         allocate (rai(x, y, z3))
         allocate (sno(x, y, z3))
         allocate (pic(x, y, z3))
         allocate (ice(x, y, z3))
         if (moment_uphys .eq. 2) then
            allocate (mvd(x, y, z3))
         end if
         allocate (snocov(x, y))
         allocate (lcp(x, y, z3))
         ! the following variables are only used
         ! for converting non-mandatory cloud variables
         ! to mixing ratio values
         allocate (rho(x, y, z3))
         allocate (virtual_t(x, y, z3))
         allocate (sh(x, y, z3))
         allocate (mr(x, y, z3 + 1))
         allocate (soilt(x, y, 10))
         allocate (soilm(x, y, 10))

         ! initialize the non-mandatory values

         lwc(:, :, :) = -999.
         rai(:, :, :) = -999.
         sno(:, :, :) = -999.
         pic(:, :, :) = -999.
         snocov(:, :) = -999.
         lcp(:, :, :) = 0.
      end if

      if (ext(loop) .eq. 'lh3') then

         ! loop over the number of variables for this data file.

         var_lh3: do var_loop = 1, num_cdf_var(loop)

            ! get the variable id.

            vid = ncvid(cdfid, trim(cdf_var_name(var_loop, loop)), rcode)
            start = (/1, 1, 1, 1/)
            count = (/x, y, z, 1/)
            call ncvgt(cdfid, vid, start, count, rh, rcode)

            !  do this just once for pressure.

            vid = ncvid(cdfid, 'level', rcode)

            call ncvgt(cdfid, vid, 1, z3, p, rcode)

            ! set the pressure level of the lowest level of our
            ! pressure array as 2001 mb to flag the surface
            p(z3 + 1) = 2001

         end do var_lh3

      else if (ext(loop) .eq. 'lq3') then

         !  loop over the number of variables for this data file.
         var_lq3: do var_loop = 1, num_cdf_var(loop)

            !  get the variable id.

            vid = ncvid(cdfid, trim(cdf_var_name(var_loop, loop)), rcode)
            start = (/1, 1, 1, 1/)
            count = (/x, y, z, 1/)
            call ncvgt(cdfid, vid, start, count, sh, rcode)

         end do var_lq3

      else if (ext(loop) .eq. 'lsx') then

         !  loop over the number of variables for this data file.

         var_lsx: do var_loop = 1, num_cdf_var(loop)

            !  get the variable id.

            vid = ncvid(cdfid, trim(cdf_var_name(var_loop, loop)), rcode)
            start = (/1, 1, 1, 1/)
            count = (/x, y, 1, 1/)

            if (cdf_var_name(var_loop, loop) .eq. 'u  ') then
               call ncvgt(cdfid, vid, start, count, u(1, 1, z3 + 1), rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'v  ') then
               call ncvgt(cdfid, vid, start, count, v(1, 1, z3 + 1), rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'vv ') then
               call ncvgt(cdfid, vid, start, count, vv(1, 1, z3 + 1), rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 't  ') then
               call ncvgt(cdfid, vid, start, count, t(1, 1, z3 + 1), rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'rh ') then
               call ncvgt(cdfid, vid, start, count, rh(1, 1, z3 + 1), rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'mr ') then
               call ncvgt(cdfid, vid, start, count, mr(1, 1, z3 + 1), rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'msl') then
               call ncvgt(cdfid, vid, start, count, slp, rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'ps ') then
               call ncvgt(cdfid, vid, start, count, psfc, rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'tgd') then
               call ncvgt(cdfid, vid, start, count, tskin, rcode)
            end if

         end do var_lsx

         ! convert sfc mixing ratio from g/kg to kg/kg.

         mr(:, :, z3 + 1) = mr(:, :, z3 + 1)*0.001

         ! yuanfu: add lm1 input for soil moisture and possible future temperature:
      else if (ext(loop) .eq. 'lm1') then

         var_l1s: do var_loop = 1, num_cdf_var(loop)

            !  get the variable id.

            vid = ncvid(cdfid, trim(cdf_var_name(var_loop, loop)), rcode)
            start = (/1, 1, 1, 1/)
            count = (/x, y, 1, 1/)

            if (cdf_var_name(var_loop, loop) .eq. 'lsm') then
               call ncvgt(cdfid, vid, start, count, soilm, rcode)
            end if

         end do var_l1s

      else if (ext(loop) .eq. 'lm2') then

         var_l2s: do var_loop = 1, num_cdf_var(loop)

            !  get the variable id.

            vid = ncvid(cdfid, trim(cdf_var_name(var_loop, loop)), rcode)
            start = (/1, 1, 1, 1/)
            count = (/x, y, 1, 1/)

            if (cdf_var_name(var_loop, loop) .eq. 'sc ') then
               call ncvgt(cdfid, vid, start, count, snocov, rcode)
            end if

         end do var_l2s

      else if (ext(loop) .eq. 'lt1') then

         !  loop over the number of variables for this data file.

         var_lt1: do var_loop = 1, num_cdf_var(loop)

            !  get the variable id.

            vid = ncvid(cdfid, trim(cdf_var_name(var_loop, loop)), rcode)
            start = (/1, 1, 1, 1/)
            count = (/x, y, z, 1/)
            if (cdf_var_name(var_loop, loop) .eq. 't3 ') then
               call ncvgt(cdfid, vid, start, count, t, rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'ht ') then
               call ncvgt(cdfid, vid, start, count, ht, rcode)
            end if

         end do var_lt1

      else if (ext(loop) .eq. 'lw3') then

         !  loop over the number of variables for this data file.

         var_lw3: do var_loop = 1, num_cdf_var(loop)

            !  get the variable id.

            vid = ncvid(cdfid, trim(cdf_var_name(var_loop, loop)), rcode)
            start = (/1, 1, 1, 1/)
            count = (/x, y, z, 1/)
            if (cdf_var_name(var_loop, loop) .eq. 'u3 ') then
               call ncvgt(cdfid, vid, start, count, u, rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'v3 ') then
               call ncvgt(cdfid, vid, start, count, v, rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'om ') then
               call ncvgt(cdfid, vid, start, count, vv, rcode)
            end if

         end do var_lw3

      else if ((ext(loop) .eq. 'lwc') .and. (hotstart)) then

         !  loop over the number of variables for this data file.

         var_lwc1: do var_loop = 1, num_cdf_var(loop)

            !  get the variable id.

            vid = ncvid(cdfid, trim(cdf_var_name(var_loop, loop)), rcode)
            start = (/1, 1, 1, 1/)
            count = (/x, y, z, 1/)
            if (cdf_var_name(var_loop, loop) .eq. 'lwc') then
               call ncvgt(cdfid, vid, start, count, lwc, rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'rai') then
               call ncvgt(cdfid, vid, start, count, rai, rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'sno') then
               call ncvgt(cdfid, vid, start, count, sno, rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'ice') then
               call ncvgt(cdfid, vid, start, count, ice, rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'pic') then
               call ncvgt(cdfid, vid, start, count, pic, rcode)
            else if (cdf_var_name(var_loop, loop) .eq. 'mvd' .and. moment_uphys .eq. 2) then
               call ncvgt(cdfid, vid, start, count, mvd, rcode)
            end if

         end do var_lwc1

      else if ((ext(loop) .eq. 'lcp') .and. (hotstart)) then

         !  loop over the number of variables for this data file.

         var_lvc: do var_loop = 1, num_cdf_var(loop)

            !  get the variable id.

            vid = ncvid(cdfid, trim(cdf_var_name(var_loop, loop)), rcode)
            start = (/1, 1, 1, 1/)
            count = (/x, y, z, 1/)
            if (cdf_var_name(var_loop, loop) .eq. 'lcp') then
               call ncvgt(cdfid, vid, start, count, lcp, rcode)
               print *, 'got cloud cover...min/max = ', minval(lcp), maxval(lcp)
            end if

         end do var_lvc

      end if

   end do file_loop

   ! compute mixing ratio from spec hum.
   ! fill missing values with sfc value.

!   assume pressure array goes from lower pressures to higher ones
   k300 = 0
   do k = 1, z3
      if (p(k) .le. 300.) k300 = k
   end do
   if (k300 .eq. 0) then
      print *, "could not find k300!"
      stop
   else
      print *, "k300,p(k300) = ", k300, p(k300)
   end if

   do k = 1, z3
   do j = 1, y
   do i = 1, x

      if ((jaxsbn) .and. (p(k) .lt. 300.)) then
         weight_bot = (p(k) - 50)/(250)
         weight_top = 1.0 - weight_bot
         newsh = weight_bot*sh(i, j, k300) + &
                 weight_top*tiny

         newsh = min(sh(i, j, k), newsh)

         ! make sure sh does not exceed
         ! ice saturation value
         call saturate_ice_points(t(i, j, k), &
                                  p(k), 1.0, &
                                  shmod, rhmod)
         sh(i, j, k) = min(shmod, newsh)
         rh(i, j, k) = make_rh(p(k), t(i, j, k) - 273.15, &
                               sh(i, j, k)*1000., -132.)*100.
      end if

      if (sh(i, j, k) .ge. 0. .and. sh(i, j, k) .lt. 1.) then
         mr(i, j, k) = sh(i, j, k)/(1.-sh(i, j, k))
      else
         mr(i, j, k) = mr(i, j, z3 + 1)
      end if
!      if (u(i,j,k) .eq. 1.e-30 .or. abs(u(i,j,k)) .gt. 200.) u(i,j,k)=u(i,j,z3+1)
!     if (v(i,j,k) .eq. 1.e-30 .or. abs(v(i,j,k)) .gt. 200.) v(i,j,k)=v(i,j,z3+1)
   end do
   end do
   end do

   !  set the lowest level of the geopotential height to topographic height

   ht(:, :, z3 + 1) = topo

   if (hotstart) then

      ! if this is a hot start, then we need to convert the microphysical
      ! species from mass per volume to mass per mass (mixing ratio).  this
      ! requires that we compute the air density from virtual temperature
      ! and divide each species by the air density.

      ! compute virtual temperature from mixing ratio and temperature
      virtual_t(:, :, :) = (1.+0.61*mr(:, :, 1:z3))*t(:, :, 1:z3)

      ! compute density from virtual temperature and gas constant for dry air
      do k = 1, z3
         rho(:, :, k) = p(k)*100./(rdry*virtual_t(:, :, k))
      end do

      ! for each of the species, ensure they are not "missing".  if missing
      ! then set their values to 0.000.  otherwise, divide by the density to
      ! convert from concentration to mixing ratio.

      print *, 'max lwc concentration (as read in) = ', maxval(lwc)

      if (maxval(lwc) .lt. 99999.) then

         ! scale lwc for grid spacing
         lwc = lwc*hydrometeor_scale_cld
         print *, 'max lwc concentration (after hydro scaling) = ', maxval(lwc)

         ! cap lwc to autoconversion rate for liquid to rain
         where (lwc .gt. autoconv_lwc2rai) lwc = autoconv_lwc2rai
         print *, 'max lwc concentration (after autoconversion) = ', maxval(lwc)

         ! convert lwc concentration to mixing ratio
         lwc(:, :, :) = lwc(:, :, :)/rho(:, :, :)   ! cloud liquid mixing ratio
         print *, 'max lwc mixing ratio (before saturate_lwc_points) = ', maxval(lwc)
         print *, 'max sh (before saturate_lwc_points) = ', maxval(sh)

         ! calculate tpw from initial fields
         do i = 1, x
         do j = 1, y
            sh_sfc = mr(i, j, z3 + 1)/(1.+mr(i, j, z3 + 1))
            call integrate_tpw(sh(i, j, z3:1:-1), sh_sfc, p(z3:1:-1)*100., psfc(i, j), 1, 1, z3, tpw_before(i, j))
         end do ! j
         end do ! i
         print *, 'max tpw (before saturate_lwc_points) = ', maxval(tpw_before)

         tpw_after = tpw_before ! initialize

         ! convert lwc mixing ratio to vapor mixing ratio
         if (lwc2vapor_thresh .gt. 0.) then
            do k = 1, z3
               do j = 1, y
                  do i = 1, x
                     if ((lcp(i, j, k) .ge. lcp_min) .and. &
                         (t(i, j, k) .ge. 263.0) .and. &
                         (lwc(i, j, k) .gt. lwc_min)) then
                        !if (lwc(i,j,k).gt.0.00010) then
                        !call lwc2vapor(lwc(i,j,k),sh(i,j,k),t(i,j,k), &
                        !               p(k),lwc2vapor_thresh, &
                        !               lwcmod,shmod,rhmod)

                        call saturate_lwc_points(sh(i, j, k), t(i, j, k), &
                                                 p(k), lwc2vapor_thresh, &
                                                 shmod, rhmod)
                        ! update moisture arrays

                        sh(i, j, k) = shmod
                        rh(i, j, k) = rhmod
                        mr(i, j, k) = shmod/(1.-shmod)
                     end if
                  end do
               end do
            end do
            ! calculate tpw from modified fields
            do i = 1, x
            do j = 1, y
               sh_sfc = mr(i, j, z3 + 1)/(1.+mr(i, j, z3 + 1))
               call integrate_tpw(sh(i, j, z3:1:-1), sh_sfc, p(z3:1:-1)*100., psfc(i, j), 1, 1, z3, tpw_after(i, j))
            end do ! j
            end do ! i
         end if

         print *, 'max lwc mixing ratio (after saturate_lwc_points) = ', maxval(lwc)
         print *, 'max sh (after saturate_lwc_points) = ', maxval(sh)
         print *, 'max tpw (after saturate_lwc_points) = ', maxval(tpw_after)
      else
         print *, 'missing cloud liquid, setting values to 0.0'
         lwc(:, :, :) = 0.0
      end if

      print *, 'max rain concentration = ', maxval(rai)
      if (maxval(rai) .lt. 99999.) then
         rai(:, :, :) = rai(:, :, :)*hydrometeor_scale_pcp
         rai(:, :, :) = rai(:, :, :)/rho(:, :, :)   ! rain mixing ratio
      else
         print *, 'missing rain, setting values to 0.0'
         rai(:, :, :) = 0.0
      end if
      print *, 'max rain mixing ratio = ', maxval(rai)

      print *, 'max snow concentration = ', maxval(sno)
      if (maxval(sno) .lt. 99999.) then
         sno(:, :, :) = sno(:, :, :)*hydrometeor_scale_pcp
         sno(:, :, :) = sno(:, :, :)/rho(:, :, :)   ! snow mixing ratio
      else
         print *, 'missing snow, setting values to 0.0'
         sno(:, :, :) = 0.0
      end if
      print *, 'max snow mixing ratio = ', maxval(sno)

      if (maxval(ice) .lt. 99999.) then
         ! limit ice to autoconversion threshold

         ice(:, :, :) = ice(:, :, :)*hydrometeor_scale_cld
         where (ice .gt. autoconv_ice2sno) ice = autoconv_ice2sno
         ice(:, :, :) = ice(:, :, :)/rho(:, :, :)   ! ice mixing ratio

         ! convert ice mixing ratio to vapor mixing ratio

         if (ice2vapor_thresh .gt. 0.) then
            do k = 1, z3
               do j = 1, y
                  do i = 1, x
                     if ((lcp(i, j, k) .ge. lcp_min) .and. &
                         (t(i, j, k) .lt. 263.) .and. &
                         (ice(i, j, k) .gt. ice_min)) then

                        ! update moisture arrays
                        call saturate_ice_points(t(i, j, k), &
                                                 p(k), ice2vapor_thresh, &
                                                 shmod, rhmod)
                        sh(i, j, k) = shmod
                        rh(i, j, k) = rhmod
                        mr(i, j, k) = sh(i, j, k)/(1.-sh(i, j, k))
                     end if
                  end do
               end do
            end do
            ! calculate tpw from modified fields
            do i = 1, x
            do j = 1, y
               sh_sfc = mr(i, j, z3 + 1)/(1.+mr(i, j, z3 + 1))
               call integrate_tpw(sh(i, j, z3:1:-1), sh_sfc, p(z3:1:-1)*100., psfc(i, j), 1, 1, z3, tpw_after(i, j))
            end do ! j
            end do ! i
         end if
         print *, 'max tpw (after saturate_ice_points) = ', maxval(tpw_after)
      else
         print *, 'missing ice, setting values to 0.0'
         ice(:, :, :) = 0.0
      end if

      print *, 'max pice concentration = ', maxval(ice)
      if (maxval(pic) .lt. 99999.) then
         pic(:, :, :) = pic(:, :, :)*hydrometeor_scale_pcp
         pic(:, :, :) = pic(:, :, :)/rho(:, :, :)   ! graupel (precipitating ice) mixing rat.
      else
         print *, 'missing pice, setting values to 0.0'
         pic(:, :, :) = 0.0
      end if
      print *, 'max pice mixing ratio = ', maxval(ice)

      ! keep 3d omega as pa/s, or fill with sfc value if missing.

      do k = 1, z3
      do j = 1, y
      do i = 1, x
         if (vv(i, j, k) .eq. 1.e-30 .or. vv(i, j, k) .eq. r_missing_data .or. abs(vv(i, j, k)) .gt. 100.) then
            vv(i, j, k) = vv(i, j, z3 + 1)
         end if

!       convert from omega to w
         if (vv(i, j, k) .ne. r_missing_data) then
            w(i, j, k) = -vv(i, j, k)/(rho(i, j, k)*g)
         else
            w(i, j, k) = 0.
         end if

      end do
      end do
      end do

   end if

   ! if make_sfc_uv set, then replace surface winds with
   ! winds interpolated from the 3d field.
   if (make_sfc_uv) then
      print *, 'creating surface u/v from 3d field...'
      do j = 1, y
         do i = 1, x
            kbot = 0
            get_lowest: do k = z3, 1, -1
               if (ht(i, j, k) .gt. topo(i, j)) then
                  kbot = k
                  exit get_lowest
               end if
            end do get_lowest
            if (kbot .ne. 0) then
               u(i, j, z3 + 1) = u(i, j, kbot)
               v(i, j, z3 + 1) = v(i, j, kbot)
            else
               print *, 'problem finding kbot.'
               stop
            end if
         end do
      end do

   end if

   ! loop over the each desired output format

   do out_loop = 1, num_output
      ! now it is time to output these arrays.  the arrays are ordered
      !  as (x,y,z).  the origin is the southwest corner at the top of the
      ! atmosphere for the 3d arrays, where the last layer (z3+1) contains
      ! the surface information.  this is where you would insert a call
      ! to a custom output routine.

      select_output:select case(output_format(out_loop))
      case ('mm5 ')
      call output_pregrid_format(p, t, ht, u, v, rh, slp, &
                                 lwc, rai, sno, ice, pic, snocov, tskin)

      case ('wrf ')
      call output_gribprep_format(p, t, ht, u, v, rh, slp, psfc, &
                                  lwc, rai, sno, ice, pic, snocov, tskin)
      case ('wps ')
      call output_metgrid_format(p, t, ht, u, v, w, rh, slp, psfc, &
                                 lwc, rai, sno, ice, pic, snocov, tskin, soilt, soilm)

      case ('rams')
      call output_ralph2_format(p, u, v, t, ht, rh, slp, psfc, snocov, tskin)
      case ('sfm ')
      print '(a)', 'support for sfm (rams 3b) coming soon...check back later!'

      case ('cdf ')
!lw added rh (with srh in z+1), snocov and tskin
!         call output_netcdf_format(p,ht,t,mr,u,v,w,slp,psfc,lwc,ice,rai,sno,pic)
      call output_netcdf_format(p, ht, t, mr, u, v, w, slp, psfc, lwc, ice, rai, sno, &
                                pic, rh, snocov, tskin)

      case default
      print '(2a)', 'unrecognized output format: ', output_format
      print '(a)', 'recognized formats include mm5, rams, wrf, sfm, and cdf'

      end select select_output
   end do
   print '(a)', 'lapsprep complete.'

end program lapsprep

