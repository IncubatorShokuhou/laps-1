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

program wfoprep

   ! this is a new version of the wfoprep program.  unlike the previous
   ! version of this programs distributed with laps, this version is
   ! focused on acquiring data for lateral boundary conditions to support
   ! a mesoscale nwp model run.  for initialization, the new program
   ! "lapsprep" builds the files necessary for the initial condition, assuming
   ! you wish to initialize with laps.  furthermore, the old version was
   ! set up to support only the sfm, which we are no longer supporting.  this
   ! version can output files for ingest by mm5v3/regridder, rams 4.3 in ralph2
   ! format, or wrfsi/hinterp input.   this is controlled by the output_type
   ! entry(ies) in the wfoprep.nl.  finally, this new version is cast in
   ! free-form f90 code, including the use of modules.
   !
   ! history
   ! -------
   ! sep 2001:  new f90 version for mm5v3/rams4.3/wrfsi support.
   !            b. shaw, noaa/cira/fsl frd/lapb
   !
   !
   use wfoprep_setup
   use wfo_models
   use map_utils
   use wfoprep_mm5
   use wfoprep_wrf

   implicit none
   real, parameter           :: rmissingval = -9999.
   integer, parameter         :: maxl = 75
   integer, parameter         :: maxt = 50
   logical                   :: filefound
   integer                   :: fcstsec(maxt)
   integer                   :: goodlevs
   real                      :: goodpct
   integer                   :: i, ii
   integer                   :: i4time_offset
   integer                   :: i4time_last
   integer                   :: i4time_valid, i4time_valid1, i4time_valid2
   integer                   :: i4time_cycle
   integer, allocatable      :: i4times_avail_max(:)
   integer, allocatable      :: i4times_avail(:)
   integer                   :: istatus
   integer                   :: k, kk
   character(len=13)        :: last_cycle_processed
   character(len=13)        :: latest_cycle_wfo

   integer                   :: m
   character(len=256)       :: modelfile_wfo
   integer                   :: nfid
   integer                   :: ntimes, ntimes_needed
   integer                   :: outfreq_sec
   character(len=256)       :: proclog
   integer                   :: ta
   integer                   :: t_id, ht_id, u_id, v_id, rh_id, msl_id
   integer                   :: nz_t, nz_ht, nz_u, nz_v, nz_rh, nz_msl
   integer                   :: np_t, np_ht, np_u, np_v, np_rh, np_rh_raw, np_msl
   logical                   :: havesfc_t, havesfc_ht, havesfc_u, havesfc_msl
   logical                   :: havesfc_v, havesfc_rh
   character(len=10)        :: t_levels_c(maxl)
   character(len=10)        :: ht_levels_c(maxl)
   character(len=10)        :: rh_levels_c(maxl)
   character(len=10)        :: u_levels_c(maxl)
   character(len=10)        :: v_levels_c(maxl)
   character(len=10)        :: msl_levels_c(maxl)
   real                      :: t_plevels(maxl)
   real                      :: ht_plevels(maxl)
   real                      :: rh_plevels(maxl), rh_plevels_raw(maxl)
   real                      :: u_plevels(maxl)
   real                      :: v_plevels(maxl)
   real                      :: msl_plevels(maxl)
   integer                   :: t_kbotp, t_ktopp, t_ksfc
   integer                   :: rh_kbotp, rh_ktopp, rh_ksfc
   integer                   :: ht_kbotp, ht_ktopp, ht_ksfc
   integer                   :: u_kbotp, u_ktopp, u_ksfc
   integer                   :: v_kbotp, v_ktopp, v_ksfc
   integer                   :: msl_kbotp, msl_ktopp, msl_ksfc
   character(len=13)        :: wfofname
   character(len=13), external :: cvt_i4time_wfo_fname13
   logical, allocatable      :: t_inv(:, :)
   logical, allocatable      :: u_inv(:, :)
   logical, allocatable      :: v_inv(:, :)
   logical, allocatable      :: rh_inv(:, :)
   logical, allocatable      :: rh_inv_raw(:, :)
   logical, allocatable      :: ht_inv(:, :)
   logical, allocatable      :: msl_inv(:, :)
   logical, allocatable      :: goodtime_flag(:)
   integer                   :: n_goodtimes
   character(len=10)         :: mslname
   integer, allocatable      :: time_index(:)
   integer                   :: nfssttimes
   type(proj_info)           :: proj
   real                      :: weight1
   ! data arrays
   real, allocatable         :: z3d(:, :, :), z3d1(:, :, :), z3d2(:, :, :)
   real, allocatable         :: zsf(:, :), zsf1(:, :), zsf2(:, :)
   real, allocatable         :: t3d(:, :, :), t3d1(:, :, :), t3d2(:, :, :)
   real, allocatable         :: tsf(:, :), tsf1(:, :), tsf2(:, :)
   real, allocatable         :: rh3d(:, :, :), rh3d1(:, :, :), rh3d2(:, :, :)
   real, allocatable         :: rhsf(:, :), rhsf1(:, :), rhsf2(:, :)
   real, allocatable         :: u3d(:, :, :), u3d1(:, :, :), u3d2(:, :, :)
   real, allocatable         :: usf(:, :), usf1(:, :), usf2(:, :)
   real, allocatable         :: v3d(:, :, :), v3d1(:, :, :), v3d2(:, :, :)
   real, allocatable         :: vsf(:, :), vsf1(:, :), vsf2(:, :)
   real, allocatable         :: slp(:, :), slp1(:, :), slp2(:, :)
   real, allocatable         :: fsstsum(:, :), fsst(:, :)
   real, allocatable         :: topo(:, :)
   real, allocatable         :: data3d_temp(:, :, :)
   real                      :: tdsf_c, psf_mb
   real, external            :: dwpt_laps, twet_fast
   integer                   :: ix, jy
   logical                   :: fixvar1, fixvar2

   istatus = 1
   print '(a)', 'calling setup routine...'
   call setup_wfoprep(istatus)
   if (istatus .ne. 1) then
      print '(a)', 'dprep: problem with setup configuration.'
      stop
   end if

   ! at this point, the wfoprep_setup module contains a list of models
   ! and information relating to their availability (e.g., frequency
   ! of model runs, delay time before a cycle is available, etc.).  we
   ! will loop through each valid model and process as requested and
   ! available.

   model_loop: do m = 1, num_models

      outfreq_sec = output_freq(m)*3600
      print '(2a)', '*** processing model ', trim(model_name(m))
      ! based on the namelist entries that tell us how frequently this
      ! model is run and how many hours after cycle time it becomes available,
      ! we can do some computations to see which cycle we should be trying to
      ! process.

      ! first, compute the i4time_offset value, which is the current time
      ! minus the amount of time between this model's cycle time and when
      ! the entire run is available.

      i4time_offset = i4time_now - model_delay(m)*3600

      ! round i4time_offset to nearest hour
      i4time_offset = i4time_offset - mod(i4time_offset, 3600)

      ! now, compute the closest cycle time for this model that is equal
      ! to or earlier than the i4time_offset time.

      i4time_cycle = (i4time_offset/(model_run_freq(m)*3600))* &
                     model_run_freq(m)*3600
      i4time_last = i4time_cycle + max_fcst_len(m)*3600

      ! what is the cycle time closest?
      latest_cycle_wfo = cvt_i4time_wfo_fname13(i4time_cycle)
      print *, '    latest possible cycle for this model = ', latest_cycle_wfo

      ! check the {model}.last file in the drep output directory
      ! to see if this cycle has been processed yet.

      print *, '    checking log of previously processed cycles.'
      proclog = trim(ext_data_path)//trim(output_name(m))//'.last'
      inquire (file=proclog, exist=filefound)
      if (filefound) then
         open (file=proclog, unit=11, status='old', form='formatted', &
               access='sequential')
         read (11, '(a13)') last_cycle_processed
         close (11)
         print '(2a)', '       last cycle processed for this model:', &
            last_cycle_processed
         if (last_cycle_processed .ge. latest_cycle_wfo) then
            print '(a)', '       current cycle <= last cycle processed. skipping...'
            cycle model_loop
         end if
      else
         print '(a)', '    no record of processing for any previous cycles.'
      end if

      ! see if the file required for this cycle is available
      modelfile_wfo = trim(model_path(m))//'/'//latest_cycle_wfo

      print '(2a)', '    searching for data file: ', trim(modelfile_wfo)
      inquire (file=modelfile_wfo, exist=filefound)
      if (.not. filefound) then
         print '(a)', ' '
         print '(a)', '#############################'
         print '(a)', 'skipping model: not available'
         print '(a)', '#############################'
         print '(a)', ' '
         cycle model_loop
      end if

      ! open the file
      call open_wfofile(modelfile_wfo, nfid, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '#############################'
         print '(a)', 'skipping model: i/o error'
         print '(a)', '#############################'
         print '(a)', ' '
         cycle model_loop
      end if

      ! get grid/projection info
      call get_wfomodel_proj(nfid, model_name(m), proj, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '#################################'
         print '(a)', 'skipping model: proj info problem'
         print '(a)', '#################################'
         print '(a)', ' '
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if

      ! process time info
      call get_wfomodel_fcsttimes(nfid, ntimes, fcstsec, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '#################################'
         print '(a)', 'skipping model: time info problem'
         print '(a)', '#################################'
         print '(a)', ' '
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if
      allocate (i4times_avail_max(ntimes))
      i4times_avail_max(1:ntimes) = i4time_cycle + fcstsec(1:ntimes)
      if (i4time_last .gt. i4times_avail_max(ntimes)) then
         print *, 'info: this source does not support max_fcst length of ', max_fcst_len(m)
         print *, 'info: resetting to ', (i4times_avail_max(ntimes) - i4time_cycle)/36
         i4time_last = i4times_avail_max(ntimes)
      else
         ! compute the ntimes_needed value, which reprenents the number
         ! of time periods actually needed to get us equal to or just past
         ! i4time_last
         ntimes_needed = 1
         find_ntimes_needed: do i = 2, ntimes
            ntimes_needed = ntimes_needed + 1
            if (i4times_avail_max(i) .ge. i4time_last) then
               exit find_ntimes_needed
            end if
         end do find_ntimes_needed

      end if

      ! get level information for state variables (t,ht,u,v,rh) and then
      ! get the inventory

      ! temperature
      call get_wfomodel_var_levels(nfid, 't         ', t_id, nz_t, t_levels_c, np_t, &
                                   t_plevels, t_kbotp, t_ktopp, havesfc_t, t_ksfc, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '#################################'
         print '(a)', 'skipping model: no t data        '
         print '(a)', '#################################'
         print '(a)', ' '
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if

      ! height is only mandatory for 3d data sets
      if (model_code(m) .gt. 1) then
         call get_wfomodel_var_levels(nfid, 'gh        ', ht_id, nz_ht, ht_levels_c, np_ht, &
                                      ht_plevels, ht_kbotp, ht_ktopp, havesfc_ht, ht_ksfc, istatus)
         if (istatus .ne. 1) then
            print '(a)', ' '
            print '(a)', '#################################'
            print '(a)', 'skipping model: no z data        '
            print '(a)', '#################################'
            print '(a)', ' '
            call close_wfofile(nfid, istatus)
            cycle model_loop
         end if
      end if

      ! relative humidity
      call get_wfomodel_var_levels(nfid, 'rh        ', rh_id, nz_rh, rh_levels_c, np_rh_raw, &
                                   rh_plevels_raw, rh_kbotp, rh_ktopp, havesfc_rh, rh_ksfc, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '#################################'
         print '(a)', 'skipping model: no rh data        '
         print '(a)', '#################################'
         print '(a)', ' '
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if

      ! u-wind component
      call get_wfomodel_var_levels(nfid, 'uw        ', u_id, nz_u, u_levels_c, np_u, &
                                   u_plevels, u_kbotp, u_ktopp, havesfc_u, u_ksfc, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '#################################'
         print '(a)', 'skipping model: no u data        '
         print '(a)', '#################################'
         print '(a)', ' '
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if

      ! v-wind component
      call get_wfomodel_var_levels(nfid, 'vw        ', v_id, nz_v, v_levels_c, np_v, &
                                   v_plevels, v_kbotp, v_ktopp, havesfc_v, v_ksfc, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '#################################'
         print '(a)', 'skipping model: no v data        '
         print '(a)', '#################################'
         print '(a)', ' '
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if

      ! mslp is mandatory for 3d data sets.  check for eta mslp first.  if non-existent,
      ! then check for standard mslp

      if (model_code(m) .gt. 1) then
         call get_wfomodel_var_levels(nfid, 'emsp      ', msl_id, nz_msl, &
                                      msl_levels_c, np_msl, &
                                      msl_plevels, msl_kbotp, msl_ktopp, havesfc_msl, msl_ksfc, istatus)
         if (istatus .eq. 1) then
            print *, 'using eta mslp field'
            mslname = 'emsp      '
         else
            call get_wfomodel_var_levels(nfid, 'pmsl      ', msl_id, nz_msl, &
                                         msl_levels_c, np_msl, &
                                         msl_plevels, msl_kbotp, msl_ktopp, havesfc_msl, msl_ksfc, istatus)
            if (istatus .eq. 1) then
               print *, 'using pmsl field'
               mslname = 'pmsl      '
            else
               call get_wfomodel_var_levels(nfid, 'mmsp      ', msl_id, nz_msl, &
                                            msl_levels_c, np_msl, &
                                            msl_plevels, msl_kbotp, msl_ktopp, havesfc_msl, msl_ksfc, istatus)
               if (istatus .eq. 1) then
                  print *, 'using mmsp field (maps sea-level pressure)'
                  mslname = 'mmsp      '
               else
                  print '(a)', ' '
                  print '(a)', '#################################'
                  print '(a)', 'skipping model: no mslp data        '
                  print '(a)', '#################################'
                  print '(a)', ' '
                  call close_wfofile(nfid, istatus)
                  cycle model_loop
               end if
            end if
         end if
      end if

      ! now, check to make sure we have t,ht, u, and v for all of the same
      ! pressure levels if this is a 3d data set.  this seems to be the case
      ! for all supported datasets.  however, rh is a bit trickier.  in some
      ! cases, rh is not provided above 300 mb or below 975 mb.  deal with this
      ! later.

      if (model_code(m) .gt. 1) then
         if ((np_t .ne. np_ht) .or. (np_u .ne. np_ht) .or. &
             (np_u .ne. np_ht)) then
            print '(a)', ' '
            print '(a)', '####################################'
            print '(a)', 'skipping model: dimension mismatches'
            print '(a)', '####################################'
            print '(a)', ' '
            deallocate (i4times_avail_max)
            call close_wfofile(nfid, istatus)
            cycle model_loop
         end if
      end if

      ! we have made it this far, so we must have variables we need, but we
      ! need to check their inventory to see which time periods are complete.
      ! if this is a 3d data source (model code > 1), then we need to have at
      ! least the top and bottom level values for each variable for the time
      ! to be included in the process.  also must have the first and last
      ! time periods.

      n_goodtimes = 0
      if (allocated(goodtime_flag)) deallocate (goodtime_flag)
      allocate (goodtime_flag(ntimes_needed))
      goodtime_flag(:) = .true.

      if (allocated(t_inv)) deallocate (t_inv)
      allocate (t_inv(nz_t, ntimes))
      call get_wfomodel_var_inv(nfid, 't         ', nz_t, ntimes, t_inv, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '####################################'
         print '(a)', 'skipping model: no t inventory      '
         print '(a)', '####################################'
         print '(a)', ' '
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if

      if (allocated(u_inv)) deallocate (u_inv)
      allocate (u_inv(nz_u, ntimes))
      call get_wfomodel_var_inv(nfid, 'uw        ', nz_u, ntimes, u_inv, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '####################################'
         print '(a)', 'skipping model: no u inventory      '
         print '(a)', '####################################'
         print '(a)', ' '
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if

      if (allocated(v_inv)) deallocate (v_inv)
      allocate (v_inv(nz_v, ntimes))
      call get_wfomodel_var_inv(nfid, 'vw        ', nz_v, ntimes, v_inv, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '####################################'
         print '(a)', 'skipping model: no v inventory      '
         print '(a)', '####################################'
         print '(a)', ' '
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if

      if (allocated(rh_inv_raw)) deallocate (rh_inv_raw)
      allocate (rh_inv_raw(nz_rh, ntimes))
      ! also allocate an array that matches the temperature array
      if (allocated(rh_inv)) deallocate (rh_inv)
      allocate (rh_inv(nz_t, ntimes))
      rh_inv(:, :) = .false.
      if (nz_rh .eq. nz_t) rh_inv = rh_inv_raw
      call get_wfomodel_var_inv(nfid, 'rh        ', nz_rh, ntimes, rh_inv_raw, istatus)
      if (istatus .ne. 1) then
         print '(a)', ' '
         print '(a)', '####################################'
         print '(a)', 'skipping model: no rh inventory      '
         print '(a)', '####################################'
         print '(a)', ' '
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if

      if (allocated(ht_inv)) deallocate (ht_inv)
      if (model_code(m) .gt. 1) then
         allocate (ht_inv(nz_ht, ntimes))
         call get_wfomodel_var_inv(nfid, 'gh        ', nz_ht, ntimes, ht_inv, istatus)
         if (istatus .ne. 1) then
            print '(a)', ' '
            print '(a)', '####################################'
            print '(a)', 'skipping model: no z inventory      '
            print '(a)', '####################################'
            print '(a)', ' '
            call close_wfofile(nfid, istatus)
            cycle model_loop
         end if
      end if

      if (allocated(msl_inv)) deallocate (msl_inv)
      if (model_code(m) .gt. 1) then
         allocate (msl_inv(nz_msl, ntimes))
         call get_wfomodel_var_inv(nfid, mslname, nz_msl, ntimes, msl_inv, istatus)
         if (istatus .ne. 1) then
            print '(a)', ' '
            print '(a)', '####################################'
            print '(a)', 'skipping model: no mslp inventory      '
            print '(a)', '####################################'
            print '(a)', ' '
            call close_wfofile(nfid, istatus)
            cycle model_loop
         end if
      end if

      ! here is the actual check for minimum requirements.
      ! minimum variable requirements:
      !  ht, t, u, v, rh on pressure levels + mslp
      ! for now, require that 100% of the levels be present for each
      ! variable, and that both the top and bottom level must be
      ! present.  in the future, we can add vertical interpolatoin
      ! to allow for missing levels, thus allowing us to reduce the
      ! 100% level threshold

      invloop: do i = 1, ntimes_needed

         if (model_code(m) .gt. 1) then

            ! check msl pressure
            if (.not. msl_inv(msl_ksfc, i)) then
               print *, 'warning: missing msl for this time:', i
               goodtime_flag(i) = .false.
               cycle invloop
            end if

            ! check height field, requiring most levels to be present!
            goodlevs = 0
            htinvloop: do k = ht_kbotp, ht_ktopp
               if (ht_inv(k, i)) goodlevs = goodlevs + 1
            end do htinvloop
            goodpct = float(goodlevs)/float(np_ht)
            if ((goodpct .lt. min_vert_frac) .or. (.not. ht_inv(ht_kbotp, i)) .or. &
                (.not. ht_inv(ht_ktopp, i))) then
               print *, 'warning: height inventory failed vertical check:', &
                  goodpct, i
               goodtime_flag(i) = .false.
               cycle invloop
            end if

            ! check temperatures
            goodlevs = 0
            tinvloop: do k = t_kbotp, t_ktopp
               if (t_inv(k, i)) goodlevs = goodlevs + 1
            end do tinvloop
            goodpct = float(goodlevs)/float(np_t)
            if ((goodpct .lt. min_vert_frac) .or. (.not. t_inv(t_kbotp, i)) .or. &
                (.not. t_inv(t_ktopp, i))) then
               print *, 'warning: temperature inventory failed vertical check:', &
                  goodpct, i
               goodtime_flag(i) = .false.
               cycle invloop
            end if

            ! check rh
            goodlevs = 0
            rhinvloop: do k = rh_kbotp, rh_ktopp
               if (rh_inv_raw(k, i)) goodlevs = goodlevs + 1
            end do rhinvloop
            goodpct = float(goodlevs)/float(np_rh_raw)
            if ((goodpct .lt. min_vert_frac) .or. (.not. rh_inv_raw(rh_kbotp, i)) .or. &
                (.not. rh_inv_raw(rh_ktopp, i))) then
               print *, 'warning: rh inventory failed vertical check:', &
                  goodpct, i
               goodtime_flag(i) = .false.
               cycle invloop
            end if

            ! check u
            goodlevs = 0
            uinvloop: do k = u_kbotp, u_ktopp
               if (u_inv(k, i)) goodlevs = goodlevs + 1
            end do uinvloop
            goodpct = float(goodlevs)/float(np_u)
            if ((goodpct .lt. min_vert_frac) .or. (.not. u_inv(u_kbotp, i)) .or. &
                (.not. u_inv(u_ktopp, i))) then
               print *, 'warning: u inventory failed vertical check:', &
                  goodpct, i
               goodtime_flag(i) = .false.
               cycle invloop
            end if

            ! check v
            goodlevs = 0
            vinvloop: do k = v_kbotp, v_ktopp
               if (v_inv(k, i)) goodlevs = goodlevs + 1
            end do vinvloop
            goodpct = float(goodlevs)/float(np_v)
            if ((goodpct .lt. min_vert_frac) .or. (.not. v_inv(v_kbotp, i)) .or. &
                (.not. v_inv(v_ktopp, i))) then
               print *, 'warning: v inventory failed vertical check:', &
                  goodpct, i
               goodtime_flag(i) = .false.
               cycle invloop
            end if

         end if
         if ((model_code(m) .eq. 1) .or. (model_code(m) .eq. 3)) then

            ! we need surface values

            if (.not. t_inv(t_ksfc, i)) then
               print *, 'warning: missing surface t', i
               goodtime_flag(i) = .false.
               cycle invloop
            end if

            if (.not. u_inv(u_ksfc, i)) then
               print *, 'warning: missing surface u', i
               goodtime_flag(i) = .false.
               cycle invloop
            end if

            if (.not. v_inv(u_ksfc, i)) then
               print *, 'warning: missing surface v', i
               goodtime_flag(i) = .false.
               cycle invloop
            end if

            if (.not. rh_inv_raw(rh_ksfc, i)) then
               print *, 'warning: missing surface rh', i
               goodtime_flag(i) = .false.
               cycle invloop
            end if

         end if
         ! if we made it this far without cycling the loop, then
         ! this is a good level
         n_goodtimes = n_goodtimes + 1
         goodtime_flag(i) = .true.
      end do invloop

      ! ok, now lets make sure that enough time periods
      ! passed the above check and that the first and last
      ! times are available
      goodlevs = 0
      do i = 1, ntimes_needed
         if (goodtime_flag(i)) goodlevs = goodlevs + 1
      end do
      goodpct = float(goodlevs)/float(ntimes_needed)
      if ((goodpct .lt. min_time_frac) .or. (.not. goodtime_flag(1)) .or. &
          (.not. goodtime_flag(ntimes_needed))) then
         print *, ' '
         print *, '###################################################'
         print *, 'skipping model:  time inventory check failed'
         print *, 'goodpct/min_time_frac = ', goodpct, min_time_frac
         print *, 'goodtime_flag(1) = ', goodtime_flag(1)
         print *, 'ntimes_needed = ', ntimes_needed
         print *, 'goodtime_flag(ntimes_needed) =', goodtime_flag(ntimes_needed)
         print *, '###################################################'
         call close_wfofile(nfid, istatus)
         cycle model_loop
      end if

      ! allocate the i4times_avail array to the actual number
      ! of usable time periods and populate it
      if (allocated(i4times_avail)) deallocate (i4times_avail)
      allocate (i4times_avail(n_goodtimes))
      if (allocated(time_index)) deallocate (time_index)
      allocate (time_index(n_goodtimes))
      ii = 1
      do i = 1, ntimes_needed
         if (goodtime_flag(i)) then
            i4times_avail(ii) = i4times_avail_max(i)
            time_index(ii) = i
            ii = ii + 1
         end if
      end do
      deallocate (i4times_avail_max)

      ! allocate 3 arrays for each state variable.  two will be needed for the
      ! data at two times bounding the time of interest.  the third is where
      ! the time interpolated values will be saved.

      if (model_code(m) .gt. 1) then
         allocate (z3d(proj%nx, proj%ny, np_ht))
         allocate (z3d1(proj%nx, proj%ny, np_ht))
         allocate (z3d2(proj%nx, proj%ny, np_ht))
         allocate (t3d(proj%nx, proj%ny, np_t))
         allocate (t3d1(proj%nx, proj%ny, np_t))
         allocate (t3d2(proj%nx, proj%ny, np_t))
         allocate (rh3d(proj%nx, proj%ny, np_t))
         allocate (rh3d1(proj%nx, proj%ny, np_rh_raw))
         allocate (rh3d2(proj%nx, proj%ny, np_rh_raw))
         allocate (u3d(proj%nx, proj%ny, np_u))
         allocate (u3d1(proj%nx, proj%ny, np_u))
         allocate (u3d2(proj%nx, proj%ny, np_u))
         allocate (v3d(proj%nx, proj%ny, np_v))
         allocate (v3d1(proj%nx, proj%ny, np_v))
         allocate (v3d2(proj%nx, proj%ny, np_v))
         allocate (slp(proj%nx, proj%ny))
         allocate (slp1(proj%nx, proj%ny))
         allocate (slp2(proj%nx, proj%ny))
      end if
      if ((model_code(m) .eq. 1) .or. (model_code(m) .eq. 3)) then
         allocate (tsf(proj%nx, proj%ny))
         allocate (tsf1(proj%nx, proj%ny))
         allocate (tsf2(proj%nx, proj%ny))
         allocate (rhsf(proj%nx, proj%ny))
         allocate (rhsf1(proj%nx, proj%ny))
         allocate (rhsf2(proj%nx, proj%ny))
         allocate (usf(proj%nx, proj%ny))
         allocate (usf1(proj%nx, proj%ny))
         allocate (usf2(proj%nx, proj%ny))
         allocate (vsf(proj%nx, proj%ny))
         allocate (vsf1(proj%nx, proj%ny))
         allocate (vsf2(proj%nx, proj%ny))
         allocate (fsstsum(proj%nx, proj%ny))
         allocate (fsst(proj%nx, proj%ny))
         fsstsum = 0.
         nfssttimes = 0
      end if

      ! initialize the i4time_valid1 and i4time_valid2, which are the two available
      ! times which bound the desired output time
      ta = 1
      i4time_valid1 = i4times_avail(ta)
      i4time_valid2 = i4times_avail(ta + 1)

      ! read in the data for these two time periods
      if (model_code(m) .gt. 1) then
         ! get the pressure level data arrays

         ! mslp
         call read_wfomodel_data(nfid, msl_id, proj, time_index(ta), &
                                 msl_ksfc, msl_ksfc, slp1, istatus)
         call read_wfomodel_data(nfid, msl_id, proj, time_index(ta + 1), &
                                 msl_ksfc, msl_ksfc, slp2, istatus)
         ! height
         call read_wfomodel_data(nfid, ht_id, proj, time_index(ta), &
                                 ht_kbotp, ht_ktopp, z3d1, istatus)
         call read_wfomodel_data(nfid, ht_id, proj, time_index(ta + 1), &
                                 ht_kbotp, ht_ktopp, z3d2, istatus)
         ! temperature
         call read_wfomodel_data(nfid, t_id, proj, time_index(ta), &
                                 t_kbotp, t_ktopp, t3d1, istatus)
         call read_wfomodel_data(nfid, t_id, proj, time_index(ta + 1), &
                                 t_kbotp, t_ktopp, t3d2, istatus)

         ! rh
         call read_wfomodel_data(nfid, rh_id, proj, time_index(ta), &
                                 rh_kbotp, rh_ktopp, rh3d1, istatus)
         call read_wfomodel_data(nfid, rh_id, proj, time_index(ta + 1), &
                                 rh_kbotp, rh_ktopp, rh3d2, istatus)

         ! u
         call read_wfomodel_data(nfid, u_id, proj, time_index(ta), &
                                 u_kbotp, u_ktopp, u3d1, istatus)
         call read_wfomodel_data(nfid, u_id, proj, time_index(ta + 1), &
                                 u_kbotp, u_ktopp, u3d2, istatus)

         ! v
         call read_wfomodel_data(nfid, v_id, proj, time_index(ta), &
                                 v_kbotp, v_ktopp, v3d1, istatus)
         call read_wfomodel_data(nfid, v_id, proj, time_index(ta + 1), &
                                 v_kbotp, v_ktopp, v3d2, istatus)
      end if
      if ((model_code(m) .eq. 1) .or. (model_code(m) .eq. 3)) then
         call read_wfomodel_data(nfid, t_id, proj, time_index(ta), &
                                 t_ksfc, t_ksfc, tsf1, istatus)
         call read_wfomodel_data(nfid, t_id, proj, time_index(ta + 1), &
                                 t_ksfc, t_ksfc, tsf2, istatus)
         call read_wfomodel_data(nfid, u_id, proj, time_index(ta), &
                                 u_ksfc, u_ksfc, usf1, istatus)
         call read_wfomodel_data(nfid, u_id, proj, time_index(ta + 1), &
                                 u_ksfc, u_ksfc, usf2, istatus)
         call read_wfomodel_data(nfid, v_id, proj, time_index(ta), &
                                 v_ksfc, v_ksfc, vsf1, istatus)
         call read_wfomodel_data(nfid, v_id, proj, time_index(ta + 1), &
                                 v_ksfc, v_ksfc, vsf2, istatus)
         call read_wfomodel_data(nfid, rh_id, proj, time_index(ta), &
                                 rh_ksfc, rh_ksfc, rhsf1, istatus)
         call read_wfomodel_data(nfid, rh_id, proj, time_index(ta + 1), &
                                 rh_ksfc, rh_ksfc, rhsf2, istatus)

      end if
      allocate (topo(proj%nx, proj%ny))
      call get_wfomodel_topo(nfid, proj, topo, istatus)

      ! main loop over all desired output times
      output_time_loop: do i4time_valid = i4time_cycle, i4time_last, outfreq_sec
         ! print *, 'getting data for i4time = ', i4time_valid
         if (i4time_valid .gt. i4time_valid2) then
            ! we need to get new bounding data times
            read_avail_loop: do while (i4time_valid .gt. i4time_valid2)
               ta = ta + 1
               i4time_valid1 = i4times_avail(ta)
               i4time_valid2 = i4times_avail(ta + 1)
            end do read_avail_loop
            ! read the data for these 2 times
            if (model_code(m) .gt. 1) then
               ! get the pressure level data arrays

               ! mslp
               call read_wfomodel_data(nfid, msl_id, proj, time_index(ta), &
                                       msl_ksfc, msl_ksfc, slp1, istatus)
               call read_wfomodel_data(nfid, msl_id, proj, time_index(ta + 1), &
                                       msl_ksfc, msl_ksfc, slp2, istatus)
               ! height
               call read_wfomodel_data(nfid, ht_id, proj, time_index(ta), &
                                       ht_kbotp, ht_ktopp, z3d1, istatus)
               call read_wfomodel_data(nfid, ht_id, proj, time_index(ta + 1), &
                                       ht_kbotp, ht_ktopp, z3d2, istatus)
               ! temperature
               call read_wfomodel_data(nfid, t_id, proj, time_index(ta), &
                                       t_kbotp, t_ktopp, t3d1, istatus)
               call read_wfomodel_data(nfid, t_id, proj, time_index(ta + 1), &
                                       t_kbotp, t_ktopp, t3d2, istatus)

               ! rh
               call read_wfomodel_data(nfid, rh_id, proj, time_index(ta), &
                                       rh_kbotp, rh_ktopp, rh3d1, istatus)
               call read_wfomodel_data(nfid, rh_id, proj, time_index(ta + 1), &
                                       rh_kbotp, rh_ktopp, rh3d2, istatus)

               ! u
               call read_wfomodel_data(nfid, u_id, proj, time_index(ta), &
                                       u_kbotp, u_ktopp, u3d1, istatus)
               call read_wfomodel_data(nfid, u_id, proj, time_index(ta + 1), &
                                       u_kbotp, u_ktopp, u3d2, istatus)

               ! v
               call read_wfomodel_data(nfid, v_id, proj, time_index(ta), &
                                       v_kbotp, v_ktopp, v3d1, istatus)
               call read_wfomodel_data(nfid, v_id, proj, time_index(ta + 1), &
                                       v_kbotp, v_ktopp, v3d2, istatus)
            end if
            if ((model_code(m) .eq. 1) .or. (model_code(m) .eq. 3)) then
               call read_wfomodel_data(nfid, t_id, proj, time_index(ta), &
                                       t_ksfc, t_ksfc, tsf1, istatus)
               call read_wfomodel_data(nfid, t_id, proj, time_index(ta + 1), &
                                       t_ksfc, t_ksfc, tsf2, istatus)
               call read_wfomodel_data(nfid, u_id, proj, time_index(ta), &
                                       u_ksfc, u_ksfc, usf1, istatus)
               call read_wfomodel_data(nfid, u_id, proj, time_index(ta + 1), &
                                       u_ksfc, u_ksfc, usf2, istatus)
               call read_wfomodel_data(nfid, v_id, proj, time_index(ta), &
                                       v_ksfc, v_ksfc, vsf1, istatus)
               call read_wfomodel_data(nfid, v_id, proj, time_index(ta + 1), &
                                       v_ksfc, v_ksfc, vsf2, istatus)
               call read_wfomodel_data(nfid, rh_id, proj, time_index(ta), &
                                       rh_ksfc, rh_ksfc, rhsf1, istatus)
               call read_wfomodel_data(nfid, rh_id, proj, time_index(ta + 1), &
                                       rh_ksfc, rh_ksfc, rhsf2, istatus)
            end if
         end if

         ! handle case where rh has fewer levels than
         ! the other state variables (only applies if
         ! we are obtaining upper air data (model_code >1)

         ! if we are processing upper air data from this model,
         ! then run through some cleanup to account for missing
         ! levels.
         if (model_code(m) .gt. 1) then

            ! handle case where rh has fewer levels than
            ! the other state variables (only applies if
            ! we are obtaining upper air data
            np_rh = np_t
            if (np_rh_raw .lt. np_t) then
               print *, 'warning: expanding rh'
               allocate (data3d_temp(proj%nx, proj%ny, np_rh))
               ! do the rh3d1 array
               data3d_temp(:, :, :) = rmissingval
               rh1_zloop: do k = 1, np_rh_raw
                  t1_zloop: do kk = 1, np_rh
                     if (rh_plevels_raw(k) .eq. t_plevels(kk)) then
                        rh_inv(kk + t_kbotp - 1, time_index(ta)) = &
                           rh_inv_raw(k + rh_kbotp - 1, time_index(ta))
                        if (rh_inv(kk + t_kbotp - 1, time_index(ta))) then
                           data3d_temp(:, :, kk) = rh3d1(:, :, k)
                        end if
                        exit t1_zloop
                     end if
                  end do t1_zloop
               end do rh1_zloop

               deallocate (rh3d1)
               allocate (rh3d1(proj%nx, proj%ny, np_rh))
               rh3d1(:, :, :) = data3d_temp(:, :, :)

               ! fill in top level if not already filled
               if (maxval(rh3d1(:, :, np_rh)) .eq. rmissingval) then
                  print *, 'setting top level rh to 5%'
                  rh3d1(:, :, np_rh) = 5.0  ! very dry
                  rh_inv(t_ktopp, time_index(ta)) = .true.
               end if
               ! fill in bottom level if not already filled (the
               ! mesoeta does not have the 1000mb value!)
               if (maxval(rh3d1(:, :, 1)) .eq. rmissingval) then
                  rh3d1(:, :, 1) = rh3d1(:, :, 2)
                  rh_inv(t_kbotp, time_index(ta)) = .true.
               end if

               ! repeat for the rh3d2 array
               data3d_temp(:, :, :) = rmissingval
               rh2_zloop: do k = 1, np_rh_raw
                  t2_zloop: do kk = 1, np_rh
                     if (rh_plevels_raw(k) .eq. t_plevels(kk)) then
                        rh_inv(kk + t_kbotp - 1, time_index(ta + 1)) = &
                           rh_inv_raw(k + rh_kbotp - 1, time_index(ta + 1))
                        if (rh_inv(kk + t_kbotp - 1, time_index(ta + 1))) then
                           data3d_temp(:, :, kk) = rh3d2(:, :, k)
                        end if
                        exit t2_zloop
                     end if
                  end do t2_zloop
               end do rh2_zloop
               deallocate (rh3d2)
               allocate (rh3d2(proj%nx, proj%ny, np_rh))
               rh3d2(:, :, :) = data3d_temp(:, :, :)

               ! fill in top level if not already filled
               if (maxval(rh3d2(:, :, np_rh)) .eq. rmissingval) then
                  rh3d2(:, :, np_rh) = 5.0  ! very dry
                  rh_inv(t_ktopp, time_index(ta + 1)) = .true.
               end if
               ! fill in bottom level if not already filled
               if (maxval(rh3d2(:, :, 1)) .eq. rmissingval) then
                  rh3d2(:, :, 1) = rh3d2(:, :, 2)
                  rh_inv(t_kbotp, time_index(ta + 1)) = .true.
               end if

               deallocate (data3d_temp)
               rh_plevels(1:np_rh) = t_plevels(1:np_t)
            else
               rh_plevels = rh_plevels_raw
            end if

            ! clean up all 3d arrays to fill in any missing levels

            ! z3d
            fixvar1 = .false.
            fixvar2 = .false.
            do k = 1, np_ht
               if (.not. ht_inv(k - 1 + ht_kbotp, time_index(ta))) then
                  z3d1(:, :, k) = rmissingval
                  fixvar1 = .true.
               end if
               if (.not. ht_inv(k - 1 + ht_kbotp, time_index(ta + 1))) then
                  z3d2(:, :, k) = rmissingval
                  fixvar2 = .true.
               end if
               if (fixvar1) then
                  print *, 'warning: filling in missing levels for z3d1'
                  ! call the fill routine to fill in for z3d1
                  call fill_missing_levs(proj%nx, proj%ny, np_ht, ht_plevels, &
                                         z3d1, rmissingval, 2)
               end if
               if (fixvar2) then
                  ! call the vinterp routine to fill in for z3d2
                  print *, 'warning: filling in missing levels for z3d2'
                  call fill_missing_levs(proj%nx, proj%ny, np_ht, ht_plevels, &
                                         z3d2, rmissingval, 2)
               end if
            end do

            ! t3d
            fixvar1 = .false.
            fixvar2 = .false.
            do k = 1, np_t
               if (.not. t_inv(k - 1 + t_kbotp, time_index(ta))) then
                  t3d1(:, :, k) = rmissingval
                  fixvar1 = .true.
               end if
               if (.not. t_inv(k - 1 + t_kbotp, time_index(ta + 1))) then
                  t3d2(:, :, k) = rmissingval
                  fixvar2 = .true.
               end if
               if (fixvar1) then
                  ! call the vinterp routine to fill in for t3d1
                  print *, 'warning: filling in missing levels for t3d1'
                  call fill_missing_levs(proj%nx, proj%ny, np_t, t_plevels, &
                                         t3d1, rmissingval, 2)
               end if
               if (fixvar2) then
                  ! call the vinterp routine to fill in for t3d2
                  print *, 'warning: filling in missing levels for t3d2'
                  call fill_missing_levs(proj%nx, proj%ny, np_t, t_plevels, &
                                         t3d2, rmissingval, 2)
               end if
            end do

            ! u3d
            fixvar1 = .false.
            fixvar2 = .false.
            do k = 1, np_u
               if (.not. u_inv(k - 1 + u_kbotp, time_index(ta))) then
                  u3d1(:, :, k) = rmissingval
                  fixvar1 = .true.
               end if
               if (.not. u_inv(k - 1 + u_kbotp, time_index(ta + 1))) then
                  u3d2(:, :, k) = rmissingval
                  fixvar2 = .true.
               end if
               if (fixvar1) then
                  ! call the vinterp routine to fill in for u3d1
                  print *, 'warning: filling in missing levels for u3d1'
                  call fill_missing_levs(proj%nx, proj%ny, np_u, u_plevels, &
                                         u3d1, rmissingval, 1)
               end if
               if (fixvar2) then
                  ! call the vinterp routine to fill in for u3d2
                  print *, 'warning: filling in missing levels for u3d2'
                  call fill_missing_levs(proj%nx, proj%ny, np_u, u_plevels, &
                                         u3d2, rmissingval, 1)
               end if
            end do

            ! v3d
            fixvar1 = .false.
            fixvar2 = .false.
            do k = 1, np_v
               if (.not. v_inv(k - 1 + v_kbotp, time_index(ta))) then
                  v3d1(:, :, k) = rmissingval
                  fixvar1 = .true.
               end if
               if (.not. v_inv(k - 1 + v_kbotp, time_index(ta + 1))) then
                  v3d2(:, :, k) = rmissingval
                  fixvar2 = .true.
               end if
               if (fixvar1) then
                  ! call the vinterp routine to fill in for v3d1
                  print *, 'warning: filling in missing levels for v3d1'
                  call fill_missing_levs(proj%nx, proj%ny, np_v, v_plevels, &
                                         v3d1, rmissingval, 1)
               end if
               if (fixvar2) then
                  ! call the vinterp routine to fill in for v3d2
                  print *, 'warning: filling in missing levels for z3d2'
                  call fill_missing_levs(proj%nx, proj%ny, np_v, v_plevels, &
                                         v3d2, rmissingval, 1)
               end if
            end do

            ! rh3d
            fixvar1 = .false.
            fixvar2 = .false.
            do k = 1, np_rh
               if (.not. rh_inv(k - 1 + t_kbotp, time_index(ta))) then
                  rh3d1(:, :, k) = rmissingval
                  fixvar1 = .true.
               end if
               if (.not. rh_inv(k - 1 + t_kbotp, time_index(ta + 1))) then
                  rh3d2(:, :, k) = rmissingval
                  fixvar2 = .true.
               end if
               if (fixvar1) then
                  ! call the vinterp routine to fill in for rh3d1
                  print *, 'warning: filling in missing levels for rh3d1'
                  call fill_missing_levs(proj%nx, proj%ny, np_rh, rh_plevels, &
                                         rh3d1, rmissingval, 1)
               end if
               if (fixvar2) then
                  ! call the vinterp routine to fill in for rh3d2
                  print *, 'warning: filling in missing levels for rh3d2'
                  call fill_missing_levs(proj%nx, proj%ny, np_rh, rh_plevels, &
                                         rh3d2, rmissingval, 1)
               end if
            end do
         end if

         ! time to do time interpolation

         ! print *, 'using bounding i4times of ', i4time_valid1,i4time_valid2
         ! at this point, we have all of the data we need for two
         ! bounding times.  so interpolate to the desired time,
         ! then output the data
         if (i4time_valid .eq. i4time_valid1) then
            if (model_code(m) .gt. 1) then
               z3d = z3d1
               t3d = t3d1
               u3d = u3d1
               v3d = v3d1
               rh3d = rh3d1
               slp = slp1
            end if
            if ((model_code(m) .eq. 1) .or. (model_code(m) .eq. 3)) then
               tsf = tsf1
               rhsf = rhsf1
               usf = usf1
               vsf = vsf1
            end if

         else if (i4time_valid .eq. i4time_valid2) then
            if (model_code(m) .gt. 1) then
               z3d = z3d2
               t3d = t3d2
               u3d = u3d2
               v3d = v3d2
               rh3d = rh3d2
               slp = slp2
            end if
            if ((model_code(m) .eq. 1) .or. (model_code(m) .eq. 3)) then
               tsf = tsf2
               rhsf = rhsf2
               usf = usf2
               vsf = vsf2
            end if

         else
            ! time interpolate
            weight1 = float(i4time_valid2 - i4time_valid)/float(i4time_valid2 - i4time_valid1)
            if (model_code(m) .gt. 1) then
               z3d = weight1*z3d1 + (1.-weight1)*z3d2
               t3d = weight1*t3d1 + (1.-weight1)*t3d2
               u3d = weight1*u3d1 + (1.-weight1)*u3d2
               v3d = weight1*v3d1 + (1.-weight1)*v3d2
               rh3d = weight1*rh3d1 + (1.-weight1)*rh3d2
               slp = weight1*slp1 + (1.-weight1)*slp2
            end if
            if ((model_code(m) .eq. 1) .or. (model_code(m) .eq. 3)) then
               tsf = weight1*tsf1 + (1.-weight1)*tsf2
               rhsf = weight1*rhsf1 + (1.-weight1)*rhsf2
               usf = weight1*usf1 + (1.-weight1)*usf2
               vsf = weight1*vsf1 + (1.-weight1)*vsf2
            end if

         end if

         ! output the data for this time
         if (output_type(1:3) .eq. 'mm5') then
            if (model_code(m) .gt. 1) then
               print *, 'writing mm5v3 3d state variables for i4time =', i4time_valid
               call output_mm5v3_basic(i4time_cycle, i4time_valid, proj, &
                                       np_ht, np_t, np_u, np_v, np_rh, &
                                       ht_plevels, t_plevels, u_plevels, v_plevels, rh_plevels, &
                                       z3d, t3d, u3d, v3d, rh3d, slp, topo, &
                                       ext_data_path, output_name(m), mm5mode_new, istatus)
               if (model_code(m) .eq. 3) then
                  print *, 'writing mm5v3 2d state variables for i4time =', i4time_valid
                  call output_mm5v3_sfc(i4time_cycle, i4time_valid, proj, &
                                        tsf, usf, vsf, rhsf, &
                                        ext_data_path, output_name(m), mm5mode_append, istatus)
               end if
            else
               print *, 'writing mm5v3 2d state variables for i4time =', i4time_valid
               call output_mm5v3_sfc(i4time_cycle, i4time_valid, proj, &
                                     tsf, usf, vsf, rhsf, &
                                     ext_data_path, output_name(m), mm5mode_new, istatus)
            end if
         else if (output_type(1:3) .eq. 'wrf') then
            if (model_code(m) .gt. 1) then
               print *, 'writing wrf 3d state variables for i4time =', i4time_valid
               call output_wrf_basic(i4time_cycle, i4time_valid, proj, &
                                     np_ht, np_t, np_u, np_v, np_rh, &
                                     ht_plevels, t_plevels, u_plevels, v_plevels, rh_plevels, &
                                     z3d, t3d, u3d, v3d, rh3d, slp, topo, &
                                     ext_data_path, output_name(m), wrfmode_new, istatus)
               if (model_code(m) .eq. 3) then
                  print *, 'writing wrf 2d state variables for i4time =', i4time_valid
                  call output_wrf_sfc(i4time_cycle, i4time_valid, proj, &
                                      tsf, usf, vsf, rhsf, &
                                      ext_data_path, output_name(m), wrfmode_append, istatus)
               end if
            else
               print *, 'writing wrf 2d state variables for i4time =', i4time_valid
               call output_wrf_sfc(i4time_cycle, i4time_valid, proj, &
                                   tsf, usf, vsf, rhsf, &
                                   ext_data_path, output_name(m), wrfmode_new, istatus)
            end if
         else if (output_type(1:4) .eq. 'rams') then
            print *, 'rams output coming soon'
            stop
         else
            print *, 'unrecognized output format requested: ', output_type
            stop
         end if

         ! if we have surface data, compute a fake sst field
         ! from the wet bulb temperature for this time.  however,
         ! we want to average wet bulb t over all time periods and write
         ! one sst field per run
         if ((model_code(m) .eq. 1) .or. (model_code(m) .eq. 3)) then
            fsst(:, :) = 300. ! initialize to something not too bad
            do jy = 1, proj%ny
               do ix = 1, proj%nx
                  ! compute tdsf (deg c)
                  tdsf_c = dwpt_laps(tsf(ix, jy) - 273.15, rhsf(ix, jy))

                  ! estimate psf from topo
                  psf_mb = 1013.25*(1.-(topo(ix, jy)*2.257e-5))**(5.259)

                  ! compute fsst (wet bulb t) from tsf/rhsf/psf
                  fsst(ix, jy) = twet_fast(tsf(ix, jy) - 273.15, tdsf_c, psf_mb) + 273.15
               end do
            end do
            ! add fsst to fsstsum and increment nfssttimes
            fsstsum(:, :) = fsstsum(:, :) + fsst(:, :)
            nfssttimes = nfssttimes + 1
         end if
      end do output_time_loop

      ! if nfssttimes > 0, compute mean fsst and output
      if (nfssttimes .gt. 0) then

         fsst(:, :) = fsstsum(:, :)/float(nfssttimes)
         print *, 'min/max fake sst: ', minval(fsst), maxval(fsst)

         ! call output routine for specific model
         if (output_type(1:3) .eq. 'mm5') then
            call output_mm5v3_fsst(i4time_cycle, i4time_cycle, proj, &
                                   fsst, ext_data_path, output_name(m), mm5mode_new, istatus)

         end if
      end if
      ! deallocate 3 arrays for each state variable.  two will be needed for the
      ! data at two times bounding the time of interest.  the third is where
      ! the time interpolated values will be saved.

      if (model_code(m) .gt. 1) then
         deallocate (z3d)
         deallocate (z3d1)
         deallocate (z3d2)
         deallocate (t3d)
         deallocate (t3d1)
         deallocate (t3d2)
         deallocate (rh3d)
         deallocate (rh3d1)
         deallocate (rh3d2)
         deallocate (u3d)
         deallocate (u3d1)
         deallocate (u3d2)
         deallocate (v3d)
         deallocate (v3d1)
         deallocate (v3d2)
         deallocate (slp)
         deallocate (slp1)
         deallocate (slp2)
      end if
      if ((model_code(m) .eq. 1) .or. (model_code(m) .eq. 3)) then
         deallocate (tsf)
         deallocate (tsf1)
         deallocate (tsf2)
         deallocate (rhsf)
         deallocate (rhsf1)
         deallocate (rhsf2)
         deallocate (usf)
         deallocate (usf1)
         deallocate (usf2)
         deallocate (vsf)
         deallocate (vsf1)
         deallocate (vsf2)
         deallocate (fsst)
         deallocate (fsstsum)
      end if
      ! show this cycle for this model as having been processed.
      print '(2a)', '    updating process log:', trim(proclog)
      open (file=proclog, unit=11, form='formatted', &
            access='sequential', status='replace')
      write (11, '(a13)') latest_cycle_wfo
      close (11)

      ! deallocate arrays specific to this model
      if (allocated(topo)) deallocate (topo)
      if (allocated(i4times_avail)) deallocate (i4times_avail)
      if (allocated(i4times_avail_max)) deallocate (i4times_avail_max)
      call close_wfofile(nfid, istatus)
   end do model_loop

end program wfoprep

