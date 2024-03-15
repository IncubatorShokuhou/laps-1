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

module wfoprep_setup

   ! this module contains data and subroutines related to the setup of
   ! the wfoprep program.
   !
   ! history
   ! -------
   ! sep 2001 - original version..............b. shaw, noaa/cira/fsl

   implicit none

   ! set default declarations to private
   private

   ! private variables available only within module
   integer, parameter           :: max_sources = 10
   ! public variables shared within and outside of this module

   integer, public              :: num_models
   integer, public              :: num_formats
   integer, public              :: model_code(max_sources)
   integer, public              :: model_run_freq(max_sources)
   integer, public              :: model_delay(max_sources)
   integer, public              :: max_fcst_len(max_sources)
   integer, public              :: output_freq(max_sources)
   real, public                 :: min_time_frac
   real, public                 :: min_vert_frac
   character(len=256)          :: fxa_data
   character(len=132), public   :: model_name(max_sources)
   character(len=256), public   :: model_path(max_sources)
   character(len=10), public   :: output_type
   character(len=32), public   :: output_name(max_sources)
   character(len=256), public   :: laps_data_root
   character(len=256), public   :: ext_data_path

   ! time stuff
   integer, public              :: i4time_now
   integer, public              :: hour_now
   integer, public              :: year4_now
   integer, public              :: year2_now
   integer, public              :: day3_now
   character(len=9), public      :: a9time_now

   ! routines publicly availble to all routines using this module
   public setup_wfoprep

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine setup_wfoprep(istatus)
      ! calls all of the wfoprep initialization routines

      implicit none
      integer :: istatus

      istatus = 1

      call get_data_roots(istatus)
      if (istatus .ne. 1) then
         print '(a)', 'setup_dprep: problem getting data roots.'
         istatus = 0
         return
      end if

      call read_wfoprep_nl(istatus)
      if (istatus .ne. 1) then
         print '(a)', 'setup_dprep: problem processing namelist.'
         istatus = 0
         return
      end if

      call get_systime(i4time_now, a9time_now, istatus)
      if (istatus .ne. 1) then
         print '(a)', 'setup_dprep: problem getting systime.'
         istatus = 0
         return
      else
         print '(a,i12,1x,a9)', 'laps i4time/a9time = ', i4time_now, a9time_now
      end if
      read (a9time_now(1:2), '(i2)') year2_now
      if (year2_now .lt. 90) then
         year4_now = year2_now + 2000
      else
         year4_now = year2_now + 1900
      end if
      read (a9time_now(3:5), '(i3)') day3_now
      read (a9time_now(6:7), '(i2)') hour_now
      return

   end subroutine setup_wfoprep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_data_roots(istatus)

      ! sets up the laps_data_root and ext_data_root.
      implicit none
      integer                    :: istatus
      character(len=256)         :: ext_data_root

      print '(a)', 'setting up data roots for wfoprep.'

      call getenv("laps_data_root", laps_data_root)
      if (laps_data_root(1:10) .eq. '          ') then
         print '(a)', 'laps_data_root not set.'
         istatus = 0
         return
      end if

      ! see if ext_data_root is set.  this allows the flexibility for
      ! users of the mm5-laps distribution to share an extdataroot between
      ! pregrid and wfoprep.  if not set, then set it to the appropriate
      ! laps directory.
      call getenv("ext_data_root", ext_data_root)
      if (ext_data_root(1:10) .eq. '          ') then
         ext_data_path = trim(laps_data_root)//'/lapsprd/dprep/'
      else
         ext_data_path = trim(ext_data_root)//'/extprd/'
      end if
      print '(2a)', 'sbnprep output will be written in ', trim(ext_data_path)

      print '(a)', 'successfully set up data roots.'
      istatus = 1
      return

   end subroutine get_data_roots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_wfoprep_nl(istatus)

      ! subroutine to read the sbnprep.nl namelist file

      implicit none

      integer, parameter          :: nlunit = 22
      integer                     :: i
      integer                     :: istatus
      integer                     :: ioerr
      character(len=256)          :: nlfile
      logical                     :: nlfound
      logical                     :: done_proc
      character(len=132)          :: model_name_tmp
      namelist /wfoprep_nl/ model_code, model_name, &
         model_run_freq, model_delay, max_fcst_len, &
         output_freq, output_type, fxa_data, &
         output_name, min_time_frac, min_vert_frac

      print '(a)', 'processing wfoprep.nl namelist information.'
      ! initialize the namelist values
      model_code(:) = -1
      model_name(:) = '                                                '
      model_path(:) = '                                                '
      model_run_freq(:) = 0
      model_delay(:) = 0
      max_fcst_len(:) = 0
      output_freq(:) = 0
      output_type = '          '
      min_time_frac = 0.8
      min_vert_frac = 0.5
      ! build the file name of the namelist and make sure it exists
      nlfile = trim(laps_data_root)//'/static/wfoprep.nl'
      inquire (file=nlfile, exist=nlfound)
      if (.not. nlfound) then
         print '(2a)', 'namelist file not found: ', trim(nlfile)
         istatus = 0
         return
      end if

      ! open the file, rewind for safety, then read the namelist data.
      open (file=nlfile, unit=nlunit, form='formatted', status='old', &
            iostat=ioerr)
      if (ioerr .ne. 0) then
         print '(a,i4)', 'error opening namelist file: ', ioerr
         istatus = 0
         return
      end if

      rewind (nlunit)
      read (nlunit, nml=wfoprep_nl, iostat=ioerr)
      if (ioerr .ne. 0) then
         print '(a,i4)', 'i/o error reading namelist: ', ioerr
         istatus = 0
         close (nlunit)
         return
      else
         close (nlunit)
      end if

      ! process the data
      num_models = 1
      done_proc = .false.

      ! get fxa_data path
      if (fxa_data(1:10) .eq. '          ') then
         call getenv('fxa_data', fxa_data)
         if (fxa_data(1:10) .eq. '          ') then
            print *, 'error, no fxa_data path set in namelist or environment.'
            stop
         end if
      end if

      do while ((.not. done_proc) .and. (num_models .lt. max_sources))

         if (model_code(num_models) .gt. 0) then
            model_name_tmp = model_name(num_models)
            if (model_name_tmp(1:5) .ne. '     ') then
               if (model_run_freq(num_models) .gt. 0) then
                  print '(a,i3)', 'will process code', model_code(num_models)
                  ! build path for this model
                  model_path(num_models) = trim(fxa_data)//'/grid/sbn/netcdf/' &
                                           //trim(model_name(num_models))
               else
                  print '(2a)', 'no model_run_freq set for ', &
                     trim(model_name(num_models))
                  done_proc = .true.
               end if
            else
               print '(a,i2)', 'no name provided for model code: ', &
                  model_code(num_models)
               done_proc = .true.
            end if
         else
            print '(a)', 'no more valid models.'
            done_proc = .true.
         end if
         if (.not. done_proc) then
            num_models = num_models + 1
         else
            num_models = num_models - 1
         end if
      end do
      if (num_models .lt. 1) then
         print '(a)', 'no valid model sources have been specified!'
         print '(a)', 'check your wfoprep.nl settings.'
         istatus = 0
         return
      else
         print '(a)', 'wfoprep will process the following sources...'
         print '(a)', 'code   name                   path'
         print '(a)', '----   --------------------   --------------------------------------'
         do i = 1, num_models
            print '(i4,3x,a20,3x,a)', model_code(i), model_name(i) (1:20), &
               trim(model_path(i))
         end do
      end if

      ! check output format for validity
      if ((output_type .ne. 'mm5       ') .and. &
          (output_type .ne. 'rams      ') .and. &
          (output_type .ne. 'wrf       ')) then
         print '(2a)', 'invalid output format requested: ', &
            output_type
         print *, 'must be one of:  mm5, rams, wrf'
         istatus = 0
         return
      else
         print '(a)', 'wfoprep will output the following format:', output_type
      end if

      print '(a)', 'successfully processed namelist.'

      istatus = 1
      return
   end subroutine read_wfoprep_nl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module wfoprep_setup
