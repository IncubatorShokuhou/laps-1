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

   !  variables for file names and laps_data_root

   character(len=256) :: laps_data_root, &
                         input_laps_file
   character(len=9)   :: laps_file_time
   integer             :: valid_yyyy, valid_jjj, &
                          valid_hh, valid_min, &
                          i4time
   real, parameter     :: missingflag = 1.e37
   ! namelist items

   logical            :: hotstart, balance, make_sfc_uv, use_sfc_bal, use_laps_skintemp, use_laps_vv
   logical            :: lapsprep_always_write
   character(len=4)  :: output_format(10)
   integer            :: num_output
   integer            :: num_soil_layers
   real               :: soil_layer_depths(10)
   real               :: snow_thresh, lwc2vapor_thresh, ice2vapor_thresh
   real               :: hydrometeor_scale_factor_pcp, hydrometeor_scale_factor_cld
   real               :: rai_frac, sno_frac
   character(len=256):: output_prefix

   !  output file info.

   integer, parameter :: output = 10
   character(len=132) :: name

   !  the items below define the number of different laps files that will
   !  be read, what their extensions are, and what variables from
   !  each file.  these are hard-coded for now, but we may put this
   !  into lapsprep.nl in the future.  mandatory extensions must be
   !  listed first

   ! yuanfu added lm1 for land soil input file:
   integer, parameter :: num_ext = 10
   character(len=3), dimension(num_ext) :: ext = (/'lt1', 'lw3', &
                                                   'lh3', 'lsx', &
                                                   'lsx', 'lq3', &
                                                   'lwc', 'lm1', &
                                                   'lm2', 'lcp'/)

   ! yuanfu added the 8th row here for lm1's lsm variable:
   character(len=3), dimension(5, num_ext) :: cdf_var_name = reshape( &
                                              (/'ht ', 't3 ', 'xxx', 'xxx', 'xxx', &
                                                'u3 ', 'v3 ', 'om ', 'xxx', 'xxx', &
                                                'rhl', 'xxx', 'xxx', 'xxx', 'xxx', &
                                                'u  ', 'v  ', 't  ', 'rh ', 'tgd', &
                                                'ps ', 'msl', 'mr ', 'vv ', 'xxx', &
                                                'sh ', 'xxx', 'xxx', 'xxx', 'xxx', &
                                                'lwc', 'rai', 'sno', 'pic', 'ice', &
                                                'lsm', 'xxx', 'xxx', 'xxx', 'xxx', &
                                                'sc ', 'xxx', 'xxx', 'xxx', 'xxx', &
                                                'lcp', 'xxx', 'xxx', 'xxx', 'xxx'/), &
                                              (/5, num_ext/))

   ! the third from the last is the number of vars of lm1:
   ! hardcoded now only reading in the soil moisture by yuanfu xie
   integer, dimension(num_ext) :: num_cdf_var = (/2, 3, 1, 5, 4, 1, 5, 1, 1, 1/)

contains

   subroutine read_namelist

      implicit none

      integer :: ioerror, nml_unit, c

      namelist /lapsprep_nl/ hotstart, &
         balance, &
         use_sfc_bal, &
         hydrometeor_scale_factor_pcp, &
         hydrometeor_scale_factor_cld, &
         make_sfc_uv, &
         output_format, &
         snow_thresh, &
         lwc2vapor_thresh, &
         ice2vapor_thresh, &
         rai_frac, &
         sno_frac, &
         use_laps_skintemp, &
         use_laps_vv, &
         lapsprep_always_write, &
         num_soil_layers, &
         soil_layer_depths

      nml_unit = 77

      output_format(:) = '    '
      ! set namelist defaults
      balance = .false.
      hydrometeor_scale_factor_pcp = 0.0
      hydrometeor_scale_factor_cld = 0.5
      snow_thresh = 0.5
      lwc2vapor_thresh = 0.
      ice2vapor_thresh = 0.
      make_sfc_uv = .false.
      rai_frac = 1.0
      sno_frac = 1.0
      use_sfc_bal = .false.
      use_laps_skintemp = .false.
      use_laps_vv = .false.
      lapsprep_always_write = .false.
      ! open the namelist

      input_laps_file = trim(laps_data_root)//'/static/lapsprep.nl'
      open (file=trim(input_laps_file), &
            unit=nml_unit, &
            status='old', &
            form='formatted', &
            iostat=ioerror)

      if (ioerror .ne. 0) then
         print '(3a,i4)', 'error opening ', trim(input_laps_file), ' code #', ioerror
         stop 'error_opening_namelist'
      end if

      ! read and print the namelist

      read (nml_unit, nml=lapsprep_nl)
      print '(a)', 'running lapsprep using the following settings:'
      write (6, nml=lapsprep_nl)

      ! close the namelist

      close (nml_unit)

      ! determine number of output formats
      num_output = 0
      find_outputs: do c = 1, 10
         if (output_format(c) .eq. '    ') then
            num_output = c - 1
            exit find_outputs
         end if
      end do find_outputs
      if (num_output .le. 0) then
         print '(a)', 'must specify at least one output format!'
         stop 'read_namelist'
      end if

   end subroutine read_namelist

end module setup
