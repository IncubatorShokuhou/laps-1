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

  subroutine output_laps_format(ht, u3, v3, w3, om, t3, sh, rh3, lwc, ice, &
                                rai, sno, pic, tke, ref, pty, u, v, w, t, td, &
                                rh, lcb, lct, msl, p, ps, lil, tpw, r01, rto, &
                                s01, sto, th, the, pbe, nbe, lcv, cce, lmt, &
                                lmr, llr, spt, lhe, li, hi, vis, terdot, &
                                lwout, swout, shflux, lhflux, pblhgt, ground_t, &
                                upb, vpb, vnt, ham, hah, fwi, &
                                press_levels, &
                                lfmprd_dir, laps_data_root, domnum, &
                                laps_reftime, laps_valtime, nx, ny, nz, &
                                realtime, write_to_lapsdir, model_name, &
                                make_donefile)

     ! creates the laps *.fua and *.fsf files for model output.  the names of the
     ! input variables correspond to their netcdf names found in the fua.cdl
     ! and fsf.cdl files.
     !
     ! use of this routine requires a valid mm5_data_root that contains
     ! a subdirectory in mm5_data_root/mm5prd/dxx/fua (and fsf), where the
     ! xx is the two digit domain number.  the fsf.cdl and fua.cdl files
     ! in laps_data_root/cdl should also have dimensions equal to nx/ny.
     !
     ! basically, this routine assumes the model output grids are identical
     ! (projection, resolution, dimensions, etc.) to the domain defined
     ! in laps_data_root/static/nest7grid.parms.
     !
     !
     !
     ! history
     ! =======
     ! initial version:  brent shaw, noaa/fsl 5 jan 01
     ! modified to output files to domain specific directories in
     ! mm5_data_root instead of laps_data_root.  2 nov 01

     ! note that the only variable that is supported in the netcdf files but
     ! not produced by this routine is the fire index.

     implicit none

     ! argument declarations
     integer, intent(in)             :: nx
     integer, intent(in)             :: ny
     integer, intent(in)             :: nz
     real, intent(in)                :: ht(nx, ny, nz)
     real, intent(in)                :: u3(nx, ny, nz)
     real, intent(in)                :: v3(nx, ny, nz)
     real, intent(in)                :: w3(nx, ny, nz)
     real, intent(in)                :: om(nx, ny, nz)
     real, intent(in)                :: t3(nx, ny, nz)
     real, intent(in)                :: sh(nx, ny, nz)
     real, intent(in)                :: rh3(nx, ny, nz)
     real, intent(in)                :: lwc(nx, ny, nz)
     real, intent(in)                :: ice(nx, ny, nz)
     real, intent(in)                :: rai(nx, ny, nz)
     real, intent(in)                :: sno(nx, ny, nz)
     real, intent(in)                :: pic(nx, ny, nz)
     real, intent(in)                :: tke(nx, ny, nz)
     real, intent(in)                :: ref(nx, ny, nz)
     real, intent(in)                :: pty(nx, ny, nz)
     real, intent(in)                :: u(nx, ny)
     real, intent(in)                :: v(nx, ny)
     real, intent(in)                :: w(nx, ny)
     real, intent(in)                :: t(nx, ny)
     real, intent(in)                :: td(nx, ny)
     real, intent(in)                :: rh(nx, ny)
     real, intent(in)                :: lcb(nx, ny)
     real, intent(in)                :: lct(nx, ny)
     real, intent(in)                :: msl(nx, ny)
     real, intent(in)                :: p(nx, ny)
     real, intent(in)                :: ps(nx, ny)
     real, intent(in)                :: lil(nx, ny)
     real, intent(in)                :: tpw(nx, ny)
     real, intent(in)                :: r01(nx, ny)
     real, intent(in)                :: rto(nx, ny)
     real, intent(in)                :: s01(nx, ny)
     real, intent(in)                :: sto(nx, ny)
     real, intent(in)                :: th(nx, ny)
     real, intent(in)                :: the(nx, ny)
     real, intent(in)                :: pbe(nx, ny)
     real, intent(in)                :: nbe(nx, ny)
     real, intent(in)                :: lcv(nx, ny)
     real, intent(in)                :: cce(nx, ny)
     real, intent(in)                :: lmt(nx, ny)
     real, intent(in)                :: lmr(nx, ny)
     real, intent(in)                :: llr(nx, ny)
     real, intent(in)                :: spt(nx, ny)
     real, intent(in)                :: lhe(nx, ny)
     real, intent(in)                :: li(nx, ny)
     real, intent(in)                :: hi(nx, ny)
     real, intent(in)                :: vis(nx, ny)
     real, intent(in)                :: terdot(nx, ny)
     real, intent(in)                :: lwout(nx, ny)
     real, intent(in)                :: swout(nx, ny)
     real, intent(in)                :: shflux(nx, ny)
     real, intent(in)                :: lhflux(nx, ny)
     real, intent(in)                :: pblhgt(nx, ny)
     real, intent(in)                :: ground_t(nx, ny)
     real, intent(in)                :: upb(nx, ny)
     real, intent(in)                :: vpb(nx, ny)
     real, intent(in)                :: vnt(nx, ny)
     real, intent(in)                :: ham(nx, ny)
     real, intent(in)                :: hah(nx, ny)
     real, intent(in)                :: fwi(nx, ny)

     real, intent(in)                :: press_levels(nz)
     character(len=*), intent(in)     :: lfmprd_dir
     character(len=*), intent(in)     :: laps_data_root
     integer, intent(in)              :: domnum
     integer, intent(in)             :: laps_reftime
     integer, intent(in)             :: laps_valtime
     logical, intent(in)             :: realtime
     logical, intent(in)             :: write_to_lapsdir
     character(len=32), intent(in)   :: model_name
     logical, intent(in)             :: make_donefile
     ! locals
     character(len=2)             :: domnum_str
     integer, parameter           :: nvar3d = 16 ! equals # of 3d arrays above!
     integer, parameter           :: nvar2d = 44 ! # of 2d arrays above!
     real, allocatable               :: laps_data(:, :, :)
     integer, allocatable            :: levels(:)
     character(len=3), allocatable    :: varname(:)
     character(len=10), allocatable   :: varunits(:)
     character(len=4), allocatable   :: varlvltype(:)
     character(len=132), allocatable  :: varcomment(:)
     integer                         :: startind
     integer                         :: stopind
     integer                         :: istatus
     character(len=255)              :: output_dir
     character(len=255)              :: cdl_dir
     real, allocatable               :: cdl_levels(:)

     integer                         :: nz_one
     real                            :: cdl_levels_one
     character(len=9)                :: asctime, gtime
     character(len=4)                :: fcst_hhmm
     character(len=256)              :: output_file
     character(len=256)              :: donefile
     integer                         :: fnlen, extlen
     integer                         :: doneunit

     write (domnum_str, '(i2.2)') domnum
     extlen = 3
     nz_one = 1
     cdl_levels_one = 0.
     allocate (cdl_levels(nz))
     cdl_levels = press_levels(nz:1:-1)
     ! lets make the fua file first (contains 3d variables).  first, allocate
     ! the big array and all of the metadata arrays

     print *, 'outputting laps format (fua/fsf) files...'
     allocate (laps_data(nx, ny, nz*nvar3d))
     allocate (varname(nvar3d*nz))
     allocate (varunits(nvar3d*nz))
     allocate (varlvltype(nvar3d*nz))
     allocate (varcomment(nvar3d*nz))
     allocate (levels(nvar3d*nz))

     ! initialize varcomment
     varcomment(:) = '                                                  '// &
                     '                                                  '// &
                     '                                '
     startind = 1
     stopind = nz

     ! heights in meters
     laps_data(:, :, startind:stopind) = ht
     varname(startind:stopind) = 'ht '
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'geopotential height                   '

     ! u-wind in m/s
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = u3
     varname(startind:stopind) = 'u3 '
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'u-component wind                '

     ! v-wind in m/s
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = v3
     varname(startind:stopind) = 'v3 '
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'v-component wind                '

     ! w in m/s
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = w3
     varname(startind:stopind) = 'w3 '
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'vertical velocity               '

     ! omega in pa/s
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = om
     varname(startind:stopind) = 'om '
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'pressure vertical velocity       '

     ! temperature in k
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = t3
     varname(startind:stopind) = 't3 '
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'temperature              '

     ! specific humidity in kg/kg
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = sh
     varname(startind:stopind) = 'sh '
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'specific humidity                  '

     ! relative humidity wrt liquid in %
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = rh3
     varname(startind:stopind) = 'rh3'
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'relative humidity                         '

     ! cloud liquid content in kg/m2
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = lwc
     varname(startind:stopind) = 'lwc'
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'cloud liquid water                    '

     ! cloud ice content in kg/m2
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = ice
     varname(startind:stopind) = 'ice'
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'cloud ice                          '

     ! rain water content in kg/m2
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = rai
     varname(startind:stopind) = 'rai'
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'rain concentration            '

     ! snow content in kg/m2
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = sno
     varname(startind:stopind) = 'sno'
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'snow concentration            '

     ! graupel content in kg/m2
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = pic
     varname(startind:stopind) = 'pic'
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'graupel concentration            '

     ! reflectivity in dbz
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = ref
     varname(startind:stopind) = 'ref'
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'sim. radar reflectivity           '

     ! coded precipitation type
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = pty
     varname(startind:stopind) = 'pty'
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'precip. type                        '

     ! turbulent kinetic energy
     startind = stopind + 1
     stopind = startind + nz - 1
     laps_data(:, :, startind:stopind) = tke
     varname(startind:stopind) = 'tke'
     levels(startind:stopind) = nint(press_levels)
     varcomment(startind:stopind) = 'turbulent kinetic energy            '

     ! write out the 3d stuff using laps library routine
     if (.not. write_to_lapsdir) then
        output_dir = trim(lfmprd_dir)//'/d'//domnum_str//'/fua/'
     else
        output_dir = trim(laps_data_root)//'/lapsprd/fua/'// &
                     trim(model_name)//'/'
     end if
     cdl_dir = trim(laps_data_root)//'/cdl/'

     ! build the output file name so we can create a "donefile" if
     ! running in realtime mode

     call make_fnam_lp(laps_reftime, gtime, istatus)
     call make_fcst_time(laps_valtime, laps_reftime, fcst_hhmm, istatus)
     call cv_i4tim_asc_lp(laps_valtime, asctime, istatus)
     call cvt_fname_v3(output_dir, gtime, fcst_hhmm, 'fua', extlen, output_file, &
                       fnlen, istatus)

     print *, 'writing 3d fields to ', trim(output_file)
     call write_laps_lfm(laps_reftime, laps_valtime, output_dir, cdl_dir, &
                         'fua', &
                         nx, ny, nz*nvar3d, nz*nvar3d, varname, levels, &
                         varlvltype, varunits, varcomment, nz, cdl_levels, &
                         laps_data, istatus)
     if (istatus .ne. 1) then
        print *, 'error writing laps 3d (fua) file.'
     end if
     print *, 'done writing 3d data.'
     if ((realtime) .and. (istatus .eq. 1) .and. (make_donefile)) then
        donefile = trim(output_file)//'.done'
        call get_file_unit(doneunit)
        open (unit=doneunit, file=donefile, status='unknown')
        close (doneunit)
     end if
     deallocate (laps_data)
     deallocate (varname)
     deallocate (levels)
     deallocate (varlvltype)
     deallocate (varunits)
     deallocate (varcomment)
     ! do 2d variables
     allocate (laps_data(nx, ny, nvar2d))
     allocate (varname(nvar2d))
     allocate (varlvltype(nvar2d))
     allocate (varunits(nvar2d))
     allocate (varcomment(nvar2d))
     allocate (levels(nvar2d))

     ! initialize varcomment

     varcomment(:) = '                                                  '// &
                     '                                                  '// &
                     '                                '

     levels(:) = 0.
     startind = 1
     laps_data(:, :, startind) = u
     varname(startind) = 'usf'
     varcomment(startind) = 'sfc u-component wind            '

     startind = startind + 1
     laps_data(:, :, startind) = v
     varname(startind) = 'vsf'
     varcomment(startind) = 'sfc v-component wind            '

     startind = startind + 1
     laps_data(:, :, startind) = w
     varname(startind) = 'wsf'
     varcomment(startind) = 'sfc vertical velocity           '

     startind = startind + 1
     laps_data(:, :, startind) = t
     varname(startind) = 'tsf'
     varcomment(startind) = 'sfc temperature          '

     startind = startind + 1
     laps_data(:, :, startind) = td
     varname(startind) = 'dsf'
     varcomment(startind) = 'sfc dewpoint temperature          '

     startind = startind + 1
     laps_data(:, :, startind) = rh
     varname(startind) = 'rh '
     varcomment(startind) = 'sfc relative humidity                     '

     startind = startind + 1
     laps_data(:, :, startind) = lcb
     varname(startind) = 'lcb'
     varcomment(startind) = 'cloud base asl                          '

     startind = startind + 1
     laps_data(:, :, startind) = lct
     varname(startind) = 'lct'
     varcomment(startind) = 'cloud top asl                          '

     startind = startind + 1
     laps_data(:, :, startind) = msl
     varname(startind) = 'slp'
     varcomment(startind) = 'sea-level pressure                    '

     startind = startind + 1
     laps_data(:, :, startind) = p
     varname(startind) = 'p  '
     varcomment(startind) = 'reduced pressure               '

     startind = startind + 1
     laps_data(:, :, startind) = ps
     varname(startind) = 'psf'
     varcomment(startind) = 'surface pressure               '

     startind = startind + 1
     laps_data(:, :, startind) = lil
     varname(startind) = 'lil'
     varcomment(startind) = 'integrated liquid water                     '

     startind = startind + 1
     laps_data(:, :, startind) = tpw
     varname(startind) = 'tpw'
     varcomment(startind) = 'total precipitable water                     '

     startind = startind + 1
     laps_data(:, :, startind) = r01
     varname(startind) = 'r01'
     varcomment(startind) = 'incremental tot. liq. precip           '

     startind = startind + 1
     laps_data(:, :, startind) = rto
     varname(startind) = 'rto'
     varcomment(startind) = 'run-total liq. precip accum            '

     startind = startind + 1
     laps_data(:, :, startind) = s01
     varname(startind) = 's01'
     varcomment(startind) = 'incremental snow depth              '

     startind = startind + 1
     laps_data(:, :, startind) = sto
     varname(startind) = 'sto'
     varcomment(startind) = 'run-total snow accum                '

     startind = startind + 1
     laps_data(:, :, startind) = th
     varname(startind) = 'th '
     varcomment(startind) = 'sfc potential temperature          '

     startind = startind + 1
     laps_data(:, :, startind) = the
     varname(startind) = 'the'
     varcomment(startind) = 'sfc equiv. potential temperature              '

     startind = startind + 1
     laps_data(:, :, startind) = pbe
     varname(startind) = 'pbe'
     varcomment(startind) = 'cape                 '

     startind = startind + 1
     laps_data(:, :, startind) = nbe
     varname(startind) = 'nbe'
     varcomment(startind) = 'cin                 '

     startind = startind + 1
     laps_data(:, :, startind) = lcv
     varname(startind) = 'lcv'
     varcomment(startind) = 'cloud fraction                          '

     startind = startind + 1
     laps_data(:, :, startind) = cce
     varname(startind) = 'cce'
     varcomment(startind) = 'cloud ceiling agl              '

     startind = startind + 1
     laps_data(:, :, startind) = lmt
     varname(startind) = 'lmt'
     varcomment(startind) = 'sim. radar echo tops                         '

     startind = startind + 1
     laps_data(:, :, startind) = lmr
     varname(startind) = 'lmr'
     varcomment(startind) = 'sim. composite reflectivity            '

     startind = startind + 1
     laps_data(:, :, startind) = llr
     varname(startind) = 'llr'
     varcomment(startind) = 'sim. sfc. reflectivity                '

     startind = startind + 1
     laps_data(:, :, startind) = spt
     varname(startind) = 'spt'
     varcomment(startind) = 'sfc precip. type                       '

     startind = startind + 1
     laps_data(:, :, startind) = lhe
     varname(startind) = 'lhe'
     varcomment(startind) = 'storm relative helicity                  '

     startind = startind + 1
     laps_data(:, :, startind) = li
     varname(startind) = 'li '
     varcomment(startind) = 'lifted index              '

     startind = startind + 1
     laps_data(:, :, startind) = hi
     varname(startind) = 'hi '
     varcomment(startind) = 'heat index              '

     startind = startind + 1
     laps_data(:, :, startind) = vis
     varname(startind) = 'vis'
     varcomment(startind) = 'sfc. visibility              '

     startind = startind + 1
     laps_data(:, :, startind) = terdot
     varname(startind) = 'ter'
     varcomment(startind) = 'model terrain           '

     startind = startind + 1
     laps_data(:, :, startind) = lwout
     varname(startind) = 'lwo'
     varcomment(startind) = 'outgoing lw radiation         '

     startind = startind + 1
     laps_data(:, :, startind) = swout
     varname(startind) = 'swo'
     varcomment(startind) = 'outgoing sw radiation         '

     startind = startind + 1
     laps_data(:, :, startind) = shflux
     varname(startind) = 'shf'
     varcomment(startind) = 'sensible heat flux            '

     startind = startind + 1
     laps_data(:, :, startind) = lhflux
     varname(startind) = 'lhf'
     varcomment(startind) = 'latent heat flux          '

     startind = startind + 1
     laps_data(:, :, startind) = pblhgt
     varname(startind) = 'blh'
     varcomment(startind) = 'boundary layer depth         '

     startind = startind + 1
     laps_data(:, :, startind) = ground_t
     varname(startind) = 'tgd'
     varcomment(startind) = 'ground temperature        '

     startind = startind + 1
     laps_data(:, :, startind) = upb
     varname(startind) = 'upb'
     varcomment(startind) = 'u-component wind in pbl    '

     startind = startind + 1
     laps_data(:, :, startind) = vpb
     varname(startind) = 'vpb'
     varcomment(startind) = 'v-component wind in pbl    '

     startind = startind + 1
     laps_data(:, :, startind) = vnt
     varname(startind) = 'vnt'
     varcomment(startind) = 'ventilation index          '

     startind = startind + 1
     laps_data(:, :, startind) = ham
     varname(startind) = 'ham'
     varcomment(startind) = 'mid-level haines index     '

     startind = startind + 1
     laps_data(:, :, startind) = hah
     varname(startind) = 'hah'
     varcomment(startind) = 'high-level haines index    '

     startind = startind + 1
     laps_data(:, :, startind) = fwi
     varname(startind) = 'fwi'
     varcomment(startind) = 'fosberg fire wx index     '

     if (.not. write_to_lapsdir) then
        output_dir = trim(lfmprd_dir)//'/d'//domnum_str//'/fsf/'
     else
        output_dir = trim(laps_data_root)//'/lapsprd/fsf/'// &
                     trim(model_name)//'/'
     end if

     call cvt_fname_v3(output_dir, gtime, fcst_hhmm, 'fsf', extlen, output_file, &
                       fnlen, istatus)

     print *, 'writing 2d fields to ', trim(output_file)

     call write_laps_lfm(laps_reftime, laps_valtime, trim(output_dir), cdl_dir, &
                         'fsf', &
                         nx, ny, nvar2d, nvar2d, varname, levels, varlvltype, &
                         varunits, varcomment, nz_one, cdl_levels_one, &
                         laps_data, istatus)
     if (istatus .ne. 1) then
        print *, 'error writing laps 2d (fsf) file.'
     end if
     print *, 'done writing 2d data.'
     if ((realtime) .and. (istatus .eq. 1) .and. (make_donefile)) then
        donefile = trim(output_file)//'.done'
        call get_file_unit(doneunit)
        open (unit=doneunit, file=donefile, status='unknown')
        close (doneunit)
     end if
     deallocate (laps_data)
     deallocate (varname)
     deallocate (levels)
     deallocate (varlvltype)
     deallocate (varunits)
     deallocate (varcomment)
     deallocate (cdl_levels)
     return

  end subroutine output_laps_format

