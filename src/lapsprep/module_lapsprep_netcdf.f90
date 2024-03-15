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

module lapsprep_netcdf

! purpose
! =======
! module to contain the various output routines needed for lapsprep.
!
! subroutines contained
! =====================
! output_netcdf_format  - used to support rams 4.x initializations
!
! remarks
! =======

!
! history
! =======
! 28 nov 2000 -- original -- brent shaw
! 14 nov 2001 -- modified -- john snook

   use setup
   use laps_static
   use date_pack

   private
   public output_netcdf_format
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine output_netcdf_format(pr, ht, tp, mr, uw, vw, ww, slp, spr, lwc, ice, rai, sno, pic, rh, snocov, tskin)

      ! subroutine to output data in netcdf format.  the skin temp from
      ! laps (lsx/tgd) is used for sst.

      ! note that the p array has a bogus 2001 mb value as the last entry
      ! to designate the surface for mm5 applications.  the surface values
      ! of the state variables are contained in the last layer of their
      ! respective 3d arrays.

      ! the u/v winds are grid-relative, and are rotated to make them true.

      ! p levels must be written as integer values

      ! rh must be written as a fraction

      ! 3d variables are stored top down, and are flipped to bottom up.

      implicit none

      include 'netcdf.inc'

      real                   :: pr(z3 + 1)      !pressure levels (mb)
      real                   :: uw(x, y, z3 + 1)  !u-component of wind wrt grid (m/s)
      real                   :: vw(x, y, z3 + 1)  !v-component of wind wrt grid (m/s)
      real                   :: ww(x, y, z3 + 1)  !w-component of wind (m/s)
      real                   :: tp(x, y, z3 + 1)  !temperature (k)
      real                   :: ht(x, y, z3 + 1)  !geopotential height (m)
      real                   :: mr(x, y, z3 + 1)  !mixing ratio (kg/kg)
      real                   :: slp(x, y)      !msl pressure (pa)
      real                   :: spr(x, y)      !surface pressure (pa)
      real                   :: lwc(x, y, z3)   !cloud water mr (kg/kg)
      real                   :: ice(x, y, z3)   !cloud ice mr (kg/kg)
      real                   :: rai(x, y, z3)   !precip rain mr (kg/kg)
      real                   :: sno(x, y, z3)   !precip snow mr (kg/kg)
      real                   :: pic(x, y, z3)   !precip ice mr (kg/kg)
      real                   :: rh(x, y, z3 + 1)  !relative humidity (%)
      real                   :: snocov(x, y)   !snow cover (fract)
      real                   :: tskin(x, y)    !skin temp (k)

      ! local variables

      real, allocatable              :: pr3(:)
      real, allocatable              :: ht3(:, :, :)
      real, allocatable              :: tp3(:, :, :)
      real, allocatable              :: mr3(:, :, :)
      real, allocatable              :: rh3(:, :, :)
      real, allocatable              :: uw3(:, :, :)
      real, allocatable              :: vw3(:, :, :)
      real, allocatable              :: ww3(:, :, :)
      real, allocatable              :: sht(:, :)
      real, allocatable              :: stp(:, :)
      real, allocatable              :: srh(:, :)
      real, allocatable              :: smr(:, :)
      real, allocatable              :: suw(:, :)
      real, allocatable              :: svw(:, :)
      real, allocatable              :: lwcf(:, :, :)
      real, allocatable              :: icef(:, :, :)
      real, allocatable              :: raif(:, :, :)
      real, allocatable              :: snof(:, :, :)
      real, allocatable              :: picf(:, :, :)
      real, allocatable              :: scvf(:, :)
      real, allocatable              :: tskf(:, :)

      integer                        :: yyyyddd, valid_mm, valid_dd
      integer                        :: icode, ncid, nid, idimid(3), start, count
      integer                        :: i, j, k
      integer*2                      :: short
      real, allocatable              :: ut(:, :, :)
      real, allocatable              :: vt(:, :, :)
      real                           :: lat0
      real*8                         :: reftime
      character(len=256)            :: output_file_name
      character(len=128)            :: agridtype
      character(len=15)             :: date_string

      ! separate surface data from upper air data.

      allocate (pr3(z3))
      pr3 = pr(1:z3)

      allocate (ht3(x, y, z3))
      ht3 = ht(1:x, 1:y, 1:z3)

      allocate (tp3(x, y, z3))
      tp3 = tp(1:x, 1:y, 1:z3)

      allocate (mr3(x, y, z3))
      mr3 = mr(1:x, 1:y, 1:z3)

      allocate (rh3(x, y, z3))
      rh3 = rh(1:x, 1:y, 1:z3)

      allocate (uw3(x, y, z3))
      uw3 = uw(1:x, 1:y, 1:z3)

      allocate (vw3(x, y, z3))
      vw3 = vw(1:x, 1:y, 1:z3)

      allocate (ww3(x, y, z3))
      ww3 = ww(1:x, 1:y, 1:z3)

      allocate (sht(x, y))
      sht = ht(1:x, 1:y, z3 + 1)

      allocate (stp(x, y))
      stp = tp(1:x, 1:y, z3 + 1)

      allocate (smr(x, y))
      smr = mr(1:x, 1:y, z3 + 1)

      allocate (srh(x, y))
      srh = rh(1:x, 1:y, z3 + 1)

      allocate (suw(x, y))
      suw = uw(1:x, 1:y, z3 + 1)

      allocate (svw(x, y))
      svw = vw(1:x, 1:y, z3 + 1)

      allocate (scvf(x, y))
      scvf = snocov(1:x, 1:y)

      allocate (tskf(x, y))
      tskf = tskin(1:x, 1:y)

      ! flip 3d arrays.

      call flip_array(1, 1, z3, pr3)
      call flip_array(x, y, z3, ht3)
      call flip_array(x, y, z3, tp3)
      call flip_array(x, y, z3, mr3)
      call flip_array(x, y, z3, rh3)
      call flip_array(x, y, z3, uw3)
      call flip_array(x, y, z3, vw3)

      if (hotstart) then
         allocate (lwcf(x, y, z3))
         lwcf = lwc
         call flip_array(x, y, z3, lwcf)

         allocate (icef(x, y, z3))
         icef = ice
         call flip_array(x, y, z3, icef)

         allocate (raif(x, y, z3))
         raif = rai
         call flip_array(x, y, z3, raif)

         allocate (snof(x, y, z3))
         snof = sno
         call flip_array(x, y, z3, snof)

         allocate (picf(x, y, z3))
         picf = pic
         call flip_array(x, y, z3, picf)
         call flip_array(x, y, z3, ww3)
      end if

      ! build the output file name

      output_prefix = trim(laps_data_root)//'/lapsprd/lapsprep/cdf/laps'
      yyyyddd = valid_yyyy*1000 + valid_jjj
      call wrf_date_to_ymd(yyyyddd, valid_yyyy, valid_mm, valid_dd)
      write (date_string, '(i4.4,"-",i2.2,"-",i2.2,"-",i2.2,i2.2)') &
         valid_yyyy, valid_mm, valid_dd, valid_hh, valid_min
      output_file_name = trim(output_prefix)//':'//date_string

      !  create laps model netcdf file.

      icode = nf_create(trim(output_file_name), nf_clobber, ncid)
      if (icode .ne. 0) then
         print *, 'could not open output file: ', trim(output_file_name)
         stop
      end if

      !  define netcdf grid dimensions.

      icode = nf_def_dim(ncid, 'nx', x, 1)
      icode = nf_def_dim(ncid, 'ny', y, 2)
      icode = nf_def_dim(ncid, 'nz', z3, 3)
      icode = nf_def_dim(ncid, 'len', 128, 4)
      icode = nf_def_dim(ncid, 'pt', 1, 5)

      !  define static and grid variables.

      icode = nf_def_var(ncid, 'grid_type', nf_char, 1, 4, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 15, 'grid projection')
      icode = nf_def_var(ncid, 'nx', nf_short, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 18, 'number of x points')
      icode = nf_def_var(ncid, 'ny', nf_short, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 18, 'number of y points')
      icode = nf_def_var(ncid, 'np', nf_short, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 25, 'number of pressure levels')
      icode = nf_def_var(ncid, 'dx', nf_real, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 16, 'x grid increment')
      icode = nf_put_att_text(ncid, nid, 'units', 10, 'kilometers')
      icode = nf_def_var(ncid, 'dy', nf_real, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 16, 'y grid increment')
      icode = nf_put_att_text(ncid, nid, 'units', 10, 'kilometers')
      icode = nf_def_var(ncid, 'lat0', nf_real, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 13, 'pole latitude')
      icode = nf_put_att_text(ncid, nid, 'units', 13, 'degrees north')
      icode = nf_def_var(ncid, 'lat1', nf_real, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 14, 'first latitude')
      icode = nf_put_att_text(ncid, nid, 'units', 13, 'degrees north')
      icode = nf_def_var(ncid, 'lat2', nf_real, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 15, 'second latitude')
      icode = nf_put_att_text(ncid, nid, 'units', 13, 'degrees north')
      icode = nf_def_var(ncid, 'lon0', nf_real, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 14, 'pole longitude')
      icode = nf_put_att_text(ncid, nid, 'units', 12, 'degrees east')
      icode = nf_def_var(ncid, 'swlat', nf_real, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 18, 'southwest latitude')
      icode = nf_put_att_text(ncid, nid, 'units', 13, 'degrees north')
      icode = nf_def_var(ncid, 'swlon', nf_real, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 19, 'southwest longitude')
      icode = nf_put_att_text(ncid, nid, 'units', 12, 'degrees east')
      icode = nf_def_var(ncid, 'nelat', nf_real, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 18, 'northeast latitude')
      icode = nf_put_att_text(ncid, nid, 'units', 13, 'degrees north')
      icode = nf_def_var(ncid, 'nelon', nf_real, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 19, 'northeast longitude')
      icode = nf_put_att_text(ncid, nid, 'units', 12, 'degrees east')
      icode = nf_def_var(ncid, 'reftime', nf_double, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 14, 'reference time')
      icode = nf_put_att_text(ncid, nid, 'units', 35, 'seconds since (1970-1-1 00:00:00.0)')
      icode = nf_def_var(ncid, 'valtime', nf_double, 1, 5, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 10, 'valid time')
      icode = nf_put_att_text(ncid, nid, 'units', 35, 'seconds since (1970-1-1 00:00:00.0)')
      icode = nf_def_var(ncid, 'level', nf_real, 1, 3, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 13, 'level of data')
      icode = nf_put_att_text(ncid, nid, 'units', 2, 'mb')
      icode = nf_def_var(ncid, 'pr', nf_real, 1, 3, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 8, 'pressure')
      icode = nf_put_att_text(ncid, nid, 'units', 7, 'pascals')

      !  define surface variables.

      idimid(1) = 1
      idimid(2) = 2
      icode = nf_def_var(ncid, 'sht', nf_real, 2, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 14, 'surface height')
      icode = nf_put_att_text(ncid, nid, 'units', 6, 'meters')
      icode = nf_def_var(ncid, 'spr', nf_real, 2, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 16, 'surface pressure')
      icode = nf_put_att_text(ncid, nid, 'units', 7, 'pascals')
      icode = nf_def_var(ncid, 'slp', nf_real, 2, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 23, 'mean sea-level pressure')
      icode = nf_put_att_text(ncid, nid, 'units', 7, 'pascals')
      icode = nf_def_var(ncid, 'stp', nf_real, 2, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 19, 'surface temperature')
      icode = nf_put_att_text(ncid, nid, 'units', 6, 'kelvin')
      icode = nf_def_var(ncid, 'smr', nf_real, 2, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 20, 'surface mixing ratio')
      icode = nf_put_att_text(ncid, nid, 'units', 5, 'kg/kg')
      icode = nf_def_var(ncid, 'srh', nf_real, 2, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 25, 'surface relative humidity')
      icode = nf_put_att_text(ncid, nid, 'units', 7, 'percent')
      icode = nf_def_var(ncid, 'suw', nf_real, 2, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 14, 'surface u-wind')
      icode = nf_put_att_text(ncid, nid, 'units', 14, 'meters/seconds')
      icode = nf_def_var(ncid, 'svw', nf_real, 2, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 14, 'surface v-wind')
      icode = nf_put_att_text(ncid, nid, 'units', 14, 'meters/seconds')
      icode = nf_def_var(ncid, 'scv', nf_real, 2, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 10, 'snow cover')
      icode = nf_put_att_text(ncid, nid, 'units', 8, 'fraction')
      icode = nf_def_var(ncid, 'tsk', nf_real, 2, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 16, 'skin temperature')
      icode = nf_put_att_text(ncid, nid, 'units', 6, 'kelvin')

      !  define upper-air variables.

      idimid(3) = 3
      icode = nf_def_var(ncid, 'ht', nf_real, 3, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 6, 'height')
      icode = nf_put_att_text(ncid, nid, 'units', 6, 'meters')
      icode = nf_def_var(ncid, 'tp', nf_real, 3, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 11, 'temperature')
      icode = nf_put_att_text(ncid, nid, 'units', 6, 'kelvin')
      icode = nf_def_var(ncid, 'mr', nf_real, 3, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 12, 'mixing ratio')
      icode = nf_put_att_text(ncid, nid, 'units', 5, 'kg/kg')
      icode = nf_def_var(ncid, 'rh', nf_real, 3, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 17, 'relative humidity')
      icode = nf_put_att_text(ncid, nid, 'units', 7, 'percent')
      icode = nf_def_var(ncid, 'uw', nf_real, 3, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 19, 'u-component of wind')
      icode = nf_put_att_text(ncid, nid, 'units', 13, 'meters/second')
      icode = nf_def_var(ncid, 'vw', nf_real, 3, idimid, nid)
      icode = nf_put_att_text(ncid, nid, 'long_name', 19, 'v-component of wind')
      icode = nf_put_att_text(ncid, nid, 'units', 13, 'meters/second')

      if (hotstart) then
         icode = nf_def_var(ncid, 'lwc', nf_real, 3, idimid, nid)
         icode = nf_put_att_text(ncid, nid, 'long_name', 31 &
                                 , 'cloud liquid water mixing ratio')
         icode = nf_put_att_text(ncid, nid, 'units', 5, 'kg/kg')
         icode = nf_def_var(ncid, 'ice', nf_real, 3, idimid, nid)
         icode = nf_put_att_text(ncid, nid, 'long_name', 22 &
                                 , 'cloud ice mixing ratio')
         icode = nf_put_att_text(ncid, nid, 'units', 5, 'kg/kg')
         icode = nf_def_var(ncid, 'sno', nf_real, 3, idimid, nid)
         icode = nf_put_att_text(ncid, nid, 'long_name', 31 &
                                 , 'precipitating snow mixing ratio')
         icode = nf_put_att_text(ncid, nid, 'units', 5, 'kg/kg')
         icode = nf_def_var(ncid, 'rai', nf_real, 3, idimid, nid)
         icode = nf_put_att_text(ncid, nid, 'long_name', 31 &
                                 , 'precipitating rain mixing ratio')
         icode = nf_put_att_text(ncid, nid, 'units', 5, 'kg/kg')
         icode = nf_def_var(ncid, 'pic', nf_real, 3, idimid, nid)
         icode = nf_put_att_text(ncid, nid, 'long_name', 30 &
                                 , 'precipitating ice mixing ratio')
         icode = nf_put_att_text(ncid, nid, 'units', 5, 'kg/kg')
         icode = nf_def_var(ncid, 'ww', nf_real, 3, idimid, nid)
         icode = nf_put_att_text(ncid, nid, 'long_name', 19 &
                                 , 'w-component of wind')
         icode = nf_put_att_text(ncid, nid, 'units', 13, 'meters/second')
      end if

      icode = nf_enddef(ncid)

      !  create reftime.

      call adate_to_i4time(laps_file_time, reftime)

      !  fill static and grid variables.

      if (grid_type(1:8) .eq. 'latlon') then
         agridtype = 'lat-lon'
      else if ((grid_type(1:24) .eq. 'secant lambert conformal') .or. &
               (grid_type(1:28) .eq. 'tangential lambert conformal')) then
         agridtype = 'lambert-conformal'
      else if (grid_type(1:19) .eq. 'polar stereographic') then
         agridtype = 'polar-stereographic'
      else
         print '(a,a,a)', 'netcdf unsupported map projection: ', &
            trim(grid_type), '.  i quit.'
         stop 'unsupported projection'
      end if

      start = 1
      count = len_trim(agridtype)
      icode = nf_inq_varid(ncid, 'grid_type', nid)
      icode = nf_put_vara_text(ncid, nid, start, count, agridtype)
      short = x
      icode = nf_inq_varid(ncid, 'nx', nid)
      icode = nf_put_var_int2(ncid, nid, short)
      short = y
      icode = nf_inq_varid(ncid, 'ny', nid)
      icode = nf_put_var_int2(ncid, nid, short)
      short = z3
      icode = nf_inq_varid(ncid, 'np', nid)
      icode = nf_put_var_int2(ncid, nid, short)
      icode = nf_inq_varid(ncid, 'dx', nid)
      icode = nf_put_var_real(ncid, nid, dx*1000.)
      icode = nf_inq_varid(ncid, 'dy', nid)
      icode = nf_put_var_real(ncid, nid, dy*1000.)
      lat0 = 90.
      icode = nf_inq_varid(ncid, 'lat0', nid)
      icode = nf_put_var_real(ncid, nid, lat0)
      icode = nf_inq_varid(ncid, 'lat1', nid)
      icode = nf_put_var_real(ncid, nid, latin1)
      icode = nf_inq_varid(ncid, 'lat2', nid)
      icode = nf_put_var_real(ncid, nid, latin2)
      icode = nf_inq_varid(ncid, 'lon0', nid)
      icode = nf_put_var_real(ncid, nid, lov)
      icode = nf_inq_varid(ncid, 'swlat', nid)
      icode = nf_put_var_real(ncid, nid, la1)
      icode = nf_inq_varid(ncid, 'swlon', nid)
      icode = nf_put_var_real(ncid, nid, lo1)
      icode = nf_inq_varid(ncid, 'nelat', nid)
      icode = nf_put_var_real(ncid, nid, la2)
      icode = nf_inq_varid(ncid, 'nelon', nid)
      icode = nf_put_var_real(ncid, nid, lo2)
      icode = nf_inq_varid(ncid, 'reftime', nid)
      icode = nf_put_var_double(ncid, nid, reftime)
      icode = nf_inq_varid(ncid, 'valtime', nid)
      icode = nf_put_var_double(ncid, nid, reftime)
      icode = nf_inq_varid(ncid, 'level', nid)
      icode = nf_put_var_real(ncid, nid, pr3)
      pr3 = pr3*100.
      icode = nf_inq_varid(ncid, 'pr', nid)
      icode = nf_put_var_real(ncid, nid, pr3)

      ! rotate the surface u and v winds to true if polar grid projection.

      allocate (ut(x, y, 1)) ! array for true u-winds
      allocate (vt(x, y, 1)) ! array for true v-winds

      if ((agridtype(1:3) .eq. 'pol') .or. (agridtype(1:3) .eq. 'lam')) then
         do j = 1, y
         do i = 1, x
            call uvgrid_to_uvtrue(suw(i, j), svw(i, j), &
                                  ut(i, j, 1), vt(i, j, 1), &
                                  lons(i, j))
         end do
         end do
      else
         ut(:, :, 1) = suw
         vt(:, :, 1) = svw
      end if

      !  fill surface variables.

      icode = nf_inq_varid(ncid, 'sht', nid)
      icode = nf_put_var_real(ncid, nid, sht)
      icode = nf_inq_varid(ncid, 'spr', nid)
      icode = nf_put_var_real(ncid, nid, spr)
      icode = nf_inq_varid(ncid, 'slp', nid)
      icode = nf_put_var_real(ncid, nid, slp)
      icode = nf_inq_varid(ncid, 'stp', nid)
      icode = nf_put_var_real(ncid, nid, stp)
      icode = nf_inq_varid(ncid, 'smr', nid)
      icode = nf_put_var_real(ncid, nid, smr)
      icode = nf_inq_varid(ncid, 'srh', nid)
      icode = nf_put_var_real(ncid, nid, srh)
      icode = nf_inq_varid(ncid, 'suw', nid)
      icode = nf_put_var_real(ncid, nid, ut)
      icode = nf_inq_varid(ncid, 'svw', nid)
      icode = nf_put_var_real(ncid, nid, vt)
      icode = nf_inq_varid(ncid, 'scv', nid)
      icode = nf_put_var_real(ncid, nid, scvf)
      icode = nf_inq_varid(ncid, 'tsk', nid)
      icode = nf_put_var_real(ncid, nid, tskf)

      deallocate (ut)
      deallocate (vt)

      ! rotate the u and v winds to true if polar grid projection.

      allocate (ut(x, y, z3)) ! array for true u-winds
      allocate (vt(x, y, z3)) ! array for true v-winds

      if ((agridtype(1:3) .eq. 'pol') .or. (agridtype(1:3) .eq. 'lam')) then
         level_loop: do k = 1, z3
         rotate_winds_j: do j = 1, y
         rotate_winds_i: do i = 1, x
            call uvgrid_to_uvtrue(uw3(i, j, k), vw3(i, j, k), &
                                  ut(i, j, k), vt(i, j, k), &
                                  lons(i, j))
         end do rotate_winds_i
         end do rotate_winds_j
         end do level_loop
      else
         ut = uw3
         vt = vw3
      end if

      !  fill upper air variables.

      icode = nf_inq_varid(ncid, 'ht', nid)
      icode = nf_put_var_real(ncid, nid, ht3)
      icode = nf_inq_varid(ncid, 'tp', nid)
      icode = nf_put_var_real(ncid, nid, tp3)
      icode = nf_inq_varid(ncid, 'mr', nid)
      icode = nf_put_var_real(ncid, nid, mr3)
      icode = nf_inq_varid(ncid, 'rh', nid)
      icode = nf_put_var_real(ncid, nid, rh3)
      icode = nf_inq_varid(ncid, 'uw', nid)
      icode = nf_put_var_real(ncid, nid, ut)
      icode = nf_inq_varid(ncid, 'vw', nid)
      icode = nf_put_var_real(ncid, nid, vt)

      if (hotstart) then
         icode = nf_inq_varid(ncid, 'lwc', nid)
         icode = nf_put_var_real(ncid, nid, lwcf)
         icode = nf_inq_varid(ncid, 'ice', nid)
         icode = nf_put_var_real(ncid, nid, icef)
         icode = nf_inq_varid(ncid, 'rai', nid)
         icode = nf_put_var_real(ncid, nid, raif)
         icode = nf_inq_varid(ncid, 'sno', nid)
         icode = nf_put_var_real(ncid, nid, snof)
         icode = nf_inq_varid(ncid, 'pic', nid)
         icode = nf_put_var_real(ncid, nid, picf)
         icode = nf_inq_varid(ncid, 'ww', nid)
         icode = nf_put_var_real(ncid, nid, ww3)
         deallocate (lwcf)
         deallocate (icef)
         deallocate (raif)
         deallocate (snof)
         deallocate (picf)
      end if

      deallocate (ut)
      deallocate (vt)
      deallocate (pr3)
      deallocate (ht3)
      deallocate (tp3)
      deallocate (mr3)
      deallocate (rh3)
      deallocate (uw3)
      deallocate (vw3)
      deallocate (ww3)
      deallocate (sht)
      deallocate (stp)
      deallocate (smr)
      deallocate (srh)
      deallocate (suw)
      deallocate (svw)
      deallocate (scvf)
      deallocate (tskf)

      ! close the netcdf file.

      icode = nf_close(ncid)

      return
   end subroutine output_netcdf_format

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine flip_array(nx, ny, nz, a)

      implicit none

      integer nx, ny, nz, i, j, k, kk

      real a(nx, ny, nz), tmp(nx, ny, nz/2)

      kk = nz + 1
      do k = 1, nz/2
         kk = kk - 1
         do j = 1, ny
         do i = 1, nx
            tmp(i, j, k) = a(i, j, k)
            a(i, j, k) = a(i, j, kk)
         end do
         end do
      end do

      kk = nz/2 + 1
      do k = (nz + 1)/2 + 1, nz
         kk = kk - 1
         do j = 1, ny
         do i = 1, nx
            a(i, j, k) = tmp(i, j, kk)
         end do
         end do
      end do

      return
   end subroutine flip_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine adate_to_i4time(adate, i4time)

      implicit none

      real*8 i4time

      integer*4 iyear, iday, ihour, imin, lp

      character*9 adate

      read (adate(1:2), '(i2)') iyear
      read (adate(3:5), '(i3)') iday
      read (adate(6:7), '(i2)') ihour
      read (adate(8:9), '(i2)') imin

      !  valid for years 1960-2060.

      if (iyear .lt. 60) iyear = iyear + 100

      lp = (iyear + 3 - 60)/4

      i4time = (iyear - 60)*31536000 &
               + (iday - 1 + lp)*86400 &
               + ihour*3600 &
               + imin*60

      return
   end subroutine adate_to_i4time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module lapsprep_netcdf
