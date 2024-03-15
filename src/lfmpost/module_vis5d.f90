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

module vis5d

! module to contain fortran routines needed by the model post processor
! to create vis5d output.  this module also depends on the v5d.c
! vis5d api.

   use setup
   use time_utils
   implicit none
   private
   include 'v5df.h'
   integer                               :: v5dtimes(maxtimes)
   integer                               :: v5ddates(maxtimes)
   integer                               :: num2dvar
   integer                               :: num3dvar
   integer                               :: totalvars
   integer                               :: ret
   character(len=10)                     :: varname(maxvars)
   character(len=255)                    :: v5dfile
   public v5dinit
   public v5dout
   public v5dend

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine v5dinit(compression)

      ! initializes the vis5d file.  only argument provided is the compression
      ! factor.  the rest of the items are input from the setup module to
      ! call various other routines.

      implicit none

      include 'v5df.h'

      ! declare arguments
      integer, intent(in)                   :: compression

      ! local variables
      integer                               :: yyyy, yy, mm, dd
      integer                               :: hh, min, doy
      integer                               :: i, j, k, t
      integer, parameter                    :: vertical_code = 3
      real                                  :: levels(maxlevels)
      integer                               :: nlevels(maxvars)
      integer                               :: varindex
      integer                               :: v5d_proj
      real                                  :: h, msf, pole_lat
      real, parameter                       :: rad_per_deg = .0174533
      real                                  :: v5d_proj_args(100)
      character(len=9)                      :: file_date_string
      integer                               :: v5d_ntimes
      character(len=3)                      :: domnum_str
    !!!!! begin code !!!!!

      ! check compression argument
      if ((compression .ne. 1) .and. &
          (compression .ne. 2) .and. &
          (compression .ne. 4)) then
         print '(a,i2)', 'v5dinit: invalid compression: ', compression
         print '(a)', '  must be 1, 2, or 4 bytes per point.'
         stop
      end if

      ! build the date/time values for vis5d
      v5d_ntimes = 0
      do t = 1, num_times_to_proc, time_index_inc
         v5d_ntimes = v5d_ntimes + 1
         read (times_to_proc(t), fmt='(i4.4)') yyyy
         read (times_to_proc(t), fmt='(2x,i2.2)') yy
         read (times_to_proc(t), fmt='(5x,i2.2)') mm
         read (times_to_proc(t), fmt='(8x,i2.2)') dd
         read (times_to_proc(t), fmt='(11x,i2.2)') hh
         read (times_to_proc(t), fmt='(14x,i2.2)') min
         doy = compute_day_of_year(yyyy, mm, dd)
         v5ddates(v5d_ntimes) = yy*1000 + doy
         v5dtimes(v5d_ntimes) = hh*10000 + min*100
         print '(a,i5.5,1x,i6.6)', 'vis5d date/time: ', v5ddates(v5d_ntimes), &
            v5dtimes(v5d_ntimes)
      end do

      ! build the vis5d output file name

      write (file_date_string, fmt='(i5.5,i4.4)') v5ddates(1), v5dtimes(1)/100
      write (domnum_str, '("d",i2.2)') domain_num
      v5dfile = trim(lfmprd_dir)//'/'//domnum_str// &
                '/v5d/'//file_date_string//'.v5d'

      ! build the levels array (descending pressures in mb)
      levels(1:kprs) = prslvl/100.

      ! build the projection arguments for vis5d

      if (trim(projection) .eq. 'lambert conformal') then
         v5d_proj = 2
         v5d_proj_args(1) = proj%truelat2
         v5d_proj_args(2) = proj%truelat1
         v5d_proj_args(3) = float(ny) - proj%polej
         v5d_proj_args(4) = proj%polei
         v5d_proj_args(5) = -proj%stdlon
         v5d_proj_args(6) = proj%dx*0.001
      else if (trim(projection) .eq. 'polar stereographic') then
         v5d_proj = 3
         if (proj%truelat1 .ge. 0.) then
            pole_lat = 90.
            h = 1.
         else
            pole_lat = -90.
            h = -1.
         end if
         v5d_proj_args(1) = pole_lat
         v5d_proj_args(2) = -proj%stdlon
         v5d_proj_args(3) = float(ny) - proj%polej
         v5d_proj_args(4) = proj%polei

         ! compute the grid spacing at the pole using
         ! map scale factor
         msf = (1.+h*sin(proj%truelat1*rad_per_deg))/ &
               (1.+h*sin(pole_lat*rad_per_deg))
         v5d_proj_args(5) = proj%dx*0.001/msf
         print *, 'dx at pole = ', v5d_proj_args(5)
      else if (trim(projection) .eq. 'mercator') then
         v5d_proj = 1
         v5d_proj_args(1) = corner_lats(2)  ! north boundary
         v5d_proj_args(2) = -corner_lons(2)  ! west boundary
         ! compute increment between rows in degrees
         v5d_proj_args(3) = (corner_lats(2) - corner_lats(1))/(float(ny) - 1.)
         ! compute increment between columns in degrees
         if (corner_lons(3) .gt. corner_lons(2)) then
            v5d_proj_args(4) = (corner_lons(3) - corner_lons(2))/(float(nx) - 1.)
         else
            v5d_proj_args(4) = (corner_lons(3) - (corner_lons(2) - 360.))/(float(nx) - 1.)
         end if
      else
         print '(2a)', 'v5dinit: unsupported map projection:', projection
         stop
      end if

!   call getproj(proj_tmp,proj_cent_lon, proj_cent_lat, grid_spacing*0.001, &
!           truelat1,truelat2,nx,ny,v5d_proj,v5d_proj_args, corner_lats(2), &
!             corner_lats(1), corner_lons(2), corner_lons(3))

      ! set up the vertical levels
      levels(1:kprs) = prslvl(:)*0.01

      ! set up the varnames and number of levels.  do the 3d variables first.
      ! these consist of all of the fua variables plus vorticity and 3d
      ! dewpoint, which will be derived by the calling routine.

      varindex = 1
      ! height (fua/ht)
      varname(varindex) = 'hgt       '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! temperature (fua/t3)
      varname(varindex) = 't         '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! dewpoint (not in fua)
      varname(varindex) = 'td        '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! rh (fua/rh3)
      varname(varindex) = 'rh        '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! u (fua/u3)
      varname(varindex) = 'u         '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! v (fua/v3)
      varname(varindex) = 'v         '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! w (fua/w3)
      varname(varindex) = 'w         '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! omega  (fua/om)
      varname(varindex) = 'om        '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! absolute vorticity (not in laps files)
      varname(varindex) = 'avort     '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! specific humidity (fua/sh)
      varname(varindex) = 'sh        '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! cloud liquid water (fua/lwc)
      varname(varindex) = 'lwc       '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! cloud ice (fua/ice)
      varname(varindex) = 'ice       '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! rain (fua/rai)
      varname(varindex) = 'rai       '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! snow (fua/rai)
      varname(varindex) = 'sno       '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! precip ice/graupel (fua/pic)
      varname(varindex) = 'pic       '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! total condensate
      varname(varindex) = 'cond      '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! radar reflectivity (fua/ref)
      varname(varindex) = 'ref       '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! precip type code (fua/pty)
      varname(varindex) = 'pty       '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! turbulent kinetic energy (fua/tke)
      varname(varindex) = 'tke       '
      nlevels(varindex) = kprs
      varindex = varindex + 1

      ! print out number of 3d variables as a sanity check
      num3dvar = varindex - 1
      print '(a,i3)', 'v5dinit: number of 3d variables to output: ', num3dvar

      ! set up 2d variables for output

      ! 1000-500 mb thickness (not in fua/fsf)
      varname(varindex) = 'thick     '
      nlevels(varindex) = 1
      varindex = varindex + 1
      ! surface temperature (fsf/t)
      varname(varindex) = 't_sfc     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! surface dewpoint (fsf/td)
      varname(varindex) = 'td_sfc    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! surface rh (fsf/rh)
      varname(varindex) = 'rh_sfc    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! surface
      ! surface u (fsf/u)
      varname(varindex) = 'u_sfc     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! surface v (fsf/v)
      varname(varindex) = 'v_sfc     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! surface w (fsf/w)
      varname(varindex) = 'w_sfc     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! mean slp (fsf/msl)
      varname(varindex) = 'mslp      '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! surface pressure (fsf/ps)
      varname(varindex) = 'p_sfc     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! cloud bases (fsf/lcb)
      varname(varindex) = 'cldbas    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! cloud tops (fsf/lct)
      varname(varindex) = 'cldtop    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! cloud cover (fsf/lcv)
      varname(varindex) = 'cldcov    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! cloud ceiling (fsf/cce)
      varname(varindex) = 'cldcei    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! integrated liquid water depth (fsf/lil)
      varname(varindex) = 'intliq    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! total precip. water (fsf/tpw)
      varname(varindex) = 'tpw       '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! incremental precip accum (fsf/r01)
      varname(varindex) = 'pcpinc    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sim. total precip accum (fsf/rto)
      varname(varindex) = 'pcptot    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! incremental convective precip accum (not in laps output)
      varname(varindex) = 'coninc    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sim. total convective precip accum (not in laps output)
      varname(varindex) = 'contot    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! incremental snow accum (fsf/s01)
      varname(varindex) = 'snoinc    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sim. total precip accum (fsf/sto)
      varname(varindex) = 'snotot    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sfc precip type code (fsf/spt)
      varname(varindex) = 'sfcpcp    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sfc potential temp (fsf_th)
      varname(varindex) = 'theta     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sfc equive pot temp (fsf/the)
      varname(varindex) = 'thetae    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sfc cape (fsf/pbe)
      varname(varindex) = 'cape      '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sfc cin (fsf/nbe)
      varname(varindex) = 'cin       '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sfc lifted index (fsf/li)
      varname(varindex) = 'li        '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! storm rel. helicity (fsf/lhe)
      varname(varindex) = 'srhel     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! low-level reflectivity (fsf/llr)
      varname(varindex) = 'llref'
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! max reflectivity (fsf/lmr)
      varname(varindex) = 'maxref    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! echo tops (fsf/lmt)
      varname(varindex) = 'echtop    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! heat index (fsf/hi)
      varname(varindex) = 'htidx     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! visibility (fsf/vis)
      varname(varindex) = 'vis       '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! model terrain (setup)
      varname(varindex) = 'topo      '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! snow cover
      varname(varindex) = 'snowcov   '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! outgoing lw radiation
      varname(varindex) = 'lwout     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! outgoing sw radiaion
      varname(varindex) = 'swout     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! incoming lw radiation
      varname(varindex) = 'lwdown    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! incoming sw radiation
      varname(varindex) = 'swdown    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sensible heat flux
      varname(varindex) = 'shflux    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! latent heat flux
      varname(varindex) = 'lhflux    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! pbl height
      varname(varindex) = 'pblhgt    '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! mean u wind in pbl
      varname(varindex) = 'u_pbl     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! mean v wind in pbl
      varname(varindex) = 'v_pbl     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! ground temp
      varname(varindex) = 't_gnd     '
      nlevels(varindex) = 1
      varindex = varindex + 1

      ! sanity check print of two-d variables
      num2dvar = varindex - 1 - num3dvar
      print '(a,i3)', 'v5dinit: number of 2d variables to output: ', num2dvar
      totalvars = num2dvar + num3dvar

      ! initialize the file with v5dcreate
      ret = v5dcreate(v5dfile, v5d_ntimes, &
                      totalvars, ny, nx, nlevels, &
                      varname, v5dtimes, v5ddates, compression, v5d_proj, &
                      v5d_proj_args, vertical_code, levels)
      if (ret .eq. 0) then
         print '(a)', '---error opening vis5d output file: v5dfile'
         call abort
      end if
      print '(a)', 'vis5d header written.'
      return
   end subroutine v5dinit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine v5dout(timestep, zprs, tprs, tdprs, rhprs, uprs, vprs, wprs, omprs, &
                     abs_vort, shprs, cldliqmr_prs, cldicemr_prs, rainmr_prs, &
                     snowmr_prs, graupelmr_prs, refl_prs, pcptype_prs, tkeprs, &
                     thick_10_5, tsfc, tdsfc, rhsfc, usfc, vsfc, wsfc, &
                     pmsl, psfc, cldbase, cldtop, cldamt, ceiling, &
                     intliqwater, totpcpwater, pcp_inc, pcp_tot, &
                     con_pcp_inc, con_pcp_tot, &
                     snow_inc, snow_tot, pcptype_sfc, thetasfc, thetaesfc, &
                     cape, cin, liftedind, srhel, refl_sfc, max_refl, echo_tops, &
                     heatind, visibility, snowcover, lwout, swout, lwdown, swdown, &
                     shflux, lhflux, &
                     pblhgt, upbl, vpbl, ground_t)

      ! subroutine to output the above variables into the already created
      ! vis5d file.  note that the variables must be output in the same order
      ! (i.e., the same index) as the names were defined in v5dinit.

      implicit none
      include 'v5df.h'

      ! declare the arguments
      integer                 :: timestep
      real                    :: zprs(nx, ny, kprs)
      real                    :: tprs(nx, ny, kprs)
      real                    :: tdprs(nx, ny, kprs)
      real                    :: rhprs(nx, ny, kprs)
      real                    :: uprs(nx, ny, kprs)
      real                    :: vprs(nx, ny, kprs)
      real                    :: wprs(nx, ny, kprs)
      real                    :: omprs(nx, ny, kprs)
      real                    :: abs_vort(nx, ny, kprs)
      real                    :: shprs(nx, ny, kprs)
      real                    :: cldliqmr_prs(nx, ny, kprs)
      real                    :: cldicemr_prs(nx, ny, kprs)
      real                    :: rainmr_prs(nx, ny, kprs)
      real                    :: snowmr_prs(nx, ny, kprs)
      real                    :: graupelmr_prs(nx, ny, kprs)
      real                    :: refl_prs(nx, ny, kprs)
      real                    :: pcptype_prs(nx, ny, kprs)
      real                    :: tkeprs(nx, ny, kprs)
      real                    :: thick_10_5(nx, ny)
      real                    :: tsfc(nx, ny)
      real                    :: tdsfc(nx, ny)
      real                    :: rhsfc(nx, ny)
      real                    :: usfc(nx, ny)
      real                    :: vsfc(nx, ny)
      real                    :: wsfc(nx, ny)
      real                    :: pmsl(nx, ny)
      real                    :: psfc(nx, ny)
      real                    :: cldbase(nx, ny)
      real                    :: cldtop(nx, ny)
      real                    :: cldamt(nx, ny)
      real                    :: ceiling(nx, ny)
      real                    :: intliqwater(nx, ny)
      real                    :: totpcpwater(nx, ny)
      real                    :: pcp_inc(nx, ny)
      real                    :: pcp_tot(nx, ny)
      real                    :: con_pcp_inc(nx, ny)
      real                    :: con_pcp_tot(nx, ny)

      real                    :: snow_inc(nx, ny)
      real                    :: snow_tot(nx, ny)
      real                    :: pcptype_sfc(nx, ny)
      real                    :: thetasfc(nx, ny)
      real                    :: thetaesfc(nx, ny)
      real                    :: cape(nx, ny)
      real                    :: cin(nx, ny)
      real                    :: liftedind(nx, ny)
      real                    :: srhel(nx, ny)
      real                    :: refl_sfc(nx, ny)
      real                    :: max_refl(nx, ny)
      real                    :: echo_tops(nx, ny)
      real                    :: heatind(nx, ny)
      real                    :: visibility(nx, ny)
      real                    :: snowcover(nx, ny)
      real                    :: lwout(nx, ny)
      real                    :: swout(nx, ny)
      real                    :: lwdown(nx, ny)
      real                    :: swdown(nx, ny)
      real                    :: shflux(nx, ny)
      real                    :: lhflux(nx, ny)
      real                    :: pblhgt(nx, ny)
      real                    :: upbl(nx, ny)
      real                    :: vpbl(nx, ny)
      real                    :: ground_t(nx, ny)

      ! local variables
      real, allocatable       :: data2d(:, :)
      real, allocatable       :: data3d(:, :, :)
      real, allocatable       :: condmr(:, :, :)
      integer                 :: varindex
      integer                 :: ret

      varindex = 1

      ! do the 3d variables first
      allocate (data3d(nx, ny, kprs))

      call post2v5d(zprs, nx, ny, kprs, data3d)
      data3d = data3d*0.1
      print '(2a)', 'v5dout: writing zprs*0.1 as varname: ', varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(tprs, nx, ny, kprs, data3d)
      data3d = data3d - 273.15
      print '(2a)', 'v5dout: writing tprs-273.15 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(tdprs, nx, ny, kprs, data3d)
      data3d = data3d - 273.15
      print '(2a)', 'v5dout: writing tdprs-273.15 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(rhprs, nx, ny, kprs, data3d)
      print '(2a)', 'v5dout: writing rhprs as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(uprs, nx, ny, kprs, data3d)
      print '(2a)', 'v5dout: writing uprs as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(vprs, nx, ny, kprs, data3d)
      print '(2a)', 'v5dout: writing vprs as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(wprs, nx, ny, kprs, data3d)
      print '(2a)', 'v5dout: writing wprs as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(omprs, nx, ny, kprs, data3d)
      print '(2a)', 'v5dout: writing omprs as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(abs_vort, nx, ny, kprs, data3d)
      data3d = data3d*100000
      print '(2a)', 'v5dout: writing abs_vort*10000 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(shprs, nx, ny, kprs, data3d)
      data3d = data3d*1000.
      print '(2a)', 'v5dout: writing shprs*1000 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(cldliqmr_prs, nx, ny, kprs, data3d)
      data3d = data3d*1000.
      print '(2a)', 'v5dout: writing cldliqmr_prs*1000 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(cldicemr_prs, nx, ny, kprs, data3d)
      data3d = data3d*1000.
      print '(2a)', 'v5dout: writing cldicemr_prs*1000 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(rainmr_prs, nx, ny, kprs, data3d)
      data3d = data3d*1000.
      print '(2a)', 'v5dout: writing rainmr_prs*1000 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(snowmr_prs, nx, ny, kprs, data3d)
      data3d = data3d*1000.
      print '(2a)', 'v5dout: writing snowmr_prs*1000 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(graupelmr_prs, nx, ny, kprs, data3d)
      data3d = data3d*1000.
      print '(2a)', 'v5dout: writing graupelmr_prs*1000 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      ! generate total condensate field
      allocate (condmr(nx, ny, kprs))
      condmr = cldliqmr_prs + cldicemr_prs + rainmr_prs + snowmr_prs &
               + graupelmr_prs
      call post2v5d(condmr, nx, ny, kprs, data3d)
      data3d = data3d*1000.
      print '(2a)', 'v5dout: writing condmr*1000 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      deallocate (condmr)
      varindex = varindex + 1

      call post2v5d(refl_prs, nx, ny, kprs, data3d)
      where (data3d .lt. 0.) data3d = 0.
      print '(2a)', 'v5dout: writing refl_prs as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(pcptype_prs, nx, ny, kprs, data3d)
      print '(2a)', 'v5dout: writing pcptype_prs as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      call post2v5d(tkeprs, nx, ny, kprs, data3d)
      print '(2a)', 'v5dout: writing tkeprs as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data3d)
      varindex = varindex + 1

      deallocate (data3d)
      ! now do the 2-d variables
      allocate (data2d(nx, ny))

      call post2v5d(thick_10_5, nx, ny, 1, data2d)
      data2d = data2d*0.1
      print '(2a)', 'v5dout: writing thick_10_5*0.1 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(tsfc, nx, ny, 1, data2d)
      data2d = data2d - 273.15
      print '(2a)', 'v5dout: writing tsfc-273.15 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(tdsfc, nx, ny, 1, data2d)
      data2d = data2d - 273.15
      print '(2a)', 'v5dout: writing tdsfc-273.15 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(rhsfc, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing rhsfc as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(usfc, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing usfc as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(vsfc, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing vsfc as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(wsfc, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing wsfc as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(pmsl, nx, ny, 1, data2d)
      data2d = data2d*0.01
      print '(2a)', 'v5dout: writing pmsl*0.01 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(psfc, nx, ny, 1, data2d)
      data2d = data2d*0.01
      print '(2a)', 'v5dout: writing psfc*0.01 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(cldbase, nx, ny, 1, data2d)
      data2d = data2d*0.1
      print '(2a)', 'v5dout: writing cldbase*0.1 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(cldtop, nx, ny, 1, data2d)
      data2d = data2d*0.1
      print '(2a)', 'v5dout: writing cldtop*0.1 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(cldamt, nx, ny, 1, data2d)
      data2d = data2d*100.
      print '(2a)', 'v5dout: writing cldamt*100 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(ceiling, nx, ny, 1, data2d)
      data2d = data2d*0.1
      print '(2a)', 'v5dout: writing ceiling*0.1 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(intliqwater, nx, ny, 1, data2d)
      data2d = data2d
      print '(2a)', 'v5dout: writing intliqwater (mm) as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(totpcpwater, nx, ny, 1, data2d)
      data2d = data2d
      print '(2a)', 'v5dout: writing totpcpwater (mm) as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(pcp_inc, nx, ny, 1, data2d)
      data2d = data2d*39.37
      print '(2a)', 'v5dout: writing pcp_inc*39.37 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(pcp_tot, nx, ny, 1, data2d)
      data2d = data2d*39.37
      print '(2a)', 'v5dout: writing pcp_tot*39.37 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(con_pcp_inc, nx, ny, 1, data2d)
      data2d = data2d*39.37
      print '(2a)', 'v5dout: writing con_pcp_inc*39.37 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(con_pcp_tot, nx, ny, 1, data2d)
      data2d = data2d*39.37
      print '(2a)', 'v5dout: writing con_pcp_tot*39.37 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(snow_inc, nx, ny, 1, data2d)
      data2d = data2d*39.37
      print '(2a)', 'v5dout: writing snow_inc*39.37 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(snow_tot, nx, ny, 1, data2d)
      data2d = data2d*39.37
      print '(2a)', 'v5dout: writing snow_tot*39.37 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(pcptype_sfc, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing pcptype_sfc as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(thetasfc, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing thetasfc as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(thetaesfc, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing thetaesfc as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(cape, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing cape as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(cin, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing cin as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(liftedind, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing liftedind as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(srhel, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing srhel as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(refl_sfc, nx, ny, 1, data2d)
      where (refl_sfc .lt. 0) refl_sfc = 0.
      print '(2a)', 'v5dout: writing refl_sfc as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(max_refl, nx, ny, 1, data2d)
      where (data2d .lt. 0) data2d = 0.
      print '(2a)', 'v5dout: writing max_refl as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(echo_tops, nx, ny, 1, data2d)
      data2d = data2d*0.1
      print '(2a)', 'v5dout: writing echo_tops*0.1 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(heatind, nx, ny, 1, data2d)
      data2d = data2d - 273.15
      print '(2a)', 'v5dout: writing heatind - 273.15 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(visibility, nx, ny, 1, data2d)
      data2d = data2d*0.001
      print '(2a)', 'v5dout: writing visibility*0.001 as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(terdot, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing terdot as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(snowcover, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing snowcover as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(lwout, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing lwout as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(swout, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing swout as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(lwdown, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing lwdown as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(swdown, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing swsdown as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(shflux, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing shflux as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(lhflux, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing lhflux as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(pblhgt, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing pblhgt as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(upbl, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing upbl as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(vpbl, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing vpbl as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      call post2v5d(ground_t, nx, ny, 1, data2d)
      print '(2a)', 'v5dout: writing ground_t as varname: ', &
         varname(varindex)
      ret = v5dwrite(timestep, varindex, data2d)
      varindex = varindex + 1

      deallocate (data2d)
      return
   end subroutine v5dout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine v5dend

      implicit none
      include 'v5df.h'
      character(len=255) donefile
      integer :: flagunit
      ret = v5dclose()
      if (ret .eq. 0) then
         print '(a)', 'problem closing the vis5d data file!'
      end if
      if (realtime) then
         donefile = trim(v5dfile)//'.done'
         call get_file_unit(flagunit)
         open (unit=flagunit, file=donefile, status='unknown')
         close (flagunit)
      end if

   end subroutine v5dend
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine post2v5d(postdata, imax, jmax, kmax, v5ddata)

!-----------------------------------------------------------------------
!  transform data in post space (i increasing east, j increasing north,
!    k increasing upward)  into vis5d space (i increasing southward,
!    j increasing east, k increasing upward)
!-----------------------------------------------------------------------

      implicit none

      integer                      :: i
      integer                      :: imax
      integer                      :: j
      integer                      :: jmax
      integer                      :: k
      integer                      :: jmaxp1
      integer                      :: kmax
      real                         :: postdata(imax, jmax, kmax)
      real                         :: v5ddata(jmax, imax, kmax)

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

      jmaxp1 = jmax + 1

      do k = 1, kmax
         do j = 1, jmax
            do i = 1, imax

               v5ddata(j, i, k) = postdata(i, jmaxp1 - j, k)

            end do
         end do
      end do

   end subroutine post2v5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module vis5d
