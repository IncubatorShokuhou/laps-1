!dis    forecast systems laboratory
!dis    noaa/oar/erl/fsl
!dis    325 broadway
!dis    boulder, co     80303
!dis
!dis    forecast research division
!dis    local analysis and prediction branch
!dis    laps
!dis
!dis    this software and its documentation are in the public domain and
!dis    are furnished "as is."  the united states government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  they assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis    technical support to users.
!dis
!dis    permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  all modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  if significant modifications or enhancements
!dis    are made to this software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

subroutine gsi_bkobs

!==========================================================
!  this routine converts laps background and observations
!  to gsi format files.
!
!        background  --> wrf_inout;
!        observation --> bufr files, e.g., prepqc.
!
!  history:
!         creation: yuanfu xie        3-2006
!        modified: yuanfu xie        10-2007 adding height to wrf.
!==========================================================

   use laps_parm

   implicit none

   ! generate wrf_inout file: background.
   call gsi_bkg(n(1), n(2), n(3), lat, lon, &
                dxy, u_wind3d, v_wind3d)

   ! generate bufr files for observation: prepqc.laps
   call gsi_obs(i4time, asctime, lat, lon, n(1), n(2), n(3), &
                nobs_point, obs_point, n_tobs, obs_temp, maxtobs)
   !call gsi_obs_test

end subroutine gsi_bkobs

subroutine gsi_bkg(imax, jmax, kmax, xlat, xlong, &
                   grid_spacing, u_laps_bkg, v_laps_bkg)

!==========================================================
!  this routine converts laps background into gsi wrf_inout
!  on the gsi mass coordinate.
!
!  history: mar. 2006 by yuanfu xie.
!==========================================================

   use laps_parm

   implicit none

   integer, intent(in) :: imax, jmax, kmax          ! dimensions
   real, intent(in) :: xlat(imax, jmax), xlong(imax, jmax)
   real, intent(in) :: grid_spacing
   real, intent(in) :: u_laps_bkg(imax, jmax, kmax)   ! u bkg
   real, intent(in) :: v_laps_bkg(imax, jmax, kmax)   ! v bkg

   ! local variables:
   real, parameter :: cp = 1004.0, rc = 287.0, t0 = 300.0 !273.15
   character varname*3, fnm*9, hr*2, mins*2, jday*5, filename*150
   integer :: i4time_sys, namelen
   integer :: istatus, i, j, k
   real :: t_mass_bkg(imax, jmax, kmax)
   real :: sh_mass_bkg(imax, jmax, kmax)
   real :: geo_mass_bkg(imax, jmax, kmax)        ! geopotential
   real :: u_mass_bkg(imax, jmax, kmax), v_mass_bkg(imax, jmax, kmax)
   real :: dam(imax, jmax), pdam(imax, jmax)
   real :: znw(kmax), znu(kmax - 1), mapfac_m(imax, jmax)

   ! times:
   character*19 :: times

   ! system time:
   call get_systime_all(i4time_sys, fnm, hr, mins, asctime, jday, istatus)
   if (i4time .ne. i4time_sys) then
      print *, 'gsibkg: error: reading background in wrong time'
      stop
   end if

   ! dry air mass in column (base state):
   call dryairmass(dam, pdam, imax, jmax, kmax, pressr1d, &
                   height3d, pressr3d, temptr3d)

   ! as default, a uniform mass vertical grid is used:
   ! eta = (p_d(k)-p_d(top))/(p_d(sfc)-p_d(top)):
   ! do k=1,kmax
   !  znw(k) = 1.0-float(k-1)/float(kmax-1)
   ! enddo
   ! do k=1,kmax-1
   !   znu(k) = 1.0-(float(k-1)+0.5)/float(kmax-1)
   ! enddo
   do k = 1, kmax
      znw(k) = (pressr1d(k) - pressr1d(kmax))/(pressr1d(1) - pressr1d(kmax))
   end do
   do k = 1, kmax - 1
      znu(k) = (0.5*(pressr1d(k) + pressr1d(k + 1)) - pressr1d(kmax))/(pressr1d(1) - pressr1d(kmax))
   end do

   ! map factor: use laps routine: get_sigma(lat,lon,fac,istatus)
   ! mapfac_m = 1.0                ! test
   do j = 1, jmax
      do i = 1, imax
         call get_sigma(xlat(i, j), xlong(i, j), mapfac_m(i, j), istatus)
      end do
   end do

   ! t background:
   ! convert t to perturbation potential temperature (theta-t0):
   do k = 1, kmax
      temptr3d(1:imax, 1:jmax, k) = temptr3d(1:imax, 1:jmax, k)* &
                                    (100000.0/pressr1d(k))**(rc/cp) - t0
   end do
   call laps2mass(temptr3d, imax, jmax, kmax, pressr1d, dam, znw, 4, 0, t_mass_bkg)

   ! qvapor background:
   call laps2mass(sphumd3d, imax, jmax, kmax, pressr1d, dam, znw, 4, 0, sh_mass_bkg)

   ! u background:
   call laps2mass(u_laps_bkg, imax, jmax, kmax, pressr1d, dam, znw, 4, 0, u_mass_bkg)

   ! v background:
   call laps2mass(v_laps_bkg, imax, jmax, kmax, pressr1d, dam, znw, 4, 0, v_mass_bkg)

   ! height background:
   call laps2mass(height3d, imax, jmax, kmax, pressr1d, dam, znw, 4, 0, geo_mass_bkg)
   ! geopotential:
   geo_mass_bkg = 9.80665*geo_mass_bkg

   ! write out state variable for post-processing:
   call get_directory('log', filename, namelen)
   filename = filename(1:namelen)//'fort.12'
   open (unit=12, file=filename(1:namelen + 7), form='unformatted')
   write (12) znw, pressr1d, dam, znu
   close (12)

   times(1:4) = asctime(8:11)
   times(5:5) = '-'
   select case (asctime(4:6))
   case ('jan')
      times(6:7) = '01'
   case ('feb')
      times(6:7) = '02'
   case ('mar')
      times(6:7) = '03'
   case ('apr')
      times(6:7) = '04'
   case ('may')
      times(6:7) = '05'
   case ('jun')
      times(6:7) = '06'
   case ('jul')
      times(6:7) = '07'
   case ('aug')
      times(6:7) = '08'
   case ('sep')
      times(6:7) = '09'
   case ('oct')
      times(6:7) = '10'
   case ('nov')
      times(6:7) = '11'
   case ('dec')
      times(6:7) = '12'
   case default
      print *, 'gsib: error: invalid month: ', asctime(4:6)
      stop
   end select
   times(8:8) = '-'
   times(9:10) = asctime(1:2)
   times(11:11) = '_'
   times(12:13) = asctime(13:14)
   times(14:14) = ':'
   times(15:16) = asctime(15:16)
   times(17:19) = ':00'

   ! write the variables into a wrf_inout netcdf file:
   call wrfbkgout(times, imax, jmax, kmax, pressr1d(kmax), &
                  znu, znw, grid_spacing, mapfac_m, xlat, &
                  xlong, dam, pdam, t_mass_bkg, geo_mass_bkg, &
                  sh_mass_bkg, u_mass_bkg, v_mass_bkg, topo)

end subroutine gsi_bkg

subroutine gsi_obs(i4time, asctime, lat, lon, imax, jmax, kmax, &
                   nobs, obs_point, n_tobs, obs_temp, maxtobs)

!==========================================================
!  this routine generates a bufr format data file for gsi.
!
!  history: mar. 2006 by yuanfu xie.
!==========================================================

   ! gsiobs reads laps ingest observation data (obs_point)
   ! and converts it into prepbufr format saved in prepqc.laps
   ! so that gsi 3dvar can use.

   implicit none

   include 'barnesob.inc'

   character*16, intent(in) :: asctime
   integer, intent(in) :: i4time, nobs, imax, jmax, kmax
   real, intent(in) :: lat(imax, jmax), lon(imax, jmax)
   type(barnesob) :: obs_point(*)
   integer, intent(in) :: n_tobs, maxtobs
   real, intent(in) :: obs_temp(maxtobs, 12)

   integer, parameter :: mxmn = 8
   integer, parameter :: mxlv = 255
   real*8, parameter :: missing = 10.0e10
   real :: rlat, rlon, ztopsa, pres_3d(imax, jmax, kmax), p
   real*8 :: r8arr(mxmn, mxlv), rval

   integer, parameter :: mxbf = 16000
   integer :: ibfmsg(mxbf/4)

   integer :: yearofreport, mnthofreport, daysofreport
   integer :: hourofreport, istatus
   integer :: idate, nlv, nlvst, jj, libf, ierr, iobs, zero

   character        :: cval*8
   equivalence(cval, rval)

   ! static access:
   character :: dirstc*256, dir*200
   integer :: dirlen

   ! ingest obs:
   print *, 'total wind obs into gsi: ', nobs
   print *, 'total temp obs into gsi: ', n_tobs

   ! read background total pressure 3d field:
   call get_pres_3d(i4time, imax, jmax, kmax, pres_3d, istatus)

   ! open the bufr messages file.
   open (unit=11, file='prepqc.laps', form='unformatted')

   ! open the bufr tables file.
   call get_directory('static', dirstc, dirlen)
   dir = dirstc(1:dirlen)//'prepobs_prep.bufrtable'
   open (unit=12, file=dir(1:dirlen + 23))

   ! associate the tables file with the messages file, and
   ! identify the latter to the bufrlib software.

   call openbf(11, 'out', 12)

   ! report date:
   yearofreport = 2005
   mnthofreport = 11
   daysofreport = 21
   hourofreport = 15
   yearofreport = 2006
   mnthofreport = 2
   daysofreport = 7
   hourofreport = 14
   zero = ichar('0')
   yearofreport = (ichar(asctime(8:8)) - zero)*1000 + &
                  (ichar(asctime(9:9)) - zero)*100 + &
                  (ichar(asctime(10:10)) - zero)*10 + &
                  (ichar(asctime(11:11)) - zero)
   select case (asctime(4:6))
   case ('jan')
      mnthofreport = 1
   case ('feb')
      mnthofreport = 2
   case ('mar')
      mnthofreport = 3
   case ('apr')
      mnthofreport = 4
   case ('may')
      mnthofreport = 5
   case ('jun')
      mnthofreport = 6
   case ('jul')
      mnthofreport = 7
   case ('aug')
      mnthofreport = 8
   case ('sep')
      mnthofreport = 9
   case ('oct')
      mnthofreport = 10
   case ('nov')
      mnthofreport = 11
   case ('dec')
      mnthofreport = 12
   case default
      print *, 'wrfbkgout: error: invalid month: ', asctime(4:6)
      stop
   end select
   daysofreport = (ichar(asctime(1:1)) - zero)*10 + &
                  ichar(asctime(2:2)) - zero
   hourofreport = (ichar(asctime(13:13)) - zero)*10 + &
                  ichar(asctime(14:14)) - zero

   ! for every wind observation data:
   do iobs = 1, nobs

      idate = ((yearofreport)*1000000) + &
              ((mnthofreport)*10000) + &
              ((daysofreport)*100) + &
              (hourofreport)

      ! open a rawinsonde bufr message in order to store the new
      ! data subset (i.e. report).

      call openmb(11, 'adpupa', idate)

      ! store the report date-time within the data subset.

      r8arr(1, 1) = (yearofreport)
      r8arr(2, 1) = (mnthofreport)
      r8arr(3, 1) = (daysofreport)
      r8arr(4, 1) = (hourofreport)

      call ufbseq(11, r8arr, mxmn, 1, nlv, 'uartm')

      !   store the station identification information within the
      !   data subset.
      !   cval = ( station id, e.g. '72403', 'dbbh', etc.)
      cval = '72403'

      ! convert to lat/lon:
      call rlapsgrid_to_latlon(obs_point(iobs)%ri, obs_point(iobs)%rj, &
                               lat, lon, imax, jmax, rlat, rlon, ierr)

      r8arr(1, 1) = rval
      r8arr(3, 1) = rlat                !( station latitude )
      r8arr(2, 1) = rlon                !( station longitude )
      r8arr(4, 1) = (obs_point(iobs)%i4time - i4time)/3600.0 !( obs time )
      r8arr(5, 1) = 210.0                !( prepbufr report type )
      r8arr(6, 1) = obs_point(iobs)%elev!( station elevation )

      call ufbint(11, r8arr, mxmn, 1, nlv, &
                  'sid xob yob dhr typ elv ')

      ! store the level data within the data subset.

      ! for laps data ingest, treat obs individually for now:
      nlvst = 1         !( number of data levels to be stored)

      do jj = 1, nlvst

         ! use trilinear_laps.f to interpolate pressure value at ri,rj,rk:
         call trilinear_laps(obs_point(iobs)%ri, &
                             obs_point(iobs)%rj, &
                             obs_point(iobs)%rk, imax, jmax, kmax, &
                             pres_3d, p)
         r8arr(1, jj) = p/100.0        ! (in mb)

         ! r8arr (1,jj) = ztopsa(obs_point(iobs)%elev)

         !899.0 !missing !( pressure, in mb, for level jj )
         r8arr(2, jj) = missing !( specific humidity observation)
         r8arr(3, jj) = missing !( temperature in c for level jj )
         r8arr(4, jj) = missing !( height, in m, for level jj )
         r8arr(5, jj) = obs_point(iobs)%value(1) !( u m/s, for level jj )
         r8arr(6, jj) = obs_point(iobs)%value(2) !( v m/s, for level jj )
         r8arr(7, jj) = missing  !( precipitable water mm, for level jj )
         r8arr(8, jj) = missing

      end do

      call ufbint(11, r8arr, mxmn, nlvst, nlv, &
                  'pob qob tob zob uob vob pwo cat')
      r8arr = 0.0
      call ufbint(11, r8arr, mxmn, nlvst, nlv, &
                  'pqm qqm tqm zqm wqm nul pwq')

      ! ( store any other available values in a similar manner )
      ! once all data values have been stored for this data subset,
      ! we are now ready to store the data subset into the message.

      call writsa(11, ibfmsg, libf)

   end do

   ! for every temperature observation data:
   do iobs = 1, n_tobs

      idate = ((yearofreport)*1000000) + &
              ((mnthofreport)*10000) + &
              ((daysofreport)*100) + &
              (hourofreport)

      ! open a rawinsonde bufr message in order to store the new
      ! data subset (i.e. report).

      call openmb(11, 'adpupa', idate)

      ! store the report date-time within the data subset.

      r8arr(1, 1) = (yearofreport)
      r8arr(2, 1) = (mnthofreport)
      r8arr(3, 1) = (daysofreport)
      r8arr(4, 1) = (hourofreport)

      call ufbseq(11, r8arr, mxmn, 1, nlv, 'uartm')

      !   store the station identification information within the
      !   data subset.
      !   cval = ( station id, e.g. '72403', 'dbbh', etc.)
      cval = '72403'

      ! convert to lat/lon:
      call rlapsgrid_to_latlon(obs_temp(iobs, 1), obs_temp(iobs, 2), &
                               lat, lon, imax, jmax, rlat, rlon, ierr)

      r8arr(1, 1) = rval
      r8arr(3, 1) = rlat                !( station latitude )
      r8arr(2, 1) = rlon                !( station longitude )
      r8arr(4, 1) = 0.0                 !( obs time; enhance this later)
      r8arr(5, 1) = 132.0                !( prepbufr report type )
      r8arr(6, 1) = missing !( station elevation )

      call ufbint(11, r8arr, mxmn, 1, nlv, &
                  'sid xob yob dhr typ elv ')

      ! store the level data within the data subset.

      ! for laps data ingest, treat obs individually for now:
      nlvst = 1         !( number of data levels to be stored)

      do jj = 1, nlvst

         ! use trilinear_laps.f to interpolate pressure value at ri,rj,rk:
         call trilinear_laps(obs_temp(iobs, 1), &
                             obs_temp(iobs, 2), &
                             obs_temp(iobs, 3), imax, jmax, kmax, &
                             pres_3d, p)
         r8arr(1, jj) = p/100.0        ! (in mb)

         ! r8arr (1,jj) = ztopsa(obs_point(iobs)%elev)

         !899.0 !missing !( pressure, in mb, for level jj )
         r8arr(2, jj) = missing !( specific humidity observation)
         r8arr(3, jj) = obs_temp(iobs, 7) !( temperature in c for level jj )
         r8arr(4, jj) = missing !( height, in m, for level jj )
         r8arr(5, jj) = missing !( u m/s, for level jj )
         r8arr(6, jj) = missing !( v m/s, for level jj )
         r8arr(7, jj) = missing  !( precipitable water mm, for level jj )
         r8arr(8, jj) = missing

      end do

      call ufbint(11, r8arr, mxmn, nlvst, nlv, &
                  'pob qob tob zob uob vob pwo cat')
      r8arr = 0.0
      call ufbint(11, r8arr, mxmn, nlvst, nlv, &
                  'pqm qqm tqm zqm wqm nul pwq')

      ! ( store any other available values in a similar manner )
      ! once all data values have been stored for this data subset,
      ! we are now ready to store the data subset into the message.

      call writsa(11, ibfmsg, libf)

   end do

   ! forcibly flush the last bufr message, if any, from the
   ! bufrlib software.

   call writsa(-11, ibfmsg, libf)
   call closbf(11)

   print *, 'wind data converted: ', nobs
   print *, 'temp data converted: ', n_tobs

end subroutine gsi_obs
