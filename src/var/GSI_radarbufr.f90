       subroutine gsi_radarbufr
!***********************************************************************
! program name : gsi_radarbufr(gsi_radarbufr.f)
!
! description: to read in laps's level ii polar radar data in netcdf format,
!              and write out ncep's gsi level iii radar bufr data.
!
! called function:
!    get_systime(${laps_src_root}/src/lib/get_systime.f)
!    get_radarbufr_parms
!    get_rfile_time
!    write_radvel_bufr
!    write_radref_bufr
!
! date :
!   original     -- may. 04, 2006 (shiow-ming, deng)
!***********************************************************************

          parameter(maxradars=99, maxtimes=999)
          character*9 a9_time
          character*150 path_lvl2(maxradars), path_table, path_output
          character*150 file_table, file_output, path
          character*9 cdf_times(maxtimes), cdfdate
          integer num_elevs(maxtimes), ielev(30, maxtimes)
          integer jelev(maxtimes)

!-----------------------------------------------------------------------
!c  to get laps's system time.

          call get_systime(i4time, a9_time, istatus)
          if (istatus .ne. 1) go to 901

!-----------------------------------------------------------------------
!c  to get paths of radar bufr table, laps netcdf level 2 radar data.

          call get_radarbufr_parms(nradars, path_table, path_output &
                                   , path_lvl2, itimeb, itimea, istatus)
          if (istatus .ne. 1) go to 902
          print *, 'get paths of radar bufr table, and '
          print *, '    the paths of laps netcdf level 2 radar data.'
          print *, 'total number of radar: ', nradars

!-----------------------------------------------------------------------
!c  to open radar bufr table file: "radar.bufrtable"

          lent = index(path_table, ' ')
          file_table = path_table(1:lent - 1)//'/radar.bufrtable'
          open (12, file=file_table, err=903)

!-----------------------------------------------------------------------
!c  to get output radar level iii bufr data file name: "radarbufr"

          len = index(path_output, ' ')
          file_output(1:len - 1) = path_output(1:len - 1)
          file_output(len:len + 9) = '/radarbufr'
          file_output(len + 10:lev + 10) = char(0)

          icheck = 0
          do 100 i = 1, nradars

!-----------------------------------------------------------------------
!c  to get laps radar polar netcdf data time.

             print *, ' '
             path = path_lvl2(i)
             print *, 'the input radar path: ', path
             call get_rfile_time(a9_time, itimeb, itimea, path, num_times &
                                 , cdf_times, num_elevs, ielev, istatus)
             if (istatus .ne. 1) go to 100
             icheck = icheck + 1

!-----------------------------------------------------------------------
!c  to open radar bufr data.

             if (icheck .eq. 1) then
!             open(11,file=file_output,form='unformatted')
                open (11, file='radarbufr', form='unformatted')
                call openbf(11, 'out', 12)
                lunin = 11
             end if

!-----------------------------------------------------------------------
!c  to write out bufr data for radial velocity field.

             len = index(path, ' ')
             len = len - 1
             do k = 1, num_times
                cdfdate(1:9) = cdf_times(k) (1:9)
                numelev = num_elevs(k)
                do j = 1, numelev
                   jelev(j) = ielev(j, k)
                end do
                call write_radvel_bufr(lunin, path, len, a9_time, cdfdate &
                                       , numelev, jelev, istatus)
             end do
100          continue
             if (icheck .eq. 0) go to 904

!-----------------------------------------------------------------------
!c  to force bufr data to be closed.

!      call writsa(-11,ibfmsg,libf)
             call closbf(11)

             return
901          continue
             print *, 'cannot get system time.'
             return
902          continue
             print *, 'cannot read in parameter of namelist file: '
             print *, 'please check file: static/radar_bufr.nl'
             return
903          continue
             print *, 'cannot open radar bufr table.'
             print *, 'radar bufr table file: ', file_table(1:lent + 15)
             return
904          continue
             print *, 'no radar data can be written to bufr.'
             return

          end

          subroutine write_radvel_bufr(lunin, directory, len, a9_time, cdfdate &
                                       , num_elev, ielev, istatus)
!***********************************************************************
! function name : write_radvel_bufr
!
! usage : call write_radvel_bufr(lunin,directory,len,a9_time,cdfdate
!        1                      ,num_elev,ielev,istatus)
!
! description      : to write out level iii radar's radial velocity data
!                    in ncep bufr format.
!
! arguments :
!  i/o/w   name,      type,       description
!    i     lunin      integer     absolute value is fortran logical unit
!                                 number for bufr file.
!    i     directory  c*150       input directory of netcdf level ii radar data.
!    i     len        integer     the character's length of directory.
!    i     a9_time    c*9         system time. (= yyjjjhhmm )
!    i     cdfdate    c*9         the date: "yyjjjhhmm".
!                                 yy: year, jjj: julian day, hh: hour, mm: minute.
!    i     num_elev   integer     number of elevation angle.
!    i ielev(num_elev) int array  elevation angle id.
!    o     istatus    integer     work or no work message,
!                                 =0, the subroutine doesnot work.
!                                 =1, the subroutine does work.
!
! sub./fun. called :
!    verjulian,
!    read_radar_polar
!    int_ppi,
!    openmb,  ufbint,  writsa  (ncep bufr api)
!
! date :
!   original     -- apr. 12, 2006 (shiow-ming, deng)
!***********************************************************************
             parameter(ngates=1600, nbeams=420, ix=401)
             parameter(mxbf=16000)

             character directory*150, a9_time*9, cdfdate*9
             integer ielev(num_elev)

             character radarfile*150, ca2*2
             character radarname*5
             real azim(nbeams), elev(nbeams)
             integer*8 itim(nbeams)

             integer iz(ngates*nbeams), iv(ngates*nbeams), iw(ngates*nbeams)
             integer it(ngates*nbeams)
             integer*8 istarttime, isec

             real xd(ix)
             integer ndatv(ix, ix), ndatw(ix, ix), ndatt(ix, ix)

             real*8 hdr(10), rwnd(7, ngates)
             character*4 sstn
             equivalence(sstn, hdr(1))
             character*8 subset
             character cout*49, time*12
             integer ibfmsg(mxbf/4)

             istatus = 0
             cout(1:49) = 'rpid clat clon selv anel year mnth days hour minu'
             i_missing = 32767
             pi = acos(-1.)
             earth = 6371.25
             rmax = 200.
             nx = ix
             do i = 1, nx
                xd(i) = -rmax + (i - 1.)
             end do

!-----------------------------------------------------------------------
!c  to get netcdf file names.

             if (directory(len:len) .eq. '/') then
                len = len - 1
             end if

!-----------------------------------------------------------------------
!c  to get idate from year,month,day,and hour.

             read (a9_time, '(i2,i3,2i2)') iyear, julian, ihour, minute
             iyear = iyear + 2000
             call verjulian(iyear, julian, month, iday, ier)
             if (ier .ne. 0) return
             idate = iyear*1000000 + month*10000 + iday*100 + ihour
             read (cdfdate, '(i2,i3,2i2)') iyear, julian, ihour, minute
             iyear = iyear + 2000
             call verjulian(iyear, julian, month, iday, ier)
             if (ier .ne. 0) return
             write (time, '(i4,4i2)') iyear, month, iday, ihour, minute

!-----------------------------------------------------------------------
!c  to open bufr message.

             subset(1:8) = 'nc006001'
             call openmb(lunin, subset, idate)

             icheck = 0
             do 100 k = 1, num_elev

!-----------------------------------------------------------------------
!c  to read in radar's data from netcdf polor-coordinated file.

                write (ca2, '(i2)') ielev(k)
                if (k .lt. 10) ca2(1:1) = '0'
                radarfile(1:len) = directory(1:len)
                radarfile(len + 1:len + 1) = '/'
                radarfile(len + 2:len + 10) = cdfdate(1:9)
                radarfile(len + 11:len + 15) = '_elev'
                radarfile(len + 16:len + 17) = ca2(1:2)
                radarfile(len + 18:len + 18) = char(0)
                call read_radar_polar(radarfile, i_missing &
                                      , radarname, slat, slon, salt, istarttime &
                                      , inumber, elevation, nobeam, nogatez, nogatev &
                                      , rangez, rangev, spacez, spacev, snyquist &
                                      , azim, elev, itim, iz, iv, iw, istat)
                if (istat .eq. 0) go to 100

!-----------------------------------------------------------------------
!c  to interpolate data for ppi scan.

                range = 1000.*rangev
                space = 1000.*spacev
                call int_ppi(nogatev, nobeam, iv, i_missing, rmax, range, space &
                             , elevation, azim, nx, xd, ndatv)
                ickw = 0
                do i = 1, nobeam*nogatev
                   if (iw(i) .ne. i_missing) ickw = 1
                end do
                if (ickw .eq. 1) then
                   call int_ppi(nogatev, nobeam, iw, i_missing, rmax, range, space &
                                , elevation, azim, nx, xd, ndatw)
                else
                   do j = 1, nx
                   do i = 1, nx
                      ndatw(i, j) = i_missing
                      if (ndatv(i, j) .ne. i_missing) ndatw(i, j) = 250
                   end do
                   end do
                end if
                ij = 0
                do j = 1, nobeam
                do i = 1, nogatev
                   ij = ij + 1
                   it(ij) = (itim(j) - istarttime)*10/6.
                end do
                end do
                call int_ppi(nogatev, nobeam, it, i_missing, rmax, range, space &
                             , elevation, azim, nx, xd, ndatt)

                idy = istarttime/86400
                isec = istarttime - 86400*idy
                ii = isec/60
                isecst = isec - 60*ii
                ih = ii/60
                imst = ii - 60*ih

                icheck = 1
                sstn(1:4) = radarname(1:4)
                hdr(2) = slat
                hdr(3) = slon
                hdr(4) = salt
                hdr(5) = elevation
                hdr(6) = iyear
                hdr(7) = month
                hdr(8) = iday
                hdr(9) = ih
                hdr(10) = imst

                ct = tan(elevation*pi/180.)
                nn = 0
                nnn = 0
                do 50 j = 1, nx
                do 50 i = 1, nx
                   if (ndatv(i, j) .eq. i_missing) go to 50
                   x = xd(i)
                   y = xd(j)
                   r = sqrt(x*x + y*y)
                   if (r .lt. 0.1) go to 50
                   plat = y/earth*180./pi
                   platm = slat + 0.5*plat
                   earth1 = earth*cos(platm*pi/180.)
                   ht = 1000*r*ct + r**2/17.
                   th = acos(y/r)*180./pi
                   if (x .lt. 0.) th = 360.-th
                   nn = nn + 1
                   nnn = nnn + 1
                   rwnd(1, nn) = 0.01*ndatt(i, j)
                   rwnd(2, nn) = slat + plat
                   rwnd(3, nn) = slon + x/earth1*180./pi
                   rwnd(4, nn) = salt + ht
                   rwnd(5, nn) = 0.01*ndatv(i, j)
                   rwnd(6, nn) = th
                   rwnd(7, nn) = 0.01*ndatw(i, j)
                   if (nn .eq. 810) then
                      call ufbint(lunin, hdr, 10, 1, levs, cout)
                      call ufbint(lunin, rwnd, 7, nn, levs &
                                  , 'stdm suplat suplon heit rwnd rwaz rstd')
                      call writsa(lunin, ibfmsg, ibf)
                      nn = 0
                   end if
50                 continue
                   if (nn .gt. 4) then
                      call ufbint(lunin, hdr, 10, 1, levs, cout)
                      call ufbint(lunin, rwnd, 7, nn, levs &
                                  , 'stdm suplat suplon heit rwnd rwaz rstd')
                      call writsa(lunin, ibfmsg, ibf)
                   end if
                   if (nnn .gt. 4) then
                      print *, 'output time: ', time(1:12) &
                         , '   sweeep number of input data: ', ca2(1:2) &
                         , '   total number of output data: ', nnn
                   end if
100                continue
                   if (icheck .eq. 0) return

                   istatus = 1
                   return
                end

                subroutine int_ppi(ng, nb, idat, miss, rmax, range, space, ang, azim &
                                   , nx, xd, ndat)
!***********************************************************************
! function name: int_ppi
!
! usage :
!    call int_ppi(ng,nb,idat,miss,rmax,range,space,ang,azim
!   1            ,nx,xd,ndat)
!
! description      : using the nearest point method intepolate
!                    the ppi dbz data  to grid point of  radars
!
! arguments :
!  i/o/w   name,      type,       description
!    i     ng         integer     gate number in one ray.
!    i     nb         integer     total ray number of input sweep scan layer.
!    i    idat(ng,nb) int array   input data.
!    i     miss       integer     missing or bad value.
!    i     rmax       real        the max. radius in kilometers.
!    i     range      real        range to first gate. (meters)
!    i     space      real        space between gates. (meters)
!    i     ang        real        sweep angle of input sweep scan layer,
!                                 (degrees)
!    i     azim(nb)   real array  azimuthal angle of a ray in one sweep.
!                                 (degrees)
!    i     nx         integer     the dimension of domain.
!    i     xd(nx)     real array  the x- or y-grid in unit of km.
!    o    ndat(nx,nx) int array   the output data.
!***********************************************************************

                   dimension idat(ng, nb)
                   dimension azim(nb), xd(nx)
                   dimension ndat(nx, nx)
                   dimension iwk(720), wk(720)

                   rspace = 0.001*space
                   rmin = 0.001*range
                   rmin1 = max(rmin, 0.1)
                   pi = acos(-1.)
                   rcos = cos(ang*pi/180.)
                   angg = max(ang, 0.1)
                   rrmax = 0.001*(range + (ng - 1)*space)

                   do i = 1, nb
                      wk(i) = azim(i)
                   end do
                   wk(nb + 1) = azim(1)
                   do 10 i = 1, 720
                      iwk(i) = -1
                      th = 0.5*i
                      do k = 1, nb
                         diff = abs(wk(k + 1) - wk(k))
                         diff1 = wk(k + 1) - th
                         diff2 = wk(k) - th
                         if (diff1*diff2 .le. 0.) then
                            if (diff .lt. 10.) then
                               if (abs(diff1) .lt. abs(diff2)) then
                                  iwk(i) = k + 1
                                  if (iwk(i) .gt. nb) iwk(i) = 1
                               else
                                  iwk(i) = k
                               end if
                               go to 10
                            end if
                         else
                            if (diff .gt. 355.) then
                               dd1 = abs(diff1)
                               dd2 = abs(diff2)
                               if (dd1 .gt. dd2) then
                                  dd1 = abs(360.-dd1)
                               else
                                  dd2 = abs(360.-dd2)
                               end if
                               if (dd1 .lt. dd2) then
                                  iwk(i) = k
                               else
                                  iwk(i) = k + 1
                                  if (iwk(i) .gt. nb) iwk(i) = 1
                               end if
                               go to 10
                            end if
                         end if
                      end do
10                    continue

                      do 100 j = 1, nx
                      do 100 i = 1, nx
                         ndat(i, j) = miss
                         rr = sqrt(xd(i)**2 + xd(j)**2)
                         rr1 = rr/rcos
                         if ((rr1 .lt. rmin1) .or. (rr1 .ge. rrmax)) go to 100
                         irad = (rr1 - rmin)/rspace + 1
                         if (irad .lt. 1) go to 100
                         th = acos(xd(j)/rr)*180./pi
                         if (xd(i) .lt. 0.) th = 360.-th
                         ith = 2.*th
                         if (ith .lt. 1) then
                            ith = 720
                            if (th .gt. 0.25) ith = 1
                         end if
                         iazi = iwk(ith)
                         if (iazi .gt. 0) ndat(i, j) = idat(irad, iazi)
100                      continue

                         return
                      end

                      subroutine verjulian(iyear, julian, month, iday, ier)
!***********************************************************************
! function name : verjulian
!
! usage : call verjulian(iyear,julian,month,iday,ier)
!
! description      : to get month, day from year, and julian day.
!
! arguments :
!  i/o/w   name,      type,       description
!    i     iyear      integer     year, ex: 1998
!    i     julian     integer     the julian day.
!                                 =1, for jan. 1st.
!    o     month      integer     month.
!    o     iday       integer     day
!    o     ier        integer     =0, success.
!                                 =1, failure.
!
! sub./fun. called : none
!
! date :
!   original     -- apr. 05, 2006 (shiow-ming, deng)
!***********************************************************************

                         dimension idate1(12), idate2(12)
                         data idate1/31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/
                         data idate2/31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366/

                         ier = 0
                         month = -1
                         iday = -1
                         if ((iyear .le. 0) .or. (julian .lt. 1) .or. (julian .gt. 366)) then
                            ier = 1
                            return
                         end if

                         ind = 0
                         i1 = iyear/4
                         i4 = iyear - 4*i1
                         if (i4 .eq. 0) ind = 1
                         i1 = iyear/400
                         i2 = iyear - 400*i1
                         if (i2 .eq. 0) ind = 1
                         i1 = iyear/100
                         i3 = iyear - 100*i1
                         if ((i3 .eq. 0) .and. (i2 .ne. 0)) ind = 0
                         if ((ind .eq. 0) .and. (julian .gt. 365)) then
                            ier = 1
                            return
                         end if

                         if (ind .eq. 0) then
                            do i = 1, 12
                               if (julian .le. idate1(i)) then
                                  month = i
                                  if (month .eq. 1) then
                                     iday = julian
                                  else
                                     iday = julian - idate1(i - 1)
                                  end if
                                  return
                               end if
                            end do
                         else
                            do i = 1, 12
                               if (julian .le. idate2(i)) then
                                  month = i
                                  if (month .eq. 1) then
                                     iday = julian
                                  else
                                     iday = julian - idate2(i - 1)
                                  end if
                                  return
                               end if
                            end do
                         end if
                         return
                      end

                      subroutine get_rfile_time(a9_time, itimeb, itimea, path, num_times &
                                                , cdf_times, num_elevs, ielev, istatus)
!***********************************************************************
! subroutine/function : get_rfile_time
!
! usage :
!      call get_rfile_time(a9_time,itimeb,itimea,path,num_times
!     1                   ,cdf_times,num_elevs,ielev,istatus)
!
! description      : to get netcdf level ii radar file names.
!
! arguments :
!  i/o/w   name,      type,       description
!    i     a9_time    c*9         system time. (= yyjjjhhmm )
!    i     itimeb     integer     time window before.
!    i     itimea     integer     time window after.
!                                 set these variables for time filtering the data.
!                                 units are seconds.  for example, if you want data with an
!                                 observation time from 60 min before to 30 min after the analysis
!                                 time to be included in the bufr file, use 3600 and 1800 for
!                                 itimeb and itimea, respectively.
!    i     path       c*150       full path to directory containing a set of input
!                                 radar tilts/volumes.
!    o     num_times  integer     number of time during system_time - itimeb and
!                                                       system_time + itimea.
!    o  cdf_times(999) c*9 array  time to write out bufr table. ( = yyjjjhhmm )
!    o  num_elevs(999) int array  number of elevation in a volue w.r.t. cdf_time.
!    o  ielevs(30,999) int array  elevation id in a volue w.r.t. cdf_time.
!    o     istatus    integer     the work or not message.
!                                 =1, read in data.
!                                 =0, didnot read any data.
!
! modules called :
!   getfilenames(getiofile.c)
!   difference_time
!
! date :
!   original     -- apr. 12, 2006 (shiow-ming, deng)
!***********************************************************************
                         parameter(maxtimes=999)
                         character a9_time*9, path*150
                         character*9 cdf_times(maxtimes)
                         integer num_elevs(maxtimes), ielev(30, maxtimes)

                         character directory*150
                         character*900000 cdat
                         character*16 file(9900), ftemp
                         integer iwork(30)
                         integer isyear, isjulian, ishour, isminute
                         integer iyear, julian, ihour, minute, jsec, ier

                         istatus = 0
                         num_times = 0
                         do j = 1, maxtimes
                            cdf_times(j) = ' '
                            num_elevs(j) = 0
                            do i = 1, 30
                               ielev(i, j) = 0
                            end do
                         end do

                         read (a9_time, '(i2,i3,2i2)', err=901) isyear, isjulian, ishour &
                            , isminute
                         len = index(path, ' ')
                         if (len .lt. 1) go to 902
                         directory(1:len - 1) = path(1:len - 1)
                         directory(len:len) = char(0)

                         isize = 9000
                         do i = 1, isize
                            cdat(i:i) = char(0)
                         end do
                         call getfilenames(directory, isize, cdat, ier)
                         if (ier .ne. 0) go to 903

                         numfile = 0
                         do i = isize, 1, -1
                            if (cdat(i:i) .ne. char(0)) then
                               ii = i
                               go to 10
                            end if
                         end do
                         return
10                       continue
                         i1 = 1
                         j = 1
20                       continue
                         j = j + 1
                         if (j .gt. ii) go to 30
                         if (cdat(j:j) .eq. ' ') then
                            i2 = j - 1
                            jj = i2 - i1 + 1
                            if ((jj .eq. 16) .and. (cdat(i1 + 9:i1 + 13) .eq. '_elev')) then
                               numfile = numfile + 1
                               file(numfile) (1:16) = cdat(i1:i2)
                            end if
                            i1 = j + 1
                         else
                            if (j .eq. ii) then
                               i2 = ii
                               jj = i2 - i1 + 1
                               if ((jj .eq. 16) .and. (cdat(i1 + 9:i1 + 13) .eq. '_elev')) then
                                  numfile = numfile + 1
                                  file(numfile) (1:16) = cdat(i1:i2)
                               end if
                               go to 30
                            end if
                         end if
                         go to 20
30                       continue
                         if (numfile .eq. 0) go to 904

                         do 50 k = 1, numfile
                            ftemp(1:16) = file(k) (1:16)
                            read (ftemp, '(i2,i3,2i2)', err=50) iyear, julian, ihour, minute
                            call difference_time(isyear, isjulian, ishour, isminute &
                                                 , iyear, julian, ihour, minute, jsec, ier)
                            if (ier .ne. 0) go to 50
                            if (jsec .gt. itimea) go to 50
                            if (abs(jsec) .gt. itimeb) go to 50
                            if (num_times .eq. 0) then
                               num_times = num_times + 1
                               cdf_times(num_times) (1:9) = ftemp(1:9)
                               go to 50
                            end if
                            do i = 1, num_times
                               if (cdf_times(i) (1:9) .eq. ftemp(1:9)) go to 50
                            end do
                            num_times = num_times + 1
                            cdf_times(num_times) (1:9) = ftemp(1:9)
50                          continue
                            if (num_times .eq. 0) go to 905

                            do 100 k = 1, num_times
                               do i = 1, 30
                                  iwork(i) = 0
                               end do
                               id = 0
                               do 60 i = 1, numfile
                                  ftemp(1:16) = file(i) (1:16)
                                  if (cdf_times(k) (1:9) .eq. ftemp(1:9)) then
                                     read (ftemp, '(14x,i2)', err=60) iva
                                     num_elevs(k) = num_elevs(k) + 1
                                     id = id + 1
                                     iwork(id) = iva
                                  end if
                                  if (num_elevs(k) .ge. 30) go to 70
60                                continue
70                                continue
                                  if (id .eq. 0) go to 100
                                  if (id .eq. 1) then
                                     ielev(1, k) = iwork(1)
                                     go to 100
                                  end if
                                  do i = 1, id - 1
                                     do j = i + 1, id
                                        if (iwork(i) .gt. iwork(j)) then
                                           iw1 = iwork(i)
                                           iwork(i) = iwork(j)
                                           iwork(j) = iw1
                                        end if
                                     end do
                                  end do
                                  do i = 1, id
                                     ielev(i, k) = iwork(i)
                                  end do
100                               continue

                                  istatus = 1
                                  return
901                               continue
                                  print *, 'cannot read system time a9_time(yyjjjhhmm): '
                                  print *, 'a9_time: ', a9_time(1:9)
                                  return
902                               continue
                                  print *, 'the path of radar level ii file does not exit.'
                                  print *, 'path: ', path(1:150)
                                  return
903                               continue
                                  print *, 'cannot read in the file names in directory: ' &
                                     , path(1:len)
                                  return
904                               continue
                                  print *, 'no any radar level ii file can be found.'
                                  print *, 'path: ', path(1:len)
                                  return
905                               continue
                                  print *, 'no any time data can be found.'
                                  print *, 'path: ', path(1:len)
                                  print *, 'system time: ', a9_time(1:9)
                                  print *, 'itimeb: ', itimeb, ' itimea: ', itimea
                                  return

                               end

                               subroutine difference_time(isyear, isjulian, ishour, isminute &
                                                          , iyear, julian, ihour, minute, jsec, ier)
!***********************************************************************
! subroutine/function : difference_time
!
! usage :
!      call difference_time(isyear,isjulian,ishour,isminute
!     1                    ,iyear,julian,ihour,minute,jsec,ier)
!
! description      : to compute diference second between
!                    target time and system time.
!
! arguments :
!  i/o/w   name,      type,       description
!    i     isyear     integer     system year.
!    i     isjulian   integer     system julian day.
!    i     ishour     integer     system hour.
!    i     isminute   integer     system minute.
!    i     iyear      integer     target year.
!    i     julian     integer     target julian day.
!    i     ihour      integer     target hour.
!    i     minute     integer     target minute.
!    o     jsec       integer     second between target time and system time.
!                                 =999999 for too large second. (ier=1)
!    o     ier        integer     error message.
!                                 =1, failure.
!                                 =0, success.
!
! modules called : none
!
! date :
!   original     -- apr. 12, 2006 (shiow-ming, deng)
!***********************************************************************

                                  integer isyear, isjulian, ishour, isminute
                                  integer iyear, julian, ihour, minute, jsec, ier

                                  ier = 1
                                  jsec = 999999
                                  jyr = iyear - isyear
                                  if (abs(jyr) .gt. 1) return
                                  if (jyr .eq. 0) then
                                     jday = julian - isjulian
                                  else
                                     if (jyr .eq. 1) then
                                        if (julian .ge. isjulian) return
                                        ind = 0
                                        i1 = isyear/4
                                        i4 = isyear - 4*i1
                                        if (i4 .eq. 0) ind = 1
                                        i1 = isyear/400
                                        i2 = isyear - 400*i1
                                        if (i2 .eq. 0) ind = 1
                                        i1 = isyear/100
                                        i3 = isyear - 100*i1
                                        if ((i3 .eq. 0) .and. (i2 .ne. 0)) ind = 0
                                        if (ind .eq. 0) then
                                           jday = julian + 365 - isjulian
                                        else
                                           jday = julian + 366 - isjulian
                                        end if
                                     else
                                        if (isjulian .ge. julian) return
                                        ind = 0
                                        i1 = iyear/4
                                        i4 = iyear - 4*i1
                                        if (i4 .eq. 0) ind = 1
                                        i1 = iyear/400
                                        i2 = iyear - 400*i1
                                        if (i2 .eq. 0) ind = 1
                                        i1 = iyear/100
                                        i3 = iyear - 100*i1
                                        if ((i3 .eq. 0) .and. (i2 .ne. 0)) ind = 0
                                        if (ind .eq. 0) then
                                           jday = julian - 365 - isjulian
                                        else
                                           jday = julian - 366 - isjulian
                                        end if
                                     end if
                                  end if
                                  jsec = 86400*jday + 3600*(ihour - ishour) + 60*(minute - isminute)
                                  ier = 0
                                  return
                               end

                               subroutine get_radarbufr_parms(nradars, path_table, path_output &
                                                              , path_lvl2, itimeb, itimea, istatus)
!***********************************************************************
! subroutine/function : get_radarbufr_parms
!
! usage :
!      call get_radarbufr_parms(nradars,path_table,path_output
!     1                        ,path_lvl2,itimeb,itimea,istatus)
!
! description      : to get parameters of namelist file:
!                    'static/radar_bufr.nl' for radar gsi bufr data.
!
! arguments :
!  i/o/w   name,      type,       description
!    o     nradars    integer     number of radars (and/or radar types) to loop through
!                                 and process.
!    o     path_table c*150       directory for input bufr table.
!    o    path_output c*150       directory for output bufr data.
!    o  path_lvl2(99) c*150 array full path to each directory containing a set of input
!                                 radar tilts/volumes.
!    o     itimeb     integer     time window before.
!    o     itimea     integer     time window after.
!                                 set these variables for time filtering the data.
!                                 units are seconds.  for example, if you want data with an
!                                 observation time from 60 min before to 30 min after the analysis
!                                 time to be included in the bufr file, use 3600 and 1800 for
!                                 itimeb and itimea, respectively.
!    o     istatus    integer     the work or not message.
!                                 =1, read in data.
!                                 =0, didnot read any data.
!
! modules called :
!   get_directory($laps_src_root/src/lib/get_directory.f)
!
! date :
!   original     -- may. 05, 2006 (shiow-ming, deng)
!***********************************************************************
                                  integer maxradars
                                  parameter(maxradars=99)
                                  character*150 path_table, path_output, path_lvl2(maxradars)

                                  character*150 path_to_radar_a(maxradars), path_to_vrc_nl
                                  character*4 laps_radar_ext_a(maxradars)
                                  logical l_line_ref_qc, l_hybrid_first_gate, l_unfold
                                  character*3 ext

                                  namelist /remap_nl/ n_radars_remap, max_times, path_to_radar_a &
                                     , laps_radar_ext_a, path_to_vrc_nl &
                                     , ref_min, min_ref_samples, min_vel_samples, dgr &
                                     , abs_vel_min, l_line_ref_qc, l_hybrid_first_gate &
                                     , l_unfold
                                  character*150 dir, filename

                                  nradars = 0
                                  istatus = 0

                                  call get_directory('nest7grid', dir, len)
                                  if (dir(len:len) .eq. '/') then
                                     path_table(1:len - 1) = dir(1:len - 1)
                                     path_table(len:len) = ' '
                                  else
                                     path_table(1:len) = dir(1:len)
                                  end if
                                  itimeb = 3600
                                  itimea = 3600

                                  filename = dir(1:len)//'/remap.nl'
                                  open (1, file=filename, status='old', err=900)
                                  read (1, remap_nl, err=901)
                                  close (1)

                                  nradars = 0
                                  do i = 1, n_radars_remap
                                     ext(1:3) = laps_radar_ext_a(i) (1:3)
                                     read (ext, '(1x,i2)', err=10) num
                                     nradars = nradars + 1
                                     if (nradars .gt. maxradars) then
                                        print *, 'stop get radar bufr data -- increase parameter maxradars.'
                                        return
                                     end if
                                     path_lvl2(nradars) = path_to_radar_a(i)
10                                   continue
                                  end do
                                  if (nradars .le. 0) return

                                  call get_directory('log', dir, len)
                                  if (dir(len:len) .eq. '/') then
                                     path_output(1:len - 1) = dir(1:len - 1)
                                     path_output(len:len) = ' '
                                  else
                                     path_output(1:len) = dir(1:len)
                                  end if

                                  istatus = 1
                                  return

900                               print *, 'error opening file ', filename
                                  istatus = 0
                                  return

901                               print *, 'error reading remap_in ', filename
                                  write (*, remap_nl)
                                  istatus = 0
                                  return

                               end

                               subroutine read_radar_polar(radarfile, i_missing &
                                                           , radarname, slat, slon, salt, istarttime &
                                                           , inumber, elevation, nobeam, nogatez, nogatev &
                                                           , rangez, rangev, spacez, spacev, snyquist &
                                                           , azim, elev, itim, iz, iv, iw, istatus)
!***********************************************************************
! subroutine/function : read_radar_polar
!
! usage :
!      call read_radar_polar(radarfile,i_missing
!    1                ,radarname,slat,slon,salt,istarttime
!    2                ,inumber,elevation,nobeam,nogatez,nogatev
!    3                ,rangez,rangev,spacez,spacev,snyquist
!    4                ,azim,elev,itim,iz,iv,iw,istatus)
!
! description      : to read radar polor-coordinate netcdf-formatted
!                    data for laps's radar importting system.
!
! parameters:
!    ngates: the maximum gate number. (=944)
!    nbeams: the maximum radial ray number. (=420)
!
! arguments :
!  i/o/w   name,      type,       description
!    i     radarfile  c*150       input radar file name.
!                                 ='yyjjjhhmm_elev01','yyjjjhhmm_elev02',...
!                                 = or 'yyjjjhhmm_elevnn'.
!                                   yy: year, jjj: julian day, hh: hour, mm: minute,
!                                   nn: the elevation number.
!    i     i_missing  integer     integer missing value.
!    o     radarname  c*5         radar name. (='rcwf ','rckt ',... )
!    o     slat       real        latitude of radar site. (degrees)
!    o     slon       real        longitude of radar site. (degrees)
!    o     salt       real        height of radar site. (meters)
!    o     istarttime integer*8      starting time, seconds since 1970-1-1 00:00:00.00
!    o     inumber    integer     number of constant ppi-slope.
!    o     elevation  real        elevation angle of ppi slope. (degrees)
!    o     nobeam     integer     total radial ray number of sweep scan.
!    o     nogatez    integer     gate number in one ray for reflectivity (z) data.
!    o     nogatev    integer     gate number in one ray for wind (v) data.
!                                 (for radial velocity or spectral width)
!    o     rangez     real        range to first gate for z field. (km)
!    o     rangev     real        range to first gate for v field. (km)
!    o     spacez     real        space between gates for z field. (km)
!    o     spacev     real        space between gates for v field. (km)
!    o     snyquist   real        nyquist velocity. (m/s)
!    o   azim(nbeams) real array  azimuthal angle of a ray in one sweep. (degrees)
!    o   elev(nbeams) real array  elevation angle of a ray in one sweep. (degrees)
!    o   itim(nbeams) int*8 array time, seconds since 1970-1-1 00:00:00.00
!    o iz(ngates*nbeams) int array the reflectivity data. (100*dbz)
!    o iv(ngates*nbeams) int array the radial velocity data. (100*m/s)
!    o iw(ngates*nbeams) int array the spectral width data. (100*m/s)
!    o     istatus    integer     work or no work message,
!                                 =0, the subroutine doesnot work.
!                                 =1, the subroutine does work.
!
! modules called :
!   netcdf fortran api.
!
! date :
!   original     -- apr. 05, 2006 (shiow-ming, deng)
!***********************************************************************
                                  parameter(ngates=944, nbeams=420)
                                  character*150 radarfile
                                  integer i_missing

                                  character radarname*5
                                  real slat, slon, salt
                                  real*8 starttime
                                  integer*8 istarttime
                                  integer inumber, nobeam, nogatez, nogatev
                                  real elevation, rangez, rangev, spacez, spacev, snyquist
                                  real azim(nbeams), elev(nbeams)
                                  integer*8 itim(nbeams)
                                  integer istatus
                                  integer iz(ngates*nbeams), iv(ngates*nbeams), iw(ngates*nbeams)

                                  integer*2 elevationnumber, numradials, numgatesz, numgatesv
                                  real*8 radialtime(nbeams)
                                  integer*1 zbyte(ngates*nbeams), vbyte(ngates*nbeams)
                                  integer*1 wbyte(ngates*nbeams)

                                  include 'netcdf.inc'

                                  istatus = 0

!-----------------------------------------------------------------------
!c  to open netcdf file.

                                  ier = nf_open(radarfile, nf_nowrite, ncid)
                                  if (ier .ne. 0) return

!-----------------------------------------------------------------------
!c  to get dimensions: z_bin, v_bin, and unlimit.

                                  ier = nf_inq_dimid(ncid, 'z_bin', id_z_bin)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_dimlen(ncid, id_z_bin, len_z_bin)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_dimid(ncid, 'v_bin', id_v_bin)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_dimlen(ncid, id_v_bin, len_v_bin)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_unlimdim(ncid, id_unlimit)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_dimlen(ncid, id_unlimit, len_unlimit)
                                  if (ier .ne. 0) return

!-----------------------------------------------------------------------
!c  to get radar name, latitude, longitude, and altitude.

                                  ier = nf_inq_varid(ncid, 'radarname', id_radarname)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_text(ncid, id_radarname, radarname)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'sitelat', id_sitelat)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_sitelat, slat)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'sitelon', id_sitelon)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_sitelon, slon)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'sitealt', id_sitealt)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_sitealt, salt)
                                  if (ier .ne. 0) return

!-----------------------------------------------------------------------
!c  to get elevationnumber, elevationangle, numradials, esstarttime,
!c         firstgaterangez, firstgaterangev, gatesizez, gatesizev,
!c         numgatesz, numgatesv, and nyquist

                                  ier = nf_inq_varid(ncid, 'elevationnumber', id_elevationnumber)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_int2(ncid, id_elevationnumber, elevationnumber)
                                  if (ier .ne. 0) return
                                  inumber = elevationnumber
                                  ier = nf_inq_varid(ncid, 'elevationangle', id_elevationangle)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_elevationangle, elevation)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'numradials', id_numradials)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_int2(ncid, id_numradials, numradials)
                                  if (ier .ne. 0) return
                                  nobeam = numradials
                                  ier = nf_inq_varid(ncid, 'esstarttime', id_esstarttime)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_double(ncid, id_esstarttime, starttime)
                                  if (ier .ne. 0) return
                                  istarttime = starttime
                                  ier = nf_inq_varid(ncid, 'firstgaterangez', id_firstgaterangez)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_firstgaterangez, rangez)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'firstgaterangev', id_firstgaterangev)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_firstgaterangev, rangev)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'gatesizez', id_gatesizez)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_gatesizez, spacez)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'gatesizev', id_gatesizev)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_gatesizev, spacev)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'numgatesz', id_numgatesz)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_int2(ncid, id_numgatesz, numgatesz)
                                  if (ier .ne. 0) return
                                  nogatez = numgatesz
                                  ier = nf_inq_varid(ncid, 'numgatesv', id_numgatesv)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_int2(ncid, id_numgatesv, numgatesv)
                                  if (ier .ne. 0) return
                                  nogatev = numgatesv
                                  ier = nf_inq_varid(ncid, 'nyquist', id_nyquist)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_nyquist, snyquist)
                                  if (ier .ne. 0) return

!-----------------------------------------------------------------------
!c  to get check (len_z_bin,numgatesz), (len_v_bin,numgatesv), and
!c               (len_unlimit,numradials)

                                  if (len_z_bin .ne. nogatez) then
!          print*,'z_bin: ',len_z_bin,' numgatesz: ',numgatesz
                                     return
                                  end if
                                  if (len_v_bin .ne. nogatev) then
!          print*,'v_bin: ',len_v_bin,' numgatesv: ',numgatesv
                                     return
                                  end if
                                  if (len_unlimit .ne. nobeam) then
!          print*,'unlimit: ',len_unlimit,' numradials: ',numradials
                                     return
                                  end if

!-----------------------------------------------------------------------
!c  to get radialazim, radialelev, and radialtime.

                                  ier = nf_inq_varid(ncid, 'radialazim', id_radialazim)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_radialazim, azim)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'radialelev', id_radialelev)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_real(ncid, id_radialelev, elev)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'radialtime', id_radialtime)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_double(ncid, id_radialtime, radialtime)
                                  if (ier .ne. 0) return
                                  do i = 1, nobeam
                                     itim(i) = radialtime(i)
                                  end do

!-----------------------------------------------------------------------
!c  to get z, v, and w.

                                  ier = nf_inq_varid(ncid, 'z', id_zbyte)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_int1(ncid, id_zbyte, zbyte)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'v', id_vbyte)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_int1(ncid, id_vbyte, vbyte)
                                  if (ier .ne. 0) return
                                  ier = nf_inq_varid(ncid, 'w', id_wbyte)
                                  if (ier .ne. 0) return
                                  ier = nf_get_var_int1(ncid, id_wbyte, wbyte)
                                  if (ier .ne. 0) return
                                  do i = 1, nobeam*nogatez
                                     if ((zbyte(i) .eq. -1) .or. (zbyte(i) .eq. 0)) then
                                        iz(i) = i_missing
                                     else
                                        ii = zbyte(i)
                                        if (ii .lt. 0) ii = ii + 256
                                        iz(i) = 50*(ii - 2) - 3200
                                     end if
                                  end do
                                  do i = 1, nobeam*nogatev
                                     if ((vbyte(i) .eq. -1) .or. (vbyte(i) .eq. 0)) then
                                        iv(i) = i_missing
                                     else
                                        ii = vbyte(i)
                                        if (ii .lt. 0) ii = ii + 256
                                        iv(i) = 50*(ii - 129)
                                     end if
                                     if ((wbyte(i) .eq. -1) .or. (wbyte(i) .eq. 0)) then
                                        iw(i) = i_missing
                                     else
                                        ii = wbyte(i)
                                        if (ii .lt. 0) ii = ii + 256
                                        iw(i) = 50*(ii - 129)
                                     end if
                                  end do

!-----------------------------------------------------------------------
!c  to close netcdf file.

                                  ier = nf_close(ncid)
                                  if (ier .ne. 0) return
                                  istatus = 1

                                  return
                               end
