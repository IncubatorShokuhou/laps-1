
      subroutine sun_eclipse_parms(i4time,rlat_r4,rlon_r4,ht,idebug
     1                            ,alt_r4,azi_r4,dist_m_r4
     1                            ,earth_radius,elgms_r4
     1                            ,r4_mag,r4_obsc,obsc_limbc)

      include 'trigd.inc'
      use mem_allsky, only: nc

      implicit real*8(a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)
      include '../../include/astparms.for'
      include 'wa.inc'

      angdif(x,y)=dmod(x-y+9.4247779607694d0,6.2831853071796d0)
     .-3.1415926535897932d0

      atan3(x,y)=dmod((datan2(x,y)+6.2831853071796d0),6.2831853071796d0)

      character blank,twi,moon,rise,set,signm,signs
      character*3 mnth(12)
      real*8 m,mag,lon,lhmsh,
     1  mxg,myg,mzg,mxr,myr,mzr,lat
      real*8 nx,ny,nz,ex,ey,ez,lst
      integer rah,decd,frame,elgmc,altdk1
      character*5 ctt,ctr,cts,btt,bmt,c5_blank
      character c8_appm*8,c8_apps*8
      data mode/1/,ifl/0/,time/0.0d0/,timeol/0.d0/,iprint/1/
      data frame/2/
      data rise/'r'/,set/'s'/,blank/' '/,c5_blank/'     '/
      data mnth/'jan','feb','mar','apr','may','jun','jul','aug','sep'
     .,'oct','nov','dec'/

      integer i4time,istatus
      integer i4time_last /0/
      real rlat_r4,rlon_r4,alt_r4,azi_r4,dist_m_r4
      real alm_r4,azm_r4,elgms_r4,r4_mag,ht,earth_radius

      real r4_ratio,r4_obsc,obsc_limbc(nc),solar_eclipse_magnitude

      save sxg,syg,szg,mxg,myg,mzg,i4time_last,et,ut,lst
      save nx,ny,nz,ex,ey,ez,znx,zny,znz,tx,ty,tz,tmag
      save topo_flag

      n_loc = 1

      if(idebug .ge. 1)then
        write(6,*)
        write(6,*)' subroutine sun_eclipse_parms:',i4time,i4time_last
      endif

      rlon = rlon_r4

      if(i4time .ne. i4time_last)then
        write(6,*)
        write(6,*)' initialize sun_eclipse_parms:',i4time,i4time_last

        i4time_last = i4time
 
        call i4time_to_jd(i4time,tb,istatus)
        tf = tb

        ti=ti/1440.d0 ! input in minutes
c
c enter loops
        l = 1
        rsn = 1.0
        iter=3

c
c enter time loops
        t = tb

        deltat=delta_t(t)
        ut = t
        et = ut + deltat

!                         day  deg rad
        call sidereal_time(ut,rlon,lst)

        do 1000 ii = 1,iter
c
c calculate position of earth (1950 coordinates - antedated)
840       continue
!         call posint(et-rsn/c,1,sxg,syg,szg)
          call posin(et-rsn/c,1,sxg,syg,szg)
          call xyz_to_polar_r(-sxg,-syg,-szg,decs,ras,rsn)
c         write(13,843)r,si
843       format(6f10.6)

1000    continue
c
c calculate coordinates of sun (coordinates of date)
        call preces(t1950,et,sxg,syg,szg,1)
c
c calculate position of moon (topocentric coordinates of date)
        call moon_brwn(et,mxg,myg,mzg)

        write(6,*)
        write(6,*)' sun coords  = ',sxg,syg,szg
        write(6,*)' moon coords = ',mxg,myg,mzg

!       calculate eclipse conditions at each grid point
        phi = rlat_r4*rpd
!       topo_flag = ht/earth_radius   ! based on mean radius
        topo_flag = ht/(r_e_km*1000.) ! based on equatorial radius
c
c calculate orientation of observer's zenith (coordinates of date)
        call topo_ff1(phi,rlon,ut,txz,tyz,tzz)
        tzmag = sqrt(txz**2 + tyz**2 + tzz**2) ! based on equ radius
c
c calculate geocentric topo position including flattening factor
        call topo(phi,rlon,ut,tx,ty,tz)
        tmag = sqrt(tx**2 + ty**2 + tz**2)

!       observer location relative to earth center
!       note that 'ht' is measured along surface normal
        tx = tx + topo_flag * txz
        ty = ty + topo_flag * tyz
        tz = tz + topo_flag * tzz

!       get direction cosines based on alt/azi (north along the horizon)
!       ra is 180 degrees from ra of meridian
!       dec is 90 - latitude
        ran = 180d0 + lst/rpd
        decn = 90d0 - rlat_r4       
        nx = cosd(decn) * cosd(ran)
        ny = cosd(decn) * sind(ran)
        nz = sind(decn)

!       direction cosines of zenith
        znx = txz / tzmag
        zny = tyz / tzmag
        znz = tzz / tzmag

!       east horizon is n horizon cross product with zenith unit vector
        call crossproduct(nx,ny,nz,znx,zny,znz,ex,ey,ez)

        write(6,*)' lst,rlon-lst ',lst/rpd,rlon-lst/rpd ! deg
        write(6,*)' decn,ran ',decn,ran                 ! deg
        write(6,*)' zenith unit vector    ',znx,zny,znz
        write(6,*)' n horizon unit vector ',nx,ny,nz
        write(6,*)' e horizon unit vector ',ex,ey,ez

      endif ! new time

      sdist_m = dist_m_r4 
!     sdist_au = (sdist_m / earth_radius) * tmag
      sdist_au = (sdist_m / 1000d0) / km_per_au 

      sinalt = sind(alt_r4)
      cosalt = cosd(alt_r4)
      sinazi = sind(azi_r4)
      cosazi = cosd(azi_r4)

!     component of ray along northern horizon
      rx = tx + sdist_au * cosazi * cosalt * nx
      ry = ty + sdist_au * cosazi * cosalt * ny
      rz = tz + sdist_au * cosazi * cosalt * nz

      if(idebug .ge. 2)then
         write(6,*)' topo (observer) coords = ',tx,ty,tz
         write(6,*)' ray coords (n) =         ',rx,ry,rz
      endif

!     component of ray along eastern horizon
      rx = rx + sdist_au * sinazi * cosalt * ex
      ry = ry + sdist_au * sinazi * cosalt * ey
      rz = rz + sdist_au * sinazi * cosalt * ez

      if(idebug .ge. 2)then
         write(6,*)' ray coords (e) = ',rx,ry,rz
      endif

!     component of ray along zenith
      rx = rx + sdist_au * sinalt * znx
      ry = ry + sdist_au * sinalt * zny
      rz = rz + sdist_au * sinalt * znz

      if(idebug .ge. 1)then

!        ray location distance from earth center (au)
         raymag = sqrt(rx**2 + ry**2 + rz**2)

         raylat = asind(rz/raymag) ! geocentric
         raylst = atan2d(ry,rx)
         raylon = raylst + (rlon-lst/rpd) ! deg
         if(raylon .gt. +180d0)raylon = raylon - 360d0
         if(raylon .lt. -180d0)raylon = raylon + 360d0
         raymsl = ((raymag-tmag)*km_per_au) * 1000d0

         write(6,*)' rlat,rlon,ut,ht,topo_flag = '
     1              ,rlat_r4,rlon,ut,ht,topo_flag
         write(6,*)' sdist_m,sdist_au = ',sdist_m,sdist_au
         write(6,*)' ray  coords = ',rx,ry,rz
         write(6,*)' ray  lat/lst/lon = ',raylat,raylst,raylon
         write(6,*)' ray  mag/tmag/msl = ',raymag,tmag,raymsl
      endif

      sxr = sxg - rx
      syr = syg - ry
      szr = szg - rz

      call xyz_to_polar_r(-sxr,-syr,-szr,decs,ras,rsn)
      has=angdif(ramr,ras)
!     write(13,*)szg,rsn,decs/rpd,ras/rpd,has/rpd

c
c calculate alt and az of sun
      call equ_to_altaz_r(decs,has,phi,als,azs)
      als = als/rpd
      call refract(als,apps,pres)

      if(als .lt. -1.0)then
          c8_apps = '        '
      else
          write(c8_apps,1002)apps
1002      format(f8.2)
      endif

      azs = azs/rpd

      mxr=mxg-rx
      myr=myg-ry
      mzr=mzg-rz
      call xyz_to_polar_r(mxr,myr,mzr,decm,ram,rmn)

      ham=angdif(ramr,ram)

      call anglevectors(-mxr,-myr,-mzr,sxr,syr,szr,elgarg)
      elgmst = elgarg/rpd

      call anglevectors(-mxg,-myg,-mzg,sxg,syg,szg,elgarg)
      elgmsg = elgarg/rpd

      if(idebug .ge. 2)then
         write(6,*)' decs / ras / has = '
     1              ,decs/rpd,ras/rpd,has/rpd
         write(6,*)' decm / ram / ham / elsmst = '
     1              ,decm/rpd,ram/rpd,ham/rpd,elgmst
      endif

!     calculate apparent positions
      call equ_to_altaz_r(decm,ham,phi,alm,azm)

      alm = alm/rpd
      call refract(alm,appm,pres)

      if(alm .lt. -1.0)then
            c8_appm = '       '
      else
          write(c8_appm,1001)appm
1001      format(f8.2)
      endif

      azm = azm/rpd

      call cjymd(ut,iyear,month,date)

      idt = int(date)
      frac = (date - int(date)) * 2.d0 * pi
      call clock(frac,ctt)

      call coords(decm,signm,dec_dg_m,dec_mn_m,dec_sc_m,
     1            ram,       ra_hr_m, ra_mn_m, ra_sc_m)

      call coords(decs,signs,dec_dg_s,dec_mn_s,dec_sc_s,
     1            ras,       ra_hr_s, ra_mn_s, ra_sc_s)


      alm_r4 = alm
      azm_r4 = azm
!     elgms_r4 = 10.
      elgms_r4 = elgmst

      phase_angle_deg = 180. - elgmst ! approximate

!     call phase_func_moon(phase_angle_deg,mode,area_rel,area_es
!    1                    ,phase_corr_moon)

      if(.false.)then
!         calculate solar eclipse magnitude
!         call magnitude(0,0,sxt,syt,szt,0.,0.,0.,amag
!    1                                          ,diam_sun)
!         call magnitude(1,1,sxt,syt,szt,
!    1      sxg+mxg,syg+myg,szg+mzg,amag,diam_moon)

      endif ! .true. 

      diam_sun = 1800.
      diam_moon = 1852.2
      elgsec = elgmst * 3600.

      overlap_sec = -(elgsec - 0.5 * (diam_sun + diam_moon))
      if(overlap_sec .gt. 0.)then
          solar_eclipse_magnitude = overlap_sec / diam_sun
!         larger term will increase exponent
!         fracmag = solar_eclipse_magnitude**30.
!         exponent = 1.35*(1.-fracmag) + 0.5*fracmag
!         r4_obsc = solar_eclipse_magnitude**exponent
!         r4_obsc = min(r4_obsc,1.0)
          r4_ratio = diam_moon/diam_sun
          call get_obscuration(solar_eclipse_magnitude,r4_ratio
     1                        ,r4_obsc,obsc_limbc)
          if(r4_obsc/r4_obsc .ne. 1.0)then
              write(6,*)' error in sun_eclipse_parms r4_obsc ='
     1                     ,r4_obsc,solar_eclipse_magnitude,r4_ratio
              stop
          endif
      else
          solar_eclipse_magnitude = 0.
          r4_obsc = 0.
          obsc_limbc(:) = 0.
      endif

      r4_mag = solar_eclipse_magnitude
  
      if(idebug .ge. 2)then
          write(6,*)' diam_sun,diam_moon,overlap_sec,r4_mag ',
     1                diam_sun,diam_moon,overlap_sec,r4_mag
          write(6,*)' return from sun_eclipse_parms'
      endif
c
      return
      end

      subroutine sun_moon(i4time,lat_2d,lon_2d,ni,nj,is,js,alm_r4,azm_r4
     1                   ,idebug                           ! i
     1                   ,ht,earth_radius                  ! i
     1                   ,elgms_r4,r4_mag,r4_rmn           ! o
     1                   ,geo_dec,geo_ra,geo_sublon,geo_dist ! o
     1                   ,solar_eclipse_magnitude,r4_obsc,obsc_limb) ! o   

      include 'trigd.inc'
      use mem_allsky, only: nc
      implicit real*8(a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)
!     include '../util/utilparms.for'
      include '../../include/astparms.for'
      include 'wa.inc'

      angdif(x,y)=dmod(x-y+9.4247779607694d0,6.2831853071796d0)
     .-3.1415926535897932d0

      atan3(x,y)=dmod((datan2(x,y)+6.2831853071796d0),6.2831853071796d0)

      character blank,twi,moon,rise,set,signm,signs
      character*3 mnth(12)
      real*8 i,m,mag,lon,lhmsh,
     1  mxg,myg,mzg,mxt,myt,mzt,lat,lst
      integer rah,decd,frame,elgmc,altdk1
      character*5 ctt,ctr,cts,btt,bmt,c5_blank
      character c8_appm*8,c8_apps*8
      character*20 c20_site
      dimension lat(9),lon(9)
      data mode/1/,ifl/0/,time/0.0d0/,timeol/0.d0/,iprint/0/
      data frame/2/
      data rise/'r'/,set/'s'/,blank/' '/,c5_blank/'     '/
      data mnth/'jan','feb','mar','apr','may','jun','jul','aug','sep'
     .,'oct','nov','dec'/

      integer i4time,istatus
      real lat_2d(ni,nj)
      real lon_2d(ni,nj)
      real alm_r4,azm_r4,elgms_r4,r4_mag,r4_rmn
      real solar_eclipse_magnitude,ht,earth_radius
      real r4_ratio,r4_obsc,obsc_limb,obsc_limbc(nc)
      real geo_dec,geo_ra,geo_sublon,geo_dist
c
      n_loc = 1
      lat = lat_2d(is,js)
      lon = lon_2d(is,js)

      call i4time_to_jd(i4time,tb,istatus)
      tf = tb
      c20_site = 'las'

!     read parameters
!     open(11,file='sun_moon.parms',status='old')
!     read(11,*)n_loc                           ! coordinates
!     do idum = 1,n_loc
!       read(11,101)c20_site
!101     format(a20)
!       read(11,*)lat(idum),lon(idum),pres
!     enddo
!     read(11,*)tb                              ! time interval
!     read(11,*)tf
!     read(11,*)ti
!     read(11,*)topo_flag

      ti=ti/1440.d0 ! input in minutes
c
c enter loops
      l = 1
      phi=lat(l)*rpd
      rsn = 1.0
      sinphi=dsin(phi)
      cosphi=dcos(phi)
      iter=3

c
c enter time loops
      t = tb
c
c write headings
      if(iprint.eq.1.and.ifl.eq.0)write(6,11)c20_site,lat(l),lon(l)
11    format(1h1,30x,a20,'  lat=',f7.2,'  lon=',f7.2/)
10    if(iprint.eq.1.and.ifl.eq.0)write(6,1)
1     format(
     1 ' year mon dt   ut                   moon         ',
     1 '     coordinates of date        sun',
     1           26x,'        elong  altdif'/
     1            21x,'alt    app     az        dec          ra',
     1      10x,'alt    app     az        dec          ra',7x,'(mag)')


      deltat=delta_t(t)
      ut = t
      et = ut + deltat
c
c calculate orientation of observer's zenith (coordinates of date)
      call topo_ff1(phi,lon(l),ut,txz,tyz,tzz)
      ramr=atan3(tyz,txz)

      if(iprint.eq.1)write(6,*)' lat,lon = ',lat(l),lon(l)
      if(iprint.eq.1)write(6,*)' txz,tyz,tzz,ramr',txz,tyz,tzz,ramr/rpd

      topo_flag = 1.0 + ht/earth_radius

c
c calculate geocentric topo position including flattening factor
      call topo(phi,lon(l),ut,tx,ty,tz)

      tx = tx * topo_flag
      ty = ty * topo_flag
      tz = tz * topo_flag

      do 1000 ii = 1,iter
c
c calculate position of earth (1950 coordinates - antedated)
840   continue
!     call posint(et-rsn/c,1,sxg,syg,szg)
      call posin(et-rsn/c,1,sxg,syg,szg)
      sxt = sxg - tx
      syt = syg - ty
      szt = szg - tz
c     write(13,843)r,si
843   format(6f10.6)

1000  continue

c
c calculate coordinates of sun (coordinates of date)
      call preces(t1950,et,sxg,syg,szg,1)
      call preces(t1950,et,sxt,syt,szt,1)

      call xyz_to_polar_r(-sxt,-syt,-szt,decs,ras,rsn)
      has=angdif(ramr,ras)
!     write(13,*)szg,rsn,decs/rpd,ras/rpd,has/rpd

c
c calculate alt and az of sun
      call equ_to_altaz_r(decs,has,phi,als,azs)
      als = als/rpd
      call refract(als,apps,pres)

      if(als .lt. -1.0)then
          c8_apps = '        '
      else
          if(iprint.eq.1)write(c8_apps,1002)apps
1002      format(f8.2)
      endif


      azs = azs/rpd

c
c calculate position of moon (topocentric coordinates of date)
      call moon_brwn(et,mxg,myg,mzg)
      mxt=mxg-tx
      myt=myg-ty
      mzt=mzg-tz
      call xyz_to_polar_r(mxt,myt,mzt,decm,ram,rmn)
      r4_rmn = rmn

!     geocentric coordinates of moon      
      call xyz_to_polar_r(mxt,myt,mzt,decmg,ramg,rmng)
      geo_dec = decmg / rpd
      geo_ra = ramg / rpd
      geo_dist = rmng * km_per_au * 1000. ! meters
      rlon = 0d0
!                       day  deg rad
      call sidereal_time(ut,rlon,lst)
      geo_sublon = modulo((ramg-lst)/rpd,360d0)
      if(idebug .ge. 2)then
         write(6,*)' ramg          (deg)',ramg/rpd
         write(6,*)' greenwich lst (deg)',lst/rpd
         write(6,*)' geo_sublon    (deg)',geo_sublon
      endif

!     insert star in lunar position
!     decm = ?
!     ram = ?
!     call preces(2433282.423357d0,et,decm,ram,rdum,2)

      ham=angdif(ramr,ram)

      call anglevectors(-mxt,-myt,-mzt,sxt,syt,szt,elgarg)
      elgmst = elgarg/rpd

      call anglevectors(-mxg,-myg,-mzg,sxg,syg,szg,elgarg)
      elgmsg = elgarg/rpd

      if(iprint.eq.1)write(6,*)' decm / ram / ham = '
     1                          ,decm/rpd,ram/rpd,ham/rpd

!     calculate apparent positions
      call equ_to_altaz_r(decm,ham,phi,alm,azm)

      alm = alm/rpd
      call refract(alm,appm,pres)

      if(alm .lt. -1.0)then
          c8_appm = '       '
      else
          if(iprint.eq.1)write(c8_appm,1001)appm
1001      format(f8.2)
      endif

      azm = azm/rpd

      call cjymd(ut,iyear,month,date)

      idt = int(date)
      frac = (date - int(date)) * 2.d0 * pi
      call clock(frac,ctt)

      call coords(decm,signm,dec_dg_m,dec_mn_m,dec_sc_m,
     1          ram,       ra_hr_m, ra_mn_m, ra_sc_m)

      call coords(decs,signs,dec_dg_s,dec_mn_s,dec_sc_s,
     1          ras,       ra_hr_s, ra_mn_s, ra_sc_s)


      alm_r4 = alm
      azm_r4 = azm
      elgms_r4 = elgmst

      phase_angle_deg = 180. - elgmst ! approximate

      call phase_func_moon(phase_angle_deg,mode,area_rel,area_es
     1                    ,phase_corr_moon)

      r4_mag = -12.74 + phase_corr_moon

      if(.true.)then

!         calculate solar eclipse magnitude
!         call magnitude(0,0,sxt,syt,szt,0.,0.,0.,amag
!    1                                          ,diam_sun)
!         call magnitude(1,1,sxt,syt,szt,
!    1      sxg+mxg,syg+myg,szg+mzg,amag,diam_moon)

          diam_sun = 1800.
          diam_moon = 1852.2
          elgsec = elgmst * 3600.

          overlap_sec = -(elgsec - 0.5 * (diam_sun + diam_moon))
          solar_eclipse_magnitude = overlap_sec / diam_sun

          if(solar_eclipse_magnitude .gt. 0.)then
!             larger term will increase exponent
              fracmag = solar_eclipse_magnitude**30.
              exponent = 1.35*(1.-fracmag) + 0.5*fracmag
!             r4_obsc = solar_eclipse_magnitude**exponent
!             r4_obsc = min(r4_obsc,1.0)
              r4_ratio = diam_moon/diam_sun
              call get_obscuration(solar_eclipse_magnitude,r4_ratio
     1                            ,r4_obsc,obsc_limbc)
          else
              r4_obsc = 0.
              obsc_limbc(:) = 0.
          endif

          obsc_limb = obsc_limbc(2)

          if(idebug .ge. 2)then
             write(6,*)
             write(6,*)' sun_moon debug:'
             write(6,*)' i4time,tb,lat,lon',i4time,tb,lat(1),lon(1)
             write(6,*)' deltat/sec',deltat,deltat*86400.
             write(6,*)' elgmst/sec',elgmst,elgmst*3600.
             write(6,*)' diam_sun,diam_moon,overlap_sec,mag,r4_obsc ',
     1diam_sun,diam_moon,overlap_sec,solar_eclipse_magnitude,r4_obsc
             write(6,*)' decs / ras / has = '
     1                  ,decs/rpd,ras/rpd,has/rpd
             write(6,*)' decm / ram / ham / elsmst = '
     1                  ,decm/rpd,ram/rpd,ham/rpd,elgmst
             write(6,*)' et / mxt / myt / mzt = '
     1                  ,et,mxt,myt,mzt
             write(6,*)
          endif

      endif ! .true.

!     correct moon's magnitude for lunar eclipse
      if(elgmsg .gt. 177.)then

!         calculate lunar eclipse magnitude
          delta_moon_au = mag(mxg,myg,mzg)
          delta_sun_au  = mag(sxg,syg,szg)
          delta_moon_km = delta_moon_au * km_per_au
          delta_sun_km  = delta_sun_au  * km_per_au

          diam_moon_km = r_m_km * 2d0

          elgsupl = 180d0 - elgmsg

          x_moon =     delta_moon_km * cosd(elgsupl)
          y_moon = abs(delta_moon_km * sind(elgsupl))

          umbral_length_km = delta_sun_km / (r_s_km / r_e_km - 1d0)

          angrad = asin(r_e_km/umbral_length_km)
          umbral_height_km = (umbral_length_km - x_moon)
     1                              * tan(angrad) * 1.02

          dist_center_km = (y_moon - umbral_height_km) * cos(angrad)
          umbral_magnitude = -(dist_center_km / diam_moon_km) + 0.5

          penumbral_length_km = delta_sun_km / (r_s_km / r_e_km + 1d0)

          angrad = asin(r_e_km/penumbral_length_km)
          penumbral_height_km  = (penumbral_length_km + x_moon)
     1                              * tan(angrad)

          dist_center_km = (y_moon - penumbral_height_km) * cos(angrad)
          penumbral_magnitude = -(dist_center_km / diam_moon_km) + 0.5

          penumbral_width = penumbral_magnitude - umbral_magnitude

          bri_2nd = .0002

          if(umbral_magnitude .gt. 1.)then          ! total    
              if(iprint.eq.1)write(6,*)' total'           
              if(iprint.eq.1)write(6,*)' umbral_magnitude = '
     1                                  ,umbral_magnitude
              frac_bri = bri_2nd / umbral_magnitude**4
          elseif(umbral_magnitude .gt. 0.)then      ! partial
              if(iprint.eq.1)then
                write(6,*)' partial'         
                write(6,*)' umbral_magnitude = ',umbral_magnitude
                write(6,*)' penumbral_magnitude = ',penumbral_magnitude
                write(6,*)' penumbral_width = ',penumbral_width
              endif
              frac_bri = ((1.0 - umbral_magnitude) * 0.5) 
     1                 / penumbral_width
!             frac_bri = max(frac_bri,bri_2nd) 
              frac_bri = frac_bri + umbral_magnitude * bri_2nd
          elseif(penumbral_magnitude .gt. 0.)then   ! penumbral
              if(iprint.eq.1)then
                write(6,*)' penumbral'       
                write(6,*)' penumbral_magnitude = ',penumbral_magnitude
                write(6,*)' penumbral_width = ',penumbral_width
              endif
              frac_bri = 1.0 - (penumbral_magnitude * 0.5)
     1                       / penumbral_width
          else
              frac_bri = 1.0
          endif

          rmag_corr = -log10(frac_bri) * 2.5

          r4_mag = r4_mag + rmag_corr

          if(iprint.eq.1)write(6,451)frac_bri,rmag_corr,r4_mag
451       format(' lunar eclipse: frac_bri/rmag_corr/r4_mag'
     1          ,f11.6,2f9.3)

      endif ! lunar eclipse magnitude correction

c
c write out data
500   continue
      if(iprint.eq.1)write(6,2)iyear,mnth(month),idt,ctt,
     1  alm,c8_appm(2:8),azm,signm,int(abs(dec_dg_m)),int(dec_mn_m),dec_
     1sc_m,
     1                        int(ra_hr_m),      int(ra_mn_m), ra_sc_m,
     1  als,c8_apps(2:8),azs,signs,int(abs(dec_dg_s)),int(dec_mn_s),dec_
     1sc_s,
     1                        int(ra_hr_s),      int(ra_mn_s), ra_sc_s,
     1                          elgmst,alm-als
2     format(1x,i4,1x,a3,1x,i2,1x,a5,
     1  f7.2,a7,f8.2,2x,a1,i2,1x,i2,f5.1,2x,i2,1x,i2,f6.2,
     1          2x,
     1  f7.2,a7,f8.2,2x,a1,i2,1x,i2,f5.1,2x,i2,1x,i2,f6.2,
     1          1x,f7.2,f7.2)

c
9999  continue

      return
      end


        subroutine coords(dec_rd,sign,dec_dg,dec_mn,dec_sc,
     1                  ra_rd,       ra_hr, ra_mn, ra_sc)

        implicit real*8(a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)
        include '../../include/astparms.for'

        character*1 sign,plus/'+'/,minus/'-'/

        if(dec_rd .ge. 0)then
            sign = plus
        else
            sign = minus
        endif

        dec_dg = abs(dec_rd / rpd)
        dec_mn = (dec_dg - int(dec_dg)) * 60d0
        dec_sc = (dec_mn - int(dec_mn)) * 60d0

        ra_hr = ra_rd / rpd / 15d0
        ra_mn = (ra_hr - int(ra_hr)) * 60d0
        ra_sc = (ra_mn - int(ra_mn)) * 60d0

        return
        end

        subroutine get_obscuration(rmag,ratio,obscuration,obsc_limbc)

!       eclipse obscuration where ratio is moon/sun angular diameter

!       limb darkening reference (hestroffer and magnan, 1997)
!       www.physics.hmc.edu/faculty/esin/a101/limbdarkening.pdf

!       optional switch to first pierce method in astrophysical quantities 
        use mem_allsky, only: nc
        include 'wa.inc'

        parameter (method = 3) ! first pierce method

        real obsc_limbc(nc)
        real a_a(nc)  / 0.78, 0.74, 0.54/
        real b_a(nc)  / 0.39, 0.43, 0.60/
        real c_a(nc)  /-0.57,-0.56,-0.44/

!       mean brightness divided by central brightness
        real cb(nc)   /.855 ,.850, .820/ ! larger increases obsc

        parameter (pi = 3.14159265)
        seg_area(theta,r) = 0.5 * r**2 * (theta - sin(theta))

        real mu
        mu(r) = sqrt(1.-r**2)
        rlimbdk1(r,u,alphal) = 1. - (u * (1. - (mu(r))**alphal))
        rlimbdk3(r,a,b,c) = a + b*mu(r) 
     1                        + c * (1.0-mu(r) * log(1.0 + 1./mu(r)))

        u = 1.00

        if(rmag .ge. 0.99999)then
            obscuration = 1.0
            obsc_limbc(:) = 1.0         
            return
        elseif(rmag .eq. 0.)then 
            obscuration = 0.0
            obsc_limbc(:) = 0.0            
        endif

        rs = 1.0
        rm = rs*ratio

        oa = -(rmag * (2.*rs) - (rs + rm))
        a = (rs**2 - rm**2 + oa**2) / (2. * oa)

        alpha = 2. * acos(a/rs)
        beta  = 2. * acos((oa-a)/rm)

        seg1 = seg_area(alpha,rs)
        seg2 = seg_area(beta,rm)  

        obscuration = (seg1 + seg2) / pi

        if(rmag .le. 0.5)then
            obsc_limbc(:) = obscuration
        else
!           calculate mean radius of "quadratic" segment
            rmean = 1. - ((1.-rmag)*(2./3.)) 
            do ic = 1,nc
!               alphal = 0.8 
                alphal = -0.023 + .292 / wa(ic)
                if(method .eq. 1)then
                    rmean_ill = rlimbdk1(rmean,u,alphal)
                else
                    rmean_ill = rlimbdk3(rmean,a_a(ic),b_a(ic),c_a(ic))
                    rmean_ill = rmean_ill / cb(ic)
                endif
                obsc_limbc(ic) = 1.0 - ((1.0 - obscuration) * rmean_ill)
            enddo ! ic
        endif

!       write(6,*)'rmag/obscuration/obsc_limbc = '
!    1            ,rmag,obscuration,obsc_limbc
!       write(6,*)'rmean/rmean_ill/mu = ',rmean,rmean_ill,mu(rmean)
!       write(6,*)'mu(r)/mu(r)**alphal',mu(rmean),mu(rmean)**alphal
!       write(6,*)'rmean_ill derived',(1.-obsc_limbc)/(1.-obscuration)
!       stop

!       write(6,1)rmag,ratio,oa,a
!1      format('rmag/ratio/oa/a',4f10.5)
!       write(6,2)alpha,beta,seg1,seg2,obscuration
!2      format('alpha/beta/seg1/seg2/obsc',5f10.5)
!       stop

        return
        end
