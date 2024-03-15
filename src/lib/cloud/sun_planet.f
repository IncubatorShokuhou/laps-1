
      subroutine sun_planet(i4time,n_planet,rlat,rlon
     1                     ,dec_r4,ra_r4,alm_r4,azm_r4,elgms_r4,r4_mag) 

      implicit real*8(a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)

      include '../../include/astparms.for'
!     include '../planets/posint_inc.for'

      angdif(x,y)=dmod(x-y+9.4247779607694d0,6.2831853071796d0)
     .-3.1415926535897932d0

      atan3(x,y)=dmod((datan2(x,y)+6.2831853071796d0),6.2831853071796d0)

      character blank,twi,moon,rise,set,sign,plus,minus
      dimension r(3),rr(3),rho(3),rri(3),p(3),q(3),w(3),pp(3),qq(3)
     +,ww(3),d1(3),d2(3),i1(3),i2(3)
      character*3 mnth(12)
      real*8 i,m,mag,lnod,mu,nlc,lon,mdegtl,lhmsh,mx,my,mz,lat
      real*8 mx_1950, my_1950, mz_1950
      integer rah,decd,frame,elgmc,altdk1
      character*5 ctt,ctr,cts,btt,bmt,c5_blank
      character c8_appm*8,c8_apps*8
      character*1 c_observe
      character*20 c20_site
      character*8  names(13)
!     character*8  names(13)/'earth','mercury','venus','mars','jupiter','saturn',' ',' ',' ',' ',' ',' ',' '/
      dimension lat(9),lon(9)

      real*8 maglimd_r8,maglimt_r8,maglimn_r8,maglim_r8
      real alm_r4,azm_r4,r4_mag,elgms_r4,rlat,rlon,dec_r4,ra_r4

!     cubert(x)=dexp(dlog(dabs(x))/3.)*x/dabs(x)

      data mode/1/,ifl/0/,time/0.0d0/,timeol/0.d0/,iprint/1/
      data frame/2/
      data rise/'r'/,set/'s'/,blank/' '/,c5_blank/'     '/
      data plus/'+'/,minus/'-'/
      data mnth/'jan','feb','mar','apr','may','jun','jul','aug','sep'
     .,'oct','nov','dec'/
      data names/'earth','mercury','venus','mars','jupiter','saturn',
     1         'uranus','neptune','pluto',4*' '/

!     saturn ring direction cosines
!     data zx,zy,zz/-.0912836831d0, -.0724471759d0, -.9931861335d0/

!     read parameters
!     open(11,file='sun_planet.parms',status='old')
!     read(11,*)n_loc,n_planet                          ! coordinates
!     do idum = 1,n_loc
!       read(11,101)c20_site
101     format(a20)
!       read(11,*)lat(idum),lon(idum),pres
!     enddo
!     read(11,*)tb                              ! time interval
!     read(11,*)tf
!     read(11,*)ti

      lat(1) = rlat
      lon(1) = rlon

      call i4time_to_jd(i4time,tb,istatus)
      tf = tb
      c20_site = 'las'

      ti=ti/1440.d0 ! input in minutes
      rph=pi/12.
      argp=argp*rpd
      lnod=lnod*rpd
      i=i*rpd
      o=obliq1950*rpd
c
c enter loops
      l = 1
      phi=lat(l)*rpd
      rsn = 1.0
      sinphi=dsin(phi)
      cosphi=dcos(phi)
      iter=3

c enter time loops
      t = tb
c
c write headings
      if(iprint.eq.1.and.ifl.eq.0)write(13,11)c20_site,lat(l),lon(l)
11    format(1h1,30x,a20,'  lat=',f7.2,'  lon=',f7.2/)
10    if(iprint.eq.1.and.ifl.eq.0)write(13,1)names(n_planet)
1     format(
     1 '  year mon dt   ut         ',a8,'        ',
     1 '          coordinates of date       sun',
     1            21x,'   elong   altdif altdif2   mag    ill    diam'
     1               ,'   mglm    vis'/
     1 '                       alt     app     az       dec    ra  ',
     1          '        alt     app     az       dec    ra')


400   if((t - tf) * ti .gt. 0.)goto9999

      deltat=delta_t(t)
      ut = t
      et = ut + deltat
c
c calculate position of observer's zenith (coordinates of date)
      call topo_ff1(phi,lon(l),ut,txz,tyz,tzz)
      ramr=atan3(tyz,txz)

      do iter = 1,3 ! light time iteration

c calculate position of earth (coordinates of date)
          tc = et - delta_planet/c
          call posint(tc,1,rri(1),rri(2),rri(3))
!         call posin(tc,1,rri(1),rri(2),rri(3))
          call preces(t1950,et,rri(1),rri(2),rri(3),1)
c         write(13,843)r,rri
843       format(6f10.6)


c calculate position of planet (coordinates of date)
          tc = et - delta_planet/c
          call posint(tc,n_planet,r(1),r(2),r(3))
!         call posin(tc,n_planet,r(1),r(2),r(3))
          call preces(t1950,et,r(1),r(2),r(3),1)
          delta_planet = sqrt((r(1)-rri(1))**2 + (r(2)-rri(2))**2
     1                                         + (r(3)-rri(3))**2 )
c
c calculate coordinates of sun (coordinates of date)
          tc = et - rsn/c
          call posint(tc,1,rri1,rri2,rri3)
!         call posin(tc,1,rri1,rri2,rri3)
          call preces(t1950,et,rri1,rri2,rri3,1)
          call xyz_to_polar_r(-rri1,-rri2,-rri3,decs,ras,rsn)

      enddo ! light time iter

      has=angdif(ramr,ras)
!     write(13,*)rri(3),rsn,decs/rpd,ras/rpd,has/rpd

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

!     relative position is planet vector minus earth vector (coords of date)
      call topo(phi,lon(l),ut,tx,ty,tz)
      mx=r(1)-(rri(1)+tx)
      my=r(2)-(rri(2)+ty)
      mz=r(3)-(rri(3)+tz)

      mx_1950 = mx
      my_1950 = my
      mz_1950 = mz

!     precess planet minus earth vector to coordinates of date
!     call preces(t1950,et,mx,my,mz,1)

      call xyz_to_polar_r(mx,my,mz,decm,ram,rmn)

!     insert star in lunar position
!     decm = ?
!     ram = ?
!     call preces(2433282.423357d0,et,decm,ram,rdum,2)

      ham=angdif(ramr,ram)
      call anglevectors(-mx,-my,-mz,rri(1),rri(2),rri(3),elgarg)
      elgms = elgarg/rpd

      call equ_to_altaz_r(decm,ham,phi,alm,azm)

      alm = alm/rpd
      call refract(alm,appm,pres)

      if(alm .lt. -1.0)then
          c8_appm = '        '
      else
          write(c8_appm,1001)appm
1001      format(f8.2)
      endif

      azm = azm/rpd

      call cjymd(ut,iyear,month,date)

      idt = int(date)
      frac = (date - int(date)) * 2.d0 * pi
      call clock(frac,ctt)

      call magnitude(n_planet,0,rri(1),rri(2),rri(3),r(1),r(2),r(3)
     1                                  ,amag,diam_sec)
      call phase(rri(1),rri(2),rri(3),r(1),r(2),r(3),phase_angle,r_ill)

      if(n_planet .eq. 6)then
!         calculate saturn's ring opening
          call anglevectors(zx,zy,zz,r(1),r(2),r(3),bprime_r)
          call anglevectors(zx,zy,zz,mx_1950,my_1950,mz_1950,b_r)
          bprime_d = bprime_r / rpd - 90.
          b_d      = b_r      / rpd - 90.
!         write(6,*)b_d,bprime_d
          if(b_d * bprime_d .lt. 0.)then ! dark side visible
              product = -b_d * bprime_d
              write(6,3)iyear,mnth(month),idt,ctt,b_d,bprime_d,product
 3            format(1x,i5,1x,a3,1x,i2,1x,a5,f9.3,f9.3,f9.4)

              write(14,3)iyear,mnth(month),idt,ctt,b_d,bprime_d,product

              write(13,2)iyear,mnth(month),idt,ctt,
     1          alm,c8_appm,azm,decm/rpd,ram/rpd,
     1          als,c8_apps,azs,decs/rpd,ras/rpd,
     1                  elgms,alm-als,amag,r_ill

          endif
      endif

!     convert to real
      alm_r4 = alm
      azm_r4 = azm
      r4_mag = amag
      elgms_r4 = elgms
      dec_r4 = decm/rpd
      ra_r4 = ram/rpd

      return

c     call visibility routine
      altdf = alm - als
      call vi(amag,elgms,altdf,al1_r8,al2_r8,0.d0,180.d0
     1       ,maglimd_r8,maglimt_r8,maglimn_r8,maglim_r8
     1       ,vis_r8,c_observe)

c
c write out data
      if(n_planet .ne. 6)then

          altdf2 = (alm-als) - (amag*2.0)

          write(13,2)iyear,mnth(month),idt,ctt,
     1          alm,c8_appm,azm,decm/rpd,ram/rpd,
     1          als,c8_apps,azs,decs/rpd,ras/rpd,
     1                  elgms,alm-als,altdf2,amag,r_ill,
     1                  diam_sec,
     1                  maglim_r8,vis_r8,c_observe
2         format(1x,i5,1x,a3,1x,i2,1x,a5,
     1          f8.2,a8,f8.2,f8.3,f8.3,
     1          4x,
     1          f8.2,a8,f8.2,2x,f6.2,f8.2,
     1          1x,f7.1,2f8.3,f8.2,f7.3,f7.2,f7.1,f7.1,a1)

      endif

!     t = t + ti
!     goto400

c
9999  return
      end
