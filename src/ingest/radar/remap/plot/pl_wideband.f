
      program pl_wideband
c
c * read and plot netcdf wsr88d level 2 data
c   19 nov 2001  dongsoo kim  original 
c
      integer    rad_max, ang_max, nz_max
      parameter (rad_max = 999, ang_max = 999, nz_max = 16)
      character  yydddhhmm*9,filename*150,atype*7(16),xxxx*6(12)
      character  path*150
      integer    dim_dim, rad_dim(16), z_dim, v_dim, vcp
c        rad_dim  !number of radials, slightly larger than 360
c        z_dim    !number reflectivity bins in each radial
c        v_dim    !number wind bins in each radial
c        vcp      !volume coverage pattern
      real  refl(460,rad_max,16)
     &       ,dwind(920,rad_max,16)
     &       ,specw(920,rad_max,16)
      real  a_ang(rad_max,16), e_ang(rad_max,16)
      integer   nsites
     &       ,vcpmode31(16),vcpmode11(16),vcpmode21(16)
      real    site_lat, site_lon, site_alt
      real    radar_lat(99), radar_lon(99), radar_alt(99)
      character radarname*5(99)
     &         ,sitename*132(99)

c * working
      character chmm*2, ch_ang*5, ch_alt*5, stname*4
      integer*2 image_z(460, rad_max)
     &       ,image_v(920, rad_max)
     &       ,image_w(920, rad_max)
      real  azim_ang(rad_max), elev_ang(rad_max), resolv
      integer   colia(460,rad_max),ibin,vlvl,nlvl,elelvl
      real      xloc(460,460), yloc(460,460)
      integer   colia_v(920,rad_max)
      real      uloc(920,920), vloc(920,920)
      real      deg2rad

c * output:
      data atype/'_elev01','_elev02','_elev03','_elev04','_elev05'
     &          ,'_elev06','_elev07','_elev08','_elev09','_elev10'
     &          ,'_elev11','_elev12','_elev13','_elev14','_elev15'
     &          ,'_elev16'/
      data  xxxx/'/kama/','/kcys/','/kddc/','/kftg/','/kfws/'
     &          ,'/kgld/','/kict/','/kinx/','/klbb/','/klzk/'
     &          ,'/ksgf/','/ktlx/'/
      data vcpmode31/1,3,5,6,8,0,0,0,0,0,0,0,0,0,0,0/
      data vcpmode11/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/
      data vcpmode21/1,3,5,6,7,8,9,10,11,0,0,0,0,0,0,0/
c          vcp = 11 !precipitation/severe weather mode
c * graphics for reflectivity
      integer   lnd2(16)
      character*3  llb2(17)
      data lnd2 /3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/
      data llb2 /' 0',' 0',' 5','10','15','20','25','30','35','40',
     &           '45','50','55','60','65','70','75'/
      character title*25, subtit1*25, subtit2*25
c * graphics for doppler winds
      integer   lnd3(16)
      character*3  llb3(17)
      data lnd3 /3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/
      data llb3 /' 0',' 0',' 5','10','15','20','25','30','35','40',
     &           '45','50','55','60','65','70','75'/
c
      deg2rad = 3.1415926/180.d0
      print *,' choose radar station (a4)'
      write(6,'(12a6)') (xxxx(k),k=1,12)
      read(5,'(a4)') stname
      
      print *,' enter yydddhhmm in /public/data/radar/wsr88d/... '
      read(5,'(a9)') yydddhhmm
c
      ncount = 1
      do nsite = 1,1
         ncount = ncount + 1
         do 10 itry=1,1 ! 10,19
!           encode(2,'(i2)',chmm) itry
            write(chmm,1)itry
1           format(i2)

!           yydddhhmm(9:9) = chmm(2:2)
!           xxxx(nsite)(2:5) = stname

            if(stname .ne. 'rcwf')then
                filename = '/public/data/radar/wsr88d/wideband/'
     1                 //stname//'/netcdf/'//yydddhhmm//atype(1)
            else
                filename = 
     1          '/data/lapb/import/lapsdat/tw/radar/rcwf/lvl2/'
     1                 //yydddhhmm//atype(1)
            endif

            call get_wideband_hdr(filename,dim_dim,z_dim,v_dim
     &                   ,site_lat,site_lon,site_alt,vcp,istat)
            if(istat .eq. -1)  then
               goto 10
            elseif(istat .ne. -1) then
               goto 20
            endif
  10     continue
  20     continue
c
        radar_lat(ncount) = site_lat
        radar_lon(ncount) = site_lon
        radar_alt(ncount) = site_alt
         write(6,*) filename
         print *,' site name=', stname
         print *,' site lat=', radar_lat(ncount)
         print *,' site lon=', radar_lon(ncount)
         print *,' site alt=', radar_alt(ncount)
         print *,' scan mode =', vcp

        vlvl = 0
        do 30 n16=1,nz_max
         if(vcp .eq. 31 .or. vcp .eq. 32) nelev=vcpmode31(n16)
         if(vcp .eq. 11) nelev=vcpmode11(n16)
         if(vcp .eq. 21) nelev=vcpmode21(n16)
         if(nelev.eq. 0) goto 40

!        filename =
c     &     'wideband'//xxxx(nsite)//'netcdf/'//yydddhhmm//atype(nelev)
!    &        'wideband/'//yydddhhmm//atype(1)

         if(stname .ne. 'rcwf')then
             filename = '/public/data/radar/wsr88d/wideband/'
     1                  //stname//'/netcdf/'//yydddhhmm//atype(nelev)
         else
             filename = 
     1          '/data/lapb/import/lapsdat/tw/radar/rcwf/lvl2/'
     1                 //yydddhhmm//atype(nelev)
         endif

         call get_wideband_hdr(filename,dim_dim,z_dim,v_dim
     &                ,site_lat,site_lon,site_alt,vcp,istat)
         call get_wideband_netcdf(filename,dim_dim,z_dim,v_dim
     &         ,image_z,image_v,image_w,azim_ang,elev_ang,resolv
     &         ,istat)
            if(istat .eq. -1) then
cdk                print *,' no file exist in this elev ', atype(nelev)
                goto 30
            endif
            vlvl = vlvl + 1
            rad_dim(vlvl) = dim_dim
            do j=1,rad_dim(vlvl)
               a_ang(j,vlvl) = azim_ang(j)
               e_ang(j,vlvl) = elev_ang(j)
               do i=1,z_dim
                  refl(i,j,vlvl) = ((image_z(i,j)-2.)/2.) - 32.
               enddo
               do i=1,v_dim
                  dwind(i,j,vlvl) = (image_v(i,j)-129)*resolv
                  specw(i,j,vlvl) = image_w(i,j)
               enddo
            enddo
cdk        write(6,*) 'azim_ang at elev_ang ', elev_ang(1) 
cdk        write(6,'(10f8.2)') (azim_ang(irad),irad=1,rad_dim(vlvl))
cdk        write(6,*) 'image_z at radial=100 '
        write(6,'(10f6.1)') (dwind(100,j,vlvl),j=1,rad_dim(vlvl))
   30   continue
   40   continue
        nlvl = vlvl   !number of vertical levels
c * 
c * convert (r,theta) |--> (x0,y0) per elevation angle
c * 
      print *,' number of vertical levels is ', nlvl
cdk      elelvl = 1 
cdk      print *,' enter vertical level 1 -',nlvl
cdk      read(5,'(i2)') nlvl

        write(6,*)' enter range of levels to plot '
        read(5,*)lvl_start,lvl_end

c
        do elelvl=lvl_start,lvl_end
c
        print *,' plotting elevation level ',e_ang(1,elelvl),' deg'      
!       encode(5,'(f5.2)',ch_ang) e_ang(1,elelvl)
        write(ch_ang,2)e_ang(1,elelvl)
 2      format(f5.2)

!       encode(5,'(i5)',ch_alt) int(radar_alt(ncount))
        write(ch_alt,3)int(radar_alt(ncount))
 3      format(i5)


c * first reflectivity
      write(6,*)' reflectivity'
      do j=1,rad_dim(elelvl)
      do i=1,z_dim
        xloc(i,j) = i*sin(a_ang(j,elelvl)*deg2rad)/1200.
        yloc(i,j) = i*cos(a_ang(j,elelvl)*deg2rad)/1200.
          if(refl(i,j,elelvl) .gt. 0. .and. 
     &       refl(i,j,elelvl).lt. 255. ) then
            colia(i,j) = ir_color(refl(i,j,elelvl))
          else
            colia(i,j) = 1
          endif
      enddo
      enddo

c * second dopple winds
      write(6,*)' velocity'
      do j=1,rad_dim(elelvl)
      do i=1,v_dim
        uloc(i,j) = i*sin(a_ang(j,elelvl)*deg2rad)/2400.
        vloc(i,j) = i*cos(a_ang(j,elelvl)*deg2rad)/2400.
          if(dwind(i,j,elelvl) .gt. -10. .and. 
     &       dwind(i,j,elelvl).lt. 30. ) then
            colia_v(i,j) = dw_color(dwind(i,j,elelvl))
          else
            colia_v(i,j) = 1
          endif
      enddo
      enddo
c
      call opngks
      call radar_color
      title = 'wsr88d wideband at '//stname
      subtit1 = 'elev angle '//ch_ang
      subtit2 = 'site elev(m) '//ch_alt
  
      call set(0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,1)
      do j=1,365,1
      do i=10,z_dim,2
         xo = xloc(i,j) + 0.5
         yo = yloc(i,j) + 0.5
         scale = 0.004*sqrt(2.35*i/rad_dim(elelvl))
         call ngdots(xo,yo,1,scale,colia(i,j))
      enddo
      enddo
      call plchlq(0.05,0.95,title,0.02,0.,-1.)
      call plchlq(0.95,0.95,yydddhhmm,0.02,0.,1.)
      call plchlq(0.05,0.92,subtit1,0.015,0.,-1.)
      call plchlq(0.05,0.89,subtit2,0.015,0.,-1.)
      call lblbar(0,0.05,0.95,0.01,0.06,16,1.,0.4,lnd2,0,llb2,17,1)
      call frame
c
      if(.false.)then

      call set(0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,1)
      do j=1,365,1
      do i=10,v_dim,4
         xo = uloc(i,j) + 0.5
         yo = vloc(i,j) + 0.5
         scale = 0.003*sqrt(2.35*i/rad_dim(elelvl))
         call ngdots(xo,yo,1,scale,colia_v(i,j))
      enddo
      enddo
      call plchlq(0.05,0.95,title,0.02,0.,-1.)
      call plchlq(0.95,0.95,yydddhhmm,0.02,0.,1.)
      call plchlq(0.05,0.92,subtit1,0.015,0.,-1.)
      call plchlq(0.05,0.89,subtit2,0.015,0.,-1.)
      call lblbar(0,0.05,0.95,0.01,0.06,16,1.,0.4,lnd2,0,llb2,17,1)
      call frame

      endif ! plot velocity
c
        enddo !elelvl
c
      enddo ! nsite
      call clsgks

      stop
      end
c
      subroutine radar_color
c
c  the color scheme is customized for nexrad image convention
c * see ofcm fm handbook no.11 part c a-7
c
      call gscr (1,0,1.,1.,1.) !white
c
c    foreground colors
c
      call gscr (1,1,0.,0.,0.) !black

      call gscr  (1,2, 0.00,0.00,0.00) !black
      call gscr  (1,3, 0.61,0.61,0.61) !med gray
      call gscr  (1,4, 0.46,0.46,0.46) !dk gray
      call gscr  (1,5, 1.00,0.67,0.67) !lt pink
      call gscr  (1,6, 0.93,0.55,0.55) !med pink
      call gscr  (1,7, 0.79,0.44,0.44) !dk pink
      call gscr  (1,8, 0.00,0.98,0.56) !lt green
      call gscr  (1,9, 0.00,0.98,0.00) !med green
      call gscr  (1,10,1.00,1.00,0.44) !lt yellow
      call gscr  (1,11,0.82,0.82,0.38) !dark yellow
      call gscr  (1,12,1.00,0.38,0.38) !lt red
      call gscr  (1,13,0.85,0.00,0.00) !med red
      call gscr  (1,14,0.70,0.00,0.00) !dk red
      call gscr  (1,15,0.00,0.00,1.00) !blue
      call gscr  (1,16,1.00,1.00,1.00) !white
      call gscr  (1,17,0.91,0.00,1.00) !purple
      call gscr  (1,18,0.91,0.00,1.00) !purple
      return
      end
c
      integer function ir_color(echo)
      real echo

      ir_color = int(echo/5)+2

      return
      end
c
      integer function dw_color(echo)
      real echo

      dw_color = int(echo/2)+2

      return
      end
