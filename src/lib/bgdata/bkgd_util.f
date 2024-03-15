      subroutine rotate_lga_winds(ldir,bgmodel,cmodel,fullname
     1,gproj,nx,ny,nz,lon,uw3d,vw3d,uw2d,vw2d,lgb_only)

      use mem_namelist
      implicit none

      character fullname*200
      character cmodel*132
      character gproj*2
      integer nx,ny,nz
      integer i,j,k
      integer bgmodel
      integer istatus,istatus_rot,ishow_timer
      logical ldir,lgb_only
      real    u_true,v_true
      real    u_grid,v_grid
      real    u_true_2d(nx,ny),v_true_2d(nx,ny)
      real    u_grid_2d(nx,ny),v_grid_2d(nx,ny)
      real    lon(nx,ny)
      real    uw3d(nx,ny,nz)
      real    vw3d(nx,ny,nz)
      real    uw2d(nx,ny)
      real    vw2d(nx,ny)
      real    angle(nx,ny,2)  !both grid to true n and true to grid n.
      real    latitude(nx,ny),projrot_latlon(nx,ny)

c
c reset or restore common projection parameters
c
      call reset_lapsparms_common(ldir,bgmodel,cmodel,fullname
     1,gproj,istatus)
c
c build look-up-tables for rotation angles in 2d grid and
c apply this to grid winds
c
!     istatus_rot=ishow_timer()

      latitude = -999. ! since lat is not yet passed in
      call projrot_latlon_2d(latitude,lon,nx,ny,projrot_latlon,istatus)

      do j = 1, ny
      do i = 1, nx
         angle(i,j,1)= -projrot_latlon(i,j)                        !grid 2 true
         angle(i,j,2)= -angle(i,j,1)                               !true 2 grid
      enddo
      enddo

      write(6,*)'rotate_lga_winds - built lookup table'
!     istatus_rot=ishow_timer()

      if(ldir)then   !from grid to true north
c 3d
        if(.not. lgb_only)then
          do k = 1, nz
          do j = 1, ny
          do i = 1, nx

            call     rotate_vec(uw3d(i,j,k),
     1                          vw3d(i,j,k),
     1                          u_true,
     1                          v_true,
     1                          angle(i,j,1))
            uw3d(i,j,k) = u_true
            vw3d(i,j,k) = v_true

          enddo
          enddo
          enddo
        endif
c 2d
        call   rotate_vec_2d(uw2d,
     1                       vw2d,
     1                       u_true_2d,
     1                       v_true_2d,
     1                       angle(1,1,1),nx,ny)

        uw2d = u_true_2d
        vw2d = v_true_2d

      else !rotate from true to grid n.
c 3d
        if(.not. lgb_only)then
          do k = 1, nz
          do j = 1, ny
          do i = 1, nx

            call     rotate_vec(uw3d(i,j,k),
     1                          vw3d(i,j,k),
     1                          u_grid,
     1                          v_grid,
     1                          angle(i,j,2))
            uw3d(i,j,k) = u_grid
            vw3d(i,j,k) = v_grid

          enddo
          enddo
          enddo
        endif
c 2d
        call rotate_vec_2d(uw2d,
     1                     vw2d,
     1                     u_grid_2d,
     1                     v_grid_2d,
     1                     angle(1,1,2),nx,ny)

        uw2d = u_grid_2d
        vw2d = v_grid_2d

      endif

      write(6,*)'end of rotate_lga_winds, lgb_only = ',lgb_only
!     istatus_rot=ishow_timer()

      return
      end
c
c=======================================================================
c
      subroutine rotate_background_uv(nx,ny,nz,lon,bgmodel,cmodel
     +,fullname,gproj,slon0,slat1,slat2,uw,vw,uw_sfc,vw_sfc,lgb_only
     +,istatus)
c
c
c

      use mem_namelist
      implicit none

      integer nx,ny,nz
      integer bgmodel
      integer istatus

      real    slon0,slat1,slat2
      real    std_lon,std_lat1,std_lat2

      real    uw(nx,ny,nz)
      real    vw(nx,ny,nz)
      real    uw_sfc(nx,ny)
      real    vw_sfc(nx,ny)
      real    lon(nx,ny)

      character  gproj*2
      character  fullname*200
      character  cmodel*132

c     character  c6_maproj*6

c     call get_c6_maproj(c6_maproj,istatus)
c     call get_standard_longitude(std_lon,istatus)
c     call get_standard_latitudes(std_lat1,std_lat2,istatus)

      logical lgb_only

      print*,'rotate u/v components'
      print*

c for the case when rotating the bkgd grid winds, we use subroutine
c reset_lapsparms_common to reset the appropriate projection parameters
c within lapsparms.cmn.

      std_lon =standard_longitude
      std_lat1=standard_latitude
      std_lat2=standard_latitude2

c ----------------------------------------------------------------
      if(c6_maproj.eq.'merctr')then
c ----------------------------------------------------------------

         if(gproj.eq.'mc'.or.gproj.eq.'ll')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
            print*,'no rotation required for background'

         else

c rotate from grid north to laps true north (this subroutine uses
c std_lon by virtue of library routines.

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

            print*,'rotate grid-north (bkgd) to true-north (laps)'

            call rotate_lga_winds(.true.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

         endif
c ----------------------------------------------------------------
      elseif(c6_maproj.eq.'plrstr')then
c ----------------------------------------------------------------

         if(gproj.eq.'mc'.or.gproj.eq.'ll'.or.gproj.eq.'le')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

            print*,'rotate true-north (bkgd) to grid-north (laps)'

            call rotate_lga_winds(.false.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

         elseif(gproj.eq.'ps')then

             if(slon0.eq.std_lon)then

c also, this only good if both domains share a common pole point.
c currently we are only considering polar stereo, not local stereo.

                call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
                print*,'no rotation required for background'

             else

                call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

                print*,'rotate grid-north (bkgd) to true-north'
                call rotate_lga_winds(.true.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

                print*,'rotate true-north to grid-north (laps)'
                call rotate_lga_winds(.false.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

             endif

         else

c background is lambert

             call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

             print*,'rotate grid-north (bkgd) to true-north'

             call rotate_lga_winds(.true.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

             print*,'rotate true-north to grid-north (laps)'

             call rotate_lga_winds(.false.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

         endif

c ----------------------------------------------------------------
      else     !laps is lambert --- if(c6_maproj.eq.'lambrt')then
c ----------------------------------------------------------------

         if(gproj.eq.'mc' .or.
     +      gproj.eq.'ll' .or.
     +      gproj.eq.'le')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

            print*,'rotate true-north (bkgd) to grid-north (laps)'   ! because mc/ll/le grids are true north

            call rotate_lga_winds(.false.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

         elseif(gproj.eq.'ps'.or.gproj.eq.'lc')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

            print*,'rotate grid-north (bkgd) to true-north'

            call rotate_lga_winds(.true.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

             print*,'rotate true-north to grid-north (laps)'

             call rotate_lga_winds(.false.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

         endif

      endif

      return
      end
c
c =======================================================================
c
      subroutine print_rotproj(gproj,c6_maproj,lon0,lat1,lat2
     +,std_lon,std_lat1,std_lat2)

      character*2 gproj
      character*6 c6_maproj
      real lon0,lat1,lat2
      real std_lon,std_lat1,std_lat2

      print*,'background: /gproj/lon0/lat1/lat2: '
      print*,'            ',gproj,' ',lon0,lat1,lat2

      print*,'analysis: /c6_maproj/std_lon/std_lat1/std_lat2: '
      print*,'          ',c6_maproj,' ',std_lon,std_lat1,std_lat2

      return
      end
c
c===============================================================================
c
      subroutine thvpc2tq(thv,pc,p,t,q)
c
c *** subprogram:  thvpc2tq - calculates temperature (k) and specific 
c                             humidity (kg/kg) given virtual potential 
c                             temperature (k), pressure (mb), and
c                             condensation pressure (mb).
c
c *** program history log:
c        93-12-20  s. benjamin - original version 
c        96-09-17  j. snook    - es calculated in a table
c
c *** usage:  call thvpc2tq(thv,pc,p,t,q)
c
c *** input argument list:
c        thv    - real  virtual potential temperature (k)
c        pc     - real  condensation pressure (mb)
c        p      - real  pressure (mb)
c
c *** output argument list:
c        t      - real  temperature (k)
c        q      - real  specific humidity (kg/kg)
c
c *** subprograms called:
c        tv2tq  - calculate temp and spec. hum. from virtual
c                    temp and relative humidity
c        es   - calculate saturation vapor pressure (from a table)
c_______________________________________________________________________________
c
      
c
      real tv,rh,p,t,q,thv,kappa,templcl,x,x1,pc
      integer it
      data kappa/0.285714/
c
      include 'bgdata.inc'
c_______________________________________________________________________________
c
      tv=thv*(p*0.001)**kappa
      templcl=thv*(pc*0.001)**kappa
      it=tv*100
      it=min(45000,max(15000,it))
      x =es(it)
      it=templcl*100
      it=min(45000,max(15000,it))
      x1=es(it)
      rh=x1/x * (p-x) / (pc-x1)
      call tv2tq(tv,rh,p,t,q)
      return
      end
c
c===============================================================================
c
      subroutine tv2tq(tv,rh,p,t,q)
c
c *** subprogram:  tv2tq - calculates temperature (k) and specific 
c                          humidity (kg/kg) given virtual temperature (k),
c                          pressure (mb), and relative humidity.
c
c *** program history log:
c        93-01-12  s. benjamin - original version
c
c *** usage:  call tv2tq(tv,rh,p,t,q)
c
c *** input argument list:
c        tv     - real  virtual temperature (k)
c        rh     - real  relative humidity (range 0.0-1.0)
c        p      - real  pressure (mb)
c
c *** output argument list:
c        t      - real  temperature (k) 
c        q      - real  specific humidity (kg/kg)
c
c *** reamrks:
c        it uses an iterative newton-raphson technique.  four iterations are
c        generally adequate to provide convergence to 5 decimal places.
c        the wobus function for saturation vapor pressure over liquid water
c        is used.
c_______________________________________________________________________________
c
      
c
      real tv,rh,p,t,q,t1,estv1,etv,t2,estv2,dt,dum
c
      integer j,it
c
      include 'bgdata.inc'
      
c_______________________________________________________________________________
c
      t1 = tv
c
c *** estv = saturation vapor pressure (mb) for tv.
c
      it=t1*100
      it=min(45000,max(15000,it))
      estv1=es(it)
c
      do j=1,3
c
c ****** etv = vapor pressure (mb) for tv and rh*
c
         etv=estv1*rh
c
c ****** q = mixing ratio for tv (kg/kg).
c
         q=0.62197*etv/(p-etv)
         t2=tv/(1.+0.608*q)
         if (abs(t2-t1) .lt. 0.001) goto 77
c
c ****** estv2 = saturation vapor pressure (mb) for estimated t (=t2).
c
         it=t2*100
         it=min(45000,max(15000,it))
         estv2=es(it)
c
c ****** etv = vapor pressure (mb) for estimated t and rh.
c
         etv=estv2*rh
c
c ****** q = mixing ratio for estimated t and rh (kg/kg).
c
         q=0.62197*etv/(p-etv)
         t=t2*(1.+0.608*q)
c
c ****** recalc. tv.
c
         dt=tv-t
         dum=(estv2-estv1)/(t2-t1)
         etv=estv2+dum*dt*rh
c
c ****** reset t1 and estv1 before next iteration.
c
         t1=t2
         estv1=estv2
c
      enddo
77    continue
      t=t2
c
      return
      end
c
c ---------------------------------
c
      subroutine reset_lapsparms_common(ldir,bgmodel,cmodel
     1,fullname,gproj,istatus)
c
c acquires background projection parameters and uses them to
c temporarily replace lapsparms parameters. it is then possible
c to use the library projection routines to calculate wind rotation
c angles for the background.
c
      use mem_namelist
      implicit none

      character*200 fullname
      character*132 cmodel
      character*2   gproj
      character*1   cgrddef
      integer       istatus
      integer       bgmodel
      integer       nxbg,nybg
      integer       nzbg_ht
      integer       nzbg_tp
      integer       nzbg_sh
      integer       nzbg_uv
      integer       nzbg_ww
      logical       ldir
      real          dlat,dlon
      real          sw(2),ne(2)
      real          dxbg,dybg
      real          la1,lo1,la2,lo2
      real          xcen,ycen
      real          xsw,ysw,xne,yne
      real          erad
      real          rlatcen,rloncen

      character*6   c6_maproj_save
      save          c6_maproj_save
      real          centrallat_save
      save          centrallat_save
      real          centrallon_save
      save          centrallon_save
      real          lat0_save
      save          lat0_save
      real          lat1_save
      save          lat1_save
      real          lon0_save
      save          lon0_save

      integer       isave
      data          isave/0/
      save          isave



      if(isave.eq.0)then
         if(ldir)then
            c6_maproj_save=c6_maproj
            centrallat_save=grid_cen_lat
            centrallon_save=grid_cen_lon
            lat0_save=standard_latitude
            lat1_save=standard_latitude2
            lon0_save=standard_longitude
            isave=1
         else
            print*,'return to rotate_background_uv: isave=0; ldir=false'
            return 
         endif
      elseif(isave.eq.1)then        
c restore original nest7grid.parms settings
         c6_maproj=c6_maproj_save
         grid_cen_lat=centrallat_save
         grid_cen_lon=centrallon_save
         standard_latitude=lat0_save
         standard_latitude2=lat1_save
         standard_longitude=lon0_save
         isave=0
         return
      endif

      print*
      print*,'reset_lapsparms_common: isave = ',isave
      print*
      print*,'bgmodel: ',bgmodel
      print*,'cmodel: ', trim(cmodel)
      print*,'fullname:',trim(fullname)
      print*

      call get_bkgd_mdl_info(bgmodel,cmodel,fullname
     &,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     &,gproj,dlat,dlon,grid_cen_lat,grid_cen_lon
     &,dxbg,dybg,standard_latitude,standard_latitude2
     &,standard_longitude,sw,ne,cgrddef,istatus)

c temporarily set c6_maproj.
      if(gproj.eq.'lc')c6_maproj='lambrt'
      if(gproj.eq.'ps')c6_maproj='plrstr'
      if(gproj.eq.'mc')c6_maproj='merctr' 
cc
c if ruc_native
      if(bgmodel.eq.5 .or. bgmodel.eq.3 .or. bgmodel.eq.13)then
         if(trim(cmodel).eq.'cwb_20fa_lambert_nf'.or.
     +      trim(cmodel).eq.'cwb_20fa_lambert_re'.or.
     +      trim(cmodel).eq.'ruc40_native'.or.
     +      trim(cmodel).eq.'ruc')then  
            la1=sw(1)
            lo1=sw(2)
            la2=ne(1)
            lo2=ne(2)
            print*,'*** sw lat/lon = ',la1,lo1,'***'
            print*,'*** ne lat/lon = ',la2,lo2,'***'
            call get_earth_radius(erad,istatus)
c get x/y for grid center
            call latlon_to_xy(la1,lo1,erad,xsw,ysw) 
            call latlon_to_xy(la2,lo2,erad,xne,yne) 
            xcen=(xne+xsw)/2.
            ycen=(yne+ysw)/2.
            call xy_to_latlon(xcen,ycen,erad,rlatcen,rloncen)
            print*,'*** center lat/lon= ',rlatcen,rloncen,' ***'
         endif
      endif

      return
      end
