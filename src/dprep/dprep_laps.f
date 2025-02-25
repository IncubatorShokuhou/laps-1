      subroutine dprep_laps(i4time,nx,ny,nz,
     .                      pr,ht,tp,uw,vw,mr,
     .			    sht,spr,slp,sth,smr,suw,svw,istatus)
c
      implicit none
c
      integer i4time
      integer*4	nx,ny,nz
      real*4 cp,kappa
      parameter (cp=1004.,kappa=287./1004.)
c
      real*4 ht(nx,ny,nz),     !laps 3d height (m)
     .       tp(nx,ny,nz),     !laps 3d temperature (k)
     .       uw(nx,ny,nz),     !laps 3d u-wind (m/s)
     .       vw(nx,ny,nz),     !laps 3d v-wind (m/s)
     .       mr(nx,ny,nz),     !laps 3d mixing ratio (kg/kg)
c     .       ex(nx,ny,nz),     !laps 3d exner function
     .       pr(nx,ny,nz),     !laps pressure levels (mb) in -> exner out
     .       sht(nx,ny),       !laps surface height (m)
     .       spr(nx,ny),       !laps surface pressure (mb)
     .       slp(nx,ny),       !laps mean sea level pressure (mb)
     .       sth(nx,ny),       !laps surface potential temperature (k)
     .       suw(nx,ny),       !laps surface u-wind (m/s)
     .       svw(nx,ny),       !laps surface v-wind (m/s)
     .       smr(nx,ny),       !laps surface mixing ratio (kg/kg)
     .       lat(nx,ny),
     .       lon(nx,ny),
c    .       lat0,lon0,        !laps polar stereo grid pole point
     .       pri(nz),
     .       factor
c    .       prbot,dpr, factor
c
      integer istatus,i,j,k
c
c
c_______________________________________________________________________________
c
c
c *** convert  3d temp to theta.
c *** compute exner function.
c

      do k=1,nz
         do j=1,ny
            do i=1,nx
               factor = (1.0e5/pr(i,j,k))**kappa
               tp(i,j,k)=tp(i,j,k)*factor
               pr(i,j,k)=cp/factor
               mr(i,j,k)=mr(i,j,k)/(1.-mr(i,j,k))
            enddo
         enddo
      enddo


      istatus = 1
      
      return
      end
c
c-------------------------------------------------------------------------------c
      subroutine get_laps_data(laps_data_root,
     .                         i4time,nx,ny,nz,
     .                         pr,ht,tp,uw,vw,mr,
     .			       sht,spr,slp,sth,smr,suw,svw,istatus)
c
      implicit none
c
      integer*4	nx,ny,nz
c
      real*4 ht(nx,ny,nz),
     .       tp(nx,ny,nz),
     .       mr(nx,ny,nz),
     .       uw(nx,ny,nz),
     .       vw(nx,ny,nz),
     .       pr(nx,ny,nz),
     .       sht(nx,ny),
     .       spr(nx,ny),
     .       slp(nx,ny),
     .       sth(nx,ny),
     .       smr(nx,ny),
     .       suw(nx,ny),
     .       svw(nx,ny),
     .       lat(nx,ny),
     .       lon(nx,ny)
c
      integer*4 i4time,error(2),i,j,k,istatus
c
      character*(*) laps_data_root
      character*256 datadir
      character*10 units
      character*125 comment
c     character*3 sfcfields(6)
c
c *** common block variables for lambert-conformal grid.
c
      integer nx_lc,ny_lc,nz_lc  !no. of lc domain grid points
      real*4 lat1,lat2,lon0,       !lambert-conformal std lat1, lat, lon
     .       sw(2),ne(2)           !sw lat, lon, ne lat, lon
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0,sw,ne
c_______________________________________________________________________________
c
      error(1)=1   !good istatus
      error(2)=0   !bad istatus

c
c *** define laps pressures.
c
      call get_pres_3d(i4time,nx,ny,nz,pr,istatus)

c
c *** get laps 3d temperature.
c
      call get_laps_3d(i4time,nx,ny,nz,'lt1','t3',units,comment,tp,
     .     istatus)
      if (istatus .ne. 1) then
         print *,'error getting laps temperature data'
         return
c         stop
      endif
c
c *** get laps 3d heights.
c
      call get_laps_3d(i4time,nx,ny,nz,'lt1','ht',units,comment,ht,
     .     istatus)
      if (istatus .ne. 1) then
         print *,'error getting laps height data'
         return
c         stop
      endif
c
c *** get laps 3d specific humidity.
c
      call get_laps_3d(i4time,nx,ny,nz,'lq3','sh',units,comment,mr,
     .     istatus)

      if (istatus .ne. 1) then
         print *,'error getting laps moisture data.'
         return
c         stop
      endif

c
c *** get laps u winds.
c
      call get_laps_3d(i4time,nx,ny,nz,'lw3','u3',units,comment,uw,
     .     istatus)

      if (istatus .ne. 1) then
         print *,'error getting laps wind data.'
         return
c         stop
      endif
c
c *** get laps v winds.
c
      call get_laps_3d(i4time,nx,ny,nz,'lw3','v3',units,comment,vw,
     .     istatus)

      if (istatus .ne. 1) then
         print *,'error getting laps wind data.'
         return
c         stop
      endif


c
c *** get laps surface data.
c
      call s_len(laps_data_root,i)
      datadir = laps_data_root(1:i)//'/lapsprd/lsx/'
      call get_laps_sfc(datadir,i4time,'lsx',nx,ny,
     .                  spr,slp,sth,smr,suw,svw,istatus)
      if (istatus .ne. 1) then
         print *,'error getting laps surface data.'
         return
c         stop
      endif


c
c *** get laps topo.
c

      call get_laps_domain(nx,ny,'nest7grid',lat,lon,sht,istatus)


      sw(1)=lat(1,1)
      sw(2)=lon(1,1)
      ne(1)=lat(nx,ny)
      ne(2)=lon(nx,ny)

c
c *** do quick qc.
c
      do k=nz/2,1,-1 
      do j=1,ny
      do i=1,nx 
         if (uw(i,j,k) .eq. 1.e-30) uw(i,j,k)=uw(i,j,k+1)
         if (vw(i,j,k) .eq. 1.e-30) vw(i,j,k)=vw(i,j,k+1)
         if (mr(i,j,k) .eq. -1.e+30 .or. mr(i,j,k) .gt. 1.e10) 
     .       mr(i,j,k)=mr(i,j,k+1)
      enddo
      enddo
      enddo
      do k=nz/2+1,nz 
      do j=1,ny
      do i=1,nx 
         if (uw(i,j,k) .eq. 1.e-30) uw(i,j,k)=uw(i,j,k-1)
         if (vw(i,j,k) .eq. 1.e-30) vw(i,j,k)=vw(i,j,k-1)
         if (mr(i,j,k) .eq. -1.e+30 .or. mr(i,j,k) .gt. 1.) 
     .       mr(i,j,k)=mr(i,j,k-1)
      enddo
      enddo
      enddo
      
      


      return
      end
c
c===============================================================================
c
      subroutine get_laps_sfc(dir,i4time,ext,nx,ny,
     .                        spr,slp,sth,smr,suw,svw,istatus)
      implicit none
c
      integer*4 nx,ny,nfld
      parameter (nfld=6)    !number of surface fields to be read.
c
      real*4 spr(nx,ny),    !laps surface pressure (mb)
     .       slp(nx,ny),    !laps mean sea level pressure (mb)
     .       sth(nx,ny),    !laps surface potential temperature (k)
     .       smr(nx,ny),    !laps surface mixing ratio (kg/kg)
     .       suw(nx,ny),    !laps surface u-wind (m/s)
     .       svw(nx,ny),    !laps surface v-wind (m/s)
     .       grid(nx,ny,nfld)
c
      integer level_req(nfld),len,i4time,i,j,istatus
c
      character*3 fldname(nfld)
      character*(*) dir,ext
      character*500 ldir
      character*1 junk1(nfld),junk2(nfld),junk3(nfld)
c
      data fldname/'u','v','msl','th','ps','mr'/
      data level_req/0,0,0,0,0,0/
c     .     13,                  !pressure id
c     .            9,        !mean sea level pressure id
c     .           11,        !theta id
c     .           15,        !mr id
c     .            1,        !u-wind id
c     .            2/        !v-wind id
c_______________________________________________________________________________
c
c *** read requested laps surface data.
c
      ldir=dir
      len=index(ldir,' ')-1

      call read_laps(i4time,i4time,ldir(1:len),ext,nx,ny,nfld,nfld, 
     +     fldname,level_req,junk1,junk2,junk3,grid,istatus)



c      call get_laps_2d(i4time,ext,nx,ny,
c     .                 nfld,fldid,grid,istatus)
      if (istatus .ne. 1) return
c
c *** fill requested fields.
c *** convert surface pressure and mslp from pa to mb.
c *** convert mixing ratio from g/kg to kg/kg.
c
      do j=1,ny
         do i=1,nx
            suw(i,j)=grid(i,j,1)
            svw(i,j)=grid(i,j,2)
            slp(i,j)=grid(i,j,3)*0.01
            sth(i,j)=grid(i,j,4)
            spr(i,j)=grid(i,j,5)*0.01
            smr(i,j)=grid(i,j,6)*0.001
         enddo
      enddo
c
      istatus=1
      return
      end
      subroutine get_laps_corners(nx,ny,sw,ne)
      integer nx,ny
      real sw(2),ne(2)
      real lat(nx,ny), lon(nx,ny), topo(nx,ny)

      call get_laps_domain(nx,ny,'nest7grid',lat,lon,topo,istatus)

      sw(1)=lat(1,1)
      sw(2)=lon(1,1)
      ne(1)=lat(nx,ny)
      ne(2)=lon(nx,ny)
      return
      end
      
