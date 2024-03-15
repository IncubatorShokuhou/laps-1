      subroutine dprep_ruc2_pub(nx,ny,nz,ht,pr,sh,uw,vw,th,gproj)
      implicit none
      include 'bgdata.inc'
      integer nx,ny,nz,i,j,k
      real ht(nx,ny,nz),pr(nx,ny,nz),sh(nx,ny,nz),uw(nx,ny,nz)
     +     ,vw(nx,ny,nz),th(nx,ny,nz)
      real cp,g,r,cpog,kappa
      parameter (cp=1004.686,g=9.80665,r=287.053,cpog=cp/g,kappa=r/cp)
      real tv, psi(nx,ny),psj(nx,ny),lat(nx,ny),lon(nx,ny),
     .       angle(nx,ny)
      character*(*) gproj
      real missing
      
      
c
c *** common block variables for lambert-conformal grid.
c
      integer nx_lc,ny_lc,nz_lc  !no. of lc domain grid points
      real*4 lat1,lat2,lon0,       !lambert-conformal std lat1, lat, lon
     .       sw(2),ne(2)           !sw lat, lon, ne lat, lon
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0,sw,ne

c *** convert pascals to mb.
c *** compute tv from thetav.
c *** compute height from msf.
c *** compute tp (returned in th) from tv.
c *** compute sh from mr.
c
      do k=1,nz
         do j=1,ny
            do i=1,nx
               missing = max(pr(i,j,k),th(i,j,k),sh(i,j,k))
               if(missing.lt.missingflag) then
                  pr(i,j,k)=pr(i,j,k)*0.01
                  tv=th(i,j,k)*(pr(i,j,k)*0.001)**kappa
                  th(i,j,k)=tv/(1.+0.61*sh(i,j,k))
                  sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))
               endif
            enddo
         enddo
      enddo
c
c *** fill lambert-conformal common block variables.
c
      gproj='lc'
      nx_lc=nx
      ny_lc=ny
      lat1=25.0
      lat2=25.0
      lon0=-95.0

      sw(1)=16.2810
      sw(2)=-126.1378
      ne(1)=55.4818
      ne(2)=-57.3794
c **** no longer needed *****
c *** convert ruc winds from grid north to true north.
c
c      do j=1,ny
c         do i=1,nx
c            psi(i,j)=float(i)
c            psj(i,j)=float(j)
c         enddo
c      enddo
c      call psij_2_latlon(nx*ny,psi,psj,lat,lon)
c
c      call uvgrid_to_uvtrue_a(uw,vw,lon,lon0,nx,ny,nz,angle)
c
      return
      end
