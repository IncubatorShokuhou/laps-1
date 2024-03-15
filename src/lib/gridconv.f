      subroutine polar_stereographic
c
c *** routines to convert from geographical lat, lon to 
c        polar stereographic grid i, j and vice-versa.
c     equations and code mostly obtained from rams.
c     snook (12/20/95)
c
ccc      implicit none
c
      integer np,n
c
      real glat(np),glon(np),    !earth lat (deg n), lon (deg +e)
     .       pslat,pslon,          !pol ste. lat, lon (deg n, deg +e)
     .       psi(np),psj(np),      !pol ste. i, j
     .       xmin,ymin,dx,dy
c
      common /pscorner/xmin,ymin,dx,dy
c
c===============================================================================
c
      entry latlon_2_psij(np,glat,glon,psi,psj)
c_______________________________________________________________________________
c
      call ps_param(xmin,ymin,dx,dy)
      do n=1,np
         call geoll_2_psll(glat(n),glon(n),pslat,pslon)

         call psll_2_psij(pslat,pslon,psi(n),psj(n))

      enddo
c
      return
c
c===============================================================================
c
      entry psij_2_latlon(np,psi,psj,glat,glon)
c_______________________________________________________________________________
c
c
      do n=1,np
         call psij_2_psll(psi(n),psj(n),pslat,pslon)
         call psll_2_geoll(pslat,pslon,glat(n),glon(n))
      enddo
c
      return
c
      end
c
c===============================================================================
c
      subroutine geoll_2_psll(glat,glon,pla,plo)
c
c     convert geographical lat/lon coordinates to polar stereographic
c     ditto with the pol.ste. pole at rlat,wlon1 (these names are
c     used to be compatible to the ramsin-parameters)
c     longitude:-180 ; 180 positive east (on input)
c              :   0 ; 360 positive east (on output)
c     latitude : -90 ;  90 posive on northern hemisphere
      include 'trigd.inc'
c     the result is rotated 270 degrees relative to 'standard pol.ste.'
c     wlon1 is defined in the same way as the input
c     approach so as to get the x-axis to point towards the east, and the
c     y-axis towards the north along 0 degrees (at np south along 180)
c
c     tsp 20/06-89
      double precision pi180,c1,c2,c3,c4,c5,c6,arg2a,bb,pla1,alpha
     +   ,plo1,pla90,argu2
c
      integer nx,ny,nz           !no. of ps domain grid points
      real rlat,wlon1,rota,       !pol ste. std lat, lon and rotation
     .       sw(2),ne(2)           !sw lat, lon, ne lat, lon
      common /psgrid/nx,ny,nz,rlat,wlon1,rota,sw,ne
c_______________________________________________________________________________
c
c     constants
c
      c1=1.
      pi180 = dasin(c1)/90.
c
c     set flag for n/s hemisphere and convert longitude to <0 ; 360> interval
c
      if(rlat.ge.0.0) then
         hsign= 1.0
      else
         hsign=-1.0
      end if
      glor=glon
      if(glor.lt.0.0) glor=360.0+glor
      rwlon1=wlon1
      if(rwlon1.lt.0.0) rwlon1=360.0+wlon1
c
c     test for a n/s pole case
c
      if(rlat.eq.90.0) then
         pla=glat
         plo=amod(glor+270.0-wlon1,360.0)
         go to 2000
      end if
      if(rlat.eq.-90.0) then
         pla=-glat
         plo=amod(glor+270.0,360.0)
         go to 2000
      end if
c
c     test for longitude on 'greenwich or date line'
c
      if(glor.eq.rwlon1) then
         if(glat.gt.rlat) then
            pla=90.0-glat+rlat
            plo=90.0
         else
            pla=90.0-rlat+glat
            plo=270.0
         end if
         go to 2000
      end if      
      if(amod(glor+180.0,360.0).eq.rwlon1) then
         pla=rlat-90.0+glat
         if(pla.lt.-90.0) then
            pla=-180.0-pla
            plo=270.0
         else
            plo= 90.0
         end if
         go to 2000         
      end if
c
c     determine longitude distance relative to rwlon1 so it belongs to
c     the absolute interval 0 - 180
c
      argu1 = glor-rwlon1
      if(argu1.gt. 180.0) argu1 = argu1-360.0
      if(argu1.lt.-180.0) argu1 = argu1+360.0
c
c     1. get the help circle bb and angle alpha (legalize arguments)
c
      c2=glat*pi180
      c3=argu1*pi180
      arg2a = dcos(c2)*dcos(c3)
      arg2a = dmax1(arg2a,-c1)
      arg2a = dmin1(arg2a, c1)         
      bb    = dacos(arg2a)
c
      c4=hsign*glat*pi180
      arg2a = dsin(c4)/dsin(bb)
      arg2a = dmax1(arg2a,-c1)
      arg2a = dmin1(arg2a, c1)
      alpha = dasin(arg2a)
c
c     2. get pla and plo (still legalizing arguments)
c
      c5=rlat*pi180
      c6=hsign*rlat*pi180
      arg2a = dcos(c5)*dcos(bb)+
     +        dsin(c6)*dsin(c4)
      arg2a = dmax1(arg2a,-c1)
      arg2a = dmin1(arg2a, c1)         
      pla1   = dasin(arg2a)
c
      arg2a = dsin(bb)*dcos(alpha)/dcos(pla1)
      arg2a = dmax1(arg2a,-c1)
      arg2a = dmin1(arg2a, c1)
      plo1   = dasin(arg2a)
c
c    test for passage of the 90 degree longitude (duallity in plo)
c         get pla for which plo=90 when glat is the latitude
c
      arg2a = dsin(c4)/dsin(c6)
      arg2a = dmax1(arg2a,-c1)
      arg2a = dmin1(arg2a, c1)         
      pla90 = dasin(arg2a)
c
c         get help arc bb and angle alpha
c
      arg2a = dcos(c5)*dsin(pla90)
      arg2a = dmax1(arg2a,-c1)
      arg2a = dmin1(arg2a, c1)
      bb    = dacos(arg2a)

      arg2a = dsin(c4)/dsin(bb)
      arg2a = dmax1(arg2a,-c1)
      arg2a = dmin1(arg2a, c1)        
      alpha = dasin(arg2a)
c
c         get glolim - it is nesc. to test for the existence of solution
c
      argu2  = dcos(c2)*dcos(bb)/
     +            (1.-dsin(c4)*dsin(bb)*dsin(alpha))
      if(dabs(argu2).gt.c1) then
      glolim = 999.0
      else
        glolim = dacos(argu2)/pi180
      end if
c
c     modify (if nesc.) the plo solution
c
      if((abs(argu1).gt.glolim.and.glat.le.rlat).or.
     +   glat.gt.rlat) then
            plo1 = pi180*180.0 - plo1
      end if
c
c     the solution is symmetric so the direction must be if'ed
c
      if(argu1.lt.0.0) then
         plo1 = -plo1
      end if
c
c     convert the radians to degrees
c
      pla = pla1/pi180        
      plo = plo1/pi180
c
c     to obtain a rotated value (ie so x-axis in pol.ste. points east)
c     add 270 to longitude
c
      plo=amod(plo+270.0,360.0)
c
 2000 continue      
      return
      end                                  
c
c===============================================================================
c
      subroutine psll_2_psij(pslat,pslon,psi,psj)
c
ccc      implicit none
c
      include 'trigd.inc'
      real pslat,pslon,      !pol ste. lat, lon (deg n, deg +e)
     .       psi,psj,          !pol ste. i,j
     .       x,y,xmin,ymin,dx,dy,
     .       mag
c
      common /pscorner/xmin,ymin,dx,dy
c_______________________________________________________________________________
c
      mag=2./(1.+sind(pslat))*cosd(pslat)
      x=mag*cosd(pslon)
      y=mag*sind(pslon)
      psi=(x-xmin)/dx+1.
      psj=(y-ymin)/dy+1.
c
      return
      end
c
c===============================================================================
c
      subroutine psij_2_psll(psi,psj,pslat,pslon)
c
ccc      implicit none
c
      include 'trigd.inc'
      real psi,psj,          !pol ste. i,j
     .       pslat,pslon,      !pol ste. lat, lon (deg n, deg +e)
     .       x,y,dist,
     .       xmin,ymin,dx,dy 
c
      common /pscorner/xmin,ymin,dx,dy
c_______________________________________________________________________________
c
      x=(psi-1.)*dx+xmin
      y=(psj-1.)*dy+ymin
      dist=sqrt(x**2+y**2)
      if (dist .eq. 0) then
         pslat=90.
         pslon=-90.
      else
         pslat=atand(dist/2.)
         pslat=90.-2.*pslat
c
         if (x .eq. 0.) then
            pslon=90.
         else
            if (x .gt. 0.) then
               pslon=atand(y/x)
            else
               pslon=atand(y/x)+180.
            endif
         endif
      endif
c
      pslon=amod(pslon+450.,360.)
c
      return
      end
c
c===============================================================================
c
      subroutine psll_2_geoll(pla,plo,glat,glon)
c
c     convert polar stereographic coordinates to geographical lat/lon
c     ditto with the pol.ste. pole at rlat,wlon1 (these names are
c     used to be compatible to the ramsin-parameters)
c     longitude:   0 ; 360 positive east (on input)
c               -180 ; 180 positive east (on output)
c     latitude : -90 ;  90 posive on northern hemisphere
c     it is assumed that the polar stereographic coordinates have been
      include 'trigd.inc'
c     rotated to the standard format with 0 degrees longitude along wlon1
c
c     tsp 21 june 89
c
      integer nx,ny,nz           !no. of ps domain grid points
      real rlat,wlon1,rota,      !pol ste. std lat, lon and rotation
     .       sw(2),ne(2)           !sw lat, lon, ne lat, lon
      common /psgrid/nx,ny,nz,rlat,wlon1,rota,sw,ne
c_______________________________________________________________________________
c
c     set flag for n/s hemisphere
c
      c1 = 1.
      pi180 = asin(c1)/90.
c      
      if(rlat.ge.0.0) then
         hsign= 1.0
      else
         hsign=-1.0
      end if
c
c     test for a n/s pole case
c
      if(rlat.eq.90.0) then
       glat=pla
         glon=mod(plo+wlon1,360.0)
         go to 2000
      end if
      if(rlat.eq.-90.0) then
         glat=-pla
         glon=mod(plo+wlon1,360.0)
         go to 2000
      end if
c
c     test for longitude on 'greenwich or date line'
c
      if(plo.eq.0) then
         glat=rlat-90.0+pla
         if(glat.lt.-90.0) then
            glat=-180.0-glat
            glon=mod(wlon1+180.0,360.0)
         else
            glon=wlon1
         end if
         go to 2000
      end if      
      if(plo.eq.180.0) then
         glat=rlat+90.0-pla
         if(glat.gt.90.0) then
            glat=180.0-glat
            glon=mod(wlon1+180.0,360.0)
         else
            glon=wlon1
         end if
         go to 2000         
      end if
c
c     determine longitude distance relative to wlon1 so it belongs to
c     the absolute interval 0 - 180
c
      argu1=plo
      if(plo.gt.180.0) argu1 = plo-360.0
c
c     get the latitude, the help circle bb and the longitude by first
c     calculating the argument and legalize it - then take the inverse fct.
c
      if(hsign.gt.0.0) then
         arg2a = sin(pla*pi180)*sin(hsign*rlat*pi180)+
     +        cos(pla*pi180)*cos(rlat*pi180)*cos((180.0-argu1)*pi180)
      else
         arg2a = sin(pla*pi180)*sin(hsign*rlat*pi180)+
     +        cos(pla*pi180)*cos(rlat*pi180)*cos(argu1*pi180)
      end if
      arg2a = min(arg2a, 1.0)
      arg2a = max(arg2a,-1.0)
      glat  = hsign*asin(arg2a)
c
      if(hsign.gt.0.0) then
         arg2a = cos(rlat*pi180)*sin(pla*pi180)+
     +        sin(rlat*pi180)*cos(pla*pi180)*cos(argu1*pi180)
      else
         arg2a = cos(rlat*pi180)*sin(pla*pi180)+
     +       sin(-rlat*pi180)*cos(pla*pi180)*cos((180.0-argu1)*pi180)
      end if
      arg2a = min(arg2a, 1.0)
      arg2a = max(arg2a,-1.0)      
      bb    = acos(arg2a)
c
      arg2a = cos(glat)*cos(bb)/(1.0-sin(glat)**2)
      arg2a = min(arg2a, 1.0)
      arg2a = max(arg2a,-1.0)      
      glon  = acos(arg2a)
c     
c     convert the radians to degrees 
c
        glat = glat/pi180
        glon = glon/pi180
c
c       the solution is symmetric so the direction must be if'ed
c
        if(argu1.lt.0.0) then
           glon = 360.0-glon
        end if
        glon=amod(glon+wlon1,360.0)
c
 2000 continue
c
c     the resultant longitude must be in the interval from -180, 180
c      
      if(glon.gt.180.0) glon=glon-360.0
      return
      end
c
c===============================================================================
c
      subroutine ps_param(xmin,ymin,dx,dy)
c
ccc      implicit none
c
      include 'trigd.inc'
      real pslat,pslon,
     .       xmin,xmax,ymin,ymax,
     .       dx,dy,mag
c
      integer nx,ny,nz           !no. of ps domain grid points
      real lat0,lon0,rota,       !pol ste. std lat, lon and rotation
     .       sw(2),ne(2)           !sw lat, lon, ne lat, lon
      common /psgrid/nx,ny,nz,lat0,lon0,rota,sw,ne
c
c_______________________________________________________________________________
c
      call geoll_2_psll(sw(1),sw(2),pslat,pslon)
      mag=2./(1.+sind(pslat))*cosd(pslat)
      xmin=mag*cosd(pslon)
      ymin=mag*sind(pslon)
      call geoll_2_psll(ne(1),ne(2),pslat,pslon)
      mag=2./(1.+sind(pslat))*cosd(pslat)
      xmax=mag*cosd(pslon)
      ymax=mag*sind(pslon)
      dx=(xmax-xmin)/float(nx-1)
      dy=(ymax-ymin)/float(ny-1)
c
      return
      end
c
c===============================================================================
c===============================================================================
c
      subroutine conical_equidistant
c
c *** routines to convert from geographical lat, lon to 
c        conical-equidistant grid i, j vice-versa.
c     equations obtained from adrian marroquin and tom black (ncep).
c     snook (12/20/95)
c
ccc      implicit none
c
      integer np,n
c
      real glat(np),glon(np),      !earth lat, lon (deg n, deg +e)
     .       celat,celon,            !con eq. lat, lon (deg n, deg +e)
     .       cei(np),cej(np)         !con eq. i,j
c
c===============================================================================
c
      entry latlon_2_coneqij(np,glat,glon,cei,cej)
c_______________________________________________________________________________
c
      do n=1,np
         call geoll_2_coneqll(glat(n),glon(n),celat,celon)
         call coneqll_2_coneqij(celat,celon,cei(n),cej(n))
      enddo
      return
c
c===============================================================================
c
      entry coneqij_2_latlon(np,cei,cej,glat,glon)
c_______________________________________________________________________________
c
      do n=1,np
         call coneqij_2_coneqll(cei(n),cej(n),celat,celon)
         call coneqll_2_geoll(celat,celon,glat(n),glon(n))
      enddo
      return
c
      end
c
c===============================================================================
c
      subroutine geoll_2_coneqll(glat,glon,celat,celon)
c
ccc      implicit none
c
      include 'trigd.inc'
      real glat,glon,      !earth lat, lon (deg n, deg +e)
     .       celat,celon,    !con eq. lat, lon (deg n, deg +e)
     .       x,y,z
c
      integer nx,ny,nz
      real lat0,lon0,dphi,dlam
      common /coneqgrid/nx,ny,nz,lat0,lon0,dphi,dlam
c_______________________________________________________________________________
c
      x=cosd(lat0)*cosd(glat)*cosd(lon0-glon)+sind(lat0)*sind(glat)
      y=-cosd(glat)*sind(lon0-glon)
      z=-sind(lat0)*cosd(glat)*cosd(lon0-glon)+cosd(lat0)*sind(glat)
c
      celat=atand(z/(x**2+y**2)**0.5)
      celon=atand(y/x)
c
      return
c
      end
c
c===============================================================================
c
      subroutine coneqll_2_coneqij(celat,celon,cei,cej)
c
ccc      implicit none
c
      real celat,celon,      !con eq. lat, lon (deg n, deg +e)
     .       cei,cej           !con eq. i,j
c
      integer nxt
c
      integer nx,ny,nz
      real lat0,lon0,dphi,dlam
      common /coneqgrid/nx,ny,nz,lat0,lon0,dphi,dlam
c_______________________________________________________________________________
c
      nxt=2*nx-1
      cei=float((nxt-1)/2+1)+celon/dlam
      cej=float((ny -1)/2+1)+celat/dphi
c
      return
      end
c
c===============================================================================
c
      subroutine coneqij_2_coneqll(cei,cej,celat,celon)
c
ccc      implicit none
c
      real cei,cej,          !con eq. i,j
     .       celat,celon       !con eq. lat, lon (deg n, deg +e)
c
      integer nyt
c
      integer nx,ny,nz
      real lat0,lon0,dphi,dlam
      common /coneqgrid/nx,ny,nz,lat0,lon0,dphi,dlam
c_______________________________________________________________________________
c
      nyt=ny/2+1
      celon=(cei-float(nx ))*dlam
      celat=(cej-float(nyt))*dphi
c
      return
      end
c
c===============================================================================
c
      subroutine coneqll_2_geoll(celat,celon,glat,glon)
c
ccc      implicit none
c
      include 'trigd.inc'
      real celat,celon,    !con eq. lat, lon (deg n, deg +e)
     .       glat,glon       !earth lat, lon (deg n, deg +e)
c
      integer nx,ny,nz
      real lat0,lon0,dphi,dlam
      common /coneqgrid/nx,ny,nz,lat0,lon0,dphi,dlam
c_______________________________________________________________________________
c
      glat=asind(sind(celat)*cosd(lat0)+
     .           cosd(celat)*sind(lat0)*cosd(celon))
      if (celon .lt. 0.) then
         glon=lon0-acosd(cosd(celat)*cosd(celon)/cosd(glat)/cosd(lat0)-
     .                   tand(glat)*tand(lat0))
      else
         glon=lon0+acosd(cosd(celat)*cosd(celon)/cosd(glat)/cosd(lat0)-
     .                   tand(glat)*tand(lat0))
      endif
c
      return
      end
c
c===============================================================================
c===============================================================================
c
      subroutine lambert_conformal
c
c *** routines to convert from geographical lat, lon to 
c        lambert-conformal grid i, j and vice-versa.
c     equations obtained from ncar graphics documentation.
c     snook (12/20/95)
c
ccc      implicit none
c
      include 'trigd.inc'
      integer np,n
c
      real glat(np),glon(np),    !earth lat (deg n), lon (deg +e)
     .       lci(np),lcj(np),      !lambert-confomal i, j
     .       s,cone,r,
     .       xmin,ymin,dx,dy,x,y
c
      real lat1,lat2,lon0,       !lambert-conformal std lat1, lat2, lon
     .       sw(2),ne(2)           !sw lat, lon, ne lat, lon
      integer nx,ny,nz           !no. of lc domain grid points
      common /lcgrid/nx,ny,nz,lat1,lat2,lon0,sw,ne
c
c===============================================================================
c
      entry latlon_2_lcij(np,glat,glon,lci,lcj)
c_______________________________________________________________________________
c
      call lc_param(s,cone,xmin,ymin,dx,dy)
c
      do n=1,np
         r=(tand(45.-s*glat(n)/2.))**cone
         x=r*sind(cone*(glon(n)-lon0))
         y=-s*r*cosd(cone*(glon(n)-lon0))
         lci(n)=(x-xmin)/dx+1
         lcj(n)=(y-ymin)/dy+1
      enddo
c
      return
c
c===============================================================================
c
      entry lcij_2_latlon(np,lci,lcj,glat,glon)
c_______________________________________________________________________________
c
      call lc_param(s,cone,xmin,ymin,dx,dy)
c
      do n=1,np
         x=(lci(n)-1)*dx+xmin
         y=(lcj(n)-1)*dy+ymin
         glon(n)=lon0+atand(-s*x/y)/cone
         glat(n)=(90.-
     .            2.*atand((x/sind(cone*(glon(n)-lon0)))**(1./cone)))/s
      enddo
c
      return
c
      end
c
c===============================================================================
c
      subroutine lc_param(s,cone,xmin,ymin,dx,dy)
c
ccc      implicit none
c
      include 'trigd.inc'
      real s,cone,r,
     .       xmin,xmax,ymin,ymax,
     .       dx,dy
c
      real lat1,lat2,lon0,       !lambert-conformal std lat1, lat2, lon
     .       sw(2),ne(2)           !sw lat, lon, ne lat, lon
      integer nx,ny,nz           !no. of lc domain grid points
      common /lcgrid/nx,ny,nz,lat1,lat2,lon0,sw,ne
c_______________________________________________________________________________
c
      if (lat1 .ge. 0.) then
         s=1.
      else
         s=-1.
      endif
      if (lat1 .ne. lat2) then
         cone=alog(cosd(lat1)/cosd(lat2))/
     .        alog(tand(45.-s*lat1/2.)/tand(45.-s*lat2/2.))
      else
         cone=cosd(90.-s*lat1)
      endif
c
      r=(tand(45.-s*sw(1)/2.))**cone
      xmin=r*sind(cone*(sw(2)-lon0))
      ymin=-s*r*cosd(cone*(sw(2)-lon0))
      r=(tand(45.-s*ne(1)/2.))**cone
      xmax=r*sind(cone*(ne(2)-lon0))
      ymax=-s*r*cosd(cone*(ne(2)-lon0))
      dx=(xmax-xmin)/float(nx-1)
      dy=(ymax-ymin)/float(ny-1)
c
      return
      end
c
c===============================================================================
c===============================================================================
c
      subroutine lat_lon
c
c *** routines to convert from geographical lat, lon to 
c        lat-lon grid i, j and vice-versa.
c     snook (11/5/96)
c     smart/pw li (04/1/03)
c
      integer np,n
      integer nx,ny,nz             !no. of ll domain grid points
 
      real   glat(np),glon(np)     !earth lat, lon (deg n, deg +e)
      real   lli(np),llj(np)       !lat-lon grid i,j
      real   diff
c     real   sw(2),ne(2)
c     common /llgrid/nx,ny,sw,ne,cgrddef

      real   lat0,lon0,dlat,dlon     !sw corner lat, lon, lat, lon spacing
      character*1 cgrddef
      common /llgrid/nx,ny,nz,lat0,lon0,dlat,dlon,cgrddef
c
c===============================================================================
c
      entry latlon_2_llij(np,glat,glon,lli,llj)
c_______________________________________________________________________________
c

c pw li code for integration later on.
c     dlat=(ne(2)-sw(2))/(nx-1)
c     dlon=(ne(1)-sw(1))/(ny-1)
c     do i=1,n
c        lli(i)=(glon(i)-sw(2))/dlon + 1.
c        llj(i)=(glat(i)-sw(1))/dlat + 1.
c     enddo
 
      dlond=dlon
      dlatd=dlat

      rlon_width = float(nx) * dlond
      write(6,*)' latlon_2_llij: rlon_width = ',rlon_width
      
      if(cgrddef.eq.'s')then
         do n=1,np
            diff=glon(n)-lon0
            if(rlon_width .lt. 170.)then ! regional
               if (diff .lt. -180.) diff=diff+360.
            else ! global (e.g. gfs)
               if (diff .lt. 0.)   diff=diff+360.
               if (diff .ge. 360.) diff=diff-360.
            endif
            lli(n)=diff/dlond+1.
            llj(n)=(glat(n)-lat0)/dlatd+1.
         enddo
      elseif(cgrddef.eq.'n')then
         do n=1,np
            diff=glon(n)-lon0
            if(rlon_width .lt. 170.)then ! regional
               if (diff .lt. -180.) diff=diff+360.
            else ! global (e.g. gfs)
               if (diff .lt. 0.)   diff=diff+360.
               if (diff .ge. 360.) diff=diff-360.
            endif
!           write(6,*)glon(n),lon0,diff
            lli(n)=diff/dlon+1.
            llj(n)=(lat0-glat(n))/dlat+1.
         enddo

      else
         print*,'you must specify whether the standard
     .           lat is southern or northern boundary'
      endif

      return
c
c===============================================================================
c
      entry llij_2_latlon(np,lli,llj,glat,glon)
c_______________________________________________________________________________
c
c note: if lat0=northern boundary of ll grid set dlat=-dlat
      do n=1,np
         glon(n)=(lli(n)-1.)*dlon+lon0
         glat(n)=(llj(n)-1.)*dlat+lat0
      enddo
      return
c
      end
c
c===============================================================================
c===============================================================================
c
      subroutine fixed_grid
c
      integer np,n,isat_fgf
      integer nx,ny,nz             !no. of fx domain grid points
 
      real   glat(np),glon(np)     !earth lat, lon (deg n, deg +e)
      real   fxi(np),fxj(np)       !lat-lon grid i,j
      real   diff
c     real   sw(2),ne(2)
c     common /llgrid/nx,ny,sw,ne,cgrddef

      real   lat0,lon0,dlat,dlon     !sw corner lat, lon, lat, lon spacing
      real*8 dlond,dlatd,glond,glatd
      real*8 scale_x,scale_y,offset_x,offset_y
      real*8 sub_lon_degrees,fgf_x,fgf_y,xmin,ymin
      character*1 cgrddef

!     x0,y0,dx,dy are in microradians
      common /fxgrid/offset_x,offset_y,scale_x,scale_y,xmin,ymin
     1              ,sub_lon_degrees
c
c===============================================================================
c
      entry latlon_2_fxij(np,glat,glon,fxi,fxj)
c_______________________________________________________________________________
c

      dlond=dlon
      dlatd=dlat
      isat_fgf = 1

      do n=1,np
         glond = glon(n)
         glatd = glat(n)
         call earth_to_fgf(isat_fgf,glond,glatd,scale_x,offset_x,scale_y
     1                    ,offset_y,sub_lon_degrees,fgf_x,fgf_y)
         fxi(n) = +fgf_x - xmin
         fxj(n) = -fgf_y - ymin
         if(n .le. 2)then
            write(6,*)
            write(6,*)'latlon_2_fxij first grid points ',n
            write(6,*)'-----------------------------------'
            write(6,*)'glatd = ',glatd
            write(6,*)'glond = ',glond
            write(6,*)'scale_x = ',scale_x
            write(6,*)'scale_y = ',scale_y
            write(6,*)'offset_x = ',offset_x
            write(6,*)'offset_y = ',offset_y
            write(6,*)'sub_lon_degrees = ',sub_lon_degrees
            write(6,*)'fgf_x = ',fgf_x
            write(6,*)'fgf_y = ',fgf_y
            write(6,*)'xmin = ',xmin
            write(6,*)'ymin = ',ymin
            write(6,*)'fxi(n) = ',fxi(n)
            write(6,*)'fxj(n) = ',fxj(n)
            write(6,*)
         endif
      enddo ! n

      return
c
c===============================================================================
c
      entry fxij_2_latlon(np,fxi,fxj,glat,glon)
c_______________________________________________________________________________
c
c note: if lat0=northern boundary of fx grid set dlat=-dlat
      do n=1,np
         glon(n)=(fxi(n)-1.)*dlon+lon0
         glat(n)=(fxj(n)-1.)*dlat+lat0
      enddo
      return
c
      end
c
c===============================================================================
c===============================================================================
c

      subroutine cylindrical_equidistant
c
c routine to convert from cylindrical equidistant to grid ri/rj
c used for mapping the wsi radar data to laps.
c
c     j. smart 11-19-98   original working version
c                         basic equation set from
c                         map projections used by the u.s. geological survey
c                         (snyder, j.p. 1983 geological survey bulletin - 1532)
c
      include 'trigd.inc'
      implicit none
      integer np,n
 
      real glat(np),glon(np),      !earth lat, lon (deg n, deg +e)
     .       lli(np),llj(np),        !lat-lon grid i,j
     .       nw(2),se(2)             !nw grid lat/lon; se grid lat/lon

      double precision     diff,x,y
     .                    ,xmin,ymin
     .                    ,xmax,ymax
     .                    ,pi,dg2rd

      double precision     coslatc
     .                    ,dx,dy
     .                    ,drlatc,drlonc
     .                    ,dglat, dglon

!     double precision     dcosd

      real  r
 
      integer nx,ny,nz             !no. of ll domain grid points
      real  rlatc,rlonc          !grid center lat, lon

c     real coslatc

      data r/6.3712e6/             !this earth radius is = to that in lapsgrid.f
      common /cegrid/nx,ny,nz,nw,se,rlatc,rlonc
c
c===============================================================================
c
      entry latlon_2_ceij(np,glat,glon,lli,llj)
c_______________________________________________________________________________
c
c note: wsi grid (1,1) is nw corner. adjustment made for rj using ny since
c       equation set assumes y-axis is on equator.
c
      pi=acos(-1.0)
      dg2rd=pi/180.
      drlatc=rlatc
      drlonc=rlonc

      xmin=r*((nw(2)-drlonc)*dg2rd)*dcosd(drlatc)
      ymin=r*se(1)*dg2rd
      xmax=r*((se(2)-drlonc)*dg2rd)*dcosd(drlatc)
      ymax=r*nw(1)*dg2rd

      diff=xmax-xmin
      if(diff.lt.0.000001)then
         xmax=dabs(xmin)
      endif

      dx=(xmax-xmin)/(nx-1)
      dy=(ymax-ymin)/(ny-1)

      coslatc=dcosd(drlatc)
      do n=1,np
         diff=glon(n)-drlonc
         if (diff .lt.-180.) diff=diff+360.
         if (diff .ge. 180.) diff=diff-360.
         diff=diff*dg2rd
         x=r*diff*coslatc
         y=r*(glat(n)*dg2rd)
         lli(n)=(x-xmin)/dx + 1.
         llj(n)=float(ny)-((y-ymin)/dy) + 1.
      enddo

      return
c
c===============================================================================
c
      entry ceij_2_latlon(np,lli,llj,glat,glon)
c
c equation set: rlat (phi) = y/r
c               rlon (lambda) = rlonc + x/(r cos(rlatc))
c        where:
c               r= earth radius
c               rlonc= central meridian
c               rlatc= central parallel

      pi=acos(-1.0)
      dg2rd=pi/180.
      drlatc=rlatc
      drlonc=rlonc

      xmin=r*((nw(2)-drlonc)*dg2rd)*dcosd(drlatc)
      ymin=r*se(1)*dg2rd
      xmax=r*((se(2)-drlonc)*dg2rd)*dcosd(drlatc)
      ymax=r*nw(1)*dg2rd

      diff=xmax-xmin
      if(diff.lt.0.000001)then
         xmax=dabs(xmin)
      endif

      dx=(xmax-xmin)/(nx-1)
      dy=(ymax-ymin)/(ny-1)

      coslatc=dcosd(drlatc)

      do n=1,np

         dglon=rlonc*dg2rd +
     &         (xmin+(lli(n)+1.0)*dx)/(r*coslatc)
         dglat=(ymin+(ny-llj(n)+1.0)*dy)/r

         glon(n)=dglon/dg2rd
         glat(n)=dglat/dg2rd

      enddo

      return
      end
c
c ===========================================================================
c
      subroutine latlon_2_mcij(n,rlat,rlon,ri,rj)

cc     implicit none
      include 'trigd.inc'

      integer n,i
      integer nx,ny
      real rlat(n),rlon(n),ri(n),rj(n)
      real rlonc,dlon,rlatc,dlat
      real sw(2),ne(2)
      real x,y
      real xmax,ymax
      real xmin,ymin
      real dx,dy
      real r,pi
      real deg2rad
      common /mcgrid/rlonc,rlatc,nx,ny,sw,ne

      call get_earth_radius(r,istatus)

      pi=acos(-1.)
      deg2rad=pi/180.

      ymax=r*(log(tan(pi/4.+0.5*(ne(1)-rlatc)*deg2rad)))
      ymin=r*(log(tan(pi/4.+0.5*(sw(1)-rlatc)*deg2rad)))

      xmax=r*(deg2rad*(ne(2)-rlonc))
      xmin=r*(deg2rad*(sw(2)-rlonc))

      dx=(xmax-xmin)/float(nx-1)
      dy=(ymax-ymin)/float(ny-1)

      do i=1,n
         dlon=rlon(i)-rlonc
         dlat=rlat(i)-rlatc
         if(dlon.gt.180.)dlon=dlon-360.
         if(dlon.lt.-180.)dlon=dlon+360.
         x=r*(deg2rad*dlon)
         y=r*(log(tan(pi/4.+0.5*(dlat*deg2rad))))
         ri(i)=(x-xmin)/dx + 1.
         rj(i)=(y-ymin)/dy + 1.
      enddo

      return
      end
c
c===============================================================================
c
      subroutine latlon_2_npij(n,rlat,rlon,ri,rj)
c
c this is conversion of lat/lon to ri/rj in a
c lat-lon grid. very similar to latlon_2_llij.
c
      integer n,i
      integer nx,ny
      real rlat(n),rlon(n),ri(n),rj(n)
      real sw(2),ne(2)
      real dx,dy
      common /npgrid/nx,ny,sw,ne

c     print *, ' inside latlon_2_npij'
c     print *, nx,ny,dx,dy
c     print *, sw(1), sw(2),ne(1),ne(2)

      dx=(ne(2)-sw(2))/(nx-1)
      dy=(ne(1)-sw(1))/(ny-1)
 
      do i=1,n
        ri(i)=(rlon(i)-sw(2))/dx + 1.
        rj(i)=(rlat(i)-sw(1))/dy + 1.
      enddo
 
      return
      end

c -------------------------------------------------------

      subroutine init_gridconv_cmn(gproj,nxbg,nybg,nzbg
     &,dlat,dlon,cenlat,cenlon,lat0,lat1,lon0
     &,sw1,sw2,ne1,ne2,cgrddef,istatus)
c
c js 4-01
c
      implicit none

      character*(*)  gproj
      character*(*)  cgrddef
      integer        istatus
      integer        nxbg,nybg,nzbg
      real           dlat,dlon
      real           lat0,lat1
      real           lon0,lon1,lon2
      real           sw1,ne1
      real           sw2,ne2
      real           cenlat,cenlon

c
c *** common block variables for lat-lon grid.
c
      integer   nx_ll,ny_ll,nz_ll
      real    lat0_ll,lon0_ll,d_lat,d_lon
      character*1 cgrddef_ll
      common /llgrid/nx_ll,ny_ll,nz_ll,lat0_ll,lon0_ll
     &,d_lat,d_lon,cgrddef_ll
c
c *** common block variables for lambert-conformal grid.
c
      integer   nx_lc,ny_lc,nz_lc
      real    lat1_lc,lat2_lc,lon0_lc,sw_lc(2),ne_lc(2)
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1_lc,lat2_lc
     &,lon0_lc,sw_lc,ne_lc
c
c *** common block variables for cyclindrical equidistant grid.
c
      integer   nx,ny,nz
      real    rlatc,rlonc,nw(2),se(2),dx,dy
      common /cegrid/nx,ny,nz,nw,se,rlatc,rlonc
c
c *** common block variables for polar stereographic grid.
c
      integer nx_ps,ny_ps,nz_ps    !no. of ps domain grid points
      real lat0_ps,lon0_ps,rota  !pol ste. std lat, lon and rotation
     .      ,sw_ps(2),ne_ps(2)     !sw lat, lon, ne lat, lon
      common /psgrid/nx_ps,ny_ps,nz_ps,lat0_ps,lon0_ps
     .              ,rota,sw_ps,ne_ps

      integer nx_np,ny_np
      real    sw_np(2),ne_np(2)
      common /npgrid/nx_np,ny_np,sw_np,ne_np
c
c *** common block variables for mercator grid.
c
      integer nx_mc,ny_mc          !no. of domain grid points
      real    sw_mc(2),ne_mc(2)    !sw lat, lon, ne lat, lon
      real    rlatc_mc,rlonc_mc    !center lat/lon of domain
      common /mcgrid/rlonc_mc,rlatc_mc,nx_mc,ny_mc,sw_mc,ne_mc

      if(gproj.eq.'lc')then
         nx_lc=nxbg
         ny_lc=nybg
         nz_lc=nzbg
         lat1_lc=lat0
         lat2_lc=lat1
         lon0_lc=lon0
         sw_lc(1)=sw1
         sw_lc(2)=sw2
         ne_lc(1)=ne1
         ne_lc(2)=ne2
      endif

      if(gproj.eq.'ll')then
         nx_ll=nxbg
         ny_ll=nybg
         nz_ll=nzbg
         lat0_ll=lat0
         lon0_ll=lon0
         d_lat=dlat
         d_lon=dlon
         cgrddef_ll=cgrddef
      endif

      if(gproj.eq.'le')then
         nx=nxbg
         ny=nybg
         nz=nzbg
         rlatc=cenlat
         rlonc=cenlon
         se(1)=sw1
         se(2)=sw2
         nw(1)=ne1
         nw(2)=ne2
      endif

      if(gproj.eq.'ps')then
         nx_ps=nxbg
         ny_ps=nybg
         nz_ps=nzbg
         lat0_ps=lat0
         lon0_ps=lon0
         sw_ps(1)=sw1
         sw_ps(2)=sw2
         ne_ps(1)=ne1
         ne_ps(2)=ne2
      endif

      if(gproj.eq.'mc')then
         nx_mc=nxbg
         ny_mc=nybg
         rlatc_mc=cenlat
         rlonc_mc=cenlon
         sw_mc(1)=sw1
         sw_mc(2)=sw2
         ne_mc(1)=ne1
         ne_mc(2)=ne2
      endif

      return
      end
c
c===============================================================================
c
      subroutine init_hinterp(nx_bg,ny_bg,nx_laps,ny_laps,gproj,
     .     lat,lon,grx,gry,bgmodel,cmodel,wrapped)

c
      implicit none
c
      integer nx_bg,ny_bg,nx_laps,ny_laps,bgmodel
c
      real lat(nx_laps,ny_laps),lon(nx_laps,ny_laps),
     .       grx(nx_laps,ny_laps),gry(nx_laps,ny_laps) 

      logical wrapped
c
      integer i,j,k
      integer istatus
c
      character*132 cmodel
      character*2   gproj

      integer lenc

      integer nxc,nyc,nzc
      integer nx,ny
      logical lprintmessage
      real sw(2),ne(2),rota,lat0,lon0
      real nw(2),se(2),rlatc,rlonc
      real tolx,toly
      real grxdifsum1,grxdifsum2
      real grydifsum1,grydifsum2
      real r_missing_data
c     parameter (tol=0.10)
      common /psgrid/nxc,nyc,nzc,lat0,lon0,rota,sw,ne
c________________________________________________________________________________
c
      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne. 1)then
         print*,'error getting r_missing_data - init_hinterp'
         return
      endif
c
c *** determine location of laps grid point in background data i,j space.
c
      if (gproj .eq. 'ps') then

c      print*,nxc,nyc,nzc,lat0,lon0,rota,sw,ne

         call latlon_2_psij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'lc') then
         call latlon_2_lcij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'ce') then
         call latlon_2_coneqij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'll') then
         call latlon_2_llij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'le') then
         call latlon_2_ceij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'np') then
         call latlon_2_npij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'mc') then
         call latlon_2_mcij(nx_laps*ny_laps,lat,lon,grx,gry)
      else
         print*,"error: unknown gproj spec in gridconv ",gproj
      endif
c
c *** check that all laps grid points are within the background data coverage.
c

c set tolerance based upon the grid spacing as a function of grx/gry
      grxdifsum1=0.0
      grxdifsum2=0.0
      grydifsum1=0.0
      grydifsum2=0.0
      do j=1,ny_laps
         grxdifsum1=grxdifsum1+(grx(2,j)-grx(1,j))
         grxdifsum2=grxdifsum2+(grx(nx_laps-1,j)-grx(nx_laps,j))
      enddo
      do i=1,nx_laps
         grydifsum1=grydifsum1+(gry(i,2)-gry(i,1))
         grydifsum2=grydifsum2+(gry(i,ny_laps-1)-gry(i,ny_laps))
      enddo

      tolx=(abs(grxdifsum1)/ny_laps+abs(grxdifsum2)/ny_laps)*0.5
      toly=(abs(grydifsum1)/nx_laps+abs(grydifsum2)/nx_laps)*0.5

      print*,'horiz mapping tolerance x/y: ',tolx,toly,wrapped
      lprintmessage=.true.
c
c *** first, check for wrapping if a global data set.
c
cwni wni comment:  this section not needed because
cwni  the latlon_2_llij already handles this
      call s_len(cmodel,lenc)

      if ( bgmodel .eq. 6 .or. 
     .     bgmodel .eq. 8) then
cwni     do j=1,ny_laps
cwni        do i=1,nx_laps
cwni           if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
cwni           if (grx(i,j) .gt. nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
cwni           if (gry(i,j) .lt. 1.) then
cwni              gry(i,j)=2.-gry(i,j)
cwni              grx(i,j)=grx(i,j)-float(nx_bg/2)
cwni              if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
cwni              if (grx(i,j).gt.nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
cwni           endif
cwni           if (gry(i,j) .gt. ny_bg) then
cwni              gry(i,j)=float(2*ny_bg)-gry(i,j)
cwni              grx(i,j)=grx(i,j)-float(nx_bg/2)
cwni              if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
cwni              if (grx(i,j).gt.nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
cwni           endif
cwni           if (grx(i,j) .lt. 1) then
cwni              grx(i,j)=grx(i,j)+1.
cwni           endif
cwni        enddo
cwni     enddo
cwni  elseif(bgmodel.eq.4.and.cmodel(1:lenc).eq.'avn_sbn_cyleq')then
      elseif( (bgmodel.eq.4.and.cmodel(1:lenc).eq.'avn_sbn_cyleq') .or. !wni
     .        (bgmodel .eq. 10 .and. cmodel .eq. 'gfs_iso')) then      !wni

cwni     do j=1,ny_laps
cwni        do i=1,nx_laps
cwni           if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
cwni           if (grx(i,j) .gt. nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
cwni           if (gry(i,j) .lt. 1.) then
cwni              gry(i,j)=2.-gry(i,j)
cwni              grx(i,j)=grx(i,j)-float(nx_bg/2)
cwni              if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
cwni              if (grx(i,j).gt.nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
cwni           endif
cwni           if (gry(i,j) .gt. ny_bg) then
cwni              gry(i,j)=float(2*ny_bg)-gry(i,j)
cwni              grx(i,j)=grx(i,j)-float(nx_bg/2)
cwni              if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
cwni              if (grx(i,j).gt.nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
cwni           endif
cwni        enddo
cwni     enddo

      elseif( wrapped ) then ! e.g. global dataset
         if (grx(i,j) .lt. 1 .or. grx(i,j) .gt. nx_bg .or.
     .       gry(i,j) .lt. 1 .or. gry(i,j) .gt. ny_bg) then
           do j=1,ny_laps
             do i=1,nx_laps
               if(lprintmessage)then
                  print*,'domain gridpt outside of bkgd data coverage.'
                  print*,'   data i,j,lat,lon - ',i,j,lat(i,j),lon(i,j)
                  print*,'   grx, gry:',grx(i,j),gry(i,j)
                  lprintmessage=.false.
c                 stop 'init_hinterp'
               endif
             enddo
           enddo
         endif

c ****** if not a global data set, then check that laps domain is fully
c           within background domain.
c
      else
         do j=1,ny_laps
            do i=1,nx_laps
c
c laps must fit into model grid which must also fit into laps grid thus we
c introduce a small fudge factor on the grid boundaries.
c               

               if(grx(i,j).gt.1.-tolx) grx(i,j) = max(1.,grx(i,j))
               if(gry(i,j).gt.1.-toly) gry(i,j) = max(1.,gry(i,j))

               if(grx(i,j).lt.float(nx_bg)+tolx) 
     +              grx(i,j) = min(float(nx_bg)-tolx,grx(i,j))
               if(gry(i,j).lt.float(ny_bg)+toly) 
     +              gry(i,j) = min(float(ny_bg)-toly,gry(i,j))

               if (grx(i,j) .lt. 1 .or. grx(i,j) .gt. nx_bg .or.
     .             gry(i,j) .lt. 1 .or. gry(i,j) .gt. ny_bg) then

                  grx(i,j) = r_missing_data
                  gry(i,j) = r_missing_data

           if(lprintmessage)then
              print*,'domain gridpt outside of bkgd data coverage.'
              print*,'   data i,j,lat,lon - ',i,j,lat(i,j),lon(i,j)
              print*,'   grx, gry:',grx(i,j),gry(i,j)
              lprintmessage=.false.
c                 stop 'init_hinterp'
           endif

               endif
            enddo
         enddo
      endif
c
      return
      end

      subroutine hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,nz,
     .     grx,gry,fvi,flaps,wrapped)
c
      implicit none
c
      integer nx_bg,ny_bg,nx_laps,ny_laps,nz        

      logical wrapped
c
c *** input vertically interpolated field - 2d should be slightly faster than
c     'hinterp_field_3d'
c *** output laps field
c
      real fvi(nx_bg,ny_bg,nz),
     .       flaps(nx_laps,ny_laps,nz),
     .       grx(nx_laps,ny_laps),gry(nx_laps,ny_laps)
c
      integer i,j,k
      integer istatus
      real    r_missing_data
c

      write(6,*)' subroutine hinterp_field'

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'error getting r_missing_data - hinterp_field'
         return
      endif
c
c *** horizontally interpolate variable.
c
      do k=1,nz
         do j=1,ny_laps
         do i=1,nx_laps
            if(grx(i,j).lt.r_missing_data.and. 
     .         gry(i,j).lt.r_missing_data)then
               call gdtost_i(fvi(1,1,k),nx_bg,ny_bg,
     .              grx(i,j),gry(i,j),flaps(i,j,k),wrapped)
            else
               flaps(i,j,k)=r_missing_data
            endif
         enddo
         enddo
      enddo
c
      return
      end

      subroutine hinterp_field_3d(nx_bg,ny_bg,nx_laps,ny_laps,nz,
     .     grx,gry,fvi,flaps,wrapped)
c
!     implicit none
c
      include 'bgdata.inc'

      integer nx_bg,ny_bg,nx_laps,ny_laps,nz        

      dimension r(4,nz),scr(4,nz),staval(nz)

      logical wrapped
c
c *** input vertically interpolated field - optimized for 3d.
c *** output laps field
c
      real fvi(nx_bg,ny_bg,nz),
     .       flaps(nx_laps,ny_laps,nz),
     .       grx(nx_laps,ny_laps),gry(nx_laps,ny_laps)
c
      integer i,j,k
      integer istatus
      real    r_missing_data
c
      write(6,*)' subroutine hinterp_field_3d ',wrapped

      iwrite = 0

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'error getting r_missing_data - hinterp_field'
         return
      endif

      write(6,*)' center input column',fvi(nx_bg/2,ny_bg/2,:)

      write(6,*)' center grx,gry ',grx(nx_laps/2,ny_laps/2)
     1                            ,gry(nx_laps/2,ny_laps/2)
c
c *** horizontally interpolate variable.
c
      ix = nx_bg
      iy = ny_bg

      do jl=1,ny_laps
      do il=1,nx_laps
        if(grx(il,jl).lt.r_missing_data.and. 
     .     gry(il,jl).lt.r_missing_data)then
!         call gdtost_i(fvi(1,1,k),nx_bg,ny_bg,
!    .         grx(i,j),gry(i,j),flaps(i,j,k),wrapped)
!         subroutine gdtost_i(ab,ix,iy,stax,stay,staval,wrapped)

          stax = grx(il,jl)
          stay = gry(il,jl)
            
!         dimension ab(ix,iy),r(4),scr(4)
!         logical wrapped ! wni added
!         include 'bgdata.inc'
c_______________________________________________________________________________
c
          iy1=int(stay)-1
          if(stay.lt.1.0)iy1=1.0
          iy2=iy1+3
          ix1=int(stax)-1
          if(stax.lt.1.0) then 
            if (wrapped) then   ! wni
              ix1 = ix1 + ix
            else ! wni
              ix1=1.0
            endif ! wni
          endif ! wni
          ix2=ix1+3
          fiym2=float(iy1)-1
          fixm2=float(ix1)-1

            ii=0

            do iii=ix1,ix2
               i=iii
c
c ****** account for wrapping around effect of global data at greenwich.
c
               if (wrapped) then    ! wni
                  if (i .lt. 1) i=i+ix  
                  if (i .gt. ix) i=i-ix
               endif

               ii=ii+1
               if (i .ge. 1 .and. i .le. ix) then 
                  jj=0
                  do jjj=iy1,iy2
                     j=jjj
c
c ************ account for n-s wrapping effect of global data.
c
                     if (wrapped) then  !wni
                        if (j .lt. 1) then
                           j=2-j
                           i=i-ix/2
                           if (i .lt. 1) i=i+ix
                           if (i .gt. ix) i=i-ix
                        endif
                        if (j .gt. iy) then
                           j=2*iy-j
                           i=i-ix/2
                           if (i .lt. 1) i=i+ix
                           if (i .gt. ix) i=i-ix
                        endif
                     endif

                     jj=jj+1

                     if (j .ge. 1 .and. j .le. iy) then
                        r(jj,:)=fvi(i,j,:) ! ab(i,j)
                     else
                        r(jj,:)=missingflag
                     endif
                  enddo
                  yy=stay-fiym2

                  if (yy .eq. 2.0) then
                    scr(ii,:)=r(2,:)
                  else
      
                    staval(:)=missingflag
                    xxx=yy
                    wt1=(xxx-3.)/(-1.)
                    wt2=1.0-wt1
                    yz22_temp=wt1*(xxx-4.)/(-2.)
                    yz23_temp=wt2*(xxx-4.)/(-1.)
                    yz24_temp=(xxx-2.)*(xxx-3.)/2.

                    yz11_temp=(xxx-2.)*(xxx-3.)/(2.)
                    yz12_temp=wt1*(xxx-1.)/(1.)
                    yz13_temp=wt2*(xxx-1.)/(2.)

                    do k=1,nz

!                    call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii)) (now inlined)
!                    subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
                     y1 = r(1,k)
                     y2 = r(2,k)
                     y3 = r(3,k)
                     y4 = r(4,k)
c
                     yyy=missingflag
                     if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
c
                       if (y4 .lt. 1.e19) then
                         yz22=yz22_temp
                         yz23=yz23_temp
                         yz24=yz24_temp
                       else
                         yz22=wt1
                         yz23=wt2
                         yz24=0.0
                       endif
c
                       if (y1 .lt. 1.e19) then
                         yz11=yz11_temp
                         yz12=yz12_temp
                         yz13=yz13_temp
                       else
                         yz11=0.0
                         yz12=wt1
                         yz13=wt2
                       endif
c
                       if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
                         yyy=wt1*y2+wt2*y3
                       else
                         yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1                      +wt2*(yz22*y2+yz23*y3+yz24*y4)
                       endif
                     endif

                     scr(ii,k) = yyy
!                    return
!                    end
                    enddo ! k
                  endif
               else 
                  scr(ii,:)=missingflag
               endif ! i is valid
            enddo ! iii

            xx=stax-fixm2
            if (xx .eq. 2.0) then
              staval(:)=scr(2,:)
            else

!             call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval) (now inlined)
!             subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)

              xxx=xx
              wt1=(xxx-3.)/(-1.)
              wt2=1.0-wt1
              yz22_temp=wt1*(xxx-4.)/(-2.)
              yz23_temp=wt2*(xxx-4.)/(-1.)
              yz24_temp=(xxx-2.)*(xxx-3.)/2.

              yz11_temp=(xxx-2.)*(xxx-3.)/(2.)
              yz12_temp=wt1*(xxx-1.)/(1.)
              yz13_temp=wt2*(xxx-1.)/(2.)

              do k=1,nz

               y1 = scr(1,k)
               y2 = scr(2,k)
               y3 = scr(3,k)
               y4 = scr(4,k)
c
               yyy=missingflag
               if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
                 wt1=(xxx-3.)/(-1.)
                 wt2=1.0-wt1
c
                 if (y4 .lt. 1.e19) then
                   yz22=yz22_temp
                   yz23=yz23_temp
                   yz24=yz24_temp
                 else
                   yz22=wt1
                   yz23=wt2
                   yz24=0.0
                 endif
c
                 if (y1 .lt. 1.e19) then
                   yz11=yz11_temp
                   yz12=yz12_temp
                   yz13=yz13_temp
                 else
                   yz11=0.0
                   yz12=wt1
                   yz13=wt2
                 endif
c
                 if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
                   yyy=wt1*y2+wt2*y3
                 else
                   yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1                +wt2*(yz22*y2+yz23*y3+yz24*y4)
                 endif
               endif

               staval(k) = yyy
!              return
!              end

               if(k .eq. 1)then
                 if(staval(1).eq.missingflag)then
                   iwrite = iwrite + 1
                   if(iwrite .le. 100)then
                     write(6,11)il,jl,i,j,stax,stay,ix1,ix2,iy1,iy2
     1                         ,fvi(il,jl,1),xx,yy,r(:,k),y1,y2,y3,y4      
     1                         ,xxx,yyy
11                   format(' hinterp_field_3d: staval = missingflag'  
     1                     ,4i5,2f9.2,4i5,13f9.2)
                   endif
                 endif
               endif

              enddo ! k

            endif ! xx = 2
c
            flaps(il,jl,:) = staval(:)
!           return
!           end

        else
          flaps(il,jl,:)=r_missing_data
        endif

      enddo ! il
      enddo ! jl
c
      return
      end


      subroutine hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,nz,
     .     grx,gry,fvi,flaps,wrapped)
c
!     implicit none
c
      include 'bgdata.inc'

      integer nx_bg,ny_bg,nx_laps,ny_laps,nz        

      dimension r(4),scr(4)

      logical wrapped
c
c *** input vertically interpolated field.
c *** output laps field
c
      real fvi(nx_bg,ny_bg,nz),
     .       flaps(nx_laps,ny_laps,nz),
     .       grx(nx_laps,ny_laps),gry(nx_laps,ny_laps)
c
      integer i,j,k
      integer istatus
      real    r_missing_data
c
      write(6,*)' subroutine hinterp_field_2d'

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'error getting r_missing_data - hinterp_field'
         return
      endif
c
c *** horizontally interpolate variable.
c
      ix = nx_bg
      iy = ny_bg

      do k=1,nz
        do jl=1,ny_laps
        do il=1,nx_laps
          if(grx(il,jl).lt.r_missing_data.and. 
     .       gry(il,jl).lt.r_missing_data)then
!           call gdtost_i(fvi(1,1,k),nx_bg,ny_bg,
!    .           grx(i,j),gry(i,j),flaps(i,j,k),wrapped)
!           subroutine gdtost_i(ab,ix,iy,stax,stay,staval,wrapped)

            stax = grx(il,jl)
            stay = gry(il,jl)
            
!           dimension ab(ix,iy),r(4),scr(4)
!           logical wrapped ! wni added
!           include 'bgdata.inc'
c_______________________________________________________________________________
c
            iy1=int(stay)-1
            if(stay.lt.1.0)iy1=1.0
            iy2=iy1+3
            ix1=int(stax)-1
            if(stax.lt.1.0) then 
              if (wrapped) then   ! wni
                ix1 = ix1 + ix
              else ! wni
                ix1=1.0
              endif ! wni
            endif ! wni
            ix2=ix1+3
            staval=missingflag
            fiym2=float(iy1)-1
            fixm2=float(ix1)-1
            ii=0
            do iii=ix1,ix2
               i=iii
c
c ****** account for wrapping around effect of global data at greenwich.
c
               if (wrapped) then    ! wni
                  if (i .lt. 1) i=i+ix  
                  if (i .gt. ix) i=i-ix
               endif

               ii=ii+1
               if (i .ge. 1 .and. i .le. ix) then 
                  jj=0
                  do jjj=iy1,iy2
                     j=jjj
c
c ************ account for n-s wrapping effect of global data.
c
                     if (wrapped) then  !wni
                        if (j .lt. 1) then
                           j=2-j
                           i=i-ix/2
                           if (i .lt. 1) i=i+ix
                           if (i .gt. ix) i=i-ix
                        endif
                        if (j .gt. iy) then
                           j=2*iy-j
                           i=i-ix/2
                           if (i .lt. 1) i=i+ix
                           if (i .gt. ix) i=i-ix
                        endif
                     endif

                     jj=jj+1
                     if (j .ge. 1 .and. j .le. iy) then
                        r(jj)=fvi(i,j,k) ! ab(i,j)
                     else
                        r(jj)=missingflag
                     endif
                  enddo
                  yy=stay-fiym2
                  if (yy .eq. 2.0) then
                     scr(ii)=r(2)
                  else
      
!                    call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii)) (now inlined)
!                    subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
                     y1 = r(1)
                     y2 = r(2)
                     y3 = r(3)
                     y4 = r(4)
                     xxx=yy
c
                     yyy=missingflag
                     if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
                       wt1=(xxx-3.)/(-1.)
                       wt2=1.0-wt1
c
                       if (y4 .lt. 1.e19) then
                         yz22=wt1*(xxx-4.)/(-2.)
                         yz23=wt2*(xxx-4.)/(-1.)
                         yz24=(xxx-2.)*(xxx-3.)/2.
                       else
                         yz22=wt1
                         yz23=wt2
                         yz24=0.0
                       endif
c
                       if (y1 .lt. 1.e19) then
                         yz11=(xxx-2.)*(xxx-3.)/(2.)
                         yz12=wt1*(xxx-1.)/(1.)
                         yz13=wt2*(xxx-1.)/(2.)
                       else
                         yz11=0.0
                         yz12=wt1
                         yz13=wt2
                       endif
c
                       if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
                         yyy=wt1*y2+wt2*y3
                       else
                         yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1                      +wt2*(yz22*y2+yz23*y3+yz24*y4)
                       endif
                     endif

                     scr(ii) = yyy
!                    return
!                    end

                  endif
               else 
                  scr(ii)=missingflag
               endif
            enddo
            xx=stax-fixm2
            if (xx .eq. 2.0) then
               staval=scr(2)
            else
!              call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval) (now inlined)
!              subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
               y1 = scr(1)
               y2 = scr(2)
               y3 = scr(3)
               y4 = scr(4)
               xxx=xx
c
               yyy=missingflag
               if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
                 wt1=(xxx-3.)/(-1.)
                 wt2=1.0-wt1
c
                 if (y4 .lt. 1.e19) then
                   yz22=wt1*(xxx-4.)/(-2.)
                   yz23=wt2*(xxx-4.)/(-1.)
                   yz24=(xxx-2.)*(xxx-3.)/(+2.)
                 else
                   yz22=wt1
                   yz23=wt2
                   yz24=0.0
                 endif
c
                 if (y1 .lt. 1.e19) then
                   yz11=(xxx-2.)*(xxx-3.)/(+2.)
                   yz12=wt1*(xxx-1.)/(1.)
                   yz13=wt2*(xxx-1.)/(2.)
                 else
                   yz11=0.0
                   yz12=wt1
                   yz13=wt2
                 endif
c
                 if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
                   yyy=wt1*y2+wt2*y3
                 else
                   yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1                +wt2*(yz22*y2+yz23*y3+yz24*y4)
                 endif
               endif

               staval = yyy
!              return
!              end

            endif
c
!           if(staval.eq.missingflag)then
!              print*,'hinterp_field_2d: staval = missingflag',il,jl
!           endif

            flaps(il,jl,k) = staval
!           return
!           end

          else
            flaps(il,jl,k)=r_missing_data
          endif
        enddo
        enddo
      enddo
c
      return
      end
c

c
c===============================================================================
c
      subroutine gdtost(ab,ix,iy,stax,stay,staval,wrapped)
c
c *** subroutine to return stations back-interpolated values(staval)
c        from uniform grid points using overlapping-quadratics.
c        gridded values of input array a dimensioned ab(ix,iy), where
c        ix = grid points in x, iy = grid points in y.  station
c        location given in terms of grid relative station x (stax)
c        and station column.
c *** values greater than 1.0e30 indicate missing data.
c
      dimension ab(ix,iy),r(4),scr(4)
      logical wrapped ! wni added
      include 'bgdata.inc'
c_______________________________________________________________________________
c
      iy1=int(stay)-1
      if(stay.lt.1.0)iy1=1.0
      iy2=iy1+3
      ix1=int(stax)-1
      if(stax.lt.1.0) then 
        if (wrapped) then   ! wni
          ix1 = ix1 + ix
        else ! wni
          ix1=1.0
        endif ! wni
      endif ! wni
      ix2=ix1+3
      staval=missingflag
      fiym2=float(iy1)-1
      fixm2=float(ix1)-1
      ii=0
      do iii=ix1,ix2
         i=iii
c
c ****** account for wrapping around effect of global data at greenwich.
c
         if (wrapped) then    ! wni
            if (i .lt. 1) i=i+ix  
            if (i .gt. ix) i=i-ix
         endif

         ii=ii+1
         if (i .ge. 1 .and. i .le. ix) then 
            jj=0
            do jjj=iy1,iy2
               j=jjj
c
c ************ account for n-s wrapping effect of global data.
c
               if (wrapped) then  !wni
                  if (j .lt. 1) then
                     j=2-j
                     i=i-ix/2
                     if (i .lt. 1) i=i+ix
                     if (i .gt. ix) i=i-ix
                  endif
                  if (j .gt. iy) then
                     j=2*iy-j
                     i=i-ix/2
                     if (i .lt. 1) i=i+ix
                     if (i .gt. ix) i=i-ix
                  endif
               endif

               jj=jj+1
               if (j .ge. 1 .and. j .le. iy) then
                  r(jj)=ab(i,j)
               else
                  r(jj)=missingflag
               endif
            enddo
            yy=stay-fiym2
            if (yy .eq. 2.0) then
               scr(ii)=r(2)
            else
               call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii))
            endif
         else 
            scr(ii)=missingflag
         endif
      enddo
      xx=stax-fixm2
      if (xx .eq. 2.0) then
         staval=scr(2)
      else
         call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval)
      endif
c
      if(staval.eq.missingflag)then
         print*,'gdtost: ',stax,stay,' staval = missingflag' 
      endif
      return
      end

c
c===============================================================================
c
      subroutine gdtost_i(ab,ix,iy,stax,stay,staval,wrapped)
c
c *** subroutine to return stations back-interpolated values(staval)
c        from uniform grid points using overlapping-quadratics.
c        gridded values of input array a dimensioned ab(ix,iy), where
c        ix = grid points in x, iy = grid points in y.  station
c        location given in terms of grid relative station x (stax)
c        and station column.
c *** values greater than 1.0e30 indicate missing data.
c
      dimension ab(ix,iy),r(4),scr(4)
      logical wrapped ! wni added
      include 'bgdata.inc'
c_______________________________________________________________________________
c
      iy1=int(stay)-1
      if(stay.lt.1.0)iy1=1.0
      iy2=iy1+3
      ix1=int(stax)-1
      if(stax.lt.1.0) then 
        if (wrapped) then   ! wni
          ix1 = ix1 + ix
        else ! wni
          ix1=1.0
        endif ! wni
      endif ! wni
      ix2=ix1+3
      staval=missingflag
      fiym2=float(iy1)-1
      fixm2=float(ix1)-1
      ii=0
      do iii=ix1,ix2
         i=iii
c
c ****** account for wrapping around effect of global data at greenwich.
c
         if (wrapped) then    ! wni
            if (i .lt. 1) i=i+ix  
            if (i .gt. ix) i=i-ix
         endif

         ii=ii+1
         if (i .ge. 1 .and. i .le. ix) then 
            jj=0
            do jjj=iy1,iy2
               j=jjj
c
c ************ account for n-s wrapping effect of global data.
c
               if (wrapped) then  !wni
                  if (j .lt. 1) then
                     j=2-j
                     i=i-ix/2
                     if (i .lt. 1) i=i+ix
                     if (i .gt. ix) i=i-ix
                  endif
                  if (j .gt. iy) then
                     j=2*iy-j
                     i=i-ix/2
                     if (i .lt. 1) i=i+ix
                     if (i .gt. ix) i=i-ix
                  endif
               endif

               jj=jj+1
               if (j .ge. 1 .and. j .le. iy) then
                  r(jj)=ab(i,j)
               else
                  r(jj)=missingflag
               endif
            enddo
            yy=stay-fiym2
            if (yy .eq. 2.0) then
               scr(ii)=r(2)
            else

!              call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii)) (now inlined)
!              subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
               y1 = r(1)
               y2 = r(2)
               y3 = r(3)
               y4 = r(4)
               xxx=yy
c
               yyy=missingflag
               if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
                 wt1=(xxx-3.)/(-1.)
                 wt2=1.0-wt1
c
                 if (y4 .lt. 1.e19) then
                   yz22=wt1*(xxx-4.)/(-2.)
                   yz23=wt2*(xxx-4.)/(-1.)
                   yz24=(xxx-2.)*(xxx-3.)/2.
                 else
                   yz22=wt1
                   yz23=wt2
                   yz24=0.0
                 endif
c
                 if (y1 .lt. 1.e19) then
                   yz11=(xxx-2.)*(xxx-3.)/(2.)
                   yz12=wt1*(xxx-1.)/(1.)
                   yz13=wt2*(xxx-1.)/(2.)
                 else
                   yz11=0.0
                   yz12=wt1
                   yz13=wt2
                 endif
c
                 if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
                   yyy=wt1*y2+wt2*y3
                 else
                   yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1                +wt2*(yz22*y2+yz23*y3+yz24*y4)
                 endif
               endif

               scr(ii) = yyy
!              return
!              end

            endif
         else 
            scr(ii)=missingflag
         endif
      enddo
      xx=stax-fixm2
      if (xx .eq. 2.0) then
         staval=scr(2)
      else
!        call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval) (now inlined)
!        subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
         y1 = scr(1)
         y2 = scr(2)
         y3 = scr(3)
         y4 = scr(4)
         xxx=xx
c
         yyy=missingflag
         if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
           wt1=(xxx-3.)/(-1.)
           wt2=1.0-wt1
c
           if (y4 .lt. 1.e19) then
             yz22=wt1*(xxx-4.)/(-2.)
             yz23=wt2*(xxx-4.)/(-1.)
             yz24=(xxx-2.)*(xxx-3.)/(+2.)
           else
             yz22=wt1
             yz23=wt2
             yz24=0.0
           endif
c
           if (y1 .lt. 1.e19) then
             yz11=(xxx-2.)*(xxx-3.)/(+2.)
             yz12=wt1*(xxx-1.)/(1.)
             yz13=wt2*(xxx-1.)/(2.)
           else
             yz11=0.0
             yz12=wt1
             yz13=wt2
           endif
c
           if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
             yyy=wt1*y2+wt2*y3
           else
             yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1          +wt2*(yz22*y2+yz23*y3+yz24*y4)
           endif
         endif

         staval = yyy
!        return
!        end

      endif
c
c      if(staval.eq.missingflag)then
c         print*,'gdtost_i: staval = missingflag'
c      endif
      return
      end

c
c===============================================================================
c
      subroutine gdtost2(a,ix,iy,stax,stay,staval)
*  subroutine to return stations back-interpolated values(staval)
*  from uniform grid points using overlapping-quadratics.
*  gridded values of input array a dimensioned a(ix,iy),where
*  ix=grid points in x, iy = grid points in y .  station
*  location given in terms of grid relative station x (stax)
*  and station column.
*  values greater than 1.0e30 indicate missing data.
*
      real a(ix,iy),r(4),scr(4),stax,stay,staval
     +  ,fixm2,fiym2,yy,xx
      iy1=int(stay)-1
      iy2=iy1+3
      ix1=int(stax)-1
      ix2=ix1+3
      staval=1e30
      fiym2=float(iy1)-1
      fixm2=float(ix1)-1
      ii=0
      do 100 i=ix1,ix2
      ii=ii+1
      if(i.ge.1.and.i.le.ix) go to 101
      scr(ii)=1e30
      go to 100
101   jj=0
      do 111 j=iy1,iy2
      jj=jj+1
      if(j.ge.1.and.j.le.iy) go to 112
      r(jj)=1e30
      go to 111
112   r(jj)=a(i,j)
111   continue
      yy=stay-fiym2
      call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii))
100   continue
      xx=stax-fixm2
      call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval)
      return
      end
c
cc ------------------------------------------------------------------
c

      subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
c
      include 'bgdata.inc'
      yyy=missingflag
      if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
         wt1=(xxx-x3)/(x2-x3)
         wt2=1.0-wt1
c
         if (y4 .lt. 1.e19) then
c           yz22=(xxx-x3)*(xxx-x4)/((x2-x3)*(x2-x4))
            yz22=wt1*(xxx-x4)/(x2-x4)
c           yz23=(xxx-x2)*(xxx-x4)/((x3-x2)*(x3-x4))
            yz23=wt2*(xxx-x4)/(x3-x4)
            yz24=(xxx-x2)*(xxx-x3)/((x4-x2)*(x4-x3))
         else
            yz22=wt1
            yz23=wt2
            yz24=0.0
         endif
c
         if (y1 .lt. 1.e19) then
            yz11=(xxx-x2)*(xxx-x3)/((x1-x2)*(x1-x3))
c           yz12=(xxx-x1)*(xxx-x3)/((x2-x1)*(x2-x3))
            yz12=wt1*(xxx-x1)/(x2-x1)
c           yz13=(xxx-x1)*(xxx-x2)/((x3-x1)*(x3-x2))
            yz13=wt2*(xxx-x1)/(x3-x1)
         else
            yz11=0.0
            yz12=wt1
            yz13=wt2
         endif
c
         if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
            yyy=wt1*y2+wt2*y3
         else
            yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1         +wt2*(yz22*y2+yz23*y3+yz24*y4)
         endif
c
      endif

      return
      end
