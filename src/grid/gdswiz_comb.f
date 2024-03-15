cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gdswiz03(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,
     &                    lrot,crot,srot)
c$$$  subprogram documentation block
c
c subprogram:  gdswiz03   gds wizard for lambert conformal conical
c   prgmmr: iredell       org: w/nmc23       date: 96-04-10
c
c abstract: this subprogram decodes the grib grid description section
c           (passed in integer form as decoded by subprogram w3fi63)
c           and returns one of the following:
c             (iopt=+1) earth coordinates of selected grid coordinates
c             (iopt=-1) grid coordinates of selected earth coordinates
c           for lambert conformal conical projections.
c           if the selected coordinates are more than one gridpoint
c           beyond the the edges of the grid domain, then the relevant
c           output elements are set to fill values.
c           the actual number of valid points computed is returned too.
c
c program history log:
c   96-04-10  iredell
c   96-10-01  iredell   protected against unresolvable points
c 1999-04-27  gilbert   corrected minor error calculating variable an
c                       for the secant projection case (rlati1.ne.rlati2).
c
c usage:    call gdswiz03(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,
c     &                   lrot,crot,srot)
c
c   input argument list:
c     kgds     - integer (200) gds parameters as decoded by w3fi63
c     iopt     - integer option flag
c                (+1 to compute earth coords of selected grid coords)
c                (-1 to compute grid coords of selected earth coords)
c     npts     - integer maximum number of coordinates
c     fill     - real fill value to set invalid output data
c                (must be impossible value; suggested value: -9999.)
c     xpts     - real (npts) grid x point coordinates if iopt>0
c     ypts     - real (npts) grid y point coordinates if iopt>0
c     rlon     - real (npts) earth longitudes in degrees e if iopt<0
c                (acceptable range: -360. to 360.)
c     rlat     - real (npts) earth latitudes in degrees n if iopt<0
c                (acceptable range: -90. to 90.)
c     lrot     - integer flag to return vector rotations if 1
c
c   output argument list:
c     xpts     - real (npts) grid x point coordinates if iopt<0
c     ypts     - real (npts) grid y point coordinates if iopt<0
c     rlon     - real (npts) earth longitudes in degrees e if iopt>0
c     rlat     - real (npts) earth latitudes in degrees n if iopt>0
c     nret     - integer number of valid points computed
c     crot     - real (npts) clockwise vector rotation cosines if lrot=1
c     srot     - real (npts) clockwise vector rotation sines if lrot=1
c                (ugrid=crot*uearth-srot*vearth;
c                 vgrid=srot*uearth+crot*vearth)
c
c attributes:
c   language: fortran 77
c
c$$$
      integer kgds(200)
      real xpts(npts),ypts(npts),rlon(npts),rlat(npts)
      real crot(npts),srot(npts)
      parameter(rerth=6.3712e6)
      parameter(pi=3.14159265358979,dpr=180./pi)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(kgds(1).eq.003) then
        im=kgds(2)
        jm=kgds(3)
        rlat1=kgds(4)*1.e-3
        rlon1=kgds(5)*1.e-3
        irot=mod(kgds(6)/8,2)
        orient=kgds(7)*1.e-3
        dx=kgds(8)
        dy=kgds(9)
        write(6,*) 'im,jm,rlat1,rlon1,dx,dy: ', im,jm,rlat1,rlon1,dx,dy
        iproj=mod(kgds(10)/128,2)
        iscan=mod(kgds(11)/128,2)
        jscan=mod(kgds(11)/64,2)
        nscan=mod(kgds(11)/32,2)
        rlati1=kgds(12)*1.e-3
        rlati2=kgds(13)*1.e-3
        h=(-1.)**iproj
        hi=(-1.)**iscan
        hj=(-1.)**(1-jscan)
        dxs=dx*hi
        dys=dy*hj
        if(rlati1.eq.rlati2) then
          an=sin(h*rlati1/dpr)
        else
          an=log(cos(rlati1/dpr)/cos(rlati2/dpr))/
     &       log(tan((h*90-rlati1)/2/dpr)/tan((h*90-rlati2)/2/dpr))
        endif
        de=rerth*cos(rlati1/dpr)*tan((h*rlati1+90)/2/dpr)**an/an
        if(h*rlat1.eq.90) then
          xp=1
          yp=1
        else
          dr=de/tan((h*rlat1+90)/2/dpr)**an
          dlon1=mod(rlon1-orient+180+3600,360.)-180
          xp=1-h*sin(an*dlon1/dpr)*dr/dxs
          yp=1+cos(an*dlon1/dpr)*dr/dys
        endif
        antr=1/(2*an)
        de2=de**2
        xmin=0
        xmax=im+1
        ymin=0
        ymax=jm+1
        nret=0
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  translate grid coordinates to earth coordinates
        if(iopt.eq.0.or.iopt.eq.1) then
          do n=1,npts
            if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and.
     &         ypts(n).ge.ymin.and.ypts(n).le.ymax) then
              di=(xpts(n)-xp)*dxs
              dj=(ypts(n)-yp)*dys
              dr2=di**2+dj**2
              if(dr2.lt.de2*1.e-6) then
                rlon(n)=0.
                rlat(n)=h*90.
              else
                rlon(n)=mod(orient+h/an*dpr*atan2(di,-dj)+3600,360.)
                rlat(n)=h*(2*dpr*atan((de2/dr2)**antr)-90)
              endif
              nret=nret+1
              if(lrot.eq.1) then
                if(irot.eq.1) then
                  dlon=mod(rlon(n)-orient+180+3600,360.)-180
                  crot(n)=h*cos(an*dlon/dpr)
                  srot(n)=sin(an*dlon/dpr)
                else
                  crot(n)=1
                  srot(n)=0
                endif
              endif
            else
              rlon(n)=fill
              rlat(n)=fill
            endif
          enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  translate earth coordinates to grid coordinates
        elseif(iopt.eq.-1) then
          do n=1,npts
            if(abs(rlon(n)).le.360.and.abs(rlat(n)).le.90.and.
     &                                 h*rlat(n).ne.-90) then
              dr=de*tan((90-h*rlat(n))/2/dpr)**an
              dlon=mod(rlon(n)-orient+180+3600,360.)-180
              xpts(n)=xp+h*sin(an*dlon/dpr)*dr/dxs
              ypts(n)=yp-cos(an*dlon/dpr)*dr/dys
              if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and.
     &           ypts(n).ge.ymin.and.ypts(n).le.ymax) then
                nret=nret+1
                if(lrot.eq.1) then
                  if(irot.eq.1) then
                    crot(n)=h*cos(an*dlon/dpr)
                    srot(n)=sin(an*dlon/dpr)
                  else
                    crot(n)=1
                    srot(n)=0
                  endif
                endif
              else
                xpts(n)=fill
                ypts(n)=fill
              endif
            else
              xpts(n)=fill
              ypts(n)=fill
            endif
          enddo
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  projection unrecognized
      else
        iret=-1
        if(iopt.ge.0) then
          do n=1,npts
            rlon(n)=fill
            rlat(n)=fill
          enddo
        endif
        if(iopt.le.0) then
          do n=1,npts
            xpts(n)=fill
            ypts(n)=fill
          enddo
        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gdswiz05(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,
     &                    lrot,crot,srot)
c$$$  subprogram documentation block
c
c subprogram:  gdswiz05   gds wizard for polar stereographic azimuthal
c   prgmmr: iredell       org: w/nmc23       date: 96-04-10
c
c abstract: this subprogram decodes the grib grid description section
c           (passed in integer form as decoded by subprogram w3fi63)
c           and returns one of the following:
c             (iopt=+1) earth coordinates of selected grid coordinates
c             (iopt=-1) grid coordinates of selected earth coordinates
c           for polar stereographic azimuthal projections.
c           if the selected coordinates are more than one gridpoint
c           beyond the the edges of the grid domain, then the relevant
c           output elements are set to fill values.
c           the actual number of valid points computed is returned too.
c
c program history log:
c   96-04-10  iredell
c
c usage:    call gdswiz05(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,
c     &                   lrot,crot,srot)
c
c   input argument list:
c     kgds     - integer (200) gds parameters as decoded by w3fi63
c     iopt     - integer option flag
c                (+1 to compute earth coords of selected grid coords)
c                (-1 to compute grid coords of selected earth coords)
c     npts     - integer maximum number of coordinates
c     fill     - real fill value to set invalid output data
c                (must be impossible value; suggested value: -9999.)
c     xpts     - real (npts) grid x point coordinates if iopt>0
c     ypts     - real (npts) grid y point coordinates if iopt>0
c     rlon     - real (npts) earth longitudes in degrees e if iopt<0
c                (acceptable range: -360. to 360.)
c     rlat     - real (npts) earth latitudes in degrees n if iopt<0
c                (acceptable range: -90. to 90.)
c     lrot     - integer flag to return vector rotations if 1
c
c   output argument list:
c     xpts     - real (npts) grid x point coordinates if iopt<0
c     ypts     - real (npts) grid y point coordinates if iopt<0
c     rlon     - real (npts) earth longitudes in degrees e if iopt>0
c     rlat     - real (npts) earth latitudes in degrees n if iopt>0
c     nret     - integer number of valid points computed
c     crot     - real (npts) clockwise vector rotation cosines if lrot=1
c     srot     - real (npts) clockwise vector rotation sines if lrot=1
c                (ugrid=crot*uearth-srot*vearth;
c                 vgrid=srot*uearth+crot*vearth)
c attributes:
c   language: fortran 77
c
c$$$
      integer kgds(200)
      real xpts(npts),ypts(npts),rlon(npts),rlat(npts)
      real crot(npts),srot(npts)
      parameter(rerth=6.3712e6)
      parameter(pi=3.14159265358979,dpr=180./pi)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(kgds(1).eq.005) then
        im=kgds(2)
        jm=kgds(3)
        rlat1=kgds(4)*1.e-3
        rlon1=kgds(5)*1.e-3
        irot=mod(kgds(6)/8,2)
        orient=kgds(7)*1.e-3
        dx=kgds(8)
        dy=kgds(9)
        iproj=mod(kgds(10)/128,2)
        iscan=mod(kgds(11)/128,2)
        jscan=mod(kgds(11)/64,2)
        nscan=mod(kgds(11)/32,2)
        h=(-1.)**iproj
        hi=(-1.)**iscan
        hj=(-1.)**(1-jscan)
        dxs=dx*hi
        dys=dy*hj
        de=(1.+sin(60./dpr))*rerth
        dr=de*cos(rlat1/dpr)/(1+h*sin(rlat1/dpr))
        xp=1-h*sin((rlon1-orient)/dpr)*dr/dxs
        yp=1+cos((rlon1-orient)/dpr)*dr/dys
        de2=de**2
        xmin=0
        xmax=im+1
        ymin=0
        ymax=jm+1
        nret=0
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  translate grid coordinates to earth coordinates
        if(iopt.eq.0.or.iopt.eq.1) then
          do n=1,npts
            if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and.
     &         ypts(n).ge.ymin.and.ypts(n).le.ymax) then
              di=(xpts(n)-xp)*dxs
              dj=(ypts(n)-yp)*dys
              dr2=di**2+dj**2
              if(dr2.lt.de2*1.e-6) then
                rlon(n)=0.
                rlat(n)=h*90.
              else
                rlon(n)=mod(orient+h*dpr*atan2(di,-dj)+3600,360.)
                rlat(n)=h*dpr*asin((de2-dr2)/(de2+dr2))
              endif
              nret=nret+1
              if(lrot.eq.1) then
                if(irot.eq.1) then
                  crot(n)=h*cos((rlon(n)-orient)/dpr)
                  srot(n)=sin((rlon(n)-orient)/dpr)
                else
                  crot(n)=1
                  srot(n)=0
                endif
              endif
            else
              rlon(n)=fill
              rlat(n)=fill
            endif
          enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  translate earth coordinates to grid coordinates
        elseif(iopt.eq.-1) then
          do n=1,npts
            if(abs(rlon(n)).le.360.and.abs(rlat(n)).le.90.and.
     &                                 h*rlat(n).ne.-90) then
              dr=de*tan((90-h*rlat(n))/2/dpr)
              xpts(n)=xp+h*sin((rlon(n)-orient)/dpr)*dr/dxs
              ypts(n)=yp-cos((rlon(n)-orient)/dpr)*dr/dys
              if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and.
     &           ypts(n).ge.ymin.and.ypts(n).le.ymax) then
                nret=nret+1
                if(lrot.eq.1) then
                  if(irot.eq.1) then
                    crot(n)=h*cos((rlon(n)-orient)/dpr)
                    srot(n)=sin((rlon(n)-orient)/dpr)
                  else
                    crot(n)=1
                    srot(n)=0
                  endif
                endif
              else
                xpts(n)=fill
                ypts(n)=fill
              endif
            else
              xpts(n)=fill
              ypts(n)=fill
            endif
          enddo
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  projection unrecognized
      else
        iret=-1
        if(iopt.ge.0) then
          do n=1,npts
            rlon(n)=fill
            rlat(n)=fill
          enddo
        endif
        if(iopt.le.0) then
          do n=1,npts
            xpts(n)=fill
            ypts(n)=fill
          enddo
        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gdswiz(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,
     &                  lrot,crot,srot)
c$$$  subprogram documentation block
c
c subprogram:  gdswiz     grid description section wizard
c   prgmmr: iredell       org: w/nmc23       date: 96-04-10
c
c abstract: this subprogram decodes the grib grid description section
c           (passed in integer form as decoded by subprogram w3fi63)
c           and returns one of the following:
c             (iopt= 0) grid and earth coordinates of all grid points
c             (iopt=+1) earth coordinates of selected grid coordinates
c             (iopt=-1) grid coordinates of selected earth coordinates
c           the current code recognizes the following projections:
c             (kgds(1)=000) equidistant cylindrical
c             (kgds(1)=001) mercator cylindrical
c             (kgds(1)=003) lambert conformal conical
c             (kgds(1)=004) gaussian cylindrical
c             (kgds(1)=005) polar stereographic azimuthal
c             (kgds(1)=201) staggered rotated equidistant cylindrical
c             (kgds(1)=202) rotated equidistant cylindrical
c             (kgds(1)=203) staggered rotated equidistant cylindrical 2-d
c           if the selected coordinates are more than one gridpoint
c           beyond the the edges of the grid domain, then the relevant
c           output elements are set to fill values.  also if iopt=0,
c           if the number of grid points exceeds the number allotted,
c           then all the output elements are set to fill values.
c           the actual number of valid points computed is returned too.
c
c program history log:
c   96-04-10  iredell
c   98-08-20  baldwin  add type 203 staggered 2-d eta grids
c
c usage:    call gdswiz(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,
c     &                 lrot,crot,srot)
c
c   input argument list:
c     kgds     - integer (200) gds parameters as decoded by w3fi63
c     iopt     - integer option flag
c                ( 0 to compute earth coords of all the grid points)
c                (+1 to compute earth coords of selected grid coords)
c                (-1 to compute grid coords of selected earth coords)
c     npts     - integer maximum number of coordinates
c     fill     - real fill value to set invalid output data
c                (must be impossible value; suggested value: -9999.)
c     xpts     - real (npts) grid x point coordinates if iopt>0
c     ypts     - real (npts) grid y point coordinates if iopt>0
c     rlon     - real (npts) earth longitudes in degrees e if iopt<0
c                (acceptable range: -360. to 360.)
c     rlat     - real (npts) earth latitudes in degrees n if iopt<0
c                (acceptable range: -90. to 90.)
c     lrot     - integer flag to return vector rotations if 1
c
c   output argument list:
c     xpts     - real (npts) grid x point coordinates if iopt<=0
c     ypts     - real (npts) grid y point coordinates if iopt<=0
c     rlon     - real (npts) earth longitudes in degrees e if iopt>=0
c     rlat     - real (npts) earth latitudes in degrees n if iopt>=0
c     nret     - integer number of valid points computed
c                (-1 if projection unrecognized)
c     crot     - real (npts) clockwise vector rotation cosines if lrot=1
c     srot     - real (npts) clockwise vector rotation sines if lrot=1
c                (ugrid=crot*uearth-srot*vearth;
c                 vgrid=srot*uearth+crot*vearth)
c
c subprograms called:
c   gdswiz00     gds wizard for equidistant cylindrical
c   gdswiz01     gds wizard for mercator cylindrical
c   gdswiz03     gds wizard for lambert conformal conical
c   gdswiz04     gds wizard for gaussian cylindrical
c   gdswiz05     gds wizard for polar stereographic azimuthal
c   gdswizc9     gds wizard for rotated equidistant cylindrical
c   gdswizca     gds wizard for rotated equidistant cylindrical
c   gdswizcb     gds wizard for rotated equidistant cylindrical 2-d
c
c attributes:
c   language: fortran 77
c
c$$$
      integer kgds(200)
      real xpts(npts),ypts(npts),rlon(npts),rlat(npts)
      real crot(npts),srot(npts)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  compute grid coordinates for all grid points
      if(iopt.eq.0) then
        if(kgds(1).eq.201) then
          im=kgds(7)*2-1
          jm=kgds(8)
          kscan=mod(kgds(11)/256,2)
          if(kscan.eq.0) then
            is1=(jm+1)/2
            nm=(im/2+1)*jm-jm/2
          else
            is1=jm/2
            nm=im/2*jm+jm/2
          endif
        elseif(kgds(1).eq.202) then
          im=kgds(7)
          jm=kgds(8)
          nm=im*jm
        elseif(kgds(1).eq.203) then
          im=kgds(2)
          jm=kgds(3)
          nm=im*jm
          kscan=mod(kgds(11)/256,2)
          if(kscan.eq.0) then
            is1=(jm+1)/2
          else
            is1=jm/2
          endif
        else
          im=kgds(2)
          jm=kgds(3)
          nm=im*jm
        endif
        nscan=mod(kgds(11)/32,2)
        if(nm.le.npts) then
          if(kgds(1).eq.201) then
            do n=1,nm
              nn=2*n-1+kscan
              if(nscan.eq.0) then
                j=(nn-1)/im+1
                i=nn-im*(j-1)
              else
                i=(nn-1)/jm+1
                j=nn-jm*(i-1)
              endif
              xpts(n)=is1+(i-(j-kscan))/2
              ypts(n)=(i+(j-kscan))/2
            enddo
          elseif(kgds(1).eq.203) then
            do n=1,nm
              if(nscan.eq.0) then
                j=(n-1)/im+1
                i=(n-im*(j-1))*2-mod(j+kscan,2)
              else
                i=(n-1)/jm+1
                j=(n-jm*(i-1))*2-mod(i+kscan,2)
              endif
              xpts(n)=is1+(i-(j-kscan))/2
              ypts(n)=(i+(j-kscan))/2
            enddo
          else
            do n=1,nm
              if(nscan.eq.0) then
                j=(n-1)/im+1
                i=n-im*(j-1)
              else
                i=(n-1)/jm+1
                j=n-jm*(i-1)
              endif
              xpts(n)=i
              ypts(n)=j
            enddo
          endif
          do n=nm+1,npts
            xpts(n)=fill
            ypts(n)=fill
          enddo
        else
          do n=1,npts
            xpts(n)=fill
            ypts(n)=fill
          enddo
        endif
        iopf=1
      else
        iopf=iopt
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  equidistant cylindrical
      if(kgds(1).eq.000) then
        call gdswiz00(kgds,iopf,npts,fill,xpts,ypts,rlon,rlat,nret,
     &                lrot,crot,srot)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  mercator cylindrical
      elseif(kgds(1).eq.001) then
        call gdswiz01(kgds,iopf,npts,fill,xpts,ypts,rlon,rlat,nret,
     &                lrot,crot,srot)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  lambert conformal conical
      elseif(kgds(1).eq.003) then
        call gdswiz03(kgds,iopf,npts,fill,xpts,ypts,rlon,rlat,nret,
     &                lrot,crot,srot)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  gaussian cylindrical
      elseif(kgds(1).eq.004) then
!        call gdswiz04(kgds,iopf,npts,fill,xpts,ypts,rlon,rlat,nret,
!     &                lrot,crot,srot)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  polar stereographic azimuthal
      elseif(kgds(1).eq.005) then
        call gdswiz05(kgds,iopf,npts,fill,xpts,ypts,rlon,rlat,nret,
     &                lrot,crot,srot)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  staggered rotated equidistant cylindrical
      elseif(kgds(1).eq.201) then
!        call gdswizc9(kgds,iopf,npts,fill,xpts,ypts,rlon,rlat,nret,
!     &                lrot,crot,srot)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  rotated equidistant cylindrical
      elseif(kgds(1).eq.202) then
!        call gdswizca(kgds,iopf,npts,fill,xpts,ypts,rlon,rlat,nret,
!     &                lrot,crot,srot)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  staggered rotated equidistant cylindrical
      elseif(kgds(1).eq.203) then
!        call gdswizcb(kgds,iopf,npts,fill,xpts,ypts,rlon,rlat,nret,
!     &                lrot,crot,srot)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  projection unrecognized
      else
        iret=-1
        if(iopt.ge.0) then
          do n=1,npts
            rlon(n)=fill
            rlat(n)=fill
          enddo
        endif
        if(iopt.le.0) then
          do n=1,npts
            xpts(n)=fill
            ypts(n)=fill
          enddo
        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gdswiz01(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,
     &                    lrot,crot,srot)
c$$$  subprogram documentation block
c
c subprogram:  gdswiz01   gds wizard for mercator cylindrical
c   prgmmr: iredell       org: w/nmc23       date: 96-04-10
c
c abstract: this subprogram decodes the grib grid description section
c           (passed in integer form as decoded by subprogram w3fi63)
c           and returns one of the following:
c             (iopt=+1) earth coordinates of selected grid coordinates
c             (iopt=-1) grid coordinates of selected earth coordinates
c           for mercator cylindrical projections.
c           if the selected coordinates are more than one gridpoint
c           beyond the the edges of the grid domain, then the relevant
c           output elements are set to fill values.
c           the actual number of valid points computed is returned too.
c
c program history log:
c   96-04-10  iredell
c   96-10-01  iredell   protected against unresolvable points
c
c usage:    call gdswiz01(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,
c     &                   lrot,crot,srot)
c
c   input argument list:
c     kgds     - integer (200) gds parameters as decoded by w3fi63
c     iopt     - integer option flag
c                (+1 to compute earth coords of selected grid coords)
c                (-1 to compute grid coords of selected earth coords)
c     npts     - integer maximum number of coordinates
c     fill     - real fill value to set invalid output data
c                (must be impossible value; suggested value: -9999.)
c     xpts     - real (npts) grid x point coordinates if iopt>0
c     ypts     - real (npts) grid y point coordinates if iopt>0
c     rlon     - real (npts) earth longitudes in degrees e if iopt<0
c                (acceptable range: -360. to 360.)
c     rlat     - real (npts) earth latitudes in degrees n if iopt<0
c                (acceptable range: -90. to 90.)
c     lrot     - integer flag to return vector rotations if 1
c
c   output argument list:
c     xpts     - real (npts) grid x point coordinates if iopt<0
c     ypts     - real (npts) grid y point coordinates if iopt<0
c     rlon     - real (npts) earth longitudes in degrees e if iopt>0
c     rlat     - real (npts) earth latitudes in degrees n if iopt>0
c     nret     - integer number of valid points computed
c     crot     - real (npts) clockwise vector rotation cosines if lrot=1
c     srot     - real (npts) clockwise vector rotation sines if lrot=1
c                (ugrid=crot*uearth-srot*vearth;
c                 vgrid=srot*uearth+crot*vearth)

c
c attributes:
c   language: fortran 77
c
c$$$
      integer kgds(200)
      real xpts(npts),ypts(npts),rlon(npts),rlat(npts)
      real crot(npts),srot(npts)
      parameter(rerth=6.3712e6)
      parameter(pi=3.14159265358979,dpr=180./pi)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(kgds(1).eq.001) then
        im=kgds(2)
        jm=kgds(3)
        rlat1=kgds(4)*1.e-3
        rlon1=kgds(5)*1.e-3
        rlat2=kgds(7)*1.e-3
        rlon2=kgds(8)*1.e-3
        rlati=kgds(9)*1.e-3
        iscan=mod(kgds(11)/128,2)
        jscan=mod(kgds(11)/64,2)
        nscan=mod(kgds(11)/32,2)
        dx=kgds(12)
        dy=kgds(13)
        hi=(-1.)**iscan
        hj=(-1.)**(1-jscan)
        dlon=hi*(mod(hi*(rlon2-rlon1)-1+3600,360.)+1)/(im-1)
        dlat=hj*dy/(rerth*cos(rlati/dpr))
        ye=1-log(tan((rlat1+90)/2/dpr))/dlat
        xmin=0
        xmax=im+1
        if(im.eq.nint(360/abs(dlon))) xmax=im+2
        ymin=0
        ymax=jm+1
        nret=0
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  translate grid coordinates to earth coordinates
        if(iopt.eq.0.or.iopt.eq.1) then
          do n=1,npts
            if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and.
     &         ypts(n).ge.ymin.and.ypts(n).le.ymax) then
              rlon(n)=mod(rlon1+dlon*(xpts(n)-1)+3600,360.)
              rlat(n)=2*atan(exp(dlat*(ypts(n)-ye)))*dpr-90
              nret=nret+1
              if(lrot.eq.1) then
                crot(n)=1
                srot(n)=0
              endif
            else
              rlon(n)=fill
              rlat(n)=fill
            endif
          enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  translate earth coordinates to grid coordinates
        elseif(iopt.eq.-1) then
          do n=1,npts
            if(abs(rlon(n)).le.360.and.abs(rlat(n)).lt.90) then
              xpts(n)=1+hi*mod(hi*(rlon(n)-rlon1)+3600,360.)/dlon
              ypts(n)=ye+log(tan((rlat(n)+90)/2/dpr))/dlat
              if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and.
     &           ypts(n).ge.ymin.and.ypts(n).le.ymax) then
                nret=nret+1
                if(lrot.eq.1) then
                  crot(n)=1
                  srot(n)=0
                endif
              else
                xpts(n)=fill
                ypts(n)=fill
              endif
            else
              xpts(n)=fill
              ypts(n)=fill
            endif
          enddo
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  projection unrecognized
      else
        iret=-1
        if(iopt.ge.0) then
          do n=1,npts
            rlon(n)=fill
            rlat(n)=fill
          enddo
        endif
        if(iopt.le.0) then
          do n=1,npts
            xpts(n)=fill
            ypts(n)=fill
          enddo
        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end

ccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gdswiz00(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,
     &                    lrot,crot,srot)
c$$$  subprogram documentation block
c
c subprogram:  gdswiz00   gds wizard for equidistant cylindrical
c   prgmmr: iredell       org: w/nmc23       date: 96-04-10
c
c abstract: this subprogram decodes the grib grid description section
c           (passed in integer form as decoded by subprogram w3fi63)
c           and returns one of the following:
c             (iopt=+1) earth coordinates of selected grid coordinates
c             (iopt=-1) grid coordinates of selected earth coordinates
c           for equidistant cylindrical projections.
c           if the selected coordinates are more than one gridpoint
c           beyond the the edges of the grid domain, then the relevant
c           output elements are set to fill values.
c           the actual number of valid points computed is returned too.
c
c program history log:
c   96-04-10  iredell
c
c usage:    call gdswiz00(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,
c     &                   lrot,crot,srot)
c
c   input argument list:
c     kgds     - integer (200) gds parameters as decoded by w3fi63
c     iopt     - integer option flag
c                (+1 to compute earth coords of selected grid coords)
c                (-1 to compute grid coords of selected earth coords)
c     npts     - integer maximum number of coordinates
c     fill     - real fill value to set invalid output data
c                (must be impossible value; suggested value: -9999.)
c     xpts     - real (npts) grid x point coordinates if iopt>0
c     ypts     - real (npts) grid y point coordinates if iopt>0
c     rlon     - real (npts) earth longitudes in degrees e if iopt<0
c                (acceptable range: -360. to 360.)
c     rlat     - real (npts) earth latitudes in degrees n if iopt<0
c                (acceptable range: -90. to 90.)
c     lrot     - integer flag to return vector rotations if 1
c
c   output argument list:
c     xpts     - real (npts) grid x point coordinates if iopt<0
c     ypts     - real (npts) grid y point coordinates if iopt<0
c     rlon     - real (npts) earth longitudes in degrees e if iopt>0
c     rlat     - real (npts) earth latitudes in degrees n if iopt>0
c     nret     - integer number of valid points computed
c     crot     - real (npts) clockwise vector rotation cosines if lrot=1
c     srot     - real (npts) clockwise vector rotation sines if lrot=1
c                (ugrid=crot*uearth-srot*vearth;
c                 vgrid=srot*uearth+crot*vearth)
c attributes:
c   language: fortran 77
c
c$$$
      integer kgds(200)
      real xpts(npts),ypts(npts),rlon(npts),rlat(npts)
      real crot(npts),srot(npts)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(kgds(1).eq.000) then
        im=kgds(2)
        jm=kgds(3)
        rlat1=kgds(4)*1.e-3
        rlon1=kgds(5)*1.e-3
        rlat2=kgds(7)*1.e-3
        rlon2=kgds(8)*1.e-3
        iscan=mod(kgds(11)/128,2)
        jscan=mod(kgds(11)/64,2)
        nscan=mod(kgds(11)/32,2)
        hi=(-1.)**iscan
        hj=(-1.)**(1-jscan)
        dlon=hi*(mod(hi*(rlon2-rlon1)-1+3600,360.)+1)/(im-1)
        dlat=(rlat2-rlat1)/(jm-1)
        xmin=0
        xmax=im+1
        if(im.eq.nint(360/abs(dlon))) xmax=im+2
        ymin=0
        ymax=jm+1
        nret=0
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  translate grid coordinates to earth coordinates
        if(iopt.eq.0.or.iopt.eq.1) then
          do n=1,npts
            if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and.
     &         ypts(n).ge.ymin.and.ypts(n).le.ymax) then
              rlon(n)=mod(rlon1+dlon*(xpts(n)-1)+3600,360.)
              rlat(n)=min(max(rlat1+dlat*(ypts(n)-1),-90.),90.)
              nret=nret+1
              if(lrot.eq.1) then
                crot(n)=1
                srot(n)=0
              endif
            else
              rlon(n)=fill
              rlat(n)=fill
            endif
          enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  translate earth coordinates to grid coordinates
        elseif(iopt.eq.-1) then
          do n=1,npts
            if(abs(rlon(n)).le.360.and.abs(rlat(n)).le.90) then
              xpts(n)=1+hi*mod(hi*(rlon(n)-rlon1)+3600,360.)/dlon
              ypts(n)=1+(rlat(n)-rlat1)/dlat
              if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and.
     &           ypts(n).ge.ymin.and.ypts(n).le.ymax) then
                nret=nret+1
                if(lrot.eq.1) then
                  crot(n)=1
                  srot(n)=0
                endif
              else
                xpts(n)=fill
                ypts(n)=fill
              endif
            else
              xpts(n)=fill
              ypts(n)=fill
            endif
          enddo
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  projection unrecognized
      else
        iret=-1
        if(iopt.ge.0) then
          do n=1,npts
            rlon(n)=fill
            rlat(n)=fill
          enddo
        endif
        if(iopt.le.0) then
          do n=1,npts
            xpts(n)=fill
            ypts(n)=fill
          enddo
        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
