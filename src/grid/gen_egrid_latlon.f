      subroutine etall(im,jm,tph0d,tlm0d,dlmd,dphd,hlat,hlon,vlat,vlon)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram: etall         compute earth latitude & lonigtude of
c                           eta grid points
c   prgmmr: rogers          org: w/np22     date: 90-06-13
c
c abstract: computes the earth latitude and longitude of eta grid
c   points (both h and v points)
c
c program history log:
c   90-06-13  e.rogers
c   98-06-09  m.baldwin - convert to 2-d code
c   01-01-03  t black   - modified for mpi
c
c usage:    call etall(hlat,hlon,vlat,vlon)
c   input argument list:
c     none
c
c   output argument list:
c     hlat     - latitude of h grid points in radians (neg=s)
c     hlon     - longitude of h grid points in radians (e)
c     vlat     - latitude of v grid points in radians (neg=s)
c     vlon     - longitude of v grid points in radians (e)
c
c
c attributes:
c   language: fortran 90
c   machine:  ibm rs/6000 sp
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
                             d i m e n s i o n
     & idat  (3)
c
     &,glath (im,jm),glonh (im,jm),glatv (im,jm),glonv (im,jm)
c
        real hlat(im,jm),hlon(im,jm),vlat(im,jm),vlon(im,jm)
c-----------------------------------------------------------------------
c--------------logical file names---------------------------------------
c-----------------------------------------------------------------------
                             d a t a
     &  list  /06/
c-----------------------------------------------------------------------
c--------------universal constants--------------------------------------
c-----------------------------------------------------------------------
                             d a t a
     & pi/3.141592654/
c-----------------------------------------------------------------------
 1000 format(' i j=',2i4,' glat=',e12.5,' glon=',e12.5,' elat=',e12.5
     &,' elon=',e12.5,' wlon = ',e12.5)
c----------------------------------------------------------------------
c--------------derived geometrical constants----------------------------
c----------------------------------------------------------------------
cc
	wbd=-(im-1)*dlmd
	sbd=-(jm-1)/2*dphd
      dtr = pi / 180.0
      tph0 = tph0d * dtr
      wb = wbd * dtr
      sb = sbd * dtr
      dlm = dlmd * dtr
      dph = dphd * dtr

	write(6,*) 'tph0,wb,sb,dlm,dph,dtr: ', tph0,wb,sb,dlm,dph,dtr

      tdlm = dlm + dlm
      tdph = dph + dph
c
      stph0 = sin(tph0)
      ctph0 = cos(tph0)
c
c-----------------------------------------------------------------------
c---compute geographic lat and long of eta grid points (h & v points)---
c-----------------------------------------------------------------------
      do 200 j = 1,jm
c
         tlmh = wb - tdlm + mod(j+1,2) * dlm
         tphh = sb+(j-1)*dph
         tlmv = wb - tdlm + mod(j,2) * dlm
         tphv = tphh
         stph = sin(tphh)
         ctph = cos(tphh)
         stpv = sin(tphv)
         ctpv = cos(tphv)
c----------------------------------------------------------------------
c---------- compute earth latitude/longitude of h points --------------
c----------------------------------------------------------------------
         do 201 i = 1,im
           tlmh = tlmh + tdlm
           sphh = ctph0 * stph + stph0 * ctph * cos(tlmh)
           glath(i,j) = asin(sphh)
           clmh = ctph * cos(tlmh) / (cos(glath(i,j)) * ctph0)
     1               - tan(glath(i,j)) * tan(tph0)
           if(clmh .gt. 1.) clmh = 1.0
           if(clmh .lt. -1.) clmh = -1.0
           facth = 1.
           if(tlmh .gt. 0.) facth = -1.
c          write(6,88888) i,j, clmh
c8888      format(2x,2i6,1x,e12.5)
           glonh(i,j) = -tlm0d * dtr + facth * acos(clmh)
c          if(i .eq. 1) then
c           write(list,99995) i,j,glath(i,j),glonh(i,j)
c9995       format(2x,2(i6,1x),2(e12.5,1x))
c          end if
c
c
           hlat(i,j) = glath(i,j) / dtr
c           hlon(i,j) = 360.0 - glonh(i,j) / dtr
	    hlon(i,j)= -glonh(i,j)/dtr
	   
           if(hlon(i,j) .gt. 180.) hlon(i,j) = hlon(i,j) - 360.
           if(hlon(i,j) .lt. -180.) hlon(i,j) = hlon(i,j) + 360.
  201    continue
c----------------------------------------------------------------------
c---------- compute earth latitude/longitude of v points --------------
c----------------------------------------------------------------------
         do 202 i = 1,im
           tlmv = tlmv + tdlm
           sphv = ctph0 * stpv + stph0 * ctpv * cos(tlmv)
           glatv(i,j) = asin(sphv)
           clmv = ctpv * cos(tlmv) / (cos(glatv(i,j)) * ctph0)
     1          - tan(glatv(i,j)) * tan(tph0)
           if(clmv .gt. 1.) clmv = 1.
           if(clmv .lt. -1.) clmv = -1.
           factv = 1.
           if(tlmv .gt. 0.) factv = -1.
           glonv(i,j) = -tlm0d * dtr + factv * acos(clmv)
c          if(i.eq.1) then
c           write(list,99995) i,j,glatv(i,j),glonv(i,j)
c          end if
c
c    convert into degrees and east longitude
c
           vlat(i,j) = glatv(i,j) / dtr
c           vlon(i,j) = 360.0 - glonv(i,j) / dtr
           vlon(i,j) = -glonv(i,j) / dtr
           if(vlon(i,j) .gt. 180.) vlon(i,j) = vlon(i,j) - 360.
           if(vlon(i,j) .lt. -180.) vlon(i,j) = vlon(i,j) + 360.
  202    continue
  200 continue
c
c     do 210 j = 1,jm
c         write(list,88888) j,hlat(1,j),hlon(1,j),vlat(1,j),vlon(1,j)
c8888    format(2x,i5,1x,4(e12.5,1x))
c      end if
c 210 continue
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       subroutine vecrot_rotlat(im,jm,tph0d,tlm0d,vlat,vlon,dray,cray)
                                                                                         
c  subprogram: rotlle        rotate winds on lat/long grid to e-grid
c   prgmmr: t.black         org: w/np22     date: ??-??-??
c
c abstract: rotates winds on the lat/long grid to the eta model grid
c               (same angles can be used to do reverse function)
c
c program history log:
c   ??-??-??  t.black
c   98-06-08  m.baldwin - convert to 2-d code
c   01-01-03  t black - modified for mpi
c
c   input argument list:
c     vlat     - latitude of e-grid v points (degrees)
c     vlon     - longitude of e-grid v points (degrees)
c
c   output argument list:
c     dray     - rotation cosine
c     cray     - rotation sine
c
c attributes:
c   language: fortran 90
c   machine:  ibm rs/6000 sp
c
c***
c*** rotate the lat-lon winds to/from the e-grid
c***
c
c    n o t e : input lat/long must be converted to radians !!!
c
c----------------------------------------------------------------------
        parameter(d2rad=1.745329e-2)
                                                                                         
                         d i m e n s i o n
     1 vlat(im,jm),vlon(im,jm)
     2,cray(im,jm),dray(im,jm)
c----------------------------------------------------------------------
                                                                                         
        erphi0=tph0d*d2rad

        if (tlm0d .lt. 0) then
        erlam0=(tlm0d+360.)*d2rad
        else
        erlam0=tlm0d*d2rad
        endif
                                                                                         
      sphi0 = sin(erphi0)
      cphi0 = cos(erphi0)
c
	
	write(6,*) 'tph0d, tlm0d: ', tph0d, tlm0d
	write(6,*) 'im,jm: ', im,jm

      do j = 1, jm
      do i = 1, im

!mp dangerous test (which had no impact)
!	if (vlon(i,j) .lt. 0) vlon(i,j)=vlon(i,j)+360.

        tlat = vlat(i,j) * d2rad
        tlon = vlon(i,j) * d2rad
        relm = tlon - erlam0
        srlm = sin(relm)
        crlm = cos(relm)
        sph = sin(tlat)
        cph = cos(tlat)
        cc = cph * crlm
        tph = asin(cphi0 * sph - sphi0 * cc)
        rctph = 1.0 / cos(tph)
        cray(i,j) = sphi0 * srlm * rctph
        dray(i,j) = (cphi0 * cph + sphi0 * sph * crlm) * rctph

	if (j .eq. jm) then
	write(6,*) 'i,relm,cray,dray: ',i,relm,cray(i,j),dray(i,j)
	endif

      enddo
      enddo
	  end subroutine vecrot_rotlat

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine vecrot_rotlat_new(im,jm,tph0d,tlm0d,vlat,vlon,
     &                     cosalp,sinalp)

          
       real:: vlat(im,jm),vlon(im,jm),cosalp(im,jm),sinalp(im,jm)


      pi =3.141592654
      d2r=pi/180.
      r2d=1./d2r
c-----------------------------------------------------------------------
c***
c***  ready to compute rotation angles.
c***
c***  formulas use geodetic longitudes positive west
c***  so negate longitudes from etall (vlon and vloni)
c***  and the central longitudes of the grids
c***  since they arrive here as negative west.

!mp	is this sufficiently general?  assumption made that all will be west?

c***
      erlam0=-tlm0d*d2r
      erphi0=tph0d*d2r
      erl0_out=erlam0/d2r
      cphi0_out=cos(erphi0)
      sphi0_out=sin(erphi0)
c***
c***  compute earth rotation angles for winds relative to both
c***  the input and output grids on the output grid points.
c***
      do j=1,jm
      do i=1,im
        x=cphi0_out*cos(vlat(i,j)*d2r)*cos((-vlon(i,j)-erl0_out)*d2r)
     1   +sphi0_out*sin(vlat(i,j)*d2r)
        y=-cos(vlat(i,j)*d2r)*sin((-vlon(i,j)-erl0_out)*d2r)
        tvlon=atan(y/x)
        if(x.lt.0.)tvlon=tvlon+pi
        arg=sphi0_out*sin(tvlon)/cos(vlat(i,j)*d2r)
        arg=amin1(arg,1.)
        arg=amax1(arg,-1.)
	alpha=asin(arg)
	cosalp(i,j)=cos(alpha)
	sinalp(i,j)=sin(alpha)
c
	if (j .eq. jm) then
	write(6,*) 'i,cosalp,sinalp: ',i,cosalp(i,j),sinalp(i,j)
	endif
      enddo
      enddo

	  end subroutine vecrot_rotlat_new
