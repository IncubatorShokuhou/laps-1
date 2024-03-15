cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis
      subroutine supmap_local 
     +                (jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,jjlts,
     +                   jgrid,jus,jdot,ier)
c   supmap is used to plot map outlines according to one of nine projections.
c the origin and orientation of the projection are selected by the user.
c points on the earth defined by latitude and longitude are transformed to
c points in the u,v plane, the plane of projection.  the u and v axes are
c respectively parallel to the x and y axes of the plotter.  a rectangular
c frame parallel to the u and v axes is chosen and only material within the
c frame is plotted.
c   continental and u.s. state and county outlines are available for plotting.
c in addition, a western u.s. subset of the counties has been extracted, and
c up to three user files may be read.

c            name          date                 description
c       --------------  -----------     ---------------------------------------
c       ?? (ncar)          jan 69       "revised"
c                          may 71       "revised"
c                          oct 73       "standardized"
c                          jul 74       "revised"
c                          aug 76       "revised"
c       m. cairns       27 jan 81       added mapcol common for intensities
c       j. wakefield    19 jun 81       added mapdas common for dash patterns
c                       29 jun 81       removed trailing blanks and modified
c                                       to read unformatted files.
c                       28 sep 81       added common block supmpb.
c                       23 mar 82       added negative jdot function, to allow
c                                       suppression of supmap call printing if
c                                       jus=0.  see comments below.
c                       16 feb 83       added sysdisk to data file names and
c                                       added error code for open failure.
c                       10 apr 84       added user files option.
c                       14 jan 85       fixed bug in part < 0 option.
c       j. ramer        13 nov 85       added common blocks necessary for using
c                                       mapinv_gc with all projections and for
c                                       recreating a map from the contents of
c                                       the common blocks.  added subroutine
c                                       resup which recreates the map.
c       j. wakefield    19 nov 86       chgd sysdisk to lib_dev for files.
c
c       p. mcdonald     17 aug 93       fixed things up to run on stardent
c                                       and (hopefully) other unix machines.

c       usage
c                       call supmap (jproj,polat,polong,rrot,pl1,pl2,pl3,pl4,
c                                    jjlts,jgrid,jus,jdot,ier)

c       dimension of    pl1(2),pl2(2),pl3(2),pl4(2)
c       arguments

c       on input        jproj
c       for supmap        |jproj| defines the projection type
c                         according to the following code:
c                               1  stereographic
c                               2  orthographic
c                               3  lambert conformal conic with two standard
c                                  parallels
c                               4  lambert equal area
c                               5  gnomonic
c                               6  azimuthal equidistant
c                               7  dummy -- this code is not used
c                               8  cylindrical equidistant
c                               9  mercator
c                              10  mollweide type
c                         if jproj < 0, the map and grid lines are omitted.

c                       polat,polong,rrot
c                         if (|jproj|.ne.3)
c                         . polat and polong define in degrees the latitude and
c                           longitude of the point on the globe which is to
c                           transform to the origin of the u,v plane.
c                               -90 .le. polat .le. 90
c                              -180 .le. polong .le. 180
c                           degrees of latitude north of the equator and
c                           degrees of longitude east of the greenwich meridian
c                           are positive.  if the origin is at the north pole,
c                           "north" is considered to be in the direction of
c                           (polong+180.). if the origin is at the south pole,
c                           "north" is in the direction of polong.
c                         . rrot is the angle between the v axis and north at
c                           the origin.  it is measured in degrees and is taken
c                           to be positive if the angular movement from north
c                           to the v axis is counter-clockwise.  for the
c                           cylindrical projections (8,9,10), the axis of the
c                           projection is parallel to the v axis.
c                         if (|jproj|.eq.3) (lambert conformal conic)
c                         . polong = central meridian of projection in degrees.
c                         . polat,rrot are the two standard parallels in deg.

c                       jjlts,pl1,pl2,pl3,pl4
c                         |jjlts| can take the values 1 through 5 and specifies
c                         one of five options on the way in which the limits of
c                         the rectangular map are defined by the parameters
c                         pl1, pl2, pl3, and pl4.

c                         |jjlts| = 1
c                           the maximum useful area produced by the projection
c                           is plotted.  pl1, pl2, pl3, and pl4 are not used
c                           and may be set to zero.

c                         |jjlts| = 2
c                           in this case (pl1,pl2) and (pl3,pl4) are the
c                           latitudes and longitudes in degrees of two points
c                           which are to be at opposite corners of the map
c                           (upper right and lower left, respectively).
c                           care must be taken when using cylindrical
c                           projections and this option.

c                         |jjlts| = 3
c                           the minimum and maximum values of u and v are
c                           specified by pl1 through pl4.  pl1 = umin,
c                           pl2 = umax, pl3 = vmin, pl4 = vmax.  knowledge of
c                           the transformation equations is necessary for this
c                           option to be used (see below).

c                         |jjlts| = 4
c                           here pl1 = aumin, pl2 = aumax, pl3 = avmin,
c                           pl4 = avmax, where
c                              aumin = angular distance from origin to left
c                                      frame of map.
c                              aumax = angular distance from origin to right
c                                      frame of map.
c                              avmin = angular distance from origin to lower
c                                      frame.
c                              avmax = angular distance from origin to upper
c                                      frame.
c                           aumin, aumax, avmin, avmax must be positive and the
c                           origin must be within the rectangular limits of the
c                           map.  this option is useful for polar projections.
c                           it is not appropriate for the lambert conformal
c                           with two standard parallels.  an error message is
c                           printed if an attempt is made to use jjlts = 4 when
c                           jproj = 3, (see below).

c                         |jjlts| = 5
c                           pl1 through pl4 are two element arrays giving the
c                           latitudes and longitudes of four points which are
c                           to be on the four sides of the rectangular frame.
c                           pl1(1), pl1(2) are respectively the latitude and
c                           longitude of a point on the left frame.  similarly,
c                           pl2 lies on the right frame, pl3 lies on the lower
c                           frame and pl4 lies on the upper frame.  note that
c                           in the calling program pl1 through pl4 will be
c                           dimensioned:
c                               real  pl1(2),pl2(2),pl3(2),pl4(2)

c                          .if jjlts is positive, the supmap call is plotted
c                           below the map.  this is omitted if jjlts is < 0.

c                       jgrid
c                         |jgrid| gives in degrees the interval at which lines
c                         of latitude and longitude are to be plotted.  a value
c                         in the range 1 through 10 will usually be appropriate
c                         but higher values are acceptable.
c                         if jgrid < 0 the border around the map is omitted.
c                         if jgrid = 0 no grid lines are plotted.

c                       jus
c                         |jus|
c                         1  world outlines
c                         2  u.s. state outlines
c                         4  u.s. counties and states
c                         8  profs-region counties and states (western u.s.)
c                  16,32,64  supmap user file(s) 1,2,3
c                         to access a combination of datasets, use the sum of
c                         their values (e.g. 10 = u.s. + profs area counties).
c                         .if jus is positive, the supmap call and values of
c                          of umin, umax, vmin, vmax are printed as an aid to
c                          debugging.  this is omitted if jus is negative.

c                       jdot
c                         |jdot|
c                         0  for continuous outlines.
c                         1  for dotted outlines.
c                         .if jdot is negative, the supmap call is neither
c                          printed nor plotted.

c       on output       all arguments except ier are unchanged.
c       for supmap

c                       ier
c                         error flag with the following meanings:
c                         if ier =
c                          0  map successfully plotted.
c                         29  error opening data file.
c                         33  attempt to use non-existent projection.
c                         34  map limits inappropriate.
c                         35  angular limits too great.
c                         36  map has zero area.

c       entry points    maplot, supcon, supfst, supmap, suptrp, supvec, qcon,
c                       qvec, vecplt, resup

c                       maplot
c                         actually draws the map.

c                       supcon
c                         once the transformation has been set up by an initial
c                         call to supmap, the subroutine supcon may be called
c                         to transform a point, (latitude, longitude) to the
c                         corresponding point, (u, v) on the plane.  contours
c                         may thus be readily drawn against the map background.
c                         (see supfst and supvec below).

c                            call supcon(rlat,rlon,u,v)

c                         on input:
c                           rlat,rlon are the latitude and longitude of a point
c                           to be transformed to the u,v plane.
c                            -90. .le. rlat .le.  90.
c                           -180. .le. rlon .le. 180.

c                         on output:
c                           rlat,rlon are unchanged.
c                           u,v are the transformed coordinates of the point.

c                       qcon
c                         actually performs the above mentioned transformation.

c                       supfst
c                       supvec
c                         to facilitate drawing lines on the map these routines
c                         which act like the plotting routines frstpt and
c                         vector are included.  they are subject to the same
c                         restrictions as supcon above.

c                            call supfst (rlat,rlon)
c                            call supvec (rlat,rlon)

c                       qvec
c                         decides what lines are to be drawn and where.

c                       suptrp
c                         performs interpolation to the edges of the frame.

c                       vecplt
c                         called by qvec to draw (dot) lines on the plotter.

c       common blocks    name   length
c                       supmp1  9 fp
c                       supmp2  1 int + 204 fp
c                       supmp3  2 int + 5 fp
c                       supmp4  5 int + 2 fp
c                       supmp5  1 int + 5 fp
c                       supmp6  6 fp
c                       supmp7  3 int + 2 fp
c                       supmp8  6 fp
c                       supmp9  3 fp
c                       supmpa  1 int
c                       supmpb  4 fp
c                       supmpc  4 int + 1 fp
c                       supmpd  9 fp
c                       mapcol  4 int
c                       mapdas  4 int

c       i/o             map plotted.  outline data is read from any of several
c                       disk files.  supmap call printed.

c       precision       single

c       language        fortran

c       algorithm       the latitudes and longitudes of successive outline
c                       points are transformed to coordinates in the plane of
c                       projection and joined by a vector.

c       references      hershey, a. v., the plotting of maps on a crt printer.
c                         nwl report no. 1844, 1963.
c                       lee, tso-hwa, students summary reports, work-study
c                         program in scientific computing.  ncar 1968.
c                       parker, r. l., ucsd supermap:  world plotting package.
c                       steers, j.a., an introduction to the study of map
c                         projections.  univ. of london press, 1962.

c       accuracy        the definition of the map produced is limited by the
c                       fact that the resolution of the virtual plotter space
c                       is 1024 units in the x and y directions.

c       plotting routines       pwrt, frstpt, vector, point, dashln, perim, set
c       used

c       required resident       atan, tan, sin, cos, alog, sqrt, atan2, acos
c       routines

        common/supmp1/pi,tovpi,dtr,rtd,eps,ov90,con1,con2,part
        common/supmp3/polong,cone,rlat,rlon,jgr,ilf,sgn
        common/supmp4/ifst,igo,igold,icross,iout,uold,vold
        common/supmp5/phioc,sino,coso,sinr,cosr,iproj
        common/supmp6/umin,umax,vmin,vmax,ueps,veps
        common/supmp7/phio,phia,igrid,idot,ilts
        common/supmp8/u,v,u1,v1,u2,v2
        common/supmp9/ds,di,dsrdi
        common/supmpa/iier
        common/supmpb/x1,y1,x2,y2
        common/supmpc/lproj,rot,jlts,lgrid,ius
        common/supmpd/aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii

      dimension pl1(2),pl2(2),pl3(2),pl4(2),laba(20),labb(18)

      real lat1,lat2,aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii

      equivalence (pla1,aumin),(pla2,aumax),(pla3,avmin),(pla4,avmax),
     1            (phia,lat1),(rot,lat2)

      data pltres,reslim/ 1024.,10./

!     external supmbd

      squ(x) = x*x
c
c       supmbd subroutine changed to block data
c     igdflt = %loc(supmbd)             !force load of block data from library

      write(6,*)' subroutine supmap...',jproj

      rot = rrot
      ilts = iabs(jjlts)
      jlts = jjlts
      lgrid = jgrid
      igrid = iabs(lgrid)
c     jgr = isign(1,lgrid)
      jgr = lgrid
      ius = jus
      iusgn = isign(1,ius)
      pla1 = pl1(1)
      pla2 = pl2(1)
      pla3 = pl3(1)
      pla4 = pl4(1)
      labl = 77
      idot = iabs(jdot)

c initialization
      lproj = jproj
      iproj = iabs(lproj)
      jpr = isign(1,lproj)
c     jgr=jpr
      iout = abs(ius)           !use iout for choosing data file(s)
      phia = polat
      polong = polon
      phio = polon
      phioc = 540.-phio
      icross = 0
      iier = 0
      ilf = 0

c compute constants appropriate to each projection
      if (iproj .ne. 3) go to 30

c lambert conformal conic
      sgn = sign(1.,0.5*(lat1+lat2))
      chi1 = (90.-sgn*lat1)*dtr
      if (lat1 .eq. lat2) go to 20
      chi2 = (90.-sgn*lat2)*dtr
      cone = alog(sin(chi1)/sin(chi2))/alog(tan(0.5*chi1)/tan(0.5*chi2))
      go to 60
   20 cone = cos(chi1)
      go to 60

c the others
   30 xt = rot*dtr
      sinr = sin(xt)
      cosr = cos(xt)
      xt = phia*dtr
      sino = sin(xt)
      coso = cos(xt)

c  calculate coefficients necessary to convert displacement vector on earth's
c  surface from cartesian system with origin at center of earth, unit vector k
c  k at the projection center, and unit vector j in the positive u direction to
c  a system with unit vector k at north pole and unit vector i at intersection
c  of equator and greenwich meridian.  needed later for mapinv_gc.
      aaa=cos(rot*dtr)*cos(polong*dtr)*sin(phia*dtr)-
     &    sin(rot*dtr)*sin(polong*dtr)
      bbb=cos(rot*dtr)*sin(polong*dtr)*sin(phia*dtr)+
     &    sin(rot*dtr)*cos(polong*dtr)
      ccc=-cos(rot*dtr)*cos(phia*dtr)
      ddd=-sin(rot*dtr)*cos(polong*dtr)*sin(phia*dtr)-
     &    cos(rot*dtr)*sin(polong*dtr)
      eee=-sin(rot*dtr)*sin(polong*dtr)*sin(phia*dtr)+
     &    cos(rot*dtr)*cos(polong*dtr)
      fff=sin(rot*dtr)*cos(phia*dtr)
      ggg=cos(polong*dtr)*cos(phia*dtr)
      hhh=sin(polong*dtr)*cos(phia*dtr)
      iii=sin(phia*dtr)

      if (iproj-7) 60,67,40

c cylindrical projections               [8,9,10]
   40 if (phia .ne. 0.0) go to 42
      if (rot .eq. 0.0) go to 45
      if (abs(rot) .eq. 180.) go to 50
   42 sino1 = coso*cosr
      coso1 = sqrt(con1-sino1*sino1)
      ovc1 = 1./coso1
      phio = phio-atan2(sinr*ovc1,-cosr*sino*ovc1)*rtd
      phioc = 540.-phio
      sinr = sinr*coso*ovc1
      cosr = -sino*ovc1
      sino = sino1
      coso = coso1
      go to 60

c use simple transforms for cylindrical projections if rot = polat = .0
c   i.e. iproj = 11, 12, 13
   45 sino = 1.0
      iproj = iproj+3
      go to 55

   50 sino = -1.0
      phio = phio+180.
      phioc = phioc+180.
   55 coso = 0.0
      sinr = 0.0
      cosr = 1.0
      ilf = 1

c ilts = 1         the maximum useful area is plotted.
c ---------
   60 go to (61 ,62 ,62 ,61 ,61 ,66 ,67 ,68 ,66 ,70 ,
     1       68 ,66 ,70 ),iproj

c stereographic                         [1]
c lambert equal area                    [4]
c gnomonic                              [5]
   61 umin = -2.0
      vmin = -2.0
      umax = 2.0
      vmax = 2.0
      go to 80

c orthographic                          [2]
c lambert conformal conic               [3]
   62 umin = -1.0
      vmin = -1.0
      umax = 1.0
      vmax = 1.0
      go to 80

c azimuthal equidistant                 [6]
c mercator with arbitrary pole          [9]
c mercator                              [12]
   66 umax = pi
      vmax = pi
      umin = -pi
      vmin = -pi
      go to 80

c dummy  --  error exit                 [7]
   67 iier = 33
      call uliber2 (iier
     1           ,' supmap-attempt to use non-existent projection')     
      go to 700

c cylindrical equidistant               [8,11]
   68 umax = 180.
      umin = -180.
      vmax = 90.
      vmin = -90.
      go to 80

c mollweide type                        [10,13]
   70 umax = 2.0
      umin = -2.0
      vmax = 1.0
      vmin = -1.0

   80 ueps = 0.5*(umax-umin)
      veps = 0.5*(vmax-vmin)
      if (iproj .eq. 3) ueps = 180.

c compute the appropriate map boundaries.
      go to (600,200,300,400,500),ilts

c ilts = 2         point (pl1,pl2) in upper right corner , (pl3,pl4) in
c ---------        lower left corner of plot.
  200 rlat = pla1
      rlon = pla2
      call qcon
      u1 = u
      v1 = v
      rlat = pla3
      rlon = pla4
      call qcon
      umax = amax1(u1,u)
      umin = amin1(u1,u)
      vmax = amax1(v1,v)
      vmin = amin1(v1,v)
      go to 600

c ilts = 3         set plot limits directly.
c ----------
  300 umax = pla2
      umin = pla1
      vmax = pla4
      vmin = pla3
      go to 600

c ilts = 4         use angular distances to set plot limits.
c ----------
  400 cosumi = cos(aumin*dtr)
      sinumi = sqrt(con1-cosumi*cosumi)
      cosuma = cos(aumax*dtr)
      sinuma = sqrt(con1-cosuma*cosuma)
      cosvmi = cos(avmin*dtr)
      sinvmi = sqrt(con1-cosvmi*cosvmi)
      cosvma = cos(avmax*dtr)
      sinvma = sqrt(con1-cosvma*cosvma)

      go to (401,402,403,404,405,406,407,408,409,410,
     1       408,409,410),iproj

c stereographic                         [1]
  401 umax = (1.-cosuma)/sinuma
      umin = -(1.-cosumi)/sinumi
      vmax = (1.-cosvma)/sinvma
      vmin = -(1.-cosvmi)/sinvmi
      go to 600

c orthographic                          [2]
  402 if (amax1(aumin,aumax,avmin,avmax) .gt. 90.) go to 900
      umax = sinuma
      umin = -sinumi
      vmax = sinvma
      vmin = -sinvmi
      go to 600

c lambert conformal conic               [3]
  403 iier = 34
      call uliber2 (iier,' supmap-map limits inappropriate')
      go to 700

c lambert equal area                    [4]
  404 umax = (1.+cosuma)/sinuma
      umin = (1.+cosumi)/sinumi
      vmax = (1.+cosvma)/sinvma
      vmin = (1.+cosvmi)/sinvmi
      umax = 2./sqrt(1.+umax*umax)
      umin = -2./sqrt(1.+umin*umin)
      vmax = 2./sqrt(1.+vmax*vmax)
      vmin = -2./sqrt(1.+vmin*vmin)
      go to 600

c gnomonic                              [5]
  405 if (amax1(aumin,aumax,avmin,avmax) .ge. 90.) go to 900
      umax = sinuma/cosuma
      umin = -sinumi/cosumi
      vmax = sinvma/cosvma
      vmin = -sinvmi/cosvmi
      go to 600

c azimuthal equidistant                 [6]
  406 umax = aumax*dtr
      umin = -aumin*dtr
      vmax = avmax*dtr
      vmin = -avmin*dtr
      go to 600

c dummy  --  error exit                 [7]
  407 go to 67

c cylindrical equidistant               [8,11]
  408 umax = aumax
      umin = -aumin
      vmax = avmax
      vmin = -avmin
      go to 600

c mercator                              [9,12]
  409 if (amax1(avmin,avmax) .ge. 90.) go to 900
      umax = aumax*dtr
      umin = -aumin*dtr
      vmax = alog((1.+sinvma)/cosvma)
      vmin = -alog((1.+sinvmi)/cosvmi)
      go to 600

c mollweide type                        [10,13]
  410 umax = aumax*ov90
      umin = -aumin*ov90
      vmax = sinvma
      vmin = -sinvmi
      go to 600

c ilts = 5         use four edge points to set limits.
c ----------
  500 plb1 = pl1(2)
      rlat = pla1
      rlon = plb1+eps
      call qcon
      umin = u
      plb2 = pl2(2)
      rlat = pla2
      rlon = plb2-eps
      call qcon
      umax = u
      plb3 = pl3(2)
      rlat = pla3
      rlon = plb3
      call qcon
      vmin = v
      plb4 = pl4(2)
      rlat = pla4
      rlon = plb4
      call qcon
      vmax = v

c compute map limits for plot
  600 du = umax-umin
      dv = vmax-vmin

c error if map has zero area
      if (du.eq.0.0 .or. dv.eq.0.0) go to 905
        if(part.le..0)goto 620          !user-supplied x1,y1,x2,y2
      if (du .gt. dv) go to 610
      y1 = 0.5*(1.-part)
      y2 = 1.-y1
      x1 = 0.5*(1.-part*du/dv)
      x2 = 1.-x1
      go to 620

  610 x1 = 0.5*(1.-part)
      x2 = 1.-x1
      y1 = 0.5*(1.-part*dv/du)
      y2 = 1.-y1

c error if map has essentially zero area
  620 if (amin1(x2-x1,y2-y1)*pltres .lt. reslim) go to 905
      call set (x1,x2,y1,y2,umin,umax,vmin,vmax,1)
      call line(umin,vmin,umin,vmax)
      call line(umin,vmax,umax,vmax)
      call line(umax,vmax,umax,vmin)
      call line(umax,vmin,umin,vmin)
      ds = squ(((x2-x1)*pltres)/du)
      dsrdi = sqrt(ds/di)
      write(6,*)' dsrdi = ',dsrdi

c do we write anything?
      if (jlts.lt.0 .and. iusgn.lt.0) go to 640
        if(jdot.lt.0)goto 640

c create the label
c     encode (labl,7000,laba(1)) lproj,phia,polong,rot,pla1,pla2,pla3,
c    1                           pla4,jlts,lgrid,ius,idot
c7000 format (8h supmap(,i3,7(1h,,f6.1),4(1h,,i3),1h))

c     if (ilts .eq. 5) encode (61,7010,labb(1)) plb1,plb2,plb3,plb4
c7010 format (7x,1h(,24x,4(1h,,f6.1),1h))

      if (jlts .lt. 0) go to 630

c write supmap call beneath the map
      call pwrt (240,17,laba(1),labl,0,0)
      if (ilts .eq. 5) call pwrt (240,1,labb(1),61,0,0)

  630 if (iusgn .lt. 0) go to 640

c print out the call et al.
      k = (labl+3)/4
      write (6,6000) (laba(i),i=1,k)
 6000 format (30a4)
      if (ilts .eq. 5) write (6,6000) (labb(i),i=1,18)
      write (6,6010) umin,umax,vmin,vmax
 6010 format (8h umin = ,f11.6,9h  umax = ,f11.6,9h  vmin = ,
     +        f11.6,9h  vmax = ,f11.6)

c draw the map
  640 if (iout.ne.0 .or. igrid.ne.0 .or. jgr.ge.0) call maplot_local
      idot = 0

c return ier
  700 ier = iier
      return

c error returns
  900 iier = 35
      call uliber2 (iier,' supmap-angular limits too great')
      go to 700
  905 iier = 36
      call uliber2 (iier,' supmap-map has zero area')
      go to 700

      end
c-------------------------------------------------------------------------------
      subroutine maplot_local
c this subroutine plots the continental and u.s. state outlines,
c meridians, parallels, limbs where appropriate. it labels key meridians
c and poles, and it draws a border.

        common/supmp1/pi,tovpi,dtr,rtd,eps,ov90,con1,con2,part
        common/supmp2/npts,maxlat,minlat,maxlon,minlon,pts(1200000)
        common/supmp3/polong,cone,rlat,rlon,jgr,ilf,sgn
        common/supmp4/ifst,igo,igold,icross,iout,uold,vold
        common/supmp5/phioc,sino,coso,sinr,cosr,iproj
        common/supmp6/umin,umax,vmin,vmax,ueps,veps
        common/supmp7/phio,phia,igrid,idot,ilts
        common/supmp8/u,v,u1,v1,u2,v2
        common/supmpa/iier
        common/mapcol/mpcol1,mpcol2,mpcol3,mpcol4
        common/mapdas/ldash1,ldash2,ldash3,ldash4
c       data                    !default line intensities and dash patterns
c    1          mpcol1,ldash1   /255,'1777'o/,  !map lines
c    2          mpcol2,ldash2   /128,'1756'o/,  !grid lines
c    3          mpcol3,ldash3   /192,'1777'o/,  !limb lines
c    4          mpcol4,ldash4   /255,'1777'o/   !perimeter

        character       supmap_dir*150
        integer       lsdir

        data iwrite /0/

        character*1000 c_line

      dimension splat(2)
      real maxlat,minlat,maxlon,minlon,midlat,midlon
      character*180 namfil
      data   sinlmb,coslmb /0.017452406, 0.99984765/
      data floorc / 10000. /

      floor(x) = aint(x+floorc)-floorc
      cling(x) = floor(x)+1.


      call get_directory('static',supmap_dir,lsdir)

      supmap_dir = supmap_dir(1:lsdir)//'/ncarg/'
c        supmap_dir = '../static/ncarg/'
      lsdir = index (supmap_dir, ' ') - 1

      call dashln(ldash1)
        call optn(2hin,mpcol1)
      go to (10,10,10,5,10,5,905,5,5,5,5,5,5),iproj
    5 icf=1
   10 if(iout.eq.0)goto 100

c***select appropriate file(s) according to bits set in iout
        ibs=1                           !start with 1st bit
        mask=1
   11   do 12 ibit=ibs,16               !check lower 16 bits
         ibs=ibs+1
c        if((iout.and.mask).ne.0)goto 13        !exit loop
         if(iand(iout,mask).ne.0)goto 13        !exit loop
         mask=mask*2
   12   continue
        go to 100                       !fell out of loop so done
   13   mask=mask*2

        if(ibit .eq. 1)goto14
        if(ibit .eq. 2)goto15
        if(ibit .eq. 3)goto16
        if(ibit .eq. 4)goto17
        if(ibit .eq. 5)goto18
        if(ibit .eq. 6)goto19
        if(ibit .eq. 7)goto20
        go to 100                       !out of range -- return

c  14   if((iout.and.2).ne.0)then       ! states with continents?
   14   if(iand(iout,2).ne.0)then       ! states with continents?
         ibs=ibs+1      !skip bit later
         mask=mask*2
c        namfil='lib_dev:[gudat]conandsta.dat'  !3 (special case)
c       namfil = '/home/star1/b/mcdonald/data/supmap/conandsta.dat'
        namfil = supmap_dir (1:lsdir) // 'conandsta.dat'
        else
c        namfil='lib_dev:[gudat]continent.dat'  !1
c       namfil = '/home/star1/b/mcdonald/data/supmap/continent.dat'
        namfil = supmap_dir (1:lsdir) // 'continent_minus_us.dat'
        endif
        go to 25
c  15   namfil='lib_dev:[gudat]state.dat'       !2
c  15 namfil = '/home/star1/b/mcdonald/data/supmap/state.dat'
   15   namfil = supmap_dir (1:lsdir) // 'state.dat'
        go to 25
c  16   namfil='lib_dev:[gudat]uscounty.dat'    !4
c  16 namfil = '/home/star1/b/mcdonald/data/supmap/uscounty.dat'
   16   namfil = supmap_dir (1:lsdir) // 'uscounty.dat'
        go to 25
c  17   namfil='lib_dev:[gudat]county.dat'      !8
c  17 namfil = '/home/star1/b/mcdonald/data/supmap/county.dat'
   17   namfil = supmap_dir (1:lsdir) // 'state_from_counties.dat'
        go to 25
c  18   namfil='supmap_userfile1'               !16
   18   namfil = supmap_dir (1:lsdir) // 'userfile1.dat'
        goto 25
c  19   namfil='supmap_userfile2'               !32
   19   namfil = supmap_dir (1:lsdir) // 'userfile2.dat'
        goto 25
c  20   namfil='supmap_userfile3'               !64
   20   namfil = supmap_dir (1:lsdir) // 'userfile3.dat'
        goto 25

c***open file
c  25   open(3,name=namfil,type='old',form='unformatted',readonly,err=26)
   25   continue
        write(6,*)namfil(1:60)
        open(3,file=namfil,status='old',form='unformatted',err=26)
        goto 30

c***error opening file -- return
   26   iier=29
        write(6,*)'supmap - error opening: ',namfil
        return

c***read next line
   30   continue

        if(.true.)then
            if(iwrite .eq. 0)
     1        write(6,*)' using simple read to read binary map info...'
            read(3,end=99)npts,maxlat,minlat,maxlon,minlon
     1                   ,(pts(m),m=1,npts)

        else ! this may be needed for dec alpha but will not work on linux
            if(iwrite .eq. 0)
     1        write(6,*)' using cio.c to read binary map info...'
            read(3,end=99,err=41)npts,maxlat,minlat,maxlon,minlon
     1                   ,(pts(m),m=1,200)
 41         continue

!           convert between bigendian and littleendian (under construction)
            call in_to_im(4,4,npts,1)
            call in_to_im(4,4,pts,npts)
            call in_to_im(4,4,maxlat,1)
            call in_to_im(4,4,minlat,1)
            call in_to_im(4,4,maxlon,1)
            call in_to_im(4,4,minlon,1)

        endif

        iwrite = iwrite + 1

   99   npts=npts/2
      if (npts .eq. 0)then
        close(3)
        goto 11                 !check next bit
      endif
      if (npts .le. 16) go to 70
      if (icf .ne. 0) go to 70

c does this line intersect the screen?
c       1  2  3
c       4     5
c       6  7  8
      midlat = (maxlat+minlat)*0.5
      midlon = (maxlon+minlon)*0.5
      rlat = maxlat
      rlon = maxlon
      call qcon
      x3 = u
      y3 = v
      rlon = midlon
      call qcon
      x2 = u
      y2 = v
      rlon = minlon
      call qcon
      x1 = u
      y1 = v
      rlat = midlat
      call qcon
      x4 = u
      y4 = v
      rlon = maxlon
      call qcon
      x5 = u
      y5 = v
      rlat = minlat
      call qcon
      x8 = u
      y8 = v
      rlon = midlon
      call qcon
      x7 = u
      y7 = v
      rlon = minlon
      call qcon
      x6 = u
      y6 = v
      xmn = amin1(x1,x2,x3,x4,x5,x6,x7,x8)
      xmx = amax1(x1,x2,x3,x4,x5,x6,x7,x8)
      ymn = amin1(y1,y2,y3,y4,y5,y6,y7,y8)
      ymx = amax1(y1,y2,y3,y4,y5,y6,y7,y8)

      dx = amin1(xmx-xmn,180.)
      dy = amin1(ymx-ymn,180.)
      xmx = xmx+0.01*dx
      xmn = xmn-0.01*dx
      ymx = ymx+0.01*dy
      ymn = ymn-0.01*dy

      if (xmn.gt.umax .or. xmx.lt.umin .or. ymn.gt.vmax .or.
     1    ymx.lt.vmin) go to 30

   70 rlat = pts(1)
      rlon = pts(2)
      ifst = 1
      igo = 0
      call qvec
      do 75 j=2,npts
         rlat = pts(2*j-1)
         rlon = pts(2*j)
         call qvec
   75 continue

      go to 30

c***done plotting -- draw and label grid, if desired
  100 splat(2) = 90.
      splat(1) = -90.
      if (igrid .eq. 0) go to 300
      idots = idot
      idot = 0

c     if (idot .eq. 0) go to 200

c letter key meridians and poles

c       south pole
      ispf = 0
      rlat = -90.
      rlon = 0.0
      call qcon
      if ((u .gt. umax) .or. (u .lt. umin) .or. (v .gt. vmax) .or.
     1    (v .lt. vmin)) go to 110
      usp = u
      vsp = v
      ispf = 1
      ipf = 1
      if (phia .lt. 0) call pwrt (u,v,2hsp,2,1,0)

c       north pole
110   inpf = 0
      ipf = 0
      rlat = 90.
      call qcon
      if ((u .gt. umax) .or. (u .lt. umin) .or. (v .gt. vmax) .or.
     1    (v .lt. vmin)) go to 120
      unp = u
      vnp = v
      inpf = 1
      ipf = 1
      if (phia .gt. 0) call pwrt (u,v,2hnp,2,1,0)

c       equator
  120 rlon = phio-10.
      rlat = 0.0
      do 125 i=1,36
         rlon = rlon+10.
         call qcon
         if (u.le.umax .and. u.ge.umin .and. v.le.vmax .and. v.ge.vmin)
     1       go to 130
  125 continue
      go to 140

  130 call pwrt (u,v,2heq,2,1,0)

c       greenwich meridian
  140 rlat = 85.
      rlon = 0.0
      do 145 i=1,16
         rlat = rlat-10.
         call qcon
         if (u.le.umax .and. u.ge.umin .and. v.le.vmax .and. v.ge.vmin)
     1       go to 150
  145 continue
      go to 160

  150 call pwrt (u,v,2hgm,2,1,0)

c       date line
  160 rlat = 85.
      rlon = 180.
      do 165 i=1,16
         rlat = rlat-10.
         call qcon
         if (u.le.umax .and. u.ge.umin .and. v.le.vmax .and. v.ge.vmin)
     1       go to 170
  165 continue
      go to 200

  170 call pwrt (u,v,1hi,1,1,0)

  200 rgrid = igrid
c      call dashln (ldash2)
c       call optn(2hin,mpcol2)

c should we bother limiting grid points transformed?
      if (icf .ne. 0) go to 270
      if (iproj.ge.8 .and. iproj.le.10) go to 270

c set up to find extrema
      dlon = rgrid
      stlon = floor(polong/rgrid)*rgrid
      if (ispf.ne.0 .and. inpf.eq.0) stlon = stlon+180.
      rlon = stlon-dlon
      splon = stlon+360.
      j = 0
      psign = 1.

c check for south pole
      if (ispf .ne. 0) psign = -1.

c do we grid poles specially?
      splat(2) = 90.*psign
      splat(1) = splat(2)

c if both poles within frame jump.
      if (inpf.ne.0 .and. ispf.ne.0) go to 270

c if either in frame use as base
      if (inpf.ne.0 .or. ispf.ne.0) go to 230

c no pole is close to the window
      j = -1
      splat(2) = floor(phia/rgrid)*rgrid
      if (abs(splat(2)) .eq. 90.) splat(2) = 0.0

c search for first point within frame.
  210 rlon = rlon+dlon
      dlat = rgrid
      rlat = splat(2)-dlat
  215 rlat = rlat+dlat
      call qcon
      if ((u .le. umax) .and. (u .ge. umin) .and. (v .le. vmax) .and.
     1    (v .ge. vmin)) go to 225
      if (abs(rlat) .lt. 90.) go to 215
      if (dlat .lt. 0.0) go to 220

c reverse latitude search direction
      rlat = splat(2)+dlat
      dlat = -dlat
      go to 215

c update longitude ! quit.
  220 j = 0
      if (rlon-splon) 210,300,300

c set up for limit search
  225 j = j+1
      stlon = rlon
      rlon = stlon-dlon
      if (rlat .eq. 0.0) rlat = sign(rlat,-psign)
      splat(2) = rlat
      splat(1) = splat(2)

      splat(1) = 90.0
      splat(2) = -90.0
      stlon = stlon - dlon
      splon = splon + dlon
      rlon  = stlon
  226 if (rlon .le. splon) then
          rlat2 = 80.0
          if (amod(rlon,90.0) .eq. 0.0) rlat2 = 90.0
          rlat  = -rlat2
          ifst  = 1
          igo   = 0
          call qvec
  227     rlat = rlat + 1.0
          if (rlat .le. rlat2) then
              if (iproj .eq. 1) rlat = rlat2
              call qvec
              if (igo .ne. 0) then
                  if (rlat .lt. splat(1)) splat(1) = rlat
                  if (rlat .gt. splat(2)) splat(2) = rlat
              endif
              go to 227
          endif
          rlon  = rlon + dlon
          go to 226
      endif
      if (dlon .ne. 0.0) go to 285

c longitude loop
c       igf     flag to signal no points within window.
c       ipf     flag signals whether a pole lies within the frame.
c       ilf     flag signals whether to plot complete longitudes
c               (i.e. to pole for all latitudes.)
  230 rlon = rlon+dlon
      if (rlon.ge.splon .or. rlon.lt.stlon) go to 285
      i1 = ipf
      i2 = mod(i1+1,2)
      tsa = psign
      dlat = -psign
      dx = amod(90.,rgrid)
      if (dx .eq. 0.0) dx = rgrid
      xlat = 90.-dx
      if (ilf.ne.0 .or. amod(rlon,90.).eq.0.0) xlat = 90.
      olat = sign(amin1(abs(splat(i2+1)),xlat),splat(i2+1))
      igf = 0
  235 ifst = 1
      igo = 0
      rlat = olat
      call qvec

c latitude loop.
  240 rlat = rlat+dlat
      igf = max0(igo,igf)
      call qvec
      if (igo .ne. 0) go to 245

c this point outside the frame
      if (rlat*tsa .le. splat(i1+1)*tsa) go to 250
  245 if (abs(rlat) .lt. xlat) go to 240
      rlat = sign(amax1(abs(splat(i1+1)),xlat),splat(i1+1))

c possible new latitude extreme.
  250 splat(i1+1) = rlat

c reverse latitude search direction
      i1 = i2
      i2 = mod(i1+1,2)
      tsa = -psign
      dlat = psign
      if (i1 .ne. 0) go to 235

c latitude loop finished.
      if (abs(splat(i2+1)) .lt. 90.) go to 255
      ipf = 1
      psign = sign(1.,splat(i2+1))
      splat(i2+1) = splat(i1+1)
      splat(i1+1) = 90.*psign
  255 if (igf .ne. 0) go to 230

c longitude extreme reached.
      if (j .ne. 0) go to 260

c change longitude direction.
      j = 1
      splon = rlon
      rlon = stlon
      dlon = -dlon
      stlon = splon-360.
      go to 230

c set up last longitude extreme
  260 if (dlon .lt. 0.0) go to 265
      splon = rlon
      go to 285
  265 stlon = rlon
      go to 285

c draw all meridians.
  270 dlon = rgrid
      stlon = 0.0
      splon = 360.
      rlon = 0.0
      splat(2) = 90.
      splat(1) = -90.
      dx = amod(90.,rgrid)
      if (dx .eq. 0.0) dx = rgrid
      olat = 90.-dx

  275 rlon = rlon+dlon
      igo = 0
      ifst = 1
      xlat = olat

      if (ilf.ne.0 .or. amod(rlon,90.).eq.0.0) then
         if (phia.lt.0) then
            rlat = 0.0
            xlat = 90.0

         else if (phia.ge.0) then
            rlat = 90.0
            xlat = 0.0

         end if

      else
         rlat = xlat

      end if

      call qvec
  280 rlat = rlat-1.
      call qvec
      if (rlat .gt. -xlat) go to 280
      if (rlon .lt. splon) go to 275

c draw parallels
  285 dlat = rgrid
      rlat = amin1(splat(2),splat(1))
      olat = amax1(splat(2),splat(1))
      splat(2) = floor(rlat/rgrid)*rgrid
      splat(1) = amin1(cling(olat/rgrid)*rgrid,90.)
      rlat = amax1(dlat-90.,splat(2))-dlat
      olat = amin1(90.-dlat,splat(1))
c      if (dlon .gt. 0.0) then
c          stlon = stlon - dlon
c      else
c          stlon = stlon + dlon
c      endif
      dlon = 1.
      if (ilf .ne. 0) ipf = 0

  290 rlat = rlat+dlat
      if (ipf .ne. 0) dlon = 1./cos(dtr*rlat)
      igo = 0
      ifst = 1
      rlon = stlon
      call qvec

  295 rlon = rlon+dlon
      call qvec
      if (rlon .le. splon) go to 295
      if (rlat .lt. olat) go to 290

      idot = idots
c      call dashln (ldash3)
c       call optn(2hin,mpcol3)

c draw limb lines
  300 idots = idot
      idot = 0
      go to (400,330,305,335,400,340,400,400,400,345,
     1       400,400,345),iproj

c lambert conformal conic           [3]
  305 continue
      go to 400 ! test to eliminate spurious line
      dlat = 1.
      rlon = phio+con2
      olat = amax1(-90.,splat(2)-dlat)
      k = cling(splat(1)-splat(2))
      do 320 i=1,2
         igo = 0
         ifst = 1
         rlat = olat
         call qvec
         do 310 j=1,k
            rlat = rlat+dlat
            call qvec
  310    continue
         rlon = phio-con2
  320 continue
      go to 400

c orthographic                  [2]
  330 radius = 1.
      axis = 1.
      go to 350

c lambert equal area            [4]
  335 radius = 2.
      axis = 1.
      go to 350

c azimuthal equdistant          [6]
  340 radius = pi
      axis = 1.
      go to 350

c mollweide                     [10,13]
  345 radius = 2.
      axis = 0.5

  350 u = radius
      v = 0.0
      w = 0.0
      igo = 0
      ifst = 1
      do 370 i=1,361
         v = axis*v
         if (u.le.umax .and. u.ge.umin .and. v.le.vmax .and. v.ge.vmin)
     1       go to 355
         igo = 0
         go to 365
  355    if (igo .ne. 0) go to 360
         call frstpt (u,v)
         igo = 1
         go to 365

  360    call vector (u,v)
  365    v = u*sinlmb+w*coslmb
         u = u*coslmb-w*sinlmb
         w = v
  370 continue

c draw border
  400 if (jgr .gt. 0)then
c       call dashln(ldash4)
c       call optn(2hin,mpcol4)
c       call perim(1,1,1,1)
c        call frstpt (umin,vmin)
c        call vector (umax,vmin)
c        call vector (umax,vmax)
c        call vector (umin,vmax)
c        call vector (umin,vmin)
      endif
      idot = idots
      return

  905 return
      end

c-------------------------------------------------------------------------------
      subroutine qcon
c this subroutine transforms the point (rlat,rlon), in degrees,
c to (u,v) on the map plane dependent upon the projection, iproj.

        common/supmp1/pi,tovpi,dtr,rtd,eps,ov90,con1,con2,part
        common/supmp3/polong,cone,rlat,rlon,jgr,ilf,sgn
        common/supmp5/phioc,sino,coso,sinr,cosr,iproj
        common/supmp6/umin,umax,vmin,vmax,ueps,veps
        common/supmp8/u,v,u1,v1,u2,v2
        common/supmpa/iier
c
c
      common /supmp7/   phio,phia,igrid,idot,ilts
c
      data      resl    / 2.9765625 /
c
c
      data oldu,oldv / 0.,0./

      u = amod(rlon+phioc,360.)-180.

      go to (50 ,50 ,130,50 ,50 ,50 ,170,50 ,50 ,50 ,
     1       210,220,230),iproj

   50 t1 = u*dtr
      t2 = rlat*dtr
      sinph = sin(t1)
      sinla = sin(t2)
      cosph = cos(t1)
      cosla = sqrt(1.-sinla*sinla)
      tcos = cosla*cosph
      cosa = sinla*sino+tcos*coso
      sina = sqrt(con1-cosa*cosa)

c patch to avoid divide by zero
      ovsina = 1.e12
      if(sina.ne.0.0)ovsina = 1./sina
c end patch

      sinb = cosla*sinph*ovsina
      cosb = (sinla*coso-tcos*sino)*ovsina

c perform transformation appropriate to the projection
      go to (110,120,130,140,150,160,170,180,190,200),iproj

c stereographic                         [1]
  110 r = (1.-cosa)*ovsina
        a = atan2 (sina, cosa) * 0.5
        r = tan (a)
      go to 300
c 110 ihemi = 1
c     if (phia .lt. 0.0) ihemi = -1
c     call maproj_poster (rlat, rlon, u, v, resl, polong, ihemi, 1)
c     v     = 8193 - v
c     go to 305

c orthographic                          [2]
  120 r = sina
      if (cosa) 320,320,300

c lambert conformal conic               [3]
  130 udif = abs(u-oldu)
      oldu = u
      chi = 90.-sgn*rlat
      if (chi .ge. con2) go to 320
      r = tan(0.5*dtr*chi)**cone
      u = u*cone*dtr
      v = -r*sgn*cos(u)
      u = r*sin(u)
      go to 310

c lambert equal area                    [4]
  140 if (abs(cosa+1.) .lt. 1.e-6) go to 320
      r = (1.+cosa)*ovsina
      r = 2./sqrt(1.+r*r)
      go to 300

c gnomonic                              [5]
  150 if (cosa .le. 0.0) go to 320
      r = sina/cosa
      go to 300

c azimuthal equididsant                 [6]
  160 if (abs(cosa+1.) .lt. 1.e-6) go to 320
      r = acos(cosa)
      go to 300

c dummy   --  error                     [7]
  170 iier = 33
      call uliber2 (iier,
     1             ' supmap-attempt to use non-existent projection')
      go to 320

c cylindrical equidistant,  arbitrary pole and orientation.
  180 if (abs(1.-cosa*cosa) .lt. 1.e-4) go to 320
      u = atan2(sinb*cosr+cosb*sinr,sinb*sinr-cosb*cosr)*rtd
      v = 90.-acos(cosa)*rtd
      go to 305

c mercator, arbitrary pole and orientation.
  190 if ((1.-cosa*cosa) .lt. 2.e-6) go to 320
      u = atan2(sinb*cosr+cosb*sinr,sinb*sinr-cosb*cosr)
      v = alog((1.+cosa)*ovsina)
      go to 305

c mollweide, arbitrary pole and orientation.
  200 if (abs(1.-cosa*cosa) .lt. 2.e-6) go to 320
      u = atan2(sinb*cosr+cosb*sinr,sinb*sinr-cosb*cosr)*tovpi
      udif = abs(u-oldu)
      oldu = u
      v = cosa
c     u = u*sqrt(1.-v*v)
      u = u*sqrt(abs(1.-v*v))
      go to 310

c cylindrical equidistant for polat = rot = 0.  [11]
  210 v = rlat
      go to 305

c mercator                              [12]
  220 u = u*dtr
      v = alog(tan(0.00872664*(rlat+90.0001)))
      go to 305

c mollweide                             [13]
  230 u = u*ov90
      v = sin(rlat*dtr)
      udif = abs(u-oldu)
      oldu = u
      u = u*sqrt(1.-v*v)
      go to 310

c terminal phase    [1,2,4,5,6]
  300 u = r*(sinb*cosr+cosb*sinr)
      v = r*(cosb*cosr-sinb*sinr)

c check for crossover
  305 udif = abs(u-oldu)
      oldu = u
  310 vdif = abs(v-oldv)
      oldv = v
      icross = 0
      if (udif.gt.ueps .or. vdif.gt.veps) icross = 1
      return

c dispense with undefined points
  320 u = 1.e12
      icross = 0
      if (abs(u-oldu) .gt. ueps) icross = 1
      oldu = u
      return

      end
c-------------------------------------------------------------------------------
      subroutine qvec
c this subroutine transforms and plots line segments for supmap and others

c inputs (passed through common.)

c  (rlat,rlon)  next point to be plotted
c  ifst         a flag used to signal the first point of a line segment
c               = 0  -  start a new line
c               = 1  -  continuation of a line

c other variables

c  (u,v)        next point transformed to the virtual screen by supconq
c  icross       a flag returned by supconq for cylindrical projections
c  igo          = 0  -  last point not plotted
c               = 1  -  last point was plotted.
c  (u1,v1),(u2,v2)  parameters passed to suptrp.

        common/supmp4/ifst,igo,igold,icross,iout,uold,vold
        common/supmp6/umin,umax,vmin,vmax,ueps,veps
        common/supmp8/u,v,u1,v1,u2,v2
        common/supmp9/ds,di,dsrdi

c     squ(x) = (x)*(x)

c transform the point
      call qcon

c have we flipped to other side of frame?
      if (icross .ne. 0) igo = 0

c are we within the frame?
      if (u.gt.umax .or. u.lt.umin .or. v.gt.vmax .or. v.lt.vmin)
     1    go to 20
      if (igo .eq. 0) go to 30

c continue line
c check proximity to previous point.
c   5 if ((squ(u-uold)+squ(v-vold))*ds .le. di) return
    5 continue
      call vecplt
   10 uold = u
      vold = v
      igold = igo
      return

c this point lies outside the frame
   20 igo = 0
      if (ifst .ne. 0) go to 65

      if (igold .eq. 0) go to 10

c it was inside - interpolate to edge of frame
c status of last point.   if not inside frame, go on

c if uninterpolatable
      if (icross .ne. 0) go to 70

      u1 = uold
      v1 = vold
      u2 = u
      v2 = v
      call suptrp

c check proximity to previous point.
c     if ((squ(u-uold)+squ(v-vold))*ds .le. di) go to 25
      call vecplt
   25 uold = u2
      vold = v2
      igold = 0
      return

c this point is within the frame

c is it the first point of a line?
   30 if (ifst .ne. 0) go to 60
      if (igold .eq. 0) go to 50

c the previous point was inside the frame on the other side.
c start a new line
   40 call frstpt (u,v)
      igo = 1
      go to 10

c last point not in frame - this one is
   50 if (icross .ne. 0) go to 40

c interpolate back to edge
      u1 = u
      v1 = v
      u2 = uold
      v2 = vold
      call suptrp
      call frstpt (u,v)
      igo = 1
      igold = 1
      uold = u
      vold = v
      u = u1
      v = v1
      go to 5

c first point on line segment  -  check for duplication of end point
   60 if (u.ne.uold .or. v.ne.vold) call frstpt (u,v)
      igo = 1
   65 ifst = 0
      go to 10

c ignore undefined point
   70 ifst = 1
      go to 10

      end
c-------------------------------------------------------------------------------
      subroutine suptrp
c the interpolation routine.  finds (u,v) on the edge of the frame nearest
c (u1,v1).  (u1,v1) must lie within the frame, (u2,v2) without.

        common/supmp6/umin,umax,vmin,vmax,ueps,veps
        common/supmp8/u,v,u1,v1,u2,v2

      f(v) = (v-v2)*du/dv+u2
      g(u) = (u-u2)*dv/du+v2

c find index to (u2,v2)
c       5 | 4 | 6
c       --+---+--
c       2 | 1 | 3
c       --+---+--
c       8 | 7 | 9
      i = 1
      du = u1-u2
      dv = v1-v2
      a = u2-umin
      b = u2-umax
      c = v2-vmin
      d = v2-vmax
      if (a) 110,140,120
  110 i = i+1
      go to 140
  120 if (b) 140,140,130
  130 i = i+2
  140 if (c) 150,200,160
  150 i = i+6
      go to 200
  160 if (d) 200,200,170
  170 i = i+3

  200 go to (900,210,220,230,240,250,260,270,280),i

  210 u = umin
      go to 300
  220 u = umax
      go to 300
  230 v = vmax
      go to 350
  240 if (f(vmax)-umin) 210,230,230
  250 if (f(vmax)-umax) 230,230,220
  260 v = vmin
      go to 350
  270 if (f(vmin)-umin) 210,260,260
  280 if (f(vmin)-umax) 260,260,220

c interpolate
  300 v = g(u)
      return
  350 u = f(v)
      return

c error exit
  900 u = u2
      v = v2
      return
      end
c-------------------------------------------------------------------------------
      subroutine vecplt
c plots the line segment from (uold,vold) to (u,v)
c inputs (passed through common)
c  (uold,vold)  the last point plotted
c  (u,v)        the next point
c  idot         control flag  [ dot vs plot ]

        common/supmp4/ifst,igo,igold,icross,iout,uold,vold
        common/supmp7/phio,phia,igrid,idot,ilts
        common/supmp8/u,v,u1,v1,u2,v2
        common/supmp9/ds,di,dsrdi

c do we dot or plot?
      if (idot .ne. 0) go to 10

c plot
      call vector (u,v)
      return

c dot
   10 du = u-uold
      dv = v-vold

      call line(uold,vold,u,v)
      return

      i = (abs(du)+abs(dv))*dsrdi
      if (i .le. 1) go to 30
      a = 1./float(i)
      i = i-1
      du = du*a
      dv = dv*a
      uo = u
      vo = v
      u = uold
      v = vold
      do 20 k=1,i
         u = u+du
         v = v+dv
         call point (u,v)
   20 continue
      u = uo
      v = vo
   30 call point (u,v)
      return
      end
c-------------------------------------------------------------------------------
      subroutine supvec (xlat,xlon)
c this subroutine allows the user to draw lines on the virtual screen set up by
c supmap, unencumbered by decisions as to whether it will be visible through
c the window.
c use supfst and supvec in exactly the same manner as frstpt and vector.

        common/supmp3/polong,cone,rlat,rlon,jgr,ilf,sgn

      rlat = xlat
      rlon = xlon
      call qvec
      return
      end
c-------------------------------------------------------------------------------
      subroutine supfst (xlat,xlon)

        common/supmp3/polong,cone,rlat,rlon,jgr,ilf,sgn
        common/supmp4/ifst,igo,igold,icross,iout,uold,vold

      igo = 0
      ifst = 1
      rlat = xlat
      rlon = xlon
      call qvec
      return

      end
c-------------------------------------------------------------------------------
      subroutine supcon (xlat,xlon,xu,xv)
c this subroutine is provided to retain compatibility with user's programs.

        common/supmp3/polong,cone,rlat,rlon,jgr,ilf,sgn
        common/supmp8/u,v,u1,v1,u2,v2

      rlat = xlat
      rlon = xlon
      call qcon
      xu = u
      xv = v

      return
      end
c-------------------------------------------------------------------------------
!       moved to laps library (steve albers 1997)
!       block data
c     subroutine supmbd
c force-load block data
!       common/supmp1/pi,tovpi,dtr,rtd,eps,ov90,con1,con2,part
!       common/supmp2/npts,maxlat,minlat,maxlon,minlon,pts(1000)
!       common/supmp3/polong,cone,rlat,rlon,jgr,ilf,sgn
!       common/supmp4/ifst,igo,igold,icross,iout,uold,vold
!       common/supmp5/phioc,sino,coso,sinr,cosr,iproj
!       common/supmp6/umin,umax,vmin,vmax,ueps,veps
!       common/supmp7/phio,phia,igrid,idot,ilts
!       common/supmp8/u,v,u1,v1,u2,v2
!       common/supmpa/iier
!       common/mapcol/mpcol1,mpcol2,mpcol3,mpcol4
!       common/mapdas/ldash1,ldash2,ldash3,ldash4
!       data                    !default line intensities and dash patterns
!    1          mpcol1,ldash1   /255,'1777'o/,  !map lines
!    2          mpcol2,ldash2   /128,'1756'o/,  !grid lines
!    3          mpcol3,ldash3   /192,'1777'o/,  !limb lines
!    4          mpcol4,ldash4   /255,'1777'o/   !perimeter


c       common/supmp1/pi,tovpi,dtr,rtd,eps,ov90,con1,con2,part
c       common/supmp4/ifst,igo,igold,icross,iout,uold,vold
!       common/supmp9/ds,di,dsrdi

!     data   con1 / 1.00001/
!     data   con2 / 179.99999/
!     data     di / 16./
!     data    dtr / 1.7453292519943e-2/
!     data    eps / 1.e-6/
!     data   ov90 / 1.11111111111111e-2/
!     data     pi / 3.1415926535898/
!     data    rtd / 57.295779513082/
!     data  tovpi / 0.63661977236758/
!     data   uold / 0.0 /
!     data   vold / 0.0 /
!     data part/1.0/          !size of picture (90% of screen)

c       return
!     end
c-------------------------------------------------------------------------------
        subroutine resup(ier)
c recall supmap using contents of common blocks.

        common/supmp1/pi,tovpi,dtr,rtd,eps,ov90,con1,con2,part
        common/supmp3/polong,cone,rlat,rlon,jgr,ilf,sgn
        common/supmp6/umin,umax,vmin,vmax,ueps,veps
        common/supmp7/phio,phia,igrid,idot,ilts
        common/supmpb/x1,y1,x2,y2
        common/supmpc/lproj,rot,jlts,lgrid,ius

        lprj=lproj
        pha=phia
        plong=polong
        rt=rot
        umn=umin
        umx=umax
        vmn=vmin
        vmx=vmax
        lgrd=lgrid
        iu=ius
        idt=idot
        part=amax1(x2-x1,y2-y1)

        call supmap(lprj,pha,plong,rt,umn,umx,vmn,vmx,-3,lgrd,iu,idt,ier
     1)

        return

        end
c
c
c
        subroutine supset (jproj, polat, polon, rot,
     +                     pl1, pl2, pl3, pl4, jlts,
     +                     jgrid, jout, jdot, jerr)
c
c
        call supmap (jproj, polat, polon, rot,
     +               pl1, pl2, pl3, pl4, jlts,
     +               0, 0, 0, jerr)
c
c
        return
        end



c!!!!!!
c*******************************************************************************
c
c       the following routines are replacements for ncar graphics routines
c       that are called by supmap.  they are included to supercede the
c       original routines.
c
c*******************************************************************************




        subroutine ichar_to_str (chars, nchars, str)

        integer       chars(*), nchars
        character       str*(*)


        integer       ic, ib, iw

        integer       iword, ichr
        character       ibyte(4)
        equivalence     (ibyte(1), iword)


        iw = 0
        ib = 4
        ic = 0
        do while (ic .lt. nchars)

            ib = ib + 1
            if (ib .gt. 4) then
                iw    = iw + 1
                iword = chars(iw)
                ib    = 1
            end if

            ic       = ic + 1
            ichr     = byte_to_i4(ibyte(ib))
            str(ic:) = char (ichr)

        end do

        return
        end




        subroutine dashln (idash_pat)

        integer       idash_pat


        return
        end




        subroutine optn (iopt, ival)

        integer       iopt, ival


        return
        end




        subroutine pwrt (u, v, chars, nchars, isize, iorient)

        real          u, v
        integer       chars(*), nchars, isize, iorient


        character       str*256


        call ichar_to_str (chars, nchars, str)

        write (*,'(2f10.5,3i5,5x,a)') u, v, isize, iorient, nchars, str

        return
        end




c       subroutine set (x1, x2, y1, y2, u1, u2, v1, v2, itype)
c
c       real          x1, x2, y1, y2, u1, u2, v1, v2
c       integer       itype
c
c
c       return
c       end




        subroutine uliber2 (ier, chars)

        integer       ier

        character       chars*(*)

        write(6,*)' supmap error: ',ier,chars

        return
        end


        subroutine supmap_block_data()

c     routine supmap.f

        common/supmp1/pi,tovpi,dtr,rtd,eps,ov90,con1,con2,part
        common/supmp2/npts,maxlat,minlat,maxlon,minlon,pts(1200000)
        common/supmp3/polong,cone,rlat,rlon,jgr,ilf,sgn
        common/supmp4/ifst,igo,igold,icross,iout,uold,vold
        common/supmp5/phioc,sino,coso,sinr,cosr,iproj
        common/supmp6/umin,umax,vmin,vmax,ueps,veps
        common/supmp7/phio,phia,igrid,idot,ilts
        common/supmp8/u,v,u1,v1,u2,v2
        common/supmpa/iier
        common/mapcol/mpcol1,mpcol2,mpcol3,mpcol4
        common/mapdas/ldash1,ldash2,ldash3,ldash4

!       data                    !default line intensities and dash patterns
!    1          mpcol1,ldash1   /255,1023/,  !map lines
!    2          mpcol2,ldash2   /128,1006/,  !grid lines
!    3          mpcol3,ldash3   /192,1023/,  !limb lines
!    4          mpcol4,ldash4   /255,1023/   !perimeter


c       common/supmp1/pi,tovpi,dtr,rtd,eps,ov90,con1,con2,part
c       common/supmp4/ifst,igo,igold,icross,iout,uold,vold
        common/supmp9/ds,di,dsrdi

        mpcol1 = 255
        mpcol2 = 128
        mpcol3 = 192
        mpcol4 = 255

        ldash1 = 1023 
        ldash2 = 1006
        ldash3 = 1023
        ldash4 = 1023

        con1 = 1.00001
        con2 = 179.99999
        di = 16.
        dtr = 1.7453292519943e-2
        eps = 1.e-6
        ov90 = 1.11111111111111e-2
        pi = 3.1415926535898
        rtd = 57.295779513082
        tovpi = 0.63661977236758
        uold = 0.0 
        vold = 0.0 
        part   = 1.0           !size of picture (90% of screen)
  

        return
        end
