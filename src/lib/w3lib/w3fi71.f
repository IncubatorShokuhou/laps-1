      subroutine w3fi71 (igrid, igds, ierr)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    w3fi71      make array used by grib packer for gds
c   prgmmr: r.e.jones        org: w/nmc42    date: 93-03-26
c
c abstract: w3fi71 makes a 18, 37, 55, 64, or 91 word integer array
c     used by w3fi72 grib packer to make the grid description section
c     (gds) - section 2.
c
c program history log:
c   92-02-21  r.e.jones
c   92-07-01  m. farley    added remarks for 'igds' array elements.
c                          added lambert conformal grids and enlarged
c                          idgs array from 14 to 18 words.
c   92-10-03  r.e.jones    added corrections to awips grib tables
c   92-10-16  r.e.jones    add gaussian grid 126 to tables
c   92-10-18  r.e.jones    corrections to lambert conformal tables
c                          and other tables
c   92-10-19  r.e.jones    add gaussian grid  98 to tables
c   93-01-25  r.e.jones    add on84 grids 87, 106, 107 to tables
c   93-03-10  r.e.jones    add on84 grids 1, 55, 56 to tables
c   93-03-26  r.e.jones    add grib grids 2, 3 to tables
c   93-03-29  r.e.jones    add save statement
c   93-06-15  r.e.jones    add grib grids 37 to 44 to tables
c   93-09-29  r.e.jones    gaussian grid document not correct,
c                          w3fi74 will be changed to agree with
c                          it. gaussian grid 98 table has wrong
c                          value.
c   93-10-12  r.e.jones    changes for on388 rev. oct 8,1993 for
c                          grid 204, 208.
c   93-10-13  r.e.jones    correction for grids 37-44, bytes 7-8,
c                          24-25 set to all bits 1 for missing.
c   93-11-23  r.e.jones    add grids 90-93 for eta model
c                          add grid 4 for 720*361 .5 deg. grid
c   94-04-12  r.e.jones    correction for grid 28
c   94-06-01  r.e.jones    add grid 45, 288*145 1.25 deg. grid
c   94-06-22  r.e.jones    add grids 94, 95 for eta model
c   95-04-11  r.e.jones    add grids 96, 97 for eta model
c   95-05-19  r.e.jones    add from 20 km eta model awips grid 215
c   95-10-19  r.e.jones    add from 20 km eta model alaska grid 216
c   95-10-31  iredell      removed saves and prints
c   96-05-08  iredell      correct first latitude for grids 27 and 28
c   96-07-02  r.e.jones    add from 10 km eta model olympic grid 218
c   96-07-02  r.e.jones    add 196 for eta model
c   96-08-15  r.e.jones    add o.n. 84 grid 8 and 53 as grib grid 8
c                          and 53
c   96-11-29  r.e.jones    correction to tables for grid 21-26, 61-64
c   97-01-31  iredell      correct first latitude for grid 30
c   97-10-20  iredell      correct last longitude for grid 98
c   98-07-07  gilbert      add grids 217 and 219 through 235
c   98-09-21  baldwin      add grids 190, 192 for eta model
c   99-01-20  baldwin      add grids 236, 237
c   99-08-18  iredell      add grid 170
c   01-03-08  rogers       changed eta grids 90-97, added eta grids
c                          194, 198. added awips grids 241,242,243,
c                          245, 246, 247, 248, and 250
c   01-03-19  vuong        added awips grids 238,239,240, and 244
c   01-04-02  vuong        correct last longitude for grid 225
c   01-05-03  rogers       added grid 249
c   01-10-10  rogers       redefined 218 for 12-km eta
c                          redefined grid 192 for new 32-km eta grid
c   02-03-27  vuong        added rsas grid 88 and awips grids 251 and 252
c   02-08-06  rogers       redefined grids 90-93,97,194,245-250 for the
c                          8km hi-res-window model and add awips grid 253
c 2003-06-30  gilbert      added grids 145 and 146 for cmaq
c                          and grid 175 for awips over guam.
c 2003-07-08  vuong        corrected latitude for grid 253 and 170, add grid
c                          110, 127, 171 and 172
c 2004-08-05  vuong        corrected latitude for grid 253
c 2004-09-01  gilbert      corrected the orientation and projection center flag
c                          for southern hemisphere grids 28, 172, 220 and 224
c 2004-09-02  vuong        added grids 147, 148, 173 and 254
c 2005-01-04  cooke        added grids 160, 161 and corrected longitude of orientation for grid 172  
c 2005-03-03  vuong        moved grid 170 to grid 174 and add grid 170
c 2005-03-21  vuong        added grids 130
c 2005-09-12  vuong        added grids 163
c
c usage:    call w3fi71 (igrid, igds, ierr)
c   input argument list:
c     igrid       - grib grid number, or office note 84 grid number
c
c   output argument list:
c     igds      - 18, 37, 55, 64, or 91 word integer array with
c                 information to make a grib grid description section.
c     ierr       - 0  correct exit
c                  1  grid type in igrid is not in table
c
c remarks:
c    1) office note grid type 26 is 6 in grib, 26 is an
c       international exchange grid.
c
c    2) values returned in 18, 37, 55, 64, or 91 word integer array
c        igds vary depending on grid representation type.
c
c       lat/lon grid:
c           igds( 1) = number of vertical coordinates
c           igds( 2) = pv, pl or 255
c           igds( 3) = data representation type (code table 6)
c           igds( 4) = no. of points along a latitude
c           igds( 5) = no. of points along a longitude meridian
c           igds( 6) = latitude of origin (south - ive)
c           igds( 7) = longitude of origin (west -ive)
c           igds( 8) = resolution flag (code table 7)
c           igds( 9) = latitude of extreme point (south - ive)
c           igds(10) = longitude of extreme point (west - ive)
c           igds(11) = latitude increment
c           igds(12) = longitude increment
c           igds(13) = scanning mode flags (code table 8)
c           igds(14) = ... through ...
c           igds(18) =   ... not used for this grid
c           igds(19) - igds(91) for grids 37-44, number of points
c                      in each of 73 rows.
c
c       gaussian grid:
c           igds( 1) = ... through ...
c           igds(10) =   ... same as lat/lon grid
c           igds(11) = number of latitude lines between a pole
c                      and the equator
c           igds(12) = longitude increment
c           igds(13) = scanning mode flags (code table 8)
c           igds(14) = ... through ...
c           igds(18) =   ... not used for this grid
c
c       spherical harmonics:
c           igds( 1) = number of vertical coordinates
c           igds( 2) = pv, pl or 255
c           igds( 3) = data representation type (code table 6)
c           igds( 4) = j - pentagonal resolution parameter
c           igds( 5) = k - pentagonal resolution parameter
c           igds( 6) = m - pentagonal resolution parameter
c           igds( 7) = representation type (code table 9)
c           igds( 8) = representation mode (code table 10)
c           igds( 9) = ... through ...
c           igds(18) =   ... not used for this grid
c
c       polar stereographic:
c           igds( 1) = number of vertical coordinates
c           igds( 2) = pv, pl or 255
c           igds( 3) = data representation type (code table 6)
c           igds( 4) = no. of points along x-axis
c           igds( 5) = no. of points along y-axis
c           igds( 6) = latitude of origin (south -ive)
c           igds( 7) = longitute of origin (west -ive)
c           igds( 8) = resolution flag (code table 7)
c           igds( 9) = longitude of meridian parallel to y-axis
c           igds(10) = x-direction grid length (increment)
c           igds(11) = y-direction grid length (increment)
c           igds(12) = projection center flag (0=north pole on plane,
c                                              1=south pole on plane,
c           igds(13) = scanning mode flags (code table 8)
c           igds(14) = ... through ...
c           igds(18) =   .. not used for this grid
c
c       mercator:
c           igds( 1) = ... through ...
c           igds(12) =   ... same as lat/lon grid
c           igds(13) = latitude at which projection cylinder
c                        intersects earth
c           igds(14) = scanning mode flags
c           igds(15) = ... through ...
c           igds(18) =   .. not used for this grid
c
c       lambert conformal:
c           igds( 1) = number of vertical coordinates
c           igds( 2) = pv, pl or 255
c           igds( 3) = data representation type (code table 6)
c           igds( 4) = no. of points along x-axis
c           igds( 5) = no. of points along y-axis
c           igds( 6) = latitude of origin (south -ive)
c           igds( 7) = longitute of origin (west -ive)
c           igds( 8) = resolution flag (code table 7)
c           igds( 9) = longitude of meridian parallel to y-axis
c           igds(10) = x-direction grid length (increment)
c           igds(11) = y-direction grid length (increment)
c           igds(12) = projection center flag (0=north pole on plane,
c                                              1=south pole on plane,
c           igds(13) = scanning mode flags (code table 8)
c           igds(14) = not used
c           igds(15) = first latitude from the pole at which the
c                      secant cone cuts the sperical earth
c           igds(16) = second latitude ...
c           igds(17) = latitude of south pole (millidegrees)
c           igds(18) = longitude of south pole (millidegrees)
c
c       arakawa semi-staggered e-grid on rotated lat/lon grid
c           igds( 1) = number of vertical coordinates
c           igds( 2) = pv, pl or 255
c           igds( 3) = data representation type (code table 6) [201]
c           igds( 4) = ni  - total number of actual data points
c                            included on grid
c           igds( 5) = nj  - dummy second dimension; set=1
c           igds( 6) = la1 - latitude  of first grid point
c           igds( 7) = lo1 - longitude of first grid point
c           igds( 8) = resolution and component flag (code table 7)
c           igds( 9) = la2 - number of mass points along
c                            southernmost row of grid
c           igds(10) = lo2 - number of rows in each column
c           igds(11) = di  - longitudinal direction increment
c           igds(12) = dj  - latitudinal  direction increment
c           igds(13) = scanning mode flags (code table 8)
c           igds(14) = ... through ...
c           igds(18) = ... not used for this grid (set to zero)
c
c       arakawa filled e-grid on rotated lat/lon grid
c           igds( 1) = number of vertical coordinates
c           igds( 2) = pv, pl or 255
c           igds( 3) = data representation type (code table 6) [202]
c           igds( 4) = ni  - total number of actual data points
c                            included on grid
c           igds( 5) = nj  - dummy second dimention; set=1
c           igds( 6) = la1 - latitude latitude of first grid point
c           igds( 7) = lo1 - longitude of first grid point
c           igds( 8) = resolution and component flag (code table 7)
c           igds( 9) = la2 - number of (zonal) points in each row
c           igds(10) = lo2 - number of (meridional) points in each
c                            column
c           igds(11) = di  - longitudinal direction increment
c           igds(12) = dj  - latitudinal  direction increment
c           igds(13) = scanning mode flags (code table 8)
c           igds(14) = ... through ...
c           igds(18) = ... not used for this grid
c
c       arakawa staggered e-grid on rotated lat/lon grid
c           igds( 1) = number of vertical coordinates
c           igds( 2) = pv, pl or 255
c           igds( 3) = data representation type (code table 6) [203]
c           igds( 4) = ni  - number of data points in each row
c           igds( 5) = nj  - number of rows
c           igds( 6) = la1 - latitude of first grid point
c           igds( 7) = lo1 - longitude of first grid point
c           igds( 8) = resolution and component flag (code table 7)
c           igds( 9) = la2 - central latitude
c           igds(10) = lo2 - central longtitude
c           igds(11) = di  - longitudinal direction increment
c           igds(12) = dj  - latitudinal  direction increment
c           igds(13) = scanning mode flags (code table 8)
c           igds(14) = ... through ...
c           igds(18) = ... not used for this grid
c
c   subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 90
c   machine:  ibm sp
c
c$$$
c
      integer       igrid
      integer       igds  (*)
      integer       grd1  (18)
      integer       grd2  (18)
      integer       grd3  (18)
      integer       grd4  (18)
      integer       grd5  (18)
      integer       grd6  (18)
      integer       grd8  (18)
      integer       grd21 (55)
      integer       grd22 (55)
      integer       grd23 (55)
      integer       grd24 (55)
      integer       grd25 (37)
      integer       grd26 (37)
      integer       grd27 (18)
      integer       grd28 (18)
      integer       grd29 (18)
      integer       grd30 (18)
      integer       grd33 (18)
      integer       grd34 (18)
      integer       grd37 (91)
      integer       grd38 (91)
      integer       grd39 (91)
      integer       grd40 (91)
      integer       grd41 (91)
      integer       grd42 (91)
      integer       grd43 (91)
      integer       grd44 (91)
      integer       grd45 (18)
      integer       grd53 (18)
      integer       grd55 (18)
      integer       grd56 (18)
      integer       grd61 (64)
      integer       grd62 (64)
      integer       grd63 (64)
      integer       grd64 (64)
      integer       grd85 (18)
      integer       grd86 (18)
      integer       grd87 (18)
      integer       grd88 (18)
      integer       grd90 (18)
      integer       grd91 (18)
      integer       grd92 (18)
      integer       grd93 (18)
      integer       grd94 (18)
      integer       grd95 (18)
      integer       grd96 (18)
      integer       grd97 (18)
      integer       grd98 (18)
      integer       grd100(18)
      integer       grd101(18)
      integer       grd103(18)
      integer       grd104(18)
      integer       grd105(18)
      integer       grd106(18)
      integer       grd107(18)
      integer       grd110(18)
      integer       grd126(18)
      integer       grd127(18)
      integer       grd130(18)
      integer       grd145(18)
      integer       grd146(18)
      integer       grd147(18)
      integer       grd148(18)
      integer       grd160(18)
      integer       grd161(18)
      integer       grd163(18)
      integer       grd170(18)
      integer       grd171(18)
      integer       grd172(18)
      integer       grd173(18)
      integer       grd174(18)
      integer       grd175(18)
      integer       grd190(18)
      integer       grd192(18)
      integer       grd194(18)
      integer       grd196(18)
      integer       grd198(18)
      integer       grd201(18)
      integer       grd202(18)
      integer       grd203(18)
      integer       grd204(18)
      integer       grd205(18)
      integer       grd206(18)
      integer       grd207(18)
      integer       grd208(18)
      integer       grd209(18)
      integer       grd210(18)
      integer       grd211(18)
      integer       grd212(18)
      integer       grd213(18)
      integer       grd214(18)
      integer       grd215(18)
      integer       grd216(18)
      integer       grd217(18)
      integer       grd218(18)
      integer       grd219(18)
      integer       grd220(18)
      integer       grd221(18)
      integer       grd222(18)
      integer       grd223(18)
      integer       grd224(18)
      integer       grd225(18)
      integer       grd226(18)
      integer       grd227(18)
      integer       grd228(18)
      integer       grd229(18)
      integer       grd230(18)
      integer       grd231(18)
      integer       grd232(18)
      integer       grd233(18)
      integer       grd234(18)
      integer       grd235(18)
      integer       grd236(18)
      integer       grd237(18)
      integer       grd238(18)
      integer       grd239(18)
      integer       grd240(18)
      integer       grd241(18)
      integer       grd242(18)
      integer       grd243(18)
      integer       grd244(18)
      integer       grd245(18)
      integer       grd246(18)
      integer       grd247(18)
      integer       grd248(18)
      integer       grd249(18)
      integer       grd250(18)
      integer       grd251(18)
      integer       grd252(18)
      integer       grd253(18)
      integer       grd254(18)
c
      data  grd1  / 0, 255, 1,  73, 23, -48090,       0, 128,   48090,
     &       0, 513669,513669, 22500, 64, 0, 0, 0, 0/
      data  grd2  / 0, 255, 0, 144, 73,  90000,       0, 128,  -90000,
     &   -2500,   2500, 2500,  0, 0, 0, 0, 0, 0/
      data  grd3  / 0, 255, 0, 360,181,  90000,       0, 128,  -90000,
     &   -1000,   1000, 1000,  0, 0, 0, 0, 0, 0/
      data  grd4  / 0, 255, 0, 720,361,  90000,       0, 128,  -90000,
     &    -500,    500,  500,  0, 0, 0, 0, 0, 0/
      data  grd5  / 0, 255, 5,  53, 57,   7647, -133443,   8, -105000,
     &  190500, 190500, 0, 64, 0, 0, 0, 0, 0/
      data  grd6  / 0, 255, 5,  53, 45,   7647, -133443,   8, -105000,
     &  190500, 190500, 0, 64, 0, 0, 0, 0, 0/
      data  grd8  / 0, 255, 1, 116, 44, -48670,    3104, 128,   61050,
     &       0, 318830, 318830, 22500, 64, 0, 0, 0, 0/
      data  grd21 / 0,  33, 0,65535,37,      0,       0, 128,   90000,
     &  180000,   2500, 5000, 64, 0, 0, 0, 0, 0,
     & 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
     & 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
     & 37, 37, 37, 37, 37, 37,  1/
      data  grd22 / 0,  33, 0,65535,37,      0, -180000, 128,   90000,
     &       0,   2500, 5000, 64, 0, 0, 0, 0, 0,
     & 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
     & 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
     & 37, 37, 37, 37, 37, 37,  1/
      data  grd23 / 0,  33, 0,65535, 37, -90000,       0, 128,       0,
     &  180000,   2500, 5000, 64, 0, 0, 0, 0, 0,
     &  1, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
     & 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
     & 37, 37, 37, 37, 37, 37, 37/
      data  grd24 / 0,  33, 0,65535, 37, -90000, -180000, 128,       0,
     &       0,   2500, 5000, 64, 0, 0, 0, 0, 0,
     &  1, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
     & 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
     & 37, 37, 37, 37, 37, 37, 37/
      data  grd25 / 0,  33, 0,65535, 19,      0,       0, 128,   90000,
     &  355000,   5000, 5000, 64, 0, 0, 0, 0, 0,
     & 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,
     & 72, 72, 72,  1/
      data  grd26 / 0,  33, 0,65535, 19, -90000,       0, 128,       0,
     &  355000,   5000, 5000, 64, 0, 0, 0, 0, 0,
     &  1, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72,
     & 72, 72, 72, 72/
      data  grd27 / 0, 255, 5,  65, 65, -20826, -125000,   8,  -80000,
     &  381000, 381000, 0, 64, 0, 0, 0, 0, 0/
      data  grd28 / 0, 255, 5,  65, 65,  20826,  145000,   8,  -80000,
     &  381000, 381000,128, 64, 0, 0, 0, 0, 0/
      data  grd29 / 0, 255, 0, 145, 37,      0,       0, 128,   90000,
     &  360000,   2500, 2500, 64, 0, 0, 0, 0, 0/
      data  grd30 / 0, 255, 0, 145, 37,  -90000,      0, 128,       0,
     &  360000,   2500, 2500, 64, 0, 0, 0, 0, 0/
      data  grd33 / 0, 255, 0, 181, 46,      0,       0, 128,   90000,
     &  360000,   2000, 2000, 64, 0, 0, 0, 0, 0/
      data  grd34 / 0, 255, 0, 181, 46, -90000,       0, 128,       0,
     &  360000,   2000, 2000, 64, 0, 0, 0, 0, 0/
      data  grd37 / 0,  33, 0,65535,73,      0,  -30000, 128,   90000,
     &   60000,  1250,65535, 64, 0, 0, 0, 0, 0,
     & 73, 73, 73, 73, 73, 73, 73, 73, 72, 72, 72, 71, 71, 71, 70,
     & 70, 69, 69, 68, 67, 67, 66, 65, 65, 64, 63, 62, 61, 60, 60,
     & 59, 58, 57, 56, 55, 54, 52, 51, 50, 49, 48, 47, 45, 44, 43,
     & 42, 40, 39, 38, 36, 35, 33, 32, 30, 29, 28, 26, 25, 23, 22,
     & 20, 19, 17, 16, 14, 12, 11,  9,  8,  6,  5,  3,  2/
      data  grd38 / 0,  33, 0,65535,73,      0,   60000, 128,   90000,
     &  150000,  1250,65535, 64, 0, 0, 0, 0, 0,
     & 73, 73, 73, 73, 73, 73, 73, 73, 72, 72, 72, 71, 71, 71, 70,
     & 70, 69, 69, 68, 67, 67, 66, 65, 65, 64, 63, 62, 61, 60, 60,
     & 59, 58, 57, 56, 55, 54, 52, 51, 50, 49, 48, 47, 45, 44, 43,
     & 42, 40, 39, 38, 36, 35, 33, 32, 30, 29, 28, 26, 25, 23, 22,
     & 20, 19, 17, 16, 14, 12, 11,  9,  8,  6,  5,  3,  2/
      data  grd39 / 0,  33, 0,65535,73,      0,  150000, 128,   90000,
     & -120000,  1250,65535, 64, 0, 0, 0, 0, 0,
     & 73, 73, 73, 73, 73, 73, 73, 73, 72, 72, 72, 71, 71, 71, 70,
     & 70, 69, 69, 68, 67, 67, 66, 65, 65, 64, 63, 62, 61, 60, 60,
     & 59, 58, 57, 56, 55, 54, 52, 51, 50, 49, 48, 47, 45, 44, 43,
     & 42, 40, 39, 38, 36, 35, 33, 32, 30, 29, 28, 26, 25, 23, 22,
     & 20, 19, 17, 16, 14, 12, 11,  9,  8,  6,  5,  3,  2/
      data  grd40 / 0,  33, 0,65535,73,       0, -120000, 128,   90000,
     &  -30000,  1250,65535, 64, 0, 0, 0, 0, 0,
     & 73, 73, 73, 73, 73, 73, 73, 73, 72, 72, 72, 71, 71, 71, 70,
     & 70, 69, 69, 68, 67, 67, 66, 65, 65, 64, 63, 62, 61, 60, 60,
     & 59, 58, 57, 56, 55, 54, 52, 51, 50, 49, 48, 47, 45, 44, 43,
     & 42, 40, 39, 38, 36, 35, 33, 32, 30, 29, 28, 26, 25, 23, 22,
     & 20, 19, 17, 16, 14, 12, 11,  9,  8,  6,  5,  3,  2/
      data  grd41 / 0,  33, 0,65535,73, -90000,  -30000, 128,       0,
     &   60000,  1250,65535, 64, 0, 0, 0, 0, 0,
     &  2,  3,  5,  6,  8,  9, 11, 12, 14, 16, 17, 19, 20, 22, 23,
     & 25, 26, 28, 29, 30, 32, 33, 35, 36, 38, 39, 40, 42, 43, 44,
     & 45, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, 59, 60, 60,
     & 61, 62, 63, 64, 65, 65, 66, 67, 67, 68, 69, 69, 70, 70, 71,
     & 71, 71, 72, 72, 72, 73, 73, 73, 73, 73, 73, 73, 73/
      data  grd42 / 0,  33, 0,65535,73, -90000,   60000, 128,       0,
     &  150000,  1250,65535, 64, 0, 0, 0, 0, 0,
     &  2,  3,  5,  6,  8,  9, 11, 12, 14, 16, 17, 19, 20, 22, 23,
     & 25, 26, 28, 29, 30, 32, 33, 35, 36, 38, 39, 40, 42, 43, 44,
     & 45, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, 59, 60, 60,
     & 61, 62, 63, 64, 65, 65, 66, 67, 67, 68, 69, 69, 70, 70, 71,
     & 71, 71, 72, 72, 72, 73, 73, 73, 73, 73, 73, 73, 73/
      data  grd43 / 0,  33, 0,65535,73, -90000,  150000, 128,       0,
     & -120000,  1250,65535, 64, 0, 0, 0, 0, 0,
     &  2,  3,  5,  6,  8,  9, 11, 12, 14, 16, 17, 19, 20, 22, 23,
     & 25, 26, 28, 29, 30, 32, 33, 35, 36, 38, 39, 40, 42, 43, 44,
     & 45, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, 59, 60, 60,
     & 61, 62, 63, 64, 65, 65, 66, 67, 67, 68, 69, 69, 70, 70, 71,
     & 71, 71, 72, 72, 72, 73, 73, 73, 73, 73, 73, 73, 73/
      data  grd44 / 0,  33, 0,65535,73, -90000, -120000, 128,       0,
     &  -30000,  1250,65535, 64, 0, 0, 0, 0, 0,
     &  2,  3,  5,  6,  8,  9, 11, 12, 14, 16, 17, 19, 20, 22, 23,
     & 25, 26, 28, 29, 30, 32, 33, 35, 36, 38, 39, 40, 42, 43, 44,
     & 45, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, 59, 60, 60,
     & 61, 62, 63, 64, 65, 65, 66, 67, 67, 68, 69, 69, 70, 70, 71,
     & 71, 71, 72, 72, 72, 73, 73, 73, 73, 73, 73, 73, 73/
      data  grd45 / 0, 255, 0, 288,145,  90000,       0, 128,  -90000,
     &   -1250,   1250, 1250,  0, 0, 0, 0, 0, 0/
      data  grd53 / 0, 255, 1, 117, 51, -61050,       0, 128,   61050,
     &       0,  318830, 318830, 22500, 64, 0, 0, 0, 0/
      data  grd55 / 0, 255, 5,  87, 71, -10947, -154289,   8, -105000,
     &  254000, 254000, 0, 64, 0, 0, 0, 0, 0/
      data  grd56 / 0, 255, 5,  87, 71,   7647, -133443,   8, -105000,
     &  127000, 127000, 0, 64, 0, 0, 0, 0, 0/
      data  grd61 / 0,  33, 0,65535, 46,      0,       0, 128,   90000,
     &  180000,   2000, 2000, 64, 0, 0, 0, 0, 0,
     & 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     & 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     & 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     &  1/
      data  grd62 / 0,  33, 0,65535, 46,      0, -180000, 128,   90000,
     &       0,   2000, 2000, 64, 0, 0, 0, 0, 0,
     & 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     & 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     & 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     &  1/
      data  grd63 / 0,  33, 0,65535, 46,      0,  -90000, 128,       0,
     &  180000,   2000, 2000, 64, 0, 0, 0, 0, 0,
     &  1, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     & 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     & 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     & 91/
      data  grd64 / 0,  33, 0,65535, 46, -90000, -180000, 128,       0,
     &       0,   2000, 2000, 64, 0, 0, 0, 0, 0,
     &  1, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     & 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     & 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91,
     & 91/
      data  grd85 / 0, 255, 0, 360, 90,    500,     500, 128,   89500,
     &  359500,   1000, 1000, 64, 0, 0, 0, 0, 0/
      data  grd86 / 0, 255, 0, 360, 90, -89500,     500, 128,    -500,
     &  359500,   1000, 1000, 64, 0, 0, 0, 0, 0/
      data  grd87 / 0, 255, 5,  81, 62,  22876, -120491,   8, -105000,
     &   68153,  68153, 0, 64, 0, 0, 0, 0, 0/
      data  grd88 / 0, 255, 5, 580,548,  10000, -128000,   8, -105000,
     &   15000,  15000, 0, 64, 0, 0, 0, 0, 0/
      data  grd90 / 0, 255,203,223,501,  23060,  -92570, 136,   37000,
     & -80000,     53,53,64, 0, 0, 0, 0, 0/
      data  grd91 / 0, 255,203,223,501,  23060, -110570, 136,   37000,
     & -98000,     53,53,64, 0, 0, 0, 0, 0/
      data  grd92 / 0, 255,203,223,501,  25986, -127871, 136,   40000,
     & -115000,    53,53,64, 0, 0, 0, 0, 0/
      data  grd93 / 0, 255,203,223,501,  44232, -169996, 136,   63000,
     & -150000,    67,66,64, 0, 0, 0, 0, 0/
      data  grd94 / 0, 255,203,345,569,  -3441, -148799, 136,   50000,
     & -111000,    154,141,64, 0, 0, 0, 0, 0/
      data  grd95 / 0, 255,203,146,247,  35222, -131741, 136,   44000,
     & -240000,     67, 66,64, 0, 0, 0, 0, 0/
      data  grd96 / 0, 255,203,606,1067, -3441, -148799, 136,   50000,
     & -111000,     88,75,64, 0, 0, 0, 0, 0/
      data  grd97 / 0, 255,203, 89,143,  14451,  -71347, 136,   18000,
     &  -66500,     53, 53,64, 0, 0, 0, 0, 0/
      data  grd98 / 0, 255, 4, 192, 94,  88542,       0, 128,  -88542,
     &    -1875, 47,1875, 0, 0, 0, 0, 0, 0/
      data  grd100/ 0, 255, 5,  83, 83,  17108, -129296,   8, -105000,
     &   91452,  91452, 0, 64, 0, 0, 0, 0, 0/
      data  grd101/ 0, 255, 5, 113, 91,  10528, -137146,   8, -105000,
     &   91452,  91452, 0, 64, 0, 0, 0, 0, 0/
      data  grd103/ 0, 255, 5,  65, 56,  22405, -121352,   8, -105000,
     &   91452,  91452, 0, 64, 0, 0, 0, 0, 0/
      data  grd104/ 0, 255, 5, 147,110,   -268, -139475,   8, -105000,
     &   90755,  90755, 0, 64, 0, 0, 0, 0, 0/
      data  grd105/ 0, 255, 5,  83, 83,  17529, -129296,   8, -105000,
     &   90755,  90755, 0, 64, 0, 0, 0, 0, 0/
      data  grd106/ 0, 255, 5, 165,117,  17533, -129296,   8, -105000,
     &   45373,  45373, 0, 64, 0, 0, 0, 0, 0/
      data  grd107/ 0, 255, 5, 120, 92,  23438, -120168,   8, -105000,
     &   45373,  45373, 0, 64, 0, 0, 0, 0, 0/
      data  grd110/ 0, 255, 0, 464,224,  25063, -124938, 128,   52938,
     & -67063,    125, 125, 64, 0, 0, 0, 0, 0/
      data  grd126/ 0, 255, 4, 384,190,  89277,       0, 128,  -89277,
     &    -938,    95, 938, 0, 0, 0, 0, 0, 0/
      data  grd127/ 0, 255, 4, 768,384,  89642,       0, 128,  -89642,
     &    -469,   192, 469, 0, 0, 0, 0, 0, 0/
      data  grd130/ 0, 255, 3, 451,337,  16281, -126138,  8,   -95000,
     &   13545,  13545, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd145/ 0, 255, 3, 169,145,  32174,  -90159,   8,  -79500,
     &   12000,  12000, 0, 64, 0, 36000, 46000, 0, 0/
      data  grd146/ 0, 255, 3, 166,142,  32353,  -89994,   8,  -79500,
     &   12000,  12000, 0, 64, 0, 36000, 46000, 0, 0/
      data  grd147/ 0, 255, 3, 268,259,  24595, -100998,   8,  -97000,
     &   12000,  12000, 0, 64, 0, 33000, 45000, 0, 0/
      data  grd148/ 0, 255, 3, 442,265,  21821, -120628,   8,  -97000,
     &   12000,  12000, 0, 64, 0, 33000, 45000, 0, 0/
      data  grd160/ 0, 255, 5, 180,156,  19132, -185837,   8, -150000,
     &   47500,  47500, 0, 64, 0, 0, 0, 0, 0/
      data  grd161/ 0, 255, 0, 137,102,  50750,  271750,  72,    -250,
     &  -19750,    500,500, 0, 0, 0, 0, 0, 0/ 
      data  grd163/ 0, 255, 3,1008,722,  20600, -118300,   8,  -95000,
     &    5000,   5000, 0, 64, 0, 38000, 38000, 0, 0/
      data  grd170/ 0, 255, 4, 512, 256, 89463,       0, 128,  -89463,
     &    -703,   128, 703, 0, 0, 0, 0, 0, 0/
      data  grd171/ 0, 255, 5, 770,930,  25009, -119560,  72,  -80000,
     &   12700,  12700, 0, 64, 0, 0, 0, 0, 0/
      data  grd172/ 0, 255, 5, 690,710, -36900, -220194,  72, -260000,
     &   12700,  12700, 128, 64, 0, 0, 0, 0, 0/
      data  grd173/ 0, 255, 0,4320,2160, 89958,     417, 128,   89958,
     &  359958,    83,  83, 64, 0, 0, 0, 0, 0/
      data  grd174/ 0, 255, 4,2880,1440, 89938,      62,  72,  -89938,
     &     -62,   125, 125,64, 0, 0, 0, 0, 0/
      data  grd175/ 0, 255, 0, 556,334,      0,  130000, 128,   30060,
     &  180040,    90,  90, 64, 0, 0, 0, 0, 0/
      data  grd190 / 0, 255,203, 92,141,   182, -149887, 136,   52000,
     & -111000,    577,538,64, 0, 0, 0, 0, 0/
      data  grd192 / 0, 255,203,237,387, -3441, -148799, 136,   50000,
     & -111000,    225,207,64, 0, 0, 0, 0, 0/
      data  grd194 / 0, 255,203, 89,143, 16444, -162244, 136,   20250,
     & -157350,     53, 53,64, 0, 0, 0, 0, 0/
      data  grd196/ 0, 255,201,45903,1,  23476,  -96745, 136,     151,
     &     305,     67, 66, 64, 0, 0, 0, 0, 0/
      data  grd198/ 0, 255,203,160,261,  -3441, -148799, 136,   50000,
     & -111000,    333,308,64, 0, 0, 0, 0, 0/
      data  grd201/ 0, 255, 5,  65, 65, -20826, -150000,   8, -105000,
     &  381000, 381000, 0, 64, 0, 0, 0, 0, 0/
      data  grd202/ 0, 255, 5,  65, 43,   7838, -141028,   8, -105000,
     &  190500, 190500, 0, 64, 0, 0, 0, 0, 0/
      data  grd203/ 0, 255, 5,  45, 39,  19132, -185837,   8, -150000,
     &  190500, 190500, 0, 64, 0, 0, 0, 0, 0/
      data  grd204/ 0, 255, 1,  93, 68, -25000,  110000, 128,   60644,
     & -109129, 160000, 160000, 20000, 64, 0, 0, 0, 0/
      data  grd205/ 0, 255, 5,  45, 39,    616,  -84904,   8,  -60000,
     &  190500, 190500, 0, 64, 0, 0, 0, 0, 0/
      data  grd206/ 0, 255, 3,  51, 41,  22289, -117991,   8, - 95000,
     &   81271,  81271, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd207/ 0, 255, 5,  49, 35,  42085, -175641,   8, -150000,
     &   95250,  95250, 0, 64, 0, 0, 0, 0, 0/
      data  grd208/ 0, 255, 1,  29, 27,   9343, -167315, 128,   28092,
     & -145878, 80000, 80000, 20000, 64, 0, 0, 0, 0/
      data  grd209/ 0, 255, 3, 275,223,  -4850, -151100,   8, -111000,
     &   44000,  44000, 0, 64, 0, 45000, 45000, 0, 0/
      data  grd210/ 0, 255, 1,  25, 25,   9000,  -77000, 128,   26422,
     &  -58625, 80000, 80000, 20000, 64, 0, 0, 0, 0/
      data  grd211/ 0, 255, 3,  93, 65,  12190, -133459,   8,  -95000,
     &   81271,  81271, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd212/ 0, 255, 3, 185,129,  12190, -133459,   8,  -95000,
     &   40635,  40635, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd213/ 0, 255, 5, 129, 85,   7838, -141028,   8, -105000,
     &   95250,  95250, 0, 64, 0, 0, 0, 0, 0/
      data  grd214/ 0, 255, 5,  97, 69,  42085, -175641,   8, -150000,
     &   47625,  47625, 0, 64, 0, 0, 0, 0, 0/
      data  grd215/ 0, 255, 3, 369,257,  12190, -133459,   8,  -95000,
     &   20318,  20318, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd216/ 0, 255, 5, 139,107,  30000, -173000,   8, -135000,
     &   45000,  45000, 0, 64, 0, 0, 0, 0, 0/
      data  grd217/ 0, 255, 5, 277,213,  30000, -173000,   8, -135000,
     &   22500,  22500, 0, 64, 0, 0, 0, 0, 0/
      data  grd218/ 0, 255, 3, 614,428,  12190, -133459,   8,  -95000,
     &   12191,  12191, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd219/ 0, 255, 5, 385,465,  25008, -119559,  72,  -80000,
     &   25400,  25400, 0, 64, 0, 0, 0, 0, 0/
      data  grd220/ 0, 255, 5, 345,355, -36889, -220194,  72, -80000,
     &   25400,  25400, 128, 64, 0, 0, 0, 0, 0/
      data  grd221/ 0, 255, 3, 349,277,   1000, -145500,   8, -107000,
     &   32463,  32463, 0, 64, 0, 50000, 50000, 0, 0/
      data  grd222/ 0, 255, 3, 138,112,  -4850, -151100,   8, -111000,
     &   88000,  88000, 0, 64, 0, 45000, 45000, 0, 0/
      data  grd223/ 0, 255, 5, 129,129, -20826, -150000,   8, -105000,
     &  190500, 190500, 0, 64, 0, 0, 0, 0, 0/
      data  grd224/ 0, 255, 5,  65, 65,  20826,  120000,   8, -105000,
     &  381000, 381000, 128, 64, 0, 0, 0, 0, 0/
      data  grd225/ 0, 255, 1, 185,135, -25000, -250000, 128,   60640,
     & -109129, 80000, 80000, 20000, 64, 0, 0, 0, 0/
      data  grd226/ 0, 255, 3, 737,513,  12190, -133459,   8,  -95000,
     &   10159,  10159, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd227/ 0, 255, 3,1473,1025,  12190, -133459,   8, -95000,
     &    5079,   5079, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd228/ 0, 255, 0, 144, 73,  90000,       0, 128,  -90000,
     &   -2500,   2500, 2500, 64, 0, 0, 0, 0, 0/
      data  grd229/ 0, 255, 0, 360,181,  90000,       0, 128,  -90000,
     &   -1000,   1000, 1000, 64, 0, 0, 0, 0, 0/
      data  grd230/ 0, 255, 0, 720,361,  90000,       0, 128,  -90000,
     &    -500,    500,  500, 64, 0, 0, 0, 0, 0/
      data  grd231/ 0, 255, 0, 720,181,      0,       0, 128,   90000,
     &    -500,    500,  500, 64, 0, 0, 0, 0, 0/
      data  grd232/ 0, 255, 0, 360, 91,      0,       0, 128,   90000,
     &   -1000,   1000, 1000, 64, 0, 0, 0, 0, 0/
      data  grd233/ 0, 255, 0, 288,157,  78000,       0, 128,  -78000,
     &   -1250,   1250, 1000, 64, 0, 0, 0, 0, 0/
      data  grd234/ 0, 255, 0, 133,121,  15000,  -98000, 128,  -45000,
     &  -65000,    250,  250, 64, 0, 0, 0, 0, 0/
      data  grd235/ 0, 255, 0, 720,360,  89750,     250,  72,  -89750,
     &    -250,    250, 1000, 64, 0, 0, 0, 0, 0/
      data  grd236/ 0, 255, 3, 151,113,  16281,  233862,   8,  -95000,
     &   40635,  40635, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd237/ 0, 255, 3,  54, 47,  16201,  285720,   8, -107000,
     &   32463,  32463, 0, 64, 0, 50000, 50000, 0, 0/
      data  grd238/ 0, 255, 0, 275, 203,  50750, 261750,  72,    -205,   
     &   -29750, 0,  0, 64, 0, 0, 0, 0, 0/
      data  grd239/ 0, 255, 0, 155, 123, 75750,  159500,  72,   44750, 
     &  -123500,  0, 0, 64, 0, 0, 0, 0, 0/
      data  grd240/ 0, 255, 5, 1121, 881, 23098, -119036,  8, -105000,
     &   47625,  47625, 0, 64, 0, 0, 0, 0, 0/
      data  grd241/ 0, 255, 3, 549,445,  -4850, -151100,   8, -111000,
     &   22000,  22000, 0, 64, 0, 45000, 45000, 0, 0/
      data  grd242/ 0, 255, 5, 553,425,  30000, -173000,   8, -135000,
     &   11250,  11250, 0, 64, 0, 0, 0, 0, 0/
      data  grd243/ 0, 255, 0, 126,101,  10000, -170000, 128,   50000,
     &  -120000, 400, 400, 64, 0, 0, 0, 0, 0/
      data  grd244/ 0, 255, 0, 275, 203,  50750, 261750,  72,    -205,   
     &   -29750, 0,  0, 64, 0, 0, 0, 0, 0/
      data  grd245/ 0, 255, 3, 336,372,  22980, -92840,   8,   -80000,
     &   8000,  8000, 0, 64, 0, 35000, 35000, 0, 0/
      data  grd246/ 0, 255, 3, 332,371,  25970, -127973,  8,  -115000,
     &   8000,  8000, 0, 64, 0, 40000, 40000, 0, 0/
      data  grd247/ 0, 255, 3, 336,372,  22980, -110840,   8,  -98000,
     &   8000,  8000, 0, 64, 0, 35000, 35000, 0, 0/
      data  grd248/ 0, 255, 0, 135,101,  14500,  -71500, 128,   22000,
     &  -61450,    75,  75, 64, 0, 0, 0, 0, 0/
      data  grd249/ 0, 255, 5, 367,343,  45400, -171600,   8, -150000,
     &   9868,  9868, 0, 64, 0, 0, 0, 0, 0/
      data  grd250/ 0, 255, 0, 135,101,  16500, -162000, 128,   24000,
     & -151950,    75,  75, 64, 0, 0, 0, 0, 0/
      data  grd251/ 0, 255, 0, 332,210,  26350,  -83050, 128,   47250,
     &  -49950,    100,  100, 64, 0, 0, 0, 0, 0/
      data  grd252/ 0, 255, 3, 301,225,  16281, -126138,  8,   -95000,
     &   20317,  20317, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd253/ 0, 255, 0, 373,224,  60500, -170250,  72,    4750,
     &  -77250,    75,  75, 64, 0, 0, 0, 0, 0/
      data  grd254/ 0, 255, 0, 369,300, -35000, -250000, 128,   60789,
     & -109129, 20000,   0, 64, 40000, 40000, 0, 0, 0/
c
      ierr = 0
c
        do 1 i = 1,18
          igds(i) = 0
 1      continue
c
      if (igrid.ge.37.and.igrid.le.44) then
        do 2 i = 19,91
          igds(i) = 0
 2      continue
      end if
c
      if (igrid.ge.21.and.igrid.le.24) then
        do i = 19,55
          igds(i) = 0
        end do
      end if
c
      if (igrid.ge.25.and.igrid.le.26) then
        do i = 19,37
          igds(i) = 0
        end do
      end if
c
      if (igrid.ge.61.and.igrid.le.64) then
        do i = 19,64
          igds(i) = 0
        end do
      end if
c
      if (igrid.eq.1) then
        do 3 i = 1,18
          igds(i) = grd1(i)
  3     continue
c
      else if (igrid.eq.2) then
        do 4 i = 1,18
          igds(i) = grd2(i)
  4     continue
c
      else if (igrid.eq.3) then
        do 5 i = 1,18
          igds(i) = grd3(i)
  5     continue
c
      else if (igrid.eq.4) then
        do 6 i = 1,18
          igds(i) = grd4(i)
  6     continue
c
      else if (igrid.eq.5) then
        do 10 i = 1,18
          igds(i) = grd5(i)
 10     continue
c
      else if (igrid.eq.6) then
        do 20 i = 1,18
          igds(i) = grd6(i)
 20     continue
c
      else if (igrid.eq.8) then
        do i = 1,18
          igds(i) = grd8(i)
        end do
c
      else if (igrid.eq.21) then
        do 30 i = 1,55
          igds(i) = grd21(i)
 30     continue
c
      else if (igrid.eq.22) then
        do 40 i = 1,55
          igds(i) = grd22(i)
 40     continue
c
      else if (igrid.eq.23) then
        do 50 i = 1,55
          igds(i) = grd23(i)
 50     continue
c
      else if (igrid.eq.24) then
        do 60 i = 1,55
          igds(i) = grd24(i)
 60     continue
c
      else if (igrid.eq.25) then
        do 70 i = 1,37
          igds(i) = grd25(i)
 70     continue
c
      else if (igrid.eq.26) then
        do 80 i = 1,37
          igds(i) = grd26(i)
 80     continue
c
      else if (igrid.eq.27) then
        do 90 i = 1,18
          igds(i) = grd27(i)
 90     continue
c
      else if (igrid.eq.28) then
        do 100 i = 1,18
          igds(i) = grd28(i)
 100    continue
c
      else if (igrid.eq.29) then
        do 110 i = 1,18
          igds(i) = grd29(i)
 110    continue
c
      else if (igrid.eq.30) then
        do 120 i = 1,18
          igds(i) = grd30(i)
 120    continue
c
      else if (igrid.eq.33) then
        do 130 i = 1,18
          igds(i) = grd33(i)
 130     continue
c
      else if (igrid.eq.34) then
        do 140 i = 1,18
          igds(i) = grd34(i)
 140    continue
c
      else if (igrid.eq.37) then
        do 141 i = 1,91
          igds(i) = grd37(i)
 141    continue
c
      else if (igrid.eq.38) then
        do 142 i = 1,91
          igds(i) = grd38(i)
 142    continue
c
      else if (igrid.eq.39) then
        do 143 i = 1,91
          igds(i) = grd39(i)
 143    continue
c
      else if (igrid.eq.40) then
        do 144 i = 1,91
          igds(i) = grd40(i)
 144    continue
c
      else if (igrid.eq.41) then
        do 145 i = 1,91
          igds(i) = grd41(i)
 145    continue
c
      else if (igrid.eq.42) then
        do 146 i = 1,91
          igds(i) = grd42(i)
 146    continue
c
      else if (igrid.eq.43) then
        do 147 i = 1,91
          igds(i) = grd43(i)
 147    continue
c
      else if (igrid.eq.44) then
        do 148 i = 1,91
          igds(i) = grd44(i)
 148    continue
c
      else if (igrid.eq.45) then
        do 149 i = 1,18
          igds(i) = grd45(i)
 149    continue
c
      else if (igrid.eq.53) then
        do i = 1,18
          igds(i) = grd53(i)
        end do
c
      else if (igrid.eq.55) then
        do 152 i = 1,18
          igds(i) = grd55(i)
 152    continue
c
      else if (igrid.eq.56) then
        do 154 i = 1,18
          igds(i) = grd56(i)
 154    continue
c
      else if (igrid.eq.61) then
        do 160 i = 1,64
          igds(i) = grd61(i)
 160    continue
c
      else if (igrid.eq.62) then
        do 170 i = 1,64
          igds(i) = grd62(i)
 170    continue
c
      else if (igrid.eq.63) then
        do 180 i = 1,64
          igds(i) = grd63(i)
 180    continue
c
      else if (igrid.eq.64) then
        do 190 i = 1,64
          igds(i) = grd64(i)
 190    continue
c
      else if (igrid.eq.85) then
        do 192 i = 1,18
          igds(i) = grd85(i)
 192    continue
c
      else if (igrid.eq.86) then
        do 194 i = 1,18
          igds(i) = grd86(i)
 194    continue
c
      else if (igrid.eq.87) then
        do 195 i = 1,18
          igds(i) = grd87(i)
 195    continue
c
      else if (igrid.eq.88) then
        do 2195 i = 1,18
          igds(i) = grd88(i)
2195    continue
c
      else if (igrid.eq.90) then
        do 196 i = 1,18
          igds(i) = grd90(i)
 196    continue
c
      else if (igrid.eq.91) then
        do 197 i = 1,18
          igds(i) = grd91(i)
 197    continue
c
      else if (igrid.eq.92) then
        do 198 i = 1,18
          igds(i) = grd92(i)
 198    continue
c
      else if (igrid.eq.93) then
        do 199 i = 1,18
          igds(i) = grd93(i)
 199    continue
c
      else if (igrid.eq.94) then
        do 200 i = 1,18
          igds(i) = grd94(i)
 200    continue
c
      else if (igrid.eq.95) then
        do 201 i = 1,18
          igds(i) = grd95(i)
 201    continue
c
      else if (igrid.eq.96) then
        do 202 i = 1,18
          igds(i) = grd96(i)
 202    continue
c
      else if (igrid.eq.97) then
        do 203 i = 1,18
          igds(i) = grd97(i)
 203    continue
c
      else if (igrid.eq.98) then
        do 204 i = 1,18
          igds(i) = grd98(i)
 204    continue
c
      else if (igrid.eq.100) then
        do 205 i = 1,18
          igds(i) = grd100(i)
 205    continue
c
      else if (igrid.eq.101) then
        do 210 i = 1,18
          igds(i) = grd101(i)
 210    continue
c
      else if (igrid.eq.103) then
        do 220 i = 1,18
          igds(i) = grd103(i)
 220   continue
c
      else if (igrid.eq.104) then
        do 230 i = 1,18
          igds(i) = grd104(i)
 230    continue
c
      else if (igrid.eq.105) then
        do 240 i = 1,18
          igds(i) = grd105(i)
 240    continue
c
      else if (igrid.eq.106) then
        do 242 i = 1,18
          igds(i) = grd106(i)
 242    continue
c
      else if (igrid.eq.107) then
        do 244 i = 1,18
          igds(i) = grd107(i)
 244    continue
c
      else if (igrid.eq.110) then
        do i = 1,18
          igds(i) = grd110(i)
        enddo
c
      else if (igrid.eq.126) then
        do 245 i = 1,18
          igds(i) = grd126(i)
 245    continue
c
      else if (igrid.eq.127) then
        do i = 1,18
          igds(i) = grd127(i)
        enddo
c
      else if (igrid.eq.130) then
        do i = 1,18
          igds(i) = grd130(i)
        enddo
c
      else if (igrid.eq.145) then
        do i = 1,18
          igds(i) = grd145(i)
        enddo
c
      else if (igrid.eq.146) then
        do i = 1,18
          igds(i) = grd146(i)
        enddo
c
      else if (igrid.eq.147) then
        do i = 1,18
          igds(i) = grd147(i)
        enddo
c
      else if (igrid.eq.148) then
        do i = 1,18
          igds(i) = grd148(i)
        enddo
c
      else if (igrid.eq.160) then
        do i = 1,18
          igds(i) = grd160(i)
        enddo
c
      else if (igrid.eq.161) then
        do i = 1,18
          igds(i) = grd161(i)
        enddo
      else if (igrid.eq.163) then
        do i = 1,18
          igds(i) = grd163(i)
        enddo
c
      else if (igrid.eq.170) then
        do i = 1,18
          igds(i) = grd170(i)
        enddo
c
      else if (igrid.eq.171) then
        do i = 1,18
          igds(i) = grd171(i)
        enddo
c
      else if (igrid.eq.172) then
        do i = 1,18
          igds(i) = grd172(i)
        enddo
c
      else if (igrid.eq.173) then
        do i = 1,18
          igds(i) = grd173(i)
        enddo
c
      else if (igrid.eq.174) then
        do i = 1,18
          igds(i) = grd174(i)
        enddo
c
      else if (igrid.eq.175) then
        do i = 1,18
          igds(i) = grd175(i)
        enddo
c
      else if (igrid.eq.190) then
        do 2190 i = 1,18
          igds(i) = grd190(i)
 2190   continue
c
      else if (igrid.eq.192) then
        do 2191 i = 1,18
          igds(i) = grd192(i)
 2191   continue
c
      else if (igrid.eq.194) then
        do 2192 i = 1,18
          igds(i) = grd194(i)
 2192   continue
c
      else if (igrid.eq.196) then
        do 249 i = 1,18
          igds(i) = grd196(i)
 249    continue
c
      else if (igrid.eq.198) then
        do 2490 i = 1,18
          igds(i) = grd198(i)
 2490   continue
c
      else if (igrid.eq.201) then
        do 250 i = 1,18
          igds(i) = grd201(i)
 250    continue
c
      else if (igrid.eq.202) then
        do 260 i = 1,18
          igds(i) = grd202(i)
 260    continue
c
      else if (igrid.eq.203) then
        do 270 i = 1,18
          igds(i) = grd203(i)
 270    continue
c
      else if (igrid.eq.204) then
        do 280 i = 1,18
          igds(i) = grd204(i)
 280    continue
c
      else if (igrid.eq.205) then
        do 290 i = 1,18
          igds(i) = grd205(i)
 290    continue
c
      else if (igrid.eq.206) then
        do 300 i = 1,18
          igds(i) = grd206(i)
 300    continue
c
      else if (igrid.eq.207) then
        do 310 i = 1,18
          igds(i) = grd207(i)
 310    continue
c
      else if (igrid.eq.208) then
        do 320 i = 1,18
          igds(i) = grd208(i)
 320    continue
c
      else if (igrid.eq.209) then
        do 330 i = 1,18
          igds(i) = grd209(i)
 330    continue
c
      else if (igrid.eq.210) then
        do 340 i = 1,18
          igds(i) = grd210(i)
 340    continue
c
      else if (igrid.eq.211) then
        do 350 i = 1,18
          igds(i) = grd211(i)
 350    continue
c
      else if (igrid.eq.212) then
        do 360 i = 1,18
          igds(i) = grd212(i)
 360    continue
c
      else if (igrid.eq.213) then
        do 370 i = 1,18
          igds(i) = grd213(i)
 370    continue
c
      else if (igrid.eq.214) then
        do 380 i = 1,18
          igds(i) = grd214(i)
 380    continue
c
      else if (igrid.eq.215) then
        do 390 i = 1,18
          igds(i) = grd215(i)
 390    continue
c
      else if (igrid.eq.216) then
        do 400 i = 1,18
          igds(i) = grd216(i)
 400    continue
c
      else if (igrid.eq.217) then
        do 401 i = 1,18
          igds(i) = grd217(i)
 401    continue
c
      else if (igrid.eq.218) then
        do 410 i = 1,18
          igds(i) = grd218(i)
 410    continue
c
      else if (igrid.eq.219) then
        do 411 i = 1,18
          igds(i) = grd219(i)
 411    continue
c
      else if (igrid.eq.220) then
        do 412 i = 1,18
          igds(i) = grd220(i)
 412    continue
c
      else if (igrid.eq.221) then
        do 413 i = 1,18
          igds(i) = grd221(i)
 413    continue
c
      else if (igrid.eq.222) then
        do 414 i = 1,18
          igds(i) = grd222(i)
 414    continue
c
      else if (igrid.eq.223) then
        do 415 i = 1,18
          igds(i) = grd223(i)
 415    continue
c
      else if (igrid.eq.224) then
        do 416 i = 1,18
          igds(i) = grd224(i)
 416    continue
c
      else if (igrid.eq.225) then
        do 417 i = 1,18
          igds(i) = grd225(i)
 417    continue
c
      else if (igrid.eq.226) then
        do 418 i = 1,18
          igds(i) = grd226(i)
 418    continue
c
      else if (igrid.eq.227) then
        do 419 i = 1,18
          igds(i) = grd227(i)
 419    continue
c
      else if (igrid.eq.228) then
        do 420 i = 1,18
          igds(i) = grd228(i)
 420    continue
c
      else if (igrid.eq.229) then
        do 421 i = 1,18
          igds(i) = grd229(i)
 421    continue
c
      else if (igrid.eq.230) then
        do 422 i = 1,18
          igds(i) = grd230(i)
 422    continue
c
      else if (igrid.eq.231) then
        do 423 i = 1,18
          igds(i) = grd231(i)
 423    continue
c
      else if (igrid.eq.232) then
        do 424 i = 1,18
          igds(i) = grd232(i)
 424    continue
c
      else if (igrid.eq.233) then
        do 425 i = 1,18
          igds(i) = grd233(i)
 425    continue
c
      else if (igrid.eq.234) then
        do 426 i = 1,18
          igds(i) = grd234(i)
 426    continue
c
      else if (igrid.eq.235) then
        do 427 i = 1,18
          igds(i) = grd235(i)
 427    continue
c
      else if (igrid.eq.236) then
        do 428 i = 1,18
          igds(i) = grd236(i)
 428    continue
c
      else if (igrid.eq.237) then
        do 429 i = 1,18
          igds(i) = grd237(i)
 429    continue
c
      else if (igrid.eq.238) then
        do i = 1,18
          igds(i) = grd238(i)
        end do
c
      else if (igrid.eq.239) then
        do i = 1,18
          igds(i) = grd239(i)
        end do
c
      else if (igrid.eq.240) then
        do i = 1,18
          igds(i) = grd240(i)
        end do
c
      else if (igrid.eq.241) then
        do 430 i = 1,18
          igds(i) = grd241(i)
 430    continue
c
      else if (igrid.eq.242) then
        do 431 i = 1,18
          igds(i) = grd242(i)
 431    continue
c
      else if (igrid.eq.243) then
        do 432 i = 1,18
          igds(i) = grd243(i)
 432    continue
c
      else if (igrid.eq.244) then
        do i = 1,18
          igds(i) = grd244(i)
        end do
c
      else if (igrid.eq.245) then
        do 433 i = 1,18
          igds(i) = grd245(i)
 433    continue
c
      else if (igrid.eq.246) then
        do 434 i = 1,18
          igds(i) = grd246(i)
 434    continue
c
      else if (igrid.eq.247) then
        do 435 i = 1,18
          igds(i) = grd247(i)
 435    continue
c
      else if (igrid.eq.248) then
        do 436 i = 1,18
          igds(i) = grd248(i)
 436    continue
c
      else if (igrid.eq.249) then
        do 437 i = 1,18
          igds(i) = grd249(i)
 437    continue
c
      else if (igrid.eq.250) then
        do 438 i = 1,18
          igds(i) = grd250(i)
 438    continue
c
      else if (igrid.eq.251) then
        do 439 i = 1,18
          igds(i) = grd251(i)
 439    continue
c
      else if (igrid.eq.252) then
        do 440 i = 1,18
          igds(i) = grd252(i)
 440    continue
      else if (igrid.eq.253) then
        do 441 i = 1,18
          igds(i) = grd253(i)
 441    continue
      else if (igrid.eq.254) then
        do 442 i = 1,18
          igds(i) = grd254(i)
 442    continue
c
      else
        ierr = 1
      endif
c
      return
      end
