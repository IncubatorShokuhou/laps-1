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
cc
c attributes:
c   language: silicongraphics 3.5 fortran 77
c   machine:  silicongraphics iris-4d/25, 35, indigo, indy
c   language: ibm vs fortran, cray cft77 fortran
c   machine:  hds, cray c916-128, y-mp8/864, cray y-mp el92/256
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
c     integer       grd50 (18)
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
      integer       grd126(18)
      integer       grd190(18)
      integer       grd192(18)
      integer       grd196(18)
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
      data  grd28 / 0, 255, 5,  65, 65,  20826,  145000,   8,  100000,
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
      data  grd90 / 0, 255,201,12902,1,    182, -149887, 136,      92,
     &     141,    577,538,64, 0, 0, 0, 0, 0/
      data  grd91 / 0, 255,202,25803,1,    182, -149887, 136,     183,
     &     141,    577,538,64, 0, 0, 0, 0, 0/
      data  grd92 / 0, 255,201,27071,3,    407, -144094, 136,     223,
     &     365,    222,205,64, 0, 0, 0, 0, 0/
      data  grd93 / 0, 255,202,32485,5,    407, -144094, 136,     445,
     &     365,    222,205,64, 0, 0, 0, 0, 0/
      data  grd94 / 0, 255,201,48916,1,   9678, -128826, 136,     181,
     &     271,    194,185,64, 0, 0, 0, 0, 0/
      data  grd95 / 0, 255,202,97831,1,   9678, -128826, 136,     361,
     &     271,    194,185,64, 0, 0, 0, 0, 0/
      data  grd96 / 0, 255,201,41630,1,  -3441, -148799, 136,     160,
     &     261,    333,308,64, 0, 0, 0, 0, 0/
      data  grd97 / 0, 255,202,83259,1,  -3441, -148799, 136,     319,
     &     261,    333,308,64, 0, 0, 0, 0, 0/
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
      data  grd126/ 0, 255, 4, 384,190,  89277,       0, 128,  -89277,
     &    -938,    95, 938, 0, 0, 0, 0, 0, 0/
      data  grd190 / 0, 255,203, 92,141,    182, -149887, 136,   52000,
     & -111000,    577,538,64, 0, 0, 0, 0, 0/
      data  grd192 / 0, 255,203,223,365,    407, -144094, 136,   50000,
     & -107000,    222,205,64, 0, 0, 0, 0, 0/
      data  grd196/ 0, 255,201,45903,1,  23476,  -96745, 136,     151,
     &     305,     67, 66, 64, 0, 0, 0, 0, 0/
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
      data  grd209/ 0, 255, 3, 101, 81,  22289, -117991,   8,  -95000,
     &   40635,  40635, 0, 64, 0, 25000, 25000, 0, 0/
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
      data  grd217/ 0, 255, 5, 289,205,  42085, -175641,   8, -150000,
     &   15875,  15875, 0, 64, 0, 0, 0, 0, 0/
      data  grd218/ 0, 255, 3, 737,513,  12190, -133459,   8,  -95000,
     &   10159,  10159, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd219/ 0, 255, 5, 385,465,  25008, -119559,  72,  -80000,
     &   25400,  25400, 0, 64, 0, 0, 0, 0, 0/
      data  grd220/ 0, 255, 5, 345,355, -36889, -220194,  72, -260000,
     &   25400,  25400, 1, 64, 0, 0, 0, 0, 0/
      data  grd221/ 0, 255, 3, 349,277,   1000, -145500,   8, -107000,
     &   32463,  32463, 0, 64, 0, 50000, 50000, 0, 0/
      data  grd222/ 0, 255, 3,  59, 47,   1000, -145500,   8, -107000,
     &  194780, 194780, 0, 64, 0, 50000, 50000, 0, 0/
      data  grd223/ 0, 255, 5, 129,129, -20826, -150000,   8, -105000,
     &  190500, 190500, 0, 64, 0, 0, 0, 0, 0/
      data  grd224/ 0, 255, 5,  65, 65,  20826,  120000,   8, -105000,
     &  381000, 381000, 0, 64, 0, 0, 0, 0, 0/
      data  grd225/ 0, 255, 1, 185,135, -25000, -250000, 128,   60640,
     & -250871, 80000, 80000, 20000, 64, 0, 0, 0, 0/
      data  grd226/ 0, 255, 3, 737,513,  12190, -133459,   8,  -95000,
     &   10159,  10159, 0, 64, 0, 25000, 25000, 0, 0/
      data  grd227/ 0, 255, 3,1473,1025,  12190, -133459,   8,  -95000,
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
        do 3 i = 1,14
          igds(i) = grd1(i)
  3     continue
c
      else if (igrid.eq.2) then
        do 4 i = 1,14
          igds(i) = grd2(i)
  4     continue
c
      else if (igrid.eq.3) then
        do 5 i = 1,14
          igds(i) = grd3(i)
  5     continue
c
      else if (igrid.eq.4) then
        do 6 i = 1,14
          igds(i) = grd4(i)
  6     continue
c
      else if (igrid.eq.5) then
        do 10 i = 1,14
          igds(i) = grd5(i)
 10     continue
c
      else if (igrid.eq.6) then
        do 20 i = 1,14
          igds(i) = grd6(i)
 20     continue
c
      else if (igrid.eq.8) then
        do i = 1,14
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
        do 90 i = 1,14
          igds(i) = grd27(i)
 90     continue
c
      else if (igrid.eq.28) then
        do 100 i = 1,14
          igds(i) = grd28(i)
 100    continue
c
      else if (igrid.eq.29) then
        do 110 i = 1,14
          igds(i) = grd29(i)
 110    continue
c
      else if (igrid.eq.30) then
        do 120 i = 1,14
         igds(i) = grd30(i)
 120    continue
c
      else if (igrid.eq.33) then
        do 130 i = 1,14
          igds(i) = grd33(i)
 130     continue
c
      else if (igrid.eq.34) then
        do 140 i = 1,14
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
        do 149 i = 1,14
          igds(i) = grd45(i)
 149    continue
c
c     else if (igrid.eq.50) then
c       do 150 i = 1,14
c         igds(i) = grd50(i)
c150    continue
c
      else if (igrid.eq.53) then
        do i = 1,14
          igds(i) = grd53(i)
        end do
c
      else if (igrid.eq.55) then
        do 152 i = 1,14
          igds(i) = grd55(i)
 152    continue
c
      else if (igrid.eq.56) then
        do 154 i = 1,14
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
        do 192 i = 1,14
          igds(i) = grd85(i)
 192    continue
c
      else if (igrid.eq.86) then
        do 194 i = 1,14
          igds(i) = grd86(i)
 194    continue
c
      else if (igrid.eq.87) then
        do 195 i = 1,14
          igds(i) = grd87(i)
 195    continue
c
      else if (igrid.eq.90) then
        do 196 i = 1,14
          igds(i) = grd90(i)
 196    continue
c
      else if (igrid.eq.91) then
        do 197 i = 1,14
          igds(i) = grd91(i)
 197    continue
c
      else if (igrid.eq.92) then
        do 198 i = 1,14
          igds(i) = grd92(i)
 198    continue
c
      else if (igrid.eq.93) then
        do 199 i = 1,14
          igds(i) = grd93(i)
 199    continue
c
      else if (igrid.eq.94) then
        do 200 i = 1,14
          igds(i) = grd94(i)
 200    continue
c
      else if (igrid.eq.95) then
        do 201 i = 1,14
          igds(i) = grd95(i)
 201    continue
c
      else if (igrid.eq.96) then
        do 202 i = 1,14
          igds(i) = grd96(i)
 202    continue
c
      else if (igrid.eq.97) then
        do 203 i = 1,14
          igds(i) = grd97(i)
 203    continue
c
      else if (igrid.eq.98) then
        do 204 i = 1,14
          igds(i) = grd98(i)
 204    continue
c
      else if (igrid.eq.100) then
        do 205 i = 1,14
          igds(i) = grd100(i)
 205    continue
c
      else if (igrid.eq.101) then
        do 210 i = 1,14
          igds(i) = grd101(i)
 210    continue
c
      else if (igrid.eq.103) then
        do 220 i = 1,14
          igds(i) = grd103(i)
 220   continue
c
      else if (igrid.eq.104) then
        do 230 i = 1,14
          igds(i) = grd104(i)
 230    continue
c
      else if (igrid.eq.105) then
        do 240 i = 1,14
          igds(i) = grd105(i)
 240    continue
c
      else if (igrid.eq.106) then
        do 242 i = 1,14
          igds(i) = grd106(i)
 242    continue
c
      else if (igrid.eq.107) then
        do 244 i = 1,14
          igds(i) = grd107(i)
 244    continue
c
      else if (igrid.eq.126) then
        do 245 i = 1,14
          igds(i) = grd126(i)
 245    continue
c
      else if (igrid.eq.190) then
        do 2190 i = 1,14
          igds(i) = grd190(i)
 2190   continue
c
      else if (igrid.eq.192) then
        do 2192 i = 1,14
          igds(i) = grd192(i)
 2192   continue
c
      else if (igrid.eq.196) then
        do 249 i = 1,14
          igds(i) = grd196(i)
 249    continue
c
      else if (igrid.eq.201) then
        do 250 i = 1,14
          igds(i) = grd201(i)
 250    continue
c
      else if (igrid.eq.202) then
        do 260 i = 1,14
          igds(i) = grd202(i)
 260    continue
c
      else if (igrid.eq.203) then
        do 270 i = 1,14
          igds(i) = grd203(i)
 270    continue
c
      else if (igrid.eq.204) then
        do 280 i = 1,14
          igds(i) = grd204(i)
 280    continue
c
      else if (igrid.eq.205) then
        do 290 i = 1,14
          igds(i) = grd205(i)
 290    continue
c
      else if (igrid.eq.206) then
        do 300 i = 1,18
          igds(i) = grd206(i)
 300    continue
c
      else if (igrid.eq.207) then
        do 310 i = 1,14
          igds(i) = grd207(i)
 310    continue
c
      else if (igrid.eq.208) then
        do 320 i = 1,14
          igds(i) = grd208(i)
 320    continue
c
      else if (igrid.eq.209) then
        do 330 i = 1,18
          igds(i) = grd209(i)
 330    continue
c
      else if (igrid.eq.210) then
        do 340 i = 1,14
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
        do 370 i = 1,14
          igds(i) = grd213(i)
 370    continue
c
      else if (igrid.eq.214) then
        do 380 i = 1,14
          igds(i) = grd214(i)
 380    continue
c
      else if (igrid.eq.215) then
        do 390 i = 1,18
          igds(i) = grd215(i)
 390    continue
c
      else if (igrid.eq.216) then
        do 400 i = 1,14
          igds(i) = grd216(i)
 400    continue
c
      else if (igrid.eq.217) then
        do 401 i = 1,14
          igds(i) = grd217(i)
 401    continue
c
      else if (igrid.eq.218) then
        do 410 i = 1,18
          igds(i) = grd218(i)
 410    continue
c
      else if (igrid.eq.219) then
        do 411 i = 1,14
          igds(i) = grd219(i)
 411    continue
c
      else if (igrid.eq.220) then
        do 412 i = 1,14
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
        do 415 i = 1,14
          igds(i) = grd223(i)
 415    continue
c
      else if (igrid.eq.224) then
        do 416 i = 1,14
          igds(i) = grd224(i)
 416    continue
c
      else if (igrid.eq.225) then
        do 417 i = 1,14
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
        do 420 i = 1,14
          igds(i) = grd228(i)
 420    continue
c
      else if (igrid.eq.229) then
        do 421 i = 1,14
          igds(i) = grd229(i)
 421    continue
c
      else if (igrid.eq.230) then
        do 422 i = 1,14
          igds(i) = grd230(i)
 422    continue
c
      else if (igrid.eq.231) then
        do 423 i = 1,14
          igds(i) = grd231(i)
 423    continue
c
      else if (igrid.eq.232) then
        do 424 i = 1,14
          igds(i) = grd232(i)
 424    continue
c
      else if (igrid.eq.233) then
        do 425 i = 1,14
          igds(i) = grd233(i)
 425    continue
c
c
      else if (igrid.eq.234) then
        do 426 i = 1,14
          igds(i) = grd234(i)
 426    continue
c
      else if (igrid.eq.235) then
        do 427 i = 1,14
          igds(i) = grd235(i)
 427    continue
c
      else
        ierr = 1
      endif
c
      return
      end

