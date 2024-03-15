      subroutine w3fi63(msga,kpds,kgds,kbms,data,kptr,kret)
c$$$  subprogram documentation  block
c                .      .    .                                       .
c subprogram:  w3fi63        unpk grib field to grib grid
c   prgmmr: farley           org: nmc421      date:94-11-22
c
c abstract: unpack a grib (edition 1) field to the exact grid
c   specified in the grib message, isolate the bit map, and make
c   the values of the product descripton section (pds) and the
c   grid description section (gds) available in return arrays.
c
c   when decoding is completed, data at each grid point has been
c          returned in the units specified in the grib manual.
c
c program history log:
c   91-09-13  cavanaugh
c   91-11-12  cavanaugh   modified size of ecmwf grids 5-8
c   91-12-22  cavanaugh   corrected processing of mercator projections
c                         in grid definition section (gds) in
c                         routine fi633
c   92-08-05  cavanaugh   corrected maximum grid size to allow for
c                         one degree by one degree global grids
c   92-08-27  cavanaugh   corrected typo error, added code to compare
c                         total byte size from section 0 with sum of
c                         section sizes.
c   92-10-21  cavanaugh   corrections were made (in fi634) to reduce
c                         processing time for international grids.
c                         removed a typographical error in fi635.
c   93-01-07  cavanaugh   corrections were made (in fi635) to
c                         facilitate use of these routines on a pc.
c                         a typographical error was also corrected
c   93-01-13  cavanaugh   corrections were made (in fi632) to
c                         properly handle condition when
c                         time range indicator = 10.
c                         added u.s.grid 87.
c   93-02-04  cavanaugh   added u.s.grids 85 and 86
c   93-02-26  cavanaugh   added grids 2, 3, 37 thru 44,and
c                         grids 55, 56, 90, 91, 92, and 93 to
c                         list of u.s. grids.
c   93-04-07  cavanaugh   added grids 67 thru 77 to
c                         list of u.s. grids.
c   93-04-20  cavanaugh   increased max size to accomodate
c                         gaussian grids.
c   93-05-26  cavanaugh   corrected grid range selection in fi634
c                         for ranges 67-71 & 75-77
c   93-06-08  cavanaugh   corrected fi635 to accept grib messages
c                         with second order packing. added routine fi636
c                         to process messages with second order packing.
c   93-09-22  cavanaugh   modified to extract sub-center number from
c                         pds byte 26
c   93-10-13  cavanaugh   modified fi634 to correct grid sizes for
c                         grids 204 and 208
c   93-10-14  cavanaugh   increased size of kgds to include entries for
c                         number of points in grid and number of words
c                         in each row
c   93-12-08  cavanaugh   corrected test for edition number instead
c                         of version number
c   93-12-15  cavanaugh   modified second order pointers to first order
c                         values and second order values correctly
c                         in routine fi636
c   94-03-02  cavanaugh   added call to w3fi83 within decoder.  user
c                         no longer needs to make call to this routine
c   94-04-22  cavanaugh   modified fi635, fi636 to process row by row
c                         second order packing, added scaling correction
c                         to fi635, and corrected typographical errors
c                         in comment fields in fi634
c   94-05-17  cavanaugh   corrected error in fi633 to extract resolution
c                         for lambert-conformal grids. added clarifying
c                         information to docblock entries
c   94-05-25  cavanaugh   added code to process column by column as well
c                         as row by row ordering of second order data
c   94-06-27  cavanaugh   added processing for grids 45, 94 and 95.
c                         includes construction of second order bit maps
c                         for thinned grids in fi636.
c   94-07-08  cavanaugh   commented out print outs used for debugging
c   94-09-08  cavanaugh   added grids 220, 221, 223 for fnoc
c   94-11-10  farley      increased mxsize from 72960 to 260000
c                         for .5 degree sst analysis fields
c   94-12-06  r.e.jones   changes in fi632 for pds greater than 28
c   95-02-14  r.e.jones   correct in fi633 for navy wafs grib
c   95-03-20  m.baldwin   fi633 modification to get
c                         data rep types [kgds(1)] 201 and 202 to work.
c   95-04-10  e.rogers    added grids 96 and 97 for eta model in fi634.
c   95-04-26  r.e.jones   fi636 corection for 2nd order complex
c                         unpacking. r
c   95-05-19  r.e.jones   added grid 215, 20 km awips grid
c   95-07-06  r.e.jones   added gaussian t62, t126 grid 98, 126
c   95-10-19  r.e.jones   added grid 216, 45 km eta awips alaska grid 
c   95-10-31  iredell     removed saves and prints
c   96-03-07  r.e.jones   continue unpack with kret error 9 in fi631. 
c   96-08-19  r.e.jones   added mercator grids 8 and 53, and grid 196
c   97-02-12  w bostelman corrects ecmwf us grid 2 processing
c   98-06-17  iredell     removed alternate return in fi637
c   98-08-31  iredell     eliminated need for mxsize
c   98-09-02  gilbert     corrected error in map size for u.s. grid 92
c   98-09-08  baldwin     add data rep type [kgds(1)] 203
c   01-03-08  rogers      changed eta grids 90-97, added eta grids
c                         194, 198. added awips grids 241,242,243,
c                         245, 246, 247, 248, and 250
c   01-03-19  vuong       added awips grids 238,239,240, and 244
c   01-05-03  rogers      added grid 249  (12km for alaska)
c   01-10-10  rogers      redefined grid 218 for 12 km eta
c                         redefined grid 192 for new 32-km eta grid
c   02-03-27  vuong       added rsas grid 88 and awips grids 219, 220,
c                         223, 224, 225, 226, 227, 228, 229, 230, 231,
c                         232, 233, 234, 235, 251, and 252
c   02-08-06  rogers      redefined grids 90-93,97,194,245-250 for the 
c                         8km hi-res-window model and add awips grid 253
c 2003-06-30  gilbert     set new values in array kptr to pass back additional
c                         packing info.
c                         kptr(19) - binary scale factor
c                         kptr(20) - num bits used to pack each datum
c 2003-06-30  gilbert     added grids 145 and 146 for cmaq
c                         and grid 175 for awips over guam.
c 2003-07-08  vuong       added grids 110, 127, 171, 172 and modified grid 170
c 2004-09-02  vuong       added awips grids 147, 148, 173 and 254
c 2005-01-04  cooke       added awips grids 160 and 161
c 2005-03-03  vuong       moved grid 170 to grid 174 and add grid 170
c 2005-03-21  vuong       added awips grids 130
c 2005-10-11  vuong       added awips grids 163
c
c usage:    call w3fi63(msga,kpds,kgds,kbms,data,kptr,kret)
c   input argument list:
c     msga     - grib field - "grib" thru "7777"   char*1
c                   (message can be preceded by junk chars)
c
c   output argument list:
c     data     - array containing data elements
c     kpds     - array containing pds elements.  (edition 1)
c          (1)   - id of center
c          (2)   - generating process id number
c          (3)   - grid definition
c          (4)   - gds/bms flag (right adj copy of octet 8)
c          (5)   - indicator of parameter
c          (6)   - type of level
c          (7)   - height/pressure , etc of level
c          (8)   - year including (century-1)
c          (9)   - month of year
c          (10)  - day of month
c          (11)  - hour of day
c          (12)  - minute of hour
c          (13)  - indicator of forecast time unit
c          (14)  - time range 1
c          (15)  - time range 2
c          (16)  - time range flag
c          (17)  - number included in average
c          (18)  - version nr of grib specification
c          (19)  - version nr of parameter table
c          (20)  - nr missing from average/accumulation
c          (21)  - century of reference time of data
c          (22)  - units decimal scale factor
c          (23)  - subcenter number
c          (24)  - pds byte 29, for nmc ensemble products
c                  128 if forecast field error
c                   64 if bias corrected fcst field
c                   32 if smoothed field
c                  warning: can be combination of more than 1
c          (25)  - pds byte 30, not used
c       (26-35)  - reserved
c       (36-n)   - consecutive bytes extracted from program
c                  definition section (pds) of grib message
c     kgds     - array containing gds elements.
c          (1)   - data representation type
c          (19)  - number of vertical coordinate parameters
c          (20)  - octet number of the list of vertical coordinate
c                  parameters
c                  or
c                  octet number of the list of numbers of points
c                  in each row
c                  or
c                  255 if neither are present
c          (21)  - for grids with pl, number of points in grid
c          (22)  - number of words in each row
c       latitude/longitude grids
c          (2)   - n(i) nr points on latitude circle
c          (3)   - n(j) nr points on longitude meridian
c          (4)   - la(1) latitude of origin
c          (5)   - lo(1) longitude of origin
c          (6)   - resolution flag (right adj copy of octet 17)
c          (7)   - la(2) latitude of extreme point
c          (8)   - lo(2) longitude of extreme point
c          (9)   - di longitudinal direction of increment
c          (10)  - dj latitudinal direction increment
c          (11)  - scanning mode flag (right adj copy of octet 28)
c       gaussian  grids
c          (2)   - n(i) nr points on latitude circle
c          (3)   - n(j) nr points on longitude meridian
c          (4)   - la(1) latitude of origin
c          (5)   - lo(1) longitude of origin
c          (6)   - resolution flag  (right adj copy of octet 17)
c          (7)   - la(2) latitude of extreme point
c          (8)   - lo(2) longitude of extreme point
c          (9)   - di longitudinal direction of increment
c          (10)  - n - nr of circles pole to equator
c          (11)  - scanning mode flag (right adj copy of octet 28)
c          (12)  - nv - nr of vert coord parameters
c          (13)  - pv - octet nr of list of vert coord parameters
c                             or
c                  pl - location of the list of numbers of points in
c                       each row (if no vert coord parameters
c                       are present
c                             or
c                  255 if neither are present
c       polar stereographic grids
c          (2)   - n(i) nr points along lat circle
c          (3)   - n(j) nr points along lon circle
c          (4)   - la(1) latitude of origin
c          (5)   - lo(1) longitude of origin
c          (6)   - resolution flag  (right adj copy of octet 17)
c          (7)   - lov grid orientation
c          (8)   - dx - x direction increment
c          (9)   - dy - y direction increment
c          (10)  - projection center flag
c          (11)  - scanning mode (right adj copy of octet 28)
c       spherical harmonic coefficients
c          (2)   - j pentagonal resolution parameter
c          (3)   - k      "          "         "
c          (4)   - m      "          "         "
c          (5)   - representation type
c          (6)   - coefficient storage mode
c       mercator grids
c          (2)   - n(i) nr points on latitude circle
c          (3)   - n(j) nr points on longitude meridian
c          (4)   - la(1) latitude of origin
c          (5)   - lo(1) longitude of origin
c          (6)   - resolution flag (right adj copy of octet 17)
c          (7)   - la(2) latitude of last grid point
c          (8)   - lo(2) longitude of last grid point
c          (9)   - latit - latitude of projection intersection
c          (10)  - reserved
c          (11)  - scanning mode flag (right adj copy of octet 28)
c          (12)  - longitudinal dir grid length
c          (13)  - latitudinal dir grid length
c       lambert conformal grids
c          (2)   - nx nr points along x-axis
c          (3)   - ny nr points along y-axis
c          (4)   - la1 lat of origin (lower left)
c          (5)   - lo1 lon of origin (lower left)
c          (6)   - resolution (right adj copy of octet 17)
c          (7)   - lov - orientation of grid
c          (8)   - dx - x-dir increment
c          (9)   - dy - y-dir increment
c          (10)  - projection center flag
c          (11)  - scanning mode flag (right adj copy of octet 28)
c          (12)  - latin 1 - first lat from pole of secant cone inter
c          (13)  - latin 2 - second lat from pole of secant cone inter
c       staggered arakawa rotated lat/lon grids (type 203)
c          (2)   - n(i) nr points on latitude circle
c          (3)   - n(j) nr points on longitude meridian
c          (4)   - la(1) latitude of origin
c          (5)   - lo(1) longitude of origin
c          (6)   - resolution flag (right adj copy of octet 17)
c          (7)   - la(2) latitude of center
c          (8)   - lo(2) longitude of center
c          (9)   - di longitudinal direction of increment
c          (10)  - dj latitudinal direction increment
c          (11)  - scanning mode flag (right adj copy of octet 28)
c     kbms       - bitmap describing location of output elements.
c                            (always constructed)
c     kptr       - array containing storage for following parameters
c          (1)   - total length of grib message
c          (2)   - length of indicator (section  0)
c          (3)   - length of pds       (section  1)
c          (4)   - length of gds       (section  2)
c          (5)   - length of bms       (section  3)
c          (6)   - length of bds       (section  4)
c          (7)   - value of current byte
c          (8)   - bit pointer
c          (9)   - grib start bit nr
c         (10)   - grib/grid element count
c         (11)   - nr unused bits at end of section 3
c         (12)   - bit map flag (copy of bms octets 5,6)
c         (13)   - nr unused bits at end of section 2
c         (14)   - bds flags (right adj copy of octet 4)
c         (15)   - nr unused bits at end of section 4
c         (16)   - reserved
c         (17)   - reserved
c         (18)   - reserved
c         (19)   - binary scale factor
c         (20)   - num bits used to pack each datum
c     kret       - flag indicating quality of completion
c
c remarks: when decoding is completed, data at each grid point has been
c          returned in the units specified in the grib manual.
c
c          values for return flag (kret)
c     kret = 0 - normal return, no errors
c          = 1 - 'grib' not found in first 100 chars
c          = 2 - '7777' not in correct location
c          = 3 - unpacked field is larger than 260000
c          = 4 - gds/ grid not one of currently accepted values
c          = 5 - grid not currently avail for center indicated
c          = 8 - temp gds indicated, but gds flag is off
c          = 9 - gds indicates size mismatch with std grid
c          =10 - incorrect center indicator
c          =11 - binary data section (bds) not completely processed.
c                program is not set to process flag combinations
c                shown in octets 4 and 14.
c          =12 - binary data section (bds) not completely processed.
c                program is not set to process flag combinations
c
c   subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 90
c
c$$$
c                                                         4 aug 1988
c                               w3fi63
c
c
c                       grib unpacking routine
c
c
c       this routine will unpack a 'grib' field to the exact grid
c  type specified in the message, return a bit map and make the
c  values of the product definition sec   (pds) and the grid
c  description sec   (gds) available in return arrays.
c  see "grib - the wmo format for the storage of weather product
c  information and the exchange of weather product messages in
c  gridded binary form" dated july 1, 1988 by john d. stackpole
c  doc, noaa, nws, national meteorological center.
c
c       the call to the grib unpacking routine is as follows:
c
c            call w3fi63(msga,kpds,kgds,lbms,data,kptr,kret)
c
c  input:
c
c       msga  = contains the grib message to be unpacked. characters
c               "grib" may begin anywhere within first 100 bytes.
c
c  output:
c
c       kpds(100)      integer
c               array to contain the elements of the product
c               definition sec  .
c         (version 1)
c            kpds(1)  - id of center
c            kpds(2)  - model identification (see "grib" table 1)
c            kpds(3)  - grid identification (see "grib" table 2)
c            kpds(4)  - gds/bms flag
c                           bit       definition
c                            25        0 - gds omitted
c                                      1 - gds included
c                            26        0 - bms omitted
c                                      1 - bms included
c                        note:- leftmost bit = 1,
c                               rightmost bit = 32
c            kpds(5)  - indicator of parameter (see "grib" table 5)
c            kpds(6)  - type of level (see "grib" tables 6 & 7)
c            kpds(7)  - height,pressure,etc  of level
c            kpds(8)  - year including century
c            kpds(9)  - month of year
c            kpds(10) - day of month
c            kpds(11) - hour of day
c            kpds(12) - minute of hour
c            kpds(13) - indicator of forecast time unit (see "grib"
c                       table 8)
c            kpds(14) - time 1               (see "grib" table 8a)
c            kpds(15) - time 2               (see "grib" table 8a)
c            kpds(16) - time range indicator (see "grib" table 8a)
c            kpds(17) - number included in average
c            kpds(18) - edition nr of grib specification
c            kpds(19) - version nr of parameter table
c
c       kgds(13)       integer
c             array containing gds elements.
c
c            kgds(1)  - data representation type
c
c         latitude/longitude grids (see "grib" table 10)
c            kgds(2)  - n(i) number of points on latitude
c                       circle
c            kgds(3)  - n(j) number of points on longitude
c                       circle
c            kgds(4)  - la(1) latitude of origin
c            kgds(5)  - lo(1) longitude of origin
c            kgds(6)  - resolution flag
c                           bit       meaning
c                            25       0 - direction increments not
c                                         given
c                                     1 - direction increments given
c            kgds(7)  - la(2) latitude of extreme point
c            kgds(8)  - lo(2) longitude of extreme point
c            kgds(9)  - di longitudinal direction increment
c            kgds(10) - regular lat/lon grid
c                           dj - latitudinal direction
c                                increment
c                       gaussian grid
c                           n  - number of latitude circles
c                                between a pole and the equator
c            kgds(11) - scanning mode flag
c                           bit       meaning
c                            25       0 - points along a latitude
c                                         scan from west to east
c                                     1 - points along a latitude
c                                         scan from east to west
c                            26       0 - points along a meridian
c                                         scan from north to south
c                                     1 - points along a meridian
c                                         scan from south to north
c                            27       0 - points scan first along
c                                         circles of latitude, then
c                                         along meridians
c                                         (fortran: (i,j))
c                                     1 - points scan first along
c                                         meridians then along
c                                         circles of latitude
c                                         (fortran: (j,i))
c
c         polar stereographic grids  (see grib table 12)
c            kgds(2)  - n(i) nr points along lat circle
c            kgds(3)  - n(j) nr points along lon circle
c            kgds(4)  - la(1) latitude of origin
c            kgds(5)  - lo(1) longitude of origin
c            kgds(6)  - reserved
c            kgds(7)  - lov grid orientation
c            kgds(8)  - dx - x direction increment
c            kgds(9)  - dy - y direction increment
c            kgds(10) - projection center flag
c            kgds(11) - scanning mode
c
c         spherical harmonic coefficients (see "grib" table 14)
c            kgds(2)  - j pentagonal resolution parameter
c            kgds(3)  - k pentagonal resolution parameter
c            kgds(4)  - m pentagonal resolution parameter
c            kgds(5)  - representation type
c            kgds(6)  - coefficient storage mode
c
c       mercator grids
c            kgds(2)   - n(i) nr points on latitude circle
c            kgds(3)   - n(j) nr points on longitude meridian
c            kgds(4)   - la(1) latitude of origin
c            kgds(5)   - lo(1) longitude of origin
c            kgds(6)   - resolution flag
c            kgds(7)   - la(2) latitude of last grid point
c            kgds(8)   - lo(2) longitude of last grid point
c            kgds(9)   - latin - latitude of projection intersection
c            kgds(10)  - reserved
c            kgds(11)  - scanning mode flag
c            kgds(12)  - longitudinal dir grid length
c            kgds(13)  - latitudinal dir grid length
c       lambert conformal grids
c            kgds(2)   - nx nr points along x-axis
c            kgds(3)   - ny nr points along y-axis
c            kgds(4)   - la1 lat of origin (lower left)
c            kgds(5)   - lo1 lon of origin (lower left)
c            kgds(6)   - resolution (right adj copy of octet 17)
c            kgds(7)   - lov - orientation of grid
c            kgds(8)   - dx - x-dir increment
c            kgds(9)   - dy - y-dir increment
c            kgds(10)  - projection center flag
c            kgds(11)  - scanning mode flag
c            kgds(12)  - latin 1 - first lat from pole of
c                        secant cone intersection
c            kgds(13)  - latin 2 - second lat from pole of
c                        secant cone intersection
c
c       lbms(*)    logical
c               array to contain the bit map describing the
c               placement of data in the output array.  if a
c               bit map is not included in the source message,
c               one will be generated automatically by the
c               unpacking routine.
c
c
c       data(*)    real
c               this array will contain the unpacked data points.
c
c                      note:- 65160 is maximun field size allowable
c
c       kptr(10)       integer
c               array containing storage for the following
c               parameters.
c
c                 (1)  -    unused
c                 (2)  -    unused
c                 (3)  -    length of pds (in bytes)
c                 (4)  -    length of gds (in bytes)
c                 (5)  -    length of bms (in bytes)
c                 (6)  -    length of bds (in bytes)
c                 (7)  -    used by unpacking routine
c                 (8)  -    number of data points for grid
c                 (9)  -    "grib" characters start in byte number
c                 (10) -    used by unpacking routine
c
c
c       kret      integer
c                 this variable will contain the return indicator.
c
c                 0    -    no errors detected.
c
c                 1    -    'grib' not found in first 100
c                           characters.
c
c                 2    -    '7777' not found, either missing or
c                           total of sec   counts of individual
c                           sections  is incorrect.
c
c                 3    -    unpacked field is larger than 65160.
c
c                 4    -    in gds, data representation type
c                           not one of the currently acceptable
c                           values. see "grib" table 9. value
c                           of incorrect type returned in kgds(1).
c
c                 5    -    grid indicated in kpds(3) is not
c                           available for the center indicated in
c                           kpds(1) and no gds sent.
c
c                 7    -    edition indicated in kpds(18) has not
c                           yet been included in the decoder.
c
c                 8    -    grid identification = 255 (not standard
c                           grid) but flag indicating presence of
c                           gds is turned off. no method of
c                           generating proper grid.
c
c                 9    -    product of kgds(2) and kgds(3) does not
c                           match standard number of points for this
c                           grid (for other than spectrals). this
c                           will occur only if the grid.
c                           identification, kpds(3), and a
c                           transmitted gds are inconsistent.
c
c                10    -    center indicator was not one indicated
c                           in "grib" table 1.  please contact ad
c                           production management branch (w/nmc42)
c                                     if this error is encountered.
c
c                11    -    binary data section (bds) not completely
c                           processed.  program is not set to process
c                           flag combinations as shown in
c                           octets 4 and 14.
c
c
c  list of text messages from code
c
c
c  w3fi63/fi632
c
c            'have encountered a new grid for nmc, please notify
c            automation division, production management branch
c            (w/nmc42)'
c
c            'have encountered a new grid for ecmwf, please notify
c            automation division, production management branch
c            (w/nmc42)'
c
c            'have encountered a new grid for u.k. meteorological
c            office, bracknell.  please notify automation division,
c            production management branch (w/nmc42)'
c
c            'have encountered a new grid for fnoc, please notify
c            automation division, production management branch
c            (w/nmc42)'
c
c
c  w3fi63/fi633
c
c            'polar stereo processing not available'  *
c
c  w3fi63/fi634
c
c            'warning - bit map may not be associated with spherical
c            coefficients'
c
c
c  w3fi63/fi637
c
c            'no current listing of fnoc grids'      *
c
c
c  * will be available in next update
c  ***************************************************************
c
c                       incoming message holder
      character*1   msga(*)
c                       bit map
      logical*1     kbms(*)
c
c                       elements of product description sec   (pds)
      integer       kpds(*)
c                       elements of grid description sec   (pds)
      integer       kgds(*)
c
c                       container for grib grid
      real          data(*)
c
c                       array of pointers and counters
      integer       kptr(*)
c
c  *****************************************************************
      integer       kkk,jsgn,jexp,ifr,npts
      character     kk(8)
      real          realkk,fval1,fdiff1
      equivalence   (kk(1),kkk)
c  *****************************************************************
c        1.0 locate beginning of 'grib' message
c             find 'grib' characters
c        2.0  use counts in each description sec   to determine
c             if '7777' is in proper place.
c        3.0  parse product definition section.
c        4.0  parse grid description sec   (if included)
c        5.0  parse bit map sec   (if included)
c        6.0  using information from product definition, grid
c                  description, and bit map sections.. extract
c                  data and place into proper array.
c  *******************************************************************
c
c                      main driver
c
c  *******************************************************************
      kptr(10) = 0
c                  see if proper 'grib' key exists, then
c                  using sec   counts, determine if '7777'
c                  is in the proper location
c
      call fi631(msga,kptr,kpds,kret)
      if(kret.ne.0) then
          go to 900
      end if
c     print *,'fi631 kptr',(kptr(i),i=1,16)
c
c                  parse parameters from product description section
c
      call fi632(msga,kptr,kpds,kret)
      if(kret.ne.0) then
          go to 900
      end if
c     print *,'fi632 kptr',(kptr(i),i=1,16)
c
c                  if available, extract new grid description
c
      if (iand(kpds(4),128).ne.0) then
          call fi633(msga,kptr,kgds,kret)
          if(kret.ne.0) then
              go to 900
          end if
c         print *,'fi633 kptr',(kptr(i),i=1,16)
      end if
c
c                  extract or generate bit map
c
      call fi634(msga,kptr,kpds,kgds,kbms,kret)
      if (kret.ne.0) then
        if (kret.ne.9) then
          go to 900
        end if
      end if
c     print *,'fi634 kptr',(kptr(i),i=1,16)
c
c                  using information from pds, bms and bit data sec  ,
c                  extract and save in grib grid, all data entries.
c
      if (kpds(18).eq.1) then
          call fi635(msga,kptr,kpds,kgds,kbms,data,kret)
          if (kptr(3).eq.50) then
c
c                     pds equal 50 bytes
c                        therefore something special is going on
c
c                        in this case 2nd difference packing
c                                needs to be undone.
c
c                   extract first value from byte 41-44 pds
c                              kptr(9) contains offset to start of
c                              grib message.
c                   extract first first-difference from bytes 45-48 pds
c
c                  and extract scale factor (e) to undo 2**e
c                  that was applied prior to 2nd order packing
c                  and placed in pds bytes 49-51
c                  factor is a signed two byte integer
c
c                  also need the decimal scaling from pds(27-28)
c                  (available in kpds(22) from unpacker)
c                  to undo the decimal scaling applied to the
c                  second differences during unpacking.
c                  second diffs always packed with 0 decimal scale
c                  but unpacker doesnt know that.
c
c             call gbyte  (msga,fval1,kptr(9)+384,32)
c
c         note integers, characters and equivalences
c         defined above to make this kkk extraction
c         work and line up on word boundaries
c
          call gbyte (msga,kkk,kptr(9)+384,32)
c
c       the next code will convert the ibm370 foating point
c       to the floating point used on your machine.
c
c       1st test to see in on 32 or 64 bit word machine
c       lw = 4 or 8; if 8 may be a cray
c
              call w3fi01(lw)
              if (lw.eq.4) then
                  call gbyte (kk,jsgn,0,1)
                  call gbyte (kk,jexp,1,7)
                  call gbyte (kk,ifr,8,24)
              else
                  call gbyte (kk,jsgn,32,1)
                  call gbyte (kk,jexp,33,7)
                  call gbyte (kk,ifr,40,24)
              endif
c
              if (ifr.eq.0) then
                  realkk = 0.0
              else if (jexp.eq.0.and.ifr.eq.0) then
                  realkk = 0.0
              else
                  realkk = float(ifr) * 16.0 ** (jexp - 64 - 6)
                  if (jsgn.ne.0) realkk = -realkk
              end if
              fval1 = realkk
c
c             call gbyte  (msga,fdiff1,kptr(9)+416,32)
c          (replaced by following extraction)
c
              call gbyte (msga,kkk,kptr(9)+416,32)
c
c       the next code will convert the ibm370 foating point
c       to the floating point used on your machine.
c
c       1st test to see in on 32 or 64 bit word machine
c       lw = 4 or 8; if 8 may be a cray
c
              call w3fi01(lw)
              if (lw.eq.4) then
                  call gbyte (kk,jsgn,0,1)
                  call gbyte (kk,jexp,1,7)
                  call gbyte (kk,ifr,8,24)
              else
                  call gbyte (kk,jsgn,32,1)
                  call gbyte (kk,jexp,33,7)
                  call gbyte (kk,ifr,40,24)
              endif
c
              if (ifr.eq.0) then
                  realkk = 0.0
              else if (jexp.eq.0.and.ifr.eq.0) then
                  realkk = 0.0
              else
                  realkk = float(ifr) * 16.0 ** (jexp - 64 - 6)
                  if (jsgn.ne.0) realkk = -realkk
              end if
              fdiff1 = realkk
c
              call gbyte  (msga,isign,kptr(9)+448,1)
              call gbyte  (msga,iscal2,kptr(9)+449,15)
              if(isign.gt.0) then
                  iscal2 = - iscal2
              endif
c             print *,'delta point 1-',fval1
c             print *,'delta point 2-',fdiff1
c             print *,'delta point 3-',iscal2
              npts  = kptr(10)
c             write (6,fmt='(''  2nd diff points in field = '',/,
c    &         10(3x,10f12.2,/))') (data(i),i=1,npts)
c             print *,'delta point 4-',kpds(22)
              call w3fi83 (data,npts,fval1,fdiff1,
     &                            iscal2,kpds(22),kpds,kgds)
c             write (6,fmt='(''  2nd diff expanded points in field = '',
c    &            /,10(3x,10f12.2,/))') (data(i),i=1,npts)
c             write (6,fmt='(''  end of array in field = '',/,
c    &         10(3x,10f12.2,/))') (data(i),i=npts-5,npts)
          end if
      else
c         print *,'fi635 not programmed for edition nr',kpds(18)
          kret   = 7
      end if
c
  900 return
      end
      subroutine fi631(msga,kptr,kpds,kret)
c$$$  subprogram documentation  block
c                .      .    .                                       .
c subprogram:    fi631       find 'grib' chars & reset pointers
c   prgmmr: bill cavanaugh   org: w/nmc42    date: 91-09-13
c
c abstract: find 'grib; characters and set pointers to the next
c   byte following 'grib'. if they exist extract counts from gds and
c   bms. extract count from bds. determine if sum of counts actually
c   places terminator '7777' at the correct location.
c
c program history log:
c   91-09-13  cavanaugh
c   95-10-31  iredell     removed saves and prints
c
c usage:    call fi631(msga,kptr,kpds,kret)
c   input argument list:
c     msga       - grib field - "grib" thru "7777"
c     kptr       - array containing storage for following parameters
c          (1)   - total length of grib message
c          (2)   - length of indicator (section  0)
c          (3)   - length of pds       (section  1)
c          (4)   - length of gds       (section  2)
c          (5)   - length of bms       (section  3)
c          (6)   - length of bds       (section  4)
c          (7)   - value of current byte
c          (8)   - bit pointer
c          (9)   - grib start bit nr
c         (10)   - grib/grid element count
c         (11)   - nr unused bits at end of section 3
c         (12)   - bit map flag
c         (13)   - nr unused bits at end of section 2
c         (14)   - bds flags
c         (15)   - nr unused bits at end of section 4
c
c   output argument list:      (including work arrays)
c     kpds     - array containing pds elements.
c          (1)   - id of center
c          (2)   - model identification
c          (3)   - grid identification
c          (4)   - gds/bms flag
c          (5)   - indicator of parameter
c          (6)   - type of level
c          (7)   - height/pressure , etc of level
c          (8)   - year of century
c          (9)   - month of year
c          (10)  - day of month
c          (11)  - hour of day
c          (12)  - minute of hour
c          (13)  - indicator of forecast time unit
c          (14)  - time range 1
c          (15)  - time range 2
c          (16)  - time range flag
c          (17)  - number included in average
c     kptr       - see input list
c     kret       - error return
c
c remarks:
c     error returns
c     kret  = 1  -  no 'grib'
c             2  -  no '7777' or mislocated (by counts)
c
c   subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 77
c   machine:  hds9000
c
c$$$
c
c                       incoming message holder
      character*1   msga(*)
c                       array of pointers and counters
      integer       kptr(*)
c                       product description section data.
      integer       kpds(*)
c
      integer       kret
c
c  ******************************************************************
      kret = 0
c  -------------------  find 'grib' key
      do 50 i = 0, 839, 8
          call gbyte (msga,mgrib,i,32)
          if (mgrib.eq.1196575042) then
              kptr(9)   = i
              go to 60
          end if
   50 continue
      kret  = 1
      return
   60 continue
c  -------------found 'grib'
c                        skip grib characters
c     print *,'fi631 grib at',i
      kptr(8)   = kptr(9) + 32
      call gbyte (msga,itotal,kptr(8),24)
c                    have lifted what may be a msg total byte count
      ipoint    = kptr(9) + itotal * 8 - 32
      call gbyte (msga,i7777,ipoint,32)
      if (i7777.eq.926365495) then
c                 have found end of message '7777' in proper location
c                 mark and process as grib version 1 or higher
c         print *,'fi631 7777 at',ipoint
          kptr(8)   = kptr(8) + 24
          kptr(1)   = itotal
          kptr(2)   = 8
          call gbyte (msga,kpds(18),kptr(8),8)
          kptr(8)   = kptr(8) + 8
      else
c                 cannot find end of grib edition 1 message
          kret      = 2
          return
      end if
c  -------------------  process section 1
c                   extract count from pds
c     print *,'start of pds',kptr(8)
      call gbyte (msga,kptr(3),kptr(8),24)
      look      = kptr(8) + 56
c                   extract gds/bms flag
      call gbyte (msga,kpds(4),look,8)
      kptr(8)   = kptr(8) + kptr(3) * 8
c     print *,'start of gds',kptr(8)
      if (iand(kpds(4),128).ne.0) then
c                   extract count from gds
          call gbyte (msga,kptr(4),kptr(8),24)
          kptr(8)   = kptr(8) + kptr(4) * 8
      else
          kptr(4)   = 0
      end if
c     print *,'start of bms',kptr(8)
      if (iand(kpds(4),64).ne.0) then
c                   extract count from bms
          call gbyte (msga,kptr(5),kptr(8),24)
      else
          kptr(5)   = 0
      end if
      kptr(8)   = kptr(8) + kptr(5) * 8
c     print *,'start of bds',kptr(8)
c                   extract count from bds
      call gbyte (msga,kptr(6),kptr(8),24)
c  ---------------  test for '7777'
c     print *,(kptr(kj),kj=1,10)
      kptr(8)   = kptr(8) + kptr(6) * 8
c                   extract four bytes from this location
c     print *,'fi631 looking for 7777 at',kptr(8)
      call gbyte (msga,k7777,kptr(8),32)
      match  = kptr(2) + kptr(3) + kptr(4) + kptr(5) + kptr(6) + 4
      if (k7777.ne.926365495.or.match.ne.kptr(1)) then
          kret  = 2
      else
c         print *,'fi631 7777 at',kptr(8)
          if (kpds(18).eq.0) then
              kptr(1)  = kptr(2) + kptr(3) + kptr(4) + kptr(5) +
     *                kptr(6) + 4
          end if
      end if
c     print *,'kptr',(kptr(i),i=1,16)
      return
      end
      subroutine fi632(msga,kptr,kpds,kret)
c$$$  subprogram documentation  block
c                .      .    .                                       .
c subprogram:    fi632       gather info from product definition sec
c   prgmmr: bill cavanaugh   org: w/nmc42    date: 91-09-13
c
c abstract: extract information from the product description
c   sec  , and generate label information to permit storage
c   in office note 84 format.
c
c program history log:
c   91-09-13  cavanaugh
c   93-12-08  cavanaugh   corrected test for edition number instead
c                         of version number
c   95-10-31  iredell     removed saves and prints
c   99-01-20  baldwin     modified to handle grid 237
c
c usage:    call fi632(msga,kptr,kpds,kret)
c   input argument list:
c     msga      - array containing grib message
c     kptr       - array containing storage for following parameters
c          (1)   - total length of grib message
c          (2)   - length of indicator (section  0)
c          (3)   - length of pds       (section  1)
c          (4)   - length of gds       (section  2)
c          (5)   - length of bms       (section  3)
c          (6)   - length of bds       (section  4)
c          (7)   - value of current byte
c          (8)   - bit pointer
c          (9)   - grib start bit nr
c         (10)   - grib/grid element count
c         (11)   - nr unused bits at end of section 3
c         (12)   - bit map flag
c         (13)   - nr unused bits at end of section 2
c         (14)   - bds flags
c         (15)   - nr unused bits at end of section 4
c
c   output argument list:      (including work arrays)
c     kpds     - array containing pds elements.
c          (1)   - id of center
c          (2)   - model identification
c          (3)   - grid identification
c          (4)   - gds/bms flag
c          (5)   - indicator of parameter
c          (6)   - type of level
c          (7)   - height/pressure , etc of level
c          (8)   - year of century
c          (9)   - month of year
c          (10)  - day of month
c          (11)  - hour of day
c          (12)  - minute of hour
c          (13)  - indicator of forecast time unit
c          (14)  - time range 1
c          (15)  - time range 2
c          (16)  - time range flag
c          (17)  - number included in average
c          (18)  -
c          (19)  -
c          (20)  - number missing from avgs/accumulations
c          (21)  - century
c          (22)  - units decimal scale factor
c          (23)  - subcenter
c     kptr       - array containing storage for following parameters
c                  see input list
c     kret   - error return
c
c remarks:
c        error return = 0 - no errors
c                     = 8 - temp gds indicated, but no gds
c
c   subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 77
c   machine:  hds9000
c
c$$$
c
c                       incoming message holder
      character*1   msga(*)
c
c                       array of pointers and counters
      integer       kptr(*)
c                       product description section entries
      integer       kpds(*)
c
      integer       kret
c  -------------------  process section 1
      kptr(8)  = kptr(9) + kptr(2) * 8 + 24
c  byte 4
c                   parameter table version nr
          call gbyte (msga,kpds(19),kptr(8),8)
          kptr(8)   = kptr(8) + 8
c  byte 5           identification of center
      call gbyte (msga,kpds(1),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 6
c                       get generating process id nr
      call gbyte (msga,kpds(2),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 7
c                      grid definition
      call gbyte (msga,kpds(3),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 8
c                      gds/bms flags
c     call gbyte (msga,kpds(4),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 9
c                      indicator of parameter
      call gbyte (msga,kpds(5),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 10
c                      type of level
      call gbyte (msga,kpds(6),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 11,12
c                      height/pressure
      call gbyte (msga,kpds(7),kptr(8),16)
      kptr(8)   = kptr(8) + 16
c  byte 13
c                      year of century
      call gbyte (msga,kpds(8),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 14
c                      month of year
      call gbyte (msga,kpds(9),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 15
c                      day of month
      call gbyte (msga,kpds(10),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 16
c                      hour of day
      call gbyte (msga,kpds(11),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 17
c                      minute
      call gbyte (msga,kpds(12),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 18
c                      indicator time unit range
      call gbyte (msga,kpds(13),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 19
c                      p1 - period of time
      call gbyte (msga,kpds(14),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 20
c                      p2 - period of time
      call gbyte (msga,kpds(15),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 21
c                      time range indicator
      call gbyte (msga,kpds(16),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c
c     if time range indicator is 10, p1 is packed in
c     pds bytes 19-20
c
      if (kpds(16).eq.10) then
          kpds(14)  = kpds(14) * 256 + kpds(15)
          kpds(15)  = 0
      end if
c  byte 22,23
c                      number included in average
      call gbyte (msga,kpds(17),kptr(8),16)
      kptr(8)   = kptr(8) + 16
c  byte 24
c                      number missing from averages/accumulations
      call gbyte (msga,kpds(20),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c  byte 25
c                      identification of century
      call gbyte (msga,kpds(21),kptr(8),8)
      kptr(8)   = kptr(8) + 8
      if (kptr(3).gt.25) then
c  byte 26              sub center number
          call gbyte (msga,kpds(23),kptr(8),8)
          kptr(8)   = kptr(8) + 8
          if (kptr(3).ge.28) then
c  byte 27-28
c                          units decimal scale factor
              call gbyte (msga,isign,kptr(8),1)
              kptr(8)  = kptr(8) + 1
              call gbyte (msga,idec,kptr(8),15)
              kptr(8)  = kptr(8) + 15
              if (isign.gt.0) then
                  kpds(22)  = - idec
              else
                  kpds(22)  = idec
              end if
              isiz  = kptr(3) - 28
              if (isiz.le.12) then
c  byte  29
                  call gbyte (msga,kpds(24),kptr(8)+8,8)
c  byte  30
                  call gbyte (msga,kpds(25),kptr(8)+16,8)
c  bytes 31-40                  currently reserved for future use
                  kptr(8)  = kptr(8) + isiz * 8
              else
c  byte  29
                  call gbyte (msga,kpds(24),kptr(8)+8,8)
c  byte  30
                  call gbyte (msga,kpds(25),kptr(8)+16,8)
c  bytes 31-40                  currently reserved for future use
                  kptr(8)  = kptr(8) + 12 * 8
c  bytes 41 - n                 local use data
                  call w3fi01(lw)
                  mwdbit  = lw * 8
                  isiz    = kptr(3) - 40
                  iter    = isiz / lw
                  if (mod(isiz,lw).ne.0) iter = iter + 1
                  call gbytes (msga,kpds(36),kptr(8),mwdbit,0,iter)
                  kptr(8)  = kptr(8) + isiz * 8
              end if
          end if
      end if
c  ----------- test for new grid
      if (iand(kpds(4),128).ne.0) then
          if (iand(kpds(4),64).ne.0) then
              if (kpds(3).ne.255) then
                  if (kpds(3).ge.21.and.kpds(3).le.26)then
                      return
                  else if (kpds(3).ge.37.and.kpds(3).le.44)then
                      return
                  else if (kpds(3).ge.61.and.kpds(3).le.64) then
                      return
                  end if
                  if (kpds(1).eq.7) then
                      if (kpds(3).ge.2.and.kpds(3).le.3) then
                      else if (kpds(3).ge.5.and.kpds(3).le.6) then
                      else if (kpds(3).eq.8) then
                      else if (kpds(3).ge.27.and.kpds(3).le.34) then
                      else if (kpds(3).eq.50) then
                      else if (kpds(3).eq.53) then
                      else if (kpds(3).ge.70.and.kpds(3).le.77) then
                      else if (kpds(3).eq.98) then
                      else if (kpds(3).ge.100.and.kpds(3).le.105) then
                      else if (kpds(3).eq.126) then
                      else if (kpds(3).eq.196) then
                      else if (kpds(3).ge.201.and.kpds(3).le.237) then
                      else
c                         print *,' have encountered a new grid for',
c    *                    ' nmc without a grid description section'
c                         print *,' please notify automation division'
c                         print *,' production management branch'
c                         print *,' w/nmc42)'
                      end if
                  else if (kpds(1).eq.98) then
                      if (kpds(3).ge.1.and.kpds(3).le.16) then
                      else
c                         print *,' have encountered a new grid for',
c    *                    ' ecmwf without a grid description section'
c                         print *,' please notify automation division'
c                         print *,' production management branch'
c                         print *,' w/nmc42)'
                      end if
                  else if (kpds(1).eq.74) then
                      if (kpds(3).ge.1.and.kpds(3).le.12) then
                      else if (kpds(3).ge.21.and.kpds(3).le.26)then
                      else if (kpds(3).ge.61.and.kpds(3).le.64) then
                      else if (kpds(3).ge.70.and.kpds(3).le.77) then
                      else
c                         print *,' have encountered a new grid for',
c    *                            ' u.k. met office, bracknell',
c    *                            ' without a grid description section'
c                         print *,' please notify automation division'
c                         print *,' production management branch'
c                         print *,' w/nmc42)'
                      end if
                  else if (kpds(1).eq.58) then
                      if (kpds(3).ge.1.and.kpds(3).le.12) then
                      else
c                         print *,' have encountered a new grid for',
c    *                      ' fnoc without a grid description section'
c                         print *,' please notify automation division'
c                         print *,' production management branch'
c                         print *,' w/nmc42)'
                      end if
                  end if
              end if
          end if
      end if
      return
      end
      subroutine fi633(msga,kptr,kgds,kret)
c$$$  subprogram documentation  block
c                .      .    .                                       .
c subprogram:    fi633       extract info from grib-gds
c   prgmmr: bill cavanaugh   org: w/nmc42    date: 91-09-13
c
c abstract: extract information on unlisted grid to allow
c   conversion to office note 84 format.
c
c program history log:
c   91-09-13  cavanaugh
c   95-03-20  m.baldwin   fi633 modification to get
c                         data rep types [kgds(1)] 201 and 202 to work.
c   95-10-31  iredell     removed saves and prints
c   98-09-08  baldwin     add data rep type [kgds(1)] 203
c                        
c
c usage:    call fi633(msga,kptr,kgds,kret)
c   input argument list:
c     msga      - array containing grib message
c     kptr       - array containing storage for following parameters
c          (1)   - total length of grib message
c          (2)   - length of indicator (section  0)
c          (3)   - length of pds       (section  1)
c          (4)   - length of gds       (section  2)
c          (5)   - length of bms       (section  3)
c          (6)   - length of bds       (section  4)
c          (7)   - value of current byte
c          (8)   - bit pointer
c          (9)   - grib start bit nr
c         (10)   - grib/grid element count
c         (11)   - nr unused bits at end of section 3
c         (12)   - bit map flag
c         (13)   - nr unused bits at end of section 2
c         (14)   - bds flags
c         (15)   - nr unused bits at end of section 4
c
c   output argument list:      (including work arrays)
c     kgds     - array containing gds elements.
c          (1)   - data representation type
c          (19)  - number of vertical coordinate parameters
c          (20)  - octet number of the list of vertical coordinate
c                  parameters
c                  or
c                  octet number of the list of numbers of points
c                  in each row
c                  or
c                  255 if neither are present
c          (21)  - for grids with pl, number of points in grid
c          (22)  - number of words in each row
c       latitude/longitude grids
c          (2)   - n(i) nr points on latitude circle
c          (3)   - n(j) nr points on longitude meridian
c          (4)   - la(1) latitude of origin
c          (5)   - lo(1) longitude of origin
c          (6)   - resolution flag
c          (7)   - la(2) latitude of extreme point
c          (8)   - lo(2) longitude of extreme point
c          (9)   - di longitudinal direction of increment
c          (10)  - dj latitudinal direction increment
c          (11)  - scanning mode flag
c       polar stereographic grids
c          (2)   - n(i) nr points along lat circle
c          (3)   - n(j) nr points along lon circle
c          (4)   - la(1) latitude of origin
c          (5)   - lo(1) longitude of origin
c          (6)   - reserved
c          (7)   - lov grid orientation
c          (8)   - dx - x direction increment
c          (9)   - dy - y direction increment
c          (10)  - projection center flag
c          (11)  - scanning mode
c       spherical harmonic coefficients
c          (2)   - j pentagonal resolution parameter
c          (3)   - k      "          "         "
c          (4)   - m      "          "         "
c          (5)   - representation type
c          (6)   - coefficient storage mode
c       mercator grids
c          (2)   - n(i) nr points on latitude circle
c          (3)   - n(j) nr points on longitude meridian
c          (4)   - la(1) latitude of origin
c          (5)   - lo(1) longitude of origin
c          (6)   - resolution flag
c          (7)   - la(2) latitude of last grid point
c          (8)   - lo(2) longitude of last grid point
c          (9)   - latin - latitude of projection intersection
c          (10)  - reserved
c          (11)  - scanning mode flag
c          (12)  - longitudinal dir grid length
c          (13)  - latitudinal dir grid length
c       lambert conformal grids
c          (2)   - nx nr points along x-axis
c          (3)   - ny nr points along y-axis
c          (4)   - la1 lat of origin (lower left)
c          (5)   - lo1 lon of origin (lower left)
c          (6)   - resolution (right adj copy of octet 17)
c          (7)   - lov - orientation of grid
c          (8)   - dx - x-dir increment
c          (9)   - dy - y-dir increment
c          (10)  - projection center flag
c          (11)  - scanning mode flag
c          (12)  - latin 1 - first lat from pole of secant cone inter
c          (13)  - latin 2 - second lat from pole of secant cone inter
c       staggered arakawa rotated lat/lon grids (203)
c          (2)   - n(i) nr points on rotated latitude circle
c          (3)   - n(j) nr points on rotated longitude meridian
c          (4)   - la(1) latitude of origin
c          (5)   - lo(1) longitude of origin
c          (6)   - resolution flag
c          (7)   - la(2) latitude of center
c          (8)   - lo(2) longitude of center
c          (9)   - di longitudinal direction of increment
c          (10)  - dj latitudinal direction increment
c          (11)  - scanning mode flag
c     kptr       - array containing storage for following parameters
c                  see input list
c     kret       - error return
c
c remarks:
c     kret = 0
c          = 4   - data representation type not currently acceptable
c
c   subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 77
c   machine:  hds9000
c
c$$$
c  ************************************************************
c                       incoming message holder
      character*1   msga(*)
c
c                       array gds elements
      integer       kgds(*)
c                       array of pointers and counters
      integer       kptr(*)
c
      integer       kret
c  ---------------------------------------------------------------
      kret    = 0
c                process grid definition section (if present)
c             make sure bit pointer is properly set
      kptr(8)  = kptr(9) + (kptr(2)*8) + (kptr(3)*8) + 24
      nsave    = kptr(8) - 24
c  byte 4
c                   nv - nr of vert coord parameters
      call gbyte (msga,kgds(19),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  byte 5
c                   pv - location - see fm92 manual
      call gbyte (msga,kgds(20),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  byte 6
c                      data representation type
      call gbyte (msga,kgds(1),kptr(8),8)
      kptr(8)   = kptr(8) + 8
c           bytes 7-32 are grid definition depending on
c           data representation type
      if (kgds(1).eq.0) then
          go to 1000
      else if (kgds(1).eq.1) then
          go to 4000
      else if (kgds(1).eq.2.or.kgds(1).eq.5) then
          go to 2000
      else if (kgds(1).eq.3) then
          go to 5000
      else if (kgds(1).eq.4) then
          go to 1000
c     else if (kgds(1).eq.10) then
c     else if (kgds(1).eq.14) then
c     else if (kgds(1).eq.20) then
c     else if (kgds(1).eq.24) then
c     else if (kgds(1).eq.30) then
c     else if (kgds(1).eq.34) then
      else if (kgds(1).eq.50) then
          go to 3000
c     else if (kgds(1).eq.60) then
c     else if (kgds(1).eq.70) then
c     else if (kgds(1).eq.80) then
      else if (kgds(1).eq.201.or.kgds(1).eq.202.or.kgds(1).eq.203) then
          go to 1000
      else
c                      mark as gds/ unknown data representation type
          kret     = 4
          return
      end if
c     byte 33-n   vertical coordinate parameters
c  -----------
c     bytes 33-42 extensions of grid definition for rotation
c                 or stretching of the coordinate system or
c                 lambert conformal projection.
c     byte 43-n   vertical coordinate parameters
c  -----------
c     bytes 33-52 extensions of grid definition for stretched
c                 and rotated coordinate system
c     byte 53-n   vertical coordinate parameters
c  -----------
c ************************************************************
c  ------------------- latitude/longitude grids
c  ------------------- arakawa staggered, semi-staggered, or filled
c                          rotated lat/lon grids
c
c  ------------------- byte 7-8     nr of points along latitude circle
 1000 continue
      call gbyte (msga,kgds(2),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 9-10    nr of points along long meridian
      call gbyte (msga,kgds(3),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 11-13   latitude of origin
      call gbyte (msga,kgds(4),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(4),8388608).ne.0) then
          kgds(4)  =  iand(kgds(4),8388607) * (-1)
      end if
c  ------------------- byte 14-16   longitude of origin
      call gbyte (msga,kgds(5),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(5),8388608).ne.0) then
          kgds(5)  =  - iand(kgds(5),8388607)
      end if
c  ------------------- byte 17      resolution flag
      call gbyte (msga,kgds(6),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 18-20   latitude of last grid point
      call gbyte (msga,kgds(7),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(7),8388608).ne.0) then
          kgds(7)  =  - iand(kgds(7),8388607)
      end if
c  ------------------- byte 21-23   longitude of last grid point
      call gbyte (msga,kgds(8),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(8),8388608).ne.0) then
          kgds(8)  =  - iand(kgds(8),8388607)
      end if
c  ------------------- byte 24-25   latitudinal dir increment
      call gbyte (msga,kgds(9),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 26-27   if regular lat/lon grid
c                                       have longit dir increment
c                                   else if gaussian grid
c                                       have nr of lat circles
c                                       between pole and equator
      call gbyte (msga,kgds(10),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 28      scanning mode flags
      call gbyte (msga,kgds(11),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 29-32   reserved
c                             skip to start of byte 33
      call gbyte (msga,kgds(12),kptr(8),32)
      kptr(8)  = kptr(8) + 32
c  -------------------
      go to 900
c  ******************************************************************
c            ' polar stereo processing '
c
c  ------------------- byte 7-8     nr of points along x=axis
 2000 continue
      call gbyte (msga,kgds(2),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 9-10    nr of points along y-axis
      call gbyte (msga,kgds(3),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 11-13   latitude of origin
      call gbyte (msga,kgds(4),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(4),8388608).ne.0) then
          kgds(4)  =  - iand(kgds(4),8388607)
      end if
c  ------------------- byte 14-16   longitude of origin
      call gbyte (msga,kgds(5),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(5),8388608).ne.0) then
          kgds(5)  =   - iand(kgds(5),8388607)
      end if
c  ------------------- byte 17      reserved
      call gbyte (msga,kgds(6),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 18-20   lov orientation of the grid
      call gbyte (msga,kgds(7),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(7),8388608).ne.0) then
          kgds(7)  =  - iand(kgds(7),8388607)
      end if
c  ------------------- byte 21-23   dx - the x direction increment
      call gbyte (msga,kgds(8),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(8),8388608).ne.0) then
          kgds(8)  =  - iand(kgds(8),8388607)
      end if
c  ------------------- byte 24-26   dy - the y direction increment
      call gbyte (msga,kgds(9),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(9),8388608).ne.0) then
          kgds(9)  =  - iand(kgds(9),8388607)
      end if
c  ------------------- byte 27      projection center flag
      call gbyte (msga,kgds(10),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 28      scanning mode
      call gbyte (msga,kgds(11),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 29-32   reserved
c                             skip to start of byte 33
      call gbyte (msga,kgds(12),kptr(8),32)
      kptr(8)  = kptr(8) + 32
c
c  -------------------
      go to 900
c
c  ******************************************************************
c  ------------------- grid description for spherical harmonic coeff.
c
c  ------------------- byte 7-8     j pentagonal resolution parameter
 3000 continue
      call gbyte (msga,kgds(2),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 9-10    k pentagonal resolution parameter
      call gbyte (msga,kgds(3),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 11-12   m pentagonal resolution parameter
      call gbyte (msga,kgds(4),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 13 representation type
      call gbyte (msga,kgds(5),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 14 coefficient storage mode
      call gbyte (msga,kgds(6),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  -------------------        empty fields - bytes 15 - 32
c                 set to start of byte 33
      kptr(8)  = kptr(8) + 18 * 8
      go to 900
c  ******************************************************************
c                      process mercator grids
c
c  ------------------- byte 7-8     nr of points along latitude circle
 4000 continue
      call gbyte (msga,kgds(2),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 9-10    nr of points along long meridian
      call gbyte (msga,kgds(3),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 11-13   latitue of origin
      call gbyte (msga,kgds(4),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(4),8388608).ne.0) then
          kgds(4)  =  - iand(kgds(4),8388607)
      end if
c  ------------------- byte 14-16   longitude of origin
      call gbyte (msga,kgds(5),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(5),8388608).ne.0) then
          kgds(5)  =  - iand(kgds(5),8388607)
      end if
c  ------------------- byte 17      resolution flag
      call gbyte (msga,kgds(6),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 18-20   latitude of extreme point
      call gbyte (msga,kgds(7),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(7),8388608).ne.0) then
          kgds(7)  =  - iand(kgds(7),8388607)
      end if
c  ------------------- byte 21-23   longitude of extreme point
      call gbyte (msga,kgds(8),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(8),8388608).ne.0) then
          kgds(8)  =  - iand(kgds(8),8388607)
      end if
c  ------------------- byte 24-26   latitude of projection intersection
      call gbyte (msga,kgds(9),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(9),8388608).ne.0) then
          kgds(9)  =  - iand(kgds(9),8388607)
      end if
c  ------------------- byte 27   reserved
      call gbyte (msga,kgds(10),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 28      scanning mode
      call gbyte (msga,kgds(11),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 29-31   longitudinal dir increment
      call gbyte (msga,kgds(12),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(12),8388608).ne.0) then
          kgds(12)  =  - iand(kgds(12),8388607)
      end if
c  ------------------- byte 32-34   latitudinal dir increment
      call gbyte (msga,kgds(13),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(13),8388608).ne.0) then
          kgds(13)  =  - iand(kgds(13),8388607)
      end if
c  ------------------- byte 35-42   reserved
c                        skip to start of byte 43
      kptr(8)  = kptr(8) + 8 * 8
c  -------------------
      go to 900
c  ******************************************************************
c                      process lambert conformal
c
c  ------------------- byte 7-8     nr of points along x-axis
 5000 continue
      call gbyte (msga,kgds(2),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 9-10    nr of points along y-axis
      call gbyte (msga,kgds(3),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  ------------------- byte 11-13   latitude of origin
      call gbyte (msga,kgds(4),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(4),8388608).ne.0) then
          kgds(4)  =  - iand(kgds(4),8388607)
      end if
c  ------------------- byte 14-16   longitude of origin (lower left)
      call gbyte (msga,kgds(5),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(5),8388608).ne.0) then
          kgds(5)  = - iand(kgds(5),8388607)
      end if
c  ------------------- byte 17      resolution
      call gbyte (msga,kgds(6),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 18-20   lov -orientation of grid
      call gbyte (msga,kgds(7),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(7),8388608).ne.0) then
          kgds(7)  = - iand(kgds(7),8388607)
      end if
c  ------------------- byte 21-23   dx - x-dir increment
      call gbyte (msga,kgds(8),kptr(8),24)
      kptr(8)  = kptr(8) + 24
c  ------------------- byte 24-26   dy - y-dir increment
      call gbyte (msga,kgds(9),kptr(8),24)
      kptr(8)  = kptr(8) + 24
c  ------------------- byte 27       projection center flag
      call gbyte (msga,kgds(10),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 28      scanning mode
      call gbyte (msga,kgds(11),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c  ------------------- byte 29-31   latin1 - 1st lat from pole
      call gbyte (msga,kgds(12),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(12),8388608).ne.0) then
          kgds(12)  =  - iand(kgds(12),8388607)
      end if
c  ------------------- byte 32-34   latin2 - 2nd lat from pole
      call gbyte (msga,kgds(13),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(13),8388608).ne.0) then
          kgds(13)  =  - iand(kgds(13),8388607)
      end if
c  ------------------- byte 35-37   latitude of southern pole
      call gbyte (msga,kgds(14),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(14),8388608).ne.0) then
          kgds(14)  =  - iand(kgds(14),8388607)
      end if
c  ------------------- byte 38-40   longitude of southern pole
      call gbyte (msga,kgds(15),kptr(8),24)
      kptr(8)  = kptr(8) + 24
      if (iand(kgds(15),8388608).ne.0) then
          kgds(15)  =  - iand(kgds(15),8388607)
      end if
c  ------------------- byte 41-42   reserved
      call gbyte (msga,kgds(16),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c  -------------------
  900 continue
c
c                        more code for grids with pl
c
      if (kgds(19).eq.0.or.kgds(19).eq.255) then
        if (kgds(20).ne.255) then
          isum  = 0
          kptr(8)  = nsave + (kgds(20) - 1) * 8
          call gbytes (msga,kgds(22),kptr(8),16,0,kgds(3))
          do 910 j = 1, kgds(3)
              isum  = isum + kgds(21+j)
  910     continue
          kgds(21)  = isum
        end if
      end if
      return
      end
      subroutine fi634(msga,kptr,kpds,kgds,kbms,kret)
c$$$  subprogram documentation  block
c                .      .    .                                       .
c subprogram:    fi634       extract or generate bit map for output
c   prgmmr: bill cavanaugh   org: w/nmc42    date: 91-09-13
c
c abstract: if bit map sec   is available in grib message, extract
c   for program use, otherwise generate an appropriate bit map.
c
c program history log:
c   91-09-13  cavanaugh
c   91-11-12  cavanaugh   modified size of ecmwf grids 5 - 8.
c   95-10-31  iredell     removed saves and prints
c   97-02-12  w bostelman corrects ecmwf us grid 2 processing
c   97-09-19  iredell     vectorized bitmap decoder
c   98-09-02  gilbert     corrected error in map size for u.s. grid 92
c   98-09-08  baldwin     add grids 190,192
c   99-01-20  baldwin     add grids 236,237
c   01-10-02  rogers      redefined grid #218 for 12 km eta
c                         redefined grid 192 for new 32-km eta grid
c 2003-06-30  gilbert      added grids 145 and 146 for cmaq
c                          and grid 175 for awips over guam.
c 2004-09-02  vuong       added awips grids 147, 148, 173 and 254
c
c usage:    call fi634(msga,kptr,kpds,kgds,kbms,kret)
c   input argument list:
c     msga       - bufr message
c     kptr       - array containing storage for following parameters
c          (1)   - total length of grib message
c          (2)   - length of indicator (section  0)
c          (3)   - length of pds       (section  1)
c          (4)   - length of gds       (section  2)
c          (5)   - length of bms       (section  3)
c          (6)   - length of bds       (section  4)
c          (7)   - value of current byte
c          (8)   - bit pointer
c          (9)   - grib start bit nr
c         (10)   - grib/grid element count
c         (11)   - nr unused bits at end of section 3
c         (12)   - bit map flag
c         (13)   - nr unused bits at end of section 2
c         (14)   - bds flags
c         (15)   - nr unused bits at end of section 4
c     kpds     - array containing pds elements.
c          (1)   - id of center
c          (2)   - model identification
c          (3)   - grid identification
c          (4)   - gds/bms flag
c          (5)   - indicator of parameter
c          (6)   - type of level
c          (7)   - height/pressure , etc of level
c          (8)   - year of century
c          (9)   - month of year
c          (10)  - day of month
c          (11)  - hour of day
c          (12)  - minute of hour
c          (13)  - indicator of forecast time unit
c          (14)  - time range 1
c          (15)  - time range 2
c          (16)  - time range flag
c          (17)  - number included in average
c
c   output argument list:
c     kbms       - bitmap describing location of output elements.
c     kptr       - array containing storage for following parameters
c                  see input list
c     kret       - error return
c
c remarks:
c     kret   = 0 - no error
c            = 5 - grid not avail for center indicated
c            =10 - incorrect center indicator
c            =12 - bytes 5-6 are not zero in bms, predefined bit map
c                  not provided by this center
c
c   subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 77
c   machine:  hds9000
c
c$$$
c
c                       incoming message holder
      character*1   msga(*)
c
c                       bit map
      logical*1     kbms(*)
c
c                       array of pointers and counters
      integer       kptr(*)
c                       array of pointers and counters
      integer       kpds(*)
      integer       kgds(*)
c
      integer       kret
      integer       mask(8)
c  ----------------------grid 21 and grid 22 are the same
      logical*1     grd21( 1369)
c  ----------------------grid 23 and grid 24 are the same
      logical*1     grd23( 1369)
      logical*1     grd25( 1368)
      logical*1     grd26( 1368)
c  ----------------------grid 27 and grid 28 are the same
c  ----------------------grid 29 and grid 30 are the same
c  ----------------------grid 33 and grid 34 are the same
      logical*1     grd50( 1188)
c  -----------------------grid 61 and grid 62 are the same
      logical*1     grd61( 4186)
c  -----------------------grid 63 and grid 64 are the same
      logical*1     grd63( 4186)
c     logical*1     grd70(16380)/16380*.true./
c  -------------------------------------------------------------
      data  grd21 /1333*.true.,36*.false./
      data  grd23 /.true.,36*.false.,1332*.true./
      data  grd25 /1297*.true.,71*.false./
      data  grd26 /.true.,71*.false.,1296*.true./
      data  grd50/
c line 1-4
     &  7*.false.,22*.true.,14*.false.,22*.true.,
     & 14*.false.,22*.true.,14*.false.,22*.true.,7*.false.,
c line 5-8
     &  6*.false.,24*.true.,12*.false.,24*.true.,
     & 12*.false.,24*.true.,12*.false.,24*.true.,6*.false.,
c line 9-12
     &  5*.false.,26*.true.,10*.false.,26*.true.,
     & 10*.false.,26*.true.,10*.false.,26*.true.,5*.false.,
c line 13-16
     &  4*.false.,28*.true., 8*.false.,28*.true.,
     &  8*.false.,28*.true., 8*.false.,28*.true.,4*.false.,
c line 17-20
     &  3*.false.,30*.true., 6*.false.,30*.true.,
     &  6*.false.,30*.true., 6*.false.,30*.true.,3*.false.,
c line 21-24
     &  2*.false.,32*.true., 4*.false.,32*.true.,
     &  4*.false.,32*.true., 4*.false.,32*.true.,2*.false.,
c line 25-28
     &    .false.,34*.true., 2*.false.,34*.true.,
     &  2*.false.,34*.true., 2*.false.,34*.true.,  .false.,
c line 29-33
     &           180*.true./
      data  grd61 /4096*.true.,90*.false./
      data  grd63 /.true.,90*.false.,4095*.true./
      data  mask  /128,64,32,16,8,4,2,1/
c
c     print *,'fi634'
      if (iand(kpds(4),64).eq.64) then
c
c                   set up bit pointer
c                          section 0    section 1     section 2
      kptr(8) = kptr(9) + (kptr(2)*8) + (kptr(3)*8) + (kptr(4)*8) + 24
c
c  byte 4           number of unused bits at end of section 3
c
      call gbyte (msga,kptr(11),kptr(8),8)
      kptr(8)  = kptr(8) + 8
c
c  byte 5,6         table reference if 0, bit map follows
c
      call gbyte (msga,kptr(12),kptr(8),16)
      kptr(8)  = kptr(8) + 16
c                   if table reference = 0, extract bit map
        if (kptr(12).eq.0) then
c                   calculate nr of bits in bit map
          ibits   = (kptr(5) - 6) * 8 - kptr(11)
          kptr(10)  = ibits
          if (kpds(3).eq.21.or.kpds(3).eq.22.or.kpds(3).eq.25.
     *             or.kpds(3).eq.61.or.kpds(3).eq.62) then
c                    northern hemisphere  21, 22, 25, 61, 62
              call fi634x(ibits,kptr(8),msga,kbms)
              if (kpds(3).eq.25) then
                  kadd     = 71
              else if (kpds(3).eq.61.or.kpds(3).eq.62) then
                  kadd     = 90
              else
                  kadd     = 36
              end if
              do 25 i = 1, kadd
                  kbms(i+ibits)  = .false.
   25         continue
              kptr(10)   = kptr(10) + kadd
              return
          else if (kpds(3).eq.23.or.kpds(3).eq.24.or.kpds(3).eq.26.
     *             or.kpds(3).eq.63.or.kpds(3).eq.64) then
c                    southern hemisphere  23, 24, 26, 63, 64
              call fi634x(ibits,kptr(8),msga,kbms)
              if (kpds(3).eq.26) then
                  kadd     = 72
              else if (kpds(3).eq.63.or.kpds(3).eq.64) then
                  kadd     = 91
              else
                  kadd     = 37
              end if
              do 26 i = 1, kadd
                  kbms(i+ibits)  = .false.
   26         continue
              kptr(10)   = kptr(10) + kadd - 1
              return
          else if (kpds(3).eq.50) then
              kpad    = 7
              kin     = 22
              kbits   = 0
              do 55 i = 1, 7
                  do 54 j = 1, 4
                      do 51 k = 1, kpad
                          kbits   = kbits + 1
                          kbms(kbits)  = .false.
   51                 continue
                      call fi634x(kin,kptr(8),msga,kbms(kbits+1))
                      kptr(8)=kptr(8)+kin
                      kbits=kbits+kin
                      do 53 k = 1, kpad
                          kbits   = kbits + 1
                          kbms(kbits)  = .false.
   53                 continue
   54             continue
                  kin    = kin + 2
                  kpad   = kpad - 1
   55         continue
              do 57 ii = 1, 5
                  call fi634x(kin,kptr(8),msga,kbms(kbits+1))
                  kptr(8)=kptr(8)+kin
                  kbits=kbits+kin
   57         continue
          else
c                        extract bit map from bms for other grids
              call fi634x(ibits,kptr(8),msga,kbms)
          end if
          return
        else
c         print *,'fi634-no predefined bit map provided by this center'
          kret = 12
          return
        end if
c
      end if
      kret = 0
c  -------------------------------------------------------
c                   process non-standard grid
c  -------------------------------------------------------
      if (kpds(3).eq.255) then
c         print *,'non standard grid, center = ',kpds(1)
          j      = kgds(2) * kgds(3)
          kptr(10) = j
          do 600 i = 1, j
              kbms(i) = .true.
  600     continue
          return
      end if
c  -------------------------------------------------------
c                   check international set
c  -------------------------------------------------------
      if (kpds(3).eq.21.or.kpds(3).eq.22) then
c                   ----- int'l grids 21, 22 - map size 1369
          j   = 1369
          kptr(10)  = j
          call fi637(j,kpds,kgds,kret)
          if(kret.ne.0) go to 820
          do 3021 i = 1, 1369
              kbms(i) = grd21(i)
 3021     continue
          return
      else if (kpds(3).eq.23.or.kpds(3).eq.24) then
c                   ----- int'l grids 23, 24 - map size 1369
          j   = 1369
          kptr(10)  = j
          call fi637(j,kpds,kgds,kret)
          if(kret.ne.0) go to 820
          do 3023 i = 1, 1369
              kbms(i) = grd23(i)
 3023     continue
          return
      else if (kpds(3).eq.25) then
c                   ----- int'l grid 25 - map size 1368
          j   = 1368
          kptr(10)  = j
          call fi637(j,kpds,kgds,kret)
          if(kret.ne.0) go to 820
          do 3025 i = 1, 1368
              kbms(i) = grd25(i)
 3025     continue
          return
      else if (kpds(3).eq.26) then
c                  ----- int'l grid  26 - map size 1368
          j   = 1368
          kptr(10)  = j
          call fi637(j,kpds,kgds,kret)
          if(kret.ne.0) go to 820
          do 3026 i = 1, 1368
              kbms(i) = grd26(i)
 3026     continue
          return
      else if (kpds(3).ge.37.and.kpds(3).le.44) then
c                  ----- int'l grid  37-44 - map size 3447
          j   = 3447
          go to 800
      else if (kpds(1).eq.7.and.kpds(3).eq.50) then
c                   ----- int'l grids 50 - map size 964
          j     = 1188
          kptr(10)  = j
          call fi637(j,kpds,kgds,kret)
          if(kret.ne.0) go to 890
          do 3050 i = 1, j
              kbms(i) = grd50(i)
 3050     continue
          return
      else if (kpds(3).eq.61.or.kpds(3).eq.62) then
c                   ----- int'l grids 61, 62 - map size 4186
          j     = 4186
          kptr(10)  = j
          call fi637(j,kpds,kgds,kret)
          if(kret.ne.0) go to 820
          do 3061 i = 1, 4186
              kbms(i) = grd61(i)
 3061     continue
          return
      else if (kpds(3).eq.63.or.kpds(3).eq.64) then
c                  ----- int'l grids 63, 64 - map size 4186
          j     = 4186
          kptr(10)  = j
          call fi637(j,kpds,kgds,kret)
          if(kret.ne.0) go to 820
          do 3063 i = 1, 4186
              kbms(i) = grd63(i)
 3063     continue
          return
      end if
c  -------------------------------------------------------
c                   check united states set
c  -------------------------------------------------------
      if (kpds(1).eq.7) then
          if (kpds(3).lt.100) then
              if (kpds(3).eq.1) then
c                       ----- u.s. grid 1 - map size 1679
                  j   = 1679
                  go to 800
              end if
              if (kpds(3).eq.2) then
c                       ----- u.s. grid 2 - map size 10512
                  j   = 10512
                  go to 800
              else if (kpds(3).eq.3) then
c                       ----- u.s. grid 3 - map size 65160
                  j   = 65160
                  go to 800
              else if (kpds(3).eq.4) then
c                       ----- u.s. grid 4 - map size 259920
                  j   = 259920
                  go to 800
              else if (kpds(3).eq.5) then
c                       ----- u.s. grid 5 - map size 3021
                  j   = 3021
                  go to 800
              else if (kpds(3).eq.6) then
c                       ----- u.s. grid 6 - map size 2385
                  j   = 2385
                  go to 800
              else if (kpds(3).eq.8) then
c                       ----- u.s. grid 8 - map size 5104
                  j   = 5104
                  go to 800
              else if (kpds(3).eq.27.or.kpds(3).eq.28) then
c                       ----- u.s. grids 27, 28 - map size 4225
                  j     = 4225
                  go to 800
              else if (kpds(3).eq.29.or.kpds(3).eq.30) then
c                       ----- u.s. grids 29,30 - map size 5365
                  j     = 5365
                  go to 800
              else if (kpds(3).eq.33.or.kpds(3).eq.34) then
c                       ----- u.s grid 33, 34 - map size 8326
                  j     = 8326
                  go to 800
              else if (kpds(3).ge.37.and.kpds(3).le.44) then
c                  -----  u.s. grid  37-44 - map size 3447
                  j   = 3447
                  go to 800
              else if (kpds(3).eq.45) then
c                  ----- u.s.  grid  45    - map size 41760
                  j   = 41760
                  go to 800
              else if (kpds(3).eq.53) then
c                  ----- u.s.  grid  53    - map size 5967
                  j   = 5967
                  go to 800
              else if (kpds(3).eq.55.or.kpds(3).eq.56) then
c                       ----- u.s grid 55, 56 - map size 6177
                  j     = 6177
                  go to 800
              else if (kpds(3).ge.67.and.kpds(3).le.71) then
c                       ----- u.s grid 67-71 - map size 13689
                  j     = 13689
                  go to 800
              else if (kpds(3).eq.72) then
c                       ----- u.s grid    72 - map size 406
                  j     = 406
                  go to 800
              else if (kpds(3).eq.73) then
c                       ----- u.s grid    73 - map size 13056
                  j     = 13056
                  go to 800
              else if (kpds(3).eq.74) then
c                       ----- u.s grid    74 - map size 10800
                  j     = 10800
                  go to 800
              else if (kpds(3).ge.75.and.kpds(3).le.77) then
c                       ----- u.s grid 75-77 - map size 12321
                  j     = 12321
                  go to 800
              else if (kpds(3).eq.85.or.kpds(3).eq.86) then
c                       ----- u.s grid 85,86 - map size 32400
                  j     = 32400
                  go to 800
              else if (kpds(3).eq.87) then
c                       ----- u.s grid 87     - map size 5022
                  j     = 5022
                  go to 800
              else if (kpds(3).eq.88) then
c                       ----- u.s grid 88     - map size 317840
                  j     = 317840
                  go to 800
              else if (kpds(3).eq.90) then
c                       ----- u.s grid 90     - map size 111723
                  j     = 111723
                  go to 800
              else if (kpds(3).eq.91) then
c                       ----- u.s grid 91     - map size 111723
                  j     = 111723
                  go to 800
              else if (kpds(3).eq.92) then
c                       ----- u.s grid 92     - map size 111723
                  j     = 111723
                  go to 800
              else if (kpds(3).eq.93) then
c                       ----- u.s grid 93     - map size 111723
                  j     = 111723
                  go to 800
              else if (kpds(3).eq.94) then
c                       ----- u.s grid 94     - map size 196305
                  j     = 196305
                  go to 800
              else if (kpds(3).eq.95) then
c                       ----- u.s grid 95     - map size 36062
                  j     = 36062
                  go to 800
              else if (kpds(3).eq.96) then
c                       ----- u.s grid 96     - map size 646602
                  j     = 646602
                  go to 800
              else if (kpds(3).eq.97) then
c                       ----- u.s grid 97     - map size 12727
                  j     = 12727
                  go to 800
              else if (kpds(3).eq.98) then
c                       ----- u.s grid 98     - map size 18048
                  j     = 18048
                  go to 800
              end if
          else if (kpds(3).ge.100.and.kpds(3).lt.200) then
              if (kpds(3).eq.100) then
c                       ----- u.s. grid 100 - map size 6889
                  j     = 6889
                  go to 800
              else if (kpds(3).eq.101) then
c                    ----- u.s. grid 101 - map size 10283
                  j     = 10283
                  go to 800
              else if (kpds(3).eq.103) then
c                     ----- u.s. grid 103 - map size 3640
                  j     = 3640
                  go to 800
              else if (kpds(3).eq.104) then
c                     ----- u.s. grid 104 - map size 16170
                  j     = 16170
                  go to 800
              else if (kpds(3).eq.105) then
c                 ----- u.s. grid 105 - map size 6889
                  j     = 6889
                  go to 800
              else if (kpds(3).eq.106) then
c                     ----- u.s. grid 106 - map size 19305
                  j     = 19305
                  go to 800
              else if (kpds(3).eq.107) then
c                 ----- u.s. grid 107 - map size 11040
                  j     = 11040
                  go to 800
              else if (kpds(3).eq.110) then
c                 ----- u.s. grid 110 - map size 103936
                  j     = 103936
                  go to 800
              else if (kpds(3).eq.126) then
c                 ----- u.s. grid 126 - map size 72960
                  j     = 72960
                  go to 800
              else if (kpds(3).eq.127) then
c                 ----- u.s. grid 127 - map size 294912
                  j     = 294912
                  go to 800
              else if (kpds(3).eq.130) then
c                 ----- u.s. grid 130 - map size 151987
                  j     = 151987
                  go to 800
              else if (kpds(3).eq.145) then
c                 ----- u.s. grid 145 - map size 24505
                  j     = 24505
                  go to 800
              else if (kpds(3).eq.146) then
c                 ----- u.s. grid 146 - map size 23572
                  j     = 23572
                  go to 800
              else if (kpds(3).eq.147) then
c                 ----- u.s. grid 147 - map size 69412
                  j     = 69412
                  go to 800
              else if (kpds(3).eq.148) then
c                 ----- u.s. grid 148 - map size 117130
                  j     = 117130
                  go to 800
              else if (kpds(3).eq.160) then
c                 ----- u.s. grid 160 - map size 28080 
                  j     = 28080
                  go to 800
              else if (kpds(3).eq.161) then
c                 ----- u.s. grid 161 - map size 13974 
                  j     = 13974 
                  go to 800
              else if (kpds(3).eq.163) then
c                 ----- u.s. grid 163 - map size 727776
                  j     = 727776
                  go to 800

              else if (kpds(3).eq.170) then
c                 ----- u.s. grid 170 - map size 131072 
                  j     =  131072
                  go to 800
              else if (kpds(3).eq.171) then
c                 ----- u.s. grid 171 - map size 716100
                  j     = 716100
                  go to 800
              else if (kpds(3).eq.172) then
c                 ----- u.s. grid 172 - map size 489900
                  j     = 489900
                  go to 800
              else if (kpds(3).eq.173) then
c                 ----- u.s. grid 173 - map size 9331200
                  j     = 9331200
                  go to 800
              else if (kpds(3).eq.174) then
c                 ----- u.s. grid 174 - map size 4147200
                  j     = 4147200
                  go to 800
              else if (kpds(3).eq.175) then
c                 ----- u.s. grid 175 - map size 185704
                  j     = 185704
                  go to 800
              else if (kpds(3).eq.190) then
c                 ----- u.s grid 190  - map size 12972
                  j     = 12972
                  go to 800
              else if (kpds(3).eq.192) then
c                 ----- u.s grid 192  - map size 91719
                  j     = 91719
                  go to 800
              else if (kpds(3).eq.194) then
c                 ----- u.s grid 194  - map size 12727
                  j     = 12727
                  go to 800
              else if (kpds(3).eq.196) then
c                 ----- u.s. grid 196 - map size 45903
                  j     = 45903
                  go to 800
              else if (kpds(3).eq.198) then
c                 ----- u.s. grid 198 - map size 41760
                  j     = 41760
                  go to 800
              else if (iand(kpds(4),128).eq.128) then
c                     ----- u.s. non-standard grid
                  go to 895
              end if
          else if (kpds(3).ge.200) then
              if (kpds(3).eq.201) then
                  j = 4225
                  go to 800
              else if (kpds(3).eq.202) then
                  j = 2795
                  go to 800
              else if (kpds(3).eq.203.or.kpds(3).eq.205) then
                  j = 1755
                  go to 800
              else if (kpds(3).eq.204) then
                  j = 6324
                  go to 800
              else if (kpds(3).eq.206) then
                  j = 2091
                  go to 800
              else if (kpds(3).eq.207) then
                  j = 1715
                  go to 800
              else if (kpds(3).eq.208) then
                  j = 783
                  go to 800
              else if (kpds(3).eq.209) then
                  j = 61325
                  go to 800
              else if (kpds(3).eq.210) then
                  j = 625
                  go to 800
              else if (kpds(3).eq.211) then
                  j = 6045
                  go to 800
              else if (kpds(3).eq.212) then
                  j = 23865
                  go to 800
              else if (kpds(3).eq.213) then
                  j = 10965
                  go to 800
              else if (kpds(3).eq.214) then
                  j = 6693
                  go to 800
              else if (kpds(3).eq.215) then
                  j = 94833
                  go to 800
              else if (kpds(3).eq.216) then
                  j = 14873
                  go to 800
              else if (kpds(3).eq.217) then
                  j = 59001
                  go to 800
              else if (kpds(3).eq.218) then
                  j = 262792
                  go to 800
              else if (kpds(3).eq.219) then
                  j = 179025
                  go to 800
              else if (kpds(3).eq.220) then
                  j = 122475
                  go to 800
              else if (kpds(3).eq.221) then
                  j = 96673
                  go to 800
              else if (kpds(3).eq.222) then
                  j = 15456
                  go to 800
              else if (kpds(3).eq.223) then
                  j = 16641
                  go to 800
              else if (kpds(3).eq.224) then
                  j = 4225
                  go to 800
              else if (kpds(3).eq.225) then
                  j = 24975
                  go to 800
              else if (kpds(3).eq.226) then
                  j = 381029
                  go to 800
              else if (kpds(3).eq.227) then
                  j = 1509825
                  go to 800
              else if (kpds(3).eq.228) then
                  j = 10512
                  go to 800
              else if (kpds(3).eq.229) then
                  j = 65160
                  go to 800
              else if (kpds(3).eq.230) then
                  j = 259920
                  go to 800
              else if (kpds(3).eq.231) then
                  j = 130320
                  go to 800
              else if (kpds(3).eq.232) then
                  j = 32760
                  go to 800
              else if (kpds(3).eq.233) then
                  j = 45216
                  go to 800
              else if (kpds(3).eq.234) then
                  j = 16093
                  go to 800
              else if (kpds(3).eq.235) then
                  j = 259200
                  go to 800
              else if (kpds(3).eq.236) then
                  j = 17063
                  go to 800
              else if (kpds(3).eq.237) then
                  j = 2538
                  go to 800
              else if (kpds(3).eq.238) then
                  j = 55825
                  go to 800
              else if (kpds(3).eq.239) then
                  j = 19065
                  go to 800
              else if (kpds(3).eq.240) then
                  j = 987601
                  go to 800
              else if (kpds(3).eq.241) then
                  j = 244305
                  go to 800
              else if (kpds(3).eq.242) then
                  j = 235025
                  go to 800
              else if (kpds(3).eq.243) then
                  j = 12726
                  go to 800
              else if (kpds(3).eq.244) then
                  j = 55825
                  go to 800
              else if (kpds(3).eq.245) then
                  j = 124992
                  go to 800
              else if (kpds(3).eq.246) then
                  j = 123172
                  go to 800
              else if (kpds(3).eq.247) then
                  j = 124992
                  go to 800
              else if (kpds(3).eq.248) then
                  j = 13635
                  go to 800
              else if (kpds(3).eq.249) then
                  j = 125881
                  go to 800
              else if (kpds(3).eq.250) then
                  j = 13635
                  go to 800
              else if (kpds(3).eq.251) then
                  j = 69720
                  go to 800
              else if (kpds(3).eq.252) then
                  j = 67725
                  go to 800
              else if (kpds(3).eq.253) then
                  j = 83552
                  go to 800
              else if (kpds(3).eq.254) then
                  j = 110700
                  go to 800
              else if (iand(kpds(4),128).eq.128) then
                  go to 895
              end if
              kret  = 5
              return
          end if
      end if
c  -------------------------------------------------------
c                   check japan meteorological agency set
c  -------------------------------------------------------
      if (kpds(1).eq.34) then
          if (iand(kpds(4),128).eq.128) then
c             print *,'jma map is not predefined, the gds will'
c             print *,'be used to unpack the data, map = ',kpds(3)
              go to 900
          end if
      end if
c  -------------------------------------------------------
c                   check canadian set
c  -------------------------------------------------------
      if (kpds(1).eq.54) then
          if (iand(kpds(4),128).eq.128) then
c             print *,'canadian map is not predefined, the gds will'
c             print *,'be used to unpack the data, map = ',kpds(3)
              go to 900
          end if
      end if
c  -------------------------------------------------------
c                   check fnoc set
c  -------------------------------------------------------
      if (kpds(1).eq.58) then
          if (kpds(3).eq.220.or.kpds(3).eq.221) then
c                      fnoc grid 220, 221 - mapsize 3969 (63 * 63)
              j  = 3969
              kptr(10)  = j
              do i = 1, j
                  kbms(i)  = .true.
              end do
              return
          end if
          if (kpds(3).eq.223) then
c                      fnoc grid 223 - mapsize 10512 (73 * 144)
              j  = 10512
              kptr(10)  = j
              do i = 1, j
                  kbms(i)  = .true.
              end do
              return
          end if
          if (iand(kpds(4),128).eq.128) then
c             print *,'fnoc map is not predefined, the gds will'
c             print *,'be used to unpack the data, map = ',kpds(3)
              go to 900
          end if
      end if
c  -------------------------------------------------------
c                   check ukmet set
c  -------------------------------------------------------
      if (kpds(1).eq.74) then
          if (iand(kpds(4),128).eq.128) then
              go to 820
          end if
      end if
c  -------------------------------------------------------
c                   check ecmwf set
c  -------------------------------------------------------
      if (kpds(1).eq.98) then
          if (kpds(3).ge.1.and.kpds(3).le.12) then
              if (kpds(3).ge.5.and.kpds(3).le.8) then
                  j     = 1073
              else
                  j     = 1369
              end if
              kptr(10)  = j
              call fi637(j,kpds,kgds,kret)
              if(kret.ne.0) go to 810
              kptr(10)  = j  ! reset for modified j
              do 1000 i = 1, j
                  kbms(i) = .true.
 1000         continue
              return
          else if (kpds(3).ge.13.and.kpds(3).le.16) then
              j         = 361
              kptr(10)  = j
              call fi637(j,kpds,kgds,kret)
              if(kret.ne.0) go to 810
              do 1013 i = 1, j
                  kbms(i) = .true.
 1013         continue
              return
          else if (iand(kpds(4),128).eq.128) then
                  go to 810
          else
              kret  = 5
              return
          end if
      else
c         print *,'center ',kpds(1),' is not defined'
          if (iand(kpds(4),128).eq.128) then
c             print *,'gds will be used to unpack the data',
c    *                        ' map = ',kpds(3)
              go to 900
          else
              kret  = 10
              return
          end if
      end if
c =======================================
c
  800 continue
      kptr(10)  = j
      call fi637 (j,kpds,kgds,kret)
      if(kret.ne.0) go to 801
      do 2201 i = 1, j
          kbms(i)  = .true.
 2201 continue
      return
  801 continue
c
c  ----- the map has a gds, byte 7 of the (pds) the grid identification
c  ----- is not 255, the size of the grid is not the same as the
c  ----- predefined sizes of the u.s. grids, or known grids of the
c  ----- of the other centers. the grid can be unknown, or from an
c  ----- unknown center, we will use the information in the gds to make
c  ----- a bit map.
c
  810 continue
c     print *,'ecmwf predefined map size does not match, i will use'
      go to 895
c
  820 continue
c     print *,'u.k. met predefined map size does not match, i will use'
      go to 895
c
  890 continue
c     print *,'predefined map size does not match, i will use'
  895 continue
c     print *,'the gds to unpack the data, map type = ',kpds(3)
c
  900 continue
        j      = kgds(2) * kgds(3)
c                    afos afos afos        special case
c                             involves next single statement only
        if (kpds(3).eq.211) kret = 0
        kptr(10) = j
        do 2203 i = 1, j
          kbms(i) = .true.
 2203   continue
c     print *,'exit fi634'
      return
      end
c-----------------------------------------------------------------------
      subroutine fi634x(npts,nskp,msga,kbms)
c$$$  subprogram documentation  block
c                .      .    .                                       .
c subprogram:    fi634x      extract bit map
c   prgmmr: iredell          org: w/np23     date: 91-09-19
c
c abstract: extract the packed bitmap into a logical array.
c
c program history log:
c   97-09-19  iredell     vectorized bitmap decoder
c
c usage:    call fi634x(npts,nskp,msga,kbms)
c   input argument list:
c     npts       - integer number of points in the bitmap field
c     nskp       - integer number of bits to skip in grib message
c     msga       - character*1 grib message
c
c   output argument list:
c     kbms       - logical*1 bitmap
c
c remarks:
c   subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 77
c   machine:  cray
c
c$$$
      character*1   msga(*)
      logical*1     kbms(npts)
      integer       ichk(npts)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call gbytes(msga,ichk,nskp,1,0,npts)
      kbms=ichk.ne.0
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
      subroutine fi635(msga,kptr,kpds,kgds,kbms,data,kret)
c$$$  subprogram documentation  block
c                .      .    .                                       .
c subprogram:    fi635         extract grib data elements from bds
c   prgmmr: bill cavanaugh   org: w/nmc42    date: 91-09-13
c
c abstract: extract grib data from binary data section and place
c           into output array in proper position.
c
c program history log:
c   91-09-13  cavanaugh
c   94-04-01  cavanaugh  modified code to include decimal scaling when
c                        calculating the value of data points specified
c                        as being equal to the reference value
c   94-11-10  farley     increased mxsize from 72960 to 260000
c                        for .5 degree sst analysis fields
c   95-10-31  iredell    removed saves and prints
c   98-08-31  iredell    eliminated need for mxsize
c
c usage:    call fi635(msga,kptr,kpds,kgds,kbms,data,kret)
c   input argument list:
c     msga       - array containing grib message
c     kptr       - array containing storage for following parameters
c          (1)   - total length of grib message
c          (2)   - length of indicator (section  0)
c          (3)   - length of pds       (section  1)
c          (4)   - length of gds       (section  2)
c          (5)   - length of bms       (section  3)
c          (6)   - length of bds       (section  4)
c          (7)   - value of current byte
c          (8)   - bit pointer
c          (9)   - grib start bit nr
c         (10)   - grib/grid element count
c         (11)   - nr unused bits at end of section 3
c         (12)   - bit map flag
c         (13)   - nr unused bits at end of section 2
c         (14)   - bds flags
c         (15)   - nr unused bits at end of section 4
c         (16)   - reserved
c         (17)   - reserved
c         (18)   - reserved
c         (19)   - binary scale factor
c         (20)   - num bits used to pack each datum
c     kpds     - array containing pds elements.
c                  see initial routine
c     kbms       - bitmap describing location of output elements.
c
c   output argument list:
c     kbds       - information extracted from binary data section
c     kbds(1)  - n1
c     kbds(2)  - n2
c     kbds(3)  - p1
c     kbds(4)  - p2
c     kbds(5)  - bit pointer to 2nd order widths
c     kbds(6)  -  "    "     "   "   "    bit maps
c     kbds(7)  -  "    "     "  first order values
c     kbds(8)  -  "    "     "  second order values
c     kbds(9)  -  "    "     start of bds
c     kbds(10) -  "    "     main bit map
c     kbds(11) - binary scaling
c     kbds(12) - decimal scaling
c     kbds(13) - bit width of first order values
c     kbds(14) - bit map flag
c                 0 = no second order bit map
c                 1 = second order bit map present
c     kbds(15) - second order bit width
c     kbds(16) - constant / different widths
c                 0 = constant widths
c                 1 = different widths
c     kbds(17) - single datum / matrix
c                 0 = single datum at each grid point
c                 1 = matrix of values at each grid point
c       (18-20)- unused
c
c     data       - real array of gridded elements in grib message.
c     kptr       - array containing storage for following parameters
c                  see input list
c     kret       - error return
c
c remarks:
c     error return
c              3 = unpacked field is larger than 65160
c              6 = does not match nr of entries for this grib/grid
c              7 = number of bits in fill too large
c
c   subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 77
c   machine:  hds9000
c
c$$$
c
      character*1   msga(*)
      character*1   kk(8)
      character*1   ckref(8)
c
      logical*1     kbms(*)
c
      integer       kpds(*)
      integer       kgds(*)
      integer       kbds(20)
      integer       kptr(*)
      integer       nrbits
      integer       kref
      integer       kkk
      integer,allocatable::  ksave(:)
      integer       kscale
c
      real          data(*)
      real          refnce
      real          scale
      real          realkk
c
      equivalence   (ckref(1),kref,refnce)
      equivalence   (kk(1),kkk,realkk)
c
c
c     changed hex values to decimal to make code more portable
c
c  *************************************************************
c     print *,'enter fi635'
c              set up bit pointer
      kptr(8) = kptr(9) + (kptr(2)*8) + (kptr(3)*8) + (kptr(4)*8)
     *                + (kptr(5)*8) + 24
c  ------------- extract flags
c            byte 4
      call gbyte(msga,kptr(14),kptr(8),4)
      kptr(8)  = kptr(8) + 4
c  --------- nr of unused bits in section 4
      call gbyte(msga,kptr(15),kptr(8),4)
      kptr(8)  = kptr(8) + 4
      kend    = kptr(9) + (kptr(2)*8) + (kptr(3)*8) + (kptr(4)*8)
     *                + (kptr(5)*8) + kptr(6) * 8 - kptr(15)
c  ------------- get scale factor
c            bytes 5,6
c                                  check sign
      call gbyte (msga,ksign,kptr(8),1)
      kptr(8)  = kptr(8) + 1
c                                  get absolute scale value
      call gbyte (msga,kscale,kptr(8),15)
      kptr(8)  = kptr(8) + 15
      if (ksign.gt.0) then
          kscale  = - kscale
      end if
      scale = 2.0**kscale
      kptr(19)=kscale
c  ------------ get reference value
c            bytes 7,10
      call gbyte (msga,kref,kptr(8),32)
      kptr(8)  = kptr(8) + 32
c
c     the next code will convert the ibm370 floating point
c     to the floating point used on your computer.
c
c     1st test to see in on 32 or 64 bit word machine
c     lw = 4 or 8;  if 8 may be a cray
c
      call w3fi01(lw)
      if (lw.eq.4) then
        call gbyte (ckref,jsgn,0,1)
        call gbyte (ckref,jexp,1,7)
        call gbyte (ckref,ifr,8,24)
      else
        call gbyte (ckref,jsgn,32,1)
        call gbyte (ckref,jexp,33,7)
        call gbyte (ckref,ifr,40,24)
      endif
c     print *,109,jsgn,jexp,ifr
c 109 format (' jsgn,jexp,ifr = ',3(1x,z8))
      if (ifr.eq.0) then
          refnce  = 0.0
      else if (jexp.eq.0.and.ifr.eq.0) then
          refnce  = 0.0
      else
          refnce  = float(ifr) * 16.0 ** (jexp - 64 - 6)
          if (jsgn.ne.0) refnce = - refnce
      end if
c     print *,'scale ',scale,' ref val ',kref,refnce
c  ------------- number of bits specified for each entry
c            byte 11
      call gbyte (msga,kbits,kptr(8),8)
      kptr(8)  = kptr(8) + 8
      kbds(4)  = kbits
c     kbds(13) = kbits
      kptr(20) = kbits
      ibyt12   = kptr(8)
c  ------------------ if there are no extended flags present
c                     this is where data begins and and the processing
c                     included in the following if...end if
c                     will be skipped
c     print *,'basic flags =',kptr(14) ,iand(kptr(14),1)
      if (iand(kptr(14),1).eq.0) then
c         print *,'no extended flags'
      else
c            bytes 12,13
          call gbyte (msga,koctet,kptr(8),16)
          kptr(8)  = kptr(8) + 16
c  --------------------------- extended flags
c            byte 14
          call gbyte (msga,kxflag,kptr(8),8)
c         print *,'have extended flags',kxflag
          kptr(8)  = kptr(8) + 8
          if (iand(kxflag,16).eq.0) then
c                          second order values constant widths
              kbds(16)  = 0
          else
c                          second order values different widths
              kbds(16)  = 1
          end if
          if (iand (kxflag,32).eq.0) then
c                         no secondary bit map
              kbds(14)  = 0
          else
c                         have secondary bit map
              kbds(14)  = 1
          end if
          if (iand (kxflag,64).eq.0) then
c                         single datum at grid point
              kbds(17)  = 0
          else
c                         matrix of values at grid point
              kbds(17)  = 1
          end if
c  ---------------------- nr - first dimension (rows) of each matrix
c            bytes 15,16
          call gbyte (msga,nr,kptr(8),16)
          kptr(8)  = kptr(8) + 16
c  ---------------------- nc - second dimension (cols) of each matrix
c            bytes 17,18
          call gbyte (msga,nc,kptr(8),16)
          kptr(8)  = kptr(8) + 16
c  ---------------------- nrv - first dim coord vals
c            byte 19
          call gbyte (msga,nrv,kptr(8),8)
          kptr(8)  = kptr(8) + 8
c  ---------------------- nc1 - nr coeff's or values
c            byte 20
          call gbyte (msga,nc1,kptr(8),8)
          kptr(8)  = kptr(8) + 8
c  ---------------------- ncv - second dim coord or value
c            byte 21
          call gbyte (msga,ncv,kptr(8),8)
          kptr(8)  = kptr(8) + 8
c  ---------------------- nc2 - nr coeff's or vals
c            byte 22
          call gbyte (msga,nc2,kptr(8),8)
          kptr(8)  = kptr(8) + 8
c  ---------------------- kphys1 - first dim physical signif
c            byte 23
          call gbyte (msga,kphys1,kptr(8),8)
          kptr(8)  = kptr(8) + 8
c  ---------------------- kphys2 - second dim physical signif
c            byte 24
          call gbyte (msga,kphys2,kptr(8),8)
          kptr(8)  = kptr(8) + 8
c            bytes 25-n
      end if
      if (kbits.eq.0) then
c                       have no bds entries, all entries = refnce
          scal10  = 10.0 ** kpds(22)
          scal10  = 1.0 / scal10
          refn10  = refnce * scal10
          kentry = kptr(10)
          do 210 i = 1, kentry
              data(i) = 0.0
              if (kbms(i)) then
                   data(i) = refn10
              end if
  210     continue
          go to 900
      end if
c     print *,'kend ',kend,' kptr(8) ',kptr(8),'kbits ',kbits
      knr     = (kend - kptr(8)) / kbits
c     print *,'number of entries in data array',knr
c  --------------------
c       cycle thru bds until have used all (specified number)
c       entries.
c  ------------- unused bits in data area
c number of bytes in data area
      nrbyte  = kptr(6) - 11
c  ------------- total nr of usable bits
      nrbits  = nrbyte * 8  - kptr(15)
c  ------------- total nr of entries
      kentry = nrbits / kbits
c                             allocate ksave
      allocate(ksave(kentry))
c
c     if (iand(kptr(14),2).eq.0) then
c        print *,'source values in floating point'
c     else
c        print *,'source values in integer'
c     end if
c
      if (iand(kptr(14),8).eq.0) then
c        print *,'processing grid point data'
         if (iand(kptr(14),4).eq.0) then
c            print *,'    with simple packing'
             if (iand(kptr(14),1).eq.0) then
c                print *,'        with no additional flags'
                 go to 4000
             else if (iand(kptr(14),1).ne.0) then
c                print *,'        with additional flags',kxflag
                 if (kbds(17).eq.0) then
c                    print *,'            single datum each grid pt'
                     if (kbds(14).eq.0) then
c                        print *,'            no sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                          ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                            ' values different widths'
                         end if
                     else if (kbds(14).ne.0) then
c                        print *,'            sec bit map'
                         if (kbds(16).eq.0) then
c                             print *,'            second order',
c    *                              ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                             ' values different widths'
                         end if
                     end if
                 else if (kbds(17).ne.0) then
c                    print *,'            matrix of vals each pt'
                     if (kbds(14).eq.0) then
c                        print *,'            no sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                          ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                              ' values different widths'
                         end if
                     else if (kbds(14).ne.0) then
c                        print *,'            sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                             ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                              ' values different widths'
                         end if
                     end if
                 end if
             end if
         else if (iand(kptr(14),4).ne.0) then
c            print *,'    with complex/second order packing'
             if (iand(kptr(14),1).eq.0) then
c                    print *,'        with no additional flags'
             else if (iand(kptr(14),1).ne.0) then
c                print *,'        with additional flags'
                 if (kbds(17).eq.0) then
c                    print *,'            single datum at each pt'
                     if (kbds(14).eq.0) then
c                            print *,'            no sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                             ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                              ' values different widths'
                         end if
c                                       row by row - col by col
                         call fi636 (data,msga,kbms,
     *                                         refnce,kptr,kpds,kgds)
                         go to 900
                     else if (kbds(14).ne.0) then
c                        print *,'            sec bit map'
                         if (kbds(16).eq.0) then
c                                print *,'            second order',
c    *                              ' values constant width'
                         else if (kbds(16).ne.0) then
c                                print *,'            second order',
c    *                              ' values different widths'
                         end if
                         call fi636 (data,msga,kbms,
     *                                         refnce,kptr,kpds,kgds)
                         go to 900
                     end if
                 else if (kbds(17).ne.0) then
c                    print *,'            matrix of vals each pt'
                     if (kbds(14).eq.0) then
c                        print *,'            no sec bit map'
                         if (kbds(16).eq.0) then
c                              print *,'            second order',
c    *                              ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                              ' values different widths'
                         end if
                     else if (kbds(14).ne.0) then
c                        print *,'            sec bit map'
                         if (kbds(16).eq.0) then
c                              print *,'            second order',
c    *                              ' values constant width'
                         else if (kbds(16).ne.0) then
c                                print *,'            second order',
c    *                              ' values different widths'
                         end if
                     end if
                 end if
             end if
         end if
      else if (iand(kptr(14),8).ne.0) then
c        print *,'processing spherical harmonic coefficients'
         if (iand(kptr(14),4).eq.0) then
c            print *,'    with simple packing'
             if (iand(kptr(14),1).eq.0) then
c                print *,'        with no additional flags'
                 go to 5000
             else if (iand(kptr(14),1).ne.0) then
c                print *,'        with additional flags'
                 if (kbds(17).eq.0) then
c                    print *,'            single datum each grid pt'
                     if (kbds(14).eq.0) then
c                        print *,'            no sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                              ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                              ' values different widths'
                         end if
                     else if (kbds(14).ne.0) then
c                        print *,'            sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                              ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                            ' values different widths'
                         end if
                     end if
                 else if (kbds(17).ne.0) then
c                    print *,'            matrix of vals each pt'
                     if (kbds(14).eq.0) then
c                        print *,'            no sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                              ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                             ' values different widths'
                         end if
                     else if (kbds(14).ne.0) then
c                        print *,'            sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                              ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                             ' values different widths'
                         end if
                     end if
                 end if
             end if
         else if (iand(kptr(14),4).ne.0) then
c                                  complex/second order packing
c            print *,'    with complex/second order packing'
             if (iand(kptr(14),1).eq.0) then
c                print *,'        with no additional flags'
             else if (iand(kptr(14),1).ne.0) then
c                print *,'        with additional flags'
                 if (kbds(17).eq.0) then
c                    print *,'            single datum each grid pt'
                     if (kbds(14).eq.0) then
c                        print *,'            no sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                             ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                              ' values different widths'
                         end if
                     else if (kbds(14).ne.0) then
c                        print *,'            sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                              ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                              ' values different widths'
                         end if
                     end if
                 else if (kbds(17).ne.0) then
c                    print *,'            matrix of vals each pt'
                     if (kbds(14).eq.0) then
c                        print *,'            no sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                            ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                              ' values different widths'
                         end if
                     else if (kbds(14).ne.0) then
c                        print *,'            sec bit map'
                         if (kbds(16).eq.0) then
c                            print *,'            second order',
c    *                              ' values constant width'
                         else if (kbds(16).ne.0) then
c                            print *,'            second order',
c    *                              ' values different widths'
                         end if
                     end if
                 end if
             end if
         end if
      end if
      if(allocated(ksave)) deallocate(ksave)
c     print *,' not processed - not processed - not processed'
      kret   = 11
      return
 4000 continue
c  ****************************************************************
c
c grid point data, simple packing, floating point, no addn'l flags
c
      scal10  = 10.0 ** kpds(22)
      scal10  = 1.0 / scal10
      if (kpds(3).eq.23.or.kpds(3).eq.24.or.kpds(3).eq.26.
     *            or.kpds(3).eq.63.or.kpds(3).eq.64) then
          if (kpds(3).eq.26) then
              kadd    = 72
          else if (kpds(3).eq.63.or.kpds(3).eq.64) then
              kadd    = 91
          else
              kadd    = 37
          end if
          call gbytes (msga,ksave,kptr(8),kbits,0,knr)
          kptr(8)   = kptr(8) + kbits * knr
          ii        = 1
          kentry    = kptr(10)
          do 4001 i = 1, kentry
              if (kbms(i)) then
                  data(i)   = (refnce+float(ksave(ii))*scale)*scal10
                  ii        = ii + 1
              else
                  data(i)   = 0.0
              end if
 4001     continue
          do 4002 i = 2, kadd
              data(i)   = data(1)
 4002     continue
      else if (kpds(3).eq.21.or.kpds(3).eq.22.or.kpds(3).eq.25.
     *            or.kpds(3).eq.61.or.kpds(3).eq.62) then
          call gbytes (msga,ksave,kptr(8),kbits,0,knr)
          ii    = 1
          kentry = kptr(10)
          do 4011 i = 1, kentry
              if (kbms(i)) then
                  data(i) = (refnce + float(ksave(ii)) * scale) * scal10
                  ii  = ii + 1
              else
                  data(i) = 0.0
              end if
 4011     continue
          if (kpds(3).eq.25) then
              kadd    = 71
          else if (kpds(3).eq.61.or.kpds(3).eq.62) then
              kadd    = 90
          else
              kadd    = 36
          end if
          lastp   = kentry - kadd
          do 4012 i = lastp+1, kentry
              data(i) = data(lastp)
 4012     continue
      else
          call gbytes (msga,ksave,kptr(8),kbits,0,knr)
          ii    = 1
          kentry = kptr(10)
          do 500 i = 1, kentry
              if (kbms(i)) then
                  data(i) = (refnce + float(ksave(ii)) * scale) * scal10
                  ii  = ii + 1
              else
                  data(i) = 0.0
              end if
  500     continue
      end if
      go to 900
c  ------------- process spherical harmonic coefficients,
c               simple packing, floating point, no addn'l flags
 5000 continue
c     print *,'check point spectral coeff'
      kptr(8)  = ibyt12
      call gbyte (msga,kkk,kptr(8),32)
      kptr(8)  = kptr(8) + 32
c
c     the next code will convert the ibm370 foating point
c     to the floating point used on your machine.
c
c     1st test to see in on 32 or 64 bit word machine
c     lw = 4 or 8;  if 8 may be a cray
c
      call w3fi01(lw)
      if (lw.eq.4) then
        call gbyte (kk,jsgn,0,1)
        call gbyte (kk,jexp,1,7)
        call gbyte (kk,ifr,8,24)
      else
        call gbyte (kk,jsgn,32,1)
        call gbyte (kk,jexp,33,7)
        call gbyte (kk,ifr,40,24)
      endif
c
      if (ifr.eq.0) then
          realkk  = 0.0
      else if (jexp.eq.0.and.ifr.eq.0) then
          realkk  = 0.0
      else
          realkk  = float(ifr) * 16.0 ** (jexp - 64 - 6)
          if (jsgn.ne.0) realkk  = -realkk
      end if
      data(1)  = realkk
      call gbytes (msga,ksave,kptr(8),kbits,0,knr)
c  --------------
      do 6000 i = 1, kentry
          data(i+1)  = refnce + float(ksave(i)) * scale
 6000 continue
  900 continue
      if(allocated(ksave)) deallocate(ksave)
c     print *,'exit fi635'
      return
      end
      subroutine fi636 (data,msga,kbms,refnce,kptr,kpds,kgds)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    fi636       process second order packing
c   prgmmr: cavanaugh        org: w/nmc42    date: 92-09-22
c
c abstract: process second order packing from the binary data section
c   (bds) for single data items grid point data
c
c program history log:
c   93-06-08  cavanaugh
c   93-12-15  cavanaugh   modified second order pointers to first order
c                         values and second order values correctly.
c   95-04-26  r.e.jones   fi636 corection for 2nd order complex
c                         unpacking.
c   95-10-31  iredell     removed saves and prints
c
c usage:    call fi636 (data,msga,kbms,refnce,kptr,kpds,kgds)
c   input argument list:
c
c     msga     - array containing grib message
c     refnce   - reference value
c     kptr     - work array
c
c   output argument list:      (including work arrays)
c     data     - location of output array
c              working array
c     kbds(1)  - n1
c     kbds(2)  - n2
c     kbds(3)  - p1
c     kbds(4)  - p2
c     kbds(5)  - bit pointer to 2nd order widths
c     kbds(6)  -  "    "     "   "   "    bit maps
c     kbds(7)  -  "    "     "  first order values
c     kbds(8)  -  "    "     "  second order values
c     kbds(9)  -  "    "     start of bds
c     kbds(10) -  "    "     main bit map
c     kbds(11) - binary scaling
c     kbds(12) - decimal scaling
c     kbds(13) - bit width of first order values
c     kbds(14) - bit map flag
c                 0 = no second order bit map
c                 1 = second order bit map present
c     kbds(15) - second order bit width
c     kbds(16) - constant / different widths
c                 0 = constant widths
c                 1 = different widths
c     kbds(17) - single datum / matrix
c                 0 = single datum at each grid point
c                 1 = matrix of values at each grid point
c       (18-20)- unused
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 77
c   machine:  hds, cray
c
c$$$
      real         data(*)
      real         refn
      real         refnce
c
      integer      kbds(20)
      integer      kptr(*)
      integer      jref,bmap2(12500)
      integer      i,ibds
      integer      kbit,ifoval,isoval
      integer      kpds(*),kgds(*)
c
      logical*1    kbms(*)
c
      character*1  msga(*)
c
      equivalence  (jref,refn)
c  *******************     setup     ******************************
c     print *,'enter fi636'
c                                start of bms (bit pointer)
      do i = 1,20
        kbds(i)  = 0
      end do
c                byte start of bds
      ibds  = kptr(2) + kptr(3) + kptr(4) + kptr(5)
c     print *,'kptr(2-5) ',kptr(2),kptr(3),kptr(4),kptr(5)
c                bit start of bds
      jptr  = ibds * 8
c     print *,'jptr ',jptr
      kbds(9) = jptr
c     print *,'start of bds         ',kbds(9)
c                    binary scale value  bds bytes 5-6
      call gbyte (msga,isign,jptr+32,1)
      call gbyte (msga,kbds(11),jptr+33,15)
      if (isign.gt.0) then
          kbds(11)  = - kbds(11)
      end if
c     print *,'binary scale value =',kbds(11)
c                  extract reference value
      call gbyte(msga,jref,jptr+48,32)
c     print *,'decoded reference value =',refn,refnce
c                f o bit width
      call gbyte(msga,kbds(13),jptr+80,8)
      jptr  = jptr + 88
c              at start of bds byte 12
c                extract n1
      call gbyte (msga,kbds(1),jptr,16)
c     print *,'n1  = ',kbds(1)
      jptr  = jptr + 16
c                 extended flags
      call gbyte (msga,kflag,jptr,8)
c                 isolate bit map flag
      if (iand(kflag,32).ne.0) then
        kbds(14)  = 1
      else
        kbds(14)  = 0
      end if
      if (iand(kflag,16).ne.0) then
        kbds(16)  = 1
      else
        kbds(16)  = 0
      end if
      if (iand(kflag,64).ne.0) then
        kbds(17)  = 1
      else
        kbds(17)  = 0
      end if
      jptr  = jptr + 8
c                extract n2
      call gbyte (msga,kbds(2),jptr,16)
c     print *,'n2  = ',kbds(2)
      jptr  = jptr + 16
c                extract p1
      call gbyte (msga,kbds(3),jptr,16)
c     print *,'p1  = ',kbds(3)
      jptr  = jptr + 16
c                extract p2
      call gbyte (msga,kbds(4),jptr,16)
c     print *,'p2  = ',kbds(4)
      jptr  = jptr + 16
c                 skip reserved byte
      jptr    = jptr + 8
c                start of second order bit widths
      kbds(5) = jptr
c                compute start of secondary bit map
      if (kbds(14).ne.0) then
c                           for included secondary bit map
          jptr    = jptr + (kbds(3) * 8)
          kbds(6) = jptr
      else
c                           for constructed secondary bit map
          kbds(6)  = 0
      end if
c                create pointer to start of first order values
      kbds(7) =  kbds(9) + kbds(1) * 8 - 8
c     print *,'bit pointer to start of fovals',kbds(7)
c                create pointer to start of second order values
      kbds(8) =  kbds(9) + kbds(2) * 8 - 8
c     print *,'bit pointer to start of sovals',kbds(8)
c     print *,'kbds( 1) - n1                         ',kbds( 1)
c     print *,'kbds( 2) - n2                         ',kbds( 2)
c     print *,'kbds( 3) - p1                         ',kbds( 3)
c     print *,'kbds( 4) - p2                         ',kbds( 4)
c     print *,'kbds( 5) - bit ptr - 2nd order widths ',kbds( 5)
c     print *,'kbds( 6) -  "   "     "   " bit maps  ',kbds( 6)
c     print *,'kbds( 7) -  "   "     f o vals        ',kbds( 7)
c     print *,'kbds( 8) -  "   "     s o vals        ',kbds( 8)
c     print *,'kbds( 9) -  "   "    start of bds     ',kbds( 9)
c     print *,'kbds(10) -  "   "    main bit map     ',kbds(10)
c     print *,'kbds(11) - binary scaling             ',kbds(11)
c     print *,'kpds(22) - decimal scaling            ',kpds(22)
c     print *,'kbds(13) - fo bit width               ',kbds(13)
c     print *,'kbds(14) - 2nd order bit map flag     ',kbds(14)
c     print *,'kbds(15) - 2nd order bit width        ',kbds(15)
c     print *,'kbds(16) - constant/different widths  ',kbds(16)
c     print *,'kbds(17) - single datum/matrix        ',kbds(17)
c     print *,'refnce val                            ',refnce
c  ************************* process data  **********************
      ij  = 0
c  ========================================================
      if (kbds(14).eq.0) then
c                           no bit map, must construct one
          if (kgds(2).eq.65535) then
              if (kgds(20).eq.255) then
c                 print *,'cannot be used here'
              else
c                               point to pl
          lp  = kptr(9) + kptr(2)*8 + kptr(3)*8 + kgds(20)*8 - 8
c                 print *,'lp = ',lp
                  jt  = 0
                  do 2000 jz = 1, kgds(3)
c                               get number in current row
                      call gbyte (msga,number,lp,16)
c                               increment to next row number
                      lp  = lp + 16
c                     print *,'number in row',jz,' = ',number
                      do 1500 jq = 1, number
                          if (jq.eq.1) then
                              call sbyte (bmap2,1,jt,1)
                          else
                              call sbyte (bmap2,0,jt,1)
                          end if
                          jt  = jt + 1
 1500                 continue
 2000             continue
              end if
          else
              if (iand(kgds(11),32).eq.0) then
c                           row by row
c                 print *,'     row by row'
                  kout  = kgds(3)
                  kin   = kgds(2)
              else
c                           col by col
c                 print *,'     col by col'
                  kin   = kgds(3)
                  kout  = kgds(2)
              end if
c             print *,'kin=',kin,' kout= ',kout
              do 200 i = 1, kout
                  do 150 j = 1, kin
                      if (j.eq.1) then
                          call sbyte (bmap2,1,ij,1)
                      else
                          call sbyte (bmap2,0,ij,1)
                      end if
                      ij  = ij + 1
  150             continue
  200         continue
          end if
      end if
c  ========================================================
c     print 99,(bmap2(j),j=1,110)
c99   format ( 10(1x,z8.8))
c     call binary (bmap2,2)
c                for each grid point entry
c
         scale2  = 2.0**kbds(11)
         scal10  = 10.0**kpds(22)
c     print *,'scale values - ',scale2,scal10
      do 1000 i = 1, kptr(10)
c                    get next master bit map bit position
c                    if next master bit map bit position is 'on' (1)
          if (kbms(i)) then
c             write(6,900)i,kbms(i)
c 900         format (1x,i4,3x,14hmain bit is on,3x,l4)
              if (kbds(14).ne.0) then
                  call gbyte (msga,kbit,kbds(6),1)
              else
                  call gbyte (bmap2,kbit,kbds(6),1)
              end if
c             print *,'kbds(6) =',kbds(6),' kbit =',kbit
              kbds(6)  = kbds(6) + 1
              if (kbit.ne.0) then
c                 print *,'          sob on'
c                                  get next first order packed value
                  call gbyte (msga,ifoval,kbds(7),kbds(13))
                  kbds(7)  = kbds(7) + kbds(13)
c                 print *,'foval =',ifoval
c                                   get second order bit width
                  call gbyte (msga,kbds(15),kbds(5),8)
                  kbds(5)  = kbds(5) + 8
c                print *,kbds(7)-kbds(13),' foval =',ifoval,' kbds(5)=',
c    *                           ,kbds(5), 'isowid =',kbds(15)
              else
c                 print *,'          sob not on'
              end if
              isoval  = 0
              if (kbds(15).eq.0) then
c                        if second order bit width = 0
c                             then second order value is 0
c                            so calculate data value for this point
c                 data(i) = (refnce + (float(ifoval) * scale2)) / scal10
              else
                  call gbyte (msga,isoval,kbds(8),kbds(15))
                  kbds(8)  = kbds(8) + kbds(15)
              end if
              data(i) = (refnce + (float(ifoval + isoval) *
     *                         scale2)) / scal10
c             print *,i,data(i),refnce,ifoval,isoval,scale2,scal10
          else
c             write(6,901) i,kbms(i)
c 901         format (1x,i4,3x,15hmain bit not on,3x,l4)
              data(i)  = 0.0
          end if
c         print *,i,data(i),ifoval,isoval,kbds(5),kbds(15)
 1000 continue
c  **************************************************************
c     print *,'exit fi636'
      return
      end
      subroutine fi637(j,kpds,kgds,kret)
c$$$  subprogram documentation  block
c                .      .    .                                       .
c subprogram:    fi637       grib grid/size test
c   prgmmr: cavanaugh        org: w/nmc42    date: 91-09-13
c
c abstract: to test when gds is available to see if size mismatch
c   on existing grids (by center) is indicated
c
c program history log:
c   91-09-13  cavanaugh
c   95-10-31  iredell     removed saves and prints
c   97-02-12  w bostelman corrects ecmwf us grid 2 processing
c   98-06-17  iredell     removed alternate return
c   99-01-20  baldwin    modify to handle grid 237
c
c usage:    call fi637(j,kpds,kgds,kret)
c   input argument list:
c     j        - size for indicated grid
c     kpds     -
c     kgds     -
c
c   output argument list:      (including work arrays)
c     j        - size for indicated grid modified for ecmwf-us 2
c     kret     - error return
c                (a mismatch was detected if kret is not zero)
c
c remarks:
c     kret     -
c          = 9 - gds indicates size mismatch with std grid
c
c   subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: fortran 77
c   machine:  hds
c
c$$$
      integer       kpds(*)
      integer       kgds(*)
      integer       j
      integer       i
c  ---------------------------------------
c  ---------------------------------------
c           if gds not indicated, return
c  ----------------------------------------
      kret=0
      if (iand(kpds(4),128).eq.0) return
c  ---------------------------------------
c            gds is indicated, proceed with testing
c  ---------------------------------------
      if (kgds(2).eq.65535) then
          return
      end if
      kret=1
      i     = kgds(2) * kgds(3)
c  ---------------------------------------
c            international set
c  ---------------------------------------
      if (kpds(3).ge.21.and.kpds(3).le.26) then
          if (i.ne.j) then
               return
          end if
      else if (kpds(3).ge.37.and.kpds(3).le.44) then
          if (i.ne.j) then
              return
          end if
      else if (kpds(3).eq.50) then
          if (i.ne.j) then
              return
          end if
      else if (kpds(3).ge.61.and.kpds(3).le.64) then
          if (i.ne.j) then
              return
          end if
c  ---------------------------------------
c            test ecmwf content
c  ---------------------------------------
      else if (kpds(1).eq.98) then
          kret  = 9
          if (kpds(3).ge.1.and.kpds(3).le.16) then
              if (i.ne.j) then
                if (kpds(3) .ne. 2) then 
                  return
                elseif (i .ne. 10512) then ! test for us grid 2
                  return
                end if
                j  = i   ! set to us grid 2, 2.5 global
              end if
          else
              kret  = 5
              return
          end if
c  ---------------------------------------
c           u.k. met office, bracknell
c  ---------------------------------------
      else if (kpds(1).eq.74) then
          kret  = 9
          if (kpds(3).ge.25.and.kpds(3).le.26) then
              if (i.ne.j) then
                  return
              end if
          else
              kret  = 5
              return
          end if
c  ---------------------------------------
c           canada
c  ---------------------------------------
      else if (kpds(1).eq.54) then
c         print *,' no current listing of canadian grids'
          return
c  ---------------------------------------
c           japan meteorological agency
c  ---------------------------------------
      else if (kpds(1).eq.34) then
c         print *,' no current listing of jma grids'
          return
c  ---------------------------------------
c           navy - fnoc
c  ---------------------------------------
      else if (kpds(1).eq.58) then
          if (kpds(3).ge.37.and.kpds(3).le.44) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.220.and.kpds(3).le.221) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.223) then
              if (i.ne.j) then
                  return
              end if
          else
              kret = 5
              return
          end if
c  ---------------------------------------
c                 u.s. grids
c  ---------------------------------------
      else if (kpds(1).eq.7) then
          kret  = 9
          if (kpds(3).ge.1.and.kpds(3).le.6) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.8) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.27.and.kpds(3).le.30) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.33.and.kpds(3).le.34) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.37.and.kpds(3).le.44) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.53) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.55.and.kpds(3).le.56) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.67.and.kpds(3).le.77) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.85.and.kpds(3).le.88) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.90.and.kpds(3).le.98) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.100.or.kpds(3).eq.101) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.103.and.kpds(3).le.107) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.110) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.126.or.kpds(3).eq.127) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.130) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.145.and.kpds(3).le.148) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.160.or.kpds(3).eq.161) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.163) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.170.and.kpds(3).le.175) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.190.or.kpds(3).eq.192) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.194.or.kpds(3).eq.196) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).eq.198) then
              if (i.ne.j) then
                  return
              end if
          else if (kpds(3).ge.201.and.kpds(3).le.254) then
              if (i.ne.j) then
                  return
              end if
          else
              kret  = 5
              return
          end if
      else
          kret  = 10
          return
      end if
c  ------------------------------------
c                    normal exit
c  ------------------------------------
      kret  = 0
      return
      end
