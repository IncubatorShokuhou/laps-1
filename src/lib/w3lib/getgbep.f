c-----------------------------------------------------------------------
      subroutine getgbep(lugb,lugi,jg,j,jpds,jgds,jens,
     &                   kg,k,kpds,kgds,kens,g,iret)
c$$$  subprogram documentation block
c
c subprogram: getgbep        finds a grib message
c   prgmmr: iredell          org: w/nmc23     date: 94-04-01
c
c abstract: find a grib message.
c   read a grib index file (or optionally the grib file itself)
c   to get the index buffer (i.e. table of contents) for the grib file.
c   (the index buffer is saved for use by future prospective calls.)
c   find in the index buffer a reference to the grib message requested.
c   the grib message request specifies the number of messages to skip
c   and the unpacked pds and gds parameters.  (a requested parameter
c   of -1 means to allow any value of this parameter to be found.)
c   if the requested grib message is found, then it is read from the
c   grib file.  its message number is returned along with the unpacked
c   pds and gds parameters and the packed grib message.  if the grib
c   message is not found, then the return code will be nonzero.
c
c program history log:
c   94-04-01  iredell
c   95-10-31  iredell     modularized portions of code into subprograms
c                         and allowed for unspecified index file
c
c usage:    call getgbep(lugb,lugi,jg,j,jpds,jgds,jens,
c    &                   kg,k,kpds,kgds,kens,g,iret)
c   input arguments:
c     lugb         integer unit of the unblocked grib data file
c     lugi         integer unit of the unblocked grib index file
c                  (=0 to get index buffer from the grib file)
c     jg           integer maximum number of bytes in the grib message
c     j            integer number of messages to skip
c                  (=0 to search from beginning)
c                  (<0 to read index buffer and skip -1-j messages)
c     jpds         integer (200) pds parameters for which to search
c                  (=-1 for wildcard)
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
c     jgds         integer (200) gds parameters for which to search
c                  (only searched if jpds(3)=255)
c                  (=-1 for wildcard)
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
c     jens         integer (200) ensemble pds parms for which to search
c                  (only searched if jpds(23)=2)
c                  (=-1 for wildcard)
c          (1)   - application identifier
c          (2)   - ensemble type
c          (3)   - ensemble identifier
c          (4)   - product identifier
c          (5)   - smoothing flag
c   output arguments:
c     kg           integer number of bytes in the grib message
c     k            integer message number unpacked
c                  (can be same as j in calling program
c                  in order to facilitate multiple searches)
c     kpds         integer (200) unpacked pds parameters
c     kgds         integer (200) unpacked gds parameters
c     kens         integer (200) unpacked ensemble pds parms
c     g            character*1 (kg) grib message
c     iret         integer return code
c                    0      all ok
c                    96     error reading index file
c                    97     error reading grib file
c                    98     number of bytes greater than jg
c                    99     request not found
c
c subprograms called:
c   getgbemp       find grib message
c
c remarks: in order to unpack grib from a multiprocessing environment
c   where each processor is attempting to read from its own pair of
c   logical units, one must directly call subprogram getgbemp as below,
c   allocating a private copy of cbuf, nlen and nnum to each processor.
c   do not engage the same logical unit from more than one processor.
c
c attributes:
c   language: fortran 77
c   machine:  cray, workstations
c
c$$$
      integer jpds(200),jgds(200),jens(200)
      integer kpds(200),kgds(200),kens(200)
      character g(jg)
      parameter(mbuf=256*1024)
      character cbuf(mbuf)
      save cbuf,nlen,nnum,mnum
      data lux/0/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  determine whether index buffer needs to be initialized
      if(lugi.gt.0.and.(j.lt.0.or.lugi.ne.lux)) then
        lux=lugi
        jj=min(j,-1-j)
      elseif(lugi.le.0.and.(j.lt.0.or.lugb.ne.lux)) then
        lux=lugb
        jj=min(j,-1-j)
      else
        jj=j
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  find and unpack grib message
      call getgbemp(lugb,lugi,jg,jj,jpds,jgds,jens,
     &              mbuf,cbuf,nlen,nnum,mnum,
     &              kg,k,kpds,kgds,kens,g,iret)
      if(iret.eq.96) lux=0
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
