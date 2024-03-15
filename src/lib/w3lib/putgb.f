c-----------------------------------------------------------------------
      subroutine putgb(lugb,kf,kpds,kgds,lb,f,iret)
c$$$  subprogram documentation block
c
c subprogram: putgb          packs and writes a grib message
c   prgmmr: iredell          org: w/nmc23     date: 94-04-01
c
c abstract: pack and write a grib message.
c   this subprogram is nearly the inverse of getgb.
c
c program history log:
c   94-04-01  iredell
c   95-10-31  iredell     removed saves and prints
c
c usage:    call putgb(lugb,kf,kpds,kgds,lb,f,iret)
c   input arguments:
c     lugb         integer unit of the unblocked grib data file
c     kf           integer number of data points
c     kpds         integer (200) pds parameters
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
c     kgds         integer (200) gds parameters
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
c     lb           logical*1 (kf) bitmap if present
c     f            real (kf) data
c   output arguments:
c     iret         integer return code
c                    0      all ok
c                    other  w3fi72 grib packer return code
c
c subprograms called:
c   r63w72         map w3fi63 parameters onto w3fi72 parameters
c   getbit         get number of bits and round data
c   w3fi72         pack grib
c   wryte          write data
c
c remarks: subprogram can be called from a multiprocessing environment.
c   do not engage the same logical unit from more than one processor.
c
c attributes:
c   language: fortran 77
c   machine:  cray, workstations
c
c$$$
      integer kpds(200),kgds(200)
      logical*1 lb(kf)
      real f(kf)
      parameter(maxbit=16)
      integer ibm(kf),ipds(200),igds(200),ibds(200)
      real fr(kf)
      character pds(400),grib(1000+kf*(maxbit+1)/8)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  get w3fi72 parameters
      call r63w72(kpds,kgds,ipds,igds)
      ibds=0
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  count valid data
      kbm=kf
      if(ipds(7).ne.0) then
        kbm=0
        do i=1,kf
          if(lb(i)) then
            ibm(i)=1
            kbm=kbm+1
          else
            ibm(i)=0
          endif
        enddo
        if(kbm.eq.kf) ipds(7)=0
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  get number of bits and round data
      if(kbm.eq.0) then
        do i=1,kf
          fr(i)=0.
        enddo
        nbit=0
      else
        call getbit(ipds(7),0,ipds(25),kf,ibm,f,fr,fmin,fmax,nbit)
        nbit=min(nbit,maxbit)
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  pack and write grib data
      call w3fi72(0,fr,0,nbit,0,ipds,pds,
     &            1,255,igds,0,0,ibm,kf,ibds,
     &            kfo,grib,lgrib,iret)
      if(iret.eq.0) call wryte(lugb,lgrib,grib)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
