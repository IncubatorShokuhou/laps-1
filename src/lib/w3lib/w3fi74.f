      subroutine w3fi74 (igds,icomp,gds,lengds,npts,igerr)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    w3fi74      construct grid definition section (gds)
c   prgmmr: farley           org: w/nmc42    date: 93-08-24
c
c abstract: this subroutine constructs a grib grid definition
c   section.
c
c program history log:
c   92-07-07  m. farley   original author
c   92-10-16  r.e.jones   add code to lat/lon section to do
c                         gaussian grids.
c   93-03-29  r.e.jones   add save statement
c   93-08-24  r.e.jones   changes for grib grids 37-44
c   93-09-29  r.e.jones   changes for gaussian grid for document
c                         change in w3fi71.
c   94-02-15  r.e.jones   changes for eta model grids 90-93
c   95-04-20  r.e.jones   change 200 and 201 to 201 and 202
c   95-10-31  iredell     removed saves and prints
c   98-08-20  baldwin     add type 203
c
c
c usage:    call w3fi74 (igds, icomp, gds, lengds, npts, igerr)
c   input argument list:
c     igds        - integer array supplied by w3fi71
c     icomp       - table 7- resolution & component flag (bit 5)
c                   for gds(17) wind components
c
c   output argument list:
c     gds       - completed grib grid definition section
c     lengds    - length of gds
c     npts      - number of points in grid
c     igerr     - 1, grid representation type not valid
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: cray cft77 fortran 77, ibm370 vs fortran
c   machine:  cray c916-128, cray y-mp8/864, cray y-mp el2/256, hds
c
c$$$
c
      integer       igds  (*)
c
      character*1   gds   (*)
c
      isum  = 0
      igerr = 0
c
c       print *,' '
c       print *,'(w3fi74-igds = )'
c       print *,(igds(i),i=1,18)
c       print *,' '
c
c     compute length of gds in octets (octets 1-3)
c       length =  32 for lat/lon, gnomic, gausian lat/lon,
c                    polar stereographic, spherical harmonics
c       length =  42 for mercator, lambert, tangent cone
c       length = 178 for mercator, lambert, tangent cone
c
      if (igds(3) .eq. 0  .or.  igds(3) .eq. 2  .or.
     &    igds(3) .eq. 4  .or.  igds(3) .eq. 5  .or.
     &    igds(3) .eq. 50 .or.  igds(3) .eq. 201.or.
     &    igds(3) .eq. 202.or.  igds(3) .eq. 203) then
          lengds = 32
c
c       correction for grids 37-44
c
        if (igds(3).eq.0.and.igds(1).eq.0.and.igds(2).ne.
     &  255) then
          lengds = igds(5) * 2 + 32
        endif
      else if (igds(3) .eq. 1  .or.  igds(3) .eq. 3  .or.
     &         igds(3) .eq. 13) then
        lengds = 42
      else
c       print *,' w3fi74 error, grid representation type not valid'
        igerr = 1
        return
      endif
c
c     put length of gds section in bytes 1,2,3
c
      gds(1) = char(mod(lengds/65536,256))
      gds(2) = char(mod(lengds/  256,256))
      gds(3) = char(mod(lengds      ,256))
c
c     octet 4 = nv, number of vertical coordinate parameters
c     octet 5 = pv, pl or 255
c     octet 6 = data representation type (table 6)
c
      gds(4) = char(igds(1))
      gds(5) = char(igds(2))
      gds(6) = char(igds(3))
c
c     fill octet the rest of the gds based on data representation
c     type (table 6)
c
c$$
c     process lat/lon grid types or gaussian grid or arakawa
c     staggered, semi-staggered, or filled e-grids
c
      if (igds(3).eq.0.or.igds(3).eq.4.or.
     &    igds(3).eq.201.or.igds(3).eq.202.or.
     &    igds(3).eq.203) then
        gds( 7) = char(mod(igds(4)/256,256))
        gds( 8) = char(mod(igds(4)    ,256))
        gds( 9) = char(mod(igds(5)/256,256))
        gds(10) = char(mod(igds(5)    ,256))
        lato    = igds(6)
        if (lato .lt. 0) then
          lato = -lato
          lato = ior(lato,8388608)
        endif
        gds(11) = char(mod(lato/65536,256))
        gds(12) = char(mod(lato/  256,256))
        gds(13) = char(mod(lato      ,256))
        lono    = igds(7)
        if (lono .lt. 0) then
          lono = -lono
          lono = ior(lono,8388608)
        endif
        gds(14) = char(mod(lono/65536,256))
        gds(15) = char(mod(lono/  256,256))
        gds(16) = char(mod(lono      ,256))
        latext  = igds(9)
        if (latext .lt. 0) then
          latext = -latext
          latext = ior(latext,8388608)
        endif
        gds(18) = char(mod(latext/65536,256))
        gds(19) = char(mod(latext/  256,256))
        gds(20) = char(mod(latext      ,256))
        lonext  = igds(10)
        if (lonext .lt. 0) then
          lonext = -lonext
          lonext = ior(lonext,8388608)
        endif
        gds(21) = char(mod(lonext/65536,256))
        gds(22) = char(mod(lonext/  256,256))
        gds(23) = char(mod(lonext      ,256))
        ires    = iand(igds(8),128)
        if (igds(3).eq.201.or.igds(3).eq.202.or.igds(3).eq.203) then
          gds(24) = char(mod(igds(11)/256,256))
          gds(25) = char(mod(igds(11)    ,256))
        else if (ires.eq.0) then
          gds(24) = char(255)
          gds(25) = char(255)
        else
          gds(24) = char(mod(igds(12)/256,256))
          gds(25) = char(mod(igds(12)    ,256))
        end if
        if (igds(3).eq.4) then
          gds(26) = char(mod(igds(11)/256,256))
          gds(27) = char(mod(igds(11)    ,256))
        else if (igds(3).eq.201.or.igds(3).eq.202.or.
     &           igds(3).eq.203) then
          gds(26) = char(mod(igds(12)/256,256))
          gds(27) = char(mod(igds(12)    ,256))
        else if (ires.eq.0) then
          gds(26) = char(255)
          gds(27) = char(255)
        else
          gds(26) = char(mod(igds(11)/256,256))
          gds(27) = char(mod(igds(11)    ,256))
        end if
        gds(28) = char(igds(13))
        gds(29) = char(0)
        gds(30) = char(0)
        gds(31) = char(0)
        gds(32) = char(0)
        if (lengds.gt.32) then
          isum = 0
          i    = 19
          do 10 j = 33,lengds,2
            isum     = isum + igds(i)
            gds(j)   = char(mod(igds(i)/256,256))
            gds(j+1) = char(mod(igds(i)    ,256))
            i        = i + 1
 10       continue
        end if
c
c$$     process mercator grid types
c
      else if (igds(3) .eq. 1) then
        gds( 7) = char(mod(igds(4)/256,256))
        gds( 8) = char(mod(igds(4)    ,256))
        gds( 9) = char(mod(igds(5)/256,256))
        gds(10) = char(mod(igds(5)    ,256))
        lato = igds(6)
        if (lato .lt. 0) then
          lato = -lato
          lato = ior(lato,8388608)
        endif
        gds(11) = char(mod(lato/65536,256))
        gds(12) = char(mod(lato/  256,256))
        gds(13) = char(mod(lato      ,256))
        lono = igds(7)
        if (lono .lt. 0) then
          lono = -lono
          lono = ior(lono,8388608)
        endif
        gds(14) = char(mod(lono/65536,256))
        gds(15) = char(mod(lono/  256,256))
        gds(16) = char(mod(lono      ,256))
        latext = igds(9)
        if (latext .lt. 0) then
          latext = -latext
          latext = ior(latext,8388608)
        endif
        gds(18) = char(mod(latext/65536,256))
        gds(19) = char(mod(latext/  256,256))
        gds(20) = char(mod(latext      ,256))
        lonext  = igds(10)
        if (lonext .lt. 0) then
          lonext = -lonext
          lonext = ior(lonext,8388608)
        endif
        gds(21) = char(mod(lonext/65536,256))
        gds(22) = char(mod(lonext/  256,256))
        gds(23) = char(mod(lonext      ,256))
        gds(24) = char(mod(igds(13)/65536,256))
        gds(25) = char(mod(igds(13)/  256,256))
        gds(26) = char(mod(igds(13)      ,256))
        gds(27) = char(0)
        gds(28) = char(igds(14))
        gds(29) = char(mod(igds(12)/65536,256))
        gds(30) = char(mod(igds(12)/  256,256))
        gds(31) = char(mod(igds(12)      ,256))
        gds(32) = char(mod(igds(11)/65536,256))
        gds(33) = char(mod(igds(11)/  256,256))
        gds(34) = char(mod(igds(11)      ,256))
        gds(35) = char(0)
        gds(36) = char(0)
        gds(37) = char(0)
        gds(38) = char(0)
        gds(39) = char(0)
        gds(40) = char(0)
        gds(41) = char(0)
        gds(42) = char(0)
c$$     process lambert conformal grid types
      else if (igds(3) .eq. 3) then
        gds( 7) = char(mod(igds(4)/256,256))
        gds( 8) = char(mod(igds(4)    ,256))
        gds( 9) = char(mod(igds(5)/256,256))
        gds(10) = char(mod(igds(5)    ,256))
        lato = igds(6)
        if (lato .lt. 0) then
          lato = -lato
          lato = ior(lato,8388608)
        endif
        gds(11) = char(mod(lato/65536,256))
        gds(12) = char(mod(lato/  256,256))
        gds(13) = char(mod(lato      ,256))
        lono = igds(7)
        if (lono .lt. 0) then
          lono = -lono
          lono = ior(lono,8388608)
        endif
        gds(14) = char(mod(lono/65536,256))
        gds(15) = char(mod(lono/  256,256))
        gds(16) = char(mod(lono      ,256))
        lonm = igds(9)
        if (lonm .lt. 0) then
          lonm = -lonm
          lonm = ior(lonm,8388608)
        endif
        gds(18) = char(mod(lonm/65536,256))
        gds(19) = char(mod(lonm/  256,256))
        gds(20) = char(mod(lonm      ,256))
        gds(21) = char(mod(igds(10)/65536,256))
        gds(22) = char(mod(igds(10)/  256,256))
        gds(23) = char(mod(igds(10)      ,256))
        gds(24) = char(mod(igds(11)/65536,256))
        gds(25) = char(mod(igds(11)/  256,256))
        gds(26) = char(mod(igds(11)      ,256))
        gds(27) = char(igds(12))
        gds(28) = char(igds(13))
        gds(29) = char(mod(igds(15)/65536,256))
        gds(30) = char(mod(igds(15)/  256,256))
        gds(31) = char(mod(igds(15)      ,256))
        gds(32) = char(mod(igds(16)/65536,256))
        gds(33) = char(mod(igds(16)/  256,256))
        gds(34) = char(mod(igds(16)      ,256))
        gds(35) = char(mod(igds(17)/65536,256))
        gds(36) = char(mod(igds(17)/  256,256))
        gds(37) = char(mod(igds(17)      ,256))
        gds(38) = char(mod(igds(18)/65536,256))
        gds(39) = char(mod(igds(18)/  256,256))
        gds(40) = char(mod(igds(18)      ,256))
        gds(41) = char(0)
        gds(42) = char(0)
c$$     process polar stereographic grid types
      else if (igds(3) .eq. 5) then
        gds( 7) = char(mod(igds(4)/256,256))
        gds( 8) = char(mod(igds(4)    ,256))
        gds( 9) = char(mod(igds(5)/256,256))
        gds(10) = char(mod(igds(5)    ,256))
        lato = igds(6)
        if (lato .lt. 0) then
          lato = -lato
          lato = ior(lato,8388608)
        endif
        gds(11) = char(mod(lato/65536,256))
        gds(12) = char(mod(lato/  256,256))
        gds(13) = char(mod(lato      ,256))
        lono = igds(7)
        if (lono .lt. 0) then
          lono = -lono
          lono = ior(lono,8388608)
        endif
        gds(14) = char(mod(lono/65536,256))
        gds(15) = char(mod(lono/  256,256))
        gds(16) = char(mod(lono      ,256))
        lonm = igds(9)
        if (lonm .lt. 0) then
          lonm = -lonm
          lonm = ior(lonm,8388608)
        endif
        gds(18) = char(mod(lonm/65536,256))
        gds(19) = char(mod(lonm/   256,256))
        gds(20) = char(mod(lonm       ,256))
        gds(21) = char(mod(igds(10)/65536,256))
        gds(22) = char(mod(igds(10)/  256,256))
        gds(23) = char(mod(igds(10)      ,256))
        gds(24) = char(mod(igds(11)/65536,256))
        gds(25) = char(mod(igds(11)/  256,256))
        gds(26) = char(mod(igds(11)      ,256))
        gds(27) = char(igds(12))
        gds(28) = char(igds(13))
        gds(29) = char(0)
        gds(30) = char(0)
        gds(31) = char(0)
        gds(32) = char(0)
      endif
c       print 10,(gds(ig),ig=1,32)
c10     format ('  gds= ',32(1x,z2.2))
c
c     compute number of points in grid by multiplying
c       igds(4) and igds(5) ... needed for packer
c
      if (igds(3).eq.0.and.igds(1).eq.0.and.igds(2).ne.
     & 255) then
        npts = isum
      else
        npts = igds(4) * igds(5)
      endif
c
c     'ior' icomp-bit 5 resolution & component flag for winds
c       with igds(8) info (rest of resolution & component flag data)
c
      icomp   = ishft(icomp,3)
      gds(17) = char(ior(igds(8),icomp))
c
      return
      end
