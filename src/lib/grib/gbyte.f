      subroutine gbyte(in,iout,iskip,nbyte)
c     
c this program written by.....
c             dr. robert c. gammill, consultant
c             national center for atmospheric research
c             may 1972
c this is the fortran version of gbyte.
c     
c             changes for fortran 90
c             august 1990  russell e. jones
c             national weather service
c             gbyte run without changes on the following compilers
c             microsoft fortran 5.0 optimizing compiler
c             svs 32 386 fortran 77 version v2.8.1b
c             sun fortran 1.3, 1.4
c             dec vax fortran
c             silicongraphics 3.3, 3.4 fortran 77
c             ibm370 vs compiler
c             intergraph green hills fortran clipper 1.8.4b
c     
c usage:    call gbyte (pckd,unpk,inofst,nbit)
c   input argument list:
c    pckd      -  the fullword in memory from which unpacking is to
c                 begin; successive fullwords will be fetched as
c                 required.
c    inofst    -  a fullword integer specifying the inital offset
c                 in bits of the first byte, counted from the
c                 leftmost bit in pckd.
c    nbits     -  a fullword integer specifying the number of bits
c                 in each byte to be unpacked.  legal byte widths
c                 are in the range 1 - 32; bytes of width .lt. 32
c                 will be right justified in the low-order positions
c                 of the unpk fullwords, with high-order zero fill.
c     
c   output argument list:
c     unpk:    -  the fullword in memory into which the initial byte
c                 of unpacked data is to be stored.
c     
      integer    in(*)
      integer    iout
      integer    masks(32)
c     
      data  nbitsw/32/
c     
c      data  masks /z'00000001',z'00000003',z'00000007',z'0000000f',
c     &             z'0000001f',z'0000003f',z'0000007f',z'000000ff',
c     &             z'000001ff',z'000003ff',z'000007ff',z'00000fff',
c     &             z'00001fff',z'00003fff',z'00007fff',z'0000ffff',
c     &             z'0001ffff',z'0003ffff',z'0007ffff',z'000fffff',
c     &             z'001fffff',z'003fffff',z'007fffff',z'00ffffff',
c     &             z'01ffffff',z'03ffffff',z'07ffffff',z'0fffffff',
c     &             z'1fffffff',z'3fffffff',z'7fffffff',z'ffffffff'/
c     
c     mask table put in decimal so it will compile on an 32 bit
c     computer
c     
      data  masks / 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047,
     & 4095, 8191, 16383, 32767, 65535, 131071, 262143, 524287,
     & 1048575, 2097151, 4194303, 8388607, 16777215, 33554431,
     & 67108863, 134217727, 268435455, 536870911, 1073741823,
     & 2147483647, -1/
c     
c nbyte must be less than or equal to nbitsw
c     
      icon   = nbitsw - nbyte
      if (icon.lt.0) return
      mask   = masks(nbyte)
c     
c index tells how many words into the array 'in' the next byte appears.
c     
      index  = iskip / nbitsw
c     
c ii tells how many bits the byte is from the left side of the word.
c     
      ii     = mod(iskip,nbitsw)
c     
c mover specifies how far to the right a byte must be moved in order
c    to be right adjusted.
c     
      mover = icon - ii
c     
c right adjust the byte.
c     
      if (mover.gt.0) then
        iout = iand(ishft(in(index+1),-mover),mask)
c     
c the byte is split across a word break.
c     
      else if (mover.lt.0) then
        movel = - mover
        mover = nbitsw - movel
        iout  = iand(ior(ishft(in(index+1),movel),
     &            ishft(in(index+2),-mover)),mask)
c     
c the byte is already right adjusted.
c     
      else
        iout = iand(in(index+1),mask)
      endif
c     
        return
      end

