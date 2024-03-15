      subroutine sbytes(iout,in,iskip,nbyte,nskip,n)
c this program written by.....
c             dr. robert c. gammill, consultant
c             national center for atmospheric research
c             july 1972
c this is the fortran versions of sbytes.
c     
c             fortran 90
c             august 1990  russell e. jones
c             national weather service
c     
c usage:    call sbytes (pckd,unpk,inofst,nbit, nskip,iter)
c     
c   input argument list:
c     unpk     -  nbits of the right side of each word of array
c                 unpk is moved to array pckd. inofst bits are
c                 skipped over before the 1st data is moved, nbits
c                 are stored, nskip bits are skipped over, the next
c                 nbits are moved,  bit are skipped over, etc. until
c                 iter groups of bits are packed.
c    inofst    -  a fullword integer specifying the inital offset
c                 in bits of the first byte, counted from the
c                 leftmost bit in pckd.
c    nbits     -  a fullword integer specifying the number of bits
c                 in each byte to be packed.  legal byte widths
c                 are in the range 1 - 32.
c    nskip     -  a fullword integer specifying the number of bits
c                 to skip between successive bytes.  all non-negative
c                 skip counts are legal.
c    iter      -  a fullword integer specifying the total number of
c                 bytes to be packed, as controlled by inofst,
c                 nbit and nskip above.   all non-negative iteration
c                 counts are legal.
c     
c   output argument list:
c    pckd      -  the fullword in memory to which packing is to
c                 begin starting at bit inofst. the inostat bits
c                 are not altered. nskip bits are not altered.
c     
      integer    in(*)
      integer    iout(*)
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
c     masks table put in decimal so it will compile on any 32 bit
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
      icon = nbitsw - nbyte
      if (icon.lt.0) return
      mask   = masks(nbyte)
c     
c index tells how many words into iout the next byte is to be stored.
c     
      index  = iskip / nbitsw
c     
c ii tells how many bits in from the left side of the word to store it.
c     
      ii     = mod(iskip,nbitsw)
c     
c istep is the distance in bits from one byte position to the next.
c     
      istep  = nbyte + nskip
c     
c iwords tells how many words to skip from one byte to the next.
c     
      iwords = istep / nbitsw
c     
c ibits tells how many bits to skip after skipping iwords.
c     
      ibits  = mod(istep,nbitsw)
c     
      do 10 i = 1,n
        j     = iand(mask,in(i))
        movel = icon - ii
c     
c byte is to be stored in middle of word.  shift left.
c     
        if (movel.gt.0) then
          msk           = ishft(mask,movel)
          iout(index+1) = ior(iand(not(msk),iout(index+1)),
     &    ishft(j,movel))
c     
c the byte is to be split across a word break.
c     
        else if (movel.lt.0) then
          msk           = masks(nbyte+movel)
          iout(index+1) = ior(iand(not(msk),iout(index+1)),
     &    ishft(j,movel))
          itemp         = iand(masks(nbitsw+movel),iout(index+2))
          iout(index+2) = ior(itemp,ishft(j,nbitsw+movel))
c     
c byte is to be stored right-adjusted.
c     
        else
          iout(index+1) = ior(iand(not(mask),iout(index+1)),j)
        endif
c     
        ii    = ii + ibits
        index = index + iwords
        if (ii.ge.nbitsw) then
          ii    = ii - nbitsw
          index = index + 1
        endif
c     
10    continue
c     
      return
      end

