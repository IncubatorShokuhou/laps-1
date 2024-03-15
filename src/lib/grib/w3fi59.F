      subroutine w3fi59(field,npts,nbits,nwork,npfld,iscale,len,rmin)
cdir$ integer=64
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    w3fi59      form and pack positive, scaled differences
c   prgmmr:  allard, r.      org: nmc41      date:  84-08-01
c
c abstract:  converts an array of single precision real numbers into
c   an array of positive scaled differences (number(s) - minimum value),
c   in integer format and packs the argument-specified number of
c   significant bits from each difference.
c
c program history log:
c   84-08-01  allard      original author
c   90-05-17  r.e.jones   convert to cray cft77 fortran
c   90-05-18  r.e.jones   change name pakmag to w3lib name w3fi59
c
c usage:    call w3fi59(field,npts,nbits,nwork,npfld,iscale,len,rmin)
c   input argument list:
c     field - array of floating point data for processing  (real)
c     npts  - number of data values to process in field (and nwork)
c             where, npts > 0
c     nbits - number of significant bits of processed data to be packed
c             where, 0 < nbits < 33
c
c   output argument list:
c     nwork - array for integer conversion  (integer)
c             if packing performed (see note below), the array will
c             contain the pre-packed, right adjusted, scaled, integer
c             differences upon return to the user.
c             (the user may equivalence field and nwork.  same size.)
c     npfld - array for packed data  (integer)
c             (dimension must be at least (nbits * npts) / 64 + 1  )
c     iscale- power of 2 for restoring data, such that
c             datum = (difference * 2**iscale) + rmin
c     len   - number of packed bytes in npfld (set to 0 if no packing)
c             where, len = (nbits * npts + 7) / 8 without remainder
c     rmin  - minimum value (reference value subtracted from input data)
c             this is a cray floating point number, it will have to be
c             converted to an ibm370 32 bit floating point number at
c             some point in your program if you are packing grib data.
c
c remarks:  len = 0 and no packing performed if
c
c        (1)  rmax = rmin  (a constant field)
c        (2)  nbits value out of range  (see input argument)
c        (3)  npts value less than 1  (see input argument)
c
c attributes:
c   language: cray cft77 fortran
c   machine:  cray y-mp8/832
c
c$$$
c
      real    field(*)
c
      integer npfld(*)
      integer nwork(*)
c
      equivalence (bigdif,ix)
c
      data kzero / 0 /
c
c / / / / / /
c
      len    = 0
      iscale = 0
      if (nbits.le.0.or.nbits.gt.32) go to 3000
      if (npts.le.0) go to 3000
c
c find the max-min values in field.
c
      rmax = field(1)
      rmin = rmax
      do 1000 k = 2,npts
        rmax = amax1(rmax,field(k))
        rmin = amin1(rmin,field(k))
 1000 continue
c
c if a constant field, return with no packing performed and 'len' = 0.
c
      if (rmax.eq.rmin) go to 3000
c
c determine largest difference in field (bigdif).
c
      bigdif = rmax - rmin
c
#ifdef cray
c--------------------------------------------------------------------------
c                       64-bit machines
c right-adjust the exponent of bigdif and save.  (32767 = 7fff hex).
c
      n = and(shift(ix,16),32767)
c
c remove cray bias of 16384 (4000 hex) from cray exponent.
c (base to its nth power gives the denominator that
c expresses the number as a proper fraction.)
c
      n = n - 16384
c
c after adjustment n is the least power of 2 such that (2**n) > bigdif.
c
#endif
#ifndef cray
c--------------------------------------------------------------------------
c                       32-bit machines
c right-adjust the exponent of bigdif and save.  (255 = ff hex).
c
      n = iand(ishft(ix,-23),255)
c
c remove bias of 127
c (base to its nth power gives the denominator that
c expresses the number as a proper fraction.)
c
      n = n - 127 + 1
c
c after adjustment n is the least power of 2 such that (2**n) > bigdif.
c
#endif
c--------------------------------------------------------------------------
c
c form differences, rescale, and convert to integer format.
c
      twon = 2.0 ** (nbits - n)
      do 2000 k = 1,npts
        nwork(k) = nint( (field(k) - rmin) * twon )
 2000 continue
c
c pack the magnitudes (rightmost nbits of each word).
c
c note: nbits = number of significant bits of the magnitude saved.
c       n = number of bits to left of binary point in scaled integer.
c       nbits-n = number of fraction bits to right of binary point in
c            scaled integer.
c
      koff  = 0
      iskip = 0
c
c     use ncar array bit packer sbytes  (gbytes package)
c
      call sbytes(npfld,nwork,koff,nbits,iskip,npts)
c
c add 7 zero-bits at end of packed data to insure byte boundary.
c     use ncar word bit packer sbyte
c
      noff = nbits * npts
      call sbyte(npfld,kzero,noff,7)
c
c determine byte length (len) of packed field (npfld).
c
      len = (noff + 7) / 8
c
c set power of 2 argument needed to restore data.
c
      iscale = n - nbits
c
 3000 continue
      return
c
      end
