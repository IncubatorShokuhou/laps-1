      subroutine w3fi58(ifield,npts,nwork,npfld,nbits,len,kmin)
cdir$ integer=64
c$$$  subprogram documentation block  ***
c                .      .    .                                       .
c subprogram:  w3fi58   - pack positive differences in least bits
c   prgmmr:  allard, r.       org:  nmc411        date:  july 1987
c
c abstract:  converts an array of integer numbers into an array of
c   positive differences (number(s) - minimum value) and packs the
c   magnitude of each difference right-adjusted into the least
c   number of bits that holds the largest difference.
c
c program history log:
c   87-09-02  allard
c   88-10-02  r.e.jones   converted to cdc cyber 205 ftn200 fortran
c   90-05-17  r.e.jones   converted to cray cft77 fortran
c   90-05-18  r.e.jones   change name vbimpk to w3lib name w3fi58
c
c usage:  call w3fi58(ifield,npts,nwork,npfld,nbits,len,kmin)
c
c   input:
c
c     ifield - array of integer data for processing
c     npts   - number of data values to process in ifield (and nwork)
c              where, npts > 0
c
c   output:
c
c     nwork  - work array with integer difference
c     npfld  - array for packed data
c              (user is responsible for an adequate dimension.)
c     nbits  - number of bits used to pack data where, 0 < nbits < 32
c              (the maximum difference without overflow is 2**31 -1)
c     len    - number of packed bytes in npfld (set to 0 if no packing)
c              where, len = (nbits * npts + 7) / 8 without remainder
c     kmin   - minimum value (subtracted from each datum). if this
c              packed data is being used for grib data, the
c              programer will have to convert the kmin value to an
c              ibm370 32 bit floating point number.
c
c   subprograms called:
c
c     w3lib:  sbytes, sbyte
c
c   exit states:  none
c
c     note:  len = 0, nbits = 0, and no packing performed if
c
c     (1) kmax = kmin  (a constant field)
c     (2) npts < 1  (see input argument)
c
c attributes:
c   language: cray cft77 fortran
c   machine:  cray y-mp8/832
c
c$$$
c
      integer  ifield(*)
      integer  npfld(*)
      integer  nwork(*)
      integer n2work(2)
c
      equivalence (bigdif,ix)
c
      data  kzero / 0 /
c
c / / / / / /
c
      len   = 0
      nbits = 0
      if (npts.le.0) go to 3000
c
c find the max-min values in integer field (ifield).
c
      kmax = ifield(1)
      kmin = kmax
      do 1000 i = 2,npts
        kmax = max(kmax,ifield(i))
        kmin = min(kmin,ifield(i))
 1000 continue
c
c if a constant field, return with no packing and 'len' and 'nbits' set
c to zero.
c
      if (kmax.eq.kmin) go to 3000
c
c determine largest difference in ifield and float (bigdif).
c
      bigdif = kmax - kmin
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
c form differences in nwork array.
c
      do 2000 k = 1,npts
        nwork(k) = ifield(k) - kmin
 2000 continue
c
c pack each magnitude in nbits (nbits = the least power of 2 or 'n')
c
      nbits  = n
      koff   = 0
      iskip  = 0
c
c     use ncar array bit packer sbytes  (gbytes package)
c
      jpts = npts
      call sbytes(npfld,nwork,koff,nbits,iskip,jpts)
c
c add 7 zero-bits at end of packed data to insure a byte boundary.
c     use ncar word bit packer sbyte
c
      noff = nbits * jpts
      call sbyte(npfld,kzero,noff,7)
c
c determine byte length (len) of packed field (npfld).
c
      len = (noff + 7) / 8
c
 3000 continue
      return
c
      end
