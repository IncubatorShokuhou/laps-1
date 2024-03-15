      subroutine w3fi76(pval,kexp,kmant,kbits)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:  w3fi76        convert to ibm370 floating point
c   prgmmr: rejones          org: nmc421      date:92-11-16
c
c abstract: converts floating point number from machine
c   representation to grib representation (ibm370 32 bit f.p.).
c
c program history log:
c   85-09-15  john hennessy  ecmwf
c   92-09-23  jones r. e.    change name, add doc block
c   93-10-27  jones,r. e.    change to agree with hennessy changes
c   95-10-31  iredell        removed saves and prints
c   98-03-10  b. vuong       remove the cdir$ integer=64 directive
c
c usage:    call w3fi76 (fval, kexp, kmant, nbits)
c   input argument list:
c     pval     - floating point number to be converted
c     kbits    - number of bits in computer word (32 or 64)
c
c   output argument list:
c     kexp     -  8 bit signed exponent
c     kmant    - 24 bit  mantissa  (fraction)
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm370 vs fortran 77, cray cft77 fortran
c   machine:  hds 9000, cray y-mp8/864< cray y-mp el2/256
c
c$$$
c
c********************************************************************
c*
c*    name      : confp3
c*
c*    function  : convert floating point number from machine
c*                representation to grib representation.
c*
c*    input     : pval  - floating point number to be converted.
c*    kbits     : kbits - number of bits in computer word
c*
c*    output    : kexp  - 8 bit signed exponent
c*                kmant - 24 bit mantissa
c*                pval  - unchanged.
c*
c*    john hennessy , ecmwf   18.06.91
c*
c********************************************************************
c
c
c     implicit none
c
      integer iexp
      integer isign
c
      integer kbits
      integer kexp
      integer kmant
c
      real pval
      real zeps
      real zref
c
c     test for floating point zero
c
      if (pval.eq.0.0) then
        kexp  = 0
        kmant = 0
        go to 900
      endif
c
c     set zeps to 1.0e-12 for 64 bit computers (cray)
c     set zeps to 1.0e-8  for 32 bit computers
c
      if (kbits.eq.32) then
        zeps = 1.0e-8
      else
        zeps = 1.0e-12
      endif
      zref = pval
c
c     sign of value
c
      isign = 0
      if (zref.lt.0.0) then
        isign =   128
        zref    = - zref
      endif
c
c     exponent
c
      iexp = int(alog(zref)*(1.0/alog(16.0))+64.0+1.0+zeps)
c
      if (iexp.lt.0  ) iexp = 0
      if (iexp.gt.127) iexp = 127
c
c     mantissa
c
c     closest number in grib format to original number
c     (equal to, greater than or less than original number).
c
      kmant = nint (zref/16.0**(iexp-70))
c
c     check that mantissa value does not exceed 24 bits
c     16777215 = 2**24 - 1
c
      if (kmant.gt.16777215) then
         iexp  = iexp + 1
c
c     closest number in grib format to original number
c     (equal to, greater than or less than original number).
c
         kmant = nint (zref/16.0**(iexp-70))
c
c        check mantissa value does not exceed 24 bits again
c
         if (kmant.gt.16777215) then
           print *,'bad mantissa value for pval = ',pval
         endif
      endif
c
c     add sign bit to exponent.
c
      kexp = iexp + isign
c
  900 continue
c
      return
      end
