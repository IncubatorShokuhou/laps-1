      subroutine w3fi83 (data,npts,fval1,fdiff1,iscal2,
     *                                isc10,kpds,kgds)
c$$$  subprogram documentation  block
c                .      .    .                                       .
c subprogram:  w3fi83        restore delta packed data to original
c   prgmmr: cavanaugh        org: nmc421      date:93-08-18
c
c abstract: restore delta packed data to original values
c           restore from boustrephedonic alignment
c
c program history log:
c   93-07-14  cavanaugh
c   93-07-22  stackpole      additions to fix scaling
c   94-01-27  cavanaugh   added reversal of even numbered rows
c                         (boustrophedonic processing) to restore
c                         data to original sequence.
c   94-03-02  cavanaugh   corrected reversal of even numbered rows
c   95-10-31  iredell     removed saves and prints
c
c usage:    call w3fi83(data,npts,fval1,fdiff1,iscal2,
c    *                                isc10,kpds,kgds)
c   input argument list:
c     data     - second order differences
c     npts     - number of points in array
c     fval1    - original first entry in array
c     fdiff1   - original first first-difference
c     iscal2   - power-of-two exponent for unscaling
c     isc10    - power-of-ten exponent for unscaling
c     kpds     - array of information for pds
c     kgds     - array of information for gds
c
c   output argument list:
c     data     - expanded original data values
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916-128, cray y-mp8/864, cray y-mp el2/256
c
c$$$
c
      real          fval1,fdiff1
      real          data(*),boust(200)
      integer       npts,nrow,ncol,kpds(*),kgds(*),isc10
c  ---------------------------------------
c
c     remove decimal un-scaling introduced during unpacking
c
      dscal = 10.0 ** isc10
      if (dscal.eq.0.0) then
          do 50 i=1,npts
              data(i) = 1.0
   50     continue
      else if (dscal.eq.1.0) then
      else
          do 51 i=1,npts
              data(i) = data(i) * dscal
   51     continue
      end if
c
      data(1)  = fval1
      data(2)  = fdiff1
      do 200 j = 3,2,-1
          do 100 k = j, npts
              data(k)  = data(k) + data(k-1)
  100     continue
  200 continue
c
c     now remove the binary scaling from the reconstructed field
c     and the decimal scaling too
c
      if (dscal.eq.0) then
          scale  = 0.0
      else
          scale =(2.0**iscal2)/dscal
      end if
      do 300 i=1,npts
        data(i) = data(i) * scale
  300 continue
c  ==========================================================
      if (iand(kpds(4),128).ne.0) then
          nrow  = kgds(3)
          ncol  = kgds(2)
c
c      data laid out boustrophedonic style
c
c
c         print*, '  reverse boustrophedon'
          do 210 i = 2, nrow, 2
c
c          reverse the even numbered rows
c
              do 201 j = 1, ncol
                  npos  = i * ncol - j + 1
                  boust(j) = data(npos)
  201         continue
              do 202 j = 1, ncol
                  npos  = ncol * (i-1) + j
                  data(npos)  = boust(j)
  202         continue
  210     continue
c
c
      end if
c  =================================================================
      return
      end
