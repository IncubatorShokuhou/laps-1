      subroutine fi7502 (iwork,istart,npts,isame)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    fi7502      second order same value collection
c   prgmmr: cavanaugh        org: w/nmc42    date: 93-06-23
c
c abstract: collect sequential same values for processing
c   as second order value for grib messages.
c
c program history log:
c   93-06-23  cavanaugh
c   95-10-31  iredell     removed saves and prints
c
c usage:    call fi7502 (iwork,istart,npts,isame)
c   input argument list:
c     iwork    - array containing source data
c     istart   - starting location for this test
c     npts     - number of points in iwork
c
c   output argument list:      (including work arrays)
c     isame    - number of sequential points having the same value
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916/256, y-mp8/64, y-mp el92/256
c
c$$$
      integer        iwork(*)
      integer        istart
      integer        isame
      integer        k
      integer        npts
c  -------------------------------------------------------------
      isame  = 0
      do 100 k = istart, npts
          if (iwork(k).ne.iwork(istart)) then
              return
          end if
          isame  = isame + 1
  100 continue
      return
      end
