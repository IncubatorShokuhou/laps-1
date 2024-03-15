      subroutine w3fi82 (ifld,fval1,fdiff1,npts)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:  w3fi82        convert to second diff array
c   prgmmr: cavanaugh        org: nmc421      date:93-08-18
c
c abstract: accept an input array, convert to array of second
c   differences.  return the original first value and the first
c   first-difference as separate values.
c
c program history log:
c   93-07-14  cavanaugh
c   93-08-18  r.e.jones   recompile for silicongraphics
c   95-10-31  iredell     removed saves and prints
c
c usage:    call w3fi82 (ifld,fval1,fdiff1,npts)
c   input argument list:
c     ifld     - integer input array
c     npts     - number of points in array
c
c   output argument list:
c     ifld     - second differenced field
c     fval1    - floating point original first value
c     fdiff1   -     "      "   first first-difference
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: silicongraphics 3.5 fortran 77
c   machine:  silicongraphics model 25, 35, indigo
c
c$$$
c
      real        fval1,fdiff1
c
      integer     ifld(*),npts
c
c  ---------------------------------------------
          do 4000 i = npts, 2, -1
              ifld(i)  = ifld(i) - ifld(i-1)
 4000     continue
          do 5000 i = npts, 3, -1
              ifld(i)  = ifld(i) - ifld(i-1)
 5000     continue
c         print *,'ifld(1) =',ifld(1),'  ifld(2) =',ifld(2)
c
c                      special for grib
c                         float output of first points to anticipate
c                         grib floating point output
c
          fval1    = ifld(1)
          fdiff1   = ifld(2)
c
c       set first two points to second diff value for better packing
c
          ifld(1)  = ifld(3)
          ifld(2)  = ifld(3)
c  -----------------------------------------------------------
      return
      end
