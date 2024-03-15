      subroutine fi7516 (iwork,npts,inrng,istart,max,min,mxval,lwidth)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    fi7516      scan number of points
c   prgmmr: cavanaugh        org: w/nmc42    date: 94-01-21
c
c abstract: scan forward from current position. collect points and
c           determine maximum and minimum values and the number
c           of points that are included. forward search is terminated
c           by encountering a set of identical values, by reaching
c           the number of points selected or by reaching the end
c           of data.
c
c program history log:
c   94-01-21  cavanaugh
c   95-10-31  iredell     removed saves and prints
c
c usage:    call fi7516 (iwork,npts,inrng,istart,max,min,mxval,lwidth)
c   input argument list:
c     *        - return address if encounter set of same values
c     iwork    - data array
c     npts     - number of points in data array
c     istart   - starting location in data
c
c   output argument list:      (including work arrays)
c     inrng    - number of points selected
c     max      - maximum value of points
c     min      - minimum value of points
c     mxval    - maximum value that can be contained in lwidth bits
c     lwidth   - number of bits to contain max diff
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916/256, y-mp8/64, y-mp el92/256
c
c$$$
      integer        iwork(*),npts,istart,inrng,max,min,lwidth,mxval
      integer        ibits(31)
c
      data           ibits/1,3,7,15,31,63,127,255,511,1023,2047,
     *               4095,8191,16383,32767,65535,131071,262143,
     *               524287,1048575,2097151,4194303,8388607,
     *               16777215,33554431,67108863,134217727,268435455,
     *               536870911,1073741823,2147483647/
c  ----------------------------------------------------------------
c
      inrng  = 1
      jq        = istart + 19
      max       = iwork(istart)
      min       = iwork(istart)
      do 1000 i = istart+1, jq
          call fi7502 (iwork,i,npts,isame)
          if (isame.ge.15) then
              go to 5000
          end if
          inrng  = inrng + 1
          if (iwork(i).gt.max) then
              max  = iwork(i)
          else if (iwork(i).lt.min) then
              min  = iwork(i)
          end if
 1000 continue
 5000 continue
      krng   = max - min
c
      do 9000 lwidth = 1, 31
          if (krng.le.ibits(lwidth)) then
c             print *,'returned',inrng,' values'
              return
          end if
 9000 continue
      return
      end
