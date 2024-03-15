      subroutine fi7517 (iret,iwork,npts,istrtb,inrnga,
     *                           maxb,minb,mxvalb,lwideb)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    fi7517      scan backward
c   prgmmr: cavanaugh        org: w/nmc42    date: 94-01-21
c
c abstract: scan backwards until a value exceeds range of group b
c           this may shorten group a
c
c program history log:
c   94-01-21  cavanaugh
c   95-10-31  iredell     removed saves and prints
c   98-06-17  iredell     removed alternate return
c
c usage:    call fi7517 (iret,iwork,npts,istrtb,inrnga,
c    *                           maxb,minb,mxvalb,lwideb)
c   input argument list:
c     iwork    -
c     istrtb   -
c     npts     -
c     inrnga   -
c
c   output argument list:      (including work arrays)
c     iret     -
c     jlast    -
c     maxb     -
c     minb     -
c     lwidth   - number of bits to contain max diff
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916/256, y-mp8/64, y-mp el92/256
c
c$$$
      integer        iwork(*),npts,istrtb,inrnga
      integer        maxb,minb,lwideb,mxvalb
      integer        ibits(31)
c
      data           ibits/1,3,7,15,31,63,127,255,511,1023,2047,
     *               4095,8191,16383,32767,65535,131071,262143,
     *               524287,1048575,2097151,4194303,8388607,
     *               16777215,33554431,67108863,134217727,268435455,
     *               536870911,1073741823,2147483647/
c  ----------------------------------------------------------------
      iret=0
c     print *,'          fi7517'
      npos  = istrtb - 1
      itst  = 0
      kset  = inrnga
c
 1000 continue
c     print *,'try npos',npos,iwork(npos),maxb,minb
      itst  = itst + 1
      if (itst.le.kset) then
          if (iwork(npos).gt.maxb) then
              if ((iwork(npos)-minb).gt.mxvalb) then
c                 print *,'went out of range at',npos
                  iret=1
                  return
              else
                  maxb    = iwork(npos)
              end if
          else if (iwork(npos).lt.minb) then
              if ((maxb-iwork(npos)).gt.mxvalb) then
c                 print *,'went out of range at',npos
                  iret=1
                  return
              else
                  minb    = iwork(npos)
              end if
          end if
          inrnga  = inrnga - 1
          npos  = npos - 1
          go to 1000
      end if
c  ----------------------------------------------------------------
c
 9000 continue
      return
      end
