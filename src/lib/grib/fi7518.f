      subroutine fi7518 (iret,iwork,npts,istrta,inrnga,inrngb,
     *                          maxa,mina,lwidea,mxvala)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    fi7518      scan forward
c   prgmmr: cavanaugh        org: w/nmc42    date: 94-01-21
c
c abstract: scan forward from start of block b towards end of block b
c           if next point under test forces a larger maxvala then
c           terminate indicating last point tested for inclusion
c           into block a.
c
c program history log:
c   94-01-21  cavanaugh
c   95-10-31  iredell     removed saves and prints
c   98-06-17  iredell     removed alternate return
c
c usage:    call fi7518 (iret,iwork,npts,istrta,inrnga,inrngb,
c     *                          maxa,mina,lwidea,mxvala)
c   input argument list:
c     ifld     -
c     jstart   -
c     npts     -
c
c   output argument list:      (including work arrays)
c     iret     -
c     jlast    -
c     max      -
c     min      -
c     lwidth   - number of bits to contain max diff
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916/256, y-mp8/64, y-mp el92/256
c
c$$$
      integer        iwork(*),npts,istrta,inrnga
      integer        maxa,mina,lwidea,mxvala
      integer        ibits(31)
c
      data           ibits/1,3,7,15,31,63,127,255,511,1023,2047,
     *               4095,8191,16383,32767,65535,131071,262143,
     *               524287,1048575,2097151,4194303,8388607,
     *               16777215,33554431,67108863,134217727,268435455,
     *               536870911,1073741823,2147483647/
c  ----------------------------------------------------------------
      iret=0
c     print *,'          fi7518'
      npos  = istrta + inrnga
      itst  = 0
c
 1000 continue
      itst  = itst + 1
      if (itst.le.inrngb) then
c         print *,'try npos',npos,iwork(npos),maxa,mina
          if (iwork(npos).gt.maxa) then
              if ((iwork(npos)-mina).gt.mxvala) then
c                 print *,'fi7518a -',itst,' range exceeds max'
                  iret=1
                  return
              else
                  maxa    = iwork(npos)
              end if
          else if (iwork(npos).lt.mina) then
              if ((maxa-iwork(npos)).gt.mxvala) then
c                 print *,'fi7518b -',itst,' range exceeds max'
                  iret=1
                  return
              else
                  mina    = iwork(npos)
              end if
          end if
          inrnga  = inrnga + 1
c         print *,'               ',itst,inrnga
          npos  = npos +1
          go to 1000
      end if
c  ----------------------------------------------------------------
 9000 continue
      return
      end
