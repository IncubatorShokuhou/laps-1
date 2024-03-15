      subroutine fi7513 (iwork,istart,npts,max,min,inrnge)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    fi7513      select block of data for packing
c   prgmmr: cavanaugh        org: w/nmc42    date: 94-01-21
c
c abstract: select a block of data for packing
c
c program history log:
c   94-01-21  cavanaugh
c   95-10-31  iredell     removed saves and prints
c
c usage:    call fi7513 (iwork,istart,npts,max,min,inrnge)
c   input argument list:
c     *        - return address if encounter set of same values
c     iwork    -
c     istart   -
c     npts     -
c
c   output argument list:      (including work arrays)
c     max      -
c     min      -
c     inrnge   -
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916/256, y-mp8/64, y-mp el92/256
c
c$$$
      integer        iwork(*),npts,istart,inrnge,inrnga,inrngb
      integer        max,min,mxval,maxb,minb,mxvalb
      integer        ibits(31)
c
      data           ibits/1,3,7,15,31,63,127,255,511,1023,2047,
     *               4095,8191,16383,32767,65535,131071,262143,
     *               524287,1048575,2097151,4194303,8388607,
     *               16777215,33554431,67108863,134217727,268435455,
     *               536870911,1073741823,2147483647/
c  ----------------------------------------------------------------
c                        identify next block of data for packing and
c                           return to caller
c  ********************************************************************
      istrta  = istart
c
c                     get block a
      call fi7516 (iwork,npts,inrnga,istrta,
     *                                  max,min,mxval,lwide)
c  ********************************************************************
c
      istrtb  = istrta + inrnga
 2000 continue
c                         if have processed all data, return
      if (istrtb.gt.npts) then
c                         no more data to look at
          inrnge  = inrnga
          return
      end if
c                     get block b
      call fi7502 (iwork,istrtb,npts,isame)
      if (isame.ge.15) then
c         print *,'block b has all identical values'
c         print *,'block a has inrnge =',inrnga
c                     block b contains all identical values
          inrnge  = inrnga
c                     exit with block a
          return
      end if
c                     get block b
c
      istrtb  = istrta + inrnga
      call fi7516 (iwork,npts,inrngb,istrtb,
     *                                  maxb,minb,mxvalb,lwideb)
c     print *,'block a',inrnga,' block b',inrngb
c  ********************************************************************
c                     perform trend analysis to determine
c                     if data collection can be improved
c
      ktrnd  = lwide - lwideb
c     print *,'trend',lwide,lwideb
      if (ktrnd.le.0) then
c         print *,'block a - smaller, should extend into block b'
          mxval   = ibits(lwide)
c
c                     if block a requires the same or fewer bits
c                             look ahead
c                        and gather those data points that can
c                        be retained in block a
c                        because this block of data
c                            uses fewer bits
c
          call fi7518 (iret,iwork,npts,istrta,inrnga,inrngb,
     *                          max,min,lwide,mxval)
          if(iret.eq.1) go to 8000
c         print *,'18 inrnga is now ',inrnga
          if (inrngb.lt.20) then
              return
          else
              go to 2000
          end if
      else
c         print *,'block a - larger, b should extend back into a'
          mxvalb  = ibits(lwideb)
c
c                     if block b requires fewer bits
c                             look back
c                            shorten block a because next block of data
c                            uses fewer bits
c
          call fi7517 (iret,iwork,npts,istrtb,inrnga,
     *                               maxb,minb,lwideb,mxvalb)
          if(iret.eq.1) go to 8000
c         print *,'17 inrnga is now ',inrnga
      end if
c
c                           pack up block a
c                           updata pointers
 8000 continue
      inrnge  = inrnga
c                           get next block a
 9000 continue
      return
      end
