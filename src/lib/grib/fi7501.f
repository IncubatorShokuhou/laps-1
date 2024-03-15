      subroutine fi7501 (iwork,ipfld,npts,ibdsfl,bds11,
     *           len,lenbds,pds,refnce,iscal2,kwide)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    fi7501      bds second order packing
c   prgmmr: cavanaugh        org: w/nmc42    date: 93-08-06
c
c abstract: perform secondary packing on grid point data,
c   generating all bds information.
c
c program history log:
c   93-08-06  cavanaugh
c   93-12-15  cavanaugh   corrected location of start of first order
c                         values and start of second order values to
c                         reflect a byte location in the bds instead
c                         of an offset.
c   95-10-31  iredell     removed saves and prints
c
c usage:    call fi7501 (iwork,ipfld,npts,ibdsfl,bds11,
c    *           len,lenbds,pds,refnce,iscal2,kwide)
c   input argument list:
c     iwork    - integer source array
c     npts     - number of points in iwork
c     ibdsfl   - flags
c
c   output argument list:      (including work arrays)
c     ipfld    - contains bds from byte 12 on
c     bds11    - contains first 11 bytes for bds
c     len      - number of bytes from 12 on
c     lenbds   - total length of bds
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916/256, y-mp8/64, y-mp el92/256
c
c$$$
      character*1     bds11(*),pds(*)
c
      real            refnce
c
      integer         iscal2,kwide
      integer         lenbds
      integer         ipfld(*)
      integer         len,kbds(22)
      integer         iwork(*)
c                        octet number in section, first order packing
c     integer         kbds(12)
c                        flags
      integer         ibdsfl(*)
c                        extended flags
c     integer         kbds(14)
c                        octet number for second order packing
c     integer         kbds(15)
c                        number of first order values
c     integer         kbds(17)
c                        number of second order packed values
c     integer         kbds(19)
c                        width of second order packing
      integer         isowid(50000)
c                        secondary bit map
      integer         isobmp(8200)
c                        first order packed values
      integer         ifoval(50000)
c                        second order packed values
      integer         isoval(100000)
c
c     integer         kbds(11)
c                        bit width table
      integer         ibits(31)
c
      data            ibits/1,3,7,15,31,63,127,255,511,1023,
     *                      2047,4095,8191,16383,32767,65535,131072,
     *                      262143,524287,1048575,2097151,4194303,
     *                      8388607,16777215,33554431,67108863,
     *                      134217727,268435455,536870911,
     *                      1073741823,2147483647/
c  ----------------------------------
c                       initialize arrays
      do 100 i = 1, 50000
          isowid(i)  = 0
          ifoval(i)  = 0
  100 continue
c
      do 101 i = 1, 8200
          isobmp(i)  = 0
  101 continue
      do 102 i = 1, 100000
          isoval(i)  = 0
  102 continue
c                      initialize pointers
c                            secondary bit width pointer
      iwdptr  = 0
c                            secondary bit map pointer
      ibmp2p  = 0
c                            first order value pointer
      ifoptr  = 0
c                            byte pointer to start of 1st order values
      kbds(12)  = 0
c                            byte pointer to start of 2nd order values
      kbds(15)  = 0
c                            to contain number of first order values
      kbds(17)  = 0
c                            to contain number of second order values
      kbds(19)  = 0
c                            second order packed value pointer
      isoptr  = 0
c  =======================================================
c
c                         data is in iwork
c
      kbds(11)  = kwide
c
c       data packing
c
      iter    = 0
      inext   = 1
      istart  = 1
c  -----------------------------------------------------------
      kount = 0
c     do 1 i = 1, npts, 10
c         print *,i,(iwork(k),k=i, i+9)
c   1 continue
 2000 continue
      iter  = iter + 1
c     print *,'next iteration starts at',istart
       if (istart.gt.npts) then
           go to 4000
       else if (istart.eq.npts) then
           kpts    = 1
           mxdiff  = 0
           go to 2200
       end if
c
c                     look for repititions of a single value
       call fi7502 (iwork,istart,npts,isame)
       if (isame.ge.15) then
           kount = kount + 1
c          print *,'fi7501 - found identical set of ',isame
           mxdiff  = 0
           kpts    = isame
       else
c
c                     look for sets of values in trend selected range
           call fi7513 (iwork,istart,npts,nmax,nmin,inrnge)
c          print *,'istart  ',istart,' inrnge',inrnge,nmax,nmin
           iend  = istart + inrnge - 1
c          do 2199 nm = istart, iend, 10
c              print *,'  ',(iwork(nm+jk),jk=0,9)
c2199      continue
           mxdiff  = nmax - nmin
           kpts    = inrnge
       end if
 2200 continue
c     print *,'                 range ',mxdiff,' max',nmax,' min',nmin
c                 increment number of first order values
      kbds(17)  = kbds(17) + 1
c                 enter first order value
      if (mxdiff.gt.0) then
          do 2220 lk = 0, kpts-1
              iwork(istart+lk)  = iwork(istart+lk) - nmin
 2220     continue
          call sbyte (ifoval,nmin,ifoptr,kbds(11))
      else
          call sbyte (ifoval,iwork(istart),ifoptr,kbds(11))
      end if
      ifoptr  = ifoptr + kbds(11)
c                  process second order bit width
      if (mxdiff.gt.0) then
          do 2330 kwide = 1, 31
              if (mxdiff.le.ibits(kwide)) then
                  go to 2331
              end if
 2330     continue
 2331     continue
      else
          kwide  = 0
      end if
      call sbyte (isowid,kwide,iwdptr,8)
      iwdptr  = iwdptr + 8
c         print *,kwide,' ifoval=',nmin,iwork(istart),kpts
c               if kwide ne 0, save second order value
      if (kwide.gt.0) then
          call sbytes (isoval,iwork(istart),isoptr,kwide,0,kpts)
          isoptr  = isoptr + kpts * kwide
          kbds(19)  = kbds(19) + kpts
c         print *,'            second order values'
c         print *,(iwork(istart+i),i=0,kpts-1)
      end if
c                 add to second order bitmap
      call sbyte (isobmp,1,ibmp2p,1)
      ibmp2p  = ibmp2p + kpts
      istart  = istart + kpts
      go to 2000
c  --------------------------------------------------------------
 4000 continue
c     print *,'there were ',iter,' second order groups'
c     print *,'there were ',kount,' strings of constants'
c                 concantenate all fields for bds
c
c                   remainder goes into ipfld
      iptr  = 0
c                               bytes 12-13
c                                          value for n1
c                                          leave space for this
      iptr   = iptr + 16
c                               byte 14
c                                          extended flags
      call sbyte (ipfld,ibdsfl(5),iptr,1)
      iptr  = iptr + 1
      call sbyte (ipfld,ibdsfl(6),iptr,1)
      iptr  = iptr + 1
      call sbyte (ipfld,ibdsfl(7),iptr,1)
      iptr  = iptr + 1
      call sbyte (ipfld,ibdsfl(8),iptr,1)
      iptr  = iptr + 1
      call sbyte (ipfld,ibdsfl(9),iptr,1)
      iptr  = iptr + 1
      call sbyte (ipfld,ibdsfl(10),iptr,1)
      iptr  = iptr + 1
      call sbyte (ipfld,ibdsfl(11),iptr,1)
      iptr  = iptr + 1
      call sbyte (ipfld,ibdsfl(12),iptr,1)
      iptr  = iptr + 1
c                               bytes 15-16
c                 skip over value  for n2
      iptr  = iptr + 16
c                               bytes 17-18
c                                     p1
      call sbyte (ipfld,kbds(17),iptr,16)
      iptr  = iptr + 16
c                               bytes 19-20
c                                   p2
      call sbyte (ipfld,kbds(19),iptr,16)
      iptr  = iptr + 16
c                               byte 21 - reserved location
      call sbyte (ipfld,0,iptr,8)
      iptr  = iptr + 8
c                               bytes 22 - ?
c                                      widths of second order packing
      ix    = (iwdptr + 32) / 32
      call sbytes (ipfld,isowid,iptr,32,0,ix)
      iptr  = iptr + iwdptr
c                                      secondary bit map
      ij    = (ibmp2p + 32) / 32
      call sbytes (ipfld,isobmp,iptr,32,0,ij)
      iptr  = iptr + ibmp2p
      if (mod(iptr,8).ne.0) then
          iptr  = iptr + 8 - mod(iptr,8)
      end if
c                                         determine location for start
c                                         of first order packed values
      kbds(12)  = iptr / 8 + 12
c                                        store location
      call sbyte (ipfld,kbds(12),0,16)
c                                     move in first order packed values
      ipass   = (ifoptr + 32) / 32
      call sbytes (ipfld,ifoval,iptr,32,0,ipass)
      iptr  = iptr + ifoptr
      if (mod(iptr,8).ne.0) then
          iptr  = iptr + 8 - mod(iptr,8)
      end if
c     print *,'ifoptr =',ifoptr,' isoptr =',isoptr
c                determine location for start
c                     of second order values
      kbds(15)  = iptr / 8 + 12
c                                   save location of second order values
      call sbyte (ipfld,kbds(15),24,16)
c                  move in second order packed values
      ix    = (isoptr + 32) / 32
      call sbytes (ipfld,isoval,iptr,32,0,ix)
      iptr  = iptr + isoptr
      nleft  = mod(iptr+88,16)
      if (nleft.ne.0) then
          nleft  = 16 - nleft
          iptr   = iptr + nleft
      end if
c                                compute length of data portion
      len     = iptr / 8
c                                    compute length of bds
      lenbds  = len + 11
c  -----------------------------------
c                               bytes 1-3
c                                   this function completed below
c                                   when length of bds is known
      call sbyte (bds11,lenbds,0,24)
c                               byte  4
      call sbyte (bds11,ibdsfl(1),24,1)
      call sbyte (bds11,ibdsfl(2),25,1)
      call sbyte (bds11,ibdsfl(3),26,1)
      call sbyte (bds11,ibdsfl(4),27,1)
c                              enter number of fill bits
      call sbyte (bds11,nleft,28,4)
c                               byte  5-6
      if (iscal2.lt.0) then
          call sbyte (bds11,1,32,1)
          iscal2 = - iscal2
      else
          call sbyte (bds11,0,32,1)
      end if
      call sbyte (bds11,iscal2,33,15)
c
c$  fill octet 7-10 with the reference value
c   convert the floating point of your machine to ibm370 32 bit
c   floating point number
c                                        reference value
c                                          first test to see if
c                                          on 32 or 64 bit computer
          call w3fi01(lw)
          if (lw.eq.4) then
              call w3fi76 (refnce,iexp,imant,32)
          else
              call w3fi76 (refnce,iexp,imant,64)
          end if
          call sbyte (bds11,iexp,48,8)
          call sbyte (bds11,imant,56,24)
c
c                               byte  11
c
      call sbyte (bds11,kbds(11),80,8)
c
      return
      end
