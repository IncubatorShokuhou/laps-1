      subroutine fi7503 (iwork,ipfld,npts,ibdsfl,bds11,
     *           len,lenbds,pds,refnce,iscal2,kwide,igds)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    fi7501      row by row, col by col packing
c   prgmmr: cavanaugh        org: w/nmc42    date: 94-05-20
c
c abstract: perform row by row or column by column packing
c   generating all bds information.
c
c program history log:
c   93-08-06  cavanaugh
c   95-10-31  iredell     removed saves and prints
c
c usage:    call fi7503 (iwork,ipfld,npts,ibdsfl,bds11,
c    *           len,lenbds,pds,refnce,iscal2,kwide,igds)
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
      integer         ipfld(*),igds(*)
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
c                         build second order bit map in either
c                         row by row or col by col format
      if (iand(igds(13),32).ne.0) then
c                              column by column
          kout  = igds(4)
          kin   = igds(5)
c         print *,'column by column',kout,kin
      else
c                              row by row
          kout  = igds(5)
          kin   = igds(4)
c         print *,'row by row',kout,kin
      end if
      kbds(17)  = kout
      kbds(19)  = npts
c
c     do 4100 j = 1, npts, 53
c         write (6,4101) (iwork(k),k=j,j+52)
 4101     format (1x,25i4)
c         print *,' '
c4100 continue
c
c                             initialize bit map pointer
      ibmp2p = 0
c                             construct working bit map
      do 2000 i = 1, kout
          do 1000 j = 1, kin
              if (j.eq.1) then
                  call sbyte (isobmp,1,ibmp2p,1)
              else
                  call sbyte (isobmp,0,ibmp2p,1)
              end if
              ibmp2p  = ibmp2p + 1
 1000     continue
 2000 continue
      len  = ibmp2p / 32 + 1
c     call binary(isobmp,len)
c
c                       process outer loop of row by row or col by col
c
      kptr  = 1
      kbds(11)  = kwide
      do 6000 i = 1, kout
c                       in current row or col
c                              find first order value
          jptr  = kptr
          lowest  = iwork(jptr)
          do 4000 j = 1, kin
              if (iwork(jptr).lt.lowest) then
                  lowest = iwork(jptr)
              end if
              jptr  = jptr + 1
 4000     continue
c                            save first order value
          call sbyte (ifoval,lowest,ifoptr,kwide)
          ifoptr  = ifoptr + kwide
c         print *,'foval',i,lowest,kwide
c                            subtract first order value from other vals
c                                         getting second order values
          jptr  = kptr
          ibig  = iwork(jptr) - lowest
          do 4200 j = 1, kin
              iwork(jptr)  = iwork(jptr) - lowest
              if (iwork(jptr).gt.ibig) then
                  ibig  = iwork(jptr)
              end if
              jptr  = jptr + 1
 4200     continue
c                            how many bits to contain largest second
c                                         order value in segment
          call fi7505 (ibig,nwide)
c                            save bit width
          call sbyte (isowid,nwide,iwdptr,8)
          iwdptr  = iwdptr + 8
c         print *,i,'soval',ibig,' in',nwide,' bits'
c         write (6,4101) (iwork(k),k=kptr,kptr+52)
c                            save second order values of this segment
          do 5000 j = 0, kin-1
              call sbyte (isoval,iwork(kptr+j),isoptr,nwide)
              isoptr  = isoptr + nwide
 5000     continue
          kptr    = kptr + kin
 6000 continue
c  =======================================================
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
c     print *,'isowid',iwdptr,ix
c     call binary (isowid,ix)
c
c                     no secondary bit map

c                                         determine location for start
c                                         of first order packed values
      kbds(12)  = iptr / 8 + 12
c                                        store location
      call sbyte (ipfld,kbds(12),0,16)
c                                     move in first order packed values
      ipass   = (ifoptr + 32) / 32
      call sbytes (ipfld,ifoval,iptr,32,0,ipass)
      iptr  = iptr + ifoptr
c     print *,'ifoval',ifoptr,ipass,kwide
c     call binary (ifoval,ipass)
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
c     print *,'isoval',isoptr,ix
c     call binary (isoval,ix)
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
      klen  = lenbds / 4 + 1
c     print *,'bds11 listing',4,lenbds
c     call binary (bds11,4)
c     print *,'ipfld listing'
c     call binary (ipfld,klen)
      return
      end
