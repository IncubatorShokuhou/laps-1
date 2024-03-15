      subroutine w3fi75 (ibitl,itype,itoss,fld,ifld,ibmap,ibdsfl,
     &  npts,bds11,ipfld,pfld,len,lenbds,iberr,pds,igds)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:  w3fi75        grib pack data and form bds octets(1-11)
c   prgmmr: farley           org: nmc421      date:94-11-22
c
c abstract: this routine packs a grib field and forms octets(1-11)
c   of the binary data section (bds).
c
c program history log:
c   92-07-10  m. farley   original author
c   92-10-01  r.e.jones   correction for field of constant data
c   92-10-16  r.e.jones   get rid of arrays fp and int
c   93-08-06  cavanaugh   added routines fi7501, fi7502, fi7503
c                         to allow second order packing in pds.
c   93-07-21 stackpole    assorted repairs to get 2nd diff pack in
c   93-10-28 cavanaugh    commented out nonoperational prints and
c                         write statements
c   93-12-15  cavanaugh   corrected location of start of first order
c                         values and start of second order values to
c                         reflect a byte location in the bds instead
c                         of an offset in subroutine fi7501.
c   94-01-27  cavanaugh   added igds as input argument to this routine
c                         and added pds and igds arrays to the call to
c                         w3fi82 to provide information needed for
c                         boustrophedonic processing.
c   94-05-25  cavanaugh   subroutine fi7503 has been added to provide
c                         for row by row or column by column second
c                         order packing.  this feature can be activated
c                         by setting ibdsfl(7) to zero.
c   94-07-08  cavanaugh   commented out print statements used for debug
c   94-11-22  farley      enlarged work arrays to handle .5degree grids
c   95-06-01  r.e.jones   correction for number of unused bits at end
c                         of section 4, in bds byte 4, bits 5-8.
c   95-10-31  iredell     removed saves and prints
c
c usage:    call w3fi75 (ibitl,itype,itoss,fld,ifld,ibmap,ibdsfl,
c    &              npts,bds11,ipfld,pfld,len,lenbds,iberr,pds,igds)
c   input argument list:
c     ibitl     - 0, computer computes packing length from power
c                    of 2 that best fits the data.
c                 8, 12, etc. computer rescales data to fit into
c                    set number of bits.
c     itype     - 0 = if input data is floating point (fld)
c                 1 = if input data is integer (ifld)
c     itoss     - 0 = no bit map is included (don't toss data)
c                 1 = toss null data according to ibmap
c     fld       - real array of data to be packed if itype=0
c     ifld      - integer array to be packed if itype=1
c     ibmap     - bit map supplied from user
c     ibdsfl    - integer array containing table 11 flag info
c                 bds octet 4:
c                 (1) 0 = grid point data
c                     1 = spherical harmonic coefficients
c                 (2) 0 = simple packing
c                     1 = second order packing
c                 (3) 0 = original data were floating point values
c                     1 = original data were integer values
c                 (4) 0 = no additional flags at octet 14
c                     1 = octet 14 contains flag bits 5-12
c                 (5) 0 = reserved - always set to 0
c                 (6) 0 = single datum at each grid point
c                     1 = matrix of values at each grid point
c                 (7) 0 = no secondary bit maps
c                     1 = secondary bit maps present
c                 (8) 0 = second order values have constant width
c                     1 = second order values have different widths
c     npts      - number of gridpoints in array to be packed
c     igds      - array of gds information
c
c   output argument list:
c     bds11     - first 11 octets of bds
c     pfld      - packed grib field
c     len       - length of pfld
c     lenbds    - length of bds
c     iberr     - 1, error converting ieee f.p. number to ibm370 f.p.
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916/256, y-mp8/64, y-mp el92/256
c
c$$$
c
      real            fld(*)
c     real            fwork(260000)
c
c     fwork can use dynamic allocation of memory on cray
c
      real            fwork(npts)
      real            rmin,refnce
c
      integer         ipfld(*)
      integer         ibdsfl(*)
      integer         ibmap(*)
      integer         ifld(*),igds(*)
c     integer         iwork(260000)
c
c     iwork can use dynamic allocation of memory on cray
c
      integer         iwork(npts)
c
      logical         const
c
      character * 1   bds11(11),pds(*)
      character * 1   pfld(*)
      character * 1   ciexp(8)
      character * 1   cimant(8)
c
      equivalence     (iexp,ciexp(1))
      equivalence     (imant,cimant(1))
c
c            1.0   pack the field.
c
c            1.1   toss data if bitmap being used,
c                  moving 'data' to work area...
c
      const = .false.
      iberr = 0
      iw    = 0
c
      if (itoss .eq. 1) then
        if (itype .eq. 0) then
          do 110 it=1,npts
            if (ibmap(it) .eq. 1) then
              iw = iw + 1
              fwork(iw) = fld(it)
            endif
  110     continue
          npts = iw
        else if (itype .eq. 1) then
          do 111 it=1,npts
            if (ibmap(it) .eq. 1) then
              iw = iw + 1
              iwork(iw) = ifld(it)
            endif
  111     continue
          npts = iw
        endif
c
c             else, just move data to work array
c
      else if (itoss .eq. 0) then
        if (itype .eq. 0) then
          do 112 it=1,npts
            fwork(it) = fld(it)
  112     continue
        else if (itype .eq. 1) then
          do 113 it=1,npts
            iwork(it) = ifld(it)
  113     continue
        endif
      endif
c
c            1.2   convert data if needed prior to packing.
c                  (integer to f.p. or f.p. to integer)
c     itype = 0...floating point data
c       ibitl = 0...pack in least # bits...convert to integer
c     itype = 1...integer data
c       ibitl > 0...pack in fixed # bits...convert to floating point
c
      if (itype .eq. 0 .and. ibitl .eq. 0) then
        do 120 if=1,npts
          iwork(if) = nint(fwork(if))
  120   continue
      else if (itype .eq. 1 .and. ibitl .ne. 0) then
        do 123 if=1,npts
          fwork(if) = float(iwork(if))
  123   continue
      endif
c
c            1.3   pack the data.
c
      if (ibdsfl(2).ne.0) then
c                                    second order packing
c
c            print*,'  doing second order packing...'
          if (ibitl.eq.0) then
c
c             print*,'    and variable bit packing'
c
c                           working with integer values
c                           since doing variable bit packing
c
              max  = iwork(1)
              min  = iwork(1)
              do 300 i = 2, npts
                  if (iwork(i).lt.min) then
                      min  = iwork(i)
                  else if (iwork(i).gt.max) then
                      max  = iwork(i)
                  end if
  300         continue
c                           extract minima
              do 400 i = 1, npts
c                 if (iwork(i).lt.0) then
c                     print *,'minima 400',i,iwork(i),npts
c                 end if
                  iwork(i)  = iwork(i) - min
  400         continue
              refnce  = min
              idiff   = max - min
c             print *,'reference value',refnce
c
c             write (6,fmt='(''  minima removed      = '',/,
c    &              10(3x,10i10,/))') (iwork(i),i=1,6)
c             write (6,fmt='(''  end of array  = '',/,
c    &              10(3x,10i10,/))') (iwork(i),i=npts-5,npts)
c
c                      find bit width of idiff
c
              call fi7505 (idiff,kwide)
c             print*,'  bit width for original data', kwide
              iscal2 = 0
c
c             multiplicative scale factor set to 1
c             in anticipation of possible use in glahn 2dn diff
c
              scal2 = 1.
c
          else
c
c             print*,'   and fixed bit packing, ibitl = ', ibitl
c                               fixed bit packing
c                               - length of field in ibitl
c                               - must be real data
c                            floating point input
c
              rmax  = fwork(1)
              rmin  = fwork(1)
              do 100 i = 2, npts
                  if (fwork(i).lt.rmin) then
                      rmin  = fwork(i)
                  else if (fwork(i).gt.rmax) then
                      rmax  = fwork(i)
                  end if
  100         continue
              refnce  = rmin
c             print *,'100 reference',refnce
c                             extract minima
              do 200 i = 1, npts
                  fwork(i)  = fwork(i) - rmin
  200         continue
c             print *,'reference value',refnce
c             write (6,fmt='(''  minima removed      = '',/,
c    &              10(3x,10f8.2,/))') (fwork(i),i=1,6)
c             write (6,fmt='(''  end of array  = '',/,
c    &              10(3x,10f8.2,/))') (fwork(i),i=npts-5,npts)
c                                find largest delta
              idelt  = nint(rmax - rmin)
c                                do binary scaling
c                                   find out what binary scale factor
c                                       permits containment of
c                                       largest delta
              call fi7505 (idelt,iwide)
c
c                                   binary scaling
c
              iscal2  = iwide - ibitl
c             print *,'scaling needed to fit =',iscal2
c             print*,'  range of  = ',idelt
c
c                                expand data with binary scaling
c                                convert to integer
              scal2  = 2.0**iscal2
              scal2  = 1./ scal2
              do 600 i = 1, npts
                  iwork(i)  = nint(fwork(i) * scal2)
  600         continue
              kwide = ibitl
          end if
c
c  *****************************************************************
c
c           following is for glahn second differencing
c           not standard grib
c
c            test for second difference packing
c            based of size of pds - size in first 3 bytes
c
          call gbyte (pds,ipdsiz,0,24)
          if (ipdsiz.eq.50) then
c             print*,'  do second difference packing '
c
c                   glahn packing to 2nd diffs
c
c             write (6,fmt='(''  call to w3fi82 with = '',/,
c    &                  10(3x,10i6,/))') (iwork(i),i=1,npts)
c
               call w3fi82 (iwork,fval1,fdiff1,npts,pds,igds)
c
c             print *,'glahn',fval1,fdiff1
c             write (6,fmt='(''  out from w3fi82 with = '',/,
c    &                  10(3x,10i6,/))') (iwork(i),i=1,npts)
c
c             must now re-remove the minimum value
c             of the second differences to assure
c             all positive numbers for second order grib packing
c
c             original reference value added to first point
c             value from the 2nd diff packer to be added
c             back in when the 2nd diff values are
c             reconstructed back to the basic values
c
c             also, the reference value is
c             power-of-two scaled to match
c             fval1.  all of this scaling
c             will be removed after the
c             glahn second differencing is undone.
c             the scaling factor needed to do that
c             is saved in the pds as a signed positive
c             two byte integer
c
c             the scaling for the 2nd dif packed
c             values is properly set to zero
c
              fval1 = fval1 + refnce*scal2
c                                          first test to see if
c                                          on 32 or 64 bit computer
              call w3fi01(lw)
              if (lw.eq.4) then
                  call w3fi76 (fval1,iexp,imant,32)
              else
                  call w3fi76 (fval1,iexp,imant,64)
              end if
              call sbyte (pds,iexp,320,8)
              call sbyte (pds,imant,328,24)
c
              if (lw.eq.4) then
                  call w3fi76 (fdiff1,iexp,imant,32)
              else
                  call w3fi76 (fdiff1,iexp,imant,64)
              end if
              call sbyte (pds,iexp,352,8)
              call sbyte (pds,imant,360,24)
c
c             turn iscal2 into signed positive integer
c             and store in two bytes
c
              if(iscal2.ge.0)  then
                call sbyte (pds,iscal2,384,16)
              else
                call sbyte (pds,1,384,1)
                iscal2 = - iscal2
                call sbyte( pds,iscal2,385,15)
              endif
c
              max  = iwork(1)
              min  = iwork(1)
              do 700 i = 2, npts
                  if (iwork(i).lt.min) then
                      min  = iwork(i)
                  else if (iwork(i).gt.max) then
                      max  = iwork(i)
                  end if
  700         continue
c                           extract minima
              do 710 i = 1, npts
                  iwork(i)  = iwork(i) - min
  710         continue
              refnce  = min
c             print *,'710 reference',refnce
              iscal2 = 0
c
c             and reset value of kwide - the bit width
c             for the range of the values
c
              idiff = max - min
              call fi7505 (idiff,kwide)
c
c             print*,'bit width (kwide) of 2nd diffs', kwide
c
c  **************************** end of glahn packing  ************
          else if (ibdsfl(2).eq.1.and.ibdsfl(7).eq.0) then
c                        have second order packing with no second order
c                        bit map. ergo row by row - col by col
              call fi7503 (iwork,ipfld,npts,ibdsfl,bds11,
     *              len,lenbds,pds,refnce,iscal2,kwide,igds)
              return
          end if
c         write (6,fmt='(''  call to fi7501 with = '',/,
c    &                  10(3x,10i6,/))') (iwork(i),i=1,npts)
c         write (6,fmt='(''  end of array = '',/,
c    &                  10(3x,10i6,/))') (iwork(i),i=npts-5,npts)
c         print*,' refnce,iscal2, kwide at call to fi7501',
c    &             refnce, iscal2,kwide
c
c                         second order packing
c
          call fi7501 (iwork,ipfld,npts,ibdsfl,bds11,
     *             len,lenbds,pds,refnce,iscal2,kwide)
c
c              bds completely assembled in fi7501 for second order
c              packing.
c
      else
c                                      simple packing
c
c                print*,'  simple first order packing...'
          if (ibitl.eq.0) then
c                print*,' with variable bit length'
c
c                  with variable bit length, adjusted
c                  to accommodate largest value
c                  binary scaling always = 0
c
              call w3fi58(iwork,npts,iwork,pfld,nbits,len,kmin)
              rmin   = kmin
              refnce  = rmin
              iscale = 0
c             print*,'  bit length came out at ...',nbits
c
c           set const .true. if all values are the same
c
              if (len.eq.0.and.nbits.eq.0) const = .true.
c
          else
c           print*,' fixed bit length, ibitl = ', ibitl
c
c             fixed bit length packing (variable precision)
c             values scaled by power of 2 (iscale) to
c             fit largest value into given bit length (ibitl)
c
              call w3fi59(fwork,npts,ibitl,iwork,pfld,iscale,len,rmin)
              refnce = rmin
c             print *,' scaling needed to fit is ...', iscale
              nbits = ibitl
c
c           set const .true. if all values are the same
c
              if (len.eq.0) then
                  const = .true.
                  nbits = 0
              end if
          end if
c
c$        compute length of bds in octets
c
          inum  = npts * nbits + 88
c         print *,'number of bits before fill added',inum
c
c                  number of fill bits
          nfill  = 0
          nleft  = mod(inum,16)
          if (nleft.ne.0) then
              inum  = inum + 16 - nleft
              nfill = 16 - nleft
          end if
c         print *,'number of bits after fill added',inum
c                  length of bds in bytes
          lenbds = inum / 8
c
c                2.0   form the binary data section (bds).
c
c                 concantenate all fields for bds
c
c                               bytes 1-3
          call sbyte (bds11,lenbds,0,24)
c
c                               byte  4
c                                       flags
          call sbyte (bds11,ibdsfl(1),24,1)
          call sbyte (bds11,ibdsfl(2),25,1)
          call sbyte (bds11,ibdsfl(3),26,1)
          call sbyte (bds11,ibdsfl(4),27,1)
c                                        nr of fill bits
          call sbyte (bds11,nfill,28,4)
c
c$      fill octets 5-6 with the scale factor.
c
c                               byte  5-6
          if (iscale.lt.0) then
              call sbyte (bds11,1,32,1)
              iscale  = - iscale
              call sbyte (bds11,iscale,33,15)
          else
              call sbyte (bds11,iscale,32,16)
          end if
c
c$  fill octet 7-10 with the reference value
c   convert the floating point of your machine to ibm370 32 bit
c   floating point number
c
c                               byte  7-10
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
c
c$                        fill octet 11 with the number of bits.
c
c                               byte  11
          call sbyte (bds11,nbits,80,8)
      end if
c
      return
      end
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
      subroutine fi7505 (n,nbits)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    fi7505      determine number of bits to contain value
c   prgmmr: cavanaugh        org: w/nmc42    date: 93-06-23
c
c abstract: calculate number of bits to contain value n, with a
c            maximum of 32 bits.
c
c program history log:
c   93-06-23  cavanaugh
c   95-10-31  iredell     removed saves and prints
c
c usage:    call fi7505 (n,nbits)
c   input argument list:
c     n        - integer value
c
c   output argument list:      (including work arrays)
c     nbits    - number of bits to contain n
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916/256, y-mp8/64, y-mp el92/256
c
c$$$
      integer        n,nbits
      integer        ibits(31)
c
      data           ibits/1,3,7,15,31,63,127,255,511,1023,2047,
     *               4095,8191,16383,32767,65535,131071,262143,
     *               524287,1048575,2097151,4194303,8388607,
     *               16777215,33554431,67108863,134217727,268435455,
     *               536870911,1073741823,2147483647/
c  ----------------------------------------------------------------
c
      do 1000 nbits = 1, 31
          if (n.le.ibits(nbits)) then
              return
          end if
 1000 continue
      return
      end
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
