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
