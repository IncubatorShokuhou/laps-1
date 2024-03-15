      subroutine pk_sect7(kfildo,ia,nxy,ipack,nd5,locn,ipos,
     1                    is5,ns5,is6,ns6,is7,ns7,inc,minpk,missp,misss,
     2                    ipkopt,locn5_20,ipos5_20,locn5_32,
     3                    ipos5_32,l3264b,ier,isevere,*)
c
c        march    2000   glahn   tdl   for grib2
c        january  2001   glahn   comments; added check on size of
c                                is5( ) and on is6(5) = 6; removed
c                                is7(5) = 7; added ier = 703 for
c                                a type of packing error
c        november 2001   glahn   a few diagnostic format changes;
c                                split test on ns7 and is7(5)
c        january  2002   glahn   set section number = 99 in packed
c                                data when severe error occurs
c        february 2002   glahn   added sections 6,7 print for simple
c        march    2002   galhn   added non fatal error returns via
c                                iersav; added ier=0
c
c        purpose
c            packs data section 7, the data section, of a grib2
c            message.
c
c        data set use
c           kfildo - unit number for output (print) file. (output)
c
c        variables
c              kfildo = unit number for output (print) file. (input)
c               ia(j) = the data to pack (j=1,nxy).  (input)
c                 nxy = the size of ia( ).  (input)
c            ipack(j) = the array that holds the actual packed message
c                       (j=1,nd5). (input/output)
c                 nd5 = the size of the array ipack( ). (input)
c                locn = the word position to place the next value.
c                       (input/output)
c                ipos = the bit position in locn to start placing
c                       the next value. (input/output)
c              is5(j) = contains the grid definition data that
c                       will be packed into ipack( ) (j=1,ns5).
c                       (input/output)
c                 ns5 = size of is5( ). (input)
c              is6(j) = contains the bit map information 
c                       corresponding to section 6 of grib2 (j=1,ns6).
c                       (input)
c                 ns6 = size of is6( ). (input)
c              is7(j) = contains the grid definition data that
c                       will be packed into ipack( ) (j=1,ns7).
c                       (input/output)
c                 ns7 = size of is7( ). (input)
c                 inc = the increment to use in defining groups.
c                       (input)
c               minpk = the minimum group size.  (input)
c               missp = the primary missing value.  (input)
c               misss = the secondary missing value.  (input)
c              ipkopt = packing indicator:
c                       0 = error, don't pack
c                       1 = pack ia( ), simple
c                       2 = pack ia( ) and ib( ), simple
c                       3 = pack complex or spatial differencing
c                       4 = pack complex.
c                      (input)
c           locn5_20 = locn for octet 20 in section 5.  (input)
c           ipos5_20 = ipos for octet 20 in section 5.  (input)
c           locn5_32 = locn for octet 32 in section 5.  (input)
c           ipos5_32 = ipos for octet 32 in section 5.  (input)
c              l3264b = the integer word length in bits of the machine
c                       being used. values of 32 and 64 are
c                       accommodated. (input)
c                 ier = return status code. (output)
c                          0 = good return.
c                          1 = packing would overflow ipack( ).
c                          2 = ipos not in range 1 to l3264b.
c                          3 = nbit not in range 0 to 32.
c                          4 = nbit eq 0, but nvalue ne 0.
c                        701 = is7(5) does not indicate section 7.
c                        702 = is7( ) has not been dimensioned large enough
c                              to continue section 7.
c                        703 = not supported type of packing.
c                        711 = lbit incorrect.  returned from pk_cmplx.
c                        712 = incorrect splitting method.  returned
c                              from pk_cmplx.
c                        713 = unrecognized missing value flag
c                              in is5(23).  returned from pk_cmplx.
c             isevere = the severity level of the error.
c                       1 = non fatal error
c                       2 = a fatal error 
c                       (output)
c                   * = alternate return when ier ne 0.
c
c             local variables
c               cfeed = contains the character representation
c                       of a printer form feed.
c               ifeed = contains the integer value of a printer
c                       form feed.
c               izero = contains the value '0'.  this is used in the
c                       code simply to emphasize that a certain 
c                       group of octets in the message are being 
c                       zeroed out.
c                   n = l3264b = the integer word length in bits of
c                       the machine being used. values of 32 and
c                       64 are accommodated.
c                ibit = the number of bits to use to pack each value
c                       for simple packing.  for complex packing, ibit
c                       is the bits to use to pack the group minima.
c                       is5(20) set to ibit.
c               ifill = number of bits necessary to fill out an octet.
c         locn7,ipos7 = word and bit position of beginnning of
c                       section 7.
c        1         2         3         4         5         6         7 x
c
c        non system subroutines called
c           pkbg, pk_smple, pk_cmplx, length
c
      character*1 cfeed
c
      dimension ia(nd5),ipack(nd5),is5(ns5),is6(ns6),is7(ns7)
c
      data ifeed/12/
      data izero/0/
c
      ier=0
      iersav=0
      n=l3264b
      cfeed=char(ifeed)
c
      locn7=locn
      ipos7=ipos
c
c        all errors generated by this routine where the alternate 
c        return is used are fatal.  non fatal errors use the 
c        normal return.
      isevere=2
c        iersav can save an error return from another subroutine.
c
c        check to make sure that data has been provided for
c        section 7.
c
      if(ns7.lt.5)then
         ier=702
         go to 900
      else
c 
         if(is7(5).ne.7)then
            ier=701
            go to 900
         endif
c
c           bytes 1-4 must be filled in later with the record length
c           in bytes; above statement holds the place.  locn7 and 
c           ipos7 hold the location.
c
         call pkbg(kfildo,ipack,nd5,locn,ipos,izero,32,n,ier,*900)
c
         locns=locn
         iposs=ipos
c           save locn and ipos in case of error, section number is
c           set to 99.
         call pkbg(kfildo,ipack,nd5,locn,ipos,is7(5),8,n,ier,*900)
c
c           determine the packing method that we will
c           use.
c
         if(ipkopt.eq.1.or.ipkopt.eq.2)then
c
c              this is the simple method.  there are is5(6) values
c              to pack.
c
            call pk_smple(kfildo,ia,is5(6),ipack,nd5,locn,ipos,ibit,
     1                       n,ier,*900)
         elseif(ipkopt.eq.3.or.ipkopt.eq.4)then
c
c              this is the complex method or spatial difference.
c
            call pk_cmplx(kfildo,ia,nxy,ipack,nd5,locn,ipos,
     1                    is5,ns5,is7,ns7,inc,minpk,missp,
     2                    misss,ibit,locn5_32,ipos5_32,n,
     3                    ier,*900)
            iersav=ier
            ier=0
c              iersav can be used to provide error tracing.
c              subroutine reduce can produce non fatal errors.
c              other errors are fatal and return is to *900.
            
         else
            ier=703
c              not recognized type of packing.
            go to 900
         endif
c
c              fill is5(20) with ibit.  this can come from
c              either pk_smple or pk_cmplx.
c
         is5(20)=ibit
         call pkbg(kfildo,ipack,nd5,locn5_20,ipos5_20,ibit,8,n,ier,
     1             *900)
c
c           pad with zeros to fill out an octet, if necessary.
c
         ifill=mod(33-ipos,8)
c
         if(ifill.ne.0)then
            call pkbg(kfildo,ipack,nd5,locn,ipos,izero,ifill,n,ier,
     1                *900)
         endif
c
         is7(1)=length(kfildo,ipack,nd5,n,locn7,ipos7,locn,
     1                 ipos,ier)
c
c           write out the contents of section 5 and the non-data
c           contents of section 7.  separate simple, complex and
c           complex with second order differences
c
         if(is5(10).eq.0)then
c           write(kfildo,5)
c5          format(/' *******************************************'
c***d           write(kfildo,5)cfeed
c***d5          format(a1,/' *******************************************'
c    1                /' data values for sections 5, 6 - simple packing'
c    2                /' *******************************************')
c           write(kfildo,6)is5(1),is5(5),is5(6),is5(10),is5(12),is5(16),
c    1                     is5(18),is5(20),is5(21),is6(6),is6(1)
c6          format(/' 1:length of section ',t56,i10,
c    1      /' 5:number of section',t60,i6,
c    2      /' 6:number of actual data points ',t56,i10,
c    3      /' 10:data representation template number ',t60,i6,
c    4      /' 12:reference value ',t54,i12,
c    5      /' 16:binary scale factor ',t60,i6,
c    6      /' 18:decimal scale factor ',t60,i6,
c    7      /' 20:bits for each packed value',t60,i6,
c    8      /' 21:type of values ',t60,i6,
c    9      /'  6:bitmap indicator (secton 6) ',t60,i6,
c    a      /'  1:length of section 6 ',t56,i10)
c
c        elseif(is5(10).eq.2)then           
c           write(kfildo,10)
c10         format(/' *******************************************'
c***d           write(kfildo,10)cfeed
c***d10         format(a1,/' *******************************************'
c    1                /' data values for section 5 - complex packing'
c    3                /' *******************************************')
c           write(kfildo,20)is5(1),is5(5),is5(6),is5(10),is5(12),
c    1                      is5(16),is5(18),is5(20),is5(21),is5(22),
c    2                      is5(23),is5(24),is5(28),is5(32),is5(36),
c    3                      is5(37),is5(38),is5(42),is5(43),is5(47)
c20         format(/' 1:length of section ',t56,i10,
c    1      /' 5:number of section',t60,i6,
c    2      /' 6:number of actual data points ',t56,i10,
c    3      /' 10:data representation template number ',t60,i6,
c    4      /' 12:reference value ',t60,i6,
c    5      /' 16:binary scale factor ',t60,i6,
c    6      /' 18:decimal scale factor ',t60,i6,
c    7      /' 20:bits to pack group references',t60,i6,
c    8      /' 21:type of values ',t60,i6,
c    9      /' 22:splitting method ',t60,i6,
c    a      /' 23:use of missing values ',t60,i6,
c    b      /' 24:primary missing value ',t60,i6,
c    c      /' 28:secondary missing value ',t60,i6,
c    d      /' 32:number of groups ',t58,i8,
c    e      /' 36:reference for group widths ',t60,i6,
c    f      /' 37:bits to pack group widths ',t60,i6,
c    g      /' 38:reference for group lengths ',t60,i6,
c    h      /' 42:length increment for group lengths ',t60,i6,
c    i      /' 43:true length of last group. ',t58,i8,
c    j      /' 47:bits to pack group lengths ',t60,i6)
c
         elseif(is5(10).eq.3)then
c           write(kfildo,10)
c15         format(/' *******************************************'
c***d           write(kfildo,10)cfeed
c***d15         format(a1,/' *******************************************'
c    1                /' data values for section 5 - complex packing'
c    2                /' with second order differences'
c    3                /' *******************************************')
c           write(kfildo,20)is5(1),is5(5),is5(6),is5(10),is5(12),
c    1                      is5(16),is5(18),is5(20),is5(21),is5(22),
c    2                      is5(23),is5(24),is5(28),is5(32),is5(36),
c    3                      is5(37),is5(38),is5(42),is5(43),is5(47)
         endif
c
c           write is5(48) and is(49) that only pertain to second order
c           differencing.
c
c        if(is5(10).eq.3)then
c           write(kfildo,25)is5(48),is5(49)
c25         format(' 48:order of spatial differencing ',t60,i6,
c    1      /' 49:field width of spatial descriptors ',t60,i6)
c        endif
c
c        if(is5(10).eq.2)then
c           write(kfildo,30)
c30         format(/' *******************************************'
c***d           write(kfildo,30)cfeed
c***d30         format(a1,/' *******************************************'
c    1                /' data values for section 7 - complex packing'
c    3                /' *******************************************')
c           write(kfildo,40)is7(1),is7(5)
c40         format(/' 1:the number of octets in this field ',t57,i9,
c    1             /' 5:the number of this section',t60,i6)
c        elseif(is5(10).eq.3)then
c           write(kfildo,45)
c45         format(/' *******************************************'
c***d           write(kfildo,45)cfeed
c***d45         format(a1,/' *******************************************'
c    1                /' data values for section 7 - complex packing'
c    2                /' with second order differences'
c    3                /' *******************************************')
c           write(kfildo,50)is7(1),is7(5),is7(6),is7(7),is7(8)
c50         format(/' 1:the number of octets in this field ',t57,i9,
c    1             /' 5:the number of this section',t60,i6,
c    2             /' 6:the first original value in the field ',t58,i8,
c    3             /' 7:the second original value in the field ',t58,i8,
c    4             /' 8:the overall min of the 2nd order diff ',t58,i8)
c        elseif(is5(10).eq.0)then
c           write(kfildo,47)
c47         format(/' *******************************************'
c***d           write(kfildo,47)cfeed
c***d47         format(a1,/' *******************************************'
c    1                /' data values for section 7 - simple packing'
c    3                /' *******************************************')
c           write(kfildo,55)is7(1),is7(5)
c55         format(/' 1:the number of octets in this field ',t57,i9,
c    1             /' 5:the number of this section',t60,i6)
c
c        endif
c
      endif
c
c         error return section
 900  if(ier.ne.0)then
c
         if(isevere.eq.2)then
c              for a severe error, the section number is set = 99
c              so that unpacking will not occur.
            iers=ier
            call pkbg(kfildo,ipack,nd5,locns,iposs,99,8,n,ier,*900)
            ier=iers
            return 1
         endif
c  
      endif
c
      if(iersav.ne.0)then
         ier=iersav
         isevere=1
c           this provides for non fatal error returns. 
      endif     
c
      return
      end
